#!/usr/bin/env python

import sys, time, os, subprocess, re
import pandas as pd


ref_file = sys.argv[1]
mut_file = sys.argv[2]
input_file = sys.argv[3]


def parse_fasta_file(fasta_file):
	
	data = open(fasta_file).read()
	seqs= [seq for seq in data.split(">") if len(seq.strip()) > 0]
	seq_records = {seq.split("\n")[0] : "".join(seq.split("\n")[1:]).upper() 
				   for seq in seqs}
	prots = set([head.split()[1] for head in seq_records.keys()])
	seqs_out = {prot : {} for prot in prots}
	
	
	for rec in seq_records.keys():
		
		prot_name = rec.split()[1].strip()
		strain = rec.split()[0]
		seq = seq_records[rec].strip()
		strain_seq = {strain : seq}
		seqs_out[prot_name].update(strain_seq)
	
	return seqs_out



def make_temp_seq_files(refs, queries, temp_dir_path):
	
	seqs = [refs, queries]
	for d in seqs:
		for prot_name in d.keys():
			for strain in d[prot_name].keys():
				seq_id = "_".join([strain.replace("/","_"), prot_name])
				head = ">" + seq_id
				seq = d[prot_name][strain]
				fname = "/".join([temp_dir_path, seq_id]) + ".fa"
				
				with open(fname, "w") as fh:
					fh.write(head + "\n")
					fh.write(seq + "\n")



def align_seqs(refs, queries, temp_dir_path):
	
	prot_combos = {prot : {} for prot in queries.keys()}
	for prot_name in queries.keys():
		combos = [(r_strain.replace("/", "_"), \
			  	  q_strain.replace("/", "_")) for \
				  r_strain in refs[prot_name].keys() for \
				  q_strain in queries[prot_name].keys()]
		prot_combos.update({prot_name : combos})
		
		for r_strain, q_strain in combos:
			r_file = temp_dir_path + "/" + "_".join([r_strain, prot_name]) + ".fa"
			q_file = temp_dir_path + "/" + "_".join([q_strain, prot_name]) + ".fa"
			outfile = temp_dir_path + "/" + "_".join([r_strain, "vs", q_strain, prot_name]) + ".aln"
			
			subprocess.check_output("needle -asequence %s -bsequence %s -outfile %s -gapopen 10 -gapextend 0.5 -aformat markx3" 
			 % (r_file, q_file, outfile), stderr = subprocess.STDOUT, shell=True)
	
	return prot_combos



def read_alignment(aln_file):
	
	data = open(aln_file).read()
	seqs = re.split(r'>.+\n', data)[1:]
	r_seq = seqs[0].replace("\n","").strip()
	q_seq = seqs[1].split("#")[0].replace("\n","").strip()
	
	return (r_seq, q_seq)



def get_alignments(refs, queries, prot_combos, temp_dir_path):
	
	q_strains = set([strain for prot in queries.keys() \
					 for strain in queries[prot].keys()])
	aln_dict = {}
	for strain in q_strains:
		s_id = strain.replace("/", "_")
		for prot_name in queries.keys():
			s_combos = [pair for pair in prot_combos[prot_name] if pair[1] == s_id]
			aln_fnames = [temp_dir_path + "/" + "_vs_".join(pair) + "_" + prot_name + ".aln" 
						  for pair in s_combos]
			
			aln_ref_ids = [r_strain for r_strain in refs[prot_name].keys() \
						   for pair in s_combos if \
						   r_strain.replace("/", "_") == pair[0]]
						   
			for ref_strain, fname in zip(aln_ref_ids, aln_fnames):
				aln_seqs = read_alignment(fname)
				aln_dict[(strain, prot_name, ref_strain)] = aln_seqs
	
	return aln_dict



def convert_motif_to_regex(motif):
	
	disam = {"B": "ND", "Z": "EQ", "J": "LI"}
	alt = []
	for i in motif:
		i = i.replace("/", "")
		for k in disam.keys():
			i = i.replace(k, disam[k])
		i = i.replace("X", "A-Z")
		if len(i) > 1:
			i = "[" + i + "]"
		alt.append(i)

	return "".join(alt)


#pos should be 0-indexed
def adjust_pos(pos, c, ref_seq):
	
	new_c = ref_seq.count("-", 0, pos + 1)

	if new_c - c == 0:
		
		return pos
	
	else:
		pos += new_c - c
		
		return adjust_pos(pos, new_c, ref_seq)



def find_mutations(aln_dict, mut_list):
	
	strains = set([k[0] for k in aln_dict.keys()])
	out_dict = {strain : [] for strain in strains}
	prots = set([k[1] for k in aln_dict.keys()])
	for strain in strains:
		for mut_info in mut_list:
			mut_info_s = [i for i in mut_info if type(i) != float]
			effect = mut_info_s[0]
			ref_name = mut_info_s[1]
			muts = mut_info_s[2:]
						
			results = []
			for mut in muts:
				f = mut.split(":")
				prot = f[0]
				pos = f[1]
				change = f[2]
				
				if prot not in prots:
					results.append(False)
					break
				
				r_seq = aln_dict[(strain, prot, ref_name)][0]
				q_seq = aln_dict[(strain, prot, ref_name)][1]
				
				if change == "mot":
					i_pos = [int(p) - 1 for p in pos.split("-")]
					pos_adj = [adjust_pos(p, 0, r_seq) for p in i_pos]
					region = q_seq[pos_adj[0] : pos_adj[1]]
					motif = convert_motif_to_regex(f[3].split("-"))
					motif_instances = re.findall(motif, region)
					if len(motif_instances) > 0:
						results.append(True)
					
				elif change == "del":
					if "-" in pos:
						i_pos = [int(p) - 1 for p in pos.split("-")]
						pos_adj = [adjust_pos(p, 0, r_seq) for p in i_pos]
						region = list(q_seq[pos_adj[0] : pos_adj[1]])
						if all(i == "-" for i in region):
							results.append(True)
					else:
						i_pos = int(pos) - 1
						pos_adj = adjust_pos(i_pos, 0, r_seq) 
						if q_seq[pos_adj] == "-":
							results.append(True)
						
				else:
					i_pos = int(pos) - 1
					pos_adj = adjust_pos(i_pos, 0, r_seq)
					if q_seq[pos_adj] == change:
						results.append(True)
					else:
						results.append(False)
					
			if False not in results and len(results) > 0:
				out_dict[strain].append((effect, muts))
	
	return out_dict
			


ref_seqs = parse_fasta_file(ref_file)
query_seqs = parse_fasta_file(input_file)

prots_removed = []
for prot in query_seqs.keys():
	if prot not in ref_seqs.keys():
		prots_removed.append(prot)
		del query_seqs[prot]

print "\nN.B: Sequences of the following proteins were not analysed as no relevant reference sequences were provided:"
print "\n".join(prots_removed)
		
temp_dir = os.getcwd() + "/temp_mut_dir_" + str(time.time())
os.mkdir(temp_dir)
make_temp_seq_files(ref_seqs, query_seqs, temp_dir)
aligned_combos = align_seqs(ref_seqs, query_seqs, temp_dir)
alignments = get_alignments(ref_seqs, query_seqs, aligned_combos, temp_dir)
mut_df = pd.read_excel(mut_file)
mutations = mut_df.values.tolist()
mutations_found = find_mutations(alignments, mutations)

'''Output'''
for strain in mutations_found.keys():
	muts = mutations_found[strain]
	print "\n" + strain
	for mut in muts:
		effect = mut[0]
		changes = mut[1]
		v_changes = []
		for change in changes:
			f = change.split(":")
			if f[2] == "mot":
				del f[2]
			v_changes.append(" ".join(f))
		str_changes = ", ".join(v_changes)
		print u'{:75.70}{:30.30}'.format(effect, str_changes)

os.system("rm -r " + temp_dir)
