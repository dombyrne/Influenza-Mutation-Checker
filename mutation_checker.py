'''
Dominic Byrne - May 2019
Animal and Plant Health Agency

This tool checks a set of protein sequences for a user-defined set of SNPs,
deletions and sequence motifs, at positions relative to a reference 
protein sequence. It was originally designed to check influenza protein
sequences for the CDC H5N1 humanising mutations inventory, however
any set of user-defined features can be checked for.
'''

import sys, time, os, subprocess, re
import pandas as pd
import numpy as np

ref_file = sys.argv[1]
mut_file = sys.argv[2]
input_file = sys.argv[3]
output_file = sys.argv[4]


def parse_fasta_file(fasta_file):
	
	data = open(fasta_file).read()
	seqs = [seq for seq in data.split(">") if len(seq.strip()) > 0]
	seq_records = {seq.split("\n")[0] : "".join(seq.split("\n")[1:]).upper() 
				   for seq in seqs}
	prots = set([head.split("|")[1].strip() for head in seq_records.keys()])
	seqs_out = {prot : {} for prot in prots}

	for rec in seq_records.keys():
		
		prot_name = rec.split("|")[1].strip()
		strain = rec.split("|")[0].strip().replace(" ","_")
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



def find_mutations(row, aln_dict):
	
	#cols = out_frame.columns
	#strains = set([k[0] for k in aln_dict.keys()])
	#out_dict = {strain : [] for strain in strains}
	prots = set([k[1] for k in aln_dict.keys()])
	strain = row["STRAIN"]
	#return strain
	ref_name = row["REF_PROTEIN"]
	mutations = [str(i) for i in row["MUTATION":] if type(i) == unicode]
	results = []
	for mut in mutations:
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
		return True
	else:
		return False


	


ref_seqs = parse_fasta_file(ref_file)
query_seqs = parse_fasta_file(input_file)

prots_removed = []
for prot in query_seqs.keys():
	if prot not in ref_seqs.keys():
		prots_removed.append(prot)
		del query_seqs[prot]

if len(prots_removed) > 0:
	print "\nN.B: Sequences of the following proteins were not analysed as no relevant reference sequences were provided:"
	print "\n".join(prots_removed)

temp_dir = os.getcwd() + "/temp_mut_dir_" + str(time.time())
os.mkdir(temp_dir)
make_temp_seq_files(ref_seqs, query_seqs, temp_dir)
aligned_combos = align_seqs(ref_seqs, query_seqs, temp_dir)
alignments = get_alignments(ref_seqs, query_seqs, aligned_combos, temp_dir)
strains = list(set([k[0] for k in alignments.keys()]))
mut_df = pd.read_excel(mut_file)
len_df = len(mut_df.index)
out_df = pd.concat([mut_df]*len(strains), ignore_index=True)
out_df.insert(0, "STRAIN", np.nan)

for i, strain in enumerate(strains):
	out_df.loc[i*len_df : (i+1)*len_df - 1, "STRAIN"] = strain

mutations_found = out_df[out_df.apply(find_mutations, aln_dict = alignments, axis = 1)].set_index("STRAIN")
mutations_found.columns = [i if "Unnamed:" not in i else "" for i in mutations_found.columns]
mutations_found.to_excel(output_file)


#Comment out the last line if you want to view the alignments created
os.system("rm -r " + temp_dir)
