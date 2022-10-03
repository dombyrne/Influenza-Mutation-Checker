'''
Dominic Byrne - May 2019
Animal and Plant Health Agency

This tool checks a set of protein sequences for a user-defined set of point mutations,
deletions and sequence motifs, at positions relative to a reference 
protein sequence. It was originally designed to check influenza protein
sequences for the CDC H5N1 humanising mutations inventory, however
any set of user-defined features can be checked for.
'''

import sys, os, subprocess, re, shutil, argparse, itertools
from pathlib import PurePath
from collections import Counter
import pandas as pd
import numpy as np

def load_mutations(f, mut_tb_cols, illegal_cols=['STRAIN']):
    #load the file
    tb = pd.read_excel(f)
    
    #check if illegal cols are present
    for col in illegal_cols:
        if col in tb.columns:
            raise Exception(f'The mutation table contained an illegal column name: {col}.\nThe following column names are not permitted in the mutation table:\n' + '\n'.join(illegal_cols))
            
    #check for missing mutations
    no_mutation = tb[mut_tb_cols['mut']].isnull()
    if no_mutation.any():
        raise Exception(f'No mutations found on the following rows of the mutations table:\n{tb[no_mutation]}')
    return tb

def load_ref_seqs(f, mutations, mut_tb_cols, ignore_missing_refs):
    #load sequences
    seqs = parse_fasta_file(f)
    
    #check we have the right reference sequences for all the mutations in the table
    if not ignore_missing_refs:
        print('Checking reference protein sequences...\n')
        check_refs_present(mutations, seqs, mut_tb_cols)
    return seqs

def load_query_seqs(f, mutations, mut_tb_cols):
    #load sequences
    seqs = parse_fasta_file(f)
    
    #check which proteins in which we have mutations to check for
    mut_prots = mutations[mut_tb_cols['mut']].str.split(':').apply(lambda x:x[0]).to_numpy()
    
    #compare list of proteins to query sequences
    will_check, wont_check = [], []
    for prot, strain in seqs.keys():
        seq_name = f'{strain}|{prot}'
        if prot in mut_prots:
            will_check.append(seq_name)
        else:
            wont_check.append(seq_name)
    print('The following protein sequences will be checked for mutations:')
    print('\n'.join(will_check) + '\n')
    if wont_check:
        print('These proteins will not be checked however, as they do not match any of the mutations provided:')
        print('\n'.join(wont_check) + '\n')
    return seqs  

def parse_fasta_file(f):
    #load the file
    with open(f) as fh:
        data = fh.read()
    
    #split data into chunks pertaining to individual sequences
    chunks = [chunk for chunk in data.split(">") if len(chunk.strip()) > 0]
    
    #loop through chunks and process individual sequences
    seqs = {}
    for chunk in chunks:
        lines = chunk.split('\n')
        header, seq = lines[0], ''.join(lines[1:]).upper()
        strain, prot = [f.strip() for f in header.split('|')]
            
        #store sequence in dict
        seqs[(prot, strain)] = seq       
    return seqs

def check_refs_present(mutations, ref_seqs, mut_cols):
    missing_refs = set()
    
    #loop through mutations
    for i in mutations.index:
        r_strain = mutations.loc[i, mut_cols['ref']]
        mutation = mutations.loc[i, mut_cols['mut']]
        prot = mutation.split(':')[0]
        
        #check if ref seq exists
        try:
            r_seq = ref_seqs[(prot, r_strain)]
        except KeyError:
            missing_refs.add(f'{r_strain}|{prot}')
    
    #raise exception if any ref seqs missing
    if missing_refs :
        excpt_str = '\n'.join((
            'Mutations present in the mutations table which require the following missing reference sequences:',
            '\n'.join(sorted(list(missing_refs))),
            'If these proteins are irrelevant and not in the input set, try --ignore_missing_refs.',
            'Otherwise, add the required reference sequences or edit the mutations file.'
        ))
        raise Exception(excpt_str)
    
def get_alignment_pairs(mutations, query_seqs, mut_tb_cols):
    
    #get all the required ref strain/protein combos we have mutations for
    r_strains = mutations[mut_tb_cols['ref']].to_numpy()
    r_prots = mutations[mut_tb_cols['mut']].str.split(':').apply(lambda x:x[0]).to_numpy()
    ref_prots_needed = set(list(zip(r_strains, r_prots)))
    
    #loop through query sequences and pair with relevant refs
    aln_pairs = []
    for q_prot, q_strain in query_seqs.keys():
            for r_strain, r_prot in ref_prots_needed:
                if r_prot == q_prot:
                    aln_pairs.append((q_prot, q_strain, r_strain))
    return aln_pairs 

def make_temp_seq_files(aln_pairs, refs, queries, temp_dir):
    
    #get unique list of sequence files that need ot be written
    all_prots = set()
    for prot, q_strain, r_strain in aln_pairs:
        all_prots.add((q_strain, prot))
        all_prots.add((r_strain, prot))
    
    #join sequences to single dict
    seqs = queries.copy()
    seqs.update(refs)
    
    #write each individual protein sequence to a fasta file
    strain_ids = {}
    for strain, prot in all_prots:
        #make filename-friendly version of the strain name
        try:
            strain_id = strain_ids[strain]
        except KeyError:
            strain_id = strain.replace('/', '_')
            strain_ids[strain] = strain_id

        #prepare strings for writing file
        header = f'>{strain_id}_{prot}'
        seq = seqs[(prot, strain)]
        file_str = f'{header}\n{seq}\n'

        #write to file
        outfile = PurePath(temp_dir, f'{strain_id}.{prot}.fa')
        with open(outfile, 'w') as fh:
            fh.write(file_str)
    return strain_ids

def get_alignments(aln_pairs, strain_ids, temp_dir, needle_go, needle_ge):
    
    aln_dict = {}
    for prot, q_strain, r_strain in aln_pairs:
            #get filename friendly strain ids
            q_strain_id, r_strain_id = strain_ids[q_strain], strain_ids[r_strain]
            
            #get the paths to the individual fasta file for each strain in the pair
            q_strain_file, r_strain_file = [PurePath(temp_dir, f'{strain_id}.{prot}.fa')
                                            for strain_id in (q_strain_id, r_strain_id)]
            
            #align the sequences
            aln_file = PurePath(temp_dir, f'{q_strain_id}_vs_{r_strain_id}.{prot}.aln')
            align_seqs(q_strain_file, r_strain_file, aln_file, needle_go, needle_ge)
                
            #read alignment
            aln_dict[(prot, q_strain, r_strain)] = read_alignment(aln_file)
    return aln_dict        
        
def align_seqs(afile, bfile, outfile, needle_go, needle_ge):
    #specify needle command
    cmd = ['needle', 
           '-asequence', str(afile), 
           '-bsequence', str(bfile), 
           '-outfile', str(outfile),
           '-gapopen', str(needle_go),
           '-gapextend', str(needle_ge),
           '-aformat', 'markx3'
          ]
    
    #run alignment
    try:
        subprocess.check_output(cmd, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise Exception(f'Needle failed, raising the following error:\n{e.stderr.decode()}')
    
def read_alignment(f):
    #load file
    with open(f) as fh:
        data = fh.read()
    
    #split on the sequence header lines
    seqs = re.split(r'>.+\n', data)
    
    #clean to get sequences only
    q_seq = seqs[1].replace('\n', '').strip()
    r_seq = seqs[2].split('#')[0].replace('\n', '').strip()
    return q_seq, r_seq

def find_mutations(mut_df, aln, mut_cols):

    strain_tbs = []
    
    #loop through each strain and check for all the mutations individually
    q_strains = list(set([k[1] for k in aln.keys()]))
    for strain in q_strains:
        strain_tbs.append(get_mutations_in_strain(strain, mut_df, aln, mut_cols))
    
    #concatenate the tables
    return pd.concat(strain_tbs)
    
def get_mutations_in_strain(strain, mut_df, aln, mut_cols):
    strain_mut_idx = []
    #loop through mutations
    for i in mut_df.index:
        #get the relevant reference strain which defines the positional numbering for the mutation in question
        r_strain = mut_df.loc[i, mut_cols['ref']]
        
        #extract the combination of mutation codes from the table
        mut_col_idx = list(mut_df.columns).index(mut_cols['mut'])
        mutations = mut_df.loc[i].to_numpy()[mut_col_idx:]
        mutations = mutations[~pd.isnull(mutations)]
        
        #check if the query sequence has all the mutations
        has_mutations = check_mutations(strain, mutations, r_strain, aln)
        if has_mutations:
            strain_mut_idx.append(i)
    
    #make output with all the mutations in this strain
    strain_tb = mut_df.loc[strain_mut_idx]
    strain_tb['STRAIN'] = strain
    
    #remove missing column names and reorder columns 
    strain_tb.columns = [i if "Unnamed:" not in i else "" for i in strain_tb.columns]
    n_cols = len(strain_tb.columns)
    new_col_idx = [n_cols-1] + [i for i in range(1, n_cols-1)]
    return strain_tb.iloc[:, new_col_idx]
    
def check_mutations(q_strain, mutations, r_strain, aln):
    
    #get the proteins for which we have sequences for the given query strain
    q_prots = set([k[0] for k in aln.keys()])
    
    #loop through all mutation in the group (row) and check individually
    for m in mutations:
        prot, pos_str, change, *motif = m.split(':')
        
        #check if we have the sequence for this protein
        if prot not in q_prots:
            return False
        
        #get the aligned protein sequences
        q_seq, r_seq = aln[(prot, q_strain, r_strain)]
        #parse position of the mutation
        pos = parse_pos(pos_str)
        try:
            start_pos, end_pos = pos
        except TypeError:
            start_pos, end_pos = pos, pos+1
        #correct positions for gaps in the ref seq
        start_pos, end_pos = [adjust_pos_for_ref_gaps(pos, r_seq) for pos in (start_pos, end_pos)]
        #slice the region of interest from the query seq
        region = q_seq[start_pos: end_pos]
        
        #deal with motif-style mutations
        if change == 'mot':
            #convert the encoded motif into a regex
            motif_regex = convert_motif_to_regex(motif[0])
            
            #use the regex to search for motif in roi
            motif_hits = re.findall(motif_regex, region)
            if motif_hits:
                continue
            else:
                return False
        #deal with deletions
        elif change =='del':
            #check if entire region is gaps - indicates deletion
            if all(aa == '-' for aa in region):
                continue
            else:
                return False
        #deal with point mutations
        else:
            #check if point mutation matches expected change
            if region == change:
                continue
            else:
                return False
    #if none of the mutations are missing then return true for all mutations present
    else:
        return True        
             
def parse_pos(pos):
    try:
        #single integers for point mutations
        return int(pos) - 1 #convert to 0-indexing
    
    except ValueError:
        #ranges for deletions/motifs
        start, end = [int(p) for p in pos.split('-')]
        return start-1, end #end stays 1-indexed for slicing
    
#this func ensures we're looking at the correct AA, as if there are gaps in the ref alignment before 
#our target AA we'll look in the wrong position in the query sequence for the mutation of interest
def adjust_pos_for_ref_gaps(pos, ref_seq, gap_count=0):
    
    #count the number of gaps in the ref alignment up to the position
    new_gap_count = ref_seq.count('-', 0, pos+1)
    count_change = new_gap_count - gap_count
    pos += count_change
    
    #if no change, return the position
    if count_change == 0:
        return pos
    
    #otherwise, add the gap_count and recursively call this function
    else:
        return adjust_pos_for_ref_gaps(pos, ref_seq, new_gap_count)
    
def convert_motif_to_regex(motif):
    #map for ambiguous amino acid codes
    disam = {'B': 'ND', 'Z': 'EQ', 'J': 'LI', 'X': 'A-Z'}
    out = []
    
    #loop through and disambiguate AA codes
    for aa in motif.split('-'):
        aa = aa.replace('/', '')
        for k, v in disam.items():
            aa = aa.replace(k, v)
            
        #add square brackets for groups of accepted AAs in a position
        if len(aa) > 1:
            aa = f'[{aa}]'
        out.append(aa)
    return ''.join(out)
 
#Get CLI
parser = argparse.ArgumentParser(description='Simple mutation checker script, originally designed to look for known mutations associated with increased virulence in Influenza A protein sequences. Takes a FASTA file of protein sequences as input, aligns them to reference protein sequences and checks user-defined positions for point mutations, deletions and sequence motifs. For more information see:\nhttps://github.com/dombyrne/Influenza-Mutation-Checker.')
parser.add_argument('-i', '--input', type=str, required=True, help='Input sequences (FASTA) to check for mutations.')
parser.add_argument('-r', '--references', type=str, required=True, 
                    help='Reference sequences (FASTA) to align query sequences to.')
parser.add_argument('-m', '--mutations', type=str, required=True, help='List of mutations to check for (Excel format).')
parser.add_argument('-o', '--output', type=str, required=True, 
                    help="Path to an output file in which to save results. If this path doesn't end in .xlsx then that extension will be added.")
parser.add_argument('--temp_dir', type=str, required=False, help='Path to a temp directory to store intermediate files. Default: ./temp.')
parser.add_argument('--keep_temp', action='store_true', help='Set this argument to retain the temp directory after the analysis finishes. Useful for inspecting alignment files. Default: False.')
parser.add_argument('--overwrite', action='store_true', help='Whether to overwrite existing output files. Will raise an error if files exist and this argument is not used. Default: False.')
parser.add_argument('--ignore_missing_refs', action='store_true', help='Will prevent an exception being raised when the reference protein sequence required to check for a mutation is not present. Useful if the missing references pertain to proteins not in the query set, but will not prevent an error being raised if the relevantr protein is included in the input. Default: False.')
parser.add_argument('--needle_ge', type=float, required=False, help='Gap extend parameter passed to needle for the sequence alignment. Default: 0.5')
parser.add_argument('--needle_go', type=float, required=False, help='Gap open parameter passed to needle for the sequence alignment. Default: 10')
parser.add_argument('--ref_col', type=str, required=False, help='Column name in the mutations table which describes the reference strain used to define the positional numbering for a given mutation. Default: REF_STRAIN.')
parser.add_argument('--mut_col', type=str, required=False, help='Name of the first mutations column in the mutations table. All columns to the right of this will be read as additional mutations regardless of column name. Mutations on the same row of the table will only be identified if all are present. Default: MUTATION.')
args = parser.parse_args()

#define default argument values
GE_DEFAULT = 0.5
GO_DEFAULT = 10
TEMP_DEFAULT = os.getcwd() + '/temp'
REF_COL_DEFAULT = 'REF_STRAIN'
MUT_COL_DEFAULT = 'MUTATION'

#set default values
if args.needle_ge is None:
    args.needle_ge = GE_DEFAULT
if args.needle_go is None:
    args.needle_go = GO_DEFAULT
if args.temp_dir is None:
    args.temp_dir = TEMP_DEFAULT
if args.ref_col is None:
    args.ref_col = REF_COL_DEFAULT
if args.mut_col is None:
    args.mut_col = MUT_COL_DEFAULT
    
#print options for logging purposes
print('Running analysis with the following parameters:')
for arg in vars(args):
    print(f'{arg}: {getattr(args, arg)}')

#update mutation table column values
mut_tb_cols = {
    'ref': args.ref_col,
    'mut': args.mut_col
}

#check if temp dir exists
args.temp_dir = PurePath(args.temp_dir)
if os.path.isdir(args.temp_dir):
    if args.overwrite:
        shutil.rmtree(args.temp_dir, ignore_errors=True)
    else:
        raise FileExistsError(f'Temp dir {args.temp_dir} already exists. Use --overwrite to ignore and overwrite.')
os.mkdir(args.temp_dir)

#check if output exists
if not args.output.endswith('.xlsx'):
    args.output += '.xlsx'
if os.path.isfile(args.output) and not args.overwrite:
    raise FileExistsError(f'Output file {args.output} already exists. Use --overwrite to ignore and overwrite.')

#read table of mutations to check for
mutations = load_mutations(args.mutations, mut_tb_cols)
print(f'\nLoaded {mutations.index.size} mutations.\n')

#load reference sequences
print('Loading reference sequences...')
ref_seqs = load_ref_seqs(args.references, mutations, mut_tb_cols, args.ignore_missing_refs)

#load query_sequences
query_seqs = load_query_seqs(args.input, mutations, mut_tb_cols)

#use the mutations table to discover which protein sequences need to be aligned
aln_pairs = get_alignment_pairs(mutations, query_seqs, mut_tb_cols)
                
#write separate fasta file for each protein sequence for seq alignment input
print('Preparing for alignment...')
strain_ids = make_temp_seq_files(aln_pairs, ref_seqs, query_seqs, args.temp_dir)

#do pairwise sequence alignments for each query/ref protein and read aligned sequences
print(f'Aligning {len(aln_pairs):,} pairs of protein sequences...\n')
alignments = get_alignments(aln_pairs, strain_ids, args.temp_dir, args.needle_go, args.needle_ge)

#check all the query protein sequences for the mutations
print('Checking sequences for mutations...')
mutations_found = find_mutations(mutations, alignments, mut_tb_cols)

#save results to file
print('Saving identified mutations to file...\n')
mutations_found.to_excel(args.output, index=False)

#remove tempfiles
if not args.keep_temp:
    print('Removing temp files...')
    shutil.rmtree(args.temp_dir, ignore_errors=True)