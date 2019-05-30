# Influenza-Mutation-Checker
A software tool designed to check influenza protein sequences for a user-defined set of mutations. This was originally designed to check for the changes found in the CDC H5N1 Genetic Changes Inventory (https://www.cdc.gov/flu/pdf/avianflu/h5n1-inventory.pdf), however this tool can be used to check for:

  - Residue substitutions
  - Single- or multi- residue deletions
  - Sequence motifs

USAGE:  
python mutation_checker.py [PATH_TO_REFERENCE_SEQUENCES] [PATH_TO_MUTATIONS_SPREADSHEET] [PATH_TO_INPUT_FILE] [PATH_FOR_OUTPUT]

The reference sequences used only define the numbering by which residue positions are defined. Reference sequences must be in FASTA format, with headers of the form:

\>[STRAIN_NAME] | [PROTEIN_NAME]
The reference sequences needed to check for the CDC mutations can be found in this repository as cdc2012_refs.fa.

The mutations spreadsheet should contain the mutations to be checked for. Multiple mutations placed on the same row of the spreadsheet will only be reported if all are found in a given strain. Some examples of how mutations can be encoded in this spreadsheet include:
  - PB2:627:K  (searches for a lysine at position 627 of PB2)
  - NA:10-30:del  (searches for full deletion of the region specified in NA, with respect to the reference provided)
  - NP:33:del  (searches for a single residue deletion at position 33 of NP)
  - HA:323-330:mot:R-X-R/K-R  (searches for the given motif anywhere within the region specified in HA)
The mutation spreadsheet should be in Excel format (.xlsx). The spreadsheet used to check for the CDC mutations can be found in this repository as cdc2012_mutations_list.xlsx.

The input file should contain protein sequences to be checked in FASTA format. Sequence headers should be in the same form as those for the reference sequences. All sequences should represent complete protein sequences, not partial.


DEPENDENCIES:

EMBOSS Needle

Pandas
