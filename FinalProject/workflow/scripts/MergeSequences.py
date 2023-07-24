import os.path

from Bio import SeqIO

fastas_input = snakemake.input["fastas"]


genes_input=snakemake.input["genes"]
with open(genes_input, 'r') as file:
    kept_genes=set([line.strip() for line in file.readlines()])

def get_merged_sequence (sample_fasta):
    kept_seqs = ""
    for record in SeqIO.parse(sample_fasta, "fasta"):
        if record.id in kept_genes:
            kept_seqs += str(record.seq)
    return (kept_seqs)

merged_entries=[]
for sample_fasta in fastas_input:
    ## get the merged sequence of all genes
    merged_seq = get_merged_sequence(sample_fasta)

    # create header for the fasta entry
    header= os.path.basename(sample_fasta) ## get file name from whole path
    header= header.split("_")[0] ## deriving the sample name
    # save it
    sample_entry= ">" + header + "\n" + merged_seq
    merged_entries.append(sample_entry)

with open(snakemake.output[0], 'w') as output:
    output.write('\n'.join(merged_entries))