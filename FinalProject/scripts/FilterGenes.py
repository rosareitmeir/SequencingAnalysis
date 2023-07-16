from Bio import SeqIO

fasta_file = snakemake.input[0]
maxNs = 1- snakemake.params["min_coverage"]
kept_seqs=[]
kept_genes=[]
for record in SeqIO.parse(fasta_file, "fasta"):
    # Access the sequence ID
    sequence_id = record.id
    # Access the sequence
    sequence = record.seq
    Ncounts= sequence.count("N")
    GeneLength= len(sequence)
    if Ncounts/GeneLength <= maxNs:
        kept_seqs.append(record)
        kept_genes.append(record.id)


SeqIO.write(kept_seqs, snakemake.output[1], "fasta")

with open(snakemake.output[0], "w") as output:
    output.write("\n".join(kept_genes))