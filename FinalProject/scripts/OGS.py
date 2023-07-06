import pysam

class Gene:
    def __init__(self, id, chr,  start, end):
        self.id = id
        self.chr= chr
        self.start = int(start)
        self.end = int(end)

    # for the binary search
    def __lt__(self, other):
        return self.start < other.start

    def __str__(self):
        return f"Gene ID: {self.id}, Start: {self.start}, End: {self.end}"

    def __eq__(self, other):
        if isinstance(other, Gene):
            return self.id == other.id
        return False
    # for the dictionary
    def __hash__(self):
        return hash(self.id)

class Read:
    def __init__(self, id, sequence,  start, end):
        self.id = id
        self.sequence= sequence
        self.start = int(start)
        self.end = int(end)


def merge_sequences(sequences, gene_start, gene_end):
    # Sort the sequences based on start positions in ascending order
    sorted_sequences = sorted(sequences, key=lambda x: x.start)

    begin_gap = sorted_sequences[0].start - gene_start
    merged_sequence = "N" * begin_gap
    last_end = float("inf")


    for read in sorted_sequences:
        start = read.start
        end = read.end
        seq = read.sequence

        # Check if there is a gap between the previous sequence and the current one
        if start > last_end:
            # Add any gap sequence between the previous end position and the current start position
            gap_length = start - last_end
            merged_sequence += "N" * gap_length

        # Append the current sequence to the merged sequence
        merged_sequence += seq

        # Update the last end position
        last_end = end

    end_gap = gene_end - sorted_sequences[len(sorted_sequences) -1 ].end
    merged_sequence += "N" * end_gap

    return merged_sequence


# read gene tsv file and save them in a list
gene_input = "/home/rosa/SequencingAnalysis/Project2/results/featurecounts/gene_list2.tsv"  # snakemake.input["gene_list"]
genes = []
with open(gene_input, 'r') as file:
    for row in file:
        data = row.strip().split('\t')
        name, chr,  start, end = data
        gene = Gene(name, chr, start, end)
        genes.append(gene)

# Sort the genes by start position for binary search approach
sorted_genes = sorted(genes, key=lambda gene: gene.start)


# output file
output_path = snakemake.output["out"]
output= open( output_path, "w")
# Open the SAM file
bam_input = snakemake.input["bam"]
bamfile = pysam.AlignmentFile(bam_input, 'rb')

# setting the minimum coverage
min_coverage = snakemake.params["min_coverage"]
gene_dict={}

for gene in genes:
    reads = []
    #print(str(gene))
    read_iters = bamfile.fetch(gene.chr, gene.start, gene.end)
    for read in read_iters:
        if read.is_proper_pair and not read.is_secondary:
            #print(read.query_name+ "\t"+ str(read.reference_start) + "\t" + str(read.reference_end))
            seq= read.query_sequence
            if read.is_reverse:
                seq = seq[::-1].translate(str.maketrans("ATGC", "TACG"))
            reads.append(Read(read.query_name,seq,read.reference_start, read.reference_start + read.query_length))

    if len(reads) > 0:
        consensus = merge_sequences(reads, gene.start, gene.end)
        if consensus.count("N")/len(consensus) <= min_coverage:
            output.write( gene.id + "\t" + consensus + "\n")

# Close the SAM file
bamfile.close()
output.close()

