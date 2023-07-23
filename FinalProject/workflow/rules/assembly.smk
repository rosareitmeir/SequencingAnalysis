
####################################### QC and preprocessing for WGS data ######################################################################
rule wgs_fastqc:
	input:
		wgs_path + "_{number}.fastq.gz"
	output:
		html="results/qc/wgs/raw/" + wgs_name + "_{number}_fastqc.html",
		zip="results/qc/wgs/raw/" + wgs_name + "_{number}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file.
	log:
		"logs/fastqc/raw/"+ wgs_name + "_{number}.log"
	wrapper:
		"v1.31.1/bio/fastqc"


rule wgs_cutadapt:
	input:
		expand(wgs_path + "_{number}.fastq.gz", number =numbers)
	output:
		fastq1="results/trimmed/wgs/"+ wgs_name + "_1.fastq",
		fastq2="results/trimmed/wgs/"+ wgs_name + "_2.fastq",
		qc="results/qc/wgs/trimmed/" + wgs_name + ".qc.txt"
	params:
		adapters= config["software"]["cutadapt"]["adapters"],
		extra= config["software"]["cutadapt"]["extra"]
	log:
		"logs/cutadapt/"+ wgs_name + ".log"
	threads:
		config["software"]["cutadapt"]["threads"]
	wrapper:
		"v1.31.1/bio/cutadapt/pe"


rule wgs_fastqc_trimmed:
	input:
		"results/trimmed/wgs/"+ wgs_name + "_{number}.fastq"
	output:
		html="results/qc/wgs/trimmed/" + wgs_name + "{number}_fastqc.html",
		zip="results/qc/wgs/trimmed/" + wgs_name + "{number}_fastqc.zip"  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	log:
		"logs/fastqc_trimmed/" + wgs_name + "_{number}.log"
	wrapper:
		"v1.31.1/bio/fastqc"


#### Screening see screen.smk

#### Subsampling seqtk

def get_files_for_norm(wildcards):
	if os.path.exists(config["contaminations"]):
		return ["results/decontamination/" + wgs_name + "_1.fastq.gz",
				"results/decontamination/" + wgs_name + "_2.fastq.gz"]
	elif (config["WGS_trimmed"]):
		return ["results/trimmed/wgs/" + wgs_name + "_1.fastq", "results/trimmed/wgs/" + wgs_name + "_2.fastq"]
	else:
		return [wgs_path + "_1.fastq.gz", wgs_path + "_2.fastq.gz"]


#Normalization by subsampling reads

rule seqtk_subsample_pe:
	input:
		fq_files=get_files_for_norm
	output:
		f1="results/seqtk/"+ wgs_name + ".1.subsampled.fastq.gz",
		f2="results/seqtk/"+ wgs_name + ".2.subsampled.fastq.gz"
	params:
		n=config["software"]["seqtk"]["reads"],
		seed=12345  #seed to make sure both paired files are subsampled the same way
	log:
		"logs/seqtk_subsample/" + wgs_name + ".log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["seqtk"]["threads"]
	shell:
		"""
        seqtk sample -s {params.seed} {input.fq_files[0]} {params.n} | pigz -9 -p {threads} > {output.f1} && \
        seqtk sample -s {params.seed} {input.fq_files[1]} {params.n} | pigz -9 -p {threads} > {output.f2} 2> {log}
        """


############################################################ Assembly Part ################################################################

rule assembly:
	input:
	#normalized samples from seqtk
		"results/seqtk/" + wgs_name + ".1.subsampled.fastq.gz",
		"results/seqtk/" + wgs_name + ".2.subsampled.fastq.gz"
	params:
		output_dir="results/assembly/megahit/" + wgs_name,
		out_prefix= wgs_name
	output:
		"results/assembly/megahit/" + wgs_name + "/checkpoints.txt",
		"results/assembly/megahit/" + wgs_name + "/done",
		directory("results/assembly/megahit/" + wgs_name + "/intermediate_contigs"),
		"results/assembly/megahit/" + wgs_name + "/options.json",
		"results/assembly/megahit/" + wgs_name + "/" + wgs_name + ".contigs.fa",
		"results/assembly/megahit/" + wgs_name + "/" + wgs_name + ".log"
	log:
		"logs/megahit/" + wgs_name + ".log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["megahit"]["threads"]
	shell:
		"""
		rm -rf {params.output_dir}
		megahit -1 {input[0]} -2 {input[1]} -t {threads} -o {params.output_dir} --out-prefix {params.out_prefix} &> {log}
		"""

### Polishing

rule bowtie_assembly_index:
	input:
		reference = "results/assembly/megahit/" + wgs_name + "/" + wgs_name + ".contigs.fa"
	output:
		expand("results/assembly/megahit/"+ wgs_name + "/" + wgs_name +  ".{index}.bt2",index=range(1,5)),
		expand("results/assembly/megahit/"+ wgs_name + "/" + wgs_name +  ".rev.{index}.bt2",index=range(1,3))
	params:
		"results/assembly/megahit/" + wgs_name + "/" + wgs_name
	log:
		"logs/bowtie2/build_assembly_index/" + wgs_name + ".log"
	conda:
		"../envs/env.yaml"
	shell:
		"bowtie2-build -f {input.reference} {params} &> {log}"


# Map WGS reads back to assembly

rule bowtie_assembly_map:
	input:
		"results/seqtk/" + wgs_name + ".1.subsampled.fastq.gz",
		"results/seqtk/" + wgs_name + ".2.subsampled.fastq.gz",
		expand("results/assembly/megahit/" + wgs_name + "/" + wgs_name + ".{index}.bt2",index=range(1,5)),
		expand("results/assembly/megahit/" + wgs_name + "/" + wgs_name + ".rev.{index}.bt2",index=range(1,3))
	params:
		"results/assembly/megahit/" + wgs_name + "/" + wgs_name
	output:
		"results/assembly/mapping/sam/"+ wgs_name + ".sam"
	log:
		"logs/bowtie2/assembly_mapping/"+ wgs_name + ".log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["bowtie2"]["threads"]
	shell:
		"bowtie2 -p {threads} -x {params} -1 {input[0]} -2 {input[1]} -S {output} &> {log}"


rule samtobam:
	input:
		"results/assembly/mapping/sam/"+ wgs_name + ".sam"
	output:
		"results/assembly/mapping/bam/"+ wgs_name + ".bam"
	log:
		"logs/samtools/samtobam/" + wgs_name + ".log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["samtools"]["threads"]
	shell:
		"samtools view -@ {threads} -Sb {input} > {output} 2> {log}"


rule bam_sort:
	input:
		"results/assembly/mapping/bam/" + wgs_name + ".bam"
	output:
		"results/assembly/mapping/bam_sorted/"+ wgs_name + "_sorted.bam"
	log:
		"logs/samtools/bam_sort/"+ wgs_name + ".log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["samtools"]["threads"]
	shell:
		"samtools sort -@ {threads} {input} -o {output} &> {log}"


rule bam_index:
	input:
		"results/assembly/mapping/bam_sorted/"+ wgs_name + "_sorted.bam"
	output:
		"results/assembly/mapping/bam_sorted/"+ wgs_name + "_sorted.bam.bai"
	log:
		"logs/samtools/bam_index/" + wgs_name + ".log"
	conda:
		"../envs/env.yaml"
	shell:
		"samtools index -b {input} &> {log}"


rule polishing:
	input:
		genome="results/assembly/megahit/" + wgs_name + "/" + wgs_name + ".contigs.fa", # Input draft genome from megahit assembly
		bam="results/assembly/mapping/bam_sorted/"+ wgs_name + "_sorted.bam", # Sorted, indexed, aligned reads
		index="results/assembly/mapping/bam_sorted/"+ wgs_name + "_sorted.bam.bai" #make sure it exists before the rule is run
	params:
		output_dir = "results/assembly/pilon",
		out_prefix = wgs_name
	output:
		"results/assembly/pilon/"+ wgs_name + ".fasta" #fasta with each contig polished
	log:
		"logs/pilon/" + wgs_name + ".log"
	conda:
		"../envs/env.yaml"
	shell:
		"pilon --genome {input.genome} --bam {input.bam} --outdir {params.output_dir} --output {params.out_prefix} &> {log}"


####################################### Genome Annotation ###########################################

rule genome_annotation:
	input:
		"results/assembly/pilon/" + wgs_name + ".fasta"
	output:
		multiext("results/assembly/annotation/" + wgs_name ,
		".err", ".faa",".ffn",".fna",".fsa", ".gbk", ".gff", ".log", ".sqn", ".tbl", ".tsv", ".txt")
	params:
		prefix= wgs_name,
		out_dir = "results/assembly/annotation/",
		cpus= config["software"]["prokka"]["cpus"]
	log:
		"logs/prokka/" + wgs_name + ".log"
	conda:
		"../envs/env.yaml"
	shell:
		"""
		rm -rf {params.out_dir}
		prokka --cpus {params.cpus} --outdir {params.out_dir} --prefix {params.prefix} {input} &> {log}
		"""

### Quality Control : Compleasm & Quast

rule quast_assembly:
	input:
		fasta = "results/assembly/megahit/" + wgs_name + "/" + wgs_name + ".contigs.fa",
	output:
		multiext("results/qc/wgs/quast/" + wgs_name + "/report.","html","tex","txt","pdf","tsv"),
		multiext("results/qc/wgs/quast/" + wgs_name + "/transposed_report.","tex","txt","tsv"),
		multiext(
			"results/qc/wgs/quast/"+ wgs_name + "/basic_stats/",
			"cumulative_plot.pdf",
			"GC_content_plot.pdf",
			#			"gc.icarus.txt",
			#			"genome_GC_content_plot.pdf",
			#			"NGx_plot.pdf",
			"Nx_plot.pdf",
		),
		# only if reference is available
		#multiext(
		#	"results/qc/quast/{sample}/contigs_reports/",
		#	"all_alignments_genome.tsv",
		#	"contigs_report_genome.mis_contigs.info",
		#	"contigs_report_genome.stderr",
		#	"contigs_report_genome.stdout",
		#),
		#	"results/qc/quast/{sample}/contigs_reports/minimap_output/genome.coords_tmp",
		"results/qc/wgs/quast/" + wgs_name+ "/icarus.html",
		"results/qc/wgs/quast/" + wgs_name+ "/icarus_viewers/contig_size_viewer.html",
		"results/qc/wgs/quast/" + wgs_name+ "/quast.log"
	log:
		"logs/quast/quast.log"
	params:
		extra=config["software"]["quast"]["extra"]
	wrapper:
		"v1.31.1/bio/quast"


rule quast_polishing:
	input:
		fasta="results/assembly/pilon/" + wgs_name + ".fasta"  #fasta with each contig polished
	#		ref= config["ref"]
	output:
		multiext("results/qc/wgs/polished/" + wgs_name + "_polished/report.","html","tex","txt","pdf","tsv"),
		multiext("results/qc/wgs/polished/" + wgs_name + "_polished/transposed_report.","tex","txt","tsv"),
		multiext(
			"results/qc/wgs/polished/" + wgs_name + "_polished/basic_stats/",
			"cumulative_plot.pdf",
			"GC_content_plot.pdf",
			#			"gc.icarus.txt",
			#			"genome_GC_content_plot.pdf",
			#			"NGx_plot.pdf",
			"Nx_plot.pdf",
		),
		# only if reference is available
		#multiext(
		#	"results/qc/quast/{sample}/contigs_reports/",
		#	"all_alignments_genome.tsv",
		#	"contigs_report_genome.mis_contigs.info",
		#	"contigs_report_genome.stderr",
		#	"contigs_report_genome.stdout",
		#),
		#	"results/qc/quast/{sample}/contigs_reports/minimap_output/genome.coords_tmp",
		"results/qc/wgs/polished/" + wgs_name + "_polished/icarus.html",
		"results/qc/wgs/polished/" + wgs_name + "_polished/icarus_viewers/contig_size_viewer.html",
		"results/qc/wgs/polished/" + wgs_name + "_polished/quast.log"
	log:
		"logs/quast/polished/"+ wgs_name + "_polished.quast.log"
	params:
		extra=config["software"]["quast"]["extra"]
	wrapper:
		"v1.31.1/bio/quast"


## Compleasm to check assembled genome completeness and annotation

rule run_compleasm:
	input:
		"results/assembly/pilon/" + wgs_name + ".fasta"
	output:
		out_dir=directory("results/qc/compleasm/txome_busco/prok")
	log:
		"logs/compleasm/proteins_compleasm_prok.log"
	conda:
		"../envs/env.yaml"
	params:
		mode="busco", #modes are busco or lite, lite: Without using hmmsearch to filtering protein alignment. busco: Using hmmsearch on all candidate predicted proteins to purify the miniprot alignment to improve accuracy.
		extra="--autolineage" #uses sepp to auto determine lineage needed
	threads:
		config["software"]["compleasm"]["threads"]
	shell:
		"compleasm run -m {params.mode} {params.extra} -a {input} -o {output.out_dir} -t {threads} &> {log}"


rule wgs_run_multiqc:
	input:
		[expand("results/qc/wgs/raw/" + wgs_name + "_{number}_fastqc.html", number=numbers),
		expand("results/qc/wgs/trimmed/" + wgs_name + "{number}_fastqc.html",number= numbers),
		"results/qc/wgs/polished/" + wgs_name + "_polished/report.html",
		"results/qc/wgs/quast/" + wgs_name + "/report.html"
		 ] if config["WGS_trimmed"] else
		[expand("results/qc/wgs/raw/" + wgs_name + "_{number}_fastqc.html", number=numbers),
		 "results/qc/wgs/polished/" + wgs_name + "_polished/report.html",
		 "results/qc/wgs/quast/" + wgs_name + "/report.html"]
	output:
		"results/qc/wgs/multiqc_report.html"
	log:
		"logs/multiqc/wgs_fasta_multiqc.log"
	conda:
		"../envs/env.yaml"
	shell:
		"multiqc results/qc/wgs  -o results/qc/wgs &> {log}"
