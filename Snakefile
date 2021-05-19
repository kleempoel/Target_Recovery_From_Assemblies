#SAMPLES = ["GCA_001375635.1", "GCA_004349715.1", "GCA_000441915.1"]
with open('AssembliesID_ls.txt') as f:
	content = f.readlines()
SAMPLES = [x.strip() for x in content] 

rule all:
	input:  expand("out_fasta/{sample}.fasta", sample=SAMPLES)

rule gunzipCK:
	input: 
		"data_gz/{sample}.fasta.gz"
	output:
		"in_fasta/{sample}.fasta"
	shell:
		"gunzip -ck {input} > {output}"

rule blastn:
	input:
		sample="in_fasta/{sample}.fasta",
		targets="/mnt/shared/scratch/kleempoe/paftol/mega_353.fasta"
	output:
		"tmp/{sample}.blast"
	threads: 32
	run:
		commands = [
				"makeblastdb -in {input.sample} -dbtype nucl -parse_seqids",
				'\n echo "running blastn for {input.sample}" \n',
				'blastn  -query {input.targets}  -db {input.sample} -outfmt "6 qseqid sseqid pident length bitscore" -out {output} -max_target_seqs 5 -num_threads {threads} -gapopen 0 -gapextend 0'
			]
		for c in commands:
			shell(c)

rule get_best_hit:
	input: 
		"tmp/{sample}.blast"
	output:
		targets="tmp/{sample}.targetstxt",
		contigs="tmp/{sample}.contigstxt",
		blast_stats="stats/{sample}_blast_stats.csv"
	shell:
		"python get_best_hits.py {input} {output.targets} {output.contigs} {output.blast_stats}"

rule subset_targets:
	input: 
		targets_fasta="/mnt/shared/scratch/kleempoe/paftol/mega_353.fasta",
		targets="tmp/{sample}.targetstxt"
	output:
		"tmp/{sample}_best_targets.fasta"
	shell:
		"seqtk subseq {input.targets_fasta} {input.targets} > {output}"

rule subset_contigs:
	input: 
		sample_fasta="in_fasta/{sample}.fasta",
		contigs="tmp/{sample}.contigstxt"
	output:
		"tmp/{sample}_best_contigs.fasta"
	shell:
		"seqtk subseq {input.sample_fasta} {input.contigs} > {output}"

rule exonerate:
	input:
		targets="tmp/{sample}_best_targets.fasta",
		contigs="tmp/{sample}_best_contigs.fasta"
	output:
		"Alignments/{sample}.aln"
	shell:
		"exonerate --model est2genome --showvulgar no --percent 10 -q {input.targets} -t {input.contigs} > {output}"

rule exonerate_to_genes:
	input: 
		aln="Alignments/{sample}.aln",
		targets="/mnt/shared/scratch/kleempoe/paftol/mega_353.fasta"
	output:
		fasta="out_fasta/{sample}.fasta",
		stats="stats/{sample}_exonerate_stats.csv"
	shell:
		"python exonerate_to_genes.py {input.aln} {input.targets} {output.fasta} {output.stats}"
