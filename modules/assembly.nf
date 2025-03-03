process ASSEMBLE_ONT {
	label "dragonflye"
	cpus 8

	publishDir "$outdir/", mode: 'copy'

	input:
	tuple val(barcode), path(reads_merged)
	val(outdir)

	output:
	tuple val(barcode), path("${barcode}.assembly.fa"), emit: assembly
	tuple val(barcode), path("${barcode}.assembly.gfa"), emit: graph

	/*
	 Filter barcodes with a file size of ~500kb or less. This avoids
	 a 'divided by 0 error' and allows the workflow to continue.
	*/

	script:
	MDL_NAME = 'r941_e81_sup_g514'
	"""
	genomeSize=\$(wc -c ${reads_merged} | cut -d " " -f 1)
	minimum_genomeSize=5120000
	if [ \$genomeSize -ge \$minimum_genomeSize ]; then
		# Usage of --model XXX leads to _very_ slow assembly times (10h<)
		# Usage of --medaka 1 leads to _very_ slow assembly times 
		# Ignoring --model and --medaka 1 does not change phylogeny
		dragonflye --reads ${reads_merged} --outdir flye_out \
			--nanohq --cpus $task.cpus --seed 101010 \
			-gsize 5M # --model $MDL_NAME # --medaka 1

		mv flye_out/contigs.fa ${barcode}.assembly.fa
		mv flye_out/*.gfa ${barcode}.assembly.gfa
	else
		touch ${barcode}.assembly.fa
		touch ${barcode}.assembly.gfa
	fi
	"""

	stub:
	"""
	touch ${barcode}.assembly.fa
	touch ${barcode}.assembly.gfa
	"""
}

process MAP_REFERENCE {
	label "reference_mapping"
	cpus=8

	publishDir "$outdir/", mode: 'copy'

	input:
	tuple val(barcode), path(reads_merged)
	path(ref_genome)
	val(outdir)

	output:
	tuple val(barcode), path("${barcode}.mapped.bam"), emit: bam

	script:
	"""
	# Index reference
	samtools faidx $ref_genome

	# Map reads to reference
	minimap2 -x lr:hq --seed 101010 -a --secondary=no -t $task.cpus -o aln.sam $ref_genome $reads_merged;
	samtools view --threads $task.cpus -F 4 -q 50 -b -o aln_filtered.bam aln.sam;
	samtools sort --threads $task.cpus -o ${barcode}.mapped.bam aln_filtered.bam;
	"""

	stub:
	"""
	touch ${barcode}.mapped.bam
	"""

}


process VARIANT_CALL_CLAIR3 {
	label "variant_call"
	cpus=6

	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), path(asmbl), env(asmbl_depth)
	path(ref_genome)
	val(outdir)
	
	output:
	tuple val(barcode), path("${barcode}.vcf.gz"), emit: vcf

	script:
	//MDL_NAME='r1041_e82_400bps_sup_v500'
	MDL_NAME='r941_prom_sup_g5014'
	DEPTH_THRESHOLD=0.5
	"""
	samtools index $asmbl
	samtools faidx $ref_genome
	
	MIN_DEPTH=\$(bc -sl <<< "\$asmbl_depth * ${DEPTH_THRESHOLD}")

	# Run clair3	
	run_clair3.sh \
		--bam_fn=$asmbl \
		--ref_fn=$ref_genome \
		--threads=${task.cpus} \
		--platform="ont" \
		--include_all_ctgs \
		--min_coverage=\$MIN_DEPTH \
		--min_mq=15 \
		--qual=25 \
		--snp_min_af=0.85 \
		--call_snp_only \
		--no_phasing_for_fa \
		--output=clair_out \
		--model_path=$params.models/clair3/${MDL_NAME}

	mv clair_out/merge_output.vcf.gz ${barcode}.vcf.gz
	"""

	stub:
	"""
	touch ${barcode}.vcf.gz
	"""
}

process VARIANT_CALL_SAMTOOLS {
	label "variant_call"
	cpus=6

	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), path(asmbl)
	path(ref_genome)
	val(outdir)
	
	output:
	tuple val(barcode), path("${barcode}.vcf.gz"), emit: vcf

	script:
	"""
	samtools index $asmbl
	samtools faidx $ref_genome

	# Run bcftools
	bcftools mpileup -X ont \
		--fasta-ref $ref_genome \
		--seed 12354 \
		--min-MQ 15 --min-BQ 15 \
		--skip-indels \
		--threads $task.cpus \
		--output-type z \
		--output ${barcode}.vcf.gz \
		$asmbl
	"""

	stub:
	"""
	touch ${barcode}.vcf.gz
	"""
}


process VARIANT_CALL_LONGSHOT {
	label "variant_call_LS"
	cpus=6

	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), path(asmbl)
	path(ref_genome)
	val(outdir)
	
	output:
	tuple val(barcode), path("${barcode}.vcf.gz"), emit: vcf

	script:
	"""
	samtools index $asmbl
	samtools faidx $ref_genome

	# Run longshot
	longshot --bam $asmbl \
			 --ref $ref_genome \
			 --min_cov 25 \
			 --min_mapq 15 \
			 --min_allele_qual 15 \
			 --hap_assignment_qual 25 \
			 --min_alt_frac=0.85 \
			 --max_cigar_indel 100 \
			 --sample_id ${barcode} \
			 --out ${barcode}.vcf

	bgzip ${barcode}.vcf
	"""

	stub:
	"""
	touch ${barcode}.vcf.gz
	"""
}
