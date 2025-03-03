process GENERATE_CONSENSUS {
	label "genome_alignments"
	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), path(barcode_variant), path(barcode_depth), env(asmbl_depth)
	path(assembly)
	path(ref_genome)
	val(outdir)

	output:
	tuple val(barcode), path("${barcode}.fa"), emit: fa

	script:
	DEPTH_THRESHOLD=0.5
	"""
	FNAME=\$(echo $barcode | cut -d'.' -f1)
	MIN_DEPTH=\$(bc -s <<< "\$asmbl_depth * $DEPTH_THRESHOLD")
	
	awk '{if (\$3<\$MIN_DEPTH) print \$1"\t"\$2}' $barcode_depth > mask.txt
	tabix -p vcf $barcode_variant
	bcftools consensus --mask mask.txt --fasta-ref $ref_genome --output ${barcode}.fa $barcode_variant
	sed -i "s/>/>\$FNAME /" "${barcode}.fa"
	"""

	stub:
	"""
	touch ${barcode}.fa
	"""
}


process ALIGN_FROM_VCF{
	label "genome_alignments"
	publishDir "$outdir/", mode: 'copy'

	input:
	val(barcode)  // DO NOT DELETE: USED TO ENSURE PREV STEP COMPLETED
	path(ref_genome)
	path(variants)
	val(outdir)

	output:
	path("msa.fa"), emit: alignment
	path("msa_nocontrols.fa"), emit: alignment_wo_controls

	script:
	"""
	for VARIANT in ${variants}/*.gz; do
		if [[ `stat -c '%s' \$VARIANT` -gt 0 ]]; then
			bcftools index -f \$VARIANT
		fi
	done

	bcftools merge --force-samples `find ${variants}/ -type f -name "*.vcf.gz" -size +1b` > multi.vcf
	python3 $params.bin/vcf2msa.py -i multi.vcf -f $ref_genome --output-prefix prefix
	mv prefix*.fasta msa.fa
	
	# Generate alignment without reference and controls
	bcftools merge --force-samples `find ${variants}/ -type f -name "*.vcf.gz" -size +1b |\
					grep -v reference | grep -v 01 | grep -v 02 | grep -v 49 |\
					grep -v 50 | grep -v 65 | grep -v 66` > multi_wo_controls.vcf
	python3 $params.bin/vcf2msa.py -i multi_wo_controls.vcf -f $ref_genome --output-prefix no_controls
	mv no_controls*.fasta msa_nocontrols.fa.temp

	# Remove reference (vcf2msa would not work without ref_genome)
	# Each entry takes 2 lines: 1 for label, 1 for seq.
	# And ref_genome is written first, so, skip it.
	tail -n +3 msa_nocontrols.fa.temp > msa_nocontrols.fa
	"""

	stub:
	"""
	touch msa.fa msa_nocontrols.fa
	"""
}


process PREALIGN_GENOMES {
	label "genome_alignments"
	publishDir "$outdir/", mode: 'copy'

	input:
	tuple val(barcode), env(READ_DEPTH)
	val(barcode_nil)  // DO NOT DELETE: USED TO ENSURE PREV STEP COMPLETED
	path(path_var)
	val(outdir)

	output:
	tuple val(barcode), path("${barcode}.vcf.gz"), emit: preprocessed
	path("${barcode}.vcf.gz.csi"), emit: index

	script:
	MIN_DEPTH=0.7
	"""
	if [[ \$READ_DEPTH > $MIN_DEPTH || \$READ_DEPTH == $MIN_DEPTH ]]; then
		if [[ `stat -c '%s' $path_var/${barcode}.vcf.gz` -gt 0 ]]; then
			echo ${barcode}.vcf.gz > tmp.txt
			bcftools reheader --samples tmp.txt $path_var/${barcode}.vcf.gz > ${barcode}.vcf.gz
			bcftools index ${barcode}.vcf.gz
		else
			touch ${barcode}.vcf.gz
			touch ${barcode}.vcf.gz.csi
		fi
	else
		touch ${barcode}.vcf.gz
		touch ${barcode}.vcf.gz.csi
	fi
	"""

	stub:
	"""
	touch ${barcode}.vcf.gz
	touch ${barcode}.vcf.gz.csi
	"""
}


process GENERATE_PHYLOGENY_GUBBINS {
	label "phylogeny_gubbins"
	cpus=16
	publishDir "$outdir", mode: 'copy'

	input:
	path("msa.fa")
	val(outdir)

	output:
	path("recomb_free.aln"), emit: aln
	path("corrected_tree.nwk"), emit: tree

	script:
	"""
	run_gubbins.py --prefix gubbins --seed 101010 --model GTRGAMMA 	--recon-model GTRGAMMA \
		--min-snps 15 --p-val 0.01 --extensive-search --min-window-size 100 \
		--sh-test --threads $task.cpus --iterations 10 --outgroup reference msa.fa

	mask_gubbins_aln.py --aln gubbins.filtered_polymorphic_sites.fasta \
		--gff gubbins.recombination_predictions.gff \
		--out recomb_free.aln

	# mv gubbins.filtered_polymorphic_sites.fasta recomb_free.aln
	mv gubbins.final_SH_support_tree.tre corrected_tree.nwk
	"""

	stub:
	"""
	touch recomb_free.aln
	touch corrected_tree.nwk
	"""
}


process GENERATE_PHYLOGENY {
	label "phylogeny"
	cpus=16
	publishDir "$outdir", mode: 'copy'

	input:
	path("msa.fa")
	val(outdir)

	output:
	path("sh_tree.nwk"), emit: tree

	script:
	"""
	# Look for best tree (no bootstraping)
	raxmlHPC-PTHREADS-SSE3 -T $task.cpus -m GTRGAMMA \
			  -o reference -s msa.fa \
			  -n BTREE -k -p 100100

	# Run SH-like test on best tree
	raxmlHPC-PTHREADS-SSE3 -T $task.cpus -m GTRGAMMA \
			  -o reference -s msa.fa \
			  -n FINAL -p 011011 \
			  -f J -t RAxML_bestTree.BTREE
	
	mv RAxML_fastTreeSH_Support.FINAL sh_tree.nwk
	"""

	stub:
	"""
	touch sh_tree.nwk
	"""
}


process CORRECT_RECOMBINATION {
	label "phylogeny"
	publishDir "$outdir", mode: 'copy'

	input:
	path(best_tree)
	path("msa.fa")
	val(outdir)

	output:
	path("recomb_free.aln"), emit: aln
	path("corrected_tree.nwk"), emit: tree

	script:
	"""
	ClonalFrameML $best_tree msa.fa corrected_tree
	mv corrected_tree.labelled_tree.newick  corrected_tree.nwk
	mv corrected_tree.ML_sequence.fasta recomb_free.aln
	"""

	stub:
	"""
	touch recomb_free.aln
	touch corrected_tree.nwk
	"""
}


process DISTANCE_MATRIX {
	label "phylogeny_gubbins"
	publishDir "$outdir", mode: 'copy'

	input:
	path(msa)
	val(outdir)

	output:
	path("distance_matrix.tsv"), emit: snp

	script:
	"""
	snp-dists $msa > distance_matrix.tsv
	"""

	stub:
	"""
	touch distance_matrix.tsv
	"""
}
