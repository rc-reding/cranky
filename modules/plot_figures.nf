process PLOT_PHYLOGENY {
	label "phylogeny_gubbins"
	publishDir "$outdir", mode: 'copy'

	input:
	val(phylogeny_dir)
	path(snp_dists)  // Only to ensure the process runs after phylogeny
	val(mlst_dir)
	val(outdir)

	output:
	path("phylogeny.pdf"), emit: phylogeny
	path("cladogram.pdf"), emit: cladogram
	path("snp_histogram.pdf"), emit: snp_hist
	path("distance_matrix.pdf"), emit: dist_matrix

	script:
	"""
	python3 $params.bin/plot_dendrogram.py $phylogeny_dir/ ./ $mlst_dir/
	python3 $params.bin/plot_snp_histogram.py $phylogeny_dir/no_controls/distance_matrix.tsv \
										   $phylogeny_dir
	"""

	stub:
	"""
	touch phylogeny.pdf cladogram.pdf snp_histogram.pdf distance_matrix.pdf
	"""
}

process PLOT_SNP_DISTANCES {
	label "phylogeny_gubbins"
	publishDir "$outdir", mode: 'copy'

	input:
	path(snp_dists) 
	val(outdir)

	output:
	path("snp_histogram.pdf"), emit: snp_hist
	path("distance_matrix.pdf"), emit: dist_matrix

	script:
	"""
	python3 $params.bin/plot_snp_histogram.py $snp_dists $outdir
	"""

	stub:
	"""
	touch snp_histogram.pdf distance_matrix.pdf
	"""
}

process PLOT_COVERAGE {
	label "phylogeny"
	publishDir "$outdir", mode: 'copy'

	input:
	val(depth_dir)
	val(coverage)  // Only to ensure process runs after depth
	val(outdir)

	output:
	path("samples_breadth_coverage.pdf"), emit: breadth
	path("samples_coverage_depth.pdf"), emit: depth

	script:
	"""
	python3 $params.bin/plot_coverage.py $depth_dir/ ./
	"""

	stub:
	"""
	touch samples_coverage_depth.pdf samples_breadth_coverage.pdf
	"""
}

process PLOT_QC {
	label "phylogeny"
	publishDir "$outdir", mode: 'copy'

	input:
	val(qc_dir)
	val(coverage)  // Only to ensure process runs after depth
	val(outdir)

	output:
	path("samples_qc.pdf"), emit: breadth
	path("samples_qc_filtered.pdf"), emit: depth

	script:
	"""
	python3 $params.bin/plot_qc.py $qc_dir/ ./
	"""

	stub:
	"""
	touch samples_qc.pdf samples_qc_filtered.pdf
	"""
}
