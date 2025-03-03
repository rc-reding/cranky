process SPECIES_PROFILING_ONT {
	label "kraken2"
	tag { "Detecting species in sample $sample" }

	publishDir "$outdir", overwrite: true, mode: 'copy'
	cpus = 4

	input:
	tuple val(sample), path(reads)
	val(outdir)
	val(results)

	output:
	tuple val(sample), path("${sample}.kraken.report"), emit: report
	tuple val(sample), path("results/${sample}.kraken.results"), emit: results

	script:
	"""
	kraken2 -db $params.db --report ${sample}.kraken.report \
		--threads $task.cpus --memory-mapping \
		--output ${sample}.kraken.results \
		$reads

	mkdir results
	mv ${sample}.kraken.results results/
	"""

	stub:
	"""
	touch $results/${sample}.kraken.results ${sample}.kraken.report
	"""
}
