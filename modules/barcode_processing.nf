process MERGE_BARCODES {
	publishDir "$params.output/", mode: 'copy'
	tag {"Processing $path_in"}

	input:
	path path_in

	output:
	path("*.fastq.gz")

	script:
	"""
	for barcode in \$(ls $path_in | awk 'BEGIN{FS="_"}{ print \$1 }' | sort | uniq); do
		cat $path_in/*\$barcode* > \$barcode.fastq.gz
	done
	"""

	stub:
	"""
	for barcode in \$(ls $path_in | awk 'BEGIN{FS="_"}{ print \$1 }' | sort | uniq); do
		touch \$barcode.fastq.gz
	done
	"""
}
