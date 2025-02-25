#!/usr/bin/sh

# Command to run. $1 allows the addition of -stub-run (needs to be there)
DIR="/home/carlos/volume1/data/outbreaks/JRH/mrsa/data"
#for DIR in /home/carlos/volume1/data/outbreaks/C_diff/data/jrh/run*; do
			echo "Processing $DIR"
			nextflow run main.nf $1 \
		  -entry preprocess_reads \
		  --input $DIR/reads/ \
		  --output $DIR/reads_merged
#done
