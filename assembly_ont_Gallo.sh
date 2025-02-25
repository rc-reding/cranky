#!/usr/bin/sh

export QT_QPA_PLATFORM=offscreen

DIR=/home/carlos/volume1/data/outbreaks/JRH/s_gallolyticus/runs_merged/
#for DIR in /users/bag/bmb039/my_space/projects/outbreaks_cdiff/data/jrh/*/; do
echo "Processing ${DIR}..."
	nextflow run main.nf $1 \
	-profile standard -resume \
	-entry assemble_reads \
	--input $DIR \
	--output $DIR \
	--reference /home/carlos/volume1/references/s_gallolyticus_FDAARGOS.fa
#done
