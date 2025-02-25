#!/usr/bin/sh

export QT_QPA_PLATFORM=offscreen

DIR=/home/carlos/volume1/data/outbreaks/JRH/mrsa/data/
#for DIR in /users/bag/bmb039/my_space/projects/outbreaks_cdiff/data/jrh/*/; do
echo "Processing ${DIR}..."
	nextflow run main.nf $1 \
	-profile standard -resume \
	-entry outbreaker \
	--input $DIR \
	--output $DIR \
	--reference /home/carlos/volume1/references/mrsa_m92.fa #\
	#--whole_genome
#done
