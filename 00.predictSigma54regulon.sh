#!/bin/bash

# Main script for Sigma54 regulon prediction in prokaryotes
#
# Anne de Jong and Tristan Achterberg, June 2023, first release 
#               
#

PROGRAMDIR=/data/ProPr54

MODEL=$PROGRAMDIR/models/model_lstm.h5
SESSIONDIR=$1
InputType=$2  
# QueryShort: Short sequences, without gff file 
# QueryLong: full genomes, with gff file
# SESSIONDIR=/tmpdrive/propr54/propr54/86.90.198.76.tj76hjst1sp8lqfh82pmo5ascj.9004

dos2unix $SESSIONDIR/query.fna


if [ "$InputType" == "QueryLong" ]; then
	echo 'Input type: full genome<br>' >> $SESSIONDIR/sessionprogress

	# Annotate genome using FACoPv2
	echo 'Annotate sequence(s)<br>' >> $SESSIONDIR/sessionprogress
	python3 /data/FACoPv2/03.FACoPv2_annotate_genomes.py -dbdir $SESSIONDIR -webserver true

	# Compare the query.g2d.diamond.tab with the Sigma UniProt regulon IDs [database/Sigma54RegulonProteins.UniProtID]
	echo 'Predict Sigma54 Regulon genes/proteins<br>' >> $SESSIONDIR/sessionprogress
	python3 $PROGRAMDIR/CompareUniProtIDs.py -sessiondir $SESSIONDIR -query query.g2d.diamond.tab -out query.regulon.tab -regulon $PROGRAMDIR/database/Sigma54RegulonProteins.UniProtID 
fi

if [ "$InputType" == "QueryShort" ]; then
	echo 'Input type: short sequences<br>' >> $SESSIONDIR/sessionprogress
	cp $SESSIONDIR/query.fna $SESSIONDIR/query.g2d.intergenic.ffn
fi	

# predict Sigma54 binding sites
echo 'Predict Sigma54 DNA binding sites<br>' >> $SESSIONDIR/sessionprogress
python3 $PROGRAMDIR/Sigma54Predict.py -sessiondir $SESSIONDIR -modelName $MODEL 


# Map Sigma54 motifs within the k-mers and remove replicates
echo 'Map Sigma54 motifs within the k-mers and remove replicates<br>' >> $SESSIONDIR/sessionprogress
python3 $PROGRAMDIR/MapMotifs.py -sessiondir $SESSIONDIR -kmers query.Sig54_motifs.table -fna query.g2d.fna -out query.mapped_motifs.tab


if [ "$InputType" == "QueryLong" ]; then
	echo 'Merge regulon with Predicted Sigma54 binding motifs<br>' >> $SESSIONDIR/sessionprogress
	python3  $PROGRAMDIR/MergeRegulonMotifs.py -sessiondir $SESSIONDIR -annotation query.g2d.table -regulon query.regulon.tab -motifs query.mapped_motifs.tab -out query.merged_regulon_motifs.tab
fi	




# Weblogo
echo 'Make WebLogo<br>' >> $SESSIONDIR/sessionprogress
weblogo -f $SESSIONDIR/query.Sig54_motifs.fna -D fasta -o $SESSIONDIR/query.sigma54.weblogo.png -F png -A dna  -c classic -s large --errorbars NO
weblogo -f $SESSIONDIR/query.merged.sigma54_motifs.fna -D fasta -o $SESSIONDIR/query.merged.sigma54.weblogo.png -F png -A dna  -c classic -s large --errorbars NO


# Processing results
echo 'DONE<br>' >> $SESSIONDIR/sessionprogress


# tell web server when run is finished
echo 'DONE' > $SESSIONDIR/sessionend

