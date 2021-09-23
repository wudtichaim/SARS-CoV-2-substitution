#!/bin/bash

#Split multiFasta file into separated single fasta file

multifasta=$1 

while read line; do 
	if [[ ${line:0:1} == '>Wuhan_Hu1_NC_045512' ]];
		then outfile=${line#>}.fa;
		echo $line > $outfile;
	else echo $line >> $outfile;
	fi;
done < $1



