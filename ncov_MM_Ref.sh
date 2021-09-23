#!/bin/bash

##For call SARS-CoV-2 nucleotide & amino acid mismatch 
##Input is multiplefasta genome without indel compared to the Reference
##Both Reference and Multi fasta have to no ambiguous base, indel
##the progarm will not compare the N position base
##Reference name have to short for bedtools
multifasta=$1
Reference=$2

#Please check the directory of nCoV-19_Mismatch_source
DIR="/media/data_storage/Hong/nCoV-19_Mismatch_Ref"

read -p "Enter output name: " output 

mkdir ./source ./Mnt ./Maa 
parallel cp *.fa* ::: Mnt Maa
workingDIR=$(pwd)

######Create Reference source
cd $workingDIR/source
cp $DIR/name_gene.txt $DIR/NC_045512.bed ../$2 .
Ref_name=${Reference%.fa*}
sed "s/NC_045512/$Ref_name/g" NC_045512.bed > $Ref_name.bed
samtools faidx $2
bedtools getfasta -fi $2 -bed $Ref_name.bed -nameOnly -fo $Ref_name.nt35.fasta
while read line; do 
	if [[ ${line:0:1} == '>' ]]; 
	then outfile=${line#>}.nt.fa; 
	echo $line > $outfile; 
	else echo $line >> $outfile;
	fi; 
done < $Ref_name.nt35.fasta

for name_gene in $(cat name_gene.txt);do
	gotranseq --sequence $name_gene.nt.fa --outseq  $name_gene.aa.tmp
	seqtk seq $name_gene.aa.tmp > $name_gene.aa.fa
	grep -v '>' $name_gene.aa.fa| grep -o -b "G\|A\|S\|T\|C\|V\|L\|I\|M\|P\|F\|Y\|W\|D\|E\|N\|Q\|H\|K\|R\|*\|X"|
	awk '{print NR,$1}'|
	column -s ":" -t|
	awk '{print $1,":", $3}'> $name_gene.aa.position;
done

grep -v '>' $2| grep -o -b "A\|T\|C\|G\|K\|M\|R\|Y\|S\|W\|B\|V\|H\|D\|N"|
awk '{print NR,$1}'|
column -s ":" -t|
awk '{print $1,":", $3}'> $Ref_name.nt.position;

######Nucleotide mismatch and Silent/Non-silent mutantion calling "callMnt_filterN_Ref.sh"
cd $workingDIR/Mnt
callMnt_filterN_Ref.sh $1 $workingDIR/source $Ref_name $output

######Amino acid mismatch calling "ncov_callMismatch_aa.sh"
cd $workingDIR/Maa
ncov_Maa_ref.sh $1 $output $workingDIR/source $Ref_name

######Making Report
cd $workingDIR
mkdir $output.MMdata
cp -r $workingDIR/Mnt/*_filterN.report ./$output.MMdata
cp $workingDIR/Maa/$output.Maa ./$output.MMdata

######Remove no need file

