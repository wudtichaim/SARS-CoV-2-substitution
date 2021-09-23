#!/bin/bash

#For nCoV-19 fill gap with N

multifasta=$1

#Please check the directory of nCoV-19_Mismatch_source
DIR="/media/data_storage/Hong/nCoV-19_Mismatch_Ref"

#Mafft
cat $DIR/NC_045512.fa $1 > $1.ref
mafft --thread -1 $1.ref > $1.ref.aln

#Convert to N
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' $1.ref.aln > $1.ref.aln.upper
sed "s/-/N/g" $1.ref.aln.upper > $1.ref.aln.upperN.fa

grep '>' $1.ref.aln | sed 's/>//g' > sample_names.tmp

#Extract multifasta
while read line; do 
	if [[ ${line:0:1} == '>' ]];
		then outfile=${line#>}.woGAP.fa;
		echo $line > $outfile;
	else echo $line >> $outfile;
	fi;
done < $1.ref.aln.upperN.fa

for sample in $(cat sample_names.tmp);do
	grep -v '>' $sample.woGAP.fa|
	grep -o -b "A\|T\|C\|G\|N\|-"|
	awk '{print NR,$1}'|
	column -s ":" -t|
	awk '{print $1, ":", $3}' > $sample.position
	rm -f $sample.woGAP.fa	
done

#find insertion position
grep 'N' NC_045512.position > insertion_position
#number_line=$(cat insertion_position) 
for sample in $(cat sample_names.tmp);do
	grep -v -xF -f insertion_position $sample.position > $sample.reposition;	
	rm -f $sample.position
done

#Gernerate woINDEL multifasta
for sample in $(cat sample_names.tmp); do
#add name
	echo ">"$sample > $sample.tr
	awk '{print $3}' $sample.reposition | tr -d '\n'>> $sample.tr
	rm -f $sample.reposition	
done
rm -f NC_045512.tr
awk 1 *.tr > $1.woINDEL.fa
##rm -f *.tr
#rm *.ref *.aln *.upper *.upperN.fa
	

