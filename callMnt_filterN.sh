#!/bin/bash

#For nCoV-19 Nucleotide Mismatch Calling
multifasta=$1
ref=$2

#Split multiFasta file into separated single fasta file
Spilt-multifasta.sh $multifasta

#call nt Ref position
grep -v '>' $2| grep -o -b "A\|T\|C\|G\|K\|M\|R\|Y\|S\|W\|B\|V\|H\|D\|N"|
awk '{print NR,$1}'|
column -s ":" -t|
awk '{print $1,":", $3}'> Ref.position

#Comparison
for sample in $(cat sample_names.tmp);do
	grep -v '>' $sample/$sample.fa| grep -o -b "A\|T\|C\|G\|K\|M\|R\|Y\|S\|W\|B\|V\|H\|D\|N"|
	awk '{print NR,$1}'|
	column -s ":" -t|
	awk '{print $1,":", $3}'> $sample/${sample%.fasta*}.position;
	diff Ref.position $sample/${sample%.fasta*}.position > $sample/${sample%.fasta*}.diff
	grep -i "<" $sample/${sample%.fasta*}.diff > $sample/${sample%.fasta*}.ref.tmb
	grep -i ">" $sample/${sample%.fasta*}.diff| awk '{print $4}' > $sample/${sample%.fasta*}.sample.tmb
	paste $sample/${sample%.fasta*}.ref.tmb $sample/${sample%.fasta*}.sample.tmb|sed 's/\t/>/g'|cut -d ' ' -f 2,3,4 > $sample/${sample%.fasta*}.nt_mismatch.tmb
	echo ${sample%.fasta*}|cat - $sample/${sample%.fasta*}.nt_mismatch.tmb > $sample/${sample%.fasta*}.nt_mismatch.tmb2
	grep -v "N>\|>N" $sample/${sample%.fasta*}.nt_mismatch.tmb2|sed  's/ : /:/g' > $sample/${sample%.fasta*}.nt_mismatch
	cp $sample/${sample%.fasta*}.nt_mismatch .
	rm -f $sample/*.tmb
	
done
paste *.nt_mismatch > $1.$2.nt
rm -f *.nt_mismatch

#count varaint position
grep -i ':' $1.$2.nt|
awk -- '{for (i = 1; i <= NF; i++) wc[$i] += 1}; END {for (w in wc) print w, wc[w]}'| sort -n > countVariant.txt

#Find sample inclue variant
mkdir FindSample
for sample in $(cat sample_names.tmp);do
	cp $sample/${sample%.fasta*}.nt_mismatch FindSample
done
cp countVariant.txt FindSample
cd FindSample
awk '{print $1}' countVariant.txt > countVaraint.tmp
for varaint in $(cat countVaraint.tmp);do
	grep -Rw ./*.nt_mismatch -e $varaint >> FindSample.txt
done
cd ..
cp FindSample/FindSample.txt .
mkdir $1.callMnt_filterN.report
mv FindSample.txt countVariant.txt $1.$2.nt $1.callMnt_filterN.report


#*****The multifasta have to no indel*****
