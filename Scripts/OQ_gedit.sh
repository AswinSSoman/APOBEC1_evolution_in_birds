#GEDIT

while read z
do
y=`echo $z | awk '{print$1}' | cut -f1,2 -d "_"`
echo " ▶ " $y
cd $(echo $z | awk '{print$2}')/aswin/$1/

#Query used
q1=`grep $y $2/gblast_query_subject_set | awk '{print$1}'`
q=`ls *.fa | grep $q1 | grep -v "flanking"`

#If the query has only one exon then don't need to look for splice sites, just pairwise align & extract the aligned region
note=`grep ">" $q -c`
if [[ $note == 1 ]]
then
start=`cat pairwise_exon_1.aln | grep -B2 "^extrac" | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":" | awk '{print$0+1}'`
#Extract the region from which the pairwise alignment start i.e. remove the regions before alignment starts (usually alignment starts with START codon ATG)
cat pairwise_exon_1.aln | grep -B2 "^extrac" | grep "|" -C1 | cut -c $start- > tmp1
end=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+1}'`
ln=`seqtk subseq $q <(echo exon_1) | grep -v ">" | awk NF | wc | awk '{print$3-$1}'`
#Extract some stats from pairwise alignment
#For each exon : Exon name, Total query length, Total query covered, Percent query covered
cat tmp1 | cut -c -$end | head -1 | tr -d "-" | wc | awk '{print$3-$1}' | awk -v p="$ln" '{print p,$1,($1/p)*100}' | sed "s/^/exon_1 /g" > gblast_edited_query_covered
#insertions
cat tmp1 | cut -c -$end | head -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exon_1 /g" > gblast_edited_insertions
#deletions
cat tmp1 | cut -c -$end | tail -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exon_1 /g" > gblast_edited_deletions
cat <(echo ">exon_1") <(cat tmp1 | cut -c -$end | tail -1) | tr -d "-" > gblast_edited_consensus.fa
rm tmp1
unset start end ln
else

#If query has more than 1 exon & exon 1 gave blast hit
#Exclude 1st and last exon from the loop to find 2 splice sites since they have only one splice site
#For the first exon : extract from start codon
if [[ -e pairwise_exon_1.aln ]]
then
start=`cat pairwise_exon_1.aln | grep -B2 "^extrac" | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":" | awk '{print$0+1}'`
#Extract the region from which the pairwise alignment start i.e. remove the regions before alignment starts (usually alignment starts with START codon ATG)
cat pairwise_exon_1.aln | grep -B2 "^extrac" | grep "|" -C1 | cut -c $start- > tmp1
#Look for nearest splice site "GT/GC" within next 5 flanking nucleotides after the query exon hit ended (the threshold flanking bases & other variable splice sites are customizable)
ssc1=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+2}'`
ssc2=`cat tmp1 | grep "|" -C1 | cut -c $ssc1- | tail -1 | cut -c -5 | grep -i "GT" -aob | head -1 | cut -f1 -d ":"`
#Allow splice site variants
if [[ $ssc2 == "" ]]; then ssc2=`cat tmp1 | grep "|" -C1 | cut -c $ssc1- | tail -1 | cut -c -5 | grep -i "GC" -aob | head -1 | cut -f1 -d ":"`; else :; fi
ssc3=`expr $ssc1 + $ssc2 | awk '{print$0-1}'`
if [[ $ssc3 == "" ]]; then er=`echo "     ↳ exon_1 - No 5 prime splice site found! "`; ssc3=`echo $ssc1 | awk '{print$0-1}'`; else :; fi
ln=`seqtk subseq $q <(echo "exon_1") | grep -v ">" | awk NF | wc | awk '{print$3-$1}'`
#Extract some stats from pairwise alignment
#For each exon : Exon name, Total query length, Total query covered, Percent query covered
cat tmp1 | cut -c -$ssc3 | head -1 | tr -d "-" | wc | awk '{print$3-$1}' | awk -v p="$ln" '{print p,$1,($1/p)*100}' | sed "s/^/exon_1 /g" > gblast_edited_query_covered
#insertions
cat tmp1 | cut -c -$ssc3 | head -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exon_1 /g" > gblast_edited_insertions
#deletions
cat tmp1 | cut -c -$ssc3 | tail -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exon_1 /g" > gblast_edited_deletions
cat <(echo ">exon_1") <(cat tmp1 | cut -c -$ssc3 | tail -1) | tr -d "-" > gblast_edited_consensus.fa
if [ -n "$er" ] ; then echo "$er"; else :; fi
rm tmp1
unset start ssc1 ssc2 ssc3 ln l p er
else
echo "     ↳ Exon 1 alignment not found"

#Create an empty file for making sblast edited output files because the first time they are created it should be direct redirection(>) not appending (>>)
echo -n > gblast_edited_query_covered
echo -n > gblast_edited_insertions
echo -n > gblast_edited_deletions
echo -n > gblast_edited_consensus.fa
fi

last_exon=`grep ">" $q | tail -1 | tr -d ">"`

#Run if total number of exons are greater than 2
if [ $note -gt 2 ]
then
#Middle exons are flanked by splice sites on both sides unlike 1st & last exons
#Extract complete exons from subject based on the pairwise alignment between query sequence & subject genome with 15bp flanking sequence
for i in `ls pairwise_exon_* | sort -V | egrep -v "exon_1.aln|$last_exon.aln"`
do
#Exon number
j=`echo $i | cut -f2,3 -d '_' | sed 's/.aln//g'`
echo ">"$j
#Splice site 1 for the respective exon
ssc4=`cat $i | grep -B2 "^extrac" | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":"`
#Look for nearest splice site "AG" within next 5 flanking nucleotides after the query exon hit ended (the threshold flanking bases & other variable splice sites are customizable)
ssc5=`cat $i | grep -B2 "^extrac" | grep "|" -C1 | cut -c -$ssc4 | tail -1 | rev | cut -c -5 | grep -i "GA" -aob | head -1 | cut -f1 -d ":"`
ssc6=`expr $ssc4 - $ssc5 | awk '{print$0+1}'`
if [[ $ssc6 == "" ]]; then echo "     ↳ " $j " - No 3 prime splice site found! " >> er1; ssc6=`echo $ssc4 | awk '{print$0+1}'`; else :; fi
cat $i | grep -B2 "^extrac" | grep "|" -C1 | cut -c $ssc6- > tmp1
ssc7=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+2}'`
#Look for nearest splice site "GT/GC" within next 5 flanking nucleotides after the query exon hit ended (the threshold flanking bases & other variable splice sites are customizable)
ssc8=`cat tmp1 | grep "|" -C1 | cut -c $ssc7- | tail -1 | cut -c -5 | grep -i "GT" -aob | head -1 | cut -f1 -d ":"`
#Allow splice site variants
if [[ $ssc8 == "" ]]; then ssc8=`cat tmp1 | grep "|" -C1 | cut -c $ssc7- | tail -1 | cut -c -5 | grep -i "GC" -aob | head -1 | cut -f1 -d ":"`; else :; fi
ssc9=`expr $ssc7 + $ssc8 | awk '{print$0-1}'`
if [[ $ssc9 == "" ]]; then echo "     ↳ " $j " - No 5 prime splice site found! " >> er2; ssc9=`echo $ssc7 | awk '{print$0-1}'`; else :; fi
ln=`seqtk subseq $q <(echo $j) | grep -v ">" | awk NF | wc | awk '{print$3-$1}'`
#Extract some stats from pairwise alignment
#For each exon : Exon name, Total query length, Total query covered, Percent query covered
cat tmp1 | cut -c -$ssc9 | head -1 | tr -d "-" | wc | awk '{print$3-$1}' | awk -v p="$ln" '{print p,$1,($1/p)*100}' | sed "s/^/$j /g" >> gblast_edited_query_covered
#insertions
cat tmp1 | cut -c -$ssc9 | head -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/$j /g" >> gblast_edited_insertions
#deletions
cat tmp1 | cut -c -$ssc9 | tail -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/$j /g" >> gblast_edited_deletions
cat tmp1 | cut -c -$ssc9 | tail -1
rm tmp1
unset j ssc4 ssc5 ssc6 ssc7 ssc8 ssc9 ln l p
done | tr -d "-" >> gblast_edited_consensus.fa
if [[ -e er1 ]] ; then cat er1; rm er1; else :; fi
if [[ -e er2 ]] ; then cat er2; rm er2; else :; fi
#If the number of exons are just 2 then don't run this,rather run the next code for last exon
else :
fi

#Extract last exon
if [[ -e $(ls pairwise_$last_exon.aln) ]] 2>/dev/null
then
ssc10=`ls pairwise_exon_* | grep $last_exon | xargs cat | grep -B2 "^extrac" | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":"`
ssc11=`ls pairwise_exon_* | grep $last_exon | xargs cat | grep -B2 "^extrac" | grep "|" -C1 | cut -c -$ssc10 | tail -1 | rev | cut -c -5 | grep -i "GA" -aob | head -1 | cut -f1 -d ":"`
ssc12=`expr $ssc10 - $ssc11 | awk '{print$0+1}'`
if [[ $ssc12 == "" ]]; then er=`echo "     ↳ " $last_exon " - No 3 prime splice site found! "`; ssc12=`echo $ssc10 | awk '{print$0+1}'`; else :; fi
ls pairwise_exon_* | grep $last_exon | xargs cat | grep -B2 "^extrac" | grep "|" -C1 | cut -c $ssc12- > tmp1
end=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+1}'`
ln=`seqtk subseq $q <(echo $last_exon) | grep -v ">" | awk NF | wc | awk '{print$3-$1}'`
cat tmp1 | cut -c -$end | head -1 | tr -d "-" | wc | awk '{print$3-$1}' | awk -v p="$ln" '{print p,$1,($1/p)*100}' | sed "s/^/$last_exon /g" >> gblast_edited_query_covered
#insertions
cat tmp1 | cut -c -$end | head -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/$last_exon /g" >> gblast_edited_insertions
#deletions
cat tmp1 | cut -c -$end | tail -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/$last_exon /g" >> gblast_edited_deletions
cat <(echo ">"$last_exon) <(cat tmp1 | cut -c -$end | tail -1) | tr -d "-" >> gblast_edited_consensus.fa
rm tmp1
if [ -n "$er" ] ; then echo "$er"; else :; fi
unset ssc10 ssc11 ssc12 end ln l p last_exon er
else
echo "     ↳ " $last_exon" alignment not found"
fi

#End of gedit
fi
#Generate ORF
transeq <(grep -v ">" gblast_edited_consensus.fa | paste -s -d "" | sed "1i \>$j") gblast_edited_consensus_orf.fa 2>/dev/null
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" gblast_edited_consensus.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" gblast_edited_consensus_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > gblast_edited_translated
#For better visualization of ORF
gnaa=`head -1 gblast_edited_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" gblast_edited_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 gblast_edited_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > gblast_edited_nucleotide_amino_acid_view
#cat <(cut -c -251 gblast_edited_translated) <(cut -c 253-504 gblast_edited_translated) <(cut -c 505- gblast_edited_translated) | GREP_COLORS="mt=41" egrep "^ M |$| \* |$| X |$" --color=always > gblast_edited_nucleotide_amino_acid_view
unset y q1 q note gnaa
cd $2
done < <(grep -if gedit_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths)
