#4. Gblast extracted/edited summary - 6.057s

#Set customized heading varaibles
exons=`echo $2`
#exons=`grep ">" $1 -c`
length=$(for i in `seq $exons`; do echo "E"$i"_L\t"; done | paste -s -d "" | sed 's/\\t$//g')
percent_similarity=$(for i in `seq $exons`; do echo "E"$i"_%sim\t"; done | paste -s -d "" | sed 's/\\t$//g')
percent_gaps=$(for i in `seq $exons`; do echo "E"$i"_%gap\t"; done | paste -s -d "" | sed 's/\\t$//g')
splice_sites=$(for i in `seq $(expr $exons - 1 | awk '{print$1*2}')`; do echo "SS_"$i"\t";done | paste -s -d "" | sed -e 's/^/SBG\\t/g' -e 's/$/SAG/g')
insertions=$(for i in `seq $exons`; do echo "E"$i"_ins\t"; done | paste -s -d "" | sed 's/\\t$//g')
deletions=$(for i in `seq $exons`; do echo "E"$i"_del\t"; done | paste -s -d "" | sed 's/\\t$//g')

#Gedit summary
while read i
do
#Species 0f interest
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
#Taxonomy order of the species
gr=`echo $j | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} ~/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $i | awk '{print$2}')/aswin/$1/
#Query used
qry1=`grep $j $4/gblast_query_subject_set | awk '{print$1}'`
qry=`ls *.fa | grep $qry1".fa$" | grep -v "flanking"`
#Exon length
k=$(for x in `grep ">" $qry | tr -d ">"`; do y=`seqtk subseq gblast_edited_consensus.fa <(echo $x) | grep -v ">" | wc | awk '{print$3-$1}'`; if [[ $y == "0" ]]; then y="-"; else :; fi; echo $y; unset y; done | paste -s -d " ")
#Similarity
sim=$(for x in `grep ">" $qry | tr -d ">"`; do y=`awk -F "(" '/Similarity/ {print$NF}' pairwise_"$x".aln | tr -d ")"`; if [[ $y == "" ]]; then y="-"; else :; fi; echo $y; unset y; done | paste -s -d " ")
#Gaps
gap=$(for x in `grep ">" $qry | tr -d ">"`; do y=`awk -F "(" '/Gaps/ {print$NF}' pairwise_"$x".aln | tr -d ")"`; if [[ $y == "" ]]; then y="-"; else :; fi; echo $y; unset y; done | paste -s -d " ")
#Insertions
l=$(for x in `grep ">" $qry | tr -d ">"`; do y=`grep -w "$x" gblast_edited_insertions | awk '{print$2}'`; if [[ $y == "" ]]; then y="-"; else :; fi; echo $y; unset y; done | paste -s -d " ")
#Deletions
m=$(for x in `grep ">" $qry | tr -d ">"`; do y=`grep -w "$x" gblast_edited_deletions | awk '{print$2}'`; if [[ $y == "" ]]; then y="-"; else :; fi; echo $y; unset y; done | paste -s -d " ")
#Check ORF presence absence
n=`grep -v ">" gblast_edited_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $n == "" ]] ; then o="❌";else o="✅"; fi
#Translated protein length
p=`grep -v ">" gblast_edited_consensus_orf.fa | wc | awk '{print$3-$1}'`
#Consensus nucleotide length
q=`grep -v ">" gblast_edited_consensus.fa | wc | awk '{print$3-$1}'`
#splice sites
for z in `grep ">" $qry | tr -d ">"`
do
z1=`seqtk subseq gblast_edited_consensus.fa <(echo $z) | grep -v ">"`
if [[ $z1 == "" ]]; then
ss1="-"
ss2="-"
else
ss1=`sed -n "s/$z1.*//p" extracted_flanking_region.fa | rev | cut -c -2 | tr [:lower:] [:upper:] | rev`
if [[ $ss1 == "" ]]; then ss1="-"; else :; fi
ss2=`sed -n "s/.*$z1//p" extracted_flanking_region.fa | cut -c -2 | tr [:lower:] [:upper:]`
if [[ $ss2 == "" ]]; then ss2="-"; else :; fi
fi
echo $ss1 $ss2
unset z1 ss1 ss2
done > gedit_splice_sites
r=`cat gedit_splice_sites | paste -s -d " " | tr [:lower:] [:upper:]`
if [[ $r == "" ]]; then r=`expr $exons - 1 | awk '{print$1*2}' | xargs -I {} sh -c 'yes " -" | head -{} | xargs'`; else :; fi
#Gaps/anonymous bases in the subject
gp=`cat pairwise_exon_* | awk '/^extracted/{print$3}' | tr -cd "N" | wc | awk '{print$3-$1}'`
#Total query % covered = total query covered / total query length * 100
tl=`cat $qry | grep -v ">" | wc | awk '{print$3-$1}'`
qc=`awk '{print$3}' gblast_edited_query_covered | awk '{sum+=$1;} END{print sum;}' | awk -v l="$tl" '{print($1/l)*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
echo -e "$j\t$gr\t$qry\t$k\t$sim\t$gap\t$r\t$l\t$m\t$q\t$p\t$gp\t$qc\t$o"
cd $4
unset j gr qry1 qry k sim gap r l m n o p q gp tl qc
done < <(grep -if <(grep -if <(cut -f1,2 -d "_" gedit_subjects |less) $3 | awk '{print$2}') /home/neo/bird_db1/aswin/database_details/all_genome_paths) | sed "1i Species\tGroup\tQuery\t$length\t$percent_similarity\t$percent_gaps\t$splice_sites\t$insertions\t$deletions\tconsensus\tprotein\tGaps\tQuery_cov\tORF\n" | column -t

#Unset heading variables
unset length percent_similarity percent_gaps splice_sites insertions deletions
