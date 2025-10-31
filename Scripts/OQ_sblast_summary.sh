#Set customized heading varaibles
exons=`echo $2`
sblast_exons=$(for i in `seq $exons`; do echo "sE"$i"_miss\t"; done | paste -s -d "" | sed 's/\\t$//g')

while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
#Taxonomy order of the species
gr=`echo $j | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} ~/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $i | awk '{print$2}')/aswin/$1/
#SRA size
a=`grep $j ~/bird_db1/aswin/database_details/all_sra_size | awk '{print$NF}'`
#Query used
qry=`ls sblast_query.fa`
#Query length
b=`grep -v ">" sblast_query.fa | wc | awk '{print$3-$1}'`
#gedit consensus length
c=`grep -v ">" gblast_edited_consensus.fa | wc | awk '{print$3-$1}'`
#sblast consensus length
d=`cat sblast_edited_consensus.fa | grep -v ">" | wc | awk '{print$3-$1}'`
#Exons with missing bases
e=`awk 'NR>1{print$2-$5}' optimal_query_sblast_summary | paste -s -d " "`
#Exons with No hits
f=`grep 'No hits found' optimal_query_sblast.outfmt3 -B5 | awk '/Query=/ {print$2}' | paste -s -d "," | sed 's/,exon_/,/g'`
if [[ $f == "" ]]; then f="-"; else :;fi
#Number of exons with almost complete (>90% covered) query hit
ewah=`awk '$NF>90{gsub(/exon_/, "");print$1}' optimal_query_sblast_summary | paste -s -d ","`
if [[ $ewah == "" ]]; then ewah="-"; else :;fi
#START of CDS
g=`grep -v ">" sblast_edited_consensus.fa | paste -s -d "" | cut -c -3 | tr '[:lower:]' '[:upper:]'`
if [[ $g == "" ]]; then g="-"; else :; fi
#Depth_range - the number of bases in the SRA blast hits to support the presence of a nuclotide (eg: some bases had 2 reads support, while some had 100 reads support)
h=`cat optimal_query_sblast.outfmt3_processed_data/*_without_query.txt_count.txt | tr -d "{}':,.[A-Z,a-z]-" | awk '{max=0;for(i=1;i<=NF;i++)if($i!~/NA/&&$i>max){max=$i;}count+=max;print max}' | sort -n | sed -n '1p;$p' | paste -s -d " " | awk '{print$1"X-"$2"X"}'`
#Poor_depths_in_SRA_hits[Depth(frequency)] - Depth means: number of reads to support a base at a particular location; frequencey means: Number of nucleotides with read support less than 3 reads
k=`cat optimal_query_sblast.outfmt3_processed_data/*_without_query.txt_count.txt | tr -d "{}':,.[A-Z,a-z]-" | awk '{max=0;for(i=1;i<=NF;i++)if($i!~/NA/&&$i>max){max=$i;}count+=max;print max}' | awk '$1<3' | wc -l`
#Variant residuals: number of positions with variant residuals below 2; Note-threshold difference can be customized
l=`cat optimal_query_sblast.outfmt3_processed_data/variants | awk '$2<2' | wc -l`
#frameshifting variants (the threshold variant residual can be customized)
m=`cat optimal_query_sblast.outfmt3_processed_data/frame_shifting_variants | awk '$2<4'| awk -F '\t' '$NF~/-/' | wc -l`
#Gaps (Ns) bases
n=`cat sblast_edited_consensus.fa  | grep -v ">" | grep -i "n" -aob | wc -l`
#Other anonymous bases other than A,T,G,C & N
oab=`awk '!/>/ {gsub(/[ATGCNatgcnx]/,"")}1' sblast_edited_consensus.fa | sed '/>exon/! s/[a-z]/\U&/g' | awk '/[A-Z]/' RS=">" | awk NF | paste -d " " - - | sed 's/exon_//g' | awk '{print length($2)"("$1")"}' | paste -s -d ","`
if [[ $oab == "" ]]; then oab="-"; else :; fi
#Total insertions
if [[ -s sblast_edited_insertions ]]; then
ins=`awk '{a+=$2} END{print a}' sblast_edited_insertions`
else ins="-"; fi
#Total deletions
if [[ -s sblast_edited_deletions ]]; then
del=`awk '{a+=$2} END{print a}' sblast_edited_deletions`
else del="-"; fi
#Total % query covered in sblast
qc1=`awk 'NR>1{print$2,$5}' optimal_query_sblast_summary | awk '{a+=$1; b+=$2} END{print b/a*100}' | awk '{$1+=0}1' CONVFMT="%.1f"`
if [[ $qc1 == "" ]]; then qc1="-"; else :; fi
#Total % query covered in sedit
tl=`cat $qry | grep -v ">" | wc | awk '{print$3-$1}'`
if [[ -f sblast_edited_query_covered ]]; then
qc2=`awk '{print$3}' sblast_edited_query_covered | awk '{sum+=$1;} END{print sum;}' | awk -v l="$tl" '{print($1/l)*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
else qc2="0"; fi
# % of the gene in which the first Stop codons appears
m2=`grep -v ">" sblast_edited_consensus_orf.fa | wc | awk '{print$3-$1}'`
sp=`tail -1 sblast_edited_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$m2" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
#Check ORF presence absence
o=`grep -v ">" sblast_edited_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $o == "" ]] ; then p="❌";else p="✅"; fi
echo $j $gr $a $b $c $d $e $f $ewah $g $h $k $l $m $n $oab $ins $del $qc1 $qc2 $sp $p
unset j gr a qry b c d e f ewah g h k l m n oab ins del qc1 qc2 m2 sp o p
cd $4
done < <(grep -if <(grep -if <(cut -f1,2 -d "_" sblast_subjects) $3 | awk '{print$2}' ) /home/neo/bird_db1/aswin/database_details/all_genome_paths) | sed "1i Species Group SRA_size Query_L gedit_L sblast_L $sblast_exons No_hits(exon) #_Almost_hits START Depth_range Poor_depth Variants Frameshifting_vars Gaps other_anonym Ins Del %_Query_cov_sblast %_Query_cov_sedit 1st_STOP_in_gene ORF" | column -t
