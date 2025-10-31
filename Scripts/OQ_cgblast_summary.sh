#2.Gblast summary table - 14.373s

#Customized headings for gblast summary
exons=`echo $2`
#exons=`grep ">" query/$2 -c`
length_and_percent_query=$(for i in `seq $exons`; do echo "E"$i"_L\tE"$i"_%q\t"; done | paste -s -d "" | sed 's/\\t$//g')
splice_sites=$(for i in `seq $(expr $exons - 1 | awk '{print$1*2}')`; do echo "SS_"$i"\t";done | paste -s -d "" | sed -e 's/^/SBG\\t/g' -e 's/$/SAG/g')
exon_start_end=$(for i in `seq $exons`; do echo "E"$i"_strt\tE"$i"_end\t"; done | paste -s -d "" | sed 's/\\t$//g')

#Gblast summary
while read i
do
#Species of interest
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
#Taxonomy order of the species
gr=`echo $j | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} /home/neo/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $i | awk '{print$2}')/aswin/$1/
#Query used
q1=`grep $j $4/gblast_query_subject_set | awk '{print$1}'`
q=`ls *.fa | grep $q1 | grep -v "flanking"`
#Exons with no hits
nh=`grep "No hits found" test.outfmt3 -B5 | grep "Query=" | awk '{print$NF}' | paste -s -d "," | sed 's/exon_//g'`
if [[ $nh == "" ]]; then nh="-"; else :; fi
#No hits in tblastx
nht=`grep -if <(grep -v "Query" tblastx.out | awk NF | awk '!a[$1]++ {print$1}') <(grep ">" $q | tr -d ">" ) -v | sed 's/exon_//g' | paste -s -d ","`
if [[ $nht == "" ]]; then nht="-"; else :; fi
#Number of exons with almost complete (>90% covered) query hit
ewah=`awk '$13>=90 {print$1}' best_hits | wc -l`
if [[ $nh == "" ]]; then nh="-"; else :; fi
#Automated consensus
a=`grep -v ">" gblast_auto_consensus.fa | paste -s -d "" | wc | awk '{print$3-$1}'`
#Exon hit length
a1=$(for z in `grep ">" $q | tr -d ">"`;do z1=`grep -w $z best_hits | awk '{print$4}'`; if [[ $z1 == "" ]]; then z1="-"; else :;fi; z2=`grep -w $z best_hits | awk '{print$13}'`; if [[ $z2 == "" ]]; then z2="-"; else :;fi; b=`echo -e "$z1\t$z2"`; echo $b; done | paste -s -d " ")
#Total query % covered = total query covered / total query length * 100
tl=`cat $q | grep -v ">" | wc | awk '{print$3-$1}'`
qc=`awk 'NR>1{print$6-$5+1}' best_hits | awk '{sum+=$1;} END{print sum;}' | awk -v l="$tl" '{print($1/l)*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
#Splice sites & sites just before & after the gene
ss=`cat splice_sites | paste -s -d " " | tr [:lower:] [:upper:]`
if [[ $ss == "" ]]; then ss=`expr $exons | awk '{print$1*2}' | xargs -I {} sh -c 'yes " -" | head -{} | xargs'`; else :; fi
#Exon start & end
a2=$(for y in `grep ">" $q | tr -d ">"`;do y1=`grep $y -A1 -w gblast_auto_consensus.fa | tail -1 | cut -c -3`; if [[ $y1 == "" ]]; then y1="-"; else :;fi; y2=`grep $y -w -A1 gblast_auto_consensus.fa | tail -1 | rev | cut -c -3 | rev`; if [[ $y2 == "" ]]; then y2="-"; else :;fi; c=`echo -e "$y1\t$y2"`; echo $c; done | paste -s -d " ")
#Check location
d=`sed 1d best_hits | awk '!seen[$1]++' | awk '{print$2}' | sort -u | wc -l`
if [ $d == "1" ]; then a3="same"; else a3="diff"; fi
#Check strandedness
e=`sed 1d best_hits | awk '!seen[$1]++' | awk '{print$NF}' | sort -u | wc -l`
if [ $e == "1" ]; then a4="same"; else a4="diff"; fi
#Paralogs
if [[ -s paralogs ]]; then pg=`cat paralogs | wc -l`; else pg="-"; fi
#duplicated exons
de=`awk '{print$1}' duplicated_exons | sort -u | paste -s -d "," | sed 's/exon_//g'`
if [[ $de == "" ]]; then de="-"; else :; fi
#Pseudogenes if GFF available
if [[ $(grep "exon" pseudogenes) == "" ]]; then psg="no"; else psg="yes"; fi
#Count deletions
a5=`cat test.out | sed 1,2d | awk '!seen[$1]++' |  awk '{print$NF}' | tr -cd '-' | wc | awk '{print$3-$1}'`
#Count Gaps
a6=`cat pairwise_exon_* | awk '/^extracted/ {print$3}' | grep -i "N" -aob | wc -l`
#Splice site disruptions
s1=`awk '{print"exon_"++i,toupper($1)}' splice_sites | sed 1d | awk '$2!="AG" {print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print $0}' | sed 's/exon_//g'`
s2=`awk '{print"exon_"++i,toupper($2)}' splice_sites | sed '$d' | awk '$2!="GT"' | awk '$2!="GC" {print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print$0}' | sed 's/exon_//g'`
#Protein consensus
pc=`grep -v ">" gblast_auto_consensus_orf.fa | paste -s -d "" | wc | awk '{print$3-$1}'`
# % of the gene in which the first Stop codons appears
sp=`tail -1 gblast_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$pc" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
o=`grep -v ">" gblast_auto_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $o == "" ]] ; then orf="❌";else orf="✅"; fi
echo -e "$j\t$gr\t$q\t$nht\t$nh\t$ewah\t$a1\t$ss\t$a2\t$a3\t$a4\t$de\t$pg\t$psg\t$a5\t$a6\t$s1\t$s2\t$a\t$pc\t$qc\t$sp\t$orf"
cd $4
unset j gr q1 q nht nh ewah a a1 a2 z z1 z2 b ss a2 y y1 y2 c d a3 e a4 de pg psg a5 a6 s1 s2 pc o tl qc sp orf
done < <(grep -if <(awk '{print$2}' $3) /home/neo/bird_db1/aswin/database_details/all_genome_paths) | sed "1i Species\tGroup\tQuery\tNo_hits_tblastx\tNo_hits_blastn\t#_Almost_hits\t$length_and_percent_query\t$splice_sites\t$exon_start_end\tLoc\tStrand\tExon_dups\tParalogs\tPseudogenes\tDel\tGaps\tSS1_disrupting_exon\tSS2_disrupting_exon\tCons\tProtein_cons\tQuery_cov\t1st_STOP_in_gene\tORF\n" | column -t

#Unset heading variables
unset length_and_percent_query splice_sites exon_start_end
