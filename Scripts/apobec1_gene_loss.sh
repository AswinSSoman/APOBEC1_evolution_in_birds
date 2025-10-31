##################################################################################################################################################################################################################################################################################################################
#IDENTIFY ALL BIRD SPECIES WITH APOBEC1 GENE LOSS USING CUSTOM PIPELINE
##################################################################################################################################################################################################################################################################################################################

#=================================================================================================================================================================================================================================================================================================================
#Create main folder for all results
mkdir APOBEC1
cd APOBEC1

##################################################################################################################################################################################################################################################################################################################
#1. DOWNLOAD & FILTER ANNOTATED QUERIES FROM NCBI
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#1.1. Get ortholog list of gene in birds

#1.1.1. From NCBI
#-------------------
#Fetch all birds with annotation
esearch -db gene -query "birds[ORGN] AND APOBEC1[Gene]" | efetch -format docsum  > apobec1_ncbi_birds_esummary_raw
cat apobec1_ncbi_birds_esummary_raw | xtract -pattern DocumentSummary -tab "\t" -def "N/A" -element Name ScientificName CommonName Description Status Summary Chromosome OtherAliases ExonCount \
| sed 's/ /_/g' | sed '1i Gene Scientific_name Common_name Description Status Summary Chromosome Othe_aliases Exon_count' | sed 's/\(DISCONTINUED\)[^\t]\+/\1/g' | column -t > APOBEC1_ncbi_birds_esummary_table

#Filter only curated genes & fetch exon-count & protein length info of them
cat APOBEC1_ncbi_birds_esummary_table | awk '$1=="APOBEC1" && $6!~"DISCONTINUED"' | awk '{print$2}' | xargs -n1 sh -c 'echo ">"$0; esearch -db gene -query "$0[ORGN] AND APOBEC1[GENE]" | efetch -format gene_table' > apobec1_ncbi_birds_basic_info_raw
egrep "^>|length" apobec1_ncbi_birds_basic_info_raw | egrep "^>|protein" | awk 'BEGIN{cmd="sort -V | head -1"} /^>/{close(cmd); print; next} {print | cmd}' | sed 's/^.*\([0-9]\+\) coding  exon.*AA length: \([0-9]\+\)/\1 \2/g' | awk -v RS=">" '{$1=$1}1' | column -t > APOBEC1_ncbi_birds_ortholog_info

#Download gene info from human as a reference
esearch -db gene -query "Homo_sapiens[ORGN] AND APOBEC1[GENE]" | efetch -format gene_table | egrep "^>|length" | egrep "^>|protein" | awk 'BEGIN{cmd="sort -V | head -1"} /^>/{close(cmd); print; next} {print | cmd}' | sed 's/^.*\([0-9]\+\) coding  exon.*AA length: \([0-9]\+\)/\1 \2/g' | awk -v RS=">" '{$1=$1}1' | sed 's/^/Homo_sapiens /g' | column -t > human_ncbi_apobec1

#1.1.2. From Ensembl
#----------------------
#Fetch all birds with annotation
wget -q --header='Content-type:text/xml' 'https://rest.ensembl.org/homology/symbol/human/APOBEC1?format=full;type=orthologues;sequence=none' -O - > apobec1_ensembl_orthologs.xml
xml2 < apobec1_ensembl_orthologs.xml | 2csv homologies @species @taxonomy_level @type @id @perc_id | sed 's/,/ /g' | column -t | grep Amniota | grep ortholog_one2one | sed 's/^./\U&/g' | sed '1i Species Group Ortholog_type Gene_ID %_identity' | column -t > APOBEC1_ensembl_aminotes_one2one_orthologs

#Fetch transcript table info of all APOBEC1 annotated birds
for sp in `awk 'NR>1{print$4}' APOBEC1_ensembl_aminotes_one2one_orthologs`
do
met=`grep $sp APOBEC1_ensembl_aminotes_one2one_orthologs | awk '{print$1}'`
echo $sp $met
wget -q --header='Content-type:text/xml' "https://rest.ensembl.org/lookup/id/${sp}?expand=1;format=full" -O - | xml2 | 2csv Transcript @Parent @display_name @biotype @length | awk -F "," '{if($2=="") print$1,"-",$3,$4; else print$0}' | sed 's/,/ /g' | awk '{print$NF,$0}' | sort -k1 -nr | awk '!($1="")' | sed 's/^/- - /g'
done | sed '1i Gene_ID Species Transcript_ID Symbol Biotype Protein_length' | column -t > apobec1_ensembl_transcript_table

#Fetch the longest transcript in each bird species
for sp in `awk 'NR>1{print$4}' APOBEC1_ensembl_aminotes_one2one_orthologs`
do
met=`grep $sp APOBEC1_ensembl_aminotes_one2one_orthologs | awk '{print$1}'`
echo $met $sp
wget -q --header='Content-type:text/xml' "https://rest.ensembl.org/lookup/id/${sp}?expand=1;format=full" -O - | xml2 | 2csv Transcript @Parent @display_name @biotype @length | sed 's/,/ /g' | awk '{print$NF,$0}' | sort -k1 -nr | awk '!($1="")' | awk 'NR==1{print$1,$NF}'
done | paste -d " " - - | sed 's/^./\U&/g' | column -t > APOBEC1_ensembl_longest_transcript_length

#1.1.3. Move all supplememtary data fetched to a seperate folder; keep only main tables in mainfolder
#-----------------------------------------------------------------------------------------------------
mkdir supplementary
mv apobec1_* supplementary/

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#1.2. Query QC

#Before downloading exon-wise sequence remove divergent queries (in terms of gene length based on IQR method)
#Filter queries that are common between NCBI annotated & genome presence in local database
grep -if <(awk '{print$1}' /home/neo/bird_db1/aswin/database_details/all_genome_paths | cut -f1,2 -d "_") <(awk '$2!=""' APOBEC1_ncbi_birds_ortholog_info) | awk '{print$1,$NF}' | column -t > common_queries_ncbi
grep -if <(awk '{print$1}' /home/neo/bird_db1/aswin/database_details/all_genome_paths | cut -f1,2 -d "_") <(awk '$4!=""' APOBEC1_ensembl_longest_transcript_length) | awk '{print$1,$NF}' | column -t > common_queries_ensembl

#Remove outliers
Rscript remove_outliers.r common_queries_ncbi T | awk 'NR>1{print$2,$3}' | column -t  > length_filtered_queries_ncbi
Rscript remove_outliers.r common_queries_ensembl T | awk 'NR>1{print$2,$3}' | column -t  > length_filtered_queries_ensembl

#Plot distribution before & after query filtering
#NCBI
paste <(awk '{print$2}' common_queries_ncbi | ministat | sed 's/x <stdin>/\n\t\tQuery-length distribution before filtering\n/g' | GREP_COLORS="mt=33;4" grep "Query-length.*[a-z]\|$" --color=always) \
	<(awk '{print$2}' common_queries_ncbi | ministat | awk 'END{print NR-1}' | xargs -n1 sh -c 'yes " ❙ " | head -$0 ' | sed '1i \\n\n' | GREP_COLORS="mt=33" grep "❙\|$" --color=always) \
	<(awk '{print$2}' length_filtered_queries_ncbi | ministat | sed 's/x <stdin>/\n\t\t\t\tQuery-length distribution after filtering\n/g' | GREP_COLORS="mt=33;4" grep "Query-length.*[a-z]\|$" --color=always) > query_length_plot_ncbi
paste <(termgraph <(awk '{print$1,$NF}' common_queries_ncbi) --format '{:.0f}' | sed 's/:/ :/g' | column -t) \
	<(termgraph <(awk '{print$1,$NF}' common_queries_ncbi) --format '{:.0f}' | awk 'END{print NR-2}' | xargs -n1 sh -c 'yes " ❙ " | head -$0' | GREP_COLORS="mt=33" grep "❙\|$" --color=always) \
	<(termgraph <(awk '{print$1,$NF}' length_filtered_queries_ncbi) --format '{:.0f}' | sed 's/:/ :/g' | column -t) | column -t | sed '1i \\n' >> query_length_plot_ncbi

#ensembl
paste <(awk '{print$2}' common_queries_ensembl | ministat | sed 's/x <stdin>/\n\t\tQuery-length distribution before filtering\n/g' | GREP_COLORS="mt=33;4" grep "Query-length.*[a-z]\|$" --color=always) \
	<(awk '{print$2}' common_queries_ensembl | ministat | awk 'END{print NR-1}' | xargs -n1 sh -c 'yes " ❙ " | head -$0 ' | sed '1i \\n\n' | GREP_COLORS="mt=33" grep "❙\|$" --color=always) \
	<(awk '{print$2}' length_filtered_queries_ensembl | ministat | sed 's/x <stdin>/\n\t\t\t\tQuery-length distribution after filtering\n/g' | GREP_COLORS="mt=33;4" grep "Query-length.*[a-z]\|$" --color=always) > query_length_plot_ensembl
paste <(termgraph <(awk '{print$1,$NF}' common_queries_ensembl) --format '{:.0f}' | sed 's/:/ :/g' | column -t) \
	<(termgraph <(awk '{print$1,$NF}' common_queries_ensembl) --format '{:.0f}' | awk 'END{print NR-2}' | xargs -n1 sh -c 'yes " ❙ " | head -$0' | GREP_COLORS="mt=33" grep "❙\|$" --color=always) \
	<(termgraph <(awk '{print$1,$NF}' length_filtered_queries_ensembl) --format '{:.0f}' | sed 's/:/ :/g' | column -t) | column -t | sed '1i \\n' >> query_length_plot_ensembl

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#1.3. Download exon_wise sequence of filtered queries

mkdir exon_wise_sequence
cd exon_wise_sequence/

#1.3.1. From NCBI
#---------------------
#For genes with very short exons use "edirect" method rather than "datasets"
awk '!/Name/ {print$1}' ../length_filtered_queries_ncbi | xargs -n1 sh -c 'exon_wise -edirect APOBEC1 $0 -ir'

#1.3.2. From Ensembl
#---------------------
while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
k=`echo $i | awk '{print$2}'`
l=`echo $i | awk '{print$3}'`
echo ">"$j
ensembl_exon_wise $k $l -ir
rename "s/$k/APOBEC1_$j/g" *report
rename "s/$k/APOBEC1_$j/g" *.fa
done < <(grep -if <(awk '{print$1}' ../length_filtered_queries_ensembl) ../APOBEC1_ensembl_longest_transcript_length)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#1.4. QC downloaded queries

#1.4.1. For NCBI downloaded genes
#-------------------------------------

#1.4.1.a. Metadata matching
for i in `ls *.fa | grep -v "ensembl"`
do
r=`echo $i | sed 's/.fa/_report/g'`
#gene symbol
j1=`echo $i | cut -f1 -d "_"`
j2=`awk 'NF' $r | awk 'NR==1{print$1}'`
#Species name
k1=`echo $i | cut -f2- -d "_" | sed 's/.fa//g'`
k2=`awk 'NF' $r | awk 'NR==1{print$0}' | sed -e 's/^.*\(\[.*\]\).*$/\1/g' -e 's/ /_/g' | tr -d "[]"`
#No of Exons from fetched metadata
ex2=`grep $k1 ../APOBEC1_ncbi_birds_ortholog_info | awk '{print$2}' | cut -f1 -d "_"`
#No of Exons from downloaded sequence
ex1=`grep ">" $i -c`
#Number of isoforms
l=`grep "isoforms present" $r | cut -f1 -d " "`
#Transcript type
m=`grep "^Longest isoform" $r | awk '{for(z=1;z<=NF;z++) if($z~/,/) print $z}' | awk -F "_" 'NR==1{print$1}'`
#CDS length
n1=`grep -v ">" $i | wc | awk '{print$3-$1}'`
n2=`grep -A1 "Longest isoform" -A1 $r | tail -1 | awk -F ":" '{print($NF+1)*3}'`
#start stop site
p=`grep -v ">" $i | paste -s -d "" | cut -c -3 | tr '[:lower:]' '[:upper:]'`
if [[ $p == "" ]]; then p="-"; else :; fi
q=`grep -v ">" $i | paste -s -d "" | rev | cut -c -3 | rev | tr '[:lower:]' '[:upper:]'`
if [[ $q == "" ]]; then q="-"; else :; fi
#check orf
s=`transeq <(grep -v ">" $i) -auto -stdout | grep -v ">" | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $s == "" ]] ; then t="❌";else t="✅"; fi
echo $i $j1 $j2 $k1 $k2 $ex1 $ex2 $l $m $n1 $n2 $p $q $t
unset r j1 j2 k1 k2 ex1 ex2 l m n1 n2 p q s t
done | sed '1i Gene Symbol_1 Symbol_2 Species_1 Species_2 #_Exons_1 #_Exons_2 #_Isoforms transcript_type CDS_length1 CDS_length2 START STOP ORF' | column -t > qc_fetched_data_ncbi

#1.4.1.b. Filter queries based on complete/incomplete ORF & human exon count
exon_wise -edirect APOBEC1 Homo_sapiens -ir
human_exon_count=`awk '{for(i=1;i<=NF;i++) if($i~/coding/) print $(i-1)}' APOBEC1_Homo_sapiens_report`
awk -v h="$human_exon_count" '{if($NF=="✅" && $6==$7 && $7==h) print$1}' qc_fetched_data_ncbi > orf_filtered_queries_ncbi

#1.4.2. For Ensembl downloaded genes
#--------------------------------------

#1.4.2.a. Metadata matching
for i in `ls *.fa | grep "ensembl"`
do
g_s=`echo $i | sed 's/_ensembl.fa//g'`
r=`echo $i | sed 's/.fa/_report/g'`
#gene symbol
j1=`echo $i | cut -f1 -d "_"`
j2=`cat $r | xml2 | grep display_name | grep -v "Transcript" | awk -F "=" '{print$NF}' | sort -u`
if [[ $j2 == "" ]]; then j2="-"; else :; fi
#Species name
k1=`echo $i | cut -f2- -d "_" | sed 's/_ensembl.fa//g'`
k2=`cat $r | xml2 | grep species | awk -F "=" '{print$NF}' | sort -u | sed 's/^./\U&/g' | cut -f1,2 -d "_"`
if [[ $k2 == "" ]]; then k2="-"; else :; fi
#Exon count
l=`grep ">" $i -c`
#Transcript length
m1=`cat $i | grep -v ">" | wc | awk '{print$3-$1}'`
m2=`cat $r | xml2 | 2csv Transcript @length | sort -nr | awk 'NR==1{print$1*3+3}'`
if [[ $m2 == "" ]]; then m2="-"; else :; fi
#start stop site
p=`grep -v ">" $i | paste -s -d "" | cut -c -3 | tr '[:lower:]' '[:upper:]'`
if [[ $p == "" ]]; then p="-"; else :; fi
q=`grep -v ">" $i | paste -s -d "" | rev | cut -c -3 | rev | tr '[:lower:]' '[:upper:]'`
if [[ $q == "" ]]; then q="-"; else :; fi
#check orf
s=`transeq <(grep -v ">" $i) -auto -stdout | grep -v ">" | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $s == "" ]] ; then t="❌";else t="✅"; fi
echo $g_s $j1 $j2 $k1 $k2 $l $m1 $m2 $p $q $t
unset g_s r j1 j2 k1 k2 l m1 m2 p q s t
done | sed '1i Gene Symbol_1 Symbol_2 Species_1 Species_2 Exon_count CDS_length1 CDS_length2 START STOP ORF' | column -t > qc_fetched_data_ensembl

#Add exon count & protein length info
awk 'NR>1{print$4,$6,$7/3-1}' qc_fetched_data_ensembl | column -t > ../APOBEC1_ensembl_birds_ortholog_info

#1.4.2.b. Filter queries based on complete/incomplete ORF & human exon count
awk -v h="$human_exon_count" '{if($NF=="✅" && $7==$8 && $6==h) print$1}' qc_fetched_data_ensembl | sed 's/$/_ensembl.fa/g' > orf_filtered_queries_ensembl

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#1.4.3. Protein alignment w.r.t humans
transeq <(ls *.fa | xargs -n1 sh -c 'awk "BEGIN{print \"\n>\"ARGV[1]} !/>/{printf\$0}" $0' | awk NF) -auto -stdout > protein.aa
clustalo -i protein.aa -t protein --wrap=250 --outfmt=clu --resno > protein.aln

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#1.4.4 Supporting evidence of queries for NCBI downloaded genes

cat orf_filtered_queries_ncbi | sed 's/.fa//g' | sed 's/_/ /1' | xargs -n2 sh -c 'echo ">"$0 $1; esearch -db nuccore -query "$0[GENE] AND $1[ORGN]" | efetch -format gb -mode xml' > raw_file_supporting_evidence
cat raw_file_supporting_evidence | egrep "^>|Supporting" | tr "<" "\n" | sed 's/.*Supporting/Supporting/g' | grep -v "/GBQualifier_value" | awk NF | sed 's/ /_/g' | sed '/>/! s/^/- /g' | column -t > supporting_evidence_ncbi
awk '/100%/{print$1}' RS=">" supporting_evidence | awk NF > queries_with_100%_rna_cov
#cat queries_with_100%_rna_cov | sed 's/$/.fa/g' > orf_filtered_queries

#Combine info of ncbi & ensembl
cd ../
cat *ortholog_info | awk '!a[$1]++' | awk '$2' | column -t > APOBEC1_ortholog_info
grep -if <(grep -if <(grep -if <(cat exon_wise_sequence/orf_filtered_queries_* | cut -f2,3 -d "_" | sed 's/.fa//g') APOBEC1_ortholog_info | awk '{print$1}') <(ls exon_wise_sequence/*.fa) | sed 's/_ensembl/ ensembl/g' | sed 's/.fa//g' | tr "/" " " | column -t | sort -k2 -r | awk '!a[$2]++' | awk '!($1="")' | column -t | sed 's/ \+/_/g') <(ls exon_wise_sequence/*.fa) | awk -F "/" '{print$2}' > qced_queries

##################################################################################################################################################################################################################################################################################################################
#2. VALIDATE THE FILTERED QUERY SET AND CREATE A CURATED SET OF QUERIES

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#2.1. gblast only in species with NCBI annotated query gene against it's own genome

while read i
do
j=`echo $i | awk '{print$1}'`
cd $(echo $i | awk '{print$2}')
#Create folders for results
mkdir -p aswin/APOBEC1
cd aswin/APOBEC1
#Query used
qs=`grep $j /home/neo/bird_db1/aswin/APOBEC1/qced_queries`
q=`find /home/neo/bird_db1/aswin/APOBEC1/exon_wise_sequence/ -name "$qs" | awk -F "/" '{print$NF}'`
echo " ▶ " $j "(G) " $q "(Q)"
#Get the query sequence from same species
find /home/neo/bird_db1/aswin/APOBEC1/exon_wise_sequence/ -name "$qs" -exec cp {} . \;
#Get the subject genome
g=`find ../../ -maxdepth 2 -name "GC*.fna"`
a=`find ../../ -maxdepth 2 -name "*genomic.gff"`
#genome blast with it's own query
if [[ $a == "" ]]; then
#blastn
gblast_short $g $q -evalue=0.01 -word_size=11 -fix_query -tblastx=no -iflank
else
#blastn
gblast_short $g $q -evalue=0.01 -word_size=11 -fix_query -tblastx=no -iflank $a
fi
#Generate ORF
transeq <(grep -v ">" gblast_auto_consensus.fa | paste -s -d "" | sed "1i \>$j") gblast_auto_consensus_orf.fa
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" gblast_auto_consensus.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" gblast_auto_consensus_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > gblast_translated
#For better visualization
gnaa=`head -1 gblast_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" gblast_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 gblast_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > gblast_nucleotide_amino_acid_view
#cat <(cut -c -251 gblast_translated) <(cut -c 253-504 gblast_translated) <(cut -c 505- gblast_translated) | GREP_COLORS="mt=41" egrep "^ M |$| \* |$| X |$" --color=always > gblast_nucleotide_amino_acid_view
#Check splice sites & save it into a file
for z in `grep ">" $q | tr -d ">"`
do
f1=`ls pairwise_exon_*aln | grep "$z\b"`
if [[ $f1 == "" ]]; then ss1="-"; ss2="-"
else
#Splice site 1
ss1c=`cat $f1 | grep "^extract" -B2 | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":"`
ss1=`cat $f1 | grep "^extract" -B2 | grep "|" -C1 | head -3 | tail -1 | cut -c -$ss1c | rev | cut -c -2 | rev | sed 's/ \+/-/g'`
if [[ $ss1 == "" ]]; then ss1="-"; else :;fi
#Splice site 2
ssc2=`cat $f1 | grep "^extract" -B2 | grep "|" -C1 | tail -3 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+2}'`
ss2=`cat $f1 | grep "^extract" -B2 | grep "|" -C1 | tail -3 | tail -1 | cut -c $ssc2- | cut -c -2 | sed 's/ \+/-/g'`
if [[ $ss2 == "" ]]; then ss2="-"; else :;fi
fi
echo -e "$ss1\t$ss2"
unset f1 ss1c ss2c ss1 ss2
done > splice_sites
#Unset temporary variables
unset j q qs s g a gnaa z
cd /home/neo/bird_db1/aswin/APOBEC1/
printf '_%.s' {1..150}; printf '\n\n'
done < <(grep -if <(cut -f2,3 -d "_" qced_queries | sed 's/.fa//g') /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#2.2. validation gblast summary
while read i
do
j=`echo $i | awk '{print$1}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
#Taxonomy order of the species
k=`echo $j | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} /home/neo/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
#Query used
qs=`grep $j /home/neo/bird_db1/aswin/APOBEC1/qced_queries`
l=`find /home/neo/bird_db1/aswin/APOBEC1/exon_wise_sequence/ -name "$qs" | awk -F "/" '{print$NF}'`
#Query length
m=`grep -v ">" $l | wc | awk '{print$3-$1}'`
#Query ORF
n=`grep -v ">" gblast_auto_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $n == "" ]] ; then n1="❌";else n1="✅"; fi
#Check location
o=`sed 1d best_hits | awk '!seen[$1]++' | awk '{print$2}' | sort -u | wc -l`
if [ $o == "1" ]; then o1="same"; else o1="diff"; fi
#Check strandedness
p=`sed 1d best_hits | awk '!seen[$1]++' | awk '{print$NF}' | sort -u | wc -l`
if [ $p == "1" ]; then p1="same"; else p1="diff"; fi
#Count Gaps
q=`cat gblast_auto_consensus.fa | awk '/^extracted/ {print$3}' | grep -i "N" -aob | wc -l`
#Splice site disruptions
s1=`awk '{print"exon_"++i,toupper($1)}' splice_sites | sed 1d | awk '$2!="AG" {print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print $0}' | sed 's/exon_//g'`
s2=`awk '{print"exon_"++i,toupper($2)}' splice_sites | sed '$d' | awk '$2!="GT"' | awk '$2!="GC" {print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print$0}' | sed 's/exon_//g'`
#Total query % covered = total query covered / total query length * 100
r=`awk 'NR>1{print$6-$5+1}' best_hits | awk '{sum+=$1;} END{print sum;}' | awk -v z="$m" '{print($1/z)*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
#Hit length
s=`grep -v ">" gblast_auto_consensus.fa | wc | awk '{print$3-$1}'`
#validate ORF
t=`grep -v ">" gblast_auto_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $t == "" ]] ; then t1="❌";else t1="✅"; fi
echo $j $k $l $m $n1 $o1 $p1 $q $s1 $s2 $r $s $t1
unset j k qs l m n n1 o o1 p p1 q s1 s2 r s t t1
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if <(cut -f2,3 -d "_" qced_queries | sed 's/.fa//g') /home/neo/bird_db1/aswin/database_details/all_genome_paths) | sed '1i Species Group Query Q_len Q_ORF Location Strand Gaps SS1_disruption SS2_disruption %_Q_cov Cons_len Cons_ORF' | column -t > validation_gblast_summary

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#2.3. Create set of validated queries & test subjects for Comprehensive gblast

#List of queries that gave complete ORF with 100% query coverage & no gaps from genome blast
awk '$5=="✅" && $8=="0" && $9=="-" && $10=="-" && $11~"100" && $NF=="✅" {print$3}' validation_gblast_summary | sed 's/_ensembl//g' > annotated_queries
#Create seperate tables for validated queries : one exclusive from each category (ncbi annotated, gblast, gedit, sblast) & other is for a combined list(total validated queries)
#Add validated queries to total validated queries (a seperate list to append)
cat annotated_queries > total_validated_queries
mkdir validated_sequences
mkdir validated_sequences_with_flanking_regions
#add sequence to folder
while read i
do
j=`echo $i | awk '{print$1}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
cat gblast_auto_consensus.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/"APOBEC1_"$j.fa
cat extracted_flanking_region.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/"APOBEC1_"$j"_flanking".fa
unset j
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if <(cut -f2,3 -d "_" annotated_queries | sed 's/.fa//g' ) /home/neo/bird_db1/aswin/database_details/all_genome_paths)

##################################################################################################################################################################################################################################################################################################################
#3. Visualize validated complete ORFs on the phylogenetic tree of birds

cp /home/neo/bird_db1/aswin/tree/species_with_genome_wo_cygnus_olor_tree gene_tree
cut -f2- -d "_" annotated_queries | sed 's/.fa//g' | xargs -n1 sh -c 'sed "/$0/ s/$/annotated/g" -i gene_tree'

##################################################################################################################################################################################################################################################################################################################
#4. COMPREHENSIVE GENOME BLAST

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Subjects for comprehensive genome blast
grep -if <(cat annotated_queries | cut -f2- -d "_" | sed 's/.fa//g') /home/neo/bird_db1/aswin/database_details/all_genome_paths -v | awk '{print$1}' > gblast_subjects

#Create all possbible pairwise distances between all species pairs - Create once & keep it in a path don't need eveytime you run the pipeline
#Query subject set
for i in `cat gblast_subjects | cut -f1,2 -d "_"`; do j=`grep -if <(sed 's/.fa//g' annotated_queries | cut -f2- -d "_") <(awk -v s="$i" '$2==s' /home/neo/bird_db1/aswin/tree/species_with_genome_data/all_possible_pairwise_distances | sort -k3 -n | awk '{print$1}') | head -1`; echo $j $i; unset j; done | column -t > gblast_query_subject_set

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#4.1. Run Comprehensive gblast for each Query-Subject set (34 min)

for sp in `awk '{print$1}' gblast_query_subject_set | sort -u`
do
echo " ▶ " $sp "(Q)"
printf '=%.s' {1..200}; printf '\n\n'
./OQ_cgblast.sh APOBEC1 $sp /home/neo/bird_db1/aswin/APOBEC1/
done
unset sp

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#4.2. gblast summary

#Summary table of hits have information in exon_wise manner - hence for proper alignment of exon summary make different summary tables if the exon count of validated quries used for gblast are varying
#Split query-subject set based on query exon-count
cp gblast_query_subject_set coding_sets
awk '$3!=""' APOBEC1_ortholog_info | awk '{print$1,$2}' | xargs -n2 sh -c 'sed "/^$0/ s/$/ $1/g" -i coding_sets'
grep -if <(cut -f2- -d "_" annotated_queries | sed 's/.fa//g') APOBEC1_ortholog_info | awk '{print$1,$1,$2}' | column -t | awk '{print >> $3"_coding_exons_set"; close($3)}'
cat coding_sets | column -t | awk '{print >> $3"_coding_exons_set"; close($3)}'
rm coding_sets

#Summary for each sets
counter=0
for set in `ls *coding_exons_set | sort -V`
do
counter=$((counter + 1))
exon_count=`echo $set | cut -f1 -d "_"`
./OQ_cgblast_summary.sh APOBEC1 $exon_count $set /home/neo/bird_db1/aswin/APOBEC1/ > gblast_summary"$counter" 2>/dev/null
unset exon_count
done
unset set counter

#Since the summary table have all information : reduce the dimensionality to increase the interpretability but at the same time minimize the information loss.
#principle compomemets(PC's) here refers to the information that contribuites the most to tthe interpretation of gene presence/loss

#From Gblast : PC's of all species : choose PC's based on univeral headings (not dependent of characteristics of genes such as number of exons)
echo -n > gblast_principal_components
for i in `ls gblast_summary[0-9]* | sort -V`
do
awk 'NR==1{for(i=1;i<=NF;i++)if($i~/Species|Group|Query|No_hits|#_Almost_hits|Loc|Strand|Exon_dups|Paralogs|Pseudogenes|Del|Gaps|Cons|Query_cov|ORF|1st_STOP_in_gene|E1_strt|Protein_cons|SS1_disrupting_exon|SS2_disrupting_exon/)f[n++]=i}{for(i=0;i<n;i++)printf"%s%s",i?" ":"",$f[i];print""}' $i | column -t >> gblast_principal_components_tmp
done
awk '!a[$1]++' gblast_principal_components_tmp | column -t > gblast_principal_components
rm gblast_principal_components_tmp

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#4.3. Add completely obtained ORFs to the total validated queries list, gene tree & add their sequences to the exon_wise folder

#Add validated queries from Comprehensive gblast to the validated list
awk '$4=="-" && $14=="0" && $15=="-" && $16=="-" && $19=="100" && $NF=="✅" {print$1}' gblast_principal_components > gblast_validated_queries
cat gblast_validated_queries | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' >> total_validated_queries

#Add to gene tree
for i in `cat gblast_validated_queries`; do sed "/$i/ s/$/"gblast_"$i/g" -i gene_tree; unset j;done

#Add validated sequence to a folder
while read i
do
j=`echo $i | awk '{print$1}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
cat gblast_auto_consensus.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/"APOBEC1_"$j.fa
cat extracted_flanking_region.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/"APOBEC1_"$j"_flanking".fa
unset j
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if gblast_validated_queries /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#4.4. QC validated sequences

./qc_validated.sh validated_sequences validated_sequences_with_flanking_regions *_coding_exons_set APOBEC1 /home/neo/bird_db1/aswin/APOBEC1 gblast_validated_qc

#Find potential false positives in validated sequences by detecting outliers in nucleotide MSA
sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" gblast_validated_qc | awk -F ":" '/Outlier/ {print$2}' | sed 's/Nothing//g' | tr "," "\n" | awk NF | sed 's/^[ ]\+//g' > outliers

#Remove outliers from the validated sequence records
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers gblast_validated_queries
cat outliers | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' > outliers2
if [[ -s outliers2 ]]; then 
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers2 total_validated_queries
else :; fi
cat outliers | xargs -n1 sh -c 'sed "/$0/ s/[ ]\+[^ ]\+$//g" gene_tree -i'
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa) | xargs rm
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/*.fa) | xargs rm
rm outliers outliers2

##################################################################################################################################################################################################################################################################################################################
#5. Synteny based on Gblast

#5.1. Print focal genes it's 3 neighbouring genes
#--------------------------------------------------------
while read i
do
j=`echo $i | awk '{print$1}'`
gr=`echo $j | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} /home/neo/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
s=`cat Synteny 2>/dev/null | grep -v "dup" | awk '!($1="")($2="")' | awk '{for(i=1;i<=NF;i++) if($i~/\[01;31m/ && i>3) print$(i-3),$(i-2),$(i-1),$i,$(i+1),$(i+2),$(i+3); else if($i~/\[01;31m/ && i==3) print"-",$(i-2),$(i-1),$i,$(i+1),$(i+2),"-"; else if($i~/\[01;31m/ && i==2) print"-","-",$(i-1),$i,$(i+1),"-","-"; else if($i~/\[01;31m/ && i==1) print"-","-","-",$i,"-","-","-" }' | uniq -c | sort -k1,1 -nr | head -1 | awk '!($1="")' | column -t`
if [[ $s == "" ]]; then s="-"; else :; fi
echo -e ">$j\t$gr\t$s"
unset j gr s
cd /home/neo/bird_db1
#Fill empty columns with dashes (-)
done < /home/neo/bird_db1/aswin/database_details/all_genome_paths | awk '{if($1~/^>/) print$0; else print"-","-",$0}' | sed 's/^>//g' | awk '{if($4=="" && $3=="-") print$0,"- - - - - -"; else print$0}' | column -t > synteny_tmp

#Convert LOC-id's to gene synbols if it is fetchable from NCBI
sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" synteny_tmp | grep -i "\bLOC[^ ]*" -o | tr -d ")(\-+" | xargs -n 1 bash -c 'echo ">"$0; efetch -db gene -id $0 -format native' > locid_efetch
awk '!/replaced/ {print">"$0}' RS=">" locid_efetch | grep -v "^>$" | egrep "^>|Designations" | awk -F ";" '{print$1}' | sed 's/Other Designations: //g' | sed '/>/! s/ /_/g' | sed 's/LOW_QUALITY_PROTEIN:_//g' | awk -F "," '{print$1}' | awk -v RS=">" '{print$1,$2}' OFS=" " | awk NF > locid_fetched_names

#Manually edit the locid fetched names to shorten the genes names. eg: homeobox_protein_NANOG-like to NANOG-LIKE
cp synteny_tmp synteny_tmp2
cat locid_fetched_names | xargs -n2 sh -c 'sed -i "s/$0/$1/g" synteny_tmp2'
cat synteny_tmp2 | column -t > synteny
rm synteny_tmp synteny_tmp2
mv locid_* supplementary/

#fill all empty columns at the end & rearrange the order of genes in one direction
#Most adjacent genes
#rg=`awk '{print$7}' synteny | sort | sed 's/(.)//g' | sed 's/^-$//g' | awk NF | sort | uniq -c | sort -nr | awk 'NR==1{print$2}'`
#lg=`awk '{print$5}' synteny | sort | sed 's/(.)//g' | sed 's/^-$//g' | awk NF | sort | uniq -c | sort -nr | awk 'NR==1{print$2}'`
#cat synteny | awk '{if($4=="")print$0,"-","-","-","-","-","-"; else print$0}' | awk '{if($5=="")print$0,"-","-","-","-","-"; else print$0}' | awk '{if($6=="")print$0,"-","-","-","-"; else print$0}' | awk '{if($7=="")print$0,"-","-","-"; else print$0}' | awk '{if($8=="")print$0,"-","-"; else print$0}'  | awk '{if($9=="")print$0,"-"; else print$0}' | awk -v l="$lg" -v r="$rg" -v IGNORECASE=1 '{if($5~l || $7~r) print$0; else if($5~r || $7~l) print$1,$2,$9,$8,$7,$6,$5,$4,$3; else print$0}' | column -t > gblast_synteny
#rm synteny
#unset rg lg

#5.2. Map synteny info on tree
#--------------------------------------------
cp /home/neo/bird_db1/aswin/tree/species_with_genome_wo_cygnus_olor_tree stree
grep "[A-Z]" stree | awk '{print$NF}' | xargs -n1 sh -c 'grep $0 synteny' | column -t > stmp
sed -z 's/\n/\n\n/g' stmp -i

#check the species in tree and species with syteny is matching
#paste stree stmp | grep '[A-Z]' | sed 's/^[^A-Z]\+//g' | tr -d "@" | awk '{if($1==$2) print"yes",$0; else print"no",$0}' | column -t | grep "^no"
grep "[A-Z]" stree | awk '{print$NF}' | xargs -n1 sh -c 'grep $0 synteny' | awk '!($1="")' | column -t > stmp
sed -z 's/\n/\n\n/g' stmp -i
paste stree stmp > synteny_tree
rm stmp

#---------------------#------------------------------#----------------------_#--------------------------------#--------------------------#--------------------------------#-----------------------_#
#5.3. Chromview of all birds

while read i
do
j=`echo $i | awk '{print$1}'`
gr=`echo $j | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} /home/neo/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
cv=`cat chromview_non_overlapping_gene | sed '/>/ s/ /_/g' | paste -s -d " "`

if [[ $cv == "" ]]
then
cv=`cat chromview | sed '/>/ s/ /_/g' | paste -s -d " "`
else :
fi

if [[ $cv == "" ]]
then
cv=`cat chromview_single_non_overlapping_gene | sed '/>/ s/ /_/g' | paste -s -d " "`
else :
fi

if [[ $cv == "" ]]
then
cv="- -"
else :
fi

echo -e ">$j\t$gr\t$cv"
unset j gr cv
cd ~/bird_db1/aswin/APOBEC1
done < /home/neo/bird_db1/aswin/database_details/all_genome_paths | column -t | sed 's/❯ \+/❯ /g' > ctmp

#Map chromview on tree
grep "[A-Z]" stree | awk '{print$NF}' | xargs -n1 sh -c 'grep $0 ctmp' | column -t | sed 's/❯ \+/❯ /g' > ctmp2
sed -z 's/\n/\n\n/g' ctmp2 -i

#check the species in tree and species with syteny is matching
#paste stree ctmp2 | grep '[A-Z]' | sed 's/^[^A-Z]\+//g' | tr -d ">" | awk '{if($1==$2) print"yes",$0; else print"no",$0}' | column -t | grep "^no"
grep "[A-Z]" stree | awk '{print$NF}' | xargs -n1 sh -c 'grep $0 ctmp' | awk '!($1="")' | column -t | sed 's/❯ \+/❯ /g' > ctmp2
sed -z 's/\n/\n\n/g' ctmp2 -i
paste stree ctmp2 > chromview_tree
rm ctmp ctmp2

##################################################################################################################################################################################################################################################################################################################
#6. GEDIT  (Genome-blast based Editing/Extraction)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.1. Filter subjects to run Gedit : Choose all species that didn't form ORF from gblast

#if the exon hits with diff location & strandedness are authentic hits don't force same location & strandedness criteria
#awk '$5=="-" && $14=="0" {print$1}' gblast_principal_components > gedit_subjects
#awk '$5=="-" && $14=="0" && $8=="same" && $9=="same" {print$1}' gblast_principal_components > gedit_subjects
grep -if gblast_validated_queries gblast_subjects -v > gedit_subjects

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run Gedit
./OQ_gedit.sh APOBEC1 /home/neo/bird_db1/aswin/APOBEC1/ 2>/dev/null

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.2. Gedit summary : For comprehensive checking

counter=0
#Summary for each sets
for set in `ls *coding_exons_set`
do
counter=$((counter + 1))
exon_count=`echo $set | cut -f1 -d "_"`
./OQ_gedit_summary.sh APOBEC1 $exon_count $set /home/neo/bird_db1/aswin/APOBEC1/ > gedit_summary"$counter" 2>/dev/null
unset exon_count
done
unset set counter

#Check if gedit ran on all gedit subjects : just for conformation in case manually intervened!
#grep -if <(awk '!/Species/ {print$1}' gedit_summary*) gedit_subjects -z

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.3. Gedit Principal components(PC's) : Components from all data that contribuites the most to the interpretation/inference

#Since gedit summary have all information : reduce the dimensionality to increasing interpretability but at the same time minimizing information loss
while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
#Taxonomy order of the species
gr=`echo $j | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} /home/neo/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
#Query used
qry1=`grep $j /home/neo/bird_db1/aswin/APOBEC1/gblast_query_subject_set | awk '{print$1}'`
qry=`ls *.fa | grep $qry1".fa$"`
#Exons with no hits
enh=`grep "No hits found" test.outfmt3 -B5 | grep Query | cut -f2 -d "=" | tr -d " " | paste -s -d ","`
if [[ $enh == "" ]]; then enh="-"; else :; fi
#Number of exons with almost (>90% covered) complete hits
if [[ -f gblast_edited_query_covered ]]; then
ewah=`awk '$NF>90 {print$1}' gblast_edited_query_covered | wc -l`
else ewah="0"; fi
#Total query % covered = total query covered / total query length * 100
tl=`cat $qry | grep -v ">" | wc | awk '{print$3-$1}'`
if [[ -f gblast_edited_query_covered ]]; then
qc=`awk '{print$3}' gblast_edited_query_covered | awk '{sum+=$1;} END{print sum;}' | awk -v l="$tl" '{print($1/l)*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
else qc="0"; fi
#Gaps within exon
#ge=`grep -v ">" gblast_edited_consensus.fa | grep -i "N" -o | wc -l`
ge=`awk '/N/{print">"$0}' RS=">" gblast_edited_consensus.fa | awk '{if($1~"^>") print$1; else print $0 | "grep -i -o N | wc -l"}' | paste -s -d " " | sed 's/>exon_//g' | awk '{print$2"("$1")"}' | paste -s -d "," | awk '{if($0~"[0-9]") print$0; else print "-"}'`
#Total Gaps : exons + flanking regions
if [[ -f pairwise_exon_* ]]; then
gf=$(grep "^extracted" pairwise_exon_* -B2 | awk '{print$3}'  | grep -i "n" -o | wc -l)
else gf="-"; fi
#Other anonymous bases (number of occurances)
oab=`awk '!/>/ {gsub(/[ATGCNatgcn]/,"")}1' gblast_auto_consensus.fa | sed '/>exon/! s/[a-z]/\U&/g' | awk '/[A-Z]/' RS=">" | awk NF | paste -d " " - - | sed 's/exon_//g' | awk '{print length($2)"("$1")"}' | paste -s -d ","`
if [[ $oab == "" ]]; then oab="-"; else :; fi
#Total insertions
if [[ -s gblast_edited_insertions ]]; then
ti=`awk '{a+=$2} END{print a}' gblast_edited_insertions`
else ti="-"; fi
#Total deletions
if [[ -s gblast_edited_deletions ]]; then
td=`awk '{a+=$2} END{print a}' gblast_edited_deletions`
else td="-"; fi
#Indel frameshifting
if [[ $ti == "-" && $td == "-" ]]; then
indf="-"
else
indf=`expr $ti - $td | tr -d "-" | xargs -I {} calc {} / 3 | awk '{if($0~/\./) print "yes"; else print "no"}'`
fi
#Calculate indels at exon boundaries for each exon alignments & save it in a table in respective folders
for pe in `ls pairwise_exon_* | sort -V`
do
pen=`echo $pe | sed 's/.aln//g' | cut -f2- -d "_"`
#Trim flanking regions outside alignment
start=`grep -B2 "^extrac" $pe | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":" | awk '{print$0+1}'`
grep -B2 "^extrac" $pe | grep "|" -C1 | cut -c $start- > tmp1
end=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+1}'`
cat tmp1 | cut -c -$end > tmp2
#insertions at exon start : counted if a dash is observed within the query sequence anywhere within the first 5% of pairwise alignment between query & subject
is=`cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/query/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $is == "" ]]; then is="-"; else :; fi
#insertions at exon end : counted if a dash is observed within the query sequence anywhere after 95% of pairwise alignment between query & subject
ie=`cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/query/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ie == "" ]]; then ie="-"; else :; fi
#deletions at exon start : counted if a dash is observed within the subject sequence anywhere within the first 10% (since subject sequence is longer due to flanking regions) of pairwise alignment between query & subject
ds=`cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/genome/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ds == "" ]]; then ds="-"; else :; fi
#deletions at exon end : counted if a dash is observed within the subject sequence anywhere after 90% of pairwise alignment between query & subject
de=`cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/genome/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
#Longest continous insertions : Position, Length & % of gene in which it appears
if [[ $(head -1 tmp2 | grep "\-") == "" ]]; then lci="- - -"; else
lci=`bedtools merge -i <(cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/query/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
fi
#Longest continous deletions : Position, Length & % of gene in which it appears
if [[ $(tail -1 tmp2 | grep "\-") == "" ]]; then lcd="- - -"; else
lcd=`bedtools merge -i <(cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/genome/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
#lci=`cat tmp2 | grep "\-\+" -o | tr -d " " | awk '{print length}' | sort -nr | head -1`
fi
if [[ $de == "" ]]; then de="-"; else :; fi
echo $pen $is $ie $ds $de $lci $lcd
unset pen start end is ie ds de lci lcd
rm tmp1 tmp2
done | sed '1i Exons Ins_in_1st_5%_exon Ins_in_last_95%_exon Del_in_1st_5%_exon Del_in_last_95%_exon Longest_continous_ins_position(lci) lci_length lci_%_range Longest_continous_del_position(lcd) lcd_length lcd_%_range' | column -t > gedit_exon_boundary_indels
#Indels at exon boundaries : format is : "number_of_indels(exon_number)"
#insertions at exon start for all exons
ais=`grep -v "Ins" gedit_exon_boundary_indels | awk '{print$2"("$1")"}' | grep -v "^-" | sed 's/exon_//g' | paste -s -d ","`
if [[ $ais == "" ]]; then ais="-"; else :; fi
#max value
aismx=`grep -v "Ins" gedit_exon_boundary_indels | awk '{print$2}' | sort -nr | head -1`
if [[ $aismx == "" ]]; then aismx="-" ; else :; fi
#insertions at exon end for all exons
aie=`grep -v "Ins" gedit_exon_boundary_indels | awk '{print$3"("$1")"}' | grep -v "^-" | sed 's/exon_//g' | paste -s -d ","`
if [[ $aie == "" ]]; then aie="-"; else :; fi
#max value
aiemx=`grep -v "Ins" gedit_exon_boundary_indels | awk '{print$3}' | sort -nr | head -1`
if [[ $aiemx == "" ]]; then aiemx="-" ; else :; fi
#deletions at exon start for all exons
ads=`grep -v "Ins" gedit_exon_boundary_indels | awk '{print$4"("$1")"}' | grep -v "^-" | sed 's/exon_//g' | paste -s -d ","`
if [[ $ads == "" ]]; then ads="-"; else :; fi
#max value
adsmx=`grep -v "Ins" gedit_exon_boundary_indels | awk '{print$4}' | sort -nr | head -1`
if [[ $adsmx == "" ]]; then adsmx="-" ; else :; fi
#deletions at exon end for all exons
ade=`grep -v "Ins" gedit_exon_boundary_indels | awk '{print$5"("$1")"}' | grep -v "^-" | sed 's/exon_//g' | paste -s -d ","`
if [[ $ade == "" ]]; then ade="-"; else :; fi
#max value
ademx=`grep -v "Ins" gedit_exon_boundary_indels | awk '{print$5}' | sort -nr | head -1`
if [[ $ademx == "" ]]; then ademx="-" ; else :; fi
#Longest continous insertion greater than 4 bases in the hit wrt to query
lci=`grep -v "Ins" gedit_exon_boundary_indels | awk '$7>4 {print$7}' | paste -s -d ","`
if [[ $lci == "" ]]; then lci="-" ; else :; fi
#Longest continous deletion greater than 4 bases in the hit wrt to query
lcd=`grep -v "Ins" gedit_exon_boundary_indels | awk '$10>4 {print$10}' | paste -s -d ","`
if [[ $lcd == "" ]]; then lcd="-" ; else :; fi
#START codon (upper case is necessary for later steps)
if [[ -s gblast_edited_consensus.fa ]]; then
stc=`grep -v ">" gblast_edited_consensus.fa | paste -s -d "" | cut -c -3 | tr '[:lower:]' '[:upper:]'`
else stc="-"; fi
#STOP codon
if [[ -s gblast_edited_consensus.fa ]]; then
enc=`grep -v ">" gblast_edited_consensus.fa | paste -s -d "" | rev | cut -c -3 | tr '[:lower:]' '[:upper:]' | rev`
else enc="-"; fi
#Splice site disruptions
s1=`awk '{print"exon_"++i,toupper($1)}' gedit_splice_sites | sed 1d | awk '$2!="AG" {print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print $0}'`
s2=`awk '{print"exon_"++i,toupper($2)}' gedit_splice_sites | sed '$d' | awk '$2!="GT"' | awk '$2!="GC"' | awk '{print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print $0}'`
#Protein consensus
if [[ -s gblast_edited_consensus_orf.fa ]]; then
pc=`grep -v ">" gblast_edited_consensus_orf.fa | paste -s -d "" | wc | awk '{print$3-$1}'`
else pc="-"; fi
# % of the gene in which the first Stop codons appears
sp=`tail -1 gblast_edited_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$pc" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
#Consensus nucleotide length
q=`grep -v ">" gblast_edited_consensus.fa | wc | awk '{print$3-$1}'`
#Check ORF presence absence
n=`grep -v ">" gblast_edited_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $n == "" ]] ; then o="❌";else o="✅"; fi
echo $j $gr $enh $ewah $qc $ge $gf $oab $ti $td $indf $ais $aismx $aie $aiemx $ads $adsmx $ade $ademx $lci $lcd $stc $enc $s1 $s2 $q $pc $sp $o
unset j gr enh ewah qry1 qry tl qc ge gf oab ti td indf ais aismx aie aiemx ads adsmx ade ademx lci lcd stc enc s1 s2 q pc sp n o
cd /home/neo/bird_db1/aswin/APOBEC1
#remove `sed 's/exon//g'` if you want exon name in the summary table (ususally not recommended beacuse adding "exon" make te table very wide for long genes )
done < <(grep -if gedit_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths) | sed 's/exon_//g' | sed '1i Species Group No_hits #_Almost_hits %_Q_cov Exon_gaps Total_gaps other_Anonyms Ins Del Frame_shifts Ins_@_exon_start[No(exon)](I1) I1_Max Ins_@_exon_end(I2) I2_Max Del_@_exon_start(D1) D1_Max Del_@_exon_end(D2) D2_Max Continous_ins_L Continous_del_L START STOP SS_1_disrupt(exon) SS_2_disrupt(exon) Cons_L Protein_con 1st_STOP_in_gene ORF' | column -t > gedit_principal_components

#(optional check points: If running manually)
#Possible losses due to partial exons
#awk '$5<90' gedit_principal_components	#query cover less than 90% at this point would suggest either partial exon loss or huge assembly error which is unlikely to occur

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.4. Add completely obtained ORFs to the total validated queries list, gene tree & add their sequences to the exon_wise folder

#List of queries that gave complete ORF with 100% query coverage & no gaps from gedit
awk '$3=="-" && $5=="100" && $6=="-" && $8=="-" && ($13<5 || $13=="-") && ($15<5 || $15=="-") && ($17<5 || $17=="-") && ($19<5 || $19=="-") && $24=="-" && $25=="-" && $NF=="✅" {print$1}' gedit_principal_components > gedit_validated_queries
cat gedit_validated_queries | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' >> total_validated_queries

#Add to gene tree
for i in `cat gedit_validated_queries`; do sed "/$i/ s/$/"gedit_"$i/g" -i gene_tree; unset j;done

#add sequence to folder
while read i
do
j=`echo $i | awk '{print$1}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
cat gblast_edited_consensus.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/"APOBEC1_"$j.fa
cat extracted_flanking_region.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/"APOBEC1_"$j"_flanking".fa
unset j
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if gedit_validated_queries /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.5. Save all exon pairwise alignment between query and subject from cgblast of all birds into a single file : later for quick reference

while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
echo ">"$j
for i in `ls pairwise_exon* | cut -f2,3 -d "_" | sed 's/.aln//g' | sort -V`; do echo "@"$i; grep "^extract" -B2 pairwise_"$i".aln | cut -c 22-; done
echo -e
unset j
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if gblast_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths) > pairwise_exon_alignments.aln

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.6. QC validated sequences

#Species with focal gene fragmented into different scaffolds or chromosomes
cut -f2- -d "_" total_validated_queries | sed 's/.fa//g' | xargs -n1 sh -c 'grep "^$0" gblast_principal_components' | column -t | awk '/ diff / {print$1}' > gene_breaking_assemblies

./qc_validated.sh validated_sequences validated_sequences_with_flanking_regions *_coding_exons_set APOBEC1 /home/neo/bird_db1/aswin/APOBEC1 gedit_validated_qc

#Find potential false positives in validated sequences by detecting outliers in nucleotide MSA
sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" gedit_validated_qc | awk -F ":" '/Outlier/ {print$2}' | sed 's/Nothing//g' | tr "," "\n" | awk NF | sed 's/^[ ]\+//g' > outliers

#Remove outliers from the validated sequence records
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers gedit_validated_queries
cat outliers | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' > outliers2
if [[ -s outliers2 ]]; then 
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers2 total_validated_queries
else :; fi
cat outliers | xargs -n1 sh -c 'sed "/$0/ s/[ ]\+[^ ]\+$//g" gene_tree -i'
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa) | xargs rm
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/*.fa) | xargs rm
rm outliers outliers2

##################################################################################################################################################################################################################################################################################################################
#7. Exonerate : To extract ORF of genes with exon boundaries (with high variability/ repeats/ false splice site ) incorrectly obtained from gedit

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#7.1. Create a list of subjects to run Exonerate : Run exonerate against extracted flanking region (not the whole genome)
cat gedit_principal_components | awk '$3=="-" && $6=="-" && $8=="-" && (($NF=="❌") || ($24!="-" && $NF=="✅") || ($25!="-" && $NF=="✅")) {print$1}' > exonerate_subjects

#Update suitable Query-Subject set
for i in `cat exonerate_subjects | cut -f1,2 -d "_"`; do j=`grep -if <(sed 's/.fa//g' total_validated_queries | cut -f2- -d "_") <(awk -v s="$i" '$2==s' /home/neo/bird_db1/aswin/tree/species_with_genome_data/all_possible_pairwise_distances | sort -k3 -n | awk '{print$1}') | head -1`; echo $j $i; unset j; done | column -t > exonerate_query_subject_set

#Compare query-subject sets
join -1 1 -2 1 <(awk '{print$2,$1}' gblast_query_subject_set | sort -k1) <(awk '{print$2,$1}' exonerate_query_subject_set | sort -k1) | cat -n | column -t | awk '{if($3==$4) print$0,"same";else print$0,"diff"}' | grep ".*diff\|$" --color=always | GREP_COLORS="mt=32" grep ".*same\|$" --color=always | sed '1i No Subject Gblast_query Exonerate_query Query_same/diff' | column -t > compare_gblast_exonerate_query_subject_sets

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#7.2. Run exonerate

while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
echo " ▶ "$j
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
#Protein as query (from nearest species)
#Find the most suitable(nearest) query species
q=`grep $j /home/neo/bird_db1/aswin/APOBEC1/exonerate_query_subject_set | awk '{print$1}'`
#Get the query sequence
pq=`find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1 -name "*$q*.fa" | awk -F "/" '{print$NF}'`
find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1 -name "*$q*.fa" -exec cp {} . \;
#Translate to amino acid ("*" is required to represent stop codon for exonerate)
transeq <(grep -v ">" $pq) -auto -stdout | sed "/>/ s/.*/>$pq/g" | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > exonerate_query.aa
#Combined extracted sequence as reference (optional: to get Splice sites & for QC)
grep -v ">" extracted_flanking_region.fa | paste -s -d "" | tr '[:lower:]' '[:upper:]' | sed '1i >extracted_flanking_region_combined' > extracted_flanking_region_combined.fa
#Run exonerate in exhaustive mode  with model protein2genome:bestfit & maximum intron allowed is 48bp
exonerate -E --model p2g:b exonerate_query.aa extracted_flanking_region_combined.fa --maxintron 48 --ryo "Rank: %r\nTarget in alignment: %ti : %tal (%tab - %tae)\nTarget in coding sequence: %ti : %tcl (%tcb - %tce)\nGene orientation: %g\nPercent identity: %pi\nPercent similarity: %ps\n>%ti (%tab - %tae)\n%tcs\n" --querytype protein --targettype dna --showtargetgff yes --fsmmemory 28000 --alignmentwidth 270 > Query_exon_combined_exonerate.out 2>/dev/null
#Create consensus sequence
#sed 1d Query_exon_combined_exonerate.out | awk '/Query range/; /^>/,/^$/' | grep -v "^>" | sed 's/.*Query range: />/g' | sed 's/ -> .*//g' | sed 's/^>/\x00&/' | sort -z -V | sort -z -V | tr -d '\0' | grep -v ">" | awk NF | paste -s -d "" | sed "1i >Query_exon_combined_exonerate" > Query_exon_combined_exonerate.out.fa
#Exons_wise sequence
bedtools getfasta -fi extracted_flanking_region_combined.fa -bed <(awk '$3=="exon" {print"extracted_flanking_region_combined",$4-1,$5,$3 | "sort -k2 -n"}' Query_exon_combined_exonerate.out | awk '{print$1,$2,$3,$4"_"++i}' OFS="\t") -name > Query_exon_combined_exonerate.out_exon_wise.fa 2>/dev/null
grep -v ">" Query_exon_combined_exonerate.out_exon_wise.fa | sed '1i >exonerate_consensus' > Query_exon_combined_exonerate.out.fa
#Check ORF
transeq <(grep -v ">" Query_exon_combined_exonerate.out.fa | paste -s -d "") Query_exon_combined_exonerate.out.fa_orf.fa 2>/dev/null
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" Query_exon_combined_exonerate.out.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" Query_exon_combined_exonerate.out.fa_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > Query_exon_combined_exonerate_translated
#For better visualization
gnaa=`head -1 Query_exon_combined_exonerate_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" Query_exon_combined_exonerate_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 Query_exon_combined_exonerate_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > Query_exon_combined_exonerate_nucleotide_amino_acid_view
unset j q pq gnaa
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if exonerate_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#7.3. Exonerate summary

while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
#Taxonomy order of the species
gr=`echo $j | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} /home/neo/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
#Query used
q=`grep $j /home/neo/bird_db1/aswin/APOBEC1/exonerate_query_subject_set | awk '{print$1}'`
pq=`find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1 -name "*$q*.fa" | awk -F "/" '{print$NF}'`
#Number of no hits
nnh=`expr $(grep ">" $pq -c) - $(grep ">" Query_exon_combined_exonerate.out_exon_wise.fa -c)`
splice_sites=`grep ">" $pq -c | awk '{print$0-1}'`
#splice site 1
a1=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$1=="GT" || $1=="GC"' | wc -l`
if [[ $a1 == $splice_sites ]]; then a2="✔"; else a2="x"; fi
#Exon with Splice site 1 disruptions
#in1=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | awk '{print++i,$0}' | awk '$2!="GT" {print$1}' | awk '{print$0+1}' | paste -s -d "," | awk '{if($1=="") print "-"; else print $0}'`
#splice site 2
a3=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$2=="AG"' | wc -l`
if [[ $a3 == $splice_sites ]]; then a4="✔"; else a4="x"; fi
#Exon with Splice site 2 disruptions
#in2=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | awk '{print++i,$0}' | awk '$3!="AG" {print$1}' | awk '{print$0+1}' | paste -s -d "," | awk '{if($1=="") print "-"; else print $0}'`
#START of the CDS
strt=`cat Query_exon_combined_exonerate.out.fa | grep -v ">" | paste -s -d "" | cut -c -3 | tr '[:lower:]' '[:upper:]'`
#% identity
b=`cat Query_exon_combined_exonerate.out | awk '/^Percent identity/{print$NF}' | paste -s -d ","`
if [[ $b == "" ]]; then b="-"; else :; fi
#% similarity
c=`cat Query_exon_combined_exonerate.out | awk '/^Percent similarity/{print$NF}' | paste -s -d ","`
if [[ $c == "" ]]; then c="-"; else :; fi
#insertions
d=`cat Query_exon_combined_exonerate.out | awk '{for(i=1;i<=NF;i++) if($i~/insertions/) print$(i+1)}' | awk '{sum+=$1;} END{print sum}' | awk '{print$1}'`
if [[ $d == "" ]]; then d="-"; else :; fi
#deletions
e=`cat Query_exon_combined_exonerate.out | awk '{for(i=1;i<=NF;i++) if($i~/deletions/) print$(i+1)}' | awk '{sum+=$1;} END{print sum}' | awk '{print$1*3}'`
if [[ $e == "" ]]; then e="-"; else :; fi
#Raw score
k=`cat Query_exon_combined_exonerate.out | awk -F ":" '/Raw score/ {print$2}' | tr -d " " | paste -s -d ","`
if [[ $k == "" ]]; then k="-"; else :; fi
#Query_length
ql=`grep -v ">" exonerate_query.aa | wc | awk '{print$3-$1}'`
#Query range
l=`cat Query_exon_combined_exonerate.out | awk -F ":" '/Query range/ {print$2}' | tr -d ">" | tr -d " " | paste -s -d ","`
if [[ $l == "" ]]; then l="-"; else :; fi
#Consensus nucleotide & protein length
m1=`grep -v ">" Query_exon_combined_exonerate.out.fa | wc | awk '{print$3-$1}'`
m2=`grep -v ">" Query_exon_combined_exonerate.out.fa_orf.fa | wc | awk '{print$3-$1}'`
# % of the gene in which the first Stop codons appears
sp=`tail -1 Query_exon_combined_exonerate_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$m2" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
#Query % covered {(Target hits - insertion + deletion - frameshifts) / query * 100}
ins=`awk '{for(i=1;i<=NF;i++) if($i~/insertions/) print$(i+1)}' Query_exon_combined_exonerate.out | awk '{s+=$1} END{print s}' | awk '{print$1}'`
del=`awk '{for(i=1;i<=NF;i++) if($i~/deletions/) print$(i+1)}' Query_exon_combined_exonerate.out | awk '{s+=$1} END{print s}' | awk '{print$1*3}'`
fs=`awk '{for(i=1;i<=NF;i++) if($i~/frameshifts/) print$(i+1)}' Query_exon_combined_exonerate.out | awk '{s+=$1} END{print s}'`
if [[ $fs == "" ]]; then fs="-"; else :; fi
tl=`grep -v ">" Query_exon_combined_exonerate.out.fa | wc | awk '{print$3-$1}'`
qnl=`expr $ql \* 3`
if [[ $fs == "-" ]]; then
qc=`calc $(calc $tl - $ins + $del) / $qnl \* 100 | tr -d "~" | awk '{$1+=0}1' CONVFMT="%.2f"`
else
qc=`calc $(calc $tl - $ins + $del - $fs) / $qnl \* 100 | tr -d "~" | awk '{$1+=0}1' CONVFMT="%.2f"`
fi
#Check ORF presence absence
n=`grep -v ">" Query_exon_combined_exonerate.out.fa_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $n == "" ]] ; then o="❌";else o="✅"; fi
echo $j $gr $nnh $a2 $a4 $strt $b $c $d $e $fs $k $ql $l $m1 $m2 $qc $sp $o
unset j gr q pq nnh splice_sites a1 a2 a3 a4 strt b c d e fs k ql l m1 m2 ins del tl qnl qc sp n o
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if exonerate_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths) | sed '1i Species Group #_no_hits SS_1 SS_2 START %_id %_sim Insertions Deletions Frameshifts Score Query_L Query_range Consensus Protein_cons %_Query_cov 1st_STOP_in_gene ORF' | column -t > exonerate_summary

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#7.4. Add completely obtained ORFs to the total validated queries list, gene tree & add their sequences to the exon_wise folder

#List of queries that gave complete ORF with no splice site discruption & 100% query coverage from exonerate
grep -v "^Species" exonerate_summary | awk '$4=="✔" && $5=="✔" && $17~"100" && $NF=="✅" {print$1}' > exonerate_validated_queries
cat exonerate_validated_queries | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' >> total_validated_queries

#Add to gene tree
for i in `cat exonerate_validated_queries`; do sed "/$i/ s/$/"exonerate_"$i/g" -i gene_tree; unset j;done

#add sequence to folder : exonerate consensus sequences are not exon-wise
while read i
do
j=`echo $i | awk '{print$1}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
cat Query_exon_combined_exonerate.out_exon_wise.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/"APOBEC1_"$j.fa
cat extracted_flanking_region.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/"APOBEC1_"$j"_flanking".fa
unset j
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if exonerate_validated_queries /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#7.5. QC validated sequences

#Species with focal gene fragmented into different scaffolds or chromosomes
cut -f2- -d "_" total_validated_queries | sed 's/.fa//g' | xargs -n1 sh -c 'grep "^$0" gblast_principal_components' | column -t | awk '/ diff / {print$1}' > gene_breaking_assemblies

./qc_validated.sh validated_sequences validated_sequences_with_flanking_regions *_coding_exons_set APOBEC1 /home/neo/bird_db1/aswin/APOBEC1 exonerate_validated_qc

#Find potential false positives in validated sequences by detecting outliers in nucleotide MSA
sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" exonerate_validated_qc | awk -F ":" '/Outlier/ {print$2}' | sed 's/Nothing//g' | tr "," "\n" | awk NF | sed 's/^[ ]\+//g' > outliers

#Remove outliers from the validated sequence records
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers exonerate_validated_queries
cat outliers | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' > outliers2
if [[ -s outliers2 ]]; then 
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers2 total_validated_queries
else :; fi
cat outliers | xargs -n1 sh -c 'sed "/$0/ s/[ ]\+[^ ]\+$//g" gene_tree -i'
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa) | xargs rm
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/*.fa) | xargs rm
rm outliers outliers2

##################################################################################################################################################################################################################################################################################################################
#8. Exoneredit : Make a consensus sequence from pairwise alignment between Exonerate query & exonerate consensus

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#8.1. Subjects for exoneredit : species that didn't form complete ORF from exonerate
grep -if exonerate_validated_queries exonerate_subjects -v > exoneredit_subjects

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#8.2. Run exoneredit : extracting final consensus based on pairwise alignment between exonerate query & exonerate target

while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
echo " ▶ "$j
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
#Find the most suitable(nearest) query species
q=`grep $j /home/neo/bird_db1/aswin/APOBEC1/exonerate_query_subject_set | awk '{print$1}'`
#Get the query sequence
pq=`find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1 -name "*$q*.fa" | awk -F "/" '{print$NF}'`
#pairwise align consensus sequences
needle -asequence <(grep -v ">" $pq | sed '1i >Query') Query_exon_combined_exonerate.out.fa -auto -stdout -awidth 280 > pairwise_query_exoneredit.aln
#find START site
start=`needle -asequence <(grep -v ">" $pq | sed '1i >Query') <(grep -v ">" Query_exon_combined_exonerate.out.fa | sed '1i >Exonerate') -auto -stdout -awidth 100000 | grep "^Query" -A2 | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":" | awk '{print$0+1}'`
#Cut everything before START site
needle -asequence <(grep -v ">" $pq | sed '1i >Query') <(grep -v ">" Query_exon_combined_exonerate.out.fa | sed '1i >Exonerate') -auto -stdout -awidth 100000 | grep "^Query" -A2 | grep "|" -C1 | cut -c $start- > tmp1
#Find STOP site
end=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+1}'`
#Calculate some basic stats about exoneredit before making consensus sequence
ln=`grep -v ">" $pq  | wc | awk '{print$3-$1}'`
#cat tmp1 | cut -c -$end | head -1 | tr -d "-" | wc | awk '{print$3-$1}' | awk -v p="$ln" '{print p,$1,($1/p)*100}' | sed "s/^/exoneredit /g" > exoneredit_query_covered
#insertions
cat tmp1 | cut -c -$end | head -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exoneredit /g" > exoneredit_insertions
#deletions
cat tmp1 | cut -c -$end | tail -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exoneredit /g" > exoneredit_deletions
#Create consensus : Delete everything after STOP site
cat <(echo ">exoneredit") <(cat tmp1 | cut -c -$end | tail -1) | tr -d "-" > exoneredit_consensus.fa
#Extract exon_wise sequence
blat extracted_flanking_region_combined.fa exoneredit_consensus.fa -out=pslx exoneredit.pslx > /dev/null 2>&1
egrep -v "psLayout|match|^---" exoneredit.pslx | awk NF | awk '{print$NF}' | tr "," "\n" | tr '[:lower:]' '[:upper:]' | awk NF | awk '{print ">exon_"++i"\n"$0}' > exoneredit_consensus_exon_wise.fa
rm tmp1
#Generate ORF
transeq <(grep -v ">" exoneredit_consensus.fa | paste -s -d "" | sed "1i \>$j") exoneredit_consensus_orf.fa 2>/dev/null
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" exoneredit_consensus.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" exoneredit_consensus_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > exoneredit_translated
#For better visualization of ORF
gnaa=`head -1 exoneredit_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" exoneredit_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 exoneredit_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > exoneredit_nucleotide_amino_acid_view
unset j q pq start ln p end gnaa
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if exoneredit_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#8.3. Exoneredit summary

while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
#Taxonomy order of the species
gr=`echo $j | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} /home/neo/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
#Query used
q=`grep $j /home/neo/bird_db1/aswin/APOBEC1/exonerate_query_subject_set | awk '{print$1}'`
pq=`find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1 -name "*$q*.fa" | awk -F "/" '{print$NF}'`
splice_sites=`grep ">" $pq -c | awk '{print$0-1}'`
#splice site 1
a1=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] |  awk '$1=="GT" || $1=="GC"' | wc -l`
if [[ $a1 == $splice_sites ]]; then a2="✔"; else a2="x"; fi
#splice site 2
a3=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$2=="AG"' | wc -l`
if [[ $a3 == $splice_sites ]]; then a4="✔"; else a4="x"; fi
#Similarity
k=`awk -F "(" '/Similarity/ {print$NF}' pairwise_query_exoneredit.aln | tr -d ")"`
if [[ $k == "" ]]; then k="-"; else :; fi
#insertions
m=`awk '{print$2}' exoneredit_insertions`
#Deletions
n=`awk '{print$2}' exoneredit_deletions`
#Query covered { (Hit length - ins + del) - query length)
ins=`grep "^Query" -w pairwise_query_exoneredit.aln | tr -c -d "-" | wc | awk '{print$3-$1}'`
del=`grep "^Query_exon" pairwise_query_exoneredit.aln | tr -c -d "-" | wc | awk '{print$3-$1}'`
tl=`grep -v ">" Query_exon_combined_exonerate.out.fa | wc | awk '{print$3-$1}'`
ql=`grep -v ">" exonerate_query.aa | wc | awk '{print($3-$1)*3}'`
o=`calc $(calc $tl - $ins + $del) / $ql \* 100 | tr -d "~" | awk '{$1+=0}1' CONVFMT="%.2f"`
#Consensus sequence length
p=`grep -v ">" exoneredit_consensus.fa  | wc | awk '{print$3-$1}'`
#Protein length
unset q
q=`grep -v ">" exoneredit_consensus_orf.fa  | wc | awk '{print$3-$1}'`
#START codon
r=`grep -v ">" exoneredit_consensus.fa  | cut -c -3 | tr '[:lower:]' '[:upper:]'`
# % of the gene in which the first Stop codons appears
s=`tail -1 exoneredit_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$q" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.1f"`
if [[ $s == "" ]]; then s="-"; else :; fi
#Check ORF presence absence
t=`grep -v ">" exoneredit_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $t == "" ]] ; then t1="❌";else t1="✅"; fi
echo $j $gr $a2 $a4 $k $m $n $o $p $q $r $s $t1
unset j gr q pq splice_sites a1 a2 a3 a4 k m n ins del tl ql o p q r s z y t t1
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if exoneredit_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths) | sed '1i Species Group SS_1 SS_2 %_Similarity Ins Del %_Query_cov Consensus_L Protein_L START 1st_STOP_in_gene ORF' | column -t > exoneredit_summary

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#8.4. Add completely obtained ORFs to the total validated queries list, gene tree & add their sequences to the exon_wise folder

#List of queries that gave complete ORF with 100% query coverage from exoneredit

#If START/STOP is within the query sequence w.r.t the subject then 100% query coverage criteria can be removed
#awk '$3=="✔" && $4=="✔" && $8~"100" && $NF=="✅" {print$1}' exoneredit_summary > exoneredit_validated_queries
awk '$3=="✔" && $4=="✔" && $8>97 && $NF=="✅" {print$1}' exoneredit_summary > exoneredit_validated_queries
cat exoneredit_validated_queries | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' >> total_validated_queries

#Add to gene tree
for i in `cat exoneredit_validated_queries`; do sed "/$i/ s/$/"exoneredit_"$i/g" -i gene_tree; unset j;done

#add sequence to folder : exoneredit consensus sequences are not exon-wise
while read i
do
j=`echo $i | awk '{print$1}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
cat exoneredit_consensus_exon_wise.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/"APOBEC1_"$j.fa
cat extracted_flanking_region.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/"APOBEC1_"$j"_flanking".fa
unset j
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if exoneredit_validated_queries /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#8.5. QC validated sequences

#Species with focal gene fragmented into different scaffolds or chromosomes
cut -f2- -d "_" total_validated_queries | sed 's/.fa//g' | xargs -n1 sh -c 'grep "^$0" gblast_principal_components' | column -t | awk '/ diff / {print$1}' > gene_breaking_assemblies

./qc_validated.sh validated_sequences validated_sequences_with_flanking_regions *_coding_exons_set APOBEC1 /home/neo/bird_db1/aswin/APOBEC1 exoneredit_validated_qc

#Find potential false positives in validated sequences by detecting outliers in nucleotide MSA
sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" exoneredit_validated_qc | awk -F ":" '/Outlier/ {print$2}' | sed 's/Nothing//g' | tr "," "\n" | awk NF | sed 's/^[ ]\+//g' > outliers

#Remove outliers from the validated sequence records
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers exoneredit_validated_queries
cat outliers | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' > outliers2
if [[ -s outliers2 ]]; then 
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers2 total_validated_queries
else :; fi
cat outliers | xargs -n1 sh -c 'sed "/$0/ s/[ ]\+[^ ]\+$//g" gene_tree -i'
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa) | xargs rm
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/*.fa) | xargs rm
rm outliers outliers2

##################################################################################################################################################################################################################################################################################################################
#9. Spaln

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#9.1. Filter out subjects that didn't form ORF from exoneredit
grep -if exoneredit_validated_queries exoneredit_subjects -v > spaln_subjects

#Update suitable Query-Subject set
for i in `cat spaln_subjects | cut -f1,2 -d "_"`; do j=`grep -if <(sed 's/.fa//g' total_validated_queries | cut -f2- -d "_") <(awk -v s="$i" '$2==s' /home/neo/bird_db1/aswin/tree/species_with_genome_data/all_possible_pairwise_distances | sort -k3 -n | awk '{print$1}') | head -1`; echo $j $i; unset j; done | column -t > spaln_query_subject_set

#Compare query-subject sets
join -1 1 -2 1 <(awk '{print$2,$1}' exonerate_query_subject_set | sort -k1) <(awk '{print$2,$1}' spaln_query_subject_set | sort -k1) | cat -n | awk '{if($3==$4) print$0,"same";else print$0,"diff"}' | column -t | grep ".*diff\|$" --color=always | GREP_COLORS="mt=32" grep ".*same\|$" --color=always > compare_gblast_2nd_gblast_tmp_query_subject_sets

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#9.2. Run spaln

while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
echo "  ▶ " $j
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1
#Query based on phylogeny
q=`grep $j /home/neo/bird_db1/aswin/APOBEC1/spaln_query_subject_set | awk '{print$1}'`
pq=`find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1 -name "*$q*.fa" | awk -F "/" '{print$NF}'`
#Get the query sequence
find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1 -name "*$q*.fa" -exec cp {} . \;
#Translate to amino acid ("*" is required to represent stop codon for Spaln)
transeq <(grep -v ">" $pq) -auto -stdout | sed "/>/ s/.*/>$pq/g" | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > spaln_query.aa
#Query exon count
pqe1=`grep ">" $pq -c`
#Query cds length
pqe2=`grep -v ">" $pq | wc | awk '{print$3-$1}'`
#Splice sites expected in target
spn=`echo $pqe1 | awk '{print$1-1}'`
#Run Spaln with different combinations of -yX(intr/cross species parameters) -yS (species-specific parameter)
for x in 0 1
do
for y in 0 50 100
do
echo -n > spaln_out_tmp
echo -n > spaln_tmp
spaln -M1 -l280 -O4 -ya0 -yX$x -yS$y extracted_flanking_region_combined.fa spaln_query.aa > spaln_out_tmp 2>/dev/null
bedtools getfasta -fi extracted_flanking_region_combined.fa -bed <(awk '!/^#|^@/{print"extracted_flanking_region_combined",$9-1,$10,"exon_"++i}' OFS="\t" spaln_out_tmp) -name > spaln_tmp 2>/dev/null
x1=`grep ">" -c spaln_tmp`
x2=`grep -v ">" spaln_tmp | wc | awk '{print$3-$1}'`
x3=`needle <(grep -v ">" spaln_tmp) <(grep -v ">" $pq) -auto -stdout -awidth 280 | awk '/Identity|Similarity|Gaps|Score/ {print$NF}' | tr -d ')(%' | paste -s -d " "` 2>/dev/null
transeq <(grep -v ">" spaln_tmp | paste -s -d "" | sed "1i \>$j") spaln_tmp_orf.fa 2>/dev/null
#Splice sites
sp1=`awk '!/^#|^@/ {print$NF}' spaln_out_tmp | tr "." " " | awk NF | tr '[:lower:]' '[:upper:]' | awk '$1=="GT" || $1=="GC"' | wc -l`
if [[ $sp1 == $spn ]]; then sp1s="✔"; else sp1s="x"; fi
sp2=`awk '!/^#|^@/ {print$NF}' spaln_out_tmp | tr "." " " | awk NF | tr '[:lower:]' '[:upper:]' | awk '$2=="AG"' | wc -l`
if [[ $sp2 == $spn ]]; then sp2s="✔"; else sp2s="x"; fi
#Check ORF presence absence
x4=`grep -v ">" spaln_tmp_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $x4 == "" ]] ; then x5="❌";else x5="✅"; fi
echo $x $y $pqe1 $x1 $pqe2 $x2 $x3 $sp1 $sp1s $sp2 $sp2s $x5
unset x1 x2 x3 x4 sp1 sp1s sp2 sp2s x5
rm spaln_out_tmp spaln_tmp spaln_tmp_orf.fa
done
done | sed '1i yX yS Query_exons Spaln_exons Query_length Spaln_length %_Identity %_Similarity %_Gaps Score No_Splice_site1 SS1_disruption No_Splice_site2 SS2_disruption ORF' | column -t > spaln_summary
unset x y

#Choose best parameters based on exon count match & alignment stats
spaln_params=`cat <(awk '$3==$4 && $12=="✔" && $14=="✔"' spaln_summary | grep "✅") <(awk '$3==$4' spaln_summary | grep "✅") <(awk '$3==$4' spaln_summary | grep "❌") | awk '!a[$0]++' | awk 'NR==1{print$1,$2}'`
if [[ $spaln_params == "" ]]
then
unset spaln_params
spaln_params=`grep -v "yX" spaln_summary | sort -k4,4 -nr | awk 'NR==1{print$1,$2}'`
echo "     ▬  No ideal Spaln parameters found! Taking the parameters that gave the highest number of exon hits : " $spaln_params
else
echo "     ▬  Spaln parameters used : yX & yS : "$spaln_params
fi
while read sp
do
x=`echo $sp | awk '{print$1}'`
y=`echo $sp | awk '{print$2}'`
#Final spaln output in multiple formats : each with different info
spaln -M1 -l280 -O4 -ya0 -yX$x -yS$y extracted_flanking_region_combined.fa spaln_query.aa > spaln_out 2>/dev/null
spaln -M1 -l280 -O1 -ya0 -yX$x -yS$y extracted_flanking_region_combined.fa spaln_query.aa > spaln_aln_out 2>/dev/null
spaln -M1 -l280 -O2 -ya0 -yX$x -yS$y extracted_flanking_region_combined.fa spaln_query.aa > spaln_gff3_out 2>/dev/null
bedtools getfasta -fi extracted_flanking_region_combined.fa -bed <(awk '!/^#|^@/{print"extracted_flanking_region_combined",$9-1,$10,"exon_"++i}' OFS="\t" spaln_out) -name > spaln_exon_wise.fa 2>/dev/null
awk '!/^#|^@/ {print$NF}' spaln_out | tr "." " " | awk NF | tr '[:lower:]' '[:upper:]' > spaln_splice_sites
done < <(echo $spaln_params)

#Generate ORF
transeq <(grep -v ">" spaln_exon_wise.fa | paste -s -d "" | sed "1i \>$j") spaln_exon_wise_orf.fa 2>/dev/null
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" spaln_exon_wise.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" spaln_exon_wise_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > spaln_translated
#For better visualization of ORF
gnaa=`head -1 spaln_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" spaln_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 spaln_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > spaln_nucleotide_amino_acid_view
#Pairwise exon-wise algnment between query & final consensus
exalign $pq spaln_exon_wise.fa > spaln_query_hit_pairwise_exon.aln 2>/dev/null
#QC spaln output : Exon boundar indels
#Exon boundary insertions
while mapfile -t -n 4 ary && ((${#ary[@]}))
do
ex=`printf '%s\n' "${ary[@]}" | head -1 | tr -d ">"`
printf '%s\n' "${ary[@]}" > tmp
cat tmp | grep "spaln_exon" | grep [0-9] -aob | head -1 | cut -f1 -d ":" | awk '{print$1+3}' | xargs -n1 sh -c 'cat tmp | cut -c $0-' 2>/dev/null | awk NF | sed 's/ \+[0-9]\+//g' > tmp1
#insertions at exon start : counted if a dash is observed within the query sequence anywhere within the first 5% of pairwise alignment between query & subject
is=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/APOBEC1/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $is == "" ]]; then is="-"; else :; fi
#insertions at exon end : counted if a dash is observed within the query sequence anywhere after 95% of pairwise alignment between query & subject
ie=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/APOBEC1/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ie == "" ]]; then ie="-"; else :; fi
#deletions at exon start : counted if a dash is observed within the subject sequence anywhere within the first 10% (since subject sequence is longer due to flanking regions) of pairwise alignment between query & subject
ds=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/spaln_exon/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ds == "" ]]; then ds="-"; else :; fi
#deletions at exon end : counted if a dash is observed within the subject sequence anywhere after 90% of pairwise alignment between query & subject
de=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/spaln_exon/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
#Longest continous insertions : Position, Length & % of gene in which it appears
if [[ $(head -1 tmp1 | grep "\-") == "" ]]; then lci="- - -"; else
lci=`bedtools merge -i <(cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/APOBEC1/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
fi
#Longest continous deletions : Position, Length & % of gene in which it appears
if [[ $(tail -1 tmp1 | grep "\-") == "" ]]; then lcd="- - -"; else
lcd=`bedtools merge -i <(cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/spaln_exon/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
#lci=`cat tmp2 | grep "\-\+" -o | tr -d " " | awk '{print length}' | sort -nr | head -1`
fi
if [[ $de == "" ]]; then de="-"; else :; fi
echo $ex $is $ie $ds $de $lci $lcd
unset ex is ie ds de lci lcd
rm tmp tmp1
done < <(cat spaln_query_hit_pairwise_exon.aln | sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" | awk NF) | sed '1i Exons Ins_in_1st_5%_exon Ins_in_last_95%_exon Del_in_1st_5%_exon Del_in_last_95%_exon Longest_continous_ins_position(lci) lci_length lci_%_range Longest_continous_del_position(lcd) lcd_length lcd_%_range' | column -t > spaln_boundary_indels_wrt_query
unset j q pq pqe1 pqe2 spn sp x y spaln_params gnaa r
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if spaln_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#9.3. Spaln summary

while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
#Taxonomy order of the species
gr=`echo $j | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} /home/neo/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
#Query based on phylogeny
q=`grep $j /home/neo/bird_db1/aswin/APOBEC1/spaln_query_subject_set | awk '{print$1}'`
pq=`find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1 -name "*$q*.fa" | awk -F "/" '{print$NF}'`
#Number of exons
qen=`grep ">" $pq -c`
ten=`grep ">" spaln_exon_wise.fa -c`
#No of exons with No hits from spaln
nh=`expr $qen - $ten`
#Splice sites
splice_sites=`grep ">" $pq -c | awk '{print$0-1}'`
#splice site 1
a1=`awk '$1=="GT" || $1=="GC"' spaln_splice_sites | wc -l`
if [[ $a1 == $splice_sites ]]; then a2="✔"; else a2="x"; fi
#splice site 2
a3=`awk '$2=="AG"' spaln_splice_sites | wc -l`
if [[ $a3 == $splice_sites ]]; then a4="✔"; else a4="x"; fi
#Splice variants - splice sites with "GC"
ssv=`cat spaln_splice_sites -n | awk '$2=="GC" {print$1}' | paste -s -d ","`
if [[ $ssv == "" ]]; then ssv="-"; else :; fi
#Exons with boundary indels
bi1=`grep -v "In" spaln_boundary_indels_wrt_query | awk '$2>0 {print$2,$1}' | sed 's/exon_//g' | sort -k1,1 -nr | awk '{print$1"("$2")"}' | paste -s -d ","`
if [[ $bi1 == "" ]]; then bi1="-"; else :; fi
bi2=`grep -v "In" spaln_boundary_indels_wrt_query | awk '$3>0 {print$3,$1}' | sed 's/exon_//g' | sort -k1,1 -nr | awk '{print$1"("$2")"}' | paste -s -d ","`
if [[ $bi2 == "" ]]; then bi2="-"; else :; fi
bi3=`grep -v "In" spaln_boundary_indels_wrt_query | awk '$4>0 {print$4,$1}' | sed 's/exon_//g' | sort -k1,1 -nr | awk '{print$1"("$2")"}' | paste -s -d ","`
if [[ $bi3 == "" ]]; then bi3="-"; else :; fi
bi4=`grep -v "In" spaln_boundary_indels_wrt_query | awk '$5>0 {print$5,$1}' | sed 's/exon_//g' | sort -k1,1 -nr | awk '{print$1"("$2")"}' | paste -s -d ","`
if [[ $bi4 == "" ]]; then bi4="-"; else :; fi
#Number of amino acid insertions spaln reported
ewi1=`egrep -v "^#|^@" spaln_gff3_out | awk '{print NR";"$0}' | awk -F ";" '{for(z=1;z<=NF;z++) if($z~/Gap.*I/) print$NF}' | sed 's/Gap=//g' | grep "I[0-9]\+" -o | sed 's/I//g' | awk '{a+=$1} END{print a}'`
if [[ $ewi1 == "" ]]; then ewi1="-"; else :; fi
#Number of amino acid deletions spaln reported including stop codon in last exon
ewi2=`egrep -v "^#|^@" spaln_gff3_out | awk '{print NR";"$0}' | awk -F ";" '{for(z=1;z<=NF;z++) if($z~/Gap.*D/) print$NF}' | sed 's/Gap=//g' | grep "D[0-9]\+" -o | sed 's/D//g' | awk '{a+=$1} END{print a}'`
if [[ $ewi2 == "" ]]; then ewi2="-"; else :; fi
#Exon boundary indels : max indel boundary observed at exon start
e1bim=`grep -v "In" spaln_boundary_indels_wrt_query | awk '$7>0 {print$1,$7}' | sort -k2,2 -nr | sed 's/exon_//g' | awk 'NR==1{print$2"("$1")"}'`
if [[ $e1bim == "" ]]; then e1bim="-"; else :; fi
#Exon boundary indels : max indel boundary observed at exon end
e2bim=`grep -v "In" spaln_boundary_indels_wrt_query | awk '$10>0 {print$1,$10}' | sort -k2,2 -nr | sed 's/exon_//g' | awk 'NR==1{print$2"("$1")"}'`
if [[ $e2bim == "" ]]; then e2bim="-"; else :; fi
#START of the CDS
strt=`grep -v ">" spaln_exon_wise.fa | paste -s -d "" | cut -c -3 | tr '[:lower:]' '[:upper:]'`
#Basic stats
bs=`needle <(grep -v ">" spaln_exon_wise.fa | sed '1i >spaln') <(grep -v ">" $pq | sed '1i >query') -auto -stdout -awidth 280 | awk '/Identity|Similarity|Score/ {print$NF}' | tr -d ')(%' | paste -s -d " "`
#insertions
ins=`needle <(grep -v ">" spaln_exon_wise.fa | sed '1i >spaln') <(grep -v ">" $pq | sed '1i >query') -auto -stdout -awidth 280 | grep "^query" | awk '{print$3}' | tr -cd "-" | wc | awk '{print$3-$1}'`
#Deletions
del=`needle <(grep -v ">" spaln_exon_wise.fa | sed '1i >spaln') <(grep -v ">" $pq | sed '1i >query') -auto -stdout -awidth 280 | grep "^spaln" | awk '{print$3}' | tr -cd "-" | wc | awk '{print$3-$1}'`
#Query length
ql=`grep -v ">" spaln_query.aa | wc | awk '{print$3-$1}'`
#Protein Consensus length
pl=`grep -v ">" spaln_exon_wise_orf.fa | wc | awk '{print$3-$1}'`
#Query range
qr=`cat spaln_aln_out | grep "^>extracted_flanking_region_combined" | grep -oP '\(\K[^\)]+' | tr -d " " | awk 'NR==2'`
#Target range
tr=`cat spaln_aln_out | grep "^>extracted_flanking_region_combined" | grep -oP '\(\K[^\)]+' | tr -d " " | awk 'NR==1'`
#DNA consensus length
dl=`grep -v ">" $pq | wc | awk '{print$3-$1}'`
# % Query covered
pqc=`echo $qr | awk -v l="$ql" -F "-" '{print($2-$1+1)/l*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
# % of the gene in which the first Stop codons appears
s=`tail -1 spaln_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$ql" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.1f"`
if [[ $s == "" ]]; then s="-"; else :; fi
#Check ORF presence absence
t=`grep -v ">" spaln_exon_wise_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $t == "" ]] ; then t1="❌";else t1="✅"; fi
echo $j $gr $qen $ten $nh $a2 $a4 $ssv $bi1 $bi2 $bi3 $bi4 $ewi1 $ewi2 $e1bim $e2bim $strt $bs $ins $del $ql $pl $qr $tr $dl $pqc $s $t1
unset j gr q pq qen ten nh splice_sites a1 a2 a3 a4 ssv bi1 bi2 bi3 bi4 ewi1 ewi2 e1bim e2bim strt bs insdel ql pl qr tr dl l pqc s y t t1
done < <(grep -if spaln_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths) | sed '1i Species Group Qry_exons Trgt_exons #_No_hits SS_1 SS_2 SS_var Exon_strt_ins Exon_end_ins Exon_strt_del Exon_end_del aa_ins aa_del lci_length lcd_length START %_id %_sim Score Ins Del Query_L Hit_L Query_range Target_range Consensus %_Query_cov 1st_STOP_in_gene ORF' | column -t > spaln_summary

#Spaln summary have some columns with indel info have exon info also in the brackets; this exon info is not required when fiding validated queries, hence another table is made removing exon info
#In case if some columns have multiple values sepered by "," then the maximum values is chosen by taking the 1st value (here the spaln summary itself is made in such a way that mulitple values are arranged in descending order)
while read i
do 
echo $i | tr " " "\n" | sed 's/([0-9]\+)//g' | cut -f1 -d "," | paste -s -d " "
done < spaln_summary | column -t > spaln_summary_for_validation

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#9.4. Add completely obtained ORFs to the total validated queries list, gene tree & add their sequences to the exon_wise folder

#List of queries that gave complete ORF with 100% query coverage from spaln (boundary indels max must be less than 7 means a hit in target extrapolated with 2 extra codons (6 bases) to make full ORF can be allowed but not more than 6 bases should be added to make a complete ORF, that would be a false positive)

#This criteria can be customized such as highest allowed exon boundary indel (here : less than 7) & longest allowed continous indel (here : less than 10)
awk '$3==$4 && $5=="0" && $6=="✔" && $7=="✔" && $9<7 && $10<7 && $11<7 && $12<7 && $15<10 && $16<10 && $28~"100" && $NF=="✅" {print$1}' spaln_summary_for_validation > spaln_validated_queries
cat spaln_validated_queries | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' >> total_validated_queries

#Add to gene tree
for i in `cat spaln_validated_queries`; do sed "/$i/ s/$/"spaln_"$i/g" -i gene_tree; unset j;done

#add sequence to folder
while read i
do
j=`echo $i | awk '{print$1}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
cat spaln_exon_wise.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/"APOBEC1_"$j.fa
cat extracted_flanking_region.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/"APOBEC1_"$j"_flanking".fa
unset j
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if spaln_validated_queries /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#9.5. QC validated sequences

#Species with focal gene fragmented into different scaffolds or chromosomes
cut -f2- -d "_" total_validated_queries | sed 's/.fa//g' | xargs -n1 sh -c 'grep "^$0" gblast_principal_components' | column -t | awk '/ diff / {print$1}' > gene_breaking_assemblies

./qc_validated.sh validated_sequences validated_sequences_with_flanking_regions *_coding_exons_set APOBEC1 /home/neo/bird_db1/aswin/APOBEC1 spaln_validated_qc

#Find potential false positives in validated sequences by detecting outliers in nucleotide MSA
sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" spaln_validated_qc | awk -F ":" '/Outlier/ {print$2}' | sed 's/Nothing//g' | tr "," "\n" | awk NF | sed 's/^[ ]\+//g' > outliers

#Remove outliers from the validated sequence records
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers spaln_validated_queries
cat outliers | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' > outliers2
if [[ -s outliers2 ]]; then 
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers2 total_validated_queries
else :; fi
cat outliers | xargs -n1 sh -c 'sed "/$0/ s/[ ]\+[^ ]\+$//g" gene_tree -i'
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa) | xargs rm
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/*.fa) | xargs rm
rm outliers outliers2

##################################################################################################################################################################################################################################################################################################################
#10. Second Genome blast & gedit on species with updated closest query

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#10.1. Filter out subjects that didn't form ORF

#grep -if <(cat total_validated_queries | cut -f2,3- -d "_" | sed 's/.fa//g') <(grep -if gedit_validated_queries gedit_subjects -v) -v > 2nd_gblast_subjects_tmp
grep -if <(cut -f2,3 -d "_" total_validated_queries|sed 's/.fa//g') gblast_subjects -v > 2nd_gblast_subjects_tmp

#Find the closest query which gives best hits for 2nd gblast subjects

#Create all possbible pairwise distances between all species pairs - Create once & keep it in a path don't need eveytime you run the pipeline
while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
#lowest distance : all species that are closest to the subject 
ld=`grep -if <(sed 's/.fa//g' total_validated_queries | cut -f2- -d "_") <(awk -v s="$j" '$2==s' /home/neo/bird_db1/aswin/tree/species_with_genome_data/all_possible_pairwise_distances | sort -k3 -n) | awk 'NR==1{print$3}'`
grep -if <(sed 's/.fa//g' total_validated_queries | cut -f2- -d "_") <(awk -v s="$j" '$2==s' /home/neo/bird_db1/aswin/tree/species_with_genome_data/all_possible_pairwise_distances | sort -k3 -n) | grep "$ld" | awk '{print$1}' > closest_query_list
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
#Loop through all the closest queries (including flanking region) & run blastn to find which gives the best hits in the subject
for x in `grep -if /home/neo/bird_db1/aswin/APOBEC1/closest_query_list <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/*.fa)`
do
blastn -task blastn -evalue 0.01 -db ../../genome/GC*.fna -query $x -num_threads 4 -outfmt "6 qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand qseq sseq"| sed '1i Query\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\tQuery_sequence\tSubject_sequence\n' > 2nd_test.out
gblast_best_hits $x 2nd_test.out 2nd_best_hits
q=`echo $x | awk -F "/" '{print$NF}'`
a=`awk 'NR>1{print$1}' 2nd_best_hits | wc -l`
b=`awk 'NR>1{print$2}' 2nd_best_hits | sort -u | wc -l`
c=`awk 'NR>1{print$3}' 2nd_best_hits | awk '{s+=$1}END{print s}'`
d=`awk 'NR>1{print$4}' 2nd_best_hits | awk '{s+=$1}END{print s}'`
e=`awk 'NR>1{print$10}' 2nd_best_hits | awk '{s+=$1}END{print s}'`
f=`awk 'NR>1{print$13}' 2nd_best_hits | awk '{s+=$1}END{print s}'`
g=`awk 'NR>1{print$14}' 2nd_best_hits | awk '{s+=$1}END{print s}'`
h=`awk 'NR>1{print$NF}' 2nd_best_hits | sort -u | wc -l`
echo $q $a $b $c $d $e $f $g $h
unset q a b c d e f g h
rm 2nd_test.out 2nd_best_hits
done | sed '1i Query No_exons No_locations Query_length Subject_length Bit_score %_Query_cov %_identity No_strands' | column -t > 2nd_gblast_query_search_table
k=`cat 2nd_gblast_query_search_table | grep -v "No_exons" | sort -k2,2 -k7,7 -nr | awk 'NR==1 {print$1}' | sed 's/APOBEC1_//g' | sed 's/_flanking.fa//g'`
echo $k $j
unset j ld x k
cd ~/bird_db1/aswin/APOBEC1
rm closest_query_list
done < <(grep -if 2nd_gblast_subjects_tmp /home/neo/bird_db1/aswin/database_details/all_genome_paths) | column -t > 2nd_gblast_query_subject_set_tmp

#Compare query-subject sets
join -1 1 -2 1 <(awk '{print$2,$1}' gblast_query_subject_set | sort -k1) <(awk '{print$2,$1}' 2nd_gblast_query_subject_set_tmp | sort -k1) | cat -n | awk '{if($3==$4) print$0,"same";else print$0,"diff"}' | column -t | grep ".*diff\|$" --color=always | GREP_COLORS="mt=32" grep ".*same\|$" --color=always > compare_gblast_2nd_gblast_tmp_query_subject_sets

#Updated squery subject list for 2nd gblast
awk '$NF~/diff/ {print$2}' compare_gblast_2nd_gblast_tmp_query_subject_sets > 2nd_gblast_subjects
awk '$NF~/diff/ {print$4,$2}' compare_gblast_2nd_gblast_tmp_query_subject_sets | column -t > 2nd_gblast_query_subject_set
#rm 2nd_gblast_subjects_tmp 2nd_gblast_query_subject_set_tmp

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#10.2. Run 2nd Genome blast, gedit, exonerate, exoneredit & spaln (18m)

while read line
do
species=`echo $line | awk '{print$1}' | cut -f1,2 -d "_"`
cd $(echo $line | awk '{print$2}')/aswin/APOBEC1
mkdir 2nd_gblast
cd 2nd_gblast
#Find the most suitable(nearest) query species
q=`grep $species /home/neo/bird_db1/aswin/APOBEC1/2nd_gblast_query_subject_set | awk '{print$1}'`
#Get the query sequence : In case the gene has very short exons : use sequences with all exons having flanking regions  (customizable)
#pq=`find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1 -name "*$q*.fa" | awk -F "/" '{print$NF}'`
pq=`find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/ -maxdepth 1 -name "*$q*.fa" | awk -F "/" '{print$NF}'`
pqwof=`find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1 -name "*$q*.fa" | awk -F "/" '{print$NF}'`
find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1  -name "*$q*.fa" -exec cp {} . \;
find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/ -name "*$q*_flanking.fa" -exec cp {} . \;

#Run gblast
echo "  ▬ " $species "(G)     ▬ " $pq "(Q)"
if [ -f ../../../genome/GC*.gff ]
then
#If the query sequence is "flanking.fa" then use -fix_query
gblast_short ../../../genome/GC*.fna $pq -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -eflankall ../../../genome/GC*.gff
else
gblast_short ../../../genome/GC*.fna $pq -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -eflankall
fi

#Create final consesnus from 2nd gblast
cp gblast_auto_consensus.fa final_consensus.fa

#------------------------------------------------------------------------------------------------------------#------------------------------------------------------------------------------------------------------------
#Run gedit if required
o1=`grep -v ">" gblast_auto_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
gaps=`grep -v ">" gblast_auto_consensus.fa | grep -i N -o | wc | awk '{print$3-$1}'`
anonym=`grep -v ">" gblast_auto_consensus.fa | tr -d 'ATGCNatgcn' | awk NF | wc | awk '{print$3-$1}'`
#splice sites
a1=`awk '{print"exon_"++i,toupper($1)}' splice_sites | sed 1d | awk '$2!="AG" {print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print $0}' | sed 's/exon_//g'`
a2=`awk '{print"exon_"++i,toupper($2)}' splice_sites | sed '$d' | awk '$2!="GT"' | awk '$2!="GC" {print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print$0}' | sed 's/exon_//g'`
if [[ ((-z $o1 && ! -z ${o1+x}) || $a1 != "-" || $a2 != "-") && $gaps == "0" && $anonym == "0" ]]
then
#Run gedit
echo "     ↳ CGblast didn't form complete ORF, running Gedit ..."
#If the query has only one exon then don't need to look for splice sites, just pairwise align & extract the aligned region
note=`grep ">" $pq -c`
if [[ $note == 1 ]]
then
start=`cat pairwise_exon_1.aln | grep -B2 "^extrac" | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":" | awk '{print$0+1}'`
#Extract the region from which the pairwise alignment start i.e. remove the regions before alignment starts (usually alignment starts with START codon ATG)
cat pairwise_exon_1.aln | grep -B2 "^extrac" | grep "|" -C1 | cut -c $start- > tmp1
end=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+1}'`
ln=`grep -v ">" $pqwof | awk NF | wc | awk '{print$3-$1}'`
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
ssc3=`expr $ssc1 + $ssc2 2>/dev/null | awk '{print$0-1}'`
if [[ $ssc3 == "" ]]; then er=`echo "        ↳ exon_1 - No 5 prime splice site found! "`; ssc3=`echo $ssc1 | awk '{print$0-1}'`; else :; fi
ln=`seqtk subseq $pqwof <(echo "exon_1") | grep -v ">" | awk NF | wc | awk '{print$3-$1}'`
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
echo "        ↳ Exon 1 alignment not found"
fi
last_exon=`grep ">" $pqwof | tail -1 | tr -d ">"`
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
ssc6=`expr $ssc4 - $ssc5 2>/dev/null | awk '{print$0+1}'`
if [[ $ssc6 == "" ]]; then echo "        ↳ "$j " - No 3 prime splice site found! " >> er1; ssc6=`echo $ssc4 | awk '{print$0+1}'`; else :; fi
cat $i | grep -B2 "^extrac" | grep "|" -C1 | cut -c $ssc6- > tmp1
ssc7=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+2}'`
#Look for nearest splice site "GT/GC" within next 5 flanking nucleotides after the query exon hit ended (the threshold flanking bases & other variable splice sites are customizable)
ssc8=`cat tmp1 | grep "|" -C1 | cut -c $ssc7- | tail -1 | cut -c -5 | grep -i "GT" -aob | head -1 | cut -f1 -d ":"`
#Allow splice site variants
if [[ $ssc8 == "" ]]; then ssc8=`cat tmp1 | grep "|" -C1 | cut -c $ssc7- | tail -1 | cut -c -5 | grep -i "GC" -aob | head -1 | cut -f1 -d ":"`; else :; fi
ssc9=`expr $ssc7 + $ssc8 2>/dev/null | awk '{print$0-1}'`
if [[ $ssc9 == "" ]]; then echo "        ↳ "$j " - No 5 prime splice site found! " >> er2; ssc9=`echo $ssc7 | awk '{print$0-1}'`; else :; fi
ln=`seqtk subseq $pqwof <(echo $j) | grep -v ">" | awk NF | wc | awk '{print$3-$1}'`
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
ssc12=`expr $ssc10 - $ssc11 2>/dev/null | awk '{print$0+1}'`
if [[ $ssc12 == "" ]]; then er=`echo "        ↳ "$last_exon " - No 3 prime splice site found! "`; ssc12=`echo $ssc10 | awk '{print$0+1}'`; else :; fi
ls pairwise_exon_* | grep $last_exon | xargs cat | grep -B2 "^extrac" | grep "|" -C1 | cut -c $ssc12- > tmp1
end=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+1}'`
ln=`seqtk subseq $pqwof <(echo $last_exon) | grep -v ">" | awk NF | wc | awk '{print$3-$1}'`
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
echo "        ↳ " $last_exon" alignment not found"
fi
#End of gedit
fi
#Gedit splice sites
#splice sites
for z in `grep ">" $pqwof | tr -d ">"`
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

#Generate ORF
transeq <(grep -v ">" gblast_edited_consensus.fa | paste -s -d "" | sed "1i \>$species") gblast_edited_consensus_orf.fa 2>/dev/null
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" gblast_edited_consensus.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" gblast_edited_consensus_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > gblast_edited_translated
#For better visualization of ORF
gnaa=`head -1 gblast_edited_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" gblast_edited_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 gblast_edited_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > gblast_edited_nucleotide_amino_acid_view
#Create final consesnus from gedit
cp gblast_edited_consensus.fa final_consensus.fa
unset o1 gaps anonym a1 a2 gnaa
#Gedit exon boundary indels
for pe in `ls pairwise_exon_* | sort -V`
do
pen=`echo $pe | sed 's/.aln//g' | cut -f2- -d "_"`
#Trim flanking regions outside alignment
start=`grep -B2 "^extrac" $pe | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":" | awk '{print$0+1}'`
grep -B2 "^extrac" $pe | grep "|" -C1 | cut -c $start- > tmp1
end=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+1}'`
cat tmp1 | cut -c -$end > tmp2
#insertions at exon start : counted if a dash is observed within the query sequence anywhere within the first 5% of pairwise alignment between query & subject
is=`cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/query/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $is == "" ]]; then is="-"; else :; fi
#insertions at exon end : counted if a dash is observed within the query sequence anywhere after 95% of pairwise alignment between query & subject
ie=`cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/query/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ie == "" ]]; then ie="-"; else :; fi
#deletions at exon start : counted if a dash is observed within the subject sequence anywhere within the first 10% (since subject sequence is longer due to flanking regions) of pairwise alignment between query & subject
ds=`cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/genome/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ds == "" ]]; then ds="-"; else :; fi
#deletions at exon end : counted if a dash is observed within the subject sequence anywhere after 90% of pairwise alignment between query & subject
de=`cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/genome/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
#Longest continous insertions : Position, Length & % of gene in which it appears
if [[ $(head -1 tmp2 | grep "\-") == "" ]]; then lci="- - -"; else
lci=`bedtools merge -i <(cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/query/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
fi
#Longest continous deletions : Position, Length & % of gene in which it appears
if [[ $(tail -1 tmp2 | grep "\-") == "" ]]; then lcd="- - -"; else
lcd=`bedtools merge -i <(cat tmp2 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"genome",$0; else print"query",$0}' | awk '/genome/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
#lci=`cat tmp2 | grep "\-\+" -o | tr -d " " | awk '{print length}' | sort -nr | head -1`
fi
if [[ $de == "" ]]; then de="-"; else :; fi
echo $pen $is $ie $ds $de $lci $lcd
unset pen start end is ie ds de lci lcd
rm tmp1 tmp2
done | sed '1i Exons Ins_in_1st_5%_exon Ins_in_last_95%_exon Del_in_1st_5%_exon Del_in_last_95%_exon Longest_continous_ins_position(lci) lci_length lci_%_range Longest_continous_del_position(lcd) lcd_length lcd_%_range' | column -t > gedit_exon_boundary_indels
unset note

#------------------------------------------------------------------------------------------------------------#------------------------------------------------------------------------------------------------------------
#Run exonerate if required : criterias : No ORF, No gaps, No anonymous bases from gedit, No no hits from blast & no splice site disruptions
o2=`grep -v ">" gblast_edited_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
gaps=`grep -v ">" gblast_edited_consensus.fa | grep -i N -o | wc | awk '{print$3-$1}'`
anonym=`grep -v ">" gblast_edited_consensus.fa | tr -d 'ATGCNatgcn' | awk NF | wc | awk '{print$3-$1}'`
nh=`grep ">" $pqwof -c | xargs -n1 bash -c 'expr $0 - $(grep ">" gblast_auto_consensus.fa -c)'`
#splice sites
a1=`awk '{print"exon_"++i,toupper($1)}' gedit_splice_sites | sed 1d | awk '$2!="AG" {print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print $0}' | sed 's/exon_//g'`
a2=`awk '{print"exon_"++i,toupper($2)}' gedit_splice_sites | sed '$d' | awk '$2!="GT"' | awk '$2!="GC"' | awk '{print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print$0}' | sed 's/exon_//g'`
#if "z" variable is set & empty not just empty
if [[ ((-z $o2 && ! -z ${o2+x}) || $a1 != "-" || $a2 != "-") && $gaps == "0" && $anonym == "0" && $nh == "0" ]]
then
#Run Exonerate
echo "     ↳ Gedit didn't form complete ORF, running Exonerate ..."
#Run exonerate
#Translate to amino acid ("*" is required to represent stop codon for exonerate)
transeq <(grep -v ">" $pqwof) -auto -stdout | sed "/>/ s/.*/>$pqof/g" | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > exonerate_query.aa 2>/dev/null
#Combined extracted sequence as reference (optional: to get Splice sites & for QC)
grep -v ">" extracted_flanking_region.fa | paste -s -d "" | tr '[:lower:]' '[:upper:]' | sed '1i >extracted_flanking_region_combined' > extracted_flanking_region_combined.fa
#Run exonerate in exhaustive mode  with model protein2genome:bestfit & maximum intron allowed is 36bp
exonerate -E --model p2g:b exonerate_query.aa extracted_flanking_region_combined.fa --maxintron 48 --ryo "Rank: %r\nTarget in alignment: %ti : %tal (%tab - %tae)\nTarget in coding sequence: %ti : %tcl (%tcb - %tce)\nGene orientation: %g\nPercent identity: %pi\nPercent similarity: %ps\n>%ti (%tab - %tae)\n%tcs\n" --querytype protein --targettype dna --showtargetgff yes --fsmmemory 28000 --alignmentwidth 270 > Query_exon_combined_exonerate.out 2>/dev/null
#Create consensus sequence
#sed 1d Query_exon_combined_exonerate.out | awk '/Query range/; /^>/,/^$/' | grep -v "^>" | sed 's/.*Query range: />/g' | sed 's/ -> .*//g' | sed 's/^>/\x00&/' | sort -z -V | sort -z -V | tr -d '\0' | grep -v ">" | awk NF | paste -s -d "" | sed "1i >Query_exon_combined_exonerate" > Query_exon_combined_exonerate.out.fa
#Exons_wise sequence
bedtools getfasta -fi extracted_flanking_region_combined.fa -bed <(awk '$3=="exon" {print"extracted_flanking_region_combined",$4-1,$5,$3 | "sort -k2 -n"}' Query_exon_combined_exonerate.out | awk '{print$1,$2,$3,$4"_"++i}' OFS="\t") -name > Query_exon_combined_exonerate.out_exon_wise.fa 2>/dev/null
grep -v ">" Query_exon_combined_exonerate.out_exon_wise.fa | sed '1i >exonerate_consensus' > Query_exon_combined_exonerate.out.fa
#Check ORF
transeq <(grep -v ">" Query_exon_combined_exonerate.out.fa | paste -s -d "") Query_exon_combined_exonerate.out.fa_orf.fa 2>/dev/null
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" Query_exon_combined_exonerate.out.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" Query_exon_combined_exonerate.out.fa_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > Query_exon_combined_exonerate_translated
#For better visualization
gnaa=`head -1 Query_exon_combined_exonerate_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" Query_exon_combined_exonerate_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 Query_exon_combined_exonerate_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > Query_exon_combined_exonerate_nucleotide_amino_acid_view
#Create consensus from exonerate
cp Query_exon_combined_exonerate.out_exon_wise.fa final_consensus.fa
unset o2 gaps anonym nh a1 a2 gnaa
#Exonerate exon boundary indels : run only if exonerate gave hits on all queries
enh=`grep ">" $pqwof -c | xargs -n1 bash -c 'expr $0 - $(grep ">" Query_exon_combined_exonerate.out_exon_wise.fa -c)'`
if [[ $enh == "0" ]]
then
while mapfile -t -n 4 ary && ((${#ary[@]}))
do
ex=`printf '%s\n' "${ary[@]}" | head -1 | tr -d ">"`
printf '%s\n' "${ary[@]}" > tmp
cat tmp | grep "Query_exon" | grep [0-9] -aob | head -1 | cut -f1 -d ":" | awk '{print$1+3}' | xargs -n1 sh -c 'cat tmp | cut -c $0-' 2>/dev/null | awk NF | sed 's/ \+[0-9]\+//g' > tmp1
#insertions at exon start : counted if a dash is observed within the query sequence anywhere within the first 5% of pairwise alignment between query & subject
is=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"Query_exon",$0; else print"APOBEC1",$0}' | awk '/APOBEC1/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $is == "" ]]; then is="-"; else :; fi
#insertions at exon end : counted if a dash is observed within the query sequence anywhere after 95% of pairwise alignment between query & subject
ie=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"Query_exon",$0; else print"APOBEC1",$0}' | awk '/APOBEC1/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ie == "" ]]; then ie="-"; else :; fi
#deletions at exon start : counted if a dash is observed within the subject sequence anywhere within the first 10% (since subject sequence is longer due to flanking regions) of pairwise alignment between query & subject
ds=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"Query_exon",$0; else print"APOBEC1",$0}' | awk '/Query_exon/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ds == "" ]]; then ds="-"; else :; fi
#deletions at exon end : counted if a dash is observed within the subject sequence anywhere after 90% of pairwise alignment between query & subject
de=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"Query_exon",$0; else print"APOBEC1",$0}' | awk '/Query_exon/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
#Longest continous insertions : Position, Length & % of gene in which it appears
if [[ $(head -1 tmp1 | grep "\-") == "" ]]; then lci="- - -"; else
lci=`bedtools merge -i <(cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"Query_exon",$0; else print"APOBEC1",$0}' | awk '/APOBEC1/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
fi
#Longest continous deletions : Position, Length & % of gene in which it appears
if [[ $(tail -1 tmp1 | grep "\-") == "" ]]; then lcd="- - -"; else
lcd=`bedtools merge -i <(cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"Query_exon",$0; else print"APOBEC1",$0}' | awk '/Query_exon/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
#lci=`cat tmp2 | grep "\-\+" -o | tr -d " " | awk '{print length}' | sort -nr | head -1`
fi
if [[ $de == "" ]]; then de="-"; else :; fi
echo $ex $is $ie $ds $de $lci $lcd
unset ex is ie ds de lci lcd
rm tmp tmp1
done < <(exalign $pqwof Query_exon_combined_exonerate.out_exon_wise.fa 2>/dev/null | sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" | awk NF) | sed '1i Exons Ins_in_1st_5%_exon Ins_in_last_95%_exon Del_in_1st_5%_exon Del_in_last_95%_exon Longest_continous_ins_position(lci) lci_length lci_%_range Longest_continous_del_position(lcd) lcd_length lcd_%_range' | column -t > exonerate_exon_boundary_indels
else
echo "     ↳ Exonerate didn't gave hit on all exons, terminating further search !"
fi

#------------------------------------------------------------------------------------------------------------#------------------------------------------------------------------------------------------------------------
#Run exoneredit if required
o3=`grep -v ">" Query_exon_combined_exonerate.out.fa_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
gaps=`grep -v ">" Query_exon_combined_exonerate.out_exon_wise.fa | grep -i N -o | wc | awk '{print$3-$1}'`
anonym=`grep -v ">" Query_exon_combined_exonerate.out_exon_wise.fa | tr -d 'ATGCNatgcn' | awk NF | wc | awk '{print$3-$1}'`
#splice sites
splice_sites=`grep ">" $pqwof -c | awk '{print$0-1}'`
a1t=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$1=="GT" || $1=="GC"' | wc -l`
if [[ $a1t == $splice_sites ]]; then a1="-"; else a1="x"; fi
a2t=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$2=="AG"' | wc -l`
if [[ $a2t == $splice_sites ]]; then a2="-"; else a2="x"; fi
if [[ ((-z $o3 && ! -z ${o3+x}) || $a1 != "-" || $a2 != "-") && $gaps == "0" && $anonym == "0" && $enh == "0" ]]
then
#Run exoneredit
echo "     ↳ Exonerate didn't form complete ORF, running Exoneredit ..."
#pairwise align consensus sequences
needle -asequence <(grep -v ">" $pqwof | sed '1i >Query') Query_exon_combined_exonerate.out.fa -auto -stdout -awidth 280 > pairwise_query_exoneredit.aln
#find START site
start=`needle -asequence <(grep -v ">" $pqwof | sed '1i >Query') <(grep -v ">" Query_exon_combined_exonerate.out.fa | sed '1i >Exonerate') -auto -stdout -awidth 100000 | grep "^Query" -A2 | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":" | awk '{print$0+1}'`
#Cut everything before START site
needle -asequence <(grep -v ">" $pqwof | sed '1i >Query') <(grep -v ">" Query_exon_combined_exonerate.out.fa | sed '1i >Exonerate') -auto -stdout -awidth 100000 | grep "^Query" -A2 | grep "|" -C1 | cut -c $start- > tmp1
#Find STOP site
end=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+1}'`
#Calculate some basic stats about exoneredit before making consensus sequence
ln=`grep -v ">" $pqwof  | wc | awk '{print$3-$1}'`
#cat tmp1 | cut -c -$end | head -1 | tr -d "-" | wc | awk '{print$3-$1}' | awk -v p="$ln" '{print p,$1,($1/p)*100}' | sed "s/^/exoneredit /g" > exoneredit_query_covered
#insertions
cat tmp1 | cut -c -$end | head -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exoneredit /g" > exoneredit_insertions
#deletions
cat tmp1 | cut -c -$end | tail -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exoneredit /g" > exoneredit_deletions
#Create consensus : Delete everything after STOP site
cat <(echo ">exoneredit") <(cat tmp1 | cut -c -$end | tail -1) | tr -d "-" > exoneredit_consensus.fa
#Extract exon_wise sequence
blat extracted_flanking_region_combined.fa exoneredit_consensus.fa -out=pslx exoneredit.pslx > /dev/null 2>&1
egrep -v "psLayout|match|^---" exoneredit.pslx | awk NF | awk '{print$NF}' | tr "," "\n" | tr '[:lower:]' '[:upper:]' | awk NF | awk '{print ">exon_"++i"\n"$0}' > exoneredit_consensus_exon_wise.fa
rm tmp1
#Generate ORF
transeq <(grep -v ">" exoneredit_consensus.fa | paste -s -d "" | sed "1i \>$species") exoneredit_consensus_orf.fa 2>/dev/null
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" exoneredit_consensus.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" exoneredit_consensus_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > exoneredit_translated
#For better visualization of ORF
gnaa=`head -1 exoneredit_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" exoneredit_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 exoneredit_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > exoneredit_nucleotide_amino_acid_view
#Create consensus from exoneredit
cp exoneredit_consensus_exon_wise.fa final_consensus.fa
unset o3 gaps anonym splice_sites a1t a1 a2t a2 start ln p end gnaa
else
echo "     ↳ Exoneredit didn't run"
unset o3 gaps anonym splice_sites a1t a1 a2t a2
fi
unset o3 enh
#------------------------------------------------------------------------------------------------------------#------------------------------------------------------------------------------------------------------------
#Run spaln if required

#Sometimes exoneredit is not run even if exonerate didn't form complete ORF; hence the conditions to run spaln will depend on whether exoneredit ran or not
if [[ -e exoneredit_consensus_orf.fa ]]
then
o4=`grep -v ">" exoneredit_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
gaps=`grep -v ">" exoneredit_consensus_exon_wise.fa | grep -i N -o | wc | awk '{print$3-$1}'`
anonym=`grep -v ">" exoneredit_consensus_exon_wise.fa | tr -d 'ATGCNatgcn' | awk NF | wc | awk '{print$3-$1}'`
else
o4=`grep -v ">" Query_exon_combined_exonerate.out.fa_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
gaps=`grep -v ">" Query_exon_combined_exonerate.out_exon_wise.fa | grep -i N -o | wc | awk '{print$3-$1}'`
anonym=`grep -v ">" Query_exon_combined_exonerate.out_exon_wise.fa | tr -d 'ATGCNatgcn' | awk NF | wc | awk '{print$3-$1}'`
fi
#no hits from gblast
nh=`grep ">" $pqwof -c | xargs -n1 bash -c 'expr $0 - $(grep ">" gblast_auto_consensus.fa -c)'`
#splice sites
splice_sites=`grep ">" $pqwof -c | awk '{print$0-1}'`
a1t=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$1=="GT" || $1=="GC"' | wc -l`
if [[ $a1t == $splice_sites ]]; then a1="-"; else a1="x"; fi
a2t=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$2=="AG"' | wc -l`
if [[ $a2t == $splice_sites ]]; then a2="-"; else a2="x"; fi
if [[ ((-z $o4 && ! -z ${o4+x}) || $a1 != "-" || $a2 != "-") && $gaps == "0" && $anonym == "0" && $nh == "0" ]]
then
#Run Spaln
echo "     ↳ Exonerate/Exoneredit didn't form complete ORF, running Spaln ..."
#Translate to amino acid ("*" is required to represent stop codon for Spaln)
transeq <(grep -v ">" $pqwof) -auto -stdout | sed "/>/ s/.*/>$pqwof/g" | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > spaln_query.aa
#Query exon count
pqe1=`grep ">" $pqwof -c`
#Query cds length
pqe2=`grep -v ">" $pqwof | wc | awk '{print$3-$1}'`
#Splice sites expected in target
spn=`echo $pqe1 | awk '{print$1-1}'`
#Run Spaln with different combinations of -yX(intr/cross species parameters) -yS (species-specific parameter)
for x in 0 1
do
for y in 0 50 100
do
echo -n > spaln_out_tmp
echo -n > spaln_tmp
spaln -M1 -l280 -O4 -ya0 -yX$x -yS$y extracted_flanking_region_combined.fa spaln_query.aa > spaln_out_tmp 2>/dev/null
bedtools getfasta -fi extracted_flanking_region_combined.fa -bed <(awk '!/^#|^@/{print"extracted_flanking_region_combined",$9-1,$10,"exon_"++i}' OFS="\t" spaln_out_tmp) -name > spaln_tmp 2>/dev/null
x1=`grep ">" -c spaln_tmp`
x2=`grep -v ">" spaln_tmp | wc | awk '{print$3-$1}'`
x3=`needle <(grep -v ">" spaln_tmp) <(grep -v ">" $pqwof) -auto -stdout -awidth 280 | awk '/Identity|Similarity|Gaps|Score/ {print$NF}' | tr -d ')(%' | paste -s -d " "` 2>/dev/null
transeq <(grep -v ">" spaln_tmp | paste -s -d "" | sed "1i \>$species") spaln_tmp_orf.fa 2>/dev/null
#Splice sites
sp1=`awk '!/^#|^@/ {print$NF}' spaln_out_tmp | tr "." " " | awk NF | tr '[:lower:]' '[:upper:]' | awk '$1=="GT" || $1=="GC"' | wc -l`
if [[ $sp1 == $spn ]]; then sp1s="✔"; else sp1s="x"; fi
sp2=`awk '!/^#|^@/ {print$NF}' spaln_out_tmp | tr "." " " | awk NF | tr '[:lower:]' '[:upper:]' | awk '$2=="AG"' | wc -l`
if [[ $sp2 == $spn ]]; then sp2s="✔"; else sp2s="x"; fi
#Check ORF presence absence
x4=`grep -v ">" spaln_tmp_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $x4 == "" ]] ; then x5="❌";else x5="✅"; fi
echo $x $y $pqe1 $x1 $pqe2 $x2 $x3 $sp1 $sp1s $sp2 $sp2s $x5
unset x1 x2 x3 x4 sp1 sp1s sp2 sp2s x5
rm spaln_out_tmp spaln_tmp spaln_tmp_orf.fa
done
done | sed '1i yX yS Query_exons Spaln_exons Query_length Spaln_length %_Identity %_Similarity %_Gaps Score No_Splice_site1 SS1_disruption No_Splice_site2 SS2_disruption ORF' | column -t > spaln_summary
unset x y

#Choose best parameters based on exon count match & alignment stats
spaln_params=`cat <(awk '$3==$4 && $12=="✔" && $14=="✔"' spaln_summary | grep "✅") <(awk '$3==$4' spaln_summary | grep "✅") <(awk '$3==$4' spaln_summary | grep "❌") | awk '!a[$0]++' | awk 'NR==1{print$1,$2}'`
if [[ $spaln_params == "" ]]
then
unset spaln_params
spaln_params=`grep -v "yX" spaln_summary | sort -k4,4 -nr | awk 'NR==1{print$1,$2}'`
echo "        ↳ No ideal Spaln parameters found! Taking the parameters that gave the highest number of exon hits : " $spaln_params
else
echo "        ↳ Spaln parameters used : yX & yS : "$spaln_params
fi
while read sp
do
x=`echo $sp | awk '{print$1}'`
y=`echo $sp | awk '{print$2}'`
#Final spaln output in multiple formats : each with different info
spaln -M1 -l280 -O4 -ya0 -yX$x -yS$y extracted_flanking_region_combined.fa spaln_query.aa > spaln_out 2>/dev/null
spaln -M1 -l280 -O1 -ya0 -yX$x -yS$y extracted_flanking_region_combined.fa spaln_query.aa > spaln_aln_out 2>/dev/null
spaln -M1 -l280 -O2 -ya0 -yX$x -yS$y extracted_flanking_region_combined.fa spaln_query.aa > spaln_gff3_out 2>/dev/null
bedtools getfasta -fi extracted_flanking_region_combined.fa -bed <(awk '!/^#|^@/{print"extracted_flanking_region_combined",$9-1,$10,"exon_"++i}' OFS="\t" spaln_out) -name > spaln_exon_wise.fa 2>/dev/null
awk '!/^#|^@/ {print$NF}' spaln_out | tr "." " " | awk NF | tr '[:lower:]' '[:upper:]' > spaln_splice_sites
done < <(echo $spaln_params)

#Generate ORF
transeq <(grep -v ">" spaln_exon_wise.fa | paste -s -d "" | sed "1i \>$species") spaln_exon_wise_orf.fa 2>/dev/null
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" spaln_exon_wise.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" spaln_exon_wise_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > spaln_translated
#For better visualization of ORF
gnaa=`head -1 spaln_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" spaln_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 spaln_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > spaln_nucleotide_amino_acid_view
#Pairwise exon-wise algnment between query & final consensus
exalign $pqwof spaln_exon_wise.fa > spaln_query_hit_pairwise_exon.aln 2>/dev/null
#QC spaln output : Exon boundar indels
#Exon boundary insertions
while mapfile -t -n 4 ary && ((${#ary[@]}))
do
ex=`printf '%s\n' "${ary[@]}" | head -1 | tr -d ">"`
printf '%s\n' "${ary[@]}" > tmp
cat tmp | grep "spaln_exon" | grep [0-9] -aob | head -1 | cut -f1 -d ":" | awk '{print$1+3}' | xargs -n1 sh -c 'cat tmp | cut -c $0-' 2>/dev/null | awk NF | sed 's/ \+[0-9]\+//g' > tmp1
#insertions at exon start : counted if a dash is observed within the query sequence anywhere within the first 5% of pairwise alignment between query & subject
is=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/APOBEC1/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $is == "" ]]; then is="-"; else :; fi
#insertions at exon end : counted if a dash is observed within the query sequence anywhere after 95% of pairwise alignment between query & subject
ie=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/APOBEC1/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ie == "" ]]; then ie="-"; else :; fi
#deletions at exon start : counted if a dash is observed within the subject sequence anywhere within the first 10% (since subject sequence is longer due to flanking regions) of pairwise alignment between query & subject
ds=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/spaln_exon/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ds == "" ]]; then ds="-"; else :; fi
#deletions at exon end : counted if a dash is observed within the subject sequence anywhere after 90% of pairwise alignment between query & subject
de=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/spaln_exon/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
#Longest continous insertions : Position, Length & % of gene in which it appears
if [[ $(head -1 tmp1 | grep "\-") == "" ]]; then lci="- - -"; else
lci=`bedtools merge -i <(cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/APOBEC1/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
fi
#Longest continous deletions : Position, Length & % of gene in which it appears
if [[ $(tail -1 tmp1 | grep "\-") == "" ]]; then lcd="- - -"; else
lcd=`bedtools merge -i <(cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"APOBEC1",$0}' | awk '/spaln_exon/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
#lci=`cat tmp2 | grep "\-\+" -o | tr -d " " | awk '{print length}' | sort -nr | head -1`
fi
if [[ $de == "" ]]; then de="-"; else :; fi
echo $ex $is $ie $ds $de $lci $lcd
unset ex is ie ds de lci lcd
rm tmp tmp1
done < <(cat spaln_query_hit_pairwise_exon.aln | sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" | awk NF) | sed '1i Exons Ins_in_1st_5%_exon Ins_in_last_95%_exon Del_in_1st_5%_exon Del_in_last_95%_exon Longest_continous_ins_position(lci) lci_length lci_%_range Longest_continous_del_position(lcd) lcd_length lcd_%_range' | column -t > spaln_boundary_indels_wrt_query
#Create consensus from spaln
cp spaln_exon_wise.fa final_consensus.fa
unset o4 gaps anonym nh splice_sites a1t a1 a2t a2 pqe1 pqe2 spn sp x y spaln_params gnaa r
else
echo "     ↳ Spaln didn't run"
unset o4 gaps anonym nh splice_sites a1t a1 a2t a2
fi

#-----------------------------------------------------------------------------------------------------------
else
#If Exonerate didn't run then it means following reasons : chose which one it is
o2=`grep -v ">" gblast_edited_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ ! -z $o2 && ! -z ${o2+x} ]]; then fo="✅"; else fo="❌"; fi
gaps=`grep -v ">" gblast_edited_consensus.fa | grep -i N -o | wc | awk '{print$3-$1}'`
anonym=`grep -v ">" gblast_edited_consensus.fa | tr -d 'ATGCNatgcn' | awk NF | wc | awk '{print$3-$1}'`
nh=`grep ">" $pqwof -c | xargs -n1 bash -c 'expr $0 - $(grep ">" gblast_auto_consensus.fa -c)'`
echo "        ↳ Terinating search due to : "
echo "           ↳ Open Reading Frame from gedit   : "$fo
echo "           ↳ Exons with No hits from gblast  : "$nh
echo "           ↳ No of gaps found in genome      : "$gaps
echo "           ↳ No of anonymous bases in genome : "$anonym
unset o2 fo gaps anonym nh
fi
else
echo "        ↳ Gblast found ORF ✅"
fi
echo -e
unset j q pq pqwof z
#Go back to home directory of Gene of interest - customize
cd /home/neo/bird_db1/aswin/APOBEC1
done < <(grep -if 2nd_gblast_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#10.3. 2nd Genome blast summary

while read i
do
#Species of interest
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
#Taxonomy order of the species
gr=`echo $j | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} /home/neo/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/2nd_gblast/
#Query used
q1=`grep $j /home/neo/bird_db1/aswin/APOBEC1/2nd_gblast_query_subject_set | awk '{print$1}'`
q=`ls *.fa | grep $q1 | grep -v "flanking"`
#Exons with no hits
#from tblastx
nht=`grep -if <(grep -v "Query" tblastx.out | awk NF | awk '!a[$1]++ {print$1}') <(grep ">" $q | tr -d ">" ) -v | sed 's/exon_//g' | paste -s -d ","`
if [[ $nht == "" ]]; then nht="-"; else :; fi
#from gblast
nh1=`grep "No hits found" test.outfmt3 -B5 | grep "Query=" | awk '{print$NF}' | paste -s -d "," | sed 's/exon_//g'`
if [[ $nh1 == "" ]]; then nh1="-"; else :; fi
#from exonerate
if [[ -s Query_exon_combined_exonerate.out_exon_wise.fa ]]
then
nh2=`exalign gblast_edited_consensus.fa Query_exon_combined_exonerate.out_exon_wise.fa -m 2>/dev/null | sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" | grep ">" | tr ")(|%" " " | column -t | awk '$3<80 || $3=="" {print$1}' | sed 's/>exon_//g' | paste -s -d ","`
if [[ $nh2 == "" ]]; then nh2="-"; else :; fi
else
nh2="-"
fi
#from spaln
if [[ -s spaln_exon_wise.fa ]]
then
nh3=`grep -if <(grep ">" spaln_exon_wise.fa) <(grep ">" $q) -v | sed 's/>exon_//g' | paste -s -d ","`
if [[ $nh3 == "" ]]; then nh3="-"; else :; fi
else
nh3="-"
fi
#Number of exons with almost complete (>90% covered) query hit
ewah=`awk '$13>90 {print$1}' best_hits | wc -l`
#START codon
sc=`grep -v ">" final_consensus.fa | paste -s -d "" | cut -c -3 | sed 's/[a-z]/\U&/g'`
if [[ $sc == "" ]]; then sc="-"; else :; fi
#Check location
d=`sed 1d best_hits | awk '!seen[$1]++' | awk '{print$2}' | sort -u | wc -l`
if [ $d == "1" ]; then a3="same"; else a3="diff"; fi
#Check strandedness
e=`sed 1d best_hits | awk '!seen[$1]++' | awk '{print$NF}' | sort -u | wc -l`
if [ $e == "1" ]; then a4="same"; else a4="diff"; fi
#duplicated exons
de=`awk '{print$1}' duplicated_exons | sort -u | paste -s -d "," | sed 's/exon_//g'`
if [[ $de == "" ]]; then de="-"; else :; fi
#Paralogs
if [[ -s paralogs ]]; then pg=`cat paralogs | wc -l`; else pg="-"; fi
#Pseudogenes if GFF available
if [[ -s pseudogenes ]]; then psg="yes"; else psg="no"; fi
#Insertions
if [[ -f spaln_exon_wise.fa ]]
then
m=`needle <(grep -v ">" spaln_exon_wise.fa | sed '1i >spaln') <(grep -v ">" $q | sed '1i >query') -auto -stdout -awidth 280 | grep "^query" | awk '{print$3}' | tr -cd "-" | wc | awk '{print$3-$1}'`
if [[ $m == "" ]]; then m="-"; else :; fi
elif [[ -f exoneredit_insertions ]]
then
m=`awk '{print$2}' exoneredit_insertions`
elif [[ -f Query_exon_combined_exonerate.out ]]
then
m=`cat Query_exon_combined_exonerate.out | awk '{for(i=1;i<=NF;i++) if($i~/insertions/) print$(i+1)}' | awk '{sum+=$1;} END{print sum}' | awk '{print$1}'`
if [[ $m == "" ]]; then m="-"; else :; fi
elif [[ -f gblast_edited_deletions ]]
then
m=`awk '{a+=$2} END{print a}' gblast_edited_insertions`
else
m="-"
fi
#Deletions
if [[ -f spaln_exon_wise.fa ]]
then
n=`needle <(grep -v ">" spaln_exon_wise.fa | sed '1i >spaln') <(grep -v ">" $q | sed '1i >query') -auto -stdout -awidth 280 | grep "^spaln" | awk '{print$3}' | tr -cd "-" | wc | awk '{print$3-$1}'`
if [[ $n == "" ]]; then n="-"; else :; fi
elif [[ -f exoneredit_deletions ]]
then
n=`awk '{print$2}' exoneredit_deletions`
elif [[ -f Query_exon_combined_exonerate.out ]]
then
n=`cat Query_exon_combined_exonerate.out | awk '{for(i=1;i<=NF;i++) if($i~/deletions/) print$(i+1)}' | awk '{sum+=$1;} END{print sum}' | awk '{print$1*3}'`
if [[ $n == "" ]]; then n="-"; else :; fi
elif [[ -f gblast_edited_deletions ]]
then
n=`awk '{a+=$2} END{print a}' gblast_edited_deletions`
else
n=`cat test.out | sed 1,2d | awk '!seen[$1]++' |  awk '{print$NF}' | tr -cd '-' | wc | awk '{print$3-$1}'`
fi
#Gaps
#Note: from exonerate output we can know whether splice site is disrupted or not but not which splice sites!
gp=`sed '/^>/! s/[a-z]/\U&/g' final_consensus.fa | grep N -B1 | tr -d [ATGC] | sed 's/^--//g' | awk NF | paste -d " " - - | sed 's/>exon_//g' | awk '{print length($2)"("$1")"}' | paste -s -d ","`
if [[ $gp == "" ]]; then gp="-"; else :; fi
#Other anonymous bases other than A,T,G,C & N in final consensus sequence : Number of anonymous bases (exon number)
ons=`cat final_consensus.fa | grep -v ">" | tr -d "[a,t,g,c,n,A,T,G,C,N]" | awk NF | xargs -n1 sh -c 'grep $0 final_consensus.fa -c; grep $0 final_consensus.fa -B1 | grep ">" | sed "s/>exon_//g"' | paste -d " " - - | awk '{print$1"("$2")"}' | paste -s -d ","`
if [[ $ons == "0()" ]]; then ons="-"; else :; fi
#Splice sites
splice_sites=`grep ">" $q -c | awk '{print$0-1}'`
if [[ -s spaln_exon_wise.fa ]]
then
#splice site 1
a1t=`awk '$1=="GT" || $1=="GC"' spaln_splice_sites | wc -l`
if [[ $a1t == $splice_sites ]]; then a1="-"; else a1="x"; fi
#splice site 2
a2t=`awk '$2=="AG"' spaln_splice_sites | wc -l`
if [[ $a2t == $splice_sites ]]; then a2="-"; else a2="x"; fi
elif (( $(cat Query_exon_combined_exonerate.out 2>/dev/null | wc -l) > 5 ))
then
a1t=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$1=="GT" || $1=="GC"' | wc -l`
if [[ $a1t == $splice_sites ]]; then a1="-"; else a1="x"; fi
a2t=`cat Query_exon_combined_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$2=="AG"' | wc -l`
if [[ $a2t == $splice_sites ]]; then a2="-"; else a2="x"; fi
elif [[ -f gedit_splice_sites ]]
then
a1=`awk '{print"exon_"++i,toupper($1)}' gedit_splice_sites | sed 1d | awk '$2!="AG" {print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print $0}' | sed 's/exon_//g'`
a2=`awk '{print"exon_"++i,toupper($2)}' gedit_splice_sites | sed '$d' | awk '$2!="GT"' | awk '$2!="GC"' | awk '{print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print$0}' | sed 's/exon_//g'`
elif [[ -f splice_sites ]]
then
a1=`awk '{print"exon_"++i,toupper($1)}' splice_sites | sed 1d | awk '$2!="AG" {print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print $0}' | sed 's/exon_//g'`
a2=`awk '{print"exon_"++i,toupper($2)}' splice_sites | sed '$d' | awk '$2!="GT"' | awk '$2!="GC" {print$1}' | paste -s -d "," | awk '{if($0=="") print"-"; else print$0}' | sed 's/exon_//g'`
else :
fi
#Consensus nucleotide length
ncon=`grep -v ">" final_consensus.fa | wc | awk '{print$3-$1}'`
#Consensus protein length
pcon=`transeq final_consensus.fa -auto -stdout | grep -v ">" | wc | awk '{print$3-$1}'`
#Query covered
if [[ -s spaln_exon_wise.fa ]]
then
ql=`grep -v ">" spaln_query.aa | wc | awk '{print$3-$1}'`
qr=`cat spaln_aln_out | grep "^>extracted_flanking_region_combined" | grep -oP '\(\K[^\)]+' | tr -d " " | awk 'NR==2'`
qc=`echo $qr | awk -v l="$ql" -F "-" '{print($2-$1+1)/l*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
elif (( $(grep -v ">" Query_exon_combined_exonerate.out.fa 2>/dev/null | wc | awk '{print$3-$1}') > 0 ))
then
ins=`awk '{for(i=1;i<=NF;i++) if($i~/insertions/) print $(i+1)}' Query_exon_combined_exonerate.out | awk '{s+=$1} END{print s}' | awk '{print$1}'`
del=`awk '{for(i=1;i<=NF;i++) if($i~/deletions/) print $(i+1)}' Query_exon_combined_exonerate.out | awk '{s+=$1} END{print s}' | awk '{print$1*3}'`
fs=`awk '{for(i=1;i<=NF;i++) if($i~/frameshifts/) print $(i+1)}' Query_exon_combined_exonerate.out | awk '{s+=$1} END{print s}'`
if [[ $fs == "" ]]; then fs="-"; else :; fi
tl=`grep -v ">" Query_exon_combined_exonerate.out.fa | wc | awk '{print$3-$1}'`
ql=`grep -v ">" $q | wc | awk '{print$3-$1}'`
if [[ $fs == "-" ]]; then
qc=`calc $(calc $tl - $ins + $del) / $ql \* 100 | tr -d "~" | awk '{$1+=0}1' CONVFMT="%.2f"`
else
qc=`calc $(calc $tl - $ins + $del - $fs) / $ql \* 100 | tr -d "~" | awk '{$1+=0}1' CONVFMT="%.2f"`
fi
elif [[ -f gblast_edited_query_covered ]]
then
tl=`cat $q | grep -v ">" | wc | awk '{print$3-$1}'`
qc=`awk '{print$3}' gblast_edited_query_covered | awk '{sum+=$1;} END{print sum;}' | awk -v l="$tl" '{print($1/l)*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
else
tl=`cat $q | grep -v ">" | wc | awk '{print$3-$1}'`
qc=`awk 'NR>1{print$6-$5+1}' best_hits | awk '{sum+=$1;} END{print sum;}' | awk -v l="$tl" '{print($1/l)*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
fi
unset tl
# % of gene in which 1st STOP codon appears
if [[ -s spaln_exon_wise.fa ]]
then
tl=`grep -v ">" spaln_query.aa | wc | awk '{print$3-$1}'`
sp=`tail -1 spaln_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$tl" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.1f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
elif (( $(grep "[A-Z]" exoneredit_translated 2>/dev/null | wc -l) > 1 ))
then
tl=`grep -v ">" exoneredit_consensus_orf.fa | wc | awk '{print$3-$1}'`
sp=`tail -1 exoneredit_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$tl" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.1f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
elif (( $(grep "[A-Z]" Query_exon_combined_exonerate_translated 2>/dev/null | wc -l) > 1 ))
then
tl=`grep -v ">" Query_exon_combined_exonerate.out.fa_orf.fa | wc | awk '{print$3-$1}'`
sp=`tail -1 Query_exon_combined_exonerate_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$tl" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.1f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
elif (( $(grep "[A-Z]" gblast_edited_consensus_orf.fa 2>/dev/null | wc -l) > 1 ))
then
tl=`grep -v ">" gblast_edited_consensus_orf.fa | wc | awk '{print$3-$1}'`
sp=`tail -1 gblast_edited_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$tl" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.1f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
else
tl=`grep -v ">" gblast_auto_consensus_orf.fa | wc | awk '{print$3-$1}'`
sp=`tail -1 gblast_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$tl" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.1f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
fi
#check ORF
if [[ -s spaln_exon_wise.fa ]]
then
co=`grep -v ">" spaln_exon_wise_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $co == "" ]] ; then orf="❌";else orf="✅"; fi
elif [[ -f exoneredit_consensus_orf.fa ]]
then
co=`grep -v ">" exoneredit_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $co == "" ]] ; then orf="❌";else orf="✅"; fi
elif [[ -f Query_exon_combined_exonerate.out.fa_orf.fa ]]
then
co=`grep -v ">" Query_exon_combined_exonerate.out.fa_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $co == "" ]] ; then orf="❌";else orf="✅"; fi
elif [[ -f gblast_edited_consensus_orf.fa ]]
then
co=`grep -v ">" gblast_edited_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $co == "" ]] ; then orf="❌";else orf="✅"; fi
else
co=`grep -v ">" gblast_auto_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $co == "" ]] ; then orf="❌";else orf="✅"; fi
fi
echo $j $gr $q $nht $nh1 $nh2 $nh3 $ewah $sc $a3 $a4 $de $pg $psg $m $n $gp $ons $a1 $a2 $ncon $pcon $qc $sp $orf
unset j gr q1 q nht nh1 nh2 nh3 ewah sc d a3 e a4 de pg psg m n gp ons splice_sites a1t a1 a2t a2 ncon pcon ins del fs tl ql qr qc sp co orf
cd /home/neo/bird_db1/aswin/APOBEC1
done < <(grep -if 2nd_gblast_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths) | sed '1i Species Group Query No_hits_tblastx No_hits_gblast No_hits_Exonerate No_hits_spaln #_Almost_hits START Loc Strand Exon_dups Paralogs Pseudogenes Ins Del Gaps(exon) Anonyms(exon) SS1_disrupting_exon SS2_disrupting_exon Cons Protein_cons Query_cov 1st_STOP_in_gene ORF' | column -t > 2nd_gblast_summary

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#10.4. Add completely obtained ORFs to the total validated queries list, gene tree & add their sequences to the exon_wise folder

#List of queries that gave complete ORF with 100% query coverage & no gaps from 2nd gblast

#If there are fragmented assemblies but with intact gene scattered in different scaffols, then location & strandedness restriction can be relieved
#awk '$5=="-" && $10"same" && $11="same" && $17"-" && $18=="-" && $19=="-" && $20=="-" && $23~"100" && $NF=="✅" {print$1}' 2nd_gblast_summary > 2nd_gblast_validated_queries
awk '$5=="-" && $17"-" && $18=="-" && $19=="-" && $20=="-" && $23~"100" && $NF=="✅" {print$1}' 2nd_gblast_summary > 2nd_gblast_validated_queries
cat 2nd_gblast_validated_queries | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' >> total_validated_queries

#Add to gene tree
for i in `cat 2nd_gblast_validated_queries`; do sed "/$i/ s/$/"2nd_gblast_"$i/g" -i gene_tree; unset j;done
#add sequence to folder
while read i
do
j=`echo $i | awk '{print$1}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/2nd_gblast
cat final_consensus.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/"APOBEC1_"$j.fa
cat extracted_flanking_region.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/"APOBEC1_"$j"_flanking".fa
unset j
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if 2nd_gblast_validated_queries /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#10.5. QC validated sequences

#Species with focal gene fragmented into different scaffolds or chromosomes
cut -f2- -d "_" total_validated_queries | sed 's/.fa//g' | xargs -n1 sh -c 'grep "^$0" gblast_principal_components' | column -t | awk '/ diff / {print$1}' > gene_breaking_assemblies

./qc_validated.sh validated_sequences validated_sequences_with_flanking_regions *_coding_exons_set APOBEC1 /home/neo/bird_db1/aswin/APOBEC1 2nd_gblast_validated_qc

#Find potential false positives in validated sequences by detecting outliers in nucleotide MSA
sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" 2nd_gblast_validated_qc | awk -F ":" '/Outlier/ {print$2}' | sed 's/Nothing//g' | tr "," "\n" | awk NF | sed 's/^[ ]\+//g' > outliers

#Remove outliers from the validated sequence records
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers 2nd_gblast_validated_queries
cat outliers | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' > outliers2
if [[ -s outliers2 ]]; then 
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers2 total_validated_queries
else :; fi
cat outliers | xargs -n1 sh -c 'sed "/$0/ s/[ ]\+[^ ]\+$//g" gene_tree -i'
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa) | xargs rm
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/*.fa) | xargs rm
rm outliers outliers2

##################################################################################################################################################################################################################################################################################################################
#11. SRA Blast

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#11.1.Filter out subjects that didn't form ORF

grep -if <(cat total_validated_queries | cut -f2,3- -d "_" | sed 's/.fa//g') gblast_subjects -v > sblast_subjects

#Find the closest query which gives best hits for sblast subjects

#Create all possbible pairwise distances between all species pairs - Create once & keep it in a path don't need eveytime you run the pipeline
while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
ld=`grep -if <(sed 's/.fa//g' total_validated_queries | cut -f2- -d "_") <(awk -v s="$j" '$2==s' /home/neo/bird_db1/aswin/tree/species_with_genome_data/all_possible_pairwise_distances | sort -k3 -n) | awk 'NR==1{print$3}'`
grep -if <(sed 's/.fa//g' total_validated_queries | cut -f2- -d "_") <(awk -v s="$j" '$2==s' /home/neo/bird_db1/aswin/tree/species_with_genome_data/all_possible_pairwise_distances | sort -k3 -n) | grep "$ld" | awk '{print$1}' > closest_query_list
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
#Loop through all the closest queries (including flanking region) & run blastn to find which gives the best hits in the subject
for x in `grep -if /home/neo/bird_db1/aswin/APOBEC1/closest_query_list <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/*.fa)`
do
blastn -task blastn -evalue 0.01 -db ../../genome/GC*.fna -query $x -num_threads 4 -outfmt "6 qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand qseq sseq"| sed '1i Query\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\tQuery_sequence\tSubject_sequence\n' > 2nd_test.out
gblast_best_hits $x 2nd_test.out 2nd_best_hits
q=`echo $x | awk -F "/" '{print$NF}'`
a=`awk 'NR>1{print$1}' 2nd_best_hits | wc -l`
b=`awk 'NR>1{print$2}' 2nd_best_hits | sort -u | wc -l`
c=`awk 'NR>1{print$3}' 2nd_best_hits | awk '{s+=$1}END{print s}'`
d=`awk 'NR>1{print$4}' 2nd_best_hits | awk '{s+=$1}END{print s}'`
e=`awk 'NR>1{print$10}' 2nd_best_hits | awk '{s+=$1}END{print s}'`
f=`awk 'NR>1{print$13}' 2nd_best_hits | awk '{s+=$1}END{print s}'`
g=`awk 'NR>1{print$14}' 2nd_best_hits | awk '{s+=$1}END{print s}'`
h=`awk 'NR>1{print$NF}' 2nd_best_hits | sort -u | wc -l`
echo $q $a $b $c $d $e $f $g $h
unset q a b c d e f g h
rm 2nd_test.out 2nd_best_hits
done | sed '1i Query No_exons No_locations Query_length Subject_length Bit_score %_Query_cov %_identity No_strands' | column -t > sblast_query_search_table
k=`cat sblast_query_search_table | grep -v "No_exons" | sort -k2,2 -k7,7 -nr | awk 'NR==1 {print$1}' | sed 's/APOBEC1_//g' | sed 's/_flanking.fa//g'`
echo $k $j
unset j ld x k
cd ~/bird_db1/aswin/APOBEC1
rm closest_query_list
done < <(grep -if sblast_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths) | column -t > sblast_query_subject_set

#Compare query-subject sets
join -1 1 -2 1 <(awk '{print$2,$1}' 2nd_gblast_query_subject_set_tmp | sort -k1) <(awk '{print$2,$1}' sblast_query_subject_set | sort -k1) | cat -n | awk '{if($3==$4) print$0,"same";else print$0,"diff"}' | column -t | grep ".*diff\|$" --color=always | GREP_COLORS="mt=32" grep ".*same\|$" --color=always > compare_2nd_gblast_sblast_query_subject_sets

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#11.2. Run SRA blast

while read z
do
y=`echo $z | awk '{print$1}' | cut -f1,2 -d "_"`
echo " ▶ " $y
cd $(echo $z | awk '{print$2}')/aswin/APOBEC1/
#Get the subject
d=`find ../../ -name "[SED]RR*.fa"`
#Find the most suitable(nearest) query species
q=`grep $y /home/neo/bird_db1/aswin/APOBEC1/sblast_query_subject_set | awk '{print$1}'`
#Get the query sequence
pqs=`find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1 -name "*$q*.fa" | awk -F "/" '{print$NF}'`
pqsf=`find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/ -maxdepth 1 -name "*$q*.fa" | awk -F "/" '{print$NF}'`
#Create optimal query with the combination of validated queries & own consensus sequence from Spaln/Gedit/Exonerate from 1st or 2nd round of genome search
#Exons to take from validated queries and own gedit consensus sequence : criterias :- Gaps, Anonyms, Splice site disruptions, Exon boundary indels, Continous indels longer than cerstain value & their range, % Query cover (not splice site should be in caps)
#---------------------------------------------------#---------------------------------------------------#---------------------------------------------------#---------------------------------------------------#---------------------------------------------------

#Create optimal query for sblast
if [[ -d "2nd_gblast" ]]
then
echo "   ↳ Taking sblast optimal query from 2nd Gblast"
cd 2nd_gblast
#Set the query for spaln (pq1) & gedit (pq)
pq1=$(grep $y /home/neo/bird_db1/aswin/APOBEC1/2nd_gblast_query_subject_set 2>/dev/null | awk '{print$1}' | xargs -I {} sh -c 'ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa 2>/dev/null | grep {} | awk -F "/" "{print\$NF}"')
#pq1=`ls APOBEC1_* | grep -v "flanking"`
pq=`echo $pq1`
#Check if the exonerate gave hits for all the exons, if exonerate was executed
enh=`grep ">" $pq1 -c | xargs -n1 bash -c 'expr $0 - $(grep ">" spaln_exon_wise.fa -c 2>/dev/null) 2>/dev/null'`
#Query for sblast
find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1  -name "*$q*.fa" -exec cp {} . \;
find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/ -maxdepth 1  -name "*$q*.fa" -exec cp {} . \;
else
echo "   ↳ Taking sblast optimal query from 1st Gblast"
#Set the query for spaln (pq1) & gedit (pq)
pq1=$(grep $y /home/neo/bird_db1/aswin/APOBEC1/spaln_query_subject_set 2>/dev/null | awk '{print$1}' | xargs -I {} sh -c 'ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa 2>/dev/null | grep {} | awk -F "/" "{print\$NF}"')
if [[ $pq1 == "" ]]; then enh="-"; else
enh=`grep ">" $pq1 -c | xargs -n1 bash -c 'expr $0 - $(grep ">" spaln_exon_wise.fa -c 2>/dev/null) 2>/dev/null'`
fi
pq=$(grep $y /home/neo/bird_db1/aswin/APOBEC1/gblast_query_subject_set 2>/dev/null | awk '{print$1}' | xargs -I {} sh -c 'ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa 2>/dev/null | grep {} | awk -F "/" "{print\$NF}"')
#Query for sblast
find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/ -maxdepth 1  -name "*$q*.fa" -exec cp {} . \;
find /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/ -maxdepth 1  -name "*$q*.fa" -exec cp {} . \;
fi

if [[ -s spaln_exon_wise.fa ]]
then
echo "      ↳ Both Spaln & Gedit are used to make sblast optimal query "
#-------------------------------------------------------------------------

#Spaln exon-wise summary for making sblast optimal query
while mapfile -t -n 2 ary && ((${#ary[@]}))
do
printf '%s\n' "${ary[@]}" > tmp
id=`grep ">" tmp | tr -d ">"`
exn=`grep ">" tmp | sed 's/>exon_//g'`
gaps=`grep -v ">" tmp | grep "N" -o | wc -l`
anonyms=`grep -v ">" tmp | tr -d '[ATGCNatgcn]' | awk NF | wc | awk '{print$3-$1}'`
#Splice sites : For 1st and last exon use only one splice site
last_exon=`grep ">" $pq1 -c`
if [[ $exn == "1" ]]
then
ss1d="-"
ss2=`grep "$id" -w <(cat spaln_splice_sites | tr " " "\n" | sed '1i -' | sed '$a -' | paste -d " " - - | awk '{print"exon_"++i,$0}') | awk '$3=="GT" || $3=="GC" {print$1}'`
if [[ $ss2 == "" ]]; then ss2d="x"; else ss2d="✔"; fi
elif [[ $exn == $last_exon ]]
then
ss1=`grep "$id" -w <(cat spaln_splice_sites | tr " " "\n" | sed '1i -' | sed '$a -' | paste -d " " - - | awk '{print"exon_"++i,$0}') | awk '$2=="AG" {print$1}'`
if [[ $ss1 == "" ]]; then ss1d="x"; else ss1d="✔"; fi
ss2d="-"
else
ss1=`grep "$id" -w <(cat spaln_splice_sites | tr " " "\n" | sed '1i -' | sed '$a -' | paste -d " " - - | awk '{print"exon_"++i,$0}') | awk '$2=="AG" {print$1}'`
if [[ $ss1 == "" ]]; then ss1d="x"; else ss1d="✔"; fi
ss2=`grep "$id" -w <(cat spaln_splice_sites | tr " " "\n" | sed '1i -' | sed '$a -' | paste -d " " - - | awk '{print"exon_"++i,$0}') | awk '$3=="GT" || $3=="GC" {print$1}'`
if [[ $ss2 == "" ]]; then ss2d="x"; else ss2d="✔"; fi
fi
#Exon boundary indels
ebi1=`grep "$id" -w spaln_boundary_indels_wrt_query | awk '{print$2 + $4}'`
ebi2=`grep "$id" -w spaln_boundary_indels_wrt_query | awk '{print$3 + $5}'`
#Longest continous indels length & range
lci=`grep "$id" -w spaln_boundary_indels_wrt_query | awk '{print$7}'`
lcir=`grep "$id" -w spaln_boundary_indels_wrt_query | awk '{print$8}'`
lcd=`grep "$id" -w spaln_boundary_indels_wrt_query | awk '{print$10}'`
lcdr=`grep "$id" -w spaln_boundary_indels_wrt_query | awk '{print$11}'`
#Total exon insertions
needle <(seqtk subseq $pq1 <(echo $id) | grep -v ">" | sed '1i >query') <(seqtk subseq spaln_exon_wise.fa <(echo $id) | grep -v ">" | sed '1i >spaln') -auto -stdout -awidth 10000 > stmp
tei=`grep "^query" stmp | awk '{print$3}' | sed 's/-\+$//g' | sed 's/^-\+//g' | tr -cd "-" | wc | awk '{print$3-$1}'`
#Total exon deletions
ted=`grep "^spaln" stmp | awk '{print$3}' | sed 's/-\+$//g' | sed 's/^-\+//g' | tr -cd "-" | wc | awk '{print$3-$1}'`
#Total exon query covered
teqc=$(cat stmp | awk '/^query/ {print$3}' | sed 's/-\+$//g' | sed 's/^-\+//g' | wc | awk '{print$3-$1}' | xargs -I {} sh -c 'grep "|" stmp | egrep -aob "\||\." | awk -F ":" "NR==1{print\$1} END{print\$1}" | paste -s -d " " | awk "{print(\$2-\$1+1)/{}*100}" | awk "{\$1+=0}1" CONVFMT="%.2f"')
echo $id $gaps $anonyms $ss1d $ss2d $ebi1 $ebi2 $lci $lcir $lcd $lcdr $tei $ted $teqc
unset id exn gaps anonyms last_exon ss1 ss1d ss2 ss2d ebi1 ebi2 lci lcir lcd lcdr tei ted teqc
rm tmp stmp
done < <(sed '/>/! s/[a-z]/\U&/g' spaln_exon_wise.fa) | sed '1i Exon Gaps Anonyms SS1 SS2 Boundary_indel_1 Boundary_indel_2 Longest_indel1 Longest_Indel1_range Longest_Indel2 Longest_Indel2_range Insertions Deletions %_query_cov' | column -t > spaln_optimal_query_summary
unset pq1 enh
else
echo "      ↳ Create optimal sblast query from Gedit since Spaln is either not executed or didn't gave hits for all exons!"
fi
unset pq1 enh
#--------------------------------------------------------

#Gedit exon-wise summary for making sblast optimal query
#For exons with no hits
exwn=`grep "No hits found" test.outfmt3 -B5 | awk '/Query=/ {print$NF}'`
if [[ ! -z $exwn && ! -z ${exwn+x} ]]
then
echo "         ↳ Gedit : Some exons didn't gave hits ! "
grep "No hits found" test.outfmt3 -B5 | awk '/Query=/ {print$NF, "- - - - - - - - - - - - -"}' > stmp1
else
echo "         ↳ Gedit : All exons gave hits "
fi
unset exwn
while mapfile -t -n 2 ary && ((${#ary[@]}))
do
printf '%s\n' "${ary[@]}" > tmp
id=`grep ">" tmp | tr -d ">"`
exn=`grep ">" tmp | sed 's/>exon_//g'`
gaps=`grep -v ">" tmp | grep "N" -o | wc -l`
anonyms=`grep -v ">" tmp | tr -d '[ATGCNatgcn]' | awk NF | wc | awk '{print$3-$1}'`
#Splice sites : For 1st and last exon use only one splice site
last_exon=`grep ">" $pq -c`
if [[ $exn == "1" ]]
then
ss1d="-"
ss2=`grep "$exn" -w <(cat gedit_splice_sites -n) | awk '$3=="GT" || $3=="GC" {print$1}'`
if [[ $ss2 == "" ]]; then ss2d="x"; else ss2d="✔"; fi
elif [[ $exn == $last_exon ]]
then
ss1=`grep "$exn" -w <(cat gedit_splice_sites -n) | awk '$2=="AG" {print$1}'`
if [[ $ss1 == "" ]]; then ss1d="x"; else ss1d="✔"; fi
ss2d="-"
else
ss1=`grep "$exn" -w <(cat gedit_splice_sites -n) | awk '$2=="AG" {print$1}'`
if [[ $ss1 == "" ]]; then ss1d="x"; else ss1d="✔"; fi
ss2=`grep "$exn" -w <(cat gedit_splice_sites -n) | awk '$3=="GT" || $3=="GC" {print$1}'`
if [[ $ss2 == "" ]]; then ss2d="x"; else ss2d="✔"; fi
fi
#Exon boundary indels
ebi1=`grep "exon_$exn" -w gedit_exon_boundary_indels | awk '{print$2 + $4}'`
ebi2=`grep "exon_$exn" -w gedit_exon_boundary_indels | awk '{print$3 + $5}'`
#Longest continous indels length & range
lci=`grep "exon_$exn" -w gedit_exon_boundary_indels | awk '{print$7}'`
lcir=`grep "exon_$exn" -w gedit_exon_boundary_indels | awk '{print$8}'`
lcd=`grep "exon_$exn" -w gedit_exon_boundary_indels | awk '{print$10}'`
lcdr=`grep "exon_$exn" -w gedit_exon_boundary_indels | awk '{print$11}'`
#Total exon insertions
tei=`grep "exon_$exn" -w gblast_edited_insertions | awk '{print$2}'`
#Total exon deletions
ted=`grep "exon_$exn" -w gblast_edited_deletions | awk '{print$2}'`
#Total exon query covered
teqc=`grep "exon_$exn" -w gblast_edited_query_covered | awk '{print$4}'`
echo $id $gaps $anonyms $ss1d $ss2d $ebi1 $ebi2 $lci $lcir $lcd $lcdr $tei $ted $teqc
unset id exn gaps anonyms last_exon ss1 ss1d ss2 ss2d ebi1 ebi2 lci lcir lcd lcdr tei ted teqc
rm tmp
done < <(sed '/>/! s/[a-z]/\U&/g' gblast_edited_consensus.fa) > stmp2
#Combine data
cat stmp1 stmp2 2>/dev/null | sort -k1,1 -V | sed '1i Exon Gaps Anonyms SS1 SS2 Boundary_indel_1 Boundary_indel_2 Longest_indel1 Longest_Indel1_range Longest_Indel2 Longest_Indel2_range Insertions Deletions %_query_cov' | column -t > gedit_optimal_query_summary
rm stmp1 stmp2 2>/dev/null
unset pq
#Create lists of exons obtained from the subject species genome that are qualified to be sblast query
if [[ -s spaln_optimal_query_summary ]]
then
awk '$2=="0" && $3=="0" && $4!="x" && $5!="x" && $6<4 && $7<4 && ($8<10 && $9!~"-9[0-9]") && ($10<10 && $11!~"-9[0-9]") && $NF~"100" {print$1}' spaln_optimal_query_summary > spaln_complete_exons
awk '$2=="0" && $3=="0" && $4!="x" && $5!="x" && $6<4 && $7<4 && ($8<10 && $9!~"-9[0-9]") && ($10<10 && $11!~"-9[0-9]") && $NF~"100" {print$1}' gedit_optimal_query_summary > gedit_complete_exons
#Extract sequence for sblast query
seqtk subseq spaln_exon_wise.fa spaln_complete_exons > sblast_query.fa_tmp
seqtk subseq extracted_flanking_region.fa spaln_complete_exons > sblast_query_flanking.fa_tmp
seqtk subseq gblast_edited_consensus.fa <(grep -if spaln_complete_exons gedit_complete_exons -w -v) >> sblast_query.fa_tmp
seqtk subseq extracted_flanking_region.fa <(grep -if spaln_complete_exons gedit_complete_exons -w -v) >> sblast_query_flanking.fa_tmp
else
awk '$2=="0" && $3=="0" && $4!="x" && $5!="x" && $6<4 && $7<4 && ($8<10 && $9!~"-9[0-9]") && ($10<10 && $11!~"-9[0-9]") && $NF~"100" {print$1}' gedit_optimal_query_summary > gedit_complete_exons
#Extract sequence for sblast query
seqtk subseq gblast_edited_consensus.fa gedit_complete_exons > sblast_query.fa_tmp
seqtk subseq extracted_flanking_region.fa gedit_complete_exons > sblast_query_flanking.fa_tmp
fi
grep -if <(grep ">" sblast_query.fa_tmp) <(grep ">" $pqs) -w -v | tr -d ">" > sblast_incomplete_exons
#Choose query sequence with flanking regions for sblast incomplete exons in case the query gene have very short exons or high diveregnce (customizable)
seqtk subseq $pqs sblast_incomplete_exons >> sblast_query.fa_tmp
seqtk subseq $pqsf sblast_incomplete_exons >> sblast_query_flanking.fa_tmp
#Sort query based on exon numberings
sed 's/^>/\x00&/' sblast_query.fa_tmp | sort -z -V | tr -d '\0' | sed '/>/! s/[a-z]/\U&/g' > sblast_query.fa
sed 's/^>/\x00&/' sblast_query_flanking.fa_tmp | sort -z -V | tr -d '\0' | sed '/>/! s/[a-z]/\U&/g' > sblast_query_flanking.fa
rm sblast_query.fa_tmp sblast_query_flanking.fa_tmp
#End of sblast query making
if [[ $(pwd | awk -F "/" '{print$NF}') == "2nd_gblast" ]]; then cp sblast_query.fa sblast_query_flanking.fa ../ ; cd ../; else :; fi
#---------------------------------------------------#---------------------------------------------------#---------------------------------------------------#---------------------------------------------------#---------------------------------------------------

#SRA blast (since the query is optimised based on the 'nearest validated query' & 'sequences from own genome', The minimum evalue for good homology matches is "0.01" & to get more accurate consensus with minimum variants in SRA blast output)
echo "   ↳ Running SRA BLAST ... "
echo $d | awk -F "/" '{print$NF}' | sed 's/^/      ▬  Subject: /g'
if [[ -s sblast_incomplete_exons ]]
then
echo "      ▬  Query: " | sed "s/$/$pqs($(find . -name "sblast_incomplete_exons" -exec sh -c 'cat {}' \; | paste -s -d " " | sed 's/ exon_/,/g')) \& ($(find . -name "*_complete_exons" -exec sh -c 'cat {}' \; | sort -uV | paste -s -d " " | sed 's/ exon_/,/g'))/g"
else
echo "      ▬  Query: " | sed "s/$/Extracted from subject Genome ($(find . -name "*_complete_exons" -exec sh -c 'cat {}' \; | sort -uV | paste -s -d " " | sed 's/ exon_/,/g'))/g"
fi
#Blast query against SRA : generate output in archive form s.t it can be reformatted to any other formats
#Run blast with query sequence with flanking regions in case the query gene have very short exons or high diveregnce (customizable)
blastn -task blastn -db $d -query sblast_query_flanking.fa -evalue 0.01 -outfmt 11 -out optimal_query_sblast.outfmt11
#Reformat to outfmt 3
blast_formatter -archive optimal_query_sblast.outfmt11 -outfmt 3 -line_length 280 -out optimal_query_sblast.outfmt3
#Reformat to sam format
blast_formatter -archive optimal_query_sblast.outfmt11 -outfmt "17 SQ" -out optimal_query_sblast.sam
#SRA blast summary
qvsra optimal_query_sblast.outfmt3 -t > optimal_query_sblast_summary
#SRA blast QC & consensus sequence making using 2 iffernt tools : to compare & quality check
sedit.sh optimal_query_sblast.outfmt3 2>/dev/null
samtools consensus optimal_query_sblast.sam -f fasta --mode simple -l 10000 | sed 's/Query/exon/g' > optimal_query_sblast.sam.fa

#QC consenus : pairwise align 2 consensus sequence
needle <(grep -v ">" optimal_query_sblast.outfmt3.fa) <(grep -v ">" optimal_query_sblast.sam.fa) -auto -stdout -awidth 280 > pairwise_sblast_consensus_seqences.aln

#---------------------------------------------------#---------------------------------------------------#---------------------------------------------------#---------------------------------------------------#---------------------------------------------------
#Create final filtered consensus sequence based on pairwise alignment between sblast hit & query exon-wise CDS sequence

#Create pairwise alignments of sblast hit & query CDS sequence
#s1=`find . -name "$pqs" | sed 's/^.\///g' | tail -1`
exalign sblast_query.fa optimal_query_sblast.outfmt3.fa 2>/dev/null | sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" | awk '/optimal/ {print">"$0}' RS=">" | awk NF | awk '/^>/ {OUT="sblast_pairwise_"substr($0,2) ".temp"}; OUT{print > OUT}'
for p in `ls sblast_pairwise_exon_*.temp | sort -V`
do
p2=`echo $p | sed 's/temp/aln/g'`
tc=`egrep -v ">|\|" $p | awk NF | grep "^ \+" -o | awk '{print length($0)}' | sort -u | awk '{print$0+1}'`
grep -v ">" $p | cut -c $tc- | awk NF > $p2
unset p2 tc
rm $p
done
#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#

#If the query has only one exon then don't need to look for splice sites, just pairwise align & extract the aligned region
note=`grep ">" sblast_query.fa -c`
if [[ $note == 1 ]]
then
start=`cat sblast_pairwise_exon_1.aln | grep -B2 "^optimal" | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":" | awk '{print$0+1}'`
#Extract the region from which the pairwise alignment start i.e. remove the regions before alignment starts (usually alignment starts with START codon ATG)
cat sblast_pairwise_exon_1.aln | grep -B2 "^optimal" | grep "|" -C1 | cut -c $start- > tmp1
end=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+1}'`
ln=`seqtk subseq sblast_query.fa <(echo exon_1) | grep -v ">" | awk NF | wc | awk '{print$3-$1}'`
#Extract some stats from pairwise alignment
#For each exon : Exon name, Total query length, Total query covered, Percent query covered
cat tmp1 | cut -c -$end | head -1 | tr -d "-" | wc | awk '{print$3-$1}' | awk -v p="$ln" '{print p,$1,($1/p)*100}' | sed "s/^/exon_1 /g" > sblast_edited_query_covered
#insertions
cat tmp1 | cut -c -$end | head -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exon_1 /g" > sblast_edited_insertions
#deletions
cat tmp1 | cut -c -$end | tail -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exon_1 /g" > sblast_edited_deletions
cat <(echo ">exon_1") <(cat tmp1 | cut -c -$end | tail -1) | tr -d "-" > sblast_edited_consensus.fa
rm tmp1
unset start end ln
else
#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#

#If query has more than 1 exon & exon 1 gave blast hit; then look for splice site "GT"
#Exclude 1st and last exon from the loop to find 2 splice sites since they have only one splice site
#For the first exon : extract from start codon
if [[ -e sblast_pairwise_exon_1.aln ]]
then
start=`cat sblast_pairwise_exon_1.aln | grep -B2 "^optimal" | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":" | awk '{print$0+1}'`
#Extract the region from which the pairwise alignment start i.e. remove the regions before alignment starts (usually alignment starts with START codon ATG)
cat sblast_pairwise_exon_1.aln | grep -B2 "^optimal" | grep "|" -C1 | cut -c $start- > tmp1
#Look for nearest splice site "GT/GC" within next 5 flanking nucleotides after the query exon hit ended (the threshold flanking bases & other variable splice sites are customizable)
ssc1=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+2}'`
ssc2=`cat tmp1 | grep "|" -C1 | cut -c $ssc1- | tail -1 | cut -c -5 | grep -i "GT" -aob | head -1 | cut -f1 -d ":"`
#Allow splice site variants
if [[ $ssc2 == "" ]]; then ssc2=`cat tmp1 | grep "|" -C1 | cut -c $ssc1- | tail -1 | cut -c -5 | grep -i "GC" -aob | head -1 | cut -f1 -d ":"`; else :; fi
ssc3=`expr $ssc1 + $ssc2 2>/dev/null | awk '{print$0-1}'`
if [[ $ssc3 == "" ]]; then er=`echo "     ↳ exon_1 - No 5 prime splice site found! "`; ssc3=`echo $ssc1 | awk '{print$0-1}'`; else :; fi
unset ln
ln=`seqtk subseq sblast_query.fa <(echo "exon_1") | grep -v ">" | awk NF | wc | awk '{print$3-$1}'`
#Extract some stats from pairwise alignment
#For each exon : Exon name, Total query length, Total query covered, Percent query covered
cat tmp1 | cut -c -$ssc3 | head -1 | tr -d "-" | wc | awk '{print$3-$1}' | awk -v p="$ln" '{print p,$1,($1/p)*100}' | sed "s/^/exon_1 /g" > sblast_edited_query_covered
#insertions
cat tmp1 | cut -c -$ssc3 | head -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exon_1 /g" > sblast_edited_insertions
#deletions
cat tmp1 | cut -c -$ssc3 | tail -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/exon_1 /g" > sblast_edited_deletions
cat <(echo ">exon_1") <(cat tmp1 | cut -c -$ssc3 | tail -1) | tr -d "-" > sblast_edited_consensus.fa
if [ -n "$er" ] ; then echo "$er"; else :; fi
rm tmp1
unset start ssc1 ssc2 ssc3 ln l p er
else
echo "     ↳ exon 1 alignment not found"
#Create an empty file for making sblast edited output files because the first time they are created it should be direct redirection(>) not appending (>>)
echo -n > sblast_edited_query_covered
echo -n > sblast_edited_insertions
echo -n > sblast_edited_deletions
echo -n > sblast_edited_consensus.fa
fi
#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#

last_exon=`grep ">" sblast_query.fa | tail -1 | tr -d ">"`

#Run if total number of exons are greater than 2
if [ $note -gt 2 ]
then
#Middle exons are flanked by splice sites on both sides unlike 1st & last exons
#Extract complete exons from subject based on the pairwise alignment between query sequence & subject genome with 15bp flanking sequence
for i in `ls sblast_pairwise_exon_* | sort -V | egrep -v "exon_1.aln|$last_exon.aln"`
do
#Exon number
j=`echo $i | awk -F "_" '{print$(NF-1)"_"$NF}' | sed 's/.aln//g'`
echo ">"$j
#Splice site 1 for the respective exon
ssc4=`cat $i | grep -B2 "^optimal" | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":"`
#Look for nearest splice site "AG" within next 5 flanking nucleotides after the query exon hit ended (the threshold flanking bases & other variable splice sites are customizable)
ssc5=`cat $i | grep -B2 "^optimal" | grep "|" -C1 | cut -c -$ssc4 | tail -1 | rev | cut -c -5 | grep -i "GA" -aob | head -1 | cut -f1 -d ":"`
ssc6=`expr $ssc4 - $ssc5 2>/dev/null | awk '{print$0+1}'`
if [[ $ssc6 == "" ]]; then echo "     ↳ "$j" - No 3 prime splice site found! " >> er1; ssc6=`echo $ssc4 | awk '{print$0+1}'`; else :; fi
cat $i | grep -B2 "^optimal" | grep "|" -C1 | cut -c $ssc6- > tmp1
ssc7=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+2}'`
#Look for nearest splice site "GT/GC" within next 5 flanking nucleotides after the query exon hit ended (the threshold flanking bases & other variable splice sites are customizable)
ssc8=`cat tmp1 | grep "|" -C1 | cut -c $ssc7- | tail -1 | cut -c -5 | grep -i "GT" -aob | head -1 | cut -f1 -d ":"`
#Allow splice site variants
if [[ $ssc8 == "" ]]; then ssc8=`cat tmp1 | grep "|" -C1 | cut -c $ssc7- | tail -1 | cut -c -5 | grep -i "GC" -aob | head -1 | cut -f1 -d ":"`; else :; fi
ssc9=`expr $ssc7 + $ssc8 2>/dev/null | awk '{print$0-1}'`
if [[ $ssc9 == "" ]]; then echo "     ↳ "$j" - No 5 prime splice site found! " >> er2; ssc9=`echo $ssc7 | awk '{print$0-1}'`; else :; fi
unset ln
ln=`seqtk subseq sblast_query.fa <(echo $j) | grep -v ">" | awk NF | wc | awk '{print$3-$1}'`
#Extract some stats from pairwise alignment
#For each exon : Exon name, Total query length, Total query covered, Percent query covered
cat tmp1 | cut -c -$ssc9 | head -1 | tr -d "-" | wc | awk '{print$3-$1}' | awk -v p="$ln" '{print p,$1,($1/p)*100}' | sed "s/^/$j /g" >> sblast_edited_query_covered
#insertions
cat tmp1 | cut -c -$ssc9 | head -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/$j /g" >> sblast_edited_insertions
#deletions
cat tmp1 | cut -c -$ssc9 | tail -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/$j /g" >> sblast_edited_deletions
cat tmp1 | cut -c -$ssc9 | tail -1
rm tmp1
unset j ssc4 ssc5 ssc6 ssc7 ssc8 ssc9 ln l p
done | tr -d "-" >> sblast_edited_consensus.fa
if [[ -e er1 ]] ; then cat er1; rm er1; else :; fi
if [[ -e er2 ]] ; then cat er2; rm er2; else :; fi
#If the number of exons are just 2 then don't run this,rather run the next code for last exon
else :
fi
#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#

#Extract last exon
if [[ -e $(ls sblast_pairwise_$last_exon.aln) ]] 2>/dev/null
then
ssc10=`ls sblast_pairwise_exon_* | grep $last_exon | xargs cat | grep -B2 "^optimal" | grep "|" -C1 | sed -n '2p' | egrep -aob "\||\." | head -1 | cut -f1 -d ":"`
ssc11=`ls sblast_pairwise_exon_* | grep $last_exon | xargs cat | grep -B2 "^optimal" | grep "|" -C1 | cut -c -$ssc10 | tail -1 | rev | cut -c -5 | grep -i "GA" -aob | head -1 | cut -f1 -d ":"`
ssc12=`expr $ssc10 - $ssc11 2>/dev/null | awk '{print$0+1}'`
if [[ $ssc12 == "" ]]; then er=`echo "     ↳ " $last_exon" - No 3 prime splice site found! "`; ssc12=`echo $ssc10 | awk '{print$0+1}'`; else :; fi
ls sblast_pairwise_exon_* | grep $last_exon | xargs cat | grep -B2 "^optimal" | grep "|" -C1 | cut -c $ssc12- > tmp1
end=`cat tmp1 | sed -n '2p' | egrep -aob "\||\." | tail -1 | cut -f1 -d ":" | awk '{print$0+1}'`
unset ln
ln=`seqtk subseq sblast_query.fa <(echo $last_exon) | grep -v ">" | awk NF | wc | awk '{print$3-$1}'`
cat tmp1 | cut -c -$end | head -1 | tr -d "-" | wc | awk '{print$3-$1}' | awk -v p="$ln" '{print p,$1,($1/p)*100}' | sed "s/^/$last_exon /g" >> sblast_edited_query_covered
#insertions
cat tmp1 | cut -c -$end | head -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/$last_exon /g" >> sblast_edited_insertions
#deletions
cat tmp1 | cut -c -$end | tail -1 | tr -cd "-" | wc | awk '{print$3-$1}' | sed "s/^/$last_exon /g" >> sblast_edited_deletions
cat <(echo ">"$last_exon) <(cat tmp1 | cut -c -$end | tail -1) | tr -d "-" >> sblast_edited_consensus.fa
rm tmp1
if [ -n "$er" ] ; then echo "$er"; else :; fi
unset ssc10 ssc11 ssc12 end ln l p last_exon er
else
echo "     ↳ "$last_exon" alignment not found"
fi
#End of sedit
fi
unset p note
#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#

#check splice sites
for spst in `grep ">" sblast_query.fa | tr -d ">"`
do
sscheck=`find . -name "sblast_incomplete_exons" -exec sh -c 'cat {}' \; | grep $spst -w`
#if the query was complete that means it already had correct splice sites
if [[ $sscheck == "" ]]
then
ss1="AG"
ss2="GT"
else
z1=`seqtk subseq sblast_edited_consensus.fa <(echo $spst) | grep -v ">"`
if [[ $z1 == "" ]]; then
ss1="-"
ss2="-"
else
ss1=`sed -n "s/$z1.*//p" optimal_query_sblast.outfmt3.fa | rev | cut -c -2 | tr [:lower:] [:upper:] | rev`
if [[ $ss1 == "" ]]; then ss1="-"; else :; fi
ss2=`sed -n "s/.*$z1//p" optimal_query_sblast.outfmt3.fa | cut -c -2 | tr [:lower:] [:upper:]`
if [[ $ss2 == "" ]]; then ss2="-"; else :; fi
fi
fi
echo $ss1 $ss2
unset z1 sscheck ss1 ss2
done > sedit_splice_sites
#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#------------------#

#Generate ORF
transeq <(grep -v ">" sblast_edited_consensus.fa | paste -s -d "" | sed "1i \>$y") sblast_edited_consensus_orf.fa 2>/dev/null
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" sblast_edited_consensus.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" sblast_edited_consensus_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > sblast_edited_translated 2>/dev/null
#For better visualization
gnaa=`head -1 sblast_edited_translated | wc | awk '{print$3-$1}'` 2>/dev/null
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" sblast_edited_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 sblast_edited_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > sblast_edited_nucleotide_amino_acid_view 2>/dev/null
unset y d q pqs pqsf s1 spst gnaa
cd /home/neo/bird_db1/aswin/APOBEC1/
printf '_%.s' {1..150}; printf "\n"
done < <(grep -if sblast_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#11.3. Sblast summary for each exon counts: For comprehensive checking

#Summary table of hits have information in exon_wise manner - hence for proper alignment of exon summary make different summary tables if the exon count of validated quries used for gblast are varying
#Split query-subject set based on query exon-count : NOTE : scientific names may have 2/3 words; ensure that there are no species missing due missing the 3rd word
#Since more queries got validated the APOBEC1_ortholog_info must be updated before making coding exon sets
for i in `find validated_sequences/ -name "*.fa"`; do j=`grep ">" $i -c`; k=`echo $i | cut -f3,4 -d "_" | sed 's/.fa//g'`; echo $k $j; unset j k; done | column -t > APOBEC1_ortholog_info_updated
cp sblast_query_subject_set coding_sets
awk '$2!=""' APOBEC1_ortholog_info_updated | awk '{print$1,$2}' | xargs -n2 sh -c 'sed "/^$0/ s/$/ $1/g" -i coding_sets'
rm *_coding_exons_sblast_set
cat coding_sets | column -t | awk '{print >> $3"_coding_exons_sblast_set"; close($3)}'
rm coding_sets

counter=0
#Summary for each sets
for set in `ls *_coding_exons_sblast_set`
do
counter=$((counter + 1))
exon_count=`echo $set | cut -f1 -d "_"`
echo -e
echo "  ▬  Generating summary tables of SRA Blast for" $exon_count "exons coding set"
./OQ_sblast_summary.sh APOBEC1 $exon_count $set /home/neo/bird_db1/aswin/APOBEC1/ > sblast_summary"$counter" 2>/dev/null
unset exon_count
echo -e
done
unset counter

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#11.4. Sblast Principal components(PC's) : Components from all data that contribuites the most to the interpretation/inference

#From Sblast : PC's of all species : choose PC's based on univeral headings (not dependent of characteristics of genes such as number of exons)
echo -n > sblast_principal_components
for i in `ls sblast_summary* | sort -V`
do
awk 'NR==1{for(i=1;i<=NF;i++)if($i~/Species|Group|SRA_size|Query_L|gedit_L|sblast_L|No_hits|#_Almost_hits|START|Depth_range|Poor_depth|Variants|Frameshifting_vars|Gaps|other_anonym|Ins|Del|%_Query_cov|1st_STOP_in_gene|ORF/)f[n++]=i}{for(i=0;i<n;i++)printf"%s%s",i?" ":"",$f[i];print""}' $i | column -t >> sblast_principal_components_tmp
done
awk '!a[$1]++' sblast_principal_components_tmp | column -t > sblast_principal_components
rm sblast_principal_components_tmp

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#11.5. Add completely obtained ORFs to the total validated queries list, gene tree & add their sequences to the exon_wise folder

#List of queries that gave complete ORF with 100% query coverage from sblast
awk '$7=="-" && $14=="0" && $19~"100" && $NF=="✅" {print$1}' sblast_principal_components > sblast_validated_queries
cat sblast_validated_queries | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' >> total_validated_queries
#Add to gene tree
for i in `cat sblast_validated_queries`; do sed "/$i/ s/$/"sblast_"$i/g" -i gene_tree; unset j;done
#add sequence to folder : sblast consensus sequences are not exon-wise
while read i
do
j=`echo $i | awk '{print$1}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/
cat sblast_edited_consensus.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/"APOBEC1_"$j.fa
cat optimal_query_sblast.outfmt3.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/"APOBEC1_"$j"_flanking".fa
unset j
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if sblast_validated_queries /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#11.6. QC validated sequences

#Species with focal gene fragmented into different scaffolds or chromosomes
cut -f2- -d "_" total_validated_queries | sed 's/.fa//g' | xargs -n1 sh -c 'grep "^$0" gblast_principal_components' | column -t | awk '/ diff / {print$1}' > gene_breaking_assemblies

./qc_validated.sh validated_sequences validated_sequences_with_flanking_regions *_coding_exons_set APOBEC1 /home/neo/bird_db1/aswin/APOBEC1 sblast_validated_qc

#Find potential false positives in validated sequences by detecting outliers in nucleotide MSA
sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" sblast_validated_qc | awk -F ":" '/Outlier/ {print$2}' | sed 's/Nothing//g' | tr "," "\n" | awk NF | sed 's/^[ ]\+//g' > outliers

#Remove outliers from the validated sequence records
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers sblast_validated_queries
cat outliers | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' > outliers2
if [[ -s outliers2 ]]; then 
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers2 total_validated_queries
else :; fi
cat outliers | xargs -n1 sh -c 'sed "/$0/ s/[ ]\+[^ ]\+$//g" gene_tree -i'
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa) | xargs rm &> /dev/null
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/*.fa) | xargs rm &> /dev/null
rm outliers outliers2

##################################################################################################################################################################################################################################################################################################################
#12. Sblast exonerate & spaln : Find exon boundaries in birds with lacking splice sites & complete ORF albeit having query cover >=90% in sblast

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#12.1. Subjects for Sblast exonerate & spaln
awk '$NF=="❌" && $7=="-" && $18>89 && $18<101 {print$1}' sblast_principal_components > sblast_exonerate_spaln_subjects

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#12.2. Run Sblast exonerate & spaln

while read z
do
species=`echo $z | awk '{print$1}' | cut -f1,2 -d "_"`
echo " ▶ " $species
cd $(echo $z | awk '{print$2}')/aswin/APOBEC1/
mkdir -p sblast_exonerate_spaln
cd sblast_exonerate_spaln

#Path of query sequence file i.e, custom made optimal sblast query (this query should not have flanking regions)
qp=`find ../ -name "sblast_query.fa" | head -1`

#Run Exonerate
echo "     ↳ Sedit didn't form complete ORF, running Exonerate to find correct exon boundaries ..."
#Run exonerate
#Translate to amino acid ("*" is required to represent stop codon for exonerate)
transeq <(grep -v ">" $qp) -auto -stdout | sed "/>/ s/.*/>sblast_query/g" | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > exonerate_query.aa 2>/dev/null
#Create subject : subject must be a single fasta file 
grep -v ">" ../optimal_query_sblast.outfmt3.fa | paste -s -d "" | sed 's/[^ATGCatgc]/N/g' | sed '1i >sblast_hit' > sblast_hit.fa
#Run exonerate in exhaustive mode  with model protein2genome:bestfit & maximum intron allowed is 36bp
exonerate -E --model p2g:b exonerate_query.aa sblast_hit.fa --ryo "Rank: %r\nTarget in alignment: %ti : %tal (%tab - %tae)\nTarget in coding sequence: %ti : %tcl (%tcb - %tce)\nGene orientation: %g\nPercent identity: %pi\nPercent similarity: %ps\n>%ti (%tab - %tae)\n%tcs\n" --querytype protein --targettype dna --showtargetgff yes --fsmmemory 28000 --alignmentwidth 270 > sblast_exonerate.out 2>/dev/null
#Create consensus sequence
#Exons_wise sequence
bedtools getfasta -fi sblast_hit.fa -bed <(awk '$3=="exon" {print"sblast_hit",$4-1,$5,$3 | "sort -k2 -n"}' sblast_exonerate.out | awk '{print$1,$2,$3,$4"_"++i}' OFS="\t") -name > sblast_exonerate.out_exon_wise.fa 2>/dev/null
grep -v ">" sblast_exonerate.out_exon_wise.fa | sed '1i >exonerate_consensus' > sblast_exonerate.out.fa
#Check ORF
transeq <(grep -v ">" sblast_exonerate.out.fa | paste -s -d "") sblast_exonerate.out.fa_orf.fa 2>/dev/null
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" sblast_exonerate.out.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" sblast_exonerate.out.fa_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > sblast_exonerate_translated
#For better visualization
gnaa=`head -1 sblast_exonerate_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" sblast_exonerate_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 sblast_exonerate_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > sblast_exonerate_nucleotide_amino_acid_view
#Create consensus from exonerate
cp sblast_exonerate.out_exon_wise.fa final_consensus.fa
unset gnaa
#Exonerate exon boundary indels : run only if exonerate gave hits on all queries
enh=`grep ">" $qp -c | xargs -n1 bash -c 'expr $0 - $(grep ">" sblast_exonerate.out_exon_wise.fa -c)'`
if [[ $enh == "0" ]]
then
while mapfile -t -n 4 ary && ((${#ary[@]}))
do 
ex=`printf '%s\n' "${ary[@]}" | head -1 | tr -d ">"`
printf '%s\n' "${ary[@]}" > tmp
cat tmp | grep "exoner" | grep [0-9] -aob | head -1 | cut -f1 -d ":" | awk '{print$1+3}' | xargs -n1 sh -c 'cat tmp | cut -c $0-' 2>/dev/null | awk NF | sed 's/ \+[0-9]\+//g' > tmp1
#insertions at exon start : counted if a dash is observed within the query sequence anywhere within the first 5% of pairwise alignment between query & subject
is=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"exoner",$0; else print"query",$0}' | awk '/query/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $is == "" ]]; then is="-"; else :; fi
#insertions at exon end : counted if a dash is observed within the query sequence anywhere after 95% of pairwise alignment between query & subject
ie=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"exoner",$0; else print"query",$0}' | awk '/query/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ie == "" ]]; then ie="-"; else :; fi
#deletions at exon start : counted if a dash is observed within the subject sequence anywhere within the first 10% (since subject sequence is longer due to flanking regions) of pairwise alignment between query & subject
ds=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"exoner",$0; else print"query",$0}' | awk '/exoner/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ds == "" ]]; then ds="-"; else :; fi
#deletions at exon end : counted if a dash is observed within the subject sequence anywhere after 90% of pairwise alignment between query & subject
de=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"exoner",$0; else print"query",$0}' | awk '/exoner/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
#Longest continous insertions : Position, Length & % of gene in which it appears
if [[ $(head -1 tmp1 | grep "\-") == "" ]]; then lci="- - -"; else
lci=`bedtools merge -i <(cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"exoner",$0; else print"query",$0}' | awk '/query/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
fi
#Longest continous deletions : Position, Length & % of gene in which it appears
if [[ $(tail -1 tmp1 | grep "\-") == "" ]]; then lcd="- - -"; else
lcd=`bedtools merge -i <(cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"exoner",$0; else print"query",$0}' | awk '/exoner/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
#lci=`cat tmp2 | grep "\-\+" -o | tr -d " " | awk '{print length}' | sort -nr | head -1`
fi
if [[ $de == "" ]]; then de="-"; else :; fi
echo $ex $is $ie $ds $de $lci $lcd
unset ex is ie ds de lci lcd
rm tmp tmp1
done < <(exalign2 $qp sblast_exonerate.out_exon_wise.fa 2>/dev/null | sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" | awk NF) | sed '1i Exons Ins_in_1st_5%_exon Ins_in_last_95%_exon Del_in_1st_5%_exon Del_in_last_95%_exon Longest_continous_ins_position(lci) lci_length lci_%_range Longest_continous_del_position(lcd) lcd_length lcd_%_range' | column -t > exonerate_exon_boundary_indels
else 
echo "     ↳ Exonerate didn't gave hit on all exons, terminating further search !"
fi

#Run spaln if required
o=`grep -v ">" sblast_exonerate.out.fa_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
splice_sites=`grep ">" $qp -c | awk '{print$0-1}'`
a1t=`cat sblast_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$1=="GT" || $1=="GC"' | wc -l`
if [[ $a1t == $splice_sites ]]; then a1="-"; else a1="x"; fi
a2t=`cat sblast_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$2=="AG"' | wc -l`
if [[ $a2t == $splice_sites ]]; then a2="-"; else a2="x"; fi

if [[ ((-z $o && ! -z ${o+x}) || $a1 != "-" || $a2 != "-") && $enh == "0" ]]
then
#Run Spaln
echo "     ↳ Exonerate/Exoneredit didn't form complete ORF, running Spaln ..."
#Translate to amino acid ("*" is required to represent stop codon for Spaln)
transeq <(grep -v ">" $qp) -auto -stdout | sed "/>/ s/.*/>sblast_query/g" | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > spaln_query.aa
#Query exon count
pqe1=`grep ">" $qp -c`
#Query cds length
pqe2=`grep -v ">" $qp | wc | awk '{print$3-$1}'`
#Splice sites expected in target
spn=`echo $pqe1 | awk '{print$1-1}'`
#Run Spaln with different combinations of -yX(intr/cross species parameters) -yS (species-specific parameter)
for x in 0 1
do
for y in 0 50 100
do
echo -n > spaln_out_tmp
echo -n > spaln_tmp
spaln -M1 -l280 -O4 -ya0 -yX$x -yS$y sblast_hit.fa spaln_query.aa > spaln_out_tmp 2>/dev/null
bedtools getfasta -fi sblast_hit.fa -bed <(awk '!/^#|^@/{print"sblast_hit",$9-1,$10,"exon_"++i}' OFS="\t" spaln_out_tmp) -name > spaln_tmp 2>/dev/null
x1=`grep ">" -c spaln_tmp`
x2=`grep -v ">" spaln_tmp | wc | awk '{print$3-$1}'`
x3=`needle <(grep -v ">" spaln_tmp) <(grep -v ">" $qp) -auto -stdout -awidth 280 | awk '/Identity|Similarity|Gaps|Score/ {print$NF}' | tr -d ')(%' | paste -s -d " "` 2>/dev/null
transeq <(grep -v ">" spaln_tmp | paste -s -d "" | sed "1i \>$species") spaln_tmp_orf.fa 2>/dev/null
#Splice sites
sp1=`awk '!/^#|^@/ {print$NF}' spaln_out_tmp | tr "." " " | awk NF | tr '[:lower:]' '[:upper:]' | awk '$1=="GT" || $1=="GC"' | wc -l`
if [[ $sp1 == $spn ]]; then sp1s="✔"; else sp1s="x"; fi
sp2=`awk '!/^#|^@/ {print$NF}' spaln_out_tmp | tr "." " " | awk NF | tr '[:lower:]' '[:upper:]' | awk '$2=="AG"' | wc -l`
if [[ $sp2 == $spn ]]; then sp2s="✔"; else sp2s="x"; fi
#Check ORF presence absence
x4=`grep -v ">" spaln_tmp_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $x4 == "" ]] ; then x5="❌";else x5="✅"; fi
echo $x $y $pqe1 $x1 $pqe2 $x2 $x3 $sp1 $sp1s $sp2 $sp2s $x5
unset x1 x2 x3 x4 sp1 sp1s sp2 sp2s x5
rm spaln_out_tmp spaln_tmp spaln_tmp_orf.fa
done
done | sed '1i yX yS Query_exons Spaln_exons Query_length Spaln_length %_Identity %_Similarity %_Gaps Score No_Splice_site1 SS1_disruption No_Splice_site2 SS2_disruption ORF' | column -t > spaln_summary
unset x y

#Choose best parameters based on exon count match & alignment stats
spaln_params=`cat <(awk '$3==$4 && $12=="✔" && $14=="✔"' spaln_summary | grep "✅") <(awk '$3==$4' spaln_summary | grep "✅") <(awk '$3==$4' spaln_summary | grep "❌") | awk '!a[$0]++' | awk 'NR==1{print$1,$2}'`
if [[ $spaln_params == "" ]]
then 
unset spaln_params
spaln_params=`grep -v "yX" spaln_summary | sort -k4,4 -nr | awk 'NR==1{print$1,$2}'`
echo "        ↳ No ideal Spaln parameters found! Taking the parameters that gave the highest number of exon hits : " $spaln_params 
else 
echo "        ↳ Spaln parameters used : yX & yS : "$spaln_params
fi
while read sp
do
x=`echo $sp | awk '{print$1}'`
y=`echo $sp | awk '{print$2}'`
#Final spaln output in multiple formats : each with different info
spaln -M1 -l280 -O4 -ya0 -yX$x -yS$y sblast_hit.fa spaln_query.aa > spaln_out 2>/dev/null
spaln -M1 -l280 -O1 -ya0 -yX$x -yS$y sblast_hit.fa spaln_query.aa > spaln_aln_out 2>/dev/null
spaln -M1 -l280 -O2 -ya0 -yX$x -yS$y sblast_hit.fa spaln_query.aa > spaln_gff3_out 2>/dev/null
bedtools getfasta -fi sblast_hit.fa -bed <(awk '!/^#|^@/{print"sblast_hit",$9-1,$10,"exon_"++i}' OFS="\t" spaln_out) -name > spaln_exon_wise.fa 2>/dev/null
awk '!/^#|^@/ {print$NF}' spaln_out | tr "." " " | awk NF | tr '[:lower:]' '[:upper:]' > spaln_splice_sites
done < <(echo $spaln_params) 

#Generate ORF
transeq <(grep -v ">" spaln_exon_wise.fa | paste -s -d "" | sed "1i \>$species") spaln_exon_wise_orf.fa 2>/dev/null
#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" spaln_exon_wise.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" spaln_exon_wise_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > spaln_translated
#For better visualization of ORF
gnaa=`head -1 spaln_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk '{print"<(cut -c "$1"-"$2" spaln_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 spaln_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > spaln_nucleotide_amino_acid_view
#Pairwise exon-wise algnment between query & final consensus 
exalign2 $qp spaln_exon_wise.fa > spaln_query_hit_pairwise_exon.aln 2>/dev/null
#QC spaln output : Exon boundar indels
#Exon boundary insertions
while mapfile -t -n 4 ary && ((${#ary[@]}))
do 
ex=`printf '%s\n' "${ary[@]}" | head -1 | tr -d ">"`
printf '%s\n' "${ary[@]}" > tmp
cat tmp | grep "spaln_exon" | grep [0-9] -aob | head -1 | cut -f1 -d ":" | awk '{print$1+3}' | xargs -n1 sh -c 'cat tmp | cut -c $0-' 2>/dev/null | awk NF | sed 's/ \+[0-9]\+//g' > tmp1
#insertions at exon start : counted if a dash is observed within the query sequence anywhere within the first 5% of pairwise alignment between query & subject
is=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"query",$0}' | awk '/query/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $is == "" ]]; then is="-"; else :; fi
#insertions at exon end : counted if a dash is observed within the query sequence anywhere after 95% of pairwise alignment between query & subject
ie=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"query",$0}' | awk '/query/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ie == "" ]]; then ie="-"; else :; fi
#deletions at exon start : counted if a dash is observed within the subject sequence anywhere within the first 10% (since subject sequence is longer due to flanking regions) of pairwise alignment between query & subject
ds=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"query",$0}' | awk '/spaln_exon/ && $NF<6' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
if [[ $ds == "" ]]; then ds="-"; else :; fi
#deletions at exon end : counted if a dash is observed within the subject sequence anywhere after 90% of pairwise alignment between query & subject
de=`cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"query",$0}' | awk '/spaln_exon/ && $NF>94' | awk '{print$2 | "sort | uniq -c "}' | awk '{print$1}'`
#Longest continous insertions : Position, Length & % of gene in which it appears
if [[ $(head -1 tmp1 | grep "\-") == "" ]]; then lci="- - -"; else
lci=`bedtools merge -i <(cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"query",$0}' | awk '/query/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
fi
#Longest continous deletions : Position, Length & % of gene in which it appears
if [[ $(tail -1 tmp1 | grep "\-") == "" ]]; then lcd="- - -"; else
lcd=`bedtools merge -i <(cat tmp1 | tr -d "|." | awk NF | sed 's/ //g' | sed 's/\t//g' | awk -v FS="" '{$1=$1};1' | awk '{for(i=1;i<=NF;i++) if($i~"-") print NR,NF,i,i/NF*100}' | awk '{if($1%2==0) print"spaln_exon",$0; else print"query",$0}' | awk '/spaln_exon/{print$4,$5}' | awk '{$2+=0}1' CONVFMT="%.1f" | awk '{print"Same",$1,$1+1,$2}' OFS="\t") -c 4 -o count,min,max | awk '{if($4>1) print$2"-"$3,$4,$5"-"$6; else print$2,$4,$5}' | sort -k2,2 -nr | awk 'NR==1{print$0}'`
#lci=`cat tmp2 | grep "\-\+" -o | tr -d " " | awk '{print length}' | sort -nr | head -1`
fi
if [[ $de == "" ]]; then de="-"; else :; fi
echo $ex $is $ie $ds $de $lci $lcd
unset ex is ie ds de lci lcd
rm tmp tmp1
done < <(cat spaln_query_hit_pairwise_exon.aln | sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" | awk NF) | sed '1i Exons Ins_in_1st_5%_exon Ins_in_last_95%_exon Del_in_1st_5%_exon Del_in_last_95%_exon Longest_continous_ins_position(lci) lci_length lci_%_range Longest_continous_del_position(lcd) lcd_length lcd_%_range' | column -t > spaln_boundary_indels_wrt_query
#Create consensus from spaln
cp spaln_exon_wise.fa final_consensus.fa
unset species qp o splice_sites a1t a1 a2t a2 enh pqe1 pqe2 spn sp x y spaln_params gnaa
else
echo "     ↳ Spaln didn't run"
fi

unset y
cd /home/neo/bird_db1/aswin/APOBEC1
done < <(grep -if sblast_exonerate_spaln_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#12.3. Summary of sblast-exonerate-spaln

while read z
do
y=`echo $z | awk '{print$1}' | cut -f1,2 -d "_"`
cd $(echo $z | awk '{print$2}')/aswin/APOBEC1/sblast_exonerate_spaln
#Taxonomy order of the species
gr=`echo $y | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} ~/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
#SRA size
a=`grep $y ~/bird_db1/aswin/database_details/all_sra_size | awk '{print$NF}'`
#Query length
sql=`grep -v ">" ../sblast_query.fa | wc | awk '{print$3-$1}'`
#Exonerate consensus length
exl=`grep -v ">" sblast_exonerate.out.fa | wc | awk '{print$3-$1}'`
#Spaln consensus length
spl=`grep -v ">" spaln_exon_wise.fa | wc | awk '{print$3-$1}'`
#START codon
stc=`grep -v ">" final_consensus.fa | paste -s -d "" | cut -c -3 | tr [:lower:] [:upper:]`
#splice sites
splice_sites=`grep ">" ../sblast_query.fa -c | awk '{print$0-1}'`
a1t=`cat sblast_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$1=="GT" || $1=="GC"' | wc -l`
if [[ $a1t == $splice_sites ]]; then a1="-"; else a1="x"; fi
a2t=`cat sblast_exonerate.out | awk '/splice_site/ {print$NF}' | tr -d '"' | paste -d " " - - | tr [:lower:] [:upper:] | awk '$2=="AG"' | wc -l`
if [[ $a2t == $splice_sites ]]; then a2="-"; else a2="x"; fi
b1t=`cat spaln_splice_sites | awk '$1=="GT" || $1=="GC"' | wc -l`
if [[ $b1t == $splice_sites ]]; then b1="-"; else b1="x"; fi
b2t=`cat spaln_splice_sites | awk '$2=="AG"' | wc -l`
if [[ $b2t == $splice_sites ]]; then b2="-"; else b2="x"; fi
#insertions
if [[ -f spaln_exon_wise.fa ]]
then
m=`needle <(grep -v ">" spaln_exon_wise.fa | sed '1i >spaln') <(grep -v ">" ../sblast_query.fa | sed '1i >query') -auto -stdout -awidth 280 | grep "^query" | awk '{print$3}' | tr -cd "-" | wc | awk '{print$3-$1}'`
if [[ $m == "" ]]; then m="-"; else :; fi
elif [[ -f sblast_exonerate.out ]]
then
m=`cat sblast_exonerate.out | awk '{for(i=1;i<=NF;i++) if($i~/insertions/) print$(i+1)}' | awk '{sum+=$1;} END{print sum}' | awk '{print$1}'`
if [[ $m == "" ]]; then m="-"; else :; fi
else
m="-"
fi
#deletions
if [[ -f spaln_exon_wise.fa ]]
then
n=`needle <(grep -v ">" spaln_exon_wise.fa | sed '1i >spaln') <(grep -v ">" ../sblast_query.fa | sed '1i >query') -auto -stdout -awidth 280 | grep "^spaln" | awk '{print$3}' | tr -cd "-" | wc | awk '{print$3-$1}'`
if [[ $n == "" ]]; then n="-"; else :; fi
elif [[ -f sblast_exonerate.out ]]
then
n=`cat sblast_exonerate.out | awk '{for(i=1;i<=NF;i++) if($i~/deletions/) print$(i+1)}' | awk '{sum+=$1;} END{print sum}' | awk '{print$1*3}'`
if [[ $n == "" ]]; then n="-"; else :; fi
else
n=`cat test.out | sed 1,2d | awk '!seen[$1]++' |  awk '{print$NF}' | tr -cd '-' | wc | awk '{print$3-$1}'`
fi
#Gaps
gp=`sed '/^>/! s/[a-z]/\U&/g' final_consensus.fa | grep N -B1 | tr -d [ATGC] | sed 's/^--//g' | awk NF | paste -d " " - - | sed 's/>exon_//g' | awk '{print length($2)"("$1")"}' | paste -s -d ","`
if [[ $gp == "" ]]; then gp="-"; else :; fi
#Query covered
if [[ -s spaln_exon_wise.fa ]]
then
ql=`grep -v ">" spaln_query.aa | wc | awk '{print$3-$1}'`
qr=`cat spaln_aln_out | grep "^>sblast_hit" | grep -oP '\(\K[^\)]+' | tr -d " " | awk 'NR==2'`
qc=`echo $qr | awk -v l="$ql" -F "-" '{print($2-$1+1)/l*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
elif (( $(grep -v ">" sblast_exonerate.out.fa 2>/dev/null | wc | awk '{print$3-$1}') > 0 ))
then
ins=`awk '{for(i=1;i<=NF;i++) if($i~/insertions/) print $(i+1)}' sblast_exonerate.out | awk '{s+=$1} END{print s}' | awk '{print$1}'`
del=`awk '{for(i=1;i<=NF;i++) if($i~/deletions/) print $(i+1)}' sblast_exonerate.out | awk '{s+=$1} END{print s}' | awk '{print$1*3}'`
fs=`awk '{for(i=1;i<=NF;i++) if($i~/frameshifts/) print $(i+1)}' sblast_exonerate.out | awk '{s+=$1} END{print s}'`
if [[ $fs == "" ]]; then fs="-"; else :; fi
tl=`grep -v ">" sblast_exonerate.out.fa | wc | awk '{print$3-$1}'`
ql=`grep -v ">" ../sblast_query.fa | wc | awk '{print$3-$1}'`
if [[ $fs == "-" ]]; then
qc=`calc $(calc $tl - $ins + $del) / $ql \* 100 | tr -d "~" | awk '{$1+=0}1' CONVFMT="%.2f"`
else
qc=`calc $(calc $tl - $ins + $del - $fs) / $ql \* 100 | tr -d "~" | awk '{$1+=0}1' CONVFMT="%.2f"`
fi
else
qc="-"
fi
unset tl
# % of gene in which 1st STOP codon appears
if [[ -s spaln_exon_wise.fa ]]
then
tl=`grep -v ">" spaln_query.aa | wc | awk '{print$3-$1}'`
sp=`tail -1 spaln_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$tl" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.1f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
elif (( $(grep "[A-Z]" sblast_exonerate_translated 2>/dev/null | wc -l) > 1 ))
then
tl=`grep -v ">" sblast_exonerate.out.fa_orf.fa | wc | awk '{print$3-$1}'`
sp=`tail -1 sblast_exonerate_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$tl" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.1f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
else 
sp="-"
fi
#check ORF
if [[ -s spaln_exon_wise.fa ]]
then
co=`grep -v ">" spaln_exon_wise_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $co == "" ]] ; then orf="❌";else orf="✅"; fi
elif [[ -f sblast_exonerate.out.fa_orf.fa ]]
then
co=`grep -v ">" sblast_exonerate.out.fa_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $co == "" ]] ; then orf="❌";else orf="✅"; fi
else
orf="-"
fi
echo $y $gr $a $sql $exl $spl $stc $m $n $gp $a1 $a2 $b1 $b2 $qc $sp $orf
unset y gr a sql exl spl stc splice_sites a1t a1 a2t a2 b1t b1 b2t b2 m n gp ql qr qc ins del fs tl sp co orf
cd /home/neo/bird_db1/aswin/APOBEC1
done < <(grep -if sblast_exonerate_spaln_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths) | sed '1i Species Group SRA_size query_L exonerate_L spaln_L START Ins Del Gaps ss1_ex ss2_ex ss1_sp ss2_sp Query_cov 1st_STOP_in_gene ORF' | column -t > sblast_exonerate_spaln_summary

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#12.4. Add completely obtained ORFs to the total validated queries list, gene tree & add their sequences to the exon_wise folder

#List of queries that gave complete ORF with 100% query coverage from sblast exonerate_spaln
awk '$7=="ATG" && $10=="-" && (($11=="-" && $12=="-") || ($13=="-" && $14=="-")) && $15~"100" && $16~"100" && $NF=="✅"' sblast_exonerate_spaln_summary > sblast_exonerate_spaln_validated_queries
cat sblast_exonerate_spaln_validated_queries | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' >> total_validated_queries
#Add to gene tree
for i in `cat sblast_exonerate_spaln_validated_queries`; do sed "/$i/ s/$/"sblast_exonerate_spaln_"$i/g" -i gene_tree; unset j;done
#add sequence to folder : sblast consensus sequences are not exon-wise
while read i
do
j=`echo $i | awk '{print$1}'`
cd $(echo $i | awk '{print$2}')/aswin/APOBEC1/sblast_exonerate_spaln/
cat final_consensus.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/"APOBEC1_"$j.fa
cat ../optimal_query_sblast.outfmt3.fa > /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/"APOBEC1_"$j"_flanking".fa
unset j
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if sblast_exonerate_spaln_validated_queries /home/neo/bird_db1/aswin/database_details/all_genome_paths)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#12.5. QC validated sequences

#Species with focal gene fragmented into different scaffolds or chromosomes
cut -f2- -d "_" total_validated_queries | sed 's/.fa//g' | xargs -n1 sh -c 'grep "^$0" gblast_principal_components' | column -t | awk '/ diff / {print$1}' > gene_breaking_assemblies

./qc_validated.sh validated_sequences validated_sequences_with_flanking_regions *_coding_exons_set APOBEC1 /home/neo/bird_db1/aswin/APOBEC1 sblast_exonerate_spaln_qc

#Find potential false positives in validated sequences by detecting outliers in nucleotide MSA
sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" sblast_exonerate_spaln_qc | awk -F ":" '/Outlier/ {print$2}' | sed 's/Nothing//g' | tr "," "\n" | awk NF | sed 's/^[ ]\+//g' > outliers

#Remove outliers from the validated sequence records
cat outliers | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' > outliers2
if [[ -s outliers2 ]]; then 
awk -iinplace 'NR==FNR{a[$0];next} !($0 in a)' outliers2 total_validated_queries
else :; fi
cat outliers | xargs -n1 sh -c 'sed "/$0/ s/[ ]\+[^ ]\+$//g" gene_tree -i'
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa) | xargs rm &> /dev/null
grep -if outliers <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/*.fa) | xargs rm &> /dev/null
rm outliers outliers2

##################################################################################################################################################################################################################################################################################################################
#13. TOGA 

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#13.1. Collect all species to run TOGA

grep -if sblast_exonerate_spaln_validated_queries <(grep -if sblast_validated_queries sblast_subjects -v) -v > toga_subjects

#Choose the best reference species : 2 criterias : 
    #1. Intact focal & flanking genes
    #2. Genomic annotation available
    
#In this case use duck as reference for all birds; since duck has good annotation
#If some birds doesn't give a good alignment then used a more closer species

#Before running TOGA : update genomes, run gblast_short & gedit on these

#Run TOGA
#Run these scripts : 
    #toga_apobec1_chicken.sh        #for testing TOGA run for a single species
    #apobec1_galliformes_toga.sh    #Run for the whole group
    #apobec1_toga_subjects_toga.sh  #Run for all species lacking ORF from cutom pipeline

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#13.2. Transfer output to neo

#In neo
mkdir ~/bird_db1/aswin/APOBEC1/TOGA
cd ~/bird_db1/aswin/APOBEC1/TOGA

#In morpheus
cd /home/morpheus/aswin/APOBEC1/TOGA
cat supplementary/galliformes_to_run_toga supplementary/toga_subjects_to_run_filtered | sed '$i Gallus_gallus' | paste -s -d " " | xargs -I {} sh -c 'scp -r {} neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/TOGA/'
scp /home/morpheus/aswin/APOBEC1/TOGA/birds_APOBEC1_XM_027449866.2_inactivating_mutations.png neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/TOGA/

#In neo
cd ~/bird_db1/aswin/APOBEC1/TOGA

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#13.3. TOGA summary (TOGA create output for each isoforms of the focal gene; In case of APOBEC1 there is a single isoform that contains all possible exons, hence create a summary for just that isoform)

for x in `grep APOBEC1 Alectura_lathami/chr1_duck_isoforms.txt | awk '{print$2}' | sed -n 2p`
do
for i in `ls -d */ | tr -d "/"`
do
gr=`echo $i | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} ~/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
#Since TOGA for Chlamydotis_macqueenii was rerun with a closer reference species due to no results when using duck as reference; hence the annotation is different henceforth the transcript id ($x) as well 
cd $i
if [[ $i == "Chlamydotis_macqueenii" ]]; then prevx=`echo $x`; unset x; x="APOBEC1_mrna"; scaffold=`ls scaffold_*.fa | grep "Chlamydotis_macqueenii"`;else scaffold=`ls scaffold_*.fa` ; fi
#Transcript quality
j=`grep "$x" toga_output/temp/transcript_quality.tsv`
if [[ $j == "" ]]; then j="-"; else :; fi
#Degree of loss / intact
k=`grep "$x" toga_output/loss_summ_data.tsv | awk '/TRANSCRIPT/{print$NF}'`
if [[ $k == "" ]]; then k="-"; else :; fi
#Detected in expected region or not
l=`grep "$x" toga_output/temp/exons_meta_data.tsv | awk '{print$6}'`
if [[ $l == "" ]]; then l="-"; else :; fi
#Exons with Gaps
exg=`grep "$x" toga_output/temp/exons_meta_data.tsv | awk '{print ++i,$9}' | grep "GAP" | awk '{print$1}' | paste -s -d ","`
if [[ $exg == "" ]]; then exg="-"; else :; fi
#Exon class
m=`grep "$x" toga_output/temp/exons_meta_data.tsv | awk '{print$10}'`
if [[ $m == "" ]]; then m="-"; else :; fi
#Orthology classification
n=`grep "$x" toga_output/orthology_classification.tsv | awk '{print$NF}'`
if [[ $n == "" ]]; then n="-"; else :; fi
#Extract consensus sequence
awk -v s="$x" '{if($0~s && $0!~"ref_") print">"$0}' RS=">" toga_output/nucleotide.fasta | tr -d "-" | awk NF > toga_"$x"_consensus.fa

#extract exon-wise consensus sequence
if [[ $(grep "$x" toga_output/query_annotation.bed | awk '{print$6}') == "+" ]]
then
bedtools getfasta -fi $scaffold -bed <(bed12ToBed6 -i <(grep "$x" toga_output/query_annotation.bed)) -s > toga_"$x"_exon_wise_consensus.fa
else
bedtools getfasta -fi $scaffold -bed <(bed12ToBed6 -i <(grep "$x" toga_output/query_annotation.bed)) -s | sed 's/^>/\x00&/' | sort -z -Vr | tr -d '\0' > toga_"$x"_exon_wise_consensus.fa
fi
exon_no=`grep "$x" toga_output/temp/exons_meta_data.tsv | awk '{if($6=="INC" && $9!="GAP" && ($10=="A+" || $10=="A" || $10=="B")) print"exon_"($2+1)}'`
while read z
do 
y1=`echo $z | awk '{print$1}'`
y2=`echo $z | awk '{print$2}'`
sed -i "s/$y2/$y1/g" toga_"$x"_exon_wise_consensus.fa
unset y1 y2
done < <(paste <(echo $exon_no | tr " " "\n") <(grep ">" toga_"$x"_exon_wise_consensus.fa | tr -d ">"))
#Generate ORF
transeq <(grep -v ">" toga_"$x"_consensus.fa | paste -s -d "" | sed "1i \>$j") toga_"$x"_consensus_orf.fa 2>/dev/null

#Align complete nucleotides & corresponding translated amino acid, highlighting the START, STOP & anonymous sites
cat <(grep -v ">" toga_"$x"_consensus.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" toga_"$x"_consensus_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > toga_"$x"_consensus_translated
#For better visualization of ORF
gnaa=`head -1 toga_"$x"_consensus_translated | wc | awk '{print$3-$1}'`
eval "$(paste <(seq 253 252 $gnaa) <(seq 504 252 $gnaa) | awk -v u="$x" '{print"<(cut -c "$1"-"$2" toga_"u"_consensus_translated)"}' | paste -s -d " " | sed 's/^/cat <(cut -c -251 toga_"$x"_consensus_translated) /g' | sed 's/$/ | GREP_COLORS="mt=41" egrep "^ M |$| \\* |$| X |$" --color=always/g')" > toga_"$x"_consensus_nucleotide_amino_acid_view
#Consensus length
cl=`grep -v ">" toga_"$x"_consensus.fa | wc | awk '{print$3-$1}'`
if [[ $cl == "" ]]; then cl="-"; else :; fi
#% Query covered
rl=`grep -v ">" ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Anas_platyrhynchos.fa | wc | awk '{print$3-$1}'`
#Start codon of consensus
start=`grep -v ">" toga_"$x"_consensus.fa | cut -c -3`
if [[ $start == "" ]]; then start="-"; else :; fi
# % of the gene in which the first Stop codons appears
ol=`grep -v ">" toga_"$x"_consensus_orf.fa | paste -s -d "" | wc | awk '{print$3-$1}'`
sp=`tail -1 toga_"$x"_consensus_translated | awk '{for(z=1;z<=NF;z++) if($z~/*/) print z}' | head -1 | awk -v y="$ol" '{print$1/y*100}' | awk '{$1+=0}1' CONVFMT="%.2f"`
if [[ $sp == "" ]]; then sp="-"; else :; fi
#Check ORF presence absence
op=`grep -v ">" toga_"$x"_consensus_orf.fa | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $op == "" ]] ; then o="❌";else o="✅"; fi

#Create an exon-wise MSA of all consensus sequences made by each methods along with the flanking genome & SRA region
cd $(grep "$i" /home/neo/bird_db1/aswin/database_details/all_genome_paths | awk '{print$2}')/aswin/APOBEC1
if [[ -d 2nd_gblast ]]
then
cd 2nd_gblast
exalign2 -msa8 ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Anas_platyrhynchos.fa gblast_edited_consensus.fa Query_exon_combined_exonerate.out_exon_wise.fa spaln_exon_wise.fa ~/bird_db1/aswin/APOBEC1/TOGA/$i/toga_*exon_wise_consensus.fa extracted_flanking_region.fa ../sblast_edited_consensus.fa ../optimal_query_sblast.outfmt3.fa > ~/bird_db1/aswin/APOBEC1/TOGA/$i/All_consensus_msa.aln 2> /dev/null
else 
exalign2 -msa8 ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Anas_platyrhynchos.fa gblast_edited_consensus.fa Query_exon_combined_exonerate.out_exon_wise.fa spaln_exon_wise.fa ~/bird_db1/aswin/APOBEC1/TOGA/$i/toga_*exon_wise_consensus.fa extracted_flanking_region.fa sblast_edited_consensus.fa optimal_query_sblast.outfmt3.fa > ~/bird_db1/aswin/APOBEC1/TOGA/$i/All_consensus_msa.aln 2> /dev/null
fi

echo $i $gr $j $l $m $exg $k $n $cl $start $sp $o
if [[ $i == "Chlamydotis_macqueenii" ]]; then unset x; x=`echo $prevx`; unset prevx; else :; fi
unset j gr scaffold k l m exg n exon_no z gnaa cl start ol sp op o
cd /home/neo/bird_db1/aswin/APOBEC1/TOGA
done | sed '1i Species Group Isoform Transcript_quality Ex1_loc Ex2_loc Ex3_loc Ex4_loc Ex5_loc Ex1_class Ex2_class Ex3_class Ex4_class Ex5_class Exons_with_gaps Presence Orthology Consensus_length START 1st_STOP_in_gene ORF' | column -t > $x"_toga_summary"
unset i
cd /home/neo/bird_db1/aswin/APOBEC1/TOGA
done
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#13.4. Add completely obtained ORFs to the total validated queries list, gene tree & add their sequences to the exon_wise folder

#List of queries that gave complete ORF with 100% query coverage from TOGA
awk '$5=="INC" && $6=="INC" && $7=="INC" && $8=="INC" && $9=="INC" && $16=="I" && $NF=="✅" {print$1}' TOGA/rna-XM_027449866.2_toga_summary > toga_validated_queries

cat toga_validated_queries | sed 's/^/APOBEC1_/g' | sed 's/$/.fa/g' >> total_validated_queries
#Add to gene tree
for i in `cat toga_validated_queries`; do sed "/$i/ s/$/"TOGA_"$i/g" -i gene_tree; unset j;done
#add sequence to folder : sblast consensus sequences are not exon-wise
cat toga_validated_queries | xargs -n1 bash -c 'paste <(echo $0) <(ls ~/bird_db1/aswin/APOBEC1/TOGA/$0/toga*exon_wise_consensus.fa)' | xargs -n2 sh -c 'cp $1 ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_"$0".fa'

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#13.5. QC validated sequences

./qc_toga_validated.sh validated_sequences *_coding_exons_set toga_validated_qc

#Find potential false positives in validated sequences by detecting outliers in nucleotide MSA
#sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" toga_validated_qc | awk -F ":" '/Outlier/ {print$2}' | sed 's/Nothing//g' | tr "," "\n" | awk NF | sed 's/^[ ]\+//g' > outliers

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Manual inspection of DNA & RNA in Uncertainly lost/intact gene & gaps containing species

awk '$16=="UL" {print$1}' TOGA/rna-XM_027449866.2_toga_summary > toga_uncertain_loss
#Species that are marked as intact and have all intact exons but lacking complete ORF
grep -if toga_validated_queries <(awk '$16=="I"' TOGA/rna-XM_027449866.2_toga_summary) -v | awk '{print$1}' > toga_uncertain_intact
#The species with gaps given by TOGA doesn't necessarily mean that the species have N's in their genome, they might or might not have N's. Hence we need to manuall inspect them
awk '$15!="-" {print$1}' TOGA/rna-XM_027449866.2_toga_summary | grep -v "Species" > toga_species_with_gaps
awk '$10=="N/A" || $11=="N/A" || $12=="N/A" || $13=="N/A" || $14=="N/A"' TOGA/rna-XM_027449866.2_toga_summary | awk '{print$1}' | sort -u > toga_species_with_missing_exons

#Inspect DNA
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Buceros_rhinoceros/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" > toga_exon_metadata.bed
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Buceros_rhinoceros/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t"  > toga_exon_metadata_exp_loc.bed 

#Inspect RNA
cd ~/bird_db1/aswin/APOBEC1	#In neo
scp toga_uncertain_loss toga_uncertain_intact sagar@172.30.1.172:/media/sagar/disk4/RNAseq_store/Birds/aswin/APOBEC1/TOGA/
cd /media/sagar/disk4/RNAseq_store/Birds/
grep -if aswin/APOBEC1/TOGA/toga_uncertain_loss <(ls -d */) | xargs -n1 sh -c 'echo ">"$0; ls -lhS $0 | grep "bam$" | grep -v "subset" | awk "{print \"- \"\$NF,\$5}"' | column -t > aswin/APOBEC1/TOGA/rna_database_toga_uncertain_loss_species

#Commands used is in the script "apobec1_manual_narrow_search.sh"

#TOGA uncertain loss :

 #Buceros_rhinoceros      - intact - exon1 found; but starting 2 bases of CDS start codon missing; but an in-frame start codon is found 123 bases upstream  
 #Cathartes_aura          - loss   - exon2(deletion) & exon4 bodaries are uncertain; Even all possible variants couldn't form complete ORF; most likely inactivation events are deletion in exon2 with stop at exon4/ stop cdon in exon4(when splice site is included) 
 #Chlamydotis_macqueenii  - loss   - exon_2 found; but exon3 has indels 
 #Gymnogyps_californianus - loss   - loss due to either deletion in exon_1 with stop codon at exon_4 or stop codon in exon_3 (if splice site is considered)
 #Numida_meleagris        - loss   - whole exons are missing; don't 
 #Phoenicopterus_ruber    - intact - assembly had error which was corrected by SRA
 #Tyto_alba               - missing data - exclude this species as it is almost(85%) intact & if it is completely inact it would be very time consuming to prove as it requires more data & even if it has a stop codon it would be in the >85% of the gene length, which is not a reliable way to establish gene loss

#TOGA uncertain intact :

 #Camarhynchus_parvulus - loss - start was missing when upstream start was added it includes a stop codon
 #Geospiza_fortis       - loss - start was missing when upstream start was added it includes a stop codon
 #Limosa_lapponica      - Intact - stop was missing; found a stop codon 66bp downstream
 #Picoides_pubescens    - missing data - The genes is at the end of chromosome' the assembly very poor to recover

#TOGA species with missing exons
 #Meleagris_gallopavo - loss - exon_5 is completely missing in genome
 #Penelope_pileata    - loss - exon_5 is completely missing in genome

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Add completely obtained ORFs to the total validated queries list, gene tree & add their sequences to the exon_wise folder
#Manually validated losses "/home/neo/bird_db1/aswin/APOBEC1/manually_validated_losses" are also saved in a folder "~/bird_db1/aswin/APOBEC1/validated_lost_genes/"

#List of queries that gave complete ORF with 100% query coverage from manual inspection

cat manually_validated_queries >> total_validated_queries
#Add to gene tree
for i in `cat manually_validated_queries`; do sed "/$i/ s/$/"manually_validated_"$i/g" -i gene_tree; unset j;done
for i in `cat manually_validated_losses`; do sed "/$i/ s/$/"manually_validated_loss_"$i/g" -i gene_tree; unset j;done
for i in `cat species_to_exclude`; do sed "/$i/ s/$/"lack_of_data_"$i/g" -i gene_tree; unset j;done

##################################################################################################################################################################################################################################################################################################################
#13.6. Synteny from TOGA output

for x in `ls -d */ | tr -d "/"`
do
gr=`echo $x | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} ~/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $x
#There are some exceptional replacements for Chlamydotis_macqueenii; since the query for it is not duck
for i in `awk '{print$4}' toga_output/query_annotation.bed | cut -f1,2 -d "." | sed -e 's/AICDA_mrna.1/AICDA_mrna/g' -e 's/APOBEC1_mrna.1/APOBEC1_mrna/g' -e 's/MFAP5_mrna.1/MFAP5_mrna/g'`
do 
j=`grep "$i" *_isoforms.txt | awk '{print$1}'`
sed -n "s/$i/$j/p" toga_output/query_annotation.bed
unset j
done | sort -k2,2 -k3,3 -n | awk '{print$4}' | awk '!a[$1]++' > tmp_synteny
y=`cat tmp_synteny | grep -i "APOBEC1" -C3 --color=always | paste -s -d " " |  awk '{for(i=1;i<=NF;i++) if($i~/\[01;31m/ && i>3) print$(i-3),$(i-2),$(i-1),$i,$(i+1),$(i+2),$(i+3); else if($i~/\[01;31m/ && i==3) print"-",$(i-2),$(i-1),$i,$(i+1),$(i+2),$(i+3); else if($i~/\[01;31m/ && i==2) print"-","-",$(i-1),$i,$(i+1),$(i+2),$(i+3); else if($i~/\[01;31m/ && i==1) print"-","-","-",$i,$(i+1),$(i+2),$(i+3)}'`
if [[ $y == "" ]]
then
#Assuming the upstrem gene is AICDA and the annotation is in 
y=`cat tmp_synteny | grep -i "AICDA" -C3 | paste -s -d " " |  awk '{for(i=1;i<=NF;i++) if($i~/AICDA/ && i>3) print$(i-2),$(i-1),$i,"<No_hits>",$(i+1),$(i+2),$(i+3); else if($i~/AICDA/ && i==3) print$(i-2),$(i-1),$i,"<No_hits>",$(i+1),$(i+2),$(i+3); else if($i~/AICDA/ && i==2) print"-",$(i-1),$i,"<No_hits>",$(i+1),$(i+2),$(i+3); else if($i~/AICDA/ && i==1) print"-","-",$i,"<No_hits>",$(i+1),$(i+2),$(i+3) }'`
else :; fi
echo $x $gr $y
unset y i
cd ~/bird_db1/aswin/APOBEC1/TOGA
done | sed 's/\.[0-9]\+//g' | sed 's/LOC101795050/NANOG_0/g' | sed 's/GeneID_101791705/NANOG_1/g' | column -t > temp_synteny

cat temp_synteny | awk '{if(NF==6) print$0, "-","-","-"; else if(NF==7) print$0,"-","-"; else if(NF==8) print$0,"-"; else if(NF==2) print$0,"-","-","-","-","-","-","-"; else print$0}' | column -t > temp_synteny1

lg="AICDA"
rg="NANOG_0"
cat temp_synteny1 | awk '{if($4=="")print$0,"-","-","-","-","-","-"; else print$0}' | awk '{if($5=="")print$0,"-","-","-","-","-"; else print$0}' | awk '{if($6=="")print$0,"-","-","-","-"; else print$0}' | awk '{if($7=="")print$0,"-","-","-"; else print$0}' | awk '{if($8=="")print$0,"-","-"; else print$0}'  | awk '{if($9=="")print$0,"-"; else print$0}' | awk -v l="$lg" -v r="$rg" -v IGNORECASE=1 '{if($5~l || $7~r) print$0; else if($5~r || $7~l) print$1,$2,$9,$8,$7,$6,$5,$4,$3; else print$0}' | column -t > TOGA_Synteny
rm temp_synteny*

##################################################################################################################################################################################################################################################################################################################
#14. Pipeline results

#14.1.  Print summary of Custom pipeline & TOGA

cd ~/bird_db1/aswin/APOBEC1

while read z
do
y=`echo $z | awk '{print$1}' | cut -f1,2 -d "_"`
#Taxonomy order of the species
gr=`echo $y | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} ~/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $z | awk '{print$2}')/aswin/APOBEC1/
#Best genome search results from custom pipeline
if [[ -d new_genome_2nd_gblast ]]
then
cd new_genome_2nd_gblast
elif [[ -d 2nd_gblast ]]
then
cd 2nd_gblast
else :
fi
cpg=$(for i in `ls APOBEC1_*.fa | grep -v "flanking" | head -1 | xargs -I {} bash -c 'grep ">" {} | tr -d ">"'`
do
excov=`grep "$i" gblast_edited_query_covered | awk '{$NF+=0}1' CONVFMT="%.1f" | awk '{print$NF}'`
if [[ $excov == "" ]]; then excov="-"; else :; fi
echo $excov
unset excov
done)
#Genome search results from TOGA pipeline
tpg=`grep "$y" /home/neo/bird_db1/aswin/APOBEC1/TOGA/rna-XM_027449866.2_toga_summary | awk '{print$10,$5,$6,$7,$8,$9,$11,$12,$13,$14,$16,$15}'`
#Best SRA search results from custom pipeline
cd $(echo $z | awk '{print$2}')/aswin/APOBEC1/
cps=`awk '{$NF+=0}1' CONVFMT="%.1f" optimal_query_sblast_summary | awk '!/Query_Length/ {print$NF}' | sed 's/^0$/-/g' | paste -s -d " "`
echo -e $y $gr $cpg $tpg $cps
unset y gr cpg i excov tpg cps
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if toga_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths | grep -v "Lyrurus_tetrix") | sed '1i Species Group OQ_gEx1 OQ_gEx2 OQ_gEx3 OQ_gEx4 OQ_gEx5 T_gEx1L T_gEx2L T_gEx3L T_gEx4L T_gEx5l T_gEx1C T_gEx2C T_gEx3C T_gEx4C T_gEx5C T_ortho T_loss OQ_sE1 OQ_sE2 OQ_sE3 OQ_sE4 OQ_sE5' | column -t > pipeline_summary_column_wise

while read z
do
y=`echo $z | awk '{print$1}' | cut -f1,2 -d "_"`
#Taxonomy order of the species
gr=`echo $y | cut -f1,2 -d "_" | xargs -I {} sh -c "echo {}; grep {} ~/bird_db1/aswin/taxonomy/orders_all_birds" | paste -s -d " " | awk '{print$3}'`
cd $(echo $z | awk '{print$2}')/aswin/APOBEC1/
#Best genome search results from custom pipeline
if [[ -d new_genome_2nd_gblast ]]
then
cd new_genome_2nd_gblast
elif [[ -d 2nd_gblast ]]
then
cd 2nd_gblast
else :
fi
cpg=$(for i in `ls APOBEC1_*.fa | grep -v "flanking" | head -1 | xargs -I {} bash -c 'grep ">" {} | tr -d ">"'`
do
excov=`grep "$i" gblast_edited_query_covered | awk '{$NF+=0}1' CONVFMT="%.1f" | awk '{print$NF}'`
if [[ $excov == "" ]]; then excov="-"; else :; fi
echo $excov
unset excov
done)
#gaps in genome
cpgg=$(for b in `ls pairwise_exon_*.aln | sort -V`
do
b1=`echo $b | sed 's/pairwise_exon_//g' | sed 's/.aln//g'`
b2=`grep -v "#" $b  | awk NF | awk '{print$3}' | grep -v "|" | tail -1 | grep -i "n" -o | wc -l`
if [[ $b2 == "" ]]; then b2="0"; else :; fi
echo $b1 $b2
done | awk '$2>0 {print$1}' | paste -s -d ",")
if [[ $cpgg == "" ]]; then cpgg="-"; else :; fi
#Genome search results from TOGA pipeline
tpg1=`grep "$y" /home/neo/bird_db1/aswin/APOBEC1/TOGA/rna-XM_027449866.2_toga_summary | awk '{print$5,$6,$7,$8,$9}'`
tpg2=`grep "$y" /home/neo/bird_db1/aswin/APOBEC1/TOGA/rna-XM_027449866.2_toga_summary | awk '{print$10,$11,$12,$13,$14,$15,$17,$16}'`
#Best SRA search results from custom pipeline
cd $(echo $z | awk '{print$2}')/aswin/APOBEC1/
cps=`awk '{$NF+=0}1' CONVFMT="%.1f" optimal_query_sblast_summary | awk '!/Query_Length/ {print$NF}' | sed 's/^0$/-/g' | paste -s -d " "`
#gaps in SRA
cpsg=`awk -v RS=">" '{$1=$1}1' sblast_edited_consensus.fa | awk NF | awk '$2~"[Nn]" {print$1}' | sed 's/exon_//g' | paste -s -d ","`
if [[ $cpsg == "" ]]; then cpsg="-"; else :; fi
echo -e ">"$y $gr " genome " $cpg $cpgg "\n . . TOGA_Loc" $tpg1 "\n . . TOGA_class" $tpg2 "\n . . SRA " $cps $cpsg"\n"
unset y gr cpg cpgg i excov tpg1 tpg2 cps cpsg
cd /home/neo/bird_db1/aswin/APOBEC1/
done < <(grep -if toga_subjects /home/neo/bird_db1/aswin/database_details/all_genome_paths | grep -v "Lyrurus_tetrix") | sed '1i Species Group Method Exon_1 Exon_2 Exon_3 Exon_4 Exon_5 Exon_with_gaps Orthology TOGA_conclusion' | column -t > pipeline_summary_row_wise

cat pipeline_summary_row_wise | sed '/^[A-Z]/ s/^/\n/g' | GREP_COLORS="mt=01;33" grep "genome.*\|$" --color=always | GREP_COLORS="mt=01;36" grep "SRA.*\|$" --color=always | GREP_COLORS="mt=01;32" grep "TOGA.*\|$" --color=always | egrep -a "\bC\b|\-" -z

#Visualize whole summary in terms of presence absence
cat pipeline_summary_row_wise | awk '{if($4~"[0-9]" || $4=="A" || $4=="B" || $4=="A+") $4="✅"; else $4="❌"}1' \
| awk '{if($5~"[0-9]" || $5=="A" || $5=="B" || $5=="A+") $5="✅"; else $5="❌"}1' \
| awk '{if($6~"[0-9]" || $6=="A" || $6=="B" || $6=="A+") $6="✅"; else $6="❌"}1' \
| awk '{if($7~"[0-9]" || $7=="A" || $7=="B" || $7=="A+") $7="✅"; else $7="❌"}1' \
| awk '{if($8~"[0-9]" || $8=="A" || $8=="B" || $8=="A+") $8="✅"; else $8="❌"}1' | column -t | sed '/>/ s/^/\n/g'


##################################################################################################################################################################################################################################################################################################################
#14.2. Visualise loss in phylogeny

#highlight losses in phylogenetic tree
GREP_COLORS="mt=41" grep -if <(cat <(awk '$16=="L" {print$1}' TOGA/rna-XM_027449866.2_toga_summary) manually_validated_losses) <(GREP_COLORS="mt=42" grep -if <(cat <(cut -f2,3 -d "_" total_validated_queries | sed 's/.fa//g') manually_validated_queries | sed 's!^!\\b!g') gene_tree -z --color=always) -az --color=always 

#Add taxonomic order info to gene tree
cp gene_tree gene_tree_with_order_names
grep '[a-z]' stree | awk '{print$NF}' | xargs -n1 sh -c 'grep $0 ~/bird_db1/aswin/taxonomy/orders_all_birds' | xargs -n2 sh -c 'sed -i "s/\b$0\b/$0 $1/g" gene_tree_with_order_names'
paste <(cat gene_tree_with_order_names | tr -d "[a-z][A-Z][0-9]_" | sed 's/[ ]\+$//g') <(cat gene_tree_with_order_names | tr -d "+|\-\\\/" | column -t | sed -z 's/\n/\n\n/g') > gene_tree_tmp
rm gene_tree_with_order_names
mv gene_tree_tmp gene_tree_with_order_names
GREP_COLORS="mt=41" grep -if <(cat <(awk '$16=="L" {print$1}' TOGA/rna-XM_027449866.2_toga_summary) manually_validated_losses | sed 's/^/\\b/g' | sed 's/$/\\b/g') <(GREP_COLORS="mt=42" grep -if <(cat <(cut -f2,3 -d "_" total_validated_queries | sed 's/.fa//g') manually_validated_queries | sed 's!^!\\b!g') gene_tree_with_order_names -z --color=always) -az --color=always > gene_tree_with_order_names_labelled

cat gene_tree_with_order_names | tr -d "+|-" | tr -d "/\\" | awk NF | awk '{if($3=="") print$0,"lost"; else print$0}' | column -t > gene_tree_table

#Update tree
sed "s/)'[0-9]\+'/)/g" all_birds.nwk -i
sed 's/Chloebia_gouldiae/Erythrura_gouldiae/g' all_birds.nwk -i
sed 's/Cyanoderma_ruficeps/Stachyris_ruficeps/g' all_birds.nwk -i
sed 's/Dryobates_pubescens/Picoides_pubescens/g' all_birds.nwk -i
ptree.sh all_birds.nwk | sed 's/\b \b/_/g' > newtree
cat gene_tree_table | xargs -n3 sh -c 'sed "/$0/ s/$/ $1 $2/g" newtree -i'
GREP_COLORS="mt=41" grep -if <(cat <(awk '$16=="L" {print$1}' TOGA/rna-XM_027449866.2_toga_summary) manually_validated_losses | sed 's/^/\\b/g' | sed 's/$/\\b/g') <(GREP_COLORS="mt=42" grep -if <(cat <(cut -f2,3 -d "_" total_validated_queries | sed 's/.fa//g') manually_validated_queries | sed 's!^!\\b!g') newtree -z --color=always) -az --color=always > new_tree_labelled

#compare 2 trees
R
`
library(phytools)
a <- phytools::read.newick("/home/neo/bird_db1/aswin/APOBEC1/all_birds.nwk")
b <- phytools::read.newick("/home/neo/bird_db1/aswin/tree/species_with_genome_wo_cygnus_olor.nwk")

ab<-cophylo(a,b)
pdf(file = "two_trees.pdf", width = 15, height = 19)
plot(ab,link.type="curved",link.lwd=4,link.lty="solid",link.col=make.transparent("red",0.25))
dev.off()
`
#highlight the APOBEC1 losses in the new and old tree and see if the placement of loss changes significantly
#2 species with loss have a change in placement in updated tree:
#Cuculus canorus
#Chlamydotis macqueeni

##################################################################################################################################################################################################################################################################################################################
#15. Protein Domain search

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#15.1. Domain identification from sequence

mkdir ~/bird_db1/aswin/APOBEC1/Domain_search
#Create multi-fasta fle
cd /home/neo/bird_db1/aswin/APOBEC1/Domain_search

#CDS & protein msa of all validated CDS 
ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | cut -f2- -d "_" | sed "s/.fa//g" | sed "s/^/>/g"; grep -v ">" $0 | paste -s -d ""' > total_validated_queries.fa
transeq <(ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/*.fa | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | cut -f2- -d "_" | sed "s/.fa//g" | sed "s/^/>/g"; grep -v ">" $0 | paste -s -d ""') -auto -stdout \
| sed '/>/! s/[a-z]/\U&/g' | sed '/>/ s/_1$//g' | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > total_validated_queries.aa 

#CDS & protein msa of all toga subjects 
find /home/neo/bird_db1/aswin/APOBEC1/TOGA -maxdepth 2 -name "toga*consensus.fa" | grep -v "exon_wise" | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$(NF-1)}" | sed "s/^/>/g"; grep -v ">" $0' | awk 'BEGIN{RS=">";FS="\n";ORS=""}$2 {print">"$0}' \
| awk '$2!~"N" {print">"$0}' RS=">" | sed '/^>$/g' | awk NF > toga_subjects.fa

transeq <(find /home/neo/bird_db1/aswin/APOBEC1/TOGA -maxdepth 2 -name "toga*consensus.fa" | grep -v "exon_wise" | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$(NF-1)}" | sed "s/^/>/g"; grep -v ">" $0' | awk 'BEGIN{RS=">";FS="\n";ORS=""}$2 {print">"$0}' \
| awk '$2!~"N" {print">"$0}' RS=">" | sed '/^>$/g' | awk NF) -auto -stdout | sed '/>/! s/[a-z]/\U&/g' | sed '/>/ s/_1$//g' | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > toga_subjects.aa
#list of species
cat total_validated_queries/total_validated_queries.fa toga_subjects/toga_subjects.fa | grep ">"| sort -u | tr -d ">" | cut -f1,2 -d "_" > all_final_species
grep '[a-z]' ../stree | awk '{print$NF}' | cut -f1,2 -d "_" | xargs -n1 sh -c 'grep $0 all_final_species' > all_final_species_in_phylo_order

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Using NCBI Batch CD search

	#Download 3 types of data in the following formats :

		  #1. Domain hits   : datat mode   : Consise, Standard, Full
		  #                 : Include domain define
		  #2. Align details : Format     : BLAST Text
		  #                 : datat mode : Full
		  #                 : Include domain define
		  #3. Features
	  
	#Download schematics of domain hits in image format : use chrome browser for better reformatting of web images

		#multi-fasta protein sequence 
			| upload
			v
		#NCBI batch CD search website
			| default parameters
			v
		#Browse results 
			| H-zoom:residue, View:full Results, select all queries
			v
		#Inspect element 
			| browser zoom level : 110%, Delete unwanted elements like header, side panels, footers, Resize display panel & image holder
			v
		#Save in image format using the chrome extension "Coco screenshot" & in pdf format using chrome's default print option in landscape mode

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Using these other tools online:


#Usin protein sequence as query:

#Motif search (https://www.genome.jp/tools/motif/)
#Motif Scan (https://myhits.sib.swiss/cgi-bin/motif_scan): one sequence at a time
#SMART: 

#Works on both DNA & protein
  #ThreaDom online- but one sequence at a time
  #InterPro/InterProScan

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#15.2. Domain alignment

#Using Cobalt 


##################################################################################################################################################################################################################################################################################################################
#16. ACTIVE SITE CHECK

#Check Catalytic sites in species with intact APOBEC1

cd ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries
egrep "H.E|PC..C" -z --color=always total_validated_queries.aa --color=always > active_sites_highlighted
cat ../total_validated_queries/total_validated_queries.aa | egrep ">[^ ]*|HAE.*PC..C" -o > active_site.fa
mafft --auto --reorder --quiet active_site.fa > active_site.aln

#Check the CDS sequence corresponding to whole active site region
GREP_COLORS="mt=07;33" grep -if <(grep "Anas_platyrhynchos" active_site.fa -A1 | tail -1 | sed 's/./&   /g') <(cat ~/soft_links/Anas_platyrhynchos/aswin/APOBEC1/gblast_translated) -z --color=always | less -SR
#CATGCTGAAGTCAACTTCTTGGAAAATTGCTTCAAGCCCTTGTCATCAGCTTCTTGCTCTATCACCTGGGTCTTATCTACTACCCCCTGTGGGAAGTGC

#align all other species
mafft --auto --reorder --quiet total_validated_queries.fa > cds_total_validated_queries.aln

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Check Catalytic sites in species with APOBEC1 loss

#collect different consensus sequences made for species with A1 loss
cd ~/bird_db1/aswin/APOBEC1
for i in `grep -if total_validated_losses_Upupa_epops_removed <(ls -d ~/soft_links/*/ | sed 's!/$!!g')`
do
cd "$i"/aswin/APOBEC1/2nd_gblast/
echo ">"$i
seqtk subseq gblast_edited_consensus.fa <(echo "exon_3") | grep -v ">" | paste -s -d ""
echo $j
unset j
cd ~/bird_db1/aswin/APOBEC1
#remove empty sequences 
done | awk NF | sed 's!/home/neo/soft_links/!!g' | myfasta -de > ~/bird_db1/aswin/APOBEC1/Domain_search/toga_subjects/exon_3_final_consensus.fa

#alignment
#grep -if <(grep ">" ../total_validated_queries/total_validated_queries.fa | tr -d ">") <(grep ">" toga_subjects.fa | tr -d ">") > species_later_validated
#cat ../total_validated_queries/total_validated_queries.fa <(awk 'BEGIN{while((getline<"species_later_validated")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' toga_subjects.fa) > total_validated_queries_toga_subjects.fa
#alignment both species with & without loss to compare catalytic sites
#mafft --auto --reorder --quiet total_validated_queries_toga_subjects.fa > total_validated_queries_toga_subjects.aln

#aligning whole seqeunce of both species with & without loss disrupts the frame of active sites in even species with intact A1
cat ../total_validated_queries/total_validated_queries.fa final_consensus.fa > total_validated_queries_final_consensus.fa
mafft --auto --reorder --quiet total_validated_queries_final_consensus.fa > total_validated_queries_final_consensus.aln

#Out of 33 losses only 22 birds have exon_3, 11 birds lost exon_3

#Align only exon_3
cd ~/bird_db1/aswin/APOBEC1/Domain_search/toga_subjects/
for i in `grep -if <(less ../total_validated_queries/total_validated_queries.fa | tr -d ">") <(ls ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*)`
do
echo $i | awk -F "/" '{print$NF}' | sed 's/APOBEC1_//g' | sed 's/.fa//g' | sed 's/^/>/g'
seqtk subseq $i <(echo "exon_3") | grep -v ">"
done > exon_3_total_validated_queries.fa

#Make the exon_3 in frame
for i in `grep -if <(less ../total_validated_queries/total_validated_queries.fa | tr -d ">") <(ls ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*)`
do
echo $i | awk -F "/" '{print$NF}' | sed 's/APOBEC1_//g' | sed 's/.fa//g' | sed 's/^/>>/g'
myfasta -codex $i -nw -col | grep -a "exon_3" -A1
done | less -RS

sed '/>/! s/^G//g' exon_3_total_validated_queries.fa -i

#Alignment
cat exon_3_total_validated_queries.fa exon_3_final_consensus.fa > exon_3_total_validated_queries_final_consensus.fa
mafft --auto --reorder --quiet exon_3_total_validated_queries_final_consensus.fa > exon_3_total_validated_queries_final_consensus.aln

aliview exon_3_total_validated_queries_final_consensus.aln
#search for this CATGCTGAAGTCAACTTCTTGGAAAATTGCTTCAAGCCCTTGTCATCAGCTTCTTGCTCTATCACCTGGGTCTTATCTACTACCCCCTGTGGGAAGTGC
#In aliview remove some highly unreliable sequences
Chlamydotis_macqueenii
Phasianus_colchicus
Pterocles_gutturalis

#Even after these changes the alignment doesn't make the active sites in frame, then it's difficult to know the protein sequence
#Since the exon_3 is incomplete (especially at exon boundaries) in many species with incomplete A1, it is difficult to make an ORF without finding a frame
#And without translaying into protein, enzyme ative sites, domains are difficult to find (especially as a batch mode)
#Best option is to translate this exon_3 remnant to all 3 forward frames & check for intact active sites in all frames in NCBIU batch CD search

#check ncbi batch CD search
mkdir ~/bird_db1/aswin/APOBEC1/Domain_search/toga_subjects/active_sites
cd ~/bird_db1/aswin/APOBEC1/Domain_search/toga_subjects/active_sites
cp ../exon_3_final_consensus.fa .

transeq exon_3_final_consensus.fa -frame F -auto -stdout > exon_3_final_consensus_3_formard_frames.aa
#Check NCBI batch CD search
scp ceglab8@172.28.65.118:/home/ceglab8/Downloads/ncbi_batch_CD_standard_results_APOBEC1_active_sites_exon_3_cds_3_forward_frames_birds_with_loss.p* .
scp ceglab8@172.28.65.118:/home/ceglab8/Downloads/APOBEC1_exon_3_active_site_birds* .

#Make alignment of exon 3 of 22 birds with A1 loss and 48 birds with intact A1
for i in $(find ~/bird_db1/aswin/APOBEC1/validated_sequences -name "APOBEC1_*.fa")
do
n=$(echo $i | awk -F "/" '{print$NF}' | sed 's/APOBEC1_//g' | sed 's/\.fa//g')
myfasta -codex $i -nw | grep ">exon_3" -A1 | sed "/>/ s/>.*/>$n\_i/g"
done | awk '{if($1~"^>") print$0; else for(i=1;i<=NF;i++) if(length($i)==3) printf $i;print""}' | awk NF > birds_with_inframe_A1_exon_3_final_consensus.fa

cat birds_with_inframe_A1_exon_3_final_consensus.fa exon_3_final_consensus.fa > birds_with_and_wo_A1_exon_3_final_consensus.fa

#Align the A1 with intact active sites
mafft --maxiterate 1000 --globalpair birds_with_inframe_A1_exon_3_final_consensus.fa > birds_with_inframe_A1_exon_3_final_consensus.aln
#View in aliview, location the active site region & manually copy paste that region into a file
#nano birds_with_inframe_A1_exon_3_final_consensus_active_site_region.fa
mafft --maxiterate 1000 --globalpair birds_with_inframe_A1_exon_3_final_consensus_active_site_region.fa > birds_with_inframe_A1_exon_3_final_consensus_active_site_region.aln

#align all A1 this time we know the exact sequence of active site region
mafft --auto birds_with_and_wo_A1_exon_3_final_consensus.fa > birds_with_and_wo_A1_exon_3_final_consensus.aln


for i in $(readlink -f ~/bird_db1/aswin/APOBEC1/validated_sequences/* | head -2)
do
myfasta -codex $i -nw | grep ">exon_3" -A1
done | less -SR

for i in $(readlink -f ~/bird_db1/aswin/APOBEC1/validated_sequences/* | head -2)
do
myfasta -orf $i 
done | less -SR

##################################################################################################################################################################################################################################################################################################################

#After finishing all the results

#Find all possible custom scripts used for APOBEC1 gene loss

#In ceglab8
ls /home/ceglab8/workspace/2.phd/9.research/9.3_gene_loss/cort/scripts | grep apobec -i | xargs -I {} sh -c 'egrep "^\./[^ ]*|[^ ]*\.sh|^bash [^ ]*|^time [^ ]*| \./[^ ]*" {} -H --color=always' | sort -V | uniq | less -NSR
ls /home/ceglab8/workspace/2.phd/9.research/9.3_gene_loss/cort/scripts | grep "apobec" -i | xargs -I {} sh -c 'egrep "^\./[^ ]*|[^ ]*\.sh|^bash [^ ]*|^time [^ ]*| \./[^ ]*" {} -o' | sort -V | uniq | less

#Summarize the apobec1_gene_loss.sh script
grep "^#[0-9]" /home/ceglab8/workspace/2.phd/9.research/9.3_gene_loss/cort/scripts/apobec1_gene_loss.sh | awk '{if(length($1)=="5" || length($1)=="6") print" - "$0; else if(length($1)==7 || length($1)==8) print"   - "$0; else if(length($1)==9 || length($1)==10) print"     - "$0; else if(length($1)==11 || length($1)==12) print"       - "$0; else print$0}' | less

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#HAE ......  PC..C
#H - Histidine - CAU, CAC
#A - Alanine - GUU, GUC, GUA, GUG
#E - Glutamate - GAA, GAG
#P - Proline - CCU, CCC, CCA, CCG, 
#C - Cysteine - UGU, UGC

##################################################################################################################################################################################################################################################################################################################
#Check number SRA datasets used for each species with A1 loss

for i in `cat total_validated_losses_Upupa_epops_removed`;do i1=`ls ~/soft_links | grep "$i"`; j=`ls ~/soft_links/"$i1"/[SED]RR*.fa | awk -F "/" '{print$NF}'`; echo $i $j; done | column -t | grep "SRR.*_\|$" 


##################################################################################################################################################################################################################################################################################################################
#Update phylogenetic tree frequently 

cd ~/bird_db1/aswin/APOBEC1/tree_update




