##############################################################################################################################################################################################################################################################################################################
#                                                                                                                             Identification of apoE gene in birds
##############################################################################################################################################################################################################################################################################################################

#Theory:
 #ApoB mRNA editing has significant functional consequences on lipoprotein metabolism. 
 #ApoB-48, the product of this editing, lacks the LDL receptor–binding domain of apoB, which is located at the carboxyl-terminal half of apoB-100, and therefore cannot bind to the LDL receptor. 
 #However, unlike those containing apoB-100, lipoproteins containing apoB-48 possess multiple copies of apoE, which serves as a ligand to mediate the clearance of the particles through both the LDL receptor and the LDL receptor–related protein.
 #Even though single copies of apoB-100 and apoE bind to the LDL receptor with similar affinities, lipoproteins with multiple copies of apoE bind to multiple LDL receptors, increasing the affinity of the lipoprotein-receptor interaction.
 #As a result, apoB-48–containing particles are cleared from the circulation more rapidly (in a matter of minutes) than are LDLs (2 to 3 days), which contain only 1 apoB-100 molecule.

##############################################################################################################################################################################################################################################################################################################
#1. Download sequences

#Main folder
cd ~/bird_db1/aswin/APOBEC1/function/APOE

#Download all APOE RNA & protein sequences from NCBI in ortholog page (these RNA sequences have UTR)
 #3 birds have apoE annotated
 #delete empty lines in fasta
sed -r '/^\s*$/d' APOE_refseq_protein.fa -i
sed -r '/^\s*$/d' APOE_refseq_transcript.fa -i

#Perform align have a quick look
mafft --auto --reorder --quiet APOE_refseq_protein.fa > APOE_refseq_protein.aln

#Download CDS sequence from NCBI using datasets
datasets download gene symbol apoE --ortholog aves --include rna,cds,product-report
datasets download gene symbol apoE --ortholog turtles --include rna,cds,product-report
datasets download gene symbol apoE --ortholog alligators --include rna,cds,product-report
datasets download gene symbol apoE --ortholog Iguania --include rna,cds,product-report
datasets download gene symbol apoE --ortholog snakes --include rna,cds,product-report
datasets download gene symbol apoE --ortholog mammals --include rna,cds,product-report
cat APOE_*cds.fa > APOE_all_cds.fa
myfasta -comb APOE_all_cds.fa | sed 's/GeneID=.*//g' | tr -d "][" | sed 's/:[^ ]\+//g' | sed 's/organism=//g' | sed 's/ /_/g' | sed 's/_$//g' > APOE_all_cds_formated.fa

mkdir ~/bird_db1/aswin/APOBEC1/function/APOE/blast
cd ~/bird_db1/aswin/APOBEC1/function/APOE/blast

##############################################################################################################################################################################################################################################################################################################
#2. Homology search in bird genomes using BLAST

#Run normal genome blast using all queries (148m58.620s)
while read i
do
j=`echo $i | awk '{print$1}'`
gn=`ls ~/soft_links/"$j"/genome/GC*.fna`
echo ">"$j $gn
blastn -task blastn -db $gn -query APOE_all_cds_formated.fa -num_threads 4 -outfmt \
"6 qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand qseq sseq" \
| sed '1i Query\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\tQuery_sequence\tSubject_sequence\n' > $j"_blast.out"
unset j gn
cd ~/bird_db1/aswin/APOBEC1/function/APOE/blast
done < <(grep -if /home/neo/bird_db1/aswin/APOBEC1/tree_update/list_of_species_to_keep_upupa_epops_removed ~/bird_db1/aswin/database_details/all_genome_paths)

#print top 20 best hits since the query is to big make a whole summary
for i in `ls *_blast.out`
do
echo ">"$i
awk '!($NF="")' $i | awk '!($NF="")' | sort -k4,4nr | sed '1i Query Subject Query_length Alignment_length Q_start Q_end S_start S_end E_value Bit_score Raw_score %_Query_covered_per_sub %_Query_covered_per_hsp %_ident Matches Mismatches Gaps Strand' | column -t | head -20
done > blast_top_20_hits

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#CGBLAST using selected exon-wise sequences

mkdir -p ~/bird_db1/aswin/APOBEC1/function/APOE/exonwise_sequences
cd ~/bird_db1/aswin/APOBEC1/function/APOE/exonwise_sequences

#Manually download exon-wise sequences from ensembl

#remove empty spaces
sed 's/ //g' *.fa -i
#make all fasta file single lined fasta
ls *.fa | xargs -n1 sh -c 'echo $0; myfasta -icomb $0'

mkdir -p ~/bird_db1/aswin/APOBEC1/function/APOE/cgblast

#Run gblast (154m52.768s)
cd ~/bird_db1/aswin/APOBEC1/function/APOE/cgblast
while read i
do
j=`echo $i | awk '{print$1}'`
mkdir $j
gn=`find ~/soft_links/"$j"/genome/ -maxdepth 1 -mindepth 1 -name "GC*.fna"`
gff=`find ~/soft_links/"$j"/genome/ -maxdepth 1 -mindepth 1 -name "GC*.gff"`
echo -e "\n>"$j | GREP_COLORS="mt=01;03;04;35" grep ".*"
#rm -r $j/*
cd $j
for k in `ls ~/bird_db1/aswin/APOBEC1/function/APOE/exonwise_sequences/APOE_*.fa`
do
k1=`echo $k | awk -F "/" '{print$NF}' | cut -f2- -d "_" | sed 's/\.fa//g'`
mkdir $k1
cd $k1
echo "  |- "$k1
if [[ -z $gff ]]
then
gblast_short $gn $k -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank
else
gblast_short $gn $k -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank $gff
fi
unset k1
cd ../
done
unset j gn gff k
cd ~/bird_db1/aswin/APOBEC1/function/APOE/cgblast
done < <(grep -if /home/neo/bird_db1/aswin/APOBEC1/tree_update/list_of_species_to_keep_upupa_epops_removed ~/bird_db1/aswin/database_details/all_genome_paths)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#summary of cgblast

cd ~/bird_db1/aswin/APOBEC1/function/APOE/cgblast
for i in `ls | grep -v "summary"`
do
#echo $i
i1=`echo $i | cut -f1,2 -d "_" | sed 's/ /_/g'`
t=`grep $i1 ~/bird_db1/aswin/taxonomy/orders_all_birds | awk '{print$2}'`
cd $i
for j in `ls`
do
l=`cat $j/best_hits | wc -l`
if [[ $l -gt 1 ]]
then
#echo "- "$j
cat $j/best_hits | grep -v "^Query" | sed "s/^/$i $t $j /g"
else :
fi
unset l
done
unset i1 t j
cd ~/bird_db1/aswin/APOBEC1/function/APOE/cgblast
done | sed '1i Subject_species Subject_Taxonomy Query_species Query Subject Query_length Alignment_length Q_start Q_end S_start S_end E_value Bit_score Raw_score %_Query_covered_per_sub %_Query_covered_per_hsp %_ident Matches Mismatches Gaps Strand' | column -t > summary

################################################################################################################################################################################################################################################################################
#2. Check Domain info

#Check Domain info of the selected APOE sequences
ls *.fa | xargs -n1 bash -c 'echo ">"$0; transeq <(grep -v ">" $0) -auto -stdout | grep -v ">"' | sed 's/*/X/g' | myfasta -comb > selected_seq.aa

scp apoE_* neo@172.30.1.174:~/bird_db1/aswin/APOBEC1/function/APOE/domain_search

cd ~/bird_db1/aswin/APOBEC1/function/APOE/domain_search
awk -F "\t" '!($NF="")' OFS="\t" apoE_hitdata.txt > apoE_hitdata_formatted.txt

for i in `ls ~/bird_db1/aswin/APOBEC1/function/APOE/exonwise_sequences/*.fa`
do
j=`echo $i | awk -F "/" '{print$NF}' | sed 's/\.fa//g'`
l=`myfasta -exdom $i apoE_hitdata_formatted.txt -specific | wc -l`
if [[ $l -gt 1 ]]
then
myfasta -exdom $i apoE_hitdata_formatted.txt -specific | grep -v "%_domain_overlap" | sed 's/^/- /g' | awk -v o="$j" '{$1=o}1'
else 
echo $j
fi
unset o j l
done | sed '1i Species Exon Exon_length Domain %_domain_overlap %_exon_of_domain_start %_exon_of_domain_end Exon_start Exon_end Domain_start Domain_end Overlap_start Overlap_end Overlap_length E-value Bit-score Completeness PSSM-ID Accession Superfamily' | column -t > domain_summary

#Statistics on summary
less summary | awk '{print$6}' | grep -v "A" | awk NF | sort | ministat

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Check how much of domain containing exon is giving hit in bird genomes 

cd ~/bird_db1/aswin/APOBEC1/function/APOE/domain_search
cp ../cgblast/summary cgblast_summary

while read i
do
i1=`echo $i | awk '{print$1}' | sed 's/APOE_//g'`
i2=`echo $i | awk '{print$2}'`
i3=`echo $i | awk '{print$14,$5,$6,$7}'`
awk -iinplace -v o="$i1" -v p="$i2" -v q="$i3" '{if($3==o && $4==p) print$0,q; else print$0}' cgblast_summary
unset i1 i2 i3
done < <(sed 1d domain_summary)

#Sort based on best hitting exons
cat cgblast_summary | awk 'NR==1 {$22="Overlap_length"; $23="%_domain_overlap"; $24="%_exon_of_domain_start"; $25="%_exon_of_domain_end"}1' | awk '{if($22=="") print$0,"- - - -"; else print$0}' \
| awk 'NR>1{print$1,$2,$3,$4,$6,$7,$8,$9,$22,$23,($8/$6)*100,($9/$6)*100,$24,$25,$12,$13,$16,$17,$18,$19,$20,$21}' | sort -k6,6 -nr \
 | sed '1i Subject_species Subject_Taxonomy Query_species Query Query_length Alignment_length Q_start Q_end Overlap_length %_domain_overlap %_exon_hit_%_start exon_hit_end %_exon_of_domain_start %_exon_of_domain_end E_value Bit_score %_Query_covered_per_hsp %_ident Matches Mismatches Gaps Strand' \
 | column -t > cgblast_domain_compare

#Exons containing domains are giving hits but their length & e-values are very poor

################################################################################################################################################################################################################################################################################
#4. Explore data

#Print single best hit per exon per species, sorted based on exon number
cd ~/bird_db1/aswin/APOBEC1/function/APOE/cgblast
while read i
do
j=$(awk -v a="$i" '$1==a {print$4,$7,$12}' summary | sed 's/exon_//g' | sort -k1,1n -k3,3g | awk '!a[$1]++' | awk '{print$1"("$2,$3")"}' OFS="," | paste -s -d " " | awk '{result=""; max_columns=4; column_index=1; for (i=1; i<=max_columns; i++) {if (column_index<=NF && $column_index~"^"i) {result=result $column_index " "; column_index++} else {result=result "- "}} print result}')
echo $i $j
unset j
done < <(grep -v "^Subject_species" summary | awk '{print$1}' | sort -u) | column -t > best_single_exon_hits_per_species

#Print single best hit per exon per species, sorted based on lowest evalue
cd ~/bird_db1/aswin/APOBEC1/function/APOE/cgblast
while read i
do
j=$(awk -v a="$i" '$1==a {print$4,$7,$12,$13}' summary | sed 's/exon_//g' | sort -k3,3g | awk '!a[$1]++' | awk '{print$1"("$2,$3,$4")"}' OFS="," | paste -s -d " ")
echo $i $j
unset j
done < <(grep -v "^Subject_species" summary | awk '{print$1}' | sort -u) | column -t > best_single_exon_hits_per_species_evalue_sorted

cat best_single_exon_hits_per_species_evalue_sorted | tr ",()" " " | sort -k5nr | column -t | less

#Map single best hit per exon per species onto phylogenetic tree
ptree.sh /home/neo/bird_db1/aswin/APOBEC1/tree_update/list_of_species_to_keep_upupa_epops_removed.nwk | sed 's/\b \b/_/g' > map_best_single_exon_hits_per_species_on_tree
while read i
do
i1=$(grep "$i" map_best_single_exon_hits_per_species_on_tree |  awk '{print$NF}')
i2=$(grep "$i" best_single_exon_hits_per_species)
sed "s/$i1/$i2/g" map_best_single_exon_hits_per_species_on_tree -i
unset i1
done < <(awk '{print$1}' best_single_exon_hits_per_species | cut -f1,2 -d "_")

################################################################################################################################################################################################################################################################################
#CLUSTERING

#Manually downloaded annotated APOE protein sequences from NCBI
cd ~/bird_db1/aswin/APOBEC1/function/APOE/Clustering

#The sequences contain a special character ("\r") at the end as windows format, so Convert Windows to Unix format
dos2unix APOE_protein_*.fa

#Combine the protein sequences an label the group name
for i in $(ls APOE_protein_*.fa)
do
n=$(echo $i | sed 's/APOE_protein_//g' | sed 's/\.fa//g')
myfasta -comb $i | sed -e 's/ /_/g' -e '/^>/ s/\[//g' -e '/^>/ s/\]//g' | sed "/^>/ s/^>/>$n\_/g"
unset n
done > APOE_tetrapoda.fa

#Clustering using CLANS online

#Open in clans.jar (Run in ceglab8, due to more cores)
java -Xmx24G -jar clans.jar -load ~/APOE.clans

#Check & compare domain info of all tetrapoda APOE
cp APOE_tetrapoda.fa /home/neo/bird_db1/aswin/APOBEC1/function/APOE/domain_search
cd /home/neo/bird_db1/aswin/APOBEC1/function/APOE/domain_search

##############################################################################################################################################################################################################################################################################################################
#TOGA

#Identify & extract the chromosome/scaffold that contains the one-to-one ortholog of APOE in alll birds for chain alignments
#This can be dont based on how reliable is the hit and synteny

#Take a quick view at how reliable are each exon hits & what are their synteny
cd ~/bird_db1/aswin/APOBEC1/function/APOE/cgblast
time while read i
do
cd $i
echo "#>"$i | GREP_COLORS="mt=01;04;07;33" grep ".*" --color=always
 find . -name "best_hits" | grep -q . && find . -name "best_hits" | xargs -n1 bash -c 'paste <(echo $0 | cut -f2 -d "/" | sed "s/^/>/g") <(grep -v "%_Query_covered_per_hsp" $0)' | sed '/^>/! s/^/- /g' \
  | sed '1i Species Query Subject Query_length Alignment_length Q_start Q_end S_start S_end E_value Bit_score Raw_score %_Query_covered_per_sub %_Query_covered_per_hsp %_ident Matches Mismatches Gaps Strand' \
  | column -t | sed 's/^/  /g' | GREP_COLORS="mt=04;38" grep -a "Species.*%_Query_covered_per_hsp.*\|$" --color=always | sed 's/^/  /g' | sed '1i @@- BEST HITS' | sed 's/@@/  /g' | GREP_COLORS="mt=01;04;07;35" grep "\- BEST HITS\|$" --color=always
 yes "-" | head -300 | paste -s -d ""
 find . -name "Exon_alignments" | grep -q . && find . -name "Exon_alignments" | xargs -n1 sh -c 'echo $0; cat $0; printf "\n"' | sed 's!^./!@!g' | awk -v RS="@" '/extracted_fla/ {print"- "$0}' \
  | sed -zE 's/([ \t\f\v\r]*\n){3,}/\n\n/g' | sed '/Exon_alignments/! s/^/   /g' | GREP_COLORS="mt=01;04;38" grep -a "[_A-Za-z/]*Exon_alignments\|$" --color=always | sed 's/\/Exon_alignments//g' \
  | sed 's/^/  /g' | sed '1i @@- EXON ALIGNMENTS' | sed 's/@@/  /g' | GREP_COLORS="mt=01;04;07;35" grep "\- EXON ALIGNMENTS\|$" --color=always | sed -e '/>exon_/ {n; /^[[:space:]]*$/d;}' -e '${/^[[:space:]]*$/d}'
 yes "-" | head -300 | paste -s -d ""
 find . -name "gblast_edited_query_covered" | grep -q . && find . -name "gblast_edited_query_covered" | xargs -n1 sh -c 'echo $0; cat $0; printf "\n"' | sed 's!^./!@!g' | awk -v RS="@" '/exon_/ {print"- "$0}' | sed -zE 's/([ \t\f\v\r]*\n){3,}/\n\n/g' \
  | sed '/gblast_edited_query_covered/! s/^/   /g' | GREP_COLORS="mt=01;04;38" grep -a "[_A-Za-z/]*gblast_edited_query_covered\|$" --color=always |  sed 's/\/consensus_qc\/gblast_edited_query_covered//g' \
  | sed 's/^/    /g' | sed '1i @@- % QUERY COVERED' | sed 's/@@/  /g' | GREP_COLORS="mt=01;04;07;35" grep "\- % QUERY COVERED\|$" --color=always
 yes "-" | head -300 | paste -s -d ""
 find . -name "Synteny" | grep -q . && find . -name "Synteny" | xargs -n1 sh -c 'echo $0; cat $0; printf "\n"' | sed 's!^./!@!g' | awk -v RS="@" '/exon/ {print"- "$0}' \
  | sed -zE 's/([ \t\f\v\r]*\n){3,}/\n\n/g' | sed '/Synteny/! s/^/   /g' | GREP_COLORS="mt=01;04;38" grep -a "[_A-Za-z/]*Synteny\|$" --color=always | sed 's/\/Synteny//g' | sed '1i @@- SYNTENY' \
  | sed 's/^/    /g' | sed 's/@@/  /g' | GREP_COLORS="mt=01;04;07;35" grep "\- SYNTENY\|$" --color=always
 yes "-" | head -300 | paste -s -d ""
 find . -name "chromview" | grep -q . && find . -name "chromview" | xargs -n1 sh -c 'echo $0; cat $0; printf "\n"' | sed 's!^./!@!g' | awk -v RS="@" '/chromview/ {print"- "$0}' | sed -zE 's/([ \t\f\v\r]*\n){3,}/\n\n/g' \
  | sed '/chromview/! s/^/   /g' | GREP_COLORS="mt=01;04;38" grep -a "[_A-Za-z/]*chromview\|$" --color=always | sed 's/\/chromview//g' | sed 's/^/  /g' | sed '1i @@- CHROMVIEW' | sed 's/@@/  /g' | GREP_COLORS="mt=01;04;07;35" grep "\- CHROMVIEW\|$" --color=always
 yes "-" | head -300 | paste -s -d ""
#find . -name "splice_sites" | xargs -n1 sh -c 'echo $0; cat $0; printf "\n"'
yes "#" | head -300 | paste -s -d ""
printf "\n"
cd ../
done < <(grep -v "^Subject_species" summary | awk '{print$1}' | sort -u) > glimpse_main_results

#Check how many species gave hits in syntenic locations
egrep -ai "#>|nectin|tomm40|tom40|apoc|cxcr3|bcl3|lip3|ceacam" glimpse_main_results | less -SR
 
##############################################################################################################################################################################################################################################################################################################
#Check ApoE manually in chicken

#New genome came in which syntenic genes are annotated (mentioned in description)
#ApoE is annotated as uncharacterized, this sequence as query nr blast gives hits to bacteria & a truncated ApoE of a bird.
#Download exon-wise sequence of ApoE, ApoC2, TOMM40 from latest chicken annotation
#Another study found 100 new genomes previously missing from chicken whcih included ApoE-like protein which is located in the  conserved synteny

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Cgblast

#Chicken syntenic genes against chicken ensembl genome
genome="/media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Gallus_gallus/genome/ensembl_genome/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa"
gff="/media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Gallus_gallus/genome/ensembl_genome/Gallus_gallus.GRCg6a.97.gff"
cd ~/soft_links/Gallus_gallus/aswin/APOE/ApoE_Gallus_gallus
gblast_short $genome ApoE_Gallus_gallus.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank $gff
cd ~/soft_links/Gallus_gallus/aswin/APOE/ApoC2_Gallus_gallus
gblast_short $genome ApoC2_Gallus_gallus.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank $gff
cd ~/soft_links/Gallus_gallus/aswin/APOE/TOMM40_Gallus_gallus
gblast_short $genome TOMM40_Gallus_gallus.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank $gff
unset genome gff

#Get bed files
for i in $(find . -name "*hit.bed"); do j=$(echo $i | cut -f3 -d "/" | sed 's/_hit.bed//g'); cat $i | awk -v a="$j" '{print$1,$2,$3,a"_"$4,$5,$6}' OFS="\t"; done > ApoE_syntenic_region.bed

cd ~/bird_db1/aswin/APOBEC1/function/APOE/cgblast/Camarhynchus_parvulus/test_Myiozetetes_cayanensis
time gblast_short ~/soft_links/Camarhynchus_parvulus/genome/GCF_901933205.1_STF_HiC_genomic.fna ~/bird_db1/aswin/APOBEC1/function/APOE/exonwise_sequences/APOE_Myiozetetes_cayanensis.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank ~/soft_links/Camarhynchus_parvulus/genome/GCF_901933205.1_STF_HiC_genomic.gff

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Long read

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#BAM files already exists for chicken: View it in IGV

#PacBio
cd /mnt/disk4/Gallus_gallus_PacBio_bam/Gal6_pacbio/aswin/APOE

#Nqnopore
cd /mnt/disk4/Gallus_gallus_PacBio_bam/Gal6_nanopore/aswin

#NOTE: Long reads shows coverage drog at the exonic region of apoE, to confirm check the genome & long read of most complete chicken genome published : huxu genome (Evolutionary analysis of a complete chicken genome)[https://doi.org/10.1073/pnas.2216641120]

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download huxu genome & respective annotation
 #NOTE: the annotation downloaded doesn't have names of genes, just some codes

#Download T2T chicken huxu genome
#link: https://www.dropbox.com/scl/fo/plq2tm2w9lzlk0ua1rzph/h?rlkey=l6z3rgmjs7ec9azun8nundnzl&e=1&dl=0
cd ~/soft_links/Gallus_gallus/genome/huxu_genome
makeblastdb -in chicken.v23.fa -out chicken.v23.fa -dbtype nucl
samtools faidx chicken.v23.fa
sortBed -i chicken.v23.gff | gff2bed > chicken.v23.gff.bed

#Cgblast
cd ~/soft_links/Gallus_gallus/aswin/APOE/huxu_genome
mkdir ApoC2_Gallus_gallus TOMM40_Gallus_gallus  ApoE_Gallus_gallus
cp ../ApoC2_Gallus_gallus/ApoC2_Gallus_gallus.fa ApoC2_Gallus_gallus/
cp ../TOMM40_Gallus_gallus/TOMM40_Gallus_gallus.fa TOMM40_Gallus_gallus/
cp ../ApoE_Gallus_gallus/ApoE_Gallus_gallus.fa ApoE_Gallus_gallus/

genome="/media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Gallus_gallus/genome/huxu_genome/chicken.v23.fa"
gff="/media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Gallus_gallus/genome/huxu_genome/chicken.v23.gff"
for i in $(ls -d */ | tr -d "/")
do
cd $i
q=$(ls | grep "Gallus_gallus.fa")
gblast_short $genome $q -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank $gff
cd ../
done
find . -name "best_hits" | xargs -n1 sh -c 'i=$(echo $0 | cut -f2 -d "/" | sed "s/$/_/g"); cat $0 | grep -v "Query" | awk -v a="$i" "{print\$2,\$7,\$8,a\$1,1,\$NF}" OFS="\t"' > ApoE_syntenic_region.bed 

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download PacBio bam files
cd /media/aswin/Long_read/Gallus_gallus_PacBio_bam

#Download from web interface in NCBI SRA Run Selector: https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15420785/SRR15420785 or using prefetch, since aws is not working 
time prefetch --max-size 100000000 SRR15420785
vdb-validate SRR15420785
fasterq-dump --fasta SRR15420785.sra

 #sam-dump --report SRR15420785
 #time sam-dump.3.0.0 -u SRR15420785 > SRR15420785.sam
 #time samtools view -bS SRR15420785.sam > SRR15420785.bam
 #time samtools index SRR15420785.bam SRR15420785.bam.bai
 #SAM file extracted shas some problem, when viewed in IGV against huxu genome, it shows no alignments! hence extract fasta & map it

time prefetch --max-size 100000000 SRR15420786
vdb-validate SRR15420786
fasterq-dump --fasta SRR15420786.sra

time prefetch --max-size 100000000 SRR15420787
vdb-validate SRR15420787
fasterq-dump --fasta SRR15420787.sra

time cat SRR15420785.fasta SRR15420786.fasta  SRR15420787.fasta > combined_pacbio.fa

#Map using minimap 2 (26m43.855s)
time minimap2 -t 32 -ax map-pb chicken.v23.fa combined_pacbio.fa > combined_pacbio.sam
#Convert sam to bam (2m55.689s)
nohup time /media/aswin/programs/samtools-1.21/samtools view -bS combined_pacbio.sam --threads 32 > combined_pacbio.bam 2> stdout_samtools_view &
#sorting (6m10.336s)
nohup time /media/aswin/programs/samtools-1.21/samtools sort combined_pacbio.bam --threads 32 -o combined_pacbio_sorted.bam 2> stdout_samtools_sort &
rm combined_pacbio.sam combined_pacbio.bam
#index the sorted bam file ()
time samtools index combined_pacbio_sorted.bam combined_pacbio_sorted.bam.bai
#View in IGV
scp neo@172.30.1.174:/media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Gallus_gallus/aswin/APOE/huxu_genome/ApoE_syntenic_region.bed .

#NOTE: PacBio reads mapped to huxu genome shows overall very poor coverage & especially at across ApoE region
#Search ApoE sequence directly in PacBio reads, since mapping depends on the reference assembly.
#How to search directly in PacBio reads
cd /media/aswin/Long_read/Gallus_gallus_PacBio_bam
time makeblastdb -in combined_pacbio.fa -out combined_pacbio.fa -dbtype nucl        #4m30.022s

#Download Nanopore
mkdir /media/aswin/Long_read/Gallus_gallus_Nanopore_db/huxu_genome
cd /media/aswin/Long_read/Gallus_gallus_Nanopore_db/huxu_genome

#time /home/ceglab25/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab25/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR154/042/SRR15421342/SRR15421342_1.fastq.gz .
#time /home/ceglab25/.aspera/connect/bin/ascp -i /home/ceglab25/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR154/042/SRR15421342/SRR15421342_1.fastq.gz .
#NOTE: In ceglab25 only prefetch.3.0.0 works, eg: time prefetch.3.0.0 --max-size 100000000 SRR15421342
time prefetch --max-size 100000000 SRR15421342
time prefetch --max-size 100000000 SRR15421343
time prefetch --max-size 100000000 SRR15421344
time prefetch --max-size 100000000 SRR15421345
time prefetch --max-size 100000000 SRR15421346

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#SRA blast

cd /home/neo/soft_links/Gallus_gallus/aswin/APOE/ApoE_Gallus_gallus/sblast
blastn -task blastn -db /home/neo/soft_links/Gallus_gallus/SRR12103809_SRR1744081_SRR3041423_SRR3954707.fa -query ../ApoE_Gallus_gallus.fa -num_threads 4 -evalue 0.01 -outfmt 11 -out sblast.outfmt11
#Reformat to outfmt 3
blast_formatter -archive sblast.outfmt11 -outfmt 3 -line_length 280 -out sblast.outfmt3
#Reformat to sam format
blast_formatter -archive sblast.outfmt11 -outfmt "17 SQ" -out sblast.sam
#SRA blast summary
qvsra sblast.outfmt3 -t > sblast_summary
#SRA blast QC & consensus sequence making using 2 iffernt tools : to compare & quality check
sedit.sh sblast.outfmt3 2>/dev/null
samtools consensus sblast.sam -f fasta --mode simple -l 10000 | sed 's/Query/exon/g' > sblast.sam.fa

##############################################################################################################################################################################################################################################################################################################
#DRAFT  SCRIPTS
##############################################################################################################################################################################################################################################################################################################

#Check cgblast results from all queries for a specific species
cd ~/bird_db1/aswin/APOBEC1/function/APOE/cgblast/Gallus_gallus
find . -name "best_hits" | xargs -n1 sh -c 'echo $0 | cut -f2 -d "/" | sed "s/^/>/g"; cat $0; printf "\n"' | column -t | awk '!a[$0]++' | sed '/^>/! s/^/- /g' | column -t > combined_best_hits
find . -name "Exon_alignments" | xargs -n1 sh -c 'echo $0; cat $0; printf "\n"' > combined_Exon_alignments
find . -name "Synteny" | xargs -n1 sh -c 'echo $0; cat $0; printf "\n"' > combined_Synteny
find . -name "splice_sites" | xargs -n1 sh -c 'echo $0; cat $0; printf "\n"' > combined splice_sites

