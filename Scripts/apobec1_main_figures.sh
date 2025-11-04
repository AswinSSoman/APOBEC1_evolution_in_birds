##################################################################################################################################################################################################################################################################################################################
#Main figures scripts
##################################################################################################################################################################################################################################################################################################################

#Update tree
mkdir ~/bird_db1/aswin/APOBEC1/main_figures
cd ~/bird_db1/aswin/APOBEC1/main_figures

cp ../species_to_exclude ../all_birds ../all_birds.nwk .
grep -ivf species_to_exclude all_birds | paste -s -d "," > keep_species
Rscript prune_tree.r all_birds.nwk keep_species final_birds.nwk

#Download species tree
less list_of_species_to_keep.nwk | tr -d ")(;:[0-9].'" | tr "," "\n"  > species_from_tree
grep -if species_from_tree list_of_species_to_keep -v | cut -f2 -d "_" > names_to_replace
ls list_of_species_to_keep species_from_tree | xargs -n1 sh -c 'echo ">"$0;grep -if names_to_replace $0 | sed "s/^/ /g"'
sed 's/Cyanoderma_ruficeps/Stachyris_ruficeps/g' list_of_species_to_keep.nwk -i
sed 's/Chloebia_gouldiae/Erythrura_gouldiae/g' list_of_species_to_keep.nwk -i

cat list_of_species_to_keep | xargs -n1 sh -c 'echo ">"$0; esearch -db taxonomy -query "$0[ORGN]" | efetch -format xml | xtract -pattern TaxaSet -element Lineage' > lineages_of_list_of_species_to_keep
sed 's/cellular organisms.*Aves; //g' lineages_of_list_of_species_to_keep | sed 's/Palaeognathae;//g' | sed 's/Neognathae;//g' | awk -v RS=">" '{$1=$1}1' | tr -d ";" | sort -k2,2 | column -t > lineages_of_list_of_species_to_keep_filtered
sed 's/\t/\n/g' lineages_of_list_of_species_to_keep | sed 's/cellular organisms.*Aves; //g' | sed 's/Palaeognathae;//g' | sed 's/Neognathae;//g' > lineages_of_list_of_species_to_keep_with_multiple_lineages

##################################################################################################################################################################################################################################################################################################################
#Figure 1 : Loss events

cd ~/bird_db1/aswin/APOBEC1/main_figures/loss_events
cp ~/bird_db1/aswin/APOBEC1/main_figures/tree/list_of_species_to_keep_figtree_edited_expanded.nwk .
cp ~/bird_db1/aswin/APOBEC1/main_figures/tree/list_of_species_to_keep.nwk .
cp /home/neo/bird_db1/aswin/APOBEC1/main_figures/list_of_species_to_keep .

#Prune the tree : remove upupa epops since the loss can be due to incomplete assembly (gene is near the end of assembly)
grep -v "Upupa_epops" list_of_species_to_keep | paste -s -d "," > list_of_species_to_keep_upupa_epops_removed.csv
Rscript /home/neo/bird_db1/aswin/APOBEC1/main_figures/tree/prune_tree.r list_of_species_to_keep.nwk list_of_species_to_keep_upupa_epops_removed.csv list_of_species_to_keep_upupa_epops_removed.nwk

#Edit the tree in figtree : edit species with loss & collpase a clade if possible
figtree list_of_species_to_keep_upupa_epops_removed.nwk

cp list_of_species_to_keep_upupa_epops_removed_collapsed_labelled.nwk list_of_species_to_keep_upupa_epops_removed_collapsed_labelled_without_underscore.nwk
for i in $(cat list_of_species_to_keep_upupa_epops_removed_collapsed_labelled | tr " " "_")
do
i1=$(echo $i | sed 's/_/ /g') 
sed "s/$i/$i1/g" list_of_species_to_keep_upupa_epops_removed_collapsed_labelled_without_underscore.nwk -i
unset i1
done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Add SRA box for specific species

#Cathartes_aura
s="Cathartes_aura"
r="XM_027449866.2"

#TOGA mutation
awk -v RS="GENE" -v a="$r" '$0~a {print"GENE"$0}' ~/bird_db1/aswin/APOBEC1/TOGA/$s/toga_output/inact_mut_data.txt | awk NF
awk -v RS=">" -v a="$r" '$0~a {print">"$0}' ~/bird_db1/aswin/APOBEC1/TOGA/$s/toga_output/codon.fasta | awk NF

#Manual alignment
exmalign2 toga_rna-XM_027449866.2_exon_wise_consensus.fa ~/soft_links/$s/aswin/APOBEC1/2nd_gblast/manual_search/final_manual_sequence.fa -m 280
exmalign2 toga_rna-XM_027449866.2_exon_wise_consensus.fa ~/soft_links/$s/aswin/APOBEC1/2nd_gblast/APOBEC1_Mesitornis_unicolor.fa -m 280
exmalign2 ~/soft_links/$s/aswin/APOBEC1/2nd_gblast/manual_search/final_manual_sequence.fa ~/soft_links/$s/aswin/APOBEC1/2nd_gblast/APOBEC1_Mesitornis_unicolor.fa -m 280

#Get validated A1 birds to compare & locate inactivating mutations
cd ~/soft_links/$s/aswin/APOBEC1/2nd_gblast/manual_search
cp ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.fa .
cp ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.aa .

#Manually aligned in aliview using muscle & refine it by removing 3 sequences with noticable gaps, especially at first half of the gene
aliview total_validated_queries.fa
awk -v RS=">" '!/Malurus_cyaneus_samueli|Pelecanus_crispus|Zonotrichia_albicollis/ {print">"$0}' total_validated_queries.fa | awk NF | grep -v "^>$" > total_validated_queries_refined.fa
awk -v RS=">" '!/Malurus_cyaneus_samueli|Pelecanus_crispus|Zonotrichia_albicollis/ {print">"$0}' total_validated_queries.aa | awk NF | grep -v "^>$" > total_validated_queries_refined.aa

#Run amino acid guided DNA alignment
#mafft --maxiterate 1000 --localpair --reorder --quiet total_validated_queries_refined.aa > total_validated_queries_refined.aa.mafft.aln
#muscle -in total_validated_queries_refined.aa -maxiters 1000 > total_validated_queries_refined.aa.muscle.aln
cat total_validated_queries_refined.fa <(grep -v ">" final_manual_sequence.fa | paste -s -d "" | sed '1i >TEST') > final_manual_sequence_with_total_validated_queries_refined.fa
cat total_validated_queries_refined.fa <(grep -v ">" ~/bird_db1/aswin/APOBEC1/TOGA/Cathartes_aura/toga_rna-XM_027449866.2_exon_wise_consensus.fa | paste -s -d "" | sed '1i >TEST') > final_manual_sequence_with_total_validated_queries_refined.fa
#Load the fasta in aliview and align using muscle not mafft, and remove all species showing gaps, then I got the same inactivation mutation given by TOGA
aliview final_manual_sequence_with_total_validated_queries_refined.fa
#Aligning in commandline muscle also shows same inactivation mutation given by TOGA
muscle -in final_manual_sequence_with_total_validated_queries_refined.fa | alv -
#Codon alignment ()Make nice alignment but not showing same inactivation mutation given by TOGA. even when using muscle as aligner
conda activate base  
translatorx -i final_manual_sequence_with_total_validated_queries_refined.fa -p M

#SRA blast
blastn -task blastn -db ~/soft_links/Cathartes_aura/SRR954275_76.fa -query ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Aquila_chrysaetos_chrysaetos.fa -evalue 0.01 -outfmt 11 -out APOBEC1_Aquila_chrysaetos_chrysaetos_against_SRR954275_76.outfmt11
blast_formatter -archive APOBEC1_Aquila_chrysaetos_chrysaetos_against_SRR954275_76.outfmt11 -outfmt 3 -line_length 280 -out APOBEC1_Aquila_chrysaetos_chrysaetos_against_SRR954275_76.outfmt3
blast_formatter -archive APOBEC1_Aquila_chrysaetos_chrysaetos_against_SRR954275_76.outfmt11 -outfmt 6 -out APOBEC1_Aquila_chrysaetos_chrysaetos_against_SRR954275_76.outfmt6
blastn -task blastn -db ~/soft_links/Cathartes_aura/SRR954275_76.fa -query /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/APOBEC1_Aquila_chrysaetos_chrysaetos_flanking.fa -evalue 0.01 -outfmt 11 -out APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR954275_76.outfmt11
blast_formatter -archive APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR954275_76.outfmt11 -outfmt 3 -line_length 280 -out APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR954275_76.outfmt3
blast_formatter -archive APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR954275_76.outfmt11 -outfmt 6 -out APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR954275_76.outfmt6
blast_formatter -archive APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR954275_76.outfmt11 -outfmt '6 qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand qseq sseq' > APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR954275_76.outfmt6

#Get region containing inactivating mutation from SRA reads for figure 1
#grep 'GCACCCAGAG' APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR954275_76.outfmt6 --color=always | less -SR
awk '/GCACCCAGAG/ {print">"$2"\n"$NF}' APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR954275_76.outfmt6 | tr -d "-" > exon_2_Cathartes_aura_sra_reads.fa
#muscle -in <(ls ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*.fa | xargs -n1 bash -c 'echo $0 | cut -f8 -d "/" | sed "s/\.fa//g" | sed "s/^/>exon_/g"; seqtk subseq $0 <(echo "exon_2") | grep -v ">"') | alv -
ls ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*.fa | xargs -n1 bash -c 'echo $0 | cut -f8 -d "/" | sed "s/\.fa//g" | sed "s/^/>exon_2_/g"; seqtk subseq $0 <(echo "exon_2") | grep -v ">"' | awk -v RS=">" '!/exon_APOBEC1_Pelecanus_crispus/ {print">"$0}' | grep -v "^>$" | awk NF > exon_2_total_validated_queries_refined.fa
mafft --localpair --maxiterate 1000 --quiet --reorder <(cat exon_2_Cathartes_aura_sra_reads.fa exon_2_total_validated_queries_refined.fa) > exon_2_sra.aln

#grep 'GGAGAAGATGATTAT' APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR954275_76.outfmt6 --color=always | less -SR
awk '/GGAGAAGATGATTAT/ {print">"$2"\n"$NF}' APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR954275_76.outfmt6 | tr -d "-" > exon_4_Cathartes_aura_sra_reads.fa
#muscle -in <(ls ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*.fa | xargs -n1 bash -c 'echo $0 | cut -f8 -d "/" | sed "s/\.fa//g" | sed "s/^/>exon_/g"; seqtk subseq $0 <(echo "exon_4") | grep -v ">"') | alv -
ls ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*.fa | xargs -n1 bash -c 'echo $0 | cut -f8 -d "/" | sed "s/\.fa//g" | sed "s/^/>exon_4_/g"; seqtk subseq $0 <(echo "exon_4") | grep -v ">"' | awk -v RS=">" '!/exon_APOBEC1_Corvus_moneduloides/ {print">"$0}' | grep -v "^>$" | awk NF > exon_4_total_validated_queries_refined.fa
mafft --localpair --maxiterate 1000 --quiet --reorder <(cat exon_4_total_validated_queries_refined.fa exon_4_Cathartes_aura_sra_reads.fa) > exon_4_sra.aln

GGAGAAGAGGATTAT
GGAGAAGATGATTAT

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Gymnogyps_californianus

s="Gymnogyps_californianus"
r="XM_027449866.2"
awk -v RS="GENE" -v a="$r" '$0~a {print"GENE"$0}' ~/bird_db1/aswin/APOBEC1/TOGA/$s/toga_output/inact_mut_data.txt | awk NF
awk -v RS=">" -v a="$r" '$0~a {print">"$0}' ~/bird_db1/aswin/APOBEC1/TOGA/$s/toga_output/codon.fasta | awk NF

#Manual alignment
exmalign2 toga_rna-XM_027449866.2_exon_wise_consensus.fa ~/soft_links/$s/aswin/APOBEC1/2nd_gblast/manual_search/final_manual_sequence.fa -m 280
exmalign2 toga_rna-XM_027449866.2_exon_wise_consensus.fa ~/soft_links/$s/aswin/APOBEC1/2nd_gblast/APOBEC1_Mesitornis_unicolor.fa -m 280
exmalign2 ~/soft_links/$s/aswin/APOBEC1/2nd_gblast/manual_search/final_manual_sequence.fa ~/soft_links/$s/aswin/APOBEC1/2nd_gblast/APOBEC1_Mesitornis_unicolor.fa -m 280

#Get validated A1 birds to compare & locate inactivating mutations
cd ~/soft_links/$s/aswin/APOBEC1/2nd_gblast/manual_search
cp ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.fa .

#Manually aligned in aliview using muscle & refine it by removing 3 sequences with noticable gaps, especially at first half of the gene
aliview total_validated_queries.fa
awk -v RS=">" '!/Malurus_cyaneus_samueli|Pelecanus_crispus|Zonotrichia_albicollis/ {print">"$0}' total_validated_queries.fa | awk NF | grep -v "^>$" > total_validated_queries_refined.fa

#Run amino acid guided DNA alignment
cat total_validated_queries_refined.fa <(grep -v ">" final_manual_sequence.fa | paste -s -d "" | sed '1i >TEST') > final_manual_sequence_with_total_validated_queries_refined.fa
cat total_validated_queries_refined.fa <(grep -v ">" ~/bird_db1/aswin/APOBEC1/TOGA/$s/toga_rna-XM_027449866.2_exon_wise_consensus.fa | paste -s -d "" | sed '1i >TEST') > toga_rna-XM_027449866.2_exon_wise_consensus_with_total_validated_queries_refined.fa
muscle -in final_manual_sequence_with_total_validated_queries_refined.fa | alv -
muscle -in toga_rna-XM_027449866.2_exon_wise_consensus_with_total_validated_queries_refined.fa | alv -
conda activate base  
translatorx -i toga_rna-XM_027449866.2_exon_wise_consensus_with_total_validated_queries_refined.fa

#SRA blast
time blastn -task blastn -db ~/soft_links/Gymnogyps_californianus/SRR11097121_22.fa -query /home/neo/bird_db1/aswin/APOBEC1/validated_sequences_with_flanking_regions/APOBEC1_Aquila_chrysaetos_chrysaetos_flanking.fa -evalue 0.01 -outfmt 11 -out APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR11097121_22.outfmt11
blast_formatter -archive APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR11097121_22.outfmt11 -outfmt 3 -line_length 280 -out APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR11097121_22.outfmt3
blast_formatter -archive APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR11097121_22.outfmt11 -outfmt '6 qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand qseq sseq' > APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR11097121_22.outfmt6

#Get region containing inactivating mutation from SRA reads for figure 1
#awk '/GCACCCAGAG/ {print">"$1"_sra\n"$NF}' APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR11097121_22.outfmt6 | tr -d "-" > exon_2_Gymnogyps_californianus_sra_reads
awk '/GCACCCAGAG/ {print">"$2"\n"$NF}' APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR11097121_22.outfmt6 | tr -d "-" > exon_2_Gymnogyps_californianus_sra_reads.fa
#muscle -in <(ls ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*.fa | xargs -n1 bash -c 'echo $0 | cut -f8 -d "/" | sed "s/\.fa//g" | sed "s/^/>exon_/g"; seqtk subseq $0 <(echo "exon_2") | grep -v ">"') | alv -
ls ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*.fa | xargs -n1 bash -c 'echo $0 | cut -f8 -d "/" | sed "s/\.fa//g" | sed "s/^/>exon_/g"; seqtk subseq $0 <(echo "exon_2") | grep -v ">"' | awk -v RS=">" '!/exon_APOBEC1_Pelecanus_crispus/ {print">"$0}' | grep -v "^>$" | awk NF > exon_2_total_validated_queries_refined.fa
mafft --localpair --maxiterate 1000 --quiet --reorder <(cat exon_2_Gymnogyps_californianus_sra_reads.fa exon_2_total_validated_queries_refined.fa) > exon_2_sra.aln

#grep 'GGAGAAGATGATTAT' APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR11097121_22.outfmt6 --color=always | less -SR
#grep 'GGAGAAGATGATTAT...' APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR11097121_22.outfmt6 --color=always -o | less -SR
awk '/GGAGAAGATGATTAT/ {print">"$2"\n"$NF}' APOBEC1_Aquila_chrysaetos_chrysaetos_flanking_against_SRR11097121_22.outfmt6 | tr -d "-" > exon_4_Gymnogyps_californianus_sra_reads.fa
#muscle -in <(ls ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*.fa | xargs -n1 bash -c 'echo $0 | cut -f8 -d "/" | sed "s/\.fa//g" | sed "s/^/>exon_/g"; seqtk subseq $0 <(echo "exon_4") | grep -v ">"') | alv -
ls ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*.fa | xargs -n1 bash -c 'echo $0 | cut -f8 -d "/" | sed "s/\.fa//g" | sed "s/^/>exon_4_/g"; seqtk subseq $0 <(echo "exon_4") | grep -v ">"' | awk -v RS=">" '!/exon_APOBEC1_Corvus_moneduloides/ {print">"$0}' | grep -v "^>$" | awk NF > exon_4_total_validated_queries_refined.fa
mafft --localpair --maxiterate 1000 --quiet --reorder <(cat exon_4_total_validated_queries_refined.fa exon_4_Gymnogyps_californianus_sra_reads.fa) > exon_4_sra.aln
myfasta -comb exon_4_sra.aln | GREP_COLORS="mt=01;32" grep -iB1 "GGAGAAGATGATTAT"  --color=always
myfasta -comb exon_4_sra.aln | GREP_COLORS="mt=01;31" grep -iB1 "GGAGAAGATGATTATTGA"  --color=always

##################################################################################################################################################################################################################################################################################################################
#Catalytic sites

mkdir ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_losses_Upupa_epops_removed

cd ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_losses_Upupa_epops_removed
while read i
do
i1=$(echo $i | awk '{print$1}')
i2=$(echo $i | awk '{print$2}')
i3=$(echo $i2 | sed 's/.\{3\}/& /g')
i4=$(transeq <(echo $i2) --auto --stdout --clean | grep -v ">" | paste -s -d "" | sed 's/.\{1\}/&  /g' | sed 's/^/ /g')
echo -e $i1 $i3"\n -"$i4
unset i1 i2 i3 i4
done < <(awk -v RS=">" '{$1=$1}1' ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.fa | awk NF) | column -t | sed '/^-/ s/^/ /g' | GREP_COLORS="mt=07;31" egrep "  H    A    E  |  P    C    .    .    C  " -z --color=always > active_site_codon_aa_view

#Get exon_3 of birds with intact A1
ls ~/bird_db1/aswin/APOBEC1/validated_sequences/*.fa | xargs -n1 bash -c 'echo $0 | cut -f8 -d "/" | cut -f2- -d "_" | sed "s/\.fa//g" | sed "s/^/>intact_exon_3_/g"; seqtk subseq $0 <(echo exon_3) | grep -v ">"' > exon_3_birds_A1_intact.fa

#exon_3 of birds with A1 loss
#check how many consensus is made for birds with A1 loss
for i in `grep -if ~/bird_db1/aswin/APOBEC1/total_validated_losses_Upupa_epops_removed <(ls -d ~/soft_links/*/ | sed 's!/$!!g')`
do
cd "$i"/aswin/APOBEC1/2nd_gblast/
j1=`ls -d */ | grep manual_search`
j2=`find . -name "final_consensus.fa"`
echo $i $j2 $j1
unset j1 j2
done | column -t

#Get exon_3 of birds with A1 loss
cd ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_losses_Upupa_epops_removed
for i in `grep -if ~/bird_db1/aswin/APOBEC1/total_validated_losses_Upupa_epops_removed <(ls -d ~/soft_links/*/ | sed 's!/$!!g') | awk -F "/" '{print$NF}'`
do
i1=$(echo $i | cut -f1,2 -d "_")
echo "#"$i | GREP_COLORS="mt=01;32" grep ".*" --color=always
seqtk subseq /home/neo/soft_links/"$i"/aswin/APOBEC1/2nd_gblast/gblast_edited_consensus.fa <(echo "exon_3") | grep -v ">" | paste -s -d "" | sed '1i >gedit_exon_3' | myfasta -de > kala_exon_3
#Some species have different folder names. eg: Struthio_camelus_australis & Struthio_camelus
i2=$(find ~/bird_db1/aswin/APOBEC1/TOGA/ -maxdepth 1 -name "$i1" -type d)
seqtk subseq "$i2"/toga_*_exon_wise_consensus.fa <(echo "exon_3") | grep -v ">" | paste -s -d "" | sed '1i >toga_exon_3' | myfasta -de >> kala_exon_3
seqtk subseq /home/neo/soft_links/"$i"/aswin/APOBEC1/2nd_gblast/final_consensus.fa <(echo "exon_3") | grep -v ">" | paste -s -d "" | sed '1i >manual_exon_3' | myfasta -de >> kala_exon_3
if [[ $(grep ">" kala_exon_3 -c) > 1 ]]
then
muscle -quiet -in kala_exon_3 -diags | myfasta -comb
else
cat kala_exon_3
fi
rm kala_exon_3
unset i1 i2
done | awk NF | grep -v "^>$" | myfasta -de > all_consensus_of_exon_3_birds_A1_loss

#It looks like consensus sequence of exon_3 from gedit & manual consensus are same, but in 13 birds with loss toga exon_3 consensus lack start portion. 
#Some times we see SRA reads have zero coverage in this region, and genome sequence have small letter for this region suggesting repeats present 

for i in `sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" all_consensus_of_exon_3_birds_A1_loss | egrep "#|^-" | awk -v RS="#" '{$1=$1}1' | awk '$2~"-" {print$1}'`
do
echo "#"$i | GREP_COLORS="mt=01;32" grep ".*" --color=always
i1=$(echo $i | cut -f1,2 -d "_")
i2=$(find ~/bird_db1/aswin/APOBEC1/TOGA/ -maxdepth 1 -name "$i1" -type d)
echo $i1 $i2
cat "$i2"/All_consensus_msa.aln
unset i1 i2
done | less -SR

#Repeat find using repeatmasker have found repeats in exon 3 of many species, but not in some birds. This could be due to lack of that repeat in repbase & require denovo type of repeat identification to identify them
#Running Repeatmodeler failed giving error : 
  #No families identified.  Perhaps the database is too small
  #or contains overly fragmented sequences.
#To quickly find repeats in this region try some online web-server: https://www.girinst.org/censor/index.php (mentioned in article: https://academic.oup.com/gbe/article/13/12/evab259/6433158#323041721)

for i in `grep -if ~/bird_db1/aswin/APOBEC1/total_validated_losses_Upupa_epops_removed <(ls -d ~/soft_links/*/ | sed 's!/$!!g') | awk -F "/" '{print$NF}'`
do
i1=$(echo $i | cut -f1,2 -d "_")
echo $i | sed 's/^/>disrupted_exon_3_/g'
seqtk subseq /home/neo/soft_links/"$i"/aswin/APOBEC1/2nd_gblast/gblast_edited_consensus.fa <(echo "exon_3") | grep -v ">" | paste -s -d ""
unset i1
done | awk NF | grep -v "^>$" | myfasta -de > exon_3_birds_A1_loss.fa

#For example
cd ~/soft_links/Centrocercus_minimus/aswin/APOBEC1/2nd_gblast
bedtools getfasta -fi ../../../genome/GC*.fna -bed <(cat test.out.bed | grep exon_3 | sed 's/plus/+/g' | sed 's/minus/-/g' | awk '{if($NF=="+") print$1,$2+21,$3-20,$4,"1",$5; else print$1,$2+20,$3-21,$4,".",$5}' OFS="\t") -s -name
bedtools getfasta -fi ../../../genome/GC*.fna -bed <(cat test.out.bed | grep exon_3 | sed 's/plus/+/g' | sed 's/minus/-/g' | awk '{if($NF=="+") print$1,$2+21-1000,$3-20,$4,"1",$5; else print$1,$2+20,$3-21+1000,$4,".",$5}' OFS="\t") -s -name

#align & find mutations
cat exon_3_birds_A1_intact.fa exon_3_birds_A1_loss.fa > exon_3_birds_A1_all.fa
#Check frame
readlink -f ~/bird_db1/aswin/APOBEC1/validated_sequences/*.fa | xargs -n1 sh -c 'myfasta -codex $0 -nw | grep "exon_3" -A1' | less -S
sed '/^>/! s/^G//g' exon_3_birds_A1_intact.fa > exon_3_birds_A1_intact_inframe.fa
transeq exon_3_birds_A1_intact_inframe.fa --auto --stdout -trim | sed '/^>/ s/_1$//g' | myfasta -comb > exon_3_birds_A1_intact_inframe.aa
muscle -in exon_3_birds_A1_intact_inframe.aa > exon_3_birds_A1_intact_inframe.aln

#Run codon alignment with 3 frames and Allow 10 codons for frameshift compensation (https://hcv.lanl.gov/content/sequence/CodonAlign/codonalign.html)
#Download nt & translated alignment in fasta & table forrat
for i in $(grep ">" exon_3_birds_A1_all_codoncode_nt.aln | tr -d ">")
do
echo ">"$i
awk -v RS=">" -v a="$i" '$0~a {print">"$0}' exon_3_birds_A1_all_codoncode_nt.aln | awk NF | grep -v "^>" | paste -s -d "" | sed 's/.\{3\}/& /g' | sed 's/^/- /g'
awk -v RS=">" -v a="$i" '$0~a {print">"$0}' exon_3_birds_A1_all_codoncode_aa.aln | awk NF | grep -v "^>" | paste -s -d "" | sed 's/.\{1\}/&  /g' | sed 's/^/- /g'
done | column -t | GREP_COLORS="mt=07;32" egrep " H    A    E | P    C    .    .    C " --color=always -z  | less -SR 


##################################################################################################################################################################################################################################################################################################################
#Figure 4: Protein alignment

cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cds_alignment
sed -e '/^>/ s/(.*)//g' birds_mammals_filtered_8.aa > birds_mammals_filtered_8_with_group_symbol.aa
muscle -in birds_mammals_filtered_8_with_group_symbol.aa > birds_mammals_filtered_8_with_group_symbol_muscle.aln
mafft --maxiterate 1000 --localpair --reorder birds_mammals_filtered_8_with_group_symbol.aa > birds_mammals_filtered_8_with_group_symbol_mafft.aln

cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cds_alignment
sed -e '/^>/ s/(.*)//g' birds_mammals_filtered_8.aa > birds_mammals_filtered_8_with_group_names.aa
sed -e '/^>/ s/(.*)//g' -e '/^>/ s/B_//g' -e '/^>/ s/M_//g' birds_mammals_filtered_8.aa > birds_mammals_filtered_8_renamed.aa
muscle -in birds_mammals_filtered_8_renamed.aa > birds_mammals_filtered_8_renamed_muscle.aln
mafft --maxiterate 1000 --localpair --reorder birds_mammals_filtered_8_renamed.aa > birds_mammals_filtered_8_renamed_mafft.aln

#Alignment & create Word compatible output
~/programmes/gismo birds_mammals_filtered_8_renamed.aa -rtf
~/programmes/gismo birds_mammals_filtered_8.aa -rtf

#Best way to get publication quality alignment images with ability to mark specific residues is "jalview"

#View alignment in jalview & manually mark by selecting the position right click > selection > Edit group > Group color > Choose colouring scheme
#When mutliple variations of residues observed within birds/mammals, check if these share same physiochemical property: 
  #For easy counting number of species supporting specific nucleotide in specific position use commandline rather than manual inspection
  #Copy the whole alignment from jalview & save it into a file "jalview_alignment"
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' jalview_alignment | awk -v RS=">" '{$1=$1}1' | awk '{print $1,substr($2,23,1)}' | column -t | less
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' jalview_alignment | awk -v RS=">" '{$1=$1}1' | awk '{print $1,substr($2,24,1)}' | awk '/^B_/{print$2}' | sort | uniq -c | sort -k1,1nr | awk '{print$2"("$1")"}' | paste -s -d ","
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' jalview_alignment | awk -v RS=">" '{$1=$1}1' | awk '{print $1,substr($2,24,1)}' | awk '/^M_/{print$2}' | sort | uniq -c | sort -k1,1nr | awk '{print$2"("$1")"}' | paste -s -d ","

#Make amino acid physiochemical properties based on imgt.org
IMGT_classes_amino_acid_properties

time for p in $(seq 1 $(awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' jalview_alignment | sed -n 2p | wc | awk '{print$3-$1}'))
do
pih=$(awk -v RS=">" '/Homo_sapiens/ {print">"$0}' jalview_alignment | awk NF | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' | awk -v RS=">" '{$1=$1}1' | awk NF | awk -v a="$p" '{print $1,substr($2,0,a)}' | awk '{print$2}' | tr -d "-" | wc | awk '{print$3-$1}')
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' jalview_alignment | awk -v RS=">" '{$1=$1}1' | awk -v a="$p" '{print $1,substr($2,a,1)}' | awk '/^B_/{print$2}' | sort | uniq -c | sort -k1,1nr > temp_count_birds
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' jalview_alignment | awk -v RS=">" '{$1=$1}1' | awk -v a="$p" '{print $1,substr($2,a,1)}' | awk '/^M_/{print$2}' | sort | uniq -c | sort -k1,1nr > temp_count_mammals
c1=$(cat temp_count_birds | awk '{print$2"("$1")"}' | paste -s -d ",")
c2=$(cat temp_count_mammals | awk '{print$2"("$1")"}' | paste -s -d ",")
o=$(sed 's/^[ ]\+//g' temp_count_birds | xargs -n2 bash -c 'paste <(grep $1 temp_count_mammals | awk "{print\$2,\$1}") <(echo $0) -d " "' | awk '$3' | awk '{print$1"("$3","$2")"}' | paste -s -d ",")
if [ -z "$o" ]; then o="-"; else :; fi
#get physiochemical properties of manio acids
	pcp=$(for group in temp_count_birds temp_count_mammals
	do
	echo -n > temp_entry
	while read tcb
	do
	tcb1=$(echo $tcb | awk '{print$1}')
	tcb2=$(echo $tcb | awk '{print$2}')
	if [ "$group" == "temp_count_birds" ]
	then
	awk -v b="$tcb2" '$3==b {print$4,$5,$6,$8,$9,$10}' OFS="\n" IMGT_classes_amino_acid_properties | sed 's/([0-9]\+)//g' | awk -v a="$tcb1" '{print $0,a/43*100}' | awk '{$NF+=0}1' CONVFMT="%.1f" | awk '{print$1"("$2")"}' | paste -s -d " " >> temp_entry
	elif [ $group == "temp_count_mammals" ]
	then
	awk -v b="$tcb2" '$3==b {print$4,$5,$6,$8,$9,$10}' OFS="\n" IMGT_classes_amino_acid_properties | sed 's/([0-9]\+)//g' | awk -v a="$tcb1" '{print $0,a/32*100}' | awk '{$NF+=0}1' CONVFMT="%.1f" | awk '{print$1"("$2")"}' | paste -s -d " " >> temp_entry
	fi
	unset tcb1 tcb2
	done < $group
	#Transpose
	awk '{for(i=1;i<=NF;i++) {a[NR,i]=$i}} NF>p{p=NF} END{for(j=1;j<=p;j++) {str=a[1,j]; for(i=2;i<=NR;i++) {str=str" "a[i,j];} print str}}' temp_entry > temp_entry2
	#cat temp_entry2 | tr "()" " " | awk '{print$2+$4+$6+$8+$10+$12+$14+$16}' | uniq -c
	#Count & sum the properties across birds or mammals
	physiochemproperty=$(while read te
	do
	echo $te | tr " " "\n" | tr -d ")" | tr "(" " " | sort -k1,1 | awk '{sum[$1]+=$2} END{for (k in sum) print k, sum[k]}' | sort -k2,2nr | sed 's/ /(/g' | sed 's/$/)/g' | paste -s -d ","
	done < temp_entry2 | sed 's/positive_charged/P/g' | sed 's/negative_charged/N/g' | sed 's/uncharged/U/g' | sed 's/verylarge/VL/g' | sed 's/large/L/g' | sed 's/medium/M/g' | sed 's/verysmall/VS/g' | sed 's/small/S/g' | sed 's/,uncharged/U/g' | sed 's/nonpolar/NP/g' | sed 's/polar/P/g' \
	| sed 's/donor_and_acceptor/DA/g' | sed 's/donor/D/g' | sed 's/acceptor/A/g' | sed 's/hydrophobic/B/g' | sed 's/hydrophilic/L/g' | sed 's/neutral/N/g' | sed 's/sulfur/Sul/g' | sed 's/hydroxyl/Hyd/g' | sed 's/basic/Bas/g' | sed 's/aliphatic/Ali/g' | sed 's/aromatic/Aro/g' \
	| sed 's/amide/Ami/g' | sed 's/acidic/Aci/g' | sed 's/none/N/g' | paste -s -d " \t")
	if [ -z "$physiochemproperty" ]; then physiochemproperty="- - - - - - "; else :; fi
	echo $physiochemproperty
	unset physiochemproperty
	rm temp_entry temp_entry2
	unset tcb te
	done)
#Major physicochemial change in each amino acid position in birds w.r.t mammals
	pcpmd=$(while read pcpe
	do
	for pcpea in $(echo $pcpe | sed 's/([0-9.]\+)//g' | tr "," " " | tr " " "\n" | awk NF | sort -u | tr -d "-" | awk NF)
	do
	pcpea1=$(echo $pcpe | awk '{print$1}' | grep "\b$pcpea([^)]\+)" -o | cut -f2 -d "(" | tr -d ")")
	pcpea2=$(echo $pcpe | awk '{print$2}' | grep "\b$pcpea([^)]\+)" -o | cut -f2 -d "(" | tr -d ")")
	if [ -z $pcpea1 ]; then pcpea3=$(echo -$pcpea2); elif [ -z $pcpea2 ]; then pcpea3=$(echo $pcpea1); else pcpea3=$(calc $pcpea1 - $pcpea2); fi
	pcpea4=$(echo $pcpea3 | sed "s/^/$pcpea /g" | awk '{if($2~"-") print$1"("$2"⬇)"; else print$1"("$2"⬆)"}' | tr -d "-")
	echo $pcpea4
	done | tr "()" " " | sort -k2,2nr | sed 's/ /(/g' | sed 's/($/)/g' | paste -s -d ","
	unset pcpea pcpea1 pcpea2 pcpea3 pcpea4
	done < <(echo "$pcp" | datamash transpose -t " ") | paste -s -d " ")
	unset pcpe
echo $p $pih $c1 $c2 $o $pcp $pcpmd
unset pih c1 c2 o pcp pcpmd
rm temp_count_birds temp_count_mammals
done | sed '1i Position_in_alignment Position_in_human Residues_in_Birds(count) Residues_in_mammals(count) Common_residues(count_in_birds,count_in_mammals) Hydropathy(B) Volume(B) Chemical(B) Charge(B) Polarity(B) Hydrogen_donor_or_acceptor_atom(B) Hydropathy(M) Volume(M) Chemical(M) Charge(M) Polarity(M) Hydrogen_donor_or_acceptor_atom(M) Hydropathy(D) Volume(D) Chemical(D) Charge(D) Polarity(D) Hydrogen_donor_or_acceptor_atom(D)' | column -t > compare_residues_43_birds_32_mammals

awk '{print$1,$2,$3,$4,$5,$12,$6,$13,$7,$14,$8,$15,$9,$16,$10,$17,$11,$18,$19,$20,$21,$22,$23}' compare_residues_43_birds_32_mammals | sed 's/[ ]\+/\t/g' > tsv_compare_residues_43_birds_32_mammals.txt

#Calculate the sum of all physicochemical changes in birds w.r.t mammals
for i in $(seq 18 23)
do
i1=$(awk -v a="$i" 'NR==1{print ">"$a}' compare_residues_43_birds_32_mammals)
i2=$(awk -v a="$i" '{print $a}' compare_residues_43_birds_32_mammals | sed 1d | awk NF | tr "," "\n" | sed 's/(/ /g' | sed 's/⬆/ ⬆/g' | sed 's/⬇/ ⬇/g' | tr -d ")" | awk '{key=$1 FS $3; sum[key]+=$2} END{for(k in sum) print k, sum[k]}' | sort -k3,3nr | awk '{print$1"("$3$2")"}' | paste -s -d ",")
i3=$(awk -v a="$i" '{print $a}' compare_residues_43_birds_32_mammals | sed 1d | awk NF | tr "," "\n" | sed 's/(/ /g' | sed 's/⬆/ ⬆/g' | sed 's/⬇/ ⬇/g' | tr -d ")" | awk '{key=$1 FS $3; sum[key]+=$2} END{for(k in sum) print k, sum[k]}' | sort -k1,1 -k2,2 | sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1 /;t1;P;D' | awk '{print$1,$3-$5}' | sort -k2,2nr | awk '{if($2~"-") print$1"("$2"⬇)"; else print$1"("$2"⬆)"}' | paste -s -d ",")
echo $i1 $i2 $i3
unset i1 i2 i3
done | column -t > major_differences_compare_residues_43_birds_32_mammals

#Birds lack ~40 aa C-terminal region & some indels are present in mammals as well. Since such a large insertion in mammals skew the major trends observed in substitutions selected for specific function
awk '$3 != "-(43)" && $3 !~ /-\(32\)/ && $3 !~ /-\(36\)/ && $3 !~ /X\(29\)/ && $4 != "-(32)"' compare_residues_43_birds_32_mammals > compare_residues_43_birds_32_mammals_indels_removed

for i in $(seq 18 23)
do
i1=$(awk -v a="$i" 'NR==1{print$a}' compare_residues_43_birds_32_mammals_indels_removed)
i2=$(awk -v a="$i" '{print $a}' compare_residues_43_birds_32_mammals_indels_removed | sed 1d | awk NF | tr "," "\n" | sed 's/(/ /g' | sed 's/⬆/ ⬆/g' | sed 's/⬇/ ⬇/g' | tr -d ")" | awk '{key=$1 FS $3; sum[key]+=$2} END{for(k in sum) print k, sum[k]}' | sort -k3,3nr | awk '{print$1"("$3$2")"}' | paste -s -d ",")
i3=$(awk -v a="$i" '{print $a}' compare_residues_43_birds_32_mammals_indels_removed | sed 1d | awk NF | tr "," "\n" | sed 's/(/ /g' | sed 's/⬆/ ⬆/g' | sed 's/⬇/ ⬇/g' | tr -d ")" | awk '{key=$1 FS $3; sum[key]+=$2} END{for(k in sum) print k, sum[k]}' | sort -k1,1 -k2,2 | sed ':1;$!N;s/^\(\(\S\+\s\+\).*\)\n\2/\1 /;t1;P;D' | awk '{print$1,$3-$5}' | sort -k2,2nr | awk '{if($2~"-") print$1"("$2"⬇)"; else print$1"("$2"⬆)"}' | paste -s -d ",")
echo $i1 $i2 $i3
unset i1 i2 i3
done | column -t > major_differences_compare_residues_43_birds_32_mammals_indels_removed

#sort the alignment based on phylogeny
#Get time calibrated phylogenetic tree from "time tree" using the list "birds_mammals_filtered_8_names"
sed 's/_/ /g' birds_mammals_filtered_8_names.nwk > birds_mammals_filtered_8_names_underscore_removed.nwk

#Make final alignment for figure with just species names 
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cds_alignment
sed '/^>/ s/^>B_/>/g' birds_mammals_filtered_8.fa | sed '/^>/ s/^>M_/>/g' | sed '/^>/ s/(.*)//g' > birds_mammals_filtered_8_only_species_names.fa
transeq birds_mammals_filtered_8_only_species_names.fa --auto --stdout --clean | myfasta -comb | sed '/^>/ s/_1//g' > birds_mammals_filtered_8_only_species_names.aa
mafft --maxiterate 1000 --localpair birds_mammals_filtered_8_only_species_names.aa > birds_mammals_filtered_8_only_species_names_mafft.aln
muscle -in birds_mammals_filtered_8_only_species_names.aa -maxiters 1000 -out birds_mammals_filtered_8_only_species_names_muscle.aln

#View in jalview & mark the sites based on table & compare species order with phylogenetic tree


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Amnino acid guided DNA aligner Tools tried:

#Need a protein guided codon or dna alignment

#---------------------------------------------------------------------------------------------------------------
#Commandline:
  #translator x : https://anaconda.org/bioconda/translatorx
    conda install bioconda::translatorx
    conda activate base  
  #coati : https://github.com/CartwrightLab/coati
    #Dependencies
    wget https://github.com/CartwrightLab/coati/archive/refs/tags/v1.0.0.zip
    unzip v1.0.0.zip
    sudo apt install meson
    meson setup builddir --buildtype=release
    meson compile -C builddir
    tar -xvzf coati*.tar.gz
    meson setup builddir --buildtype=release
    meson compile -C builddir
    meson install -C builddir
#Framealign: https://bip.weizmann.ac.il/education/materials/gcg/framealign.html#command-line
#MACSE: https://www.agap-ge2pop.org/macsee-pipelines/

#---------------------------------------------------------------------------------------------------------------
#GUI
  #PFAAT:
    #Download PFAAt aligner in mac
  #Codoncode:
    #Download codoncode aligner in mac (free trial demo)
    #Find mutation option is still not available in demo mode
#---------------------------------------------------------------------------------------------------------------
#Webserver
  #Codon alignment (HCV Sequence database): https://hcv.lanl.gov/content/sequence/CodonAlign/codonalign.html






##################################################################################################################################################################################################################################################################################################################################
#DRAFT
##################################################################################################################################################################################################################################################################################################################################

#Make physiochemical properties table
for p in $(seq 1 $(awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' jalview_alignment | sed -n 2p | wc | awk '{print$3-$1}'))
do
pih=$(awk -v RS=">" '/Homo_sapiens/ {print">"$0}' jalview_alignment | awk NF | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' | awk -v RS=">" '{$1=$1}1' | awk NF | awk -v a="$p" '{print $1,substr($2,0,a)}' | awk '{print$2}' | tr -d "-" | wc | awk '{print$3-$1}')
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' jalview_alignment | awk -v RS=">" '{$1=$1}1' | awk -v a="$p" '{print $1,substr($2,a,1)}' | awk '/^B_/{print$2}' | sort | uniq -c | sort -k1,1nr > temp_count_birds
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' jalview_alignment | awk -v RS=">" '{$1=$1}1' | awk -v a="$p" '{print $1,substr($2,a,1)}' | awk '/^M_/{print$2}' | sort | uniq -c | sort -k1,1nr > temp_count_mammals
c1=$(cat temp_count_birds | awk '{print$2"("$1")"}' | paste -s -d ",")
c2=$(cat temp_count_mammals | awk '{print$2"("$1")"}' | paste -s -d ",")
o=$(sed 's/^[ ]\+//g' temp_count_birds | xargs -n2 bash -c 'paste <(grep $1 temp_count_mammals | awk "{print\$2,\$1}") <(echo $0) -d " "' | awk '$3' | awk '{print$1"("$3","$2")"}' | paste -s -d ",")
echo $p $pih $c1 $c2 $o
unset pih c1 c2 o
rm temp_count_birds temp_count_mammals
done | sed '1i Position_in_alignment Position_in_human Residues_in_Birds(count) Residues_in_mammals(count) Common_residues(count_in_birds,count_in_mammals)' | column -t > compare_residues_43_birds_32_mammals

#
./phychem_compare.sh birds_mammals_filtered_8_with_group_symbol_mafft.aln

