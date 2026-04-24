
############################################################################################################################################################################################################################################################################################################
#TEST ADH gene

#Get example codes & input
cd /media/aswin/gene_loss/APOBEC1/selection/positive_selection/paml
git clone https://github.com/Dating-gene-loss/Mammal_ADH_IV.git

#Prepare folders
cd /media/aswin/gene_loss/APOBEC1/selection/positive_selection/paml/Mammal_ADH_IV
cd /media/aswin/gene_loss/APOBEC1/selection/positive_selection/paml/Mammal_ADH_IV/Dating/codeml/F1X4_model/hyphy

#Identify labels used
cd ~/bird_db1/aswin/APOBEC1/Dating/paml/Mammal_ADH_IV/Dating/codeml/F3X4_model
cat ADH_PAML_tree_unrooted.nwk | tr -d ")(;" | tr "," "\n" | grep -v "#" > bg.txt
cat ADH_PAML_tree_unrooted.nwk | tr -d ")(;" | tr "," "\n" | grep "#" | awk '$2=="#1"' | cut -f1 -d " " > fg1.txt
cat ADH_PAML_tree_unrooted.nwk | tr -d ")(;" | tr "," "\n" | grep "#" | awk '$2=="#2"' | cut -f1 -d " " > fg2.txt

#Get all pairwise split time
cd ~/bird_db1/aswin/APOBEC1/Dating/paml/Mammal_ADH_IV/Dating/codeml/F3X4_model
cp ~/bird_db1/aswin/APOBEC1/Dating/paml/COA1_GENE/Geneloss_timing/Galliformes/AllPseudogene_AllMix_AllFunctional/test .
cp ../../../PGLS/Mammal_ADH_tree_revisions.nwk .
Rscript ./splittime.r

#Manually create the groups 


#Calculate gene loss timings
time for sp in $(cat fg1.txt)
do
./calculate_gene_inactivation.sh $sp -f F3X4_model.out ../F1X4_model/F1X4_model.out ADH_PAML_tree_unrooted.nwk Mammal_ADH_tree_revisions.nwk -s | grep -v Mixed_branch_length
done | sed '1i Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' > fg1_inactivation.tsv

time for sp in $(cat fg2.txt)
do
gr=$(grep $sp groups | awk '{print$2}')
./calculate_gene_inactivation.sh $sp -f F3X4_model.out ../F1X4_model/F1X4_model.out ADH_PAML_tree_unrooted.nwk Mammal_ADH_tree_revisions.nwk -s | grep -v Mixed_branch_length | sed "s/^/$gr\t/g"
unset gr
done | sed '1i Group Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' > fg2_inactivation.tsv

#With newer code
~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh Equus_caballus -f ../F1X4_model/F1X4_model.out F3X4_model.out -wp=d ADH_PAML_tree_unrooted.nwk Mammal_ADH_tree_revisions.nwk 


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Using hyphy

#Remove codeml specific labels from tree
sed 's/#1//g' ../ADH_PAML_tree_unrooted.nwk | sed 's/#2//g' > ADH_PAML_tree_unrooted_unlabelled.nwk

#ABSREL (13m11.040s)
time /media/aswin/programs/hyphy-2.5.70/hyphy absrel --alignment ../ADH_PAML_alignment.phy --tree ADH_PAML_tree_unrooted_unlabelled.nwk --output absrel.json &> absrel.stdout

#FitMG94 (21m14.930s)
time /media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/FitMG94/FitMG94.bf --tree ADH_PAML_tree_unrooted_unlabelled.nwk --alignment ADH_PAML_alignment.phy --output FitMG94.json &> FitMG94.stdout
time /media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/FitMG94/FitMG94.bf --tree ADH_PAML_tree_unrooted_unlabelled.nwk --alignment ADH_PAML_alignment.phy --output FitMG94_lineage.json --type lineage &> FitMG94_lineage.stdout

#Relabel hyphy specific 
sed -E 's/#1/{FG}/g; s/#2//g' ../ADH_PAML_tree_unrooted.nwk > ADH_PAML_tree_unrooted_labelled.nwk
#FitMG94 20m5.393s
time /media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/FitMG94/FitMG94.bf --tree ADH_PAML_tree_unrooted_labelled.nwk --alignment ADH_PAML_alignment.phy --output FitMG94_lineage_labelled_tree.json --type lineage &> FitMG94_lineage_labelled_tree.stdout
time /media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/FitMG94/FitMG94.bf --tree ADH_PAML_tree_unrooted_labelled.nwk --alignment ADH_PAML_alignment.phy --output FitMG94_local_labelled_tree.json --type local &> FitMG94_local_labelled_tree.stdout

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Using codeml

#codeml 37m31.120s
cd /media/aswin/gene_loss/APOBEC1/selection/positive_selection/paml/Mammal_ADH_IV/Dating/codeml/F1X4_model/paml

############################################################################################################################################################################################################################################################################################################
#Test COA1 gene

#Prepare folders
cd /media/aswin/gene_loss/APOBEC1/selection/positive_selection/paml/COA1_GENE/Geneloss_timing/Galliformes/All_mixed_with_functional
mkdir test
cp Galliformes_F1x4.ctl Galliformes_F3x4.ctl Galliformes.aln galliformes_labelled.nwk test/
cd /media/aswin/gene_loss/APOBEC1/selection/positive_selection/paml/COA1_GENE/Geneloss_timing/Galliformes/All_mixed_with_functional/test

#Using codeml (0m28.427s)
time /media/aswin/programs/paml-4.10.10-linux-x86_64/bin/codeml Galliformes_F3x4.ctl 


############################################################################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################################################################
#APOBEC1

#Prepare folders
mkdir ~/bird_db1/aswin/APOBEC1/Dating
cd ~/bird_db1/aswin/APOBEC1/Dating

#Get input sequences
mkdir ~/bird_db1/aswin/APOBEC1/Dating/cds
cp ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*.fa ~/bird_db1/aswin/APOBEC1/Dating/cds
cd ~/bird_db1/aswin/APOBEC1/Dating/cds
ls APOBEC1_*.fa | cut -f2- -d "_" | sed 's/\.fa//g' > intact_apobec1

############################################################################################################################################################################################################################################################################################################
#QC fasta sequence

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC of each exons

	#Split exons & combine into separate sequence
	mkdir /home/neo/bird_db1/aswin/APOBEC1/Dating/cds/exon_wise
	while read -r sp; do
	echo "$sp"
	f=$(ls | grep -m1 "$sp")
	for exon in {1..5}; do
	myfasta -mfp "$f" "exon_${exon}" | sed "s/^>.*/>$sp/" >> "exon_wise/exon_${exon}.fa"
	done
	done < intact_apobec1

	#Align each exons
	for i in $(ls exon*.fa)
	do
	echo ">"$i
	mafft --maxiterate 1000 --localpair --reorder --quiet $i > $i".aln"
	done

	cp /home/neo/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.fa ~/bird_db1/aswin/APOBEC1/Dating
	myfasta -l /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Acanthisitta_chloris.fa

	#Manually create 3 lists of species groups that lost APOBEC1 
	grep galliformes -i ~/bird_db1/aswin/APOBEC1/pipeline_summary_row_wise | awk '{print$1}' | tr -d ">" | grep -v "Alectura_lathami" > lost_galliformes
	myfasta -mfl toga_subjects.fa lost_palaeognathae | myfasta -l | sort -k2nr,2
	myfasta -mfl toga_subjects.fa lost_independent | myfasta -l | sort -k2nr,2

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC exon & gap placement 

	#Manually check accuracy of consensus sequence of Struthio_camelus_australis before using for gene loss dating
	cd ~/soft_links/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/manual_search

	#Check how inframe exon_4 should start and end
	for i in $(find ~/bird_db1/aswin/APOBEC1/validated_sequences/ -name "APOBEC1_*.fa")
	do
	echo $i | awk -F "/" '{print"@"$NF}'
	myfasta -codex $i -w -col $i -w -col
	done > codon_splitted_all_validated_APOBEC1_birds

	#Compare if ostrich exon_4 is similar to the phase zero (inframe and codon start from 1st position) exon_4 from other validated birds.
	egrep -a "^@|>exon_4" -A1 codon_splitted_all_validated_APOBEC1_birds | egrep -v "exon|--" | awk -v RS="@" '{$1=$1}1' > exon_4_codon_splitted_all_validated_APOBEC1_birds
	#Align only exon 4 from ostrich with phase zero exon_4 of validated birds
	cat phase_zero_exon_4_APOBEC1_Struthio_camelus.fa <(cat ~/bird_db1/aswin/APOBEC1/Dating/cds/exon_wise/exon_4.fa | sed 's/^AT//g' | sed 's/^at//g') > exon_4_all_validated_with_Struthio_camelus.fa
	mafft --maxiterate 1000 --localpair --quiet --reorder exon_4_all_validated_with_Struthio_camelus.fa > exon_4_all_validated_with_Struthio_camelus.aln

############################################################################################################################################################################################################################################################################################################
#Collect sequence of lost species

mkdir ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species
mkdir ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA

#Collect TOGA given A1 sequnece of lost birds  
cd ~/bird_db1/aswin/APOBEC1/TOGA
for sp in $(cat all_toga_subjects_in_phylogenetic_order)
do
cp "$sp"/toga_rna-XM_027449866.2_exon_wise_consensus.fa ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA/"$sp"_toga.fa
done
#For one species the naming is slightly different
cp ~/bird_db1/aswin/APOBEC1/TOGA/Chlamydotis_macqueenii/toga_APOBEC1_mrna_exon_wise_consensus.fa ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA/Chlamydotis_macqueenii_toga.fa

############################################################################################################################################################################################################################################################################################################
#QC

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC of each exons & find oultiers

#Create list of list birds with some sequence remains
cd ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species
grep galliformes -i ~/bird_db1/aswin/APOBEC1/pipeline_summary_row_wise | awk '{print$1}' | tr -d ">" | egrep -v "Alectura_lathami|Colinus_virginianus" > lost_galliformes
cat lost_galliformes lost_palaeognathae lost_independent > all_lost_birds_with_remains

#create a list of expected exon lengths for A1
nano expected_exon_lengths

#Extract exon-wise sequence from lost species
cd ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA
awk '/^>/ {exon = substr($0,2)   # remove ">"}
exon ~ /^exon_[1-5]$/ {
    outfile = exon ".fa"
    getline seq
    print ">" FILENAME >> outfile
    print seq >> outfile}' *.fa

#Manually check for gappy alignments in each exons
	cat exon_1.fa ~/bird_db1/aswin/APOBEC1/Dating/cds/exon_wise/exon_1.fa | myfasta -vfp - Buceros_rhinoceros > all_exon_1.fa
	mafft --maxiterate 1000 --localpair --quiet --reorder all_exon_1.fa > all_exon_1.aln
cat exon_2.fa ~/bird_db1/aswin/APOBEC1/Dating/cds/exon_wise/exon_2.fa > all_exon_2.fa
mafft --maxiterate 1000 --localpair --quiet --reorder all_exon_2.fa > all_exon_2.aln
	cat exon_3.fa ~/bird_db1/aswin/APOBEC1/Dating/cds/exon_wise/exon_3.fa > all_exon_3.fa
	mafft --maxiterate 1000 --localpair --quiet --reorder all_exon_3.fa > all_exon_3.aln
cat exon_4.fa ~/bird_db1/aswin/APOBEC1/Dating/cds/exon_wise/exon_4.fa > all_exon_4.fa
mafft --maxiterate 1000 --localpair --quiet --reorder all_exon_4.fa > all_exon_4.aln
	cat exon_5.fa ~/bird_db1/aswin/APOBEC1/Dating/cds/exon_wise/exon_5.fa > all_exon_5.fa
	mafft --maxiterate 1000 --localpair --quiet --reorder all_exon_5.fa > all_exon_5.aln

#NOTE: 
	#exon_1: Buceros_rhinoceros had very long exon_1 hence remove it from analysis altogether
	#exon_2: Meleagris_gallopavo_toga.fa has shorter 5' end, Pelecanus_crispus_toga.fa & Pelecanus_crispus had shorter 3' end
	#exon_3: Galliformes had a long stretch lacking compared to others in exon_3 likely due to repeat insertion
	#exon_4: Cathartes_aura_toga.fa & Gymnogyps_californianus_toga.fa had longer 3' end, Tyto_alba_toga.fa has shorter 3' end.
	#exon_5: Limosa_lapponica has longer 3' end

#Create a file containing a list of species_to_remove
cd ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA
for i in $(grep -vf species_to_remove <(ls *_toga.fa))
do
n=$(echo $i | sed 's/_toga.fa//g')
echo ">"$n
grep -v ">" $i | paste -s -d ""
unset n
done

grep -f all_lost_birds_with_remains ~/bird_db1/aswin/APOBEC1/pipeline_summary_column_wise | colnum.sh

#Compare expected and observed lengths of each exons from species with loss
for sp in $(grep -f all_lost_birds_with_remains <(ls TOGA/*))
do
n=$(echo $sp | cut -f2 -d "/" | sed 's/_toga.fa//g')
gr=$(grep $n lost_galliformes lost_palaeognathae lost_independent | cut -f1 -d ":" | cut -f2 -d "_")
ol=$(myfasta -l $sp | awk '{split($1,a,"_"); len[a[2]]=$2}
END {for (i=1;i<=5;i++) printf "%s ", (i in len ? len[i] : "-"); print ""}')
ofl=$(myfasta -l $sp | awk '{a+=$2} END{print a}')
el=$(grep -f <(grep ">" $sp | cut -f2 -d "_") expected_exon_lengths | awk '{a+=$2} END{print a}')
dif=$(calc $el -$ofl)
echo $n $gr $ol $ofl $el $dif 
unset n gr ol ofl el dif
done | sort -k10nr,10 | sed '1i Species Group exon_1 exon_2 exon_3 exon_4 exon_5 Total_observed_length Total_expected_length Difference' | column -t > length_exon_wise_comparison_observed_with_expected

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Filter

#Extract cds of all lost birds in single multifasta
for sp in $(grep -f <(cat ../lost_galliformes ../lost_independent ../lost_palaeognathae) <(ls *_toga.fa))
do
echo ">"$sp | sed 's/_toga.fa//g'
grep -v ">" $sp | paste -s -d ""
done > all_lost_species_cds.fa

#Filter species with large gaps
myfasta -vpl all_lost_species_cds.fa species_to_remove > all_lost_species_cds_filtered.fa 

#Extract cds of all lost birds in single multifasta
cd ~/bird_db1/aswin/APOBEC1/Dating/cds
for sp in $(grep -f intact_apobec1 <(ls APOBEC1*.fa))
do
echo $sp | cut -f2- -d "_" | sed 's/.fa//g' | sed 's/^/>/g'
grep -v ">" $sp | paste -s -d ""
done > all_intact_species_cds.fa

#Filter species with large gaps
myfasta -vpl all_intact_species_cds.fa lost_species/TOGA/species_to_remove > all_intact_species_cds_filtered.fa 

#Final fasta of filtered birds 
cat all_intact_species_cds_filtered.fa lost_species/TOGA/all_lost_species_cds_filtered.fa > all_birds_filtered.fa

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC of LINE repeats in galliformes

#Only exon_3 in galliformes seems to have repeats insertions overlapping exon_3
#To use final cds for gene loss dating remove this LINE repeat insertion
#See how many species has repeats included in final consensus given by TOGA; exon_3 alignment between only galliformes & validated species
cd ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA
for sp in $(grep -f ../lost_galliformes <(ls *toga.fa)); do echo ">"$sp; myfasta -mpp "$sp" "exon_3" | grep -v ">"; done | myfasta -de > galliformes_exon_3.fa
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/cds/exon_wise/exon_3.fa intact_exon_3.fa
cat galliformes_exon_3.fa intact_exon_3.fa > galliformes_and_intact_exon_3.fa
mafft --maxiterate 1000 --localpair --reorder --quiet galliformes_exon_3.fa > galliformes_exon_3.aln
mafft --maxiterate 1000 --localpair --reorder --quiet galliformes_and_intact_exon_3.fa > galliformes_and_intact_exon_3.aln

#Check exact overlap region between repeat & exon_3
mkdir ~/bird_db1/aswin/APOBEC1/Dating/repeat_removal
cd ~/bird_db1/aswin/APOBEC1/Repeat_find2
for sp in $(cat ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/lost_galliformes)
do
bedtools getfasta -fi "$sp"/"$sp"_APOBEC1_region.fa -bed <(bedtools intersect -a "$sp"/adjusted_TOGA.bed -b "$sp"/repeat.bed) -s -name+ | sed "s/>/>${sp}_/g"
done > ~/bird_db1/aswin/APOBEC1/Dating/repeat_removal/repeat_overlaps_in_galliformes_exon_3.fa
cd ~/bird_db1/aswin/APOBEC1/Dating/repeat_removal
myfasta -de repeat_overlaps_in_galliformes_exon_3.fa > repeat_overlaps_in_galliformes_exon_3_non_empty.fa
mafft --maxiterate 1000 --localpair --reorder --quiet repeat_overlaps_in_galliformes_exon_3_non_empty.fa > repeat_overlaps_in_galliformes_exon_3_non_empty.aln

#Inspect the overlap region in each species
for sp in $(grep ">" ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA/galliformes_exon_3.fa | sed 's/_toga.fa//g' | tr -d ">")
do
echo ">"$sp
cd ~/bird_db1/aswin/APOBEC1/Repeat_find2/"$sp"
repeat=$(bedtools intersect -a repeat.bed -b adjusted_TOGA.bed -wo | awk '{print$4}')
bedtools getfasta -fi "$sp"_APOBEC1_region.fa -bed repeat.bed -s -name+ | myfasta -mpl - <(echo "$repeat") > ~/bird_db1/aswin/APOBEC1/Dating/repeat_removal/"$sp"_repeat_exon_3.fa
myfasta -mfp ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA/"$sp"_toga.fa exon_3 | sed 's/^>.*/>toga/g' >> ~/bird_db1/aswin/APOBEC1/Dating/repeat_removal/"$sp"_repeat_exon_3.fa
myfasta -mfp ~/soft_links/"$sp"/aswin/APOBEC1/2nd_gblast/extracted_flanking_region.fa exon_3 | sed 's/^>.*/>flanking_region/g' >> ~/bird_db1/aswin/APOBEC1/Dating/repeat_removal/"$sp"_repeat_exon_3.fa
mafft --maxiterate 1000 --localpair --reorder --quiet ~/bird_db1/aswin/APOBEC1/Dating/repeat_removal/"$sp"_repeat_exon_3.fa > ~/bird_db1/aswin/APOBEC1/Dating/repeat_removal/"$sp"_repeat_exon_3.aln
#Consensus
#em_cons ~/bird_db1/aswin/APOBEC1/Dating/repeat_removal/"$sp"_repeat_exon_3.aln -plurality 3 --auto --stdout | myfasta -comb | awk '/^>/ {print; next} {gsub(/[Nn]/,""); print}'
bedtools getfasta -fi "$sp"_APOBEC1_region.fa -bed <(bedtools intersect -a adjusted_TOGA.bed -b repeat.bed) -s -name+ > ~/bird_db1/aswin/APOBEC1/Dating/repeat_removal/"$sp"_repeat_exon_3_overlap.fa
unset repeat
done

cd ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA
cp galliformes_* intact_exon_3.fa ~/bird_db1/aswin/APOBEC1/Dating/repeat_removal/

#Realized coturnix have some extra repeats at starting of exon_3 that is absent in other birds
sed 's/ccaggaattGGACTTGGTGATCAATAactgtcccttccaactcag//g' galliformes_exon_3.fa > galliformes_exon_3_cotutnix_trimmed.fa
mafft --maxiterate 1000 --localpair --reorder --quiet galliformes_exon_3_cotutnix_trimmed.fa > galliformes_exon_3_cotutnix_trimmed.aln
grep -f <(grep -v ">" repeat_overlaps_in_galliformes_exon_3.fa | tr "[A-Z]" "[a-z]") <(myfasta -comb galliformes_exon_3_cotutnix_trimmed.aln) -z 

#Looks like almost all LINE repeats ends with or have "GTGAATTCT" "TGAATTCT" near end. hence use this to trim exon_3 in all galliformes, such that all galiformes will have a normalized sequence start and manually made sure this seq doesn't randomly anywhere else
sed 's/.*TGAATTCT//Ig' galliformes_exon_3.fa > galliformes_exon_3_trimmed.fa
mafft --maxiterate 1000 --localpair --reorder --quiet galliformes_exon_3_trimmed.fa > galliformes_exon_3_trimmed.aln
#This alignment looked reasonable for final codon alignment
alv galliformes_exon_3_trimmed.aln


############################################################################################################################################################################################################################################################################################################
#Codon alignment

#Extract cds of all lost birds in single multifasta
cd ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA
for sp in $(grep -f <(cat ../lost_independent ../lost_palaeognathae) <(ls *_toga.fa))
do
echo ">"$sp | sed 's/_toga.fa//g'
grep -v ">" $sp | paste -s -d ""
done > palaeognathae_independent_lost_species_cds.fa

#Trim repeat region of exon_3 in only galliformes
for sp in $(grep -f <(cat ../lost_galliformes) <(ls *_toga.fa))
do
echo ">"$sp | sed 's/_toga.fa//g'
awk '/^>/ {in_exon3 = ($0 ~ /^>exon_3/); print; next} in_exon3 {IGNORECASE=1; sub(/.*TGAATTCT/,"")} {print}' "$sp" | grep -v ">" | paste -s -d ""
#grep -v ">" $sp | paste -s -d ""
done > galliformes_cds_trimmed.fa

cat palaeognathae_independent_lost_species_cds.fa galliformes_cds_trimmed.fa > all_lost_species_cds_trimmed.fa
#Filter species with large gaps
myfasta -vpl all_lost_species_cds_trimmed.fa species_to_remove > all_lost_species_cds_trimmed_filtered.fa 

#Combine
cd ~/bird_db1/aswin/APOBEC1/Dating/cds
cat all_intact_species_cds_filtered.fa lost_species/TOGA/all_lost_species_cds_trimmed_filtered.fa > all_birds_trimmed_filtered.fa

mkdir ~/bird_db1/aswin/APOBEC1/Dating/alignment
cd ~/bird_db1/aswin/APOBEC1/Dating/alignment
cp ~/bird_db1/aswin/APOBEC1/Dating/cds/all_birds_trimmed_filtered.fa .

#Use protein alignment & translate it to codon level
#Get protein sequence
transeq all_birds_trimmed_filtered.fa --auto --stdout --clean | myfasta -comb | sed '/^>/ s/_1//g' > all_birds_trimmed_filtered.aa
#protein alignment
mafft --maxiterate 1000 --localpair --reorder --quiet all_birds_trimmed_filtered.aa > all_birds_trimmed_filtered.aa.mafft.aln
muscle -in all_birds_trimmed_filtered.aa > all_birds_trimmed_filtered.aa.muscle.aln
#mafft looks better, since species such as Cuculus_canorus have only exon_4 & 5, & mafft shows long gap at start & align at end, but muscle shows interleaved alignment from middle of alignment.
time mafft --maxiterate 1000 --localpair --reorder --quiet all_birds_filtered.fa > all_birds_filtered.aln

#NOTE: Observed pal2nal alignment failed since the trimmed fasta is not divisible by 3
#This is expected as gene loss essentially means ORF disruption or PTC.
cd ~/bird_db1/aswin/APOBEC1/Dating/alignment
myfasta -l all_birds_trimmed_filtered.fa | awk '{print$0,$2/3}'

############################################################################################################################################################################################################################################################################################################
#QC to correct frameshitfts in input sequence

	#Get all filtered intact A1 seq
	mkdir -p ~/bird_db1/aswin/APOBEC1/Dating/alignment/protein_ref/exon_wise_protein/intact && cp ~/bird_db1/aswin/APOBEC1/Dating/cds/APOBEC1_*.fa ~/bird_db1/aswin/APOBEC1/Dating/alignment/protein_ref/exon_wise_protein/intact
	grep -f ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA/species_to_remove <(ls ~/bird_db1/aswin/APOBEC1/Dating/alignment/protein_ref/exon_wise_protein/intact/APOBEC1_*.fa) | xargs rm
	#Get all filtered lost A1 seq
	mkdir ~/bird_db1/aswin/APOBEC1/Dating/alignment/protein_ref/exon_wise_protein/lost && cd ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA
	grep -f <(grep -vf species_to_remove <(cat ../lost_galliformes ../lost_palaeognathae ../lost_independent)) <(ls *_toga.fa) | xargs -n1 sh -c 'cp $0 ~/bird_db1/aswin/APOBEC1/Dating/alignment/protein_ref/exon_wise_protein/lost/'
	#Make protein of each intact seq
	cd ~/bird_db1/aswin/APOBEC1/Dating/alignment/protein_ref/exon_wise_protein/intact
	for i in $(ls APOBEC1*.fa); do n=$(echo $i | sed 's/.fa/_protein.fa/g'); myfasta -cds $i | transeq stdin --auto --stdout -clean | myfasta -comb > $n; unset n; done
	#Make cds of each intact seq
	cd ~/bird_db1/aswin/APOBEC1/Dating/alignment/protein_ref/exon_wise_protein/lost
	for i in $(ls *_toga.fa); do n=$(echo $i | sed 's/.fa/_cds.fa/g'); myfasta -cds $i > $n; unset n; done

	#Make consensus of all intact
	cd ~/bird_db1/aswin/APOBEC1/Dating/alignment/protein_ref/exon_wise_protein/intact
	transeq ~/bird_db1/aswin/APOBEC1/Dating/cds/all_intact_species_cds_filtered.fa --auto --stdout -clean | myfasta -comb | sed 's/_1//g' > all_intact_species_cds_filtered.aa
	mafft --auto --reorder --quiet all_intact_species_cds_filtered.aa > all_intact_species_cds_filtered.aln

#Clustering
	mkdir ~/bird_db1/aswin/APOBEC1/Dating/alignment/clustering
	cd ~/bird_db1/aswin/APOBEC1/Dating/alignment/clustering
	transeq ~/bird_db1/aswin/APOBEC1/Dating/cds/all_intact_species_cds_filtered.fa --auto --stdout -clean | myfasta -comb | sed 's/_1//g' > all_intact_species_cds_filtered.aa
	#Run clans
	java -Xmx4G -jar ~/programmes/clans.jar APOBEC1_intact_birds.clans

#Refine intact A1 group if possible
	mkdir ~/bird_db1/aswin/APOBEC1/Dating/alignment/filtering
	cd ~/bird_db1/aswin/APOBEC1/Dating/alignment/filtering
	mafft --auto --reorder --quiet ~/bird_db1/aswin/APOBEC1/Dating/cds/all_intact_species_cds_filtered.fa > all_intact_species_cds_filtered.aln
	mafft --auto --reorder --quiet <(myfasta -vpp ~/bird_db1/aswin/APOBEC1/Dating/cds/all_intact_species_cds_filtered.fa Malurus_cyaneus_samueli) > all_intact_species_cds_filtered_malurus_removed.aln

#Check protein-nucleotide alignment
	mkdir ~/bird_db1/aswin/APOBEC1/Dating/alignment/protein_ref/exon_wise_protein/pairwise_protein_align
	cd ~/bird_db1/aswin/APOBEC1/Dating/alignment/protein_ref/exon_wise_protein/pairwise_protein_align
	time for i in $(find ../lost/ -name "*_toga_cds.fa")
	do
	n=$(echo $i | awk -F "/" '{print$NF}' | sed 's/_toga_cds.fa//g')
	#Use 3 different reference based on total length & divergence time
	for ref in ../intact/APOBEC1_Acanthisitta_chloris_protein.fa ../intact/APOBEC1_Columba_livia_protein.fa ../intact/APOBEC1_Ficedula_albicollis_protein.fa
	do
	refn=$(echo $ref | awk -F "/" '{print$NF}' | sed 's/_protein.fa//g')
	echo ">Test:" $n "| Reference: "$refn | GREP_COLORS="mt=04;01;37" grep ">Test.*" --color=always
	#exonerate -E --model p2g:b $ref $i --ryo "Rank: %r\nTarget in alignment: %ti : %tal (%tab - %tae)\nTarget in coding sequence: %ti : %tcl (%tcb - %tce)\nGene orientation: %g\nPercent identity: %pi\nPercent similarity: %ps\n>%ti (%tab - %tae)\n%tcs\n" --querytype protein --targettype dna --showtargetgff yes --fsmmemory 28000 --alignmentwidth 270 > "$n"_"$refn"_codon_compare.out
	#spaln -M1 -l254 -O1 -ya0 -yX0 -yS50 $i $ref >> "$n"_"$refn"_codon_compare.out
	exonerate -E --model p2g:b $ref $i --ryo "Rank: %r\nTarget in alignment: %ti : %tal (%tab - %tae)\nTarget in coding sequence: %ti : %tcl (%tcb - %tce)\nGene orientation: %g\nPercent identity: %pi\nPercent similarity: %ps\n>%ti (%tab - %tae)\n%tcs\n" \
		--querytype protein --targettype dna --showtargetgff yes --fsmmemory 28000 --alignmentwidth 10000 | awk '/Target range:/{flag=1; next} /vulgar:/{flag=0} flag' | awk NF > kala_temp1
	if [ -z "$(grep "C4 Alignment" kala_temp1)" ]
	then 
	exonerate --model p2g $ref $i --ryo "Rank: %r\nTarget in alignment: %ti : %tal (%tab - %tae)\nTarget in coding sequence: %ti : %tcl (%tcb - %tce)\nGene orientation: %g\nPercent identity: %pi\nPercent similarity: %ps\n>%ti (%tab - %tae)\n%tcs\n" \
		--querytype protein --targettype dna --showtargetgff yes --fsmmemory 28000 --alignmentwidth 10000 | awk '/Target range:/{flag=1; next} /vulgar:/{flag=0} flag' | awk NF > kala_temp1
	else :
	fi
	cat kala_temp1
	spaln -M1 -l10000 -O1 -ya0 -yX0 -yS50 $i $ref | awk '/^ALIGNMENT/{flag=1; next} EOF{flag=0} flag' | sed 's/^  //g' | awk NF > kala_temp2
	if [ -z "$(grep "C4 Alignment" kala_temp2)" ]
	then
	spaln -M1 -l10000 -O1 -ya0 -yX1 -yS50 $i $ref | awk '/^ALIGNMENT/{flag=1; next} EOF{flag=0} flag' | sed 's/^  //g' | awk NF > kala_temp2
	else :
	fi
	cat kala_temp2
	echo -e
	#cat <(awk '/Target range:/{flag=1; next} /vulgar:/{flag=0} flag' "$n"_"$refn"_codon_compare.out) <(awk '/^ALIGNMENT/{flag=1; next} EOF{flag=0} flag' "$n"_"$refn"_codon_compare.out | sed 's/^  //g') > "$n"_"$refn"_codon_compare_aln.out
	unset refn
	rm kala_temp1 kala_temp2
	done | GREP_COLORS="mt=31;07" egrep "#|\*|\-" --color=always -z > $n"_pairwise_protein_nuc_align.out"
	unset ref n
	done

############################################################################################################################################################################################################################################################################################################
#Autmated filtering & Codon alignmnent 

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#HYPHY

cd ~/bird_db1/aswin/APOBEC1/Dating/alignment/hyphy
#Run these sequences through pre-msa.bf in order to correct frame-shift mutations and translate the resulting sequences to proteins.
hyphy ~/programmes/hyphy-analyses/codon-msa/pre-msa.bf --input all_birds_trimmed_filtered.fa
#Take the output of step 2 and run in through the general MSA program to generate a protein MSA
muscle -in all_birds_trimmed_filtered.fa_protein.fas -out all_birds_trimmed_filtered.fa_protein.msa
#Run the protein MSA and the frameshift corrected nucleotide sequences from step 2 through post-msa.bf to obtain a nucleotide msa. This step will also, optionally, compress all identical sequences, i.e. replace them with a single representative sequence.
hyphy ~/programmes/hyphy-analyses/codon-msa/post-msa.bf --protein-msa all_birds_trimmed_filtered.fa_protein.msa --nucleotide-sequences all_birds_trimmed_filtered.fa_nuc.fas --output all_birds_trimmed_filtered.fa_nuc.msa --compress No
#sed '/^>/ s/_[0-9]\+$//g' all_birds_trimmed_filtered.fa_nuc.msa -i

#Observations when inspecting alignment:
	#Pterocles_gutturalis:
			#only lack exon_2 rest of the exons must be aligning well, but in alignment multiple large gaps are seen in exon_3 & a gap was present in the place of whole exon_5
			#exon_3 starting has good homoogy but poor homology at 3' end.
	#Gallus_gallus:
			#exon_3 start has lower homology

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#MACSE pipeline

#Use MACSE Pipeline which automatically filters & corrects the frameshifts (Install MACSE pipeline not simple macse aligner from https://github.com/ranwez/MACSE_V2_PIPELINES.git)
#Use refined & updated pipeline called OMM_MACSE: this pipeline also produces a codon-aware alignment thanks to MACSE, which could be filtered by HmmCleaner, but it can handle larger datasets by relying on MAFFT, MUSCLE or PRANK to scale up.
#Install apptainer, since the pipeline is run using apptainer:
mkdir ~/bird_db1/aswin/APOBEC1/Dating/alignment/macse
cd ~/bird_db1/aswin/APOBEC1/Dating/alignment/macse

#Use different aligners: (Aligner had not much effect in alignment)
	time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa --out_dir mafft --out_file_prefix apobec1 --alignAA_soft MAFFT --aligner_extra_option "--localpair --maxiterate 1000 --reorder" --java_mem 4g 
	time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa --out_dir muscle --out_file_prefix apobec1 --alignAA_soft MUSCLE --out_detail_dir muscle/SAVE_DETAILS/
	time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa --out_dir prank --out_file_prefix apobec1 --alignAA_soft PRANK --out_detail_dir prank/SAVE_DETAILS/

#NOTE: the placement of gaps in alignment is compared with the exon presence absence figure 1; i.e. if Penelope pileata has only exon_5 missing & rest of the exons are present then the only gap observed should be in last 3' end of alignment, not at start or centre.
#At start when used with default parameters, the gaps are observed at centre for Penelope regardless of aligners used, hence try different filtering options.
#Check different filtering strategies
	cd ~/bird_db1/aswin/APOBEC1/Dating/alignment/macse
	time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa --out_dir mafft_no_prefiltering --out_file_prefix apobec1 --alignAA_soft MAFFT --aligner_extra_option "--localpair --maxiterate 1000 --reorder" --java_mem 4g --no_prefiltering
	time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa --out_dir mafft_no_filtering --out_file_prefix apobec1 --alignAA_soft MAFFT --aligner_extra_option "--localpair --maxiterate 1000 --reorder" --java_mem 4g --no_filtering
	time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa --out_dir mafft_no_postfiltering --out_file_prefix apobec1 --alignAA_soft MAFFT --aligner_extra_option "--localpair --maxiterate 1000 --reorder" --java_mem 4g --no_postfiltering
	time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa --out_dir mafft_replace_FS_by_gaps --out_file_prefix apobec1 --alignAA_soft MAFFT --aligner_extra_option "--localpair --maxiterate 1000 --reorder" --java_mem 4g --replace_FS_by_gaps

#Best option
time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa --out_dir mafft_no_prefiltering_no_filtering_no_postfiltering --out_file_prefix apobec1 \
 --alignAA_soft MAFFT --aligner_extra_option "--localpair --maxiterate 1000 --reorder" --java_mem 4g --no_prefiltering --no_filtering --no_postfiltering

time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa --out_dir mafft_no_prefiltering_no_filtering_no_postfiltering_replace_FS_by_gaps --out_file_prefix apobec1 \
 --alignAA_soft MAFFT --aligner_extra_option "--localpair --maxiterate 1000 --reorder" --java_mem 4g --no_prefiltering --no_filtering --no_postfiltering --replace_FS_by_gaps

#Print length of alignment
for i in $(find . -name "mafft*" -type d); do paste <(echo $i) <(myfasta -l "$i"/apobec1_final_align_NT.aln | awk '{print$NF}' | sort -u); done | column -t
#Pairwise align full seq with filtered seq of Penelope_pileata
needle <(myfasta -mfp ../all_birds_trimmed_filtered.fa Penelope_pileata) <(myfasta -mfp mafft_no_prefiltering_no_filtering_no_postfiltering/apobec1_final_align_NT.aln Penelope_pileata | sed '/^>/! s/\-//g') --auto --stdout --awidth 240
#Multi seq align Penelope_pileata
mafft --auto --reorder --quiet <(cat <(find . -name "apobec1_final_align_NT.aln" -type f | grep mafft | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$2}" | sed "s/^/>/g"; myfasta -mfp $0 Penelope_pileata | grep -v ">"' | sed '/^>/! s/\-//g') <(myfasta -mfp ../all_birds_trimmed_filtered.fa Penelope_pileata | sed 's/>.*/>Full_seq/g')) | alv -

#Based on the options tried mafft_no_prefiltering_no_filtering_no_postfiltering is the best option: 
#Try different aligners with this parameters
time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa --out_dir muscle_no_prefiltering_no_filtering_no_postfiltering --out_file_prefix apobec1 \
 --alignAA_soft MUSCLE --java_mem 4g --no_prefiltering --no_filtering --no_postfiltering

time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa --out_dir prank_no_prefiltering_no_filtering_no_postfiltering --out_file_prefix apobec1 \
 --alignAA_soft PRANK --java_mem 4g --no_prefiltering --no_filtering --no_postfiltering

#Mark lost species for better view & inspection
	tools="mafft muscle prank"
	suffix="_no_prefiltering_no_filtering_no_postfiltering"
	input="apobec1_final_align_NT.aln"
	output="labelled_apobec1_final_align_NT.aln"
	#copy
	for t in $tools; do d=${t}${suffix}; cp "$d/$input" "$d/$output"; done
	#replace
	for group in galliformes palaeognathae independent
	do
		   file=~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/lost_${group}
		   while read -r i
		   do
		       for t in $tools; do sed -i "s/$i/&_${group}/g" ${t}${suffix}/$output; done
		   done < "$file"
	done

#NOTE:
	#Using muscle: Rhea is placed wrong, None of galliformes shows exon 2 in expected location
	#Using prank: Struthio shows exon_1 but lacks exon_1, None of galliformes shows exon 2 in expected location, exon_2 of galliformes is no aligning to any region in rest of the species while it should have aligned to exon_2

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Readd some previously removed species

#IMPORTANT: 
#I noticed that earlier I removed Cathartes_aura, Gymnogyps_californianus from ipnut fasta since it's exon_4 had longer 3' end. 
#But as these 2 species lost APOBEC1, to estimate their timing these shouldn't be removed, can be carefully trimmed such that the extended 3' end is atleast resonably aligned 
#Hence add these species to alignment.

mkdir ~/bird_db1/aswin/APOBEC1/Dating/alignment/readd
alen ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA/all_exon_4.aln

#Gymnogyps_californianus
	cat ~/bird_db1/aswin/APOBEC1/TOGA/Gymnogyps_californianus/All_consensus_msa.aln | grep -i catcagtaggcacttgaagaaaacactgaga -z
	cat ~/soft_links/Gymnogyps_californianus/aswin/APOBEC1/2nd_gblast/Query_exon_combined_exonerate.out
	exmalign2 -cds mafft ~/bird_db1/aswin/APOBEC1/TOGA/Gymnogyps_californianus/toga_rna-XM_027449866.2_exon_wise_consensus.fa ~/soft_links/Gymnogyps_californianus/aswin/APOBEC1/2nd_gblast/gblast_edited_consensus.fa \
		~/soft_links/Gymnogyps_californianus/aswin/APOBEC1/2nd_gblast/final_consensus.fa ~/soft_links/Gymnogyps_californianus/aswin/APOBEC1/2nd_gblast/spaln_exon_wise.fa ~/soft_links/Gymnogyps_californianus/aswin/APOBEC1/2nd_gblast/Query_exon_combined_exonerate.out_exon_wise.fa | colnum.sh 
	#Keep all exons except exon_4 from toga consensus, take exon_4 from exonerate
	cp ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA/Gymnogyps_californianus_toga.fa .
	cp ~/soft_links/Gymnogyps_californianus/aswin/APOBEC1/2nd_gblast/Query_exon_combined_exonerate.out_exon_wise.fa Gymnogyps_californianus_exonerate.fa
	myfasta -vfp Gymnogyps_californianus_toga.fa exon_4 > Gymnogyps_californianus_manual.fa
	#Add this manually
	myfasta -mfp Gymnogyps_californianus_exonerate.fa exon_4

#Cathartes_aura
	cat ~/bird_db1/aswin/APOBEC1/TOGA/Cathartes_aura/All_consensus_msa.aln | grep -i catcagtaggcacttgaagaaaacactgaga -z
	cat ~/soft_links/Cathartes_aura/aswin/APOBEC1/2nd_gblast/Query_exon_combined_exonerate.out
	exmalign2 -cds mafft ~/bird_db1/aswin/APOBEC1/TOGA/Cathartes_aura/toga_rna-XM_027449866.2_exon_wise_consensus.fa ~/soft_links/Cathartes_aura/aswin/APOBEC1/2nd_gblast/gblast_edited_consensus.fa \
		~/soft_links/Cathartes_aura/aswin/APOBEC1/2nd_gblast/final_consensus.fa ~/soft_links/Cathartes_aura/aswin/APOBEC1/2nd_gblast/spaln_exon_wise.fa ~/soft_links/Cathartes_aura/aswin/APOBEC1/2nd_gblast/Query_exon_combined_exonerate.out_exon_wise.fa | colnum.sh 
	#Keep all exons except exon_4 from toga consensus, take exon_4 from exonerate
	cp ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA/Cathartes_aura_toga.fa .
	cp ~/soft_links/Cathartes_aura/aswin/APOBEC1/2nd_gblast/Query_exon_combined_exonerate.out_exon_wise.fa Cathartes_aura_exonerate.fa
	myfasta -vfp Cathartes_aura_toga.fa exon_4 > Cathartes_aura_manual.fa
	#Add this manually
	myfasta -mfp Cathartes_aura_exonerate.fa exon_4

cd ~/bird_db1/aswin/APOBEC1/Dating/alignment/readd
cp ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered.fa
cat all_birds_trimmed_filtered.fa <(myfasta -cds Gymnogyps_californianus_manual.fa Cathartes_aura_manual.fa | sed 's/_manual//g') > all_birds_trimmed_filtered_readd.fa

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Rerun MACSE pipeline

	cd ~/bird_db1/aswin/APOBEC1/Dating/alignment
	cp ~/bird_db1/aswin/APOBEC1/Dating/alignment/readd/all_birds_trimmed_filtered_readd.fa .
	mkdir ~/bird_db1/aswin/APOBEC1/Dating/alignment/readd/readd_macse
	cd ~/bird_db1/aswin/APOBEC1/Dating/alignment/readd/readd_macse

#Run codon alignment
	time for aln in MAFFT MUSCLE PRANK
	do
	in=$(ls /home/neo/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered_readd.fa)
	echo ">MACSE pipeline using:" $aln
	time apptainer run /home/neo/programmes/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v12.02.sif --in_seq_file $in --out_dir "$aln"_no_prefiltering_no_filtering_no_postfiltering --out_file_prefix apobec1 --alignAA_soft $aln --java_mem 4g --no_prefiltering --no_filtering --no_postfiltering
	unset in
	done

#Mark lost species for better view & inspection
	tools="MAFFT MUSCLE PRANK"
	suffix="_no_prefiltering_no_filtering_no_postfiltering"
	input="apobec1_final_align_NT.aln"
	output="labelled_apobec1_final_align_NT.aln"
	#copy
	for t in $tools; do d=${t}${suffix}; cp "$d/$input" "$d/$output"; done
	#replace
	for group in galliformes palaeognathae independent
	do
		   file=~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/lost_${group}
		   while read -r i
		   do
		       for t in $tools; do sed -i "s/$i/&_${group}/g" ${t}${suffix}/$output; done
		   done < "$file"
	done

#Confirm the NT & codon alignment follows the intact APOBEC1's frame
myfasta -mfl ~/bird_db1/aswin/APOBEC1/Dating/alignment/all_birds_trimmed_filtered_readd.fa ~/bird_db1/aswin/APOBEC1/Dating/cds/intact_apobec1 | transeq --auto --stdout -clean stdin | myfasta -comb | sed '/^>/ s/_1//g' > intact_birds_filtered.aa
mafft --maxiterate 1000 --localpair --reorder --quiet intact_birds_filtered.aa > intact_birds_filtered.aln
em_cons intact_birds_filtered.aln --auto --stdout | myfasta -comb | sed 's/^>.*/>inatct_filtered/g' > intact_birds_filtered_consensus.fa

em_cons MAFFT_no_prefiltering_no_filtering_no_postfiltering/apobec1_final_align_AA.aln --auto --stdout | myfasta -comb | sed 's/^>.*/>mafft/g' > mafft_macse_pipeline_consensus.fa
em_cons MUSCLE_no_prefiltering_no_filtering_no_postfiltering/apobec1_final_align_AA.aln --auto --stdout | myfasta -comb | sed 's/^>.*/>muscle/g' > muscle_macse_pipeline_consensus.fa
em_cons PRANK_no_prefiltering_no_filtering_no_postfiltering/apobec1_final_align_AA.aln --auto --stdout | myfasta -comb | sed 's/^>.*/>prank/g' > prank_macse_pipeline_consensus.fa
cat intact_birds_filtered_consensus.fa  mafft_macse_pipeline_consensus.fa  muscle_macse_pipeline_consensus.fa  prank_macse_pipeline_consensus.fa > intact_and_macse_pipeline_consensus.fa
mafft --maxiterate 1000 --localpair --reorder --quiet intact_and_macse_pipeline_consensus.fa > intact_and_macse_pipeline_consensus.aln
mafft --maxiterate 1000 --localpair --reorder --quiet <(myfasta -vfp intact_and_macse_pipeline_consensus.fa prank) | alv -
#NOTE: Despite some gaps, globally the amino acids are similar between consensus

#NOTE:
	#Upon inspection prank had worst placement of exons & gaps, while mafft & muscle were almost similar except in case of palaeognathae mafft was slightly better. 
	#E.g. Rhea has only exon_5, mafft place this at end of alignment, while muscle places thios exon_5 in the centre.
	#Hence choose this alignment : /home/neo/bird_db1/aswin/APOBEC1/Dating/alignment/readd/readd_macse/MAFFT_no_prefiltering_no_filtering_no_postfiltering/apobec1_final_align_NT.aln

#Convert alignment format to phylip
num=$( grep '>' MAFFT_no_prefiltering_no_filtering_no_postfiltering/apobec1_final_align_NT.aln | wc -l )
len=$( sed -n '2,2p' MAFFT_no_prefiltering_no_filtering_no_postfiltering/apobec1_final_align_NT.aln | sed 's/\r//' | sed 's/\n//' | wc -L )
perl /home/neo/programmes/paml-tutorial/positive-selection/00_data/scripts/FASTAtoPHYL.pl MAFFT_no_prefiltering_no_filtering_no_postfiltering/apobec1_final_align_NT.aln $num $len 

num=$( grep '>' data1_nuc_aln.fasta | wc -l )
len=$( sed -n '2,2p' data1_nuc_aln.fasta | sed 's/\r//' | sed 's/\n//' | wc -L )
perl ../../scripts/FASTAtoPHYL.pl data1_nuc_aln.fasta $num $len 
mv data1_nuc_aln.phy data1.phy

#Hence final alignment to use : /home/neo/bird_db1/aswin/APOBEC1/Dating/alignment/readd/readd_macse/apobec1_final_align_NT.aln.phy

############################################################################################################################################################################################################################################################################################################
#Get tree

cd ~/bird_db1/aswin/APOBEC1/Dating/cds
mkdir ~/bird_db1/aswin/APOBEC1/Dating/tree
grep ">" all_birds_filtered.fa | tr -d ">" | tr "_" " " > ~/bird_db1/aswin/APOBEC1/Dating/tree/list_all_birds_filtered
cd ~/bird_db1/aswin/APOBEC1/Dating/tree

awk -F' \\(' '{gsub(/ /,"_",$1)
    print $1 " (" $2}' replace > replace2

cat all_birds_filtered_raw.nwk | sed -e 's/Tinamus/Tinamus_guttatus/g' -e 's/Chaetura/Chaetura_pelagica/g' -e 's/Egretta novaehollandiae/Egretta_garzetta/g' -e 's/Podilymbus/Podilymbus_podiceps/g' -e 's/Phoenicopterus/Phoenicopterus_ruber/g' -e 's/Merops apiaster/Merops_nubicus/g' -e 's/Galbula/Galbula_dea/g' -e 's/Colius/Colius_striatus/g' -e 's/Aquila chrysaetos/Aquila_chrysaetos_chrysaetos/g' -e 's/Lonchura striata/Lonchura_striata_domestica/g' -e 's/Malurus cyaneus/Malurus_cyaneus_samueli/g' > all_birds_filtered.nwk 
cp ~/bird_db1/aswin/APOBEC1/Dating/cds/all_birds_filtered.fa .

#Check if IDs differ b/w alignment & tree
colordiff -y <(cat all_birds_filtered.nwk | tr -d ")(;.[0-9]:'" | tr "," "\n" | sort -k1 -V) <(grep ">" all_birds_filtered.fa | tr -d ">" | sort -k1V)
colordiff -y <(grep ">" all_birds_filtered.fa | tr -d ">" | sort -k1V) <(cat all_birds_filtered.nwk | tr -d ")(;.[0-9]:'" | tr "," "\n" | sort -k1 -V)

#rename names in tree to match alignment
sed -e 's/Chaetura_pelagica_pelagica/Chaetura_pelagica/g' -e 's/Colius_striatus_striatus/Colius_striatus/g' -e 's/Chloebia_gouldiae/Erythrura_gouldiae/g' -e 's/Cyanoderma_ruficeps/Stachyris_ruficeps/g' -e 's/Chloebia_gouldiae/Erythrura_gouldiae/g' \
 -e 's/Galbula_dea_dea/Galbula_dea/g' -e 's/Phoenicopterus_ruber_ruber/Phoenicopterus_ruber/g' -e 's/Podilymbus_podiceps_podiceps/Podilymbus_podiceps/g' -e 's/Tinamus_guttatus_guttatus/Tinamus_guttatus/g' all_birds_filtered.nwk > all_birds_filtered_renamed.nwk

colordiff -y <(grep ">" all_birds_filtered.fa | tr -d ">" | sort -k1V) <(cat all_birds_filtered_renamed.nwk | tr -d ")(;.[0-9]:'" | tr "," "\n" | sort -k1 -V)

#sed "s/'[0-9]\+'//g" all_birds_filtered_renamed.nwk > all_birds_filtered_renamed_edited.nwk

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Update tree with final alignment 

#Get species names
mkdir ~/bird_db1/aswin/APOBEC1/Dating/tree/readd
cd ~/bird_db1/aswin/APOBEC1/Dating/tree/readd
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/alignment/readd/readd_macse/MAFFT_no_prefiltering_no_filtering_no_postfiltering/apobec1_final_align_NT.aln .
grep ">" apobec1_final_align_NT.aln | tr -d ">" | tr "_" " " | cut -f1,2 -d " " > list_apobec1_final_align_NT

#Use timetree website
mv list_apobec1_final_align_NT.nwk apobec1_final_align_NT_raw.nwk
mv prunetree.jpg apobec1_final_align_NT_tree.pdf
awk -F' \\(' '{gsub(/ /,"_",$1)
    print $1 " (" $2}' replace > replace2

#Replace names changed by timetree
timetree_replace.sh replace2
cat apobec1_final_align_NT_raw.nwk | sed -e 's/Tinamus/Tinamus_guttatus/g' -e 's/Chaetura/Chaetura_pelagica/g' -e 's/Egretta novaehollandiae/Egretta_garzetta/g' -e 's/Podilymbus/Podilymbus_podiceps/g' -e 's/Phoenicopterus/Phoenicopterus_ruber/g' \
 -e 's/Merops apiaster/Merops_nubicus/g' -e 's/Galbula/Galbula_dea/g' -e 's/Colius/Colius_striatus/g' > apobec1_final_align_NT_replaced.nwk

#Compare tree & alignment
colordiff -y <(grep ">" apobec1_final_align_NT.aln | tr -d ">" | sort -k1V) <(cat apobec1_final_align_NT_replaced.nwk | tr -d ")(;.[0-9]:'" | tr "," "\n" | sort -k1 -V)

#Rename slight differences in names
cat apobec1_final_align_NT_replaced.nwk | sed -e 's/Aquila_chrysaetos/Aquila_chrysaetos_chrysaetos/g' -e 's/Chaetura_pelagica_pelagica/Chaetura_pelagica/g' -e 's/Chloebia_gouldiae/Erythrura_gouldiae/g' -e 's/Colius_striatus_striatus/Colius_striatus/g' -e 's/Cyanoderma_ruficeps/Stachyris_ruficeps/g' \
 -e 's/Galbula_dea_dea/Galbula_dea/g' -e 's/Lonchura_striata/Lonchura_striata_domestica/g' -e 's/Malurus_cyaneus/Malurus_cyaneus_samueli/g' -e 's/Phoenicopterus_ruber_ruber/Phoenicopterus_ruber/g' -e 's/Podilymbus_podiceps_podiceps/Podilymbus_podiceps/g' \
 -e 's/Tinamus_guttatus_guttatus/Tinamus_guttatus/g' > apobec1_final_align_NT.nwk

#Compare tree & alignment
colordiff -y <(grep ">" apobec1_final_align_NT.aln | tr -d ">" | sort -k1V) <(cat apobec1_final_align_NT.nwk | tr -d ")(;.[0-9]:'" | tr "," "\n" | sort -k1 -V)

sed "s/'[0-9]\+'//g" apobec1_final_align_NT.nwk > apobec1_final_align_NT_without_split_times.nwk
sed 's/[0-9]//g' apobec1_final_align_NT_without_split_times.nwk | tr -d ":." > apobec1_final_align_NT_without_split_times_branch_lengths.nwk

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Unroot tree
cd ~/bird_db1/aswin/APOBEC1/Dating/tree/readd
Rscript /home/neo/bird_db1/aswin/APOBEC1/Dating/scripts/unroot_tree.R apobec1_final_align_NT.nwk

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Label tree
cd ~/bird_db1/aswin/APOBEC1/Dating/cds
grep -vf lost_species/TOGA/species_to_remove intact_apobec1 > intact_apobec1_filtered
cd ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species
grep -v Meleagris_gallopavo lost_galliformes > lost_galliformes_Meleagris_gallopavo_removed
cp lost_galliformes_Meleagris_gallopavo_removed lost_palaeognathae lost_independent ../intact_apobec1_filtered ~/bird_db1/aswin/APOBEC1/Dating/tree/readd/

cd ~/bird_db1/aswin/APOBEC1/Dating/tree/readd
cat lost_independent > fg1.txt
cat lost_galliformes_Meleagris_gallopavo_removed lost_palaeognathae > fg2.txt 
cat intact_apobec1_filtered > bg.txt
Rscript /home/neo/bird_db1/aswin/APOBEC1/Dating/scripts/label_tree.R apobec1_final_align_NT_unroot.nwk fg1.txt fg2.txt


#Manual labelling of branches based on loss & mixed using hyphy website: https://phylotree.hyphy.org/#
#NOTE: label in the following way & export the tree & save in the name "apobec1_final_align_NT_unroot_hyphy_labelled.nwk"
	#Pseudogene branch (even sister clade has loss) as ps
	#Mixed branch (sister clade has intact gene) as mi ()

cd ~/bird_db1/aswin/APOBEC1/Dating/tree/readd
sed 's/{mi}/ #1/g' apobec1_final_align_NT_unroot_hyphy_labelled.nwk | sed 's/{ps}/ #2/g' > apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk
sed 's/:[0-9.]\+//g' apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk > apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml_branch_labels_removed.nwk


#The old final tree is : /home/neo/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_labeled.nwk
#The new final tree is : /home/neo/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk

############################################################################################################################################################################################################################################################################################################
#PAML

mkdir ~/bird_db1/aswin/APOBEC1/Dating/paml
cd ~/bird_db1/aswin/APOBEC1/Dating/paml

#Get inputs
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/alignment/readd/readd_macse/MAFFT_no_prefiltering_no_filtering_no_postfiltering/apobec1_final_align_NT.aln .
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT.nwk .

cp ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/lost_* .
sed '/Meleagris_gallopavo/d' lost_galliformes -i
cat lost_galliformes lost_palaeognathae lost_independent > all_lost

#Classify branches
grep -vf ~/bird_db1/aswin/APOBEC1/Dating/cds/lost_species/TOGA/species_to_remove ~/bird_db1/aswin/APOBEC1/Dating/cds/intact_apobec1 > funtional

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Check scripts used in COA1 gene

git clone --filter=blob:none --no-checkout https://github.com/ceglabsagarshinde/COA1_GENE.git
cd COA1_GENE
git sparse-checkout init --cone
git sparse-checkout set Geneloss_timing
git checkout 01b43da88117ed8249d073f9e1b78ddb4732ba83

cd ~/bird_db1/aswin/APOBEC1/Dating/paml/COA1_GENE
find Geneloss_timing/Galliformes/ -name "*.ctl" -type f | sort -V | paste -d " " - - | xargs -n2 sh -c 'echo ">"$0 $1; paste $0 $1 | nl' | less
find Geneloss_timing/Galliformes/ -name "*.ctl" -type f | sort -V | paste -d " " - - | xargs -n2 sh -c 'echo ">"$0 $1; paste $0 $1 -d "\t" | tr " " "_"' | sed 's!Geneloss_timing/Galliformes/!!g' | colnum.sh
find Geneloss_timing/Galliformes/ -name "*.ctl" -type f | sort -V | paste -d " " - - | xargs -n2 sh -c 'echo ">"$0 $1; paste $0 $1' | sed 's!Geneloss_timing/Galliformes/!!g' > all_galliformes_control_files
cat Geneloss_timing/Galliformes/AllPseudogene_AllMix_AllFunctional/Galliformes_F1x4.ctl| cut -f1 -d "=" | tr -d " " > ctl_params

for p in $(cat ctl_params | egrep -v "seqfile|treefile|outfile" | sort -V)
do
echo ">"$p
grep -w $p all_galliformes_control_files | sort | uniq -c
done | less

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run codeml

mkdir -p ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/F1X4
cd ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/F1X4
cat /home/neo/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_labeled.nwk | sed 's/:[0-9.]\+//g' | sed 's/#/ #/g' > apobec1_final_align_NT_unroot_labeled_branch_labels_removed.nwk

BASE=~/bird_db1/aswin/APOBEC1/Dating
TREE="$BASE/tree/readd/apobec1_final_align_NT_unroot_labeled.nwk"
ALN="$BASE/alignment/readd/readd_macse/apobec1_final_align_NT.aln.phy"
CTL="$BASE/paml/COA1_GENE/Geneloss_timing/Galliformes/AllPseudogene_AllMix_AllFunctional"
CODEML=~/programmes/paml-4.10.10-linux-x86_64/bin/codeml
OUT="$BASE/paml/Independent_loss_as_label1_Group_loss_as_label2_Intact_as_unlabelled"

# make branch-label removed tree
mkdir -p "$OUT"
sed 's/:[0-9.]\+//g; s/#/ #/g' "$TREE" > "$OUT/apobec1_final_align_NT_unroot_labeled_nobranch.nwk"

for m in F1x4 F3x4; do
  for t in "$TREE" "$OUT/apobec1_final_align_NT_unroot_labeled_nobranch.nwk"; do
    [[ "$t" == *nobranch* ]] && suf="branch_labels_removed_$m" || suf="$m"
    d="$OUT/$suf"; mkdir -p "$d"; cd "$d"
    cp "$t" "$ALN" .
    cp "$CTL/Galliformes_${m^}.ctl" "$m.ctl"
    sed -i -e 's|Galliformes.aln|apobec1_final_align_NT.aln.phy|g' -e "s|Galliformes.nwk|$(basename "$t")|g" -e "s|omega_mix_functional_${m^}|paml_out_$m|g" "$m.ctl"
    time "$CODEML" "$m.ctl" > run.stdout &
  done
done
wait


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run codeml based on tree manually lablled using hyphy tool later converted to paml

/home/neo/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk

cd ~/bird_db1/aswin/APOBEC1/Dating/paml

BASE=~/bird_db1/aswin/APOBEC1/Dating
TREE="$BASE/tree/readd/apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk"
ALN="$BASE/alignment/readd/readd_macse/apobec1_final_align_NT.aln.phy"
CTL="$BASE/paml/COA1_GENE/Geneloss_timing/Galliformes/AllPseudogene_AllMix_AllFunctional"
CODEML=~/programmes/paml-4.10.10-linux-x86_64/bin/codeml
OUT="$BASE/paml/Pseudo_as_label1_Mixed_as_label2_Inatct_as_unlabelled"

# make branch-label removed tree
mkdir -p "$OUT"
sed 's/:[0-9.]\+//g' "$TREE" > "$OUT/apobec1_final_align_NT_unroot_labeled_nobranch.nwk"

start_time=$(date +%s)
for m in F1x4 F3x4; do
  for t in "$TREE" "$OUT/apobec1_final_align_NT_unroot_labeled_nobranch.nwk"; do
    [[ "$t" == *nobranch* ]] && suf="branch_labels_removed_$m" || suf="$m"
    d="$OUT/$suf"; mkdir -p "$d"; cd "$d"
    cp "$t" "$ALN" .
    cp "$CTL/Galliformes_${m^}.ctl" "$m.ctl"
    sed -i -e 's|Galliformes.aln|apobec1_final_align_NT.aln.phy|g' -e "s|Galliformes.nwk|$(basename "$t")|g" -e "s|omega_mix_functional_${m^}|paml_out_$m|g" "$m.ctl"
    time "$CODEML" "$m.ctl" > run.stdout &
  done
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e



#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calculate gene inactivation times
/home/neo/bird_db1/aswin/APOBEC1/Dating/paml/Mammal_ADH_IV/Dating/codeml/F3X4_model/calculate_gene_inactivation.sh


cd ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix
cp ~/bird_db1/aswin/APOBEC1/Dating/tree/readd/{fg1.txt,fg2.txt,bg.txt} .

cd ~/bird_db1/aswin/APOBEC1/Dating/paml/Pseudo_as_label1_Mixed_as_label2_Inatct_as_unlabelled
for i in $(cat ../funtional)
do
~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh Gallus_gallus $i F1x4/paml_out_F1x4 F3x4/paml_out_F3x4 -wp=1 ~/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk -s | grep -v "Mixed_branch_length"
done | sed '1i Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' | column -t > gallus_gallus_gene_loss_date_wrt_diff_functional_branches.out

time for sp in $(cat ../all_lost)
do
gr=$(grep $sp ~/bird_db1/aswin/taxonomy/orders_all_birds | awk '{print$2}')
~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh $sp -f F1x4/paml_out_F1x4 F3x4/paml_out_F3x4 -wp=1 ~/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk -s | grep -v Mixed_branch_length | sed "s/^/$gr\t/g"
unset gr
done | sed '1i Group Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' > all_gene_loss_dates.tsv



#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#DRAFT SCRIPTS:

#Outgroup

#NOTE: An outgroup is not necessarily required to run gene loss dating
#The closest sequences to the A1 of birds (not A1-like from birds) is identifiable from the phylogenetic tree build based on cds of all AID-APOBEC members from vertebrates
#Some variations existed in terms of which sequences are closest to birds, as phylogeny built from different strategies or aligners gave slightly different sister clades.
#The tree & alignment used directly for the main figure is used :

mkdir ~/bird_db1/aswin/APOBEC1/Dating/cds/outgroup
cd ~/bird_db1/aswin/APOBEC1/Dating/cds/outgroup
cp /home/neo/bird_db1/aswin/APOBEC1/positive_selection/A1_A1_like_all_vertebrates/AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa .

for i in crocodylia  squamata  testudines
do
myfasta -mfl AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa $i > "$i"_cds.fa
done

#Alignment
for i in crocodylia_cds.fa  squamata_cds.fa  testudines_cds.fa
do
mafft --maxiterate 1000 --localpair --reorder --quiet $i > $i".aln"
done

cat crocodylia_cds.fa squamata_cds.fa testudines_cds.fa > all_outgroup_cds.fa
mafft --maxiterate 1000 --localpair --reorder --quiet all_outgroup_cds.fa > all_outgroup_cds.aln

mafft --maxiterate 1000 --localpair --reorder --quiet all_outgroup.aa > all_outgroup.aln
mafft --maxiterate 1000 --localpair --reorder --quiet birds_outgroup.aa > birds_outgroup.aln

#Outliers in manual inspection of aa alignment
ensembl_A1_paralog_Squamata_ENSPMRG00000022315
orthodb_A1_paralog_Squamata_Podarcis_muralis14588080
ensembl_A1_paralog_Testudines_ENSCPBG00000014331
ensembl_A1_paralog_Squamata_ENSACAG00000035333
orthodb_A1_paralog_Squamata_Python_bivittatus03048122
ensembl_A1_paralog_Crocodylia_ENSCPRG00005000694
ensembl_A1_paralog_Squamata_ENSSMRG00000004573
orthodb_A1_paralog_Squamata_Sphaerodactylus_townsendi25440208
orthodb_A1_paralog_Squamata_Anolis_carolinensis00534649
orthodb_A1_paralog_Squamata_Anolis_carolinensis03277744
orthodb_A1_paralog_Testudines_Mauremys_reevesii20372170

#Outliers in manual inspection of nucleotide alignment
orthodb_A1_paralog_Squamata_Podarcis_muralis_114588080
ensembl_A1_paralog_Squamata_ENSPMRG00000022315
ensembl_A1_paralog_Squamata_ENSACAG00000035333
ensembl_A1_paralog_Squamata_ENSSMRG00000004573
ensembl_A1_paralog_Testudines_ENSCPBG00000014331

myfasta -vfl all_outgroup.aa outliers1 | mafft --maxiterate 1000 --localpair --reorder --quiet - | alv -
myfasta -vfl all_outgroup.aa outliers2 | mafft --maxiterate 1000 --localpair --reorder --quiet - | alv -
cat <(myfasta -vfl all_outgroup.aa outliers2) all_birds_trimmed_filtered_readd.aa | mafft --maxiterate 1000 --localpair --reorder --quiet - | alv -

#NOTE: Squamata alignment has a lot of gaps

#Clustering
mkdir ~/bird_db1/aswin/APOBEC1/Dating/cds/outgroup/clustering
cd ~/bird_db1/aswin/APOBEC1/Dating/cds/outgroup
cat crocodylia_cds.fa squamata_cds.fa testudines_cds.fa | transeq --auto --stdout --clean stdin | myfasta -comb | sed 's/_1//g' > clustering/all_outgroup.aa
cat ~/bird_db1/aswin/APOBEC1/Dating/alignment/readd/all_birds_trimmed_filtered_readd.fa | transeq --auto --stdout --clean stdin | myfasta -comb | sed 's/_1//g' > all_birds_trimmed_filtered_readd.aa
cat all_outgroup.aa all_birds_trimmed_filtered_readd.aa > birds_outgroup.aa

#scp 9831869.clans.zip neo@172.28.65.224:~/bird_db1/aswin/APOBEC1/Dating/cds/outgroup/clustering/
java -Xmx4G -jar ~/programmes/clans.jar birds_outgroup.clans

#
mkdir ~/bird_db1/aswin/APOBEC1/Dating/cds/outgroup/similarity_to_bird_consensus
cp ~/bird_db1/aswin/APOBEC1/Dating/cds/all_intact_species_cds_filtered.fa .
cp ../all_outgroup_cds.fa .

em_cons MAFFT_no_prefiltering_no_filtering_no_postfiltering/apobec1_final_align_AA.aln --auto --stdout | myfasta -comb | sed 's/^>.*/>mafft/g' > mafft_macse_pipeline_consensus.fa

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Paml test run

#Input files
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_labeled.nwk /home/neo/bird_db1/aswin/APOBEC1/Dating/alignment/readd/readd_macse/apobec1_final_align_NT.aln.phy .
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/paml/COA1_GENE/Geneloss_timing/Galliformes/AllPseudogene_AllMix_AllFunctional/Galliformes_F1x4.ctl F1x4.ctl 
sed -e 's/Galliformes.aln/apobec1_final_align_NT.aln.phy/g' -e 's/Galliformes.nwk/apobec1_final_align_NT_unroot_labeled.nwk/g' -e 's/omega_mix_functional_F1x4/paml_out_F1x4/g' F1x4.ctl -i
#Run codeml
time ~/programmes/paml-4.10.10-linux-x86_64/bin/codeml F1x4.ctl > run.stdout &

mkdir -p ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/F3X4
cd ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/F3X4
#Input files
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_labeled.nwk /home/neo/bird_db1/aswin/APOBEC1/Dating/alignment/readd/readd_macse/apobec1_final_align_NT.aln.phy .
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/paml/COA1_GENE/Geneloss_timing/Galliformes/AllPseudogene_AllMix_AllFunctional/Galliformes_F3x4.ctl F3x4.ctl 
sed -e 's/Galliformes.aln/apobec1_final_align_NT.aln.phy/g' -e 's/Galliformes.nwk/apobec1_final_align_NT_unroot_labeled.nwk/g' -e 's/omega_mix_functional_F3x4/paml_out_F3x4/g' F3x4.ctl -i
#Run codeml (47m37.113s)
#delete previous paml output: rm rub rst1 rst lnf 2NG.dN  2NG.dS  2NG.t
time ~/programmes/paml-4.10.10-linux-x86_64/bin/codeml F3x4.ctl > run.stdout &


cd ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix
cat /home/neo/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_labeled.nwk | sed 's/:[0-9.]\+//g' | sed 's/#/ #/g' > apobec1_final_align_NT_unroot_labeled_branch_labels_removed.nwk

mkdir -p ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/branch_labels_removed_F1X4
cd ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/branch_labels_removed_F1X4
#Input files
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/apobec1_final_align_NT_unroot_labeled_branch_labels_removed.nwk /home/neo/bird_db1/aswin/APOBEC1/Dating/alignment/readd/readd_macse/apobec1_final_align_NT.aln.phy .
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/paml/COA1_GENE/Geneloss_timing/Galliformes/AllPseudogene_AllMix_AllFunctional/Galliformes_F1x4.ctl F1x4.ctl 
sed -e 's/Galliformes.aln/apobec1_final_align_NT.aln.phy/g' -e 's/Galliformes.nwk/apobec1_final_align_NT_unroot_labeled_branch_labels_removed.nwk/g' -e 's/omega_mix_functional_F3x4/paml_out_F1x4/g' F1x4.ctl -i
#Run codeml
#delete previous paml output: rm rub rst1 rst lnf 2NG.dN  2NG.dS  2NG.t
time ~/programmes/paml-4.10.10-linux-x86_64/bin/codeml F3x4.ctl > run.stdout &

mkdir -p ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/branch_labels_removed_F3X4
cd ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/branch_labels_removed_F3X4
#Input files
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/apobec1_final_align_NT_unroot_labeled_branch_labels_removed.nwk /home/neo/bird_db1/aswin/APOBEC1/Dating/alignment/readd/readd_macse/apobec1_final_align_NT.aln.phy .
cp /home/neo/bird_db1/aswin/APOBEC1/Dating/paml/COA1_GENE/Geneloss_timing/Galliformes/AllPseudogene_AllMix_AllFunctional/Galliformes_F3x4.ctl F3x4.ctl 
sed -e 's/Galliformes.aln/apobec1_final_align_NT.aln.phy/g' -e 's/Galliformes.nwk/apobec1_final_align_NT_unroot_labeled_branch_labels_removed.nwk/g' -e 's/omega_mix_functional_F3x4/paml_out_F3x4/g' F3x4.ctl -i
#Run codeml
#delete previous paml output: rm rub rst1 rst lnf 2NG.dN  2NG.dS  2NG.t
time ~/programmes/paml-4.10.10-linux-x86_64/bin/codeml F3x4.ctl > run.stdout &


