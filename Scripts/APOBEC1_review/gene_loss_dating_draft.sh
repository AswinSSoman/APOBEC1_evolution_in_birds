
############################################################################################################################################################################################################################################################################################################
#TEST ADH gene


#Get example codes & input
cd /media/aswin/gene_loss/APOBEC1/selection/positive_selection/paml
git clone https://github.com/Dating-gene-loss/Mammal_ADH_IV.git

#Prepare folders
cd /media/aswin/gene_loss/APOBEC1/selection/positive_selection/paml/Mammal_ADH_IV
cd /media/aswin/gene_loss/APOBEC1/selection/positive_selection/paml/Mammal_ADH_IV/Dating/codeml/F1X4_model/hyphy

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



#Prepare folders
mkdir ~/bird_db1/aswin/APOBEC1/Dating
cd ~/bird_db1/aswin/APOBEC1/Dating

#Get input sequences
mkdir ~/bird_db1/aswin/APOBEC1/Dating/cds
cp ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*.fa ~/bird_db1/aswin/APOBEC1/Dating/cds

cd ~/bird_db1/aswin/APOBEC1/Dating/cds
ls APOBEC1_*.fa | cut -f2- -d "_" | sed 's/\.fa//g' > intact_apobec1

#exon-wise separate sequence
mkdir /home/neo/bird_db1/aswin/APOBEC1/Dating/cds/exon_wise
while read -r sp; do
echo "$sp"
f=$(ls | grep -m1 "$sp")
for exon in {1..5}; do
myfasta -mfp "$f" "exon_${exon}" | sed "s/^>.*/>$sp/" >> "exon_wise/exon_${exon}.fa"
done
done < intact_apobec1

for i in $(ls exon*.fa)
do
echo ">"$i
mafft --maxiterate 1000 --localpair --reorder --quiet $i > $i".aln"
done


cp /home/neo/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.fa ~/bird_db1/aswin/APOBEC1/Dating
myfasta -l /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Acanthisitta_chloris.fa

grep galliformes -i ~/bird_db1/aswin/APOBEC1/pipeline_summary_row_wise | awk '{print$1}' | tr -d ">" | grep -v "Alectura_lathami" > lost_galliformes
myfasta -mfl toga_subjects.fa lost_palaeognathae | myfasta -l | sort -k2nr,2
myfasta -mfl toga_subjects.fa lost_independent | myfasta -l | sort -k2nr,2

############################################################################################################################################################################################################################################################################################################

#Manually check accuracy of consensus sequence of Struthio_camelus_australis before using for gene loss dating
cd ~/soft_links/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/manual_search

#Check how inframe eoxn_4 should start and end
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

grep -f all_lost_birds_with_remains ~/bird_db1/aswin/APOBEC1/pipeline_summary_column_wise | colnum.sh

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

for sp in $(find TOGA/ -name "*.fa" -type f)

############################################################################################################################################################################################################################################################################################################





