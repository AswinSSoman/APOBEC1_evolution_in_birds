
##########################################################################################################################################################################################################################################################################################################
#
##########################################################################################################################################################################################################################################################################################################


##########################################################################################################################################################################################################################################################################################################
#Comment 3 

#apoB alignment:
mkdir ~/bird_db1/aswin/APOBEC1/function/APOB/other_vertebrates
cd ~/bird_db1/aswin/APOBEC1/function/APOB/other_vertebrates

time datasets download gene symbol apob --include cds --ortholog crocodylia --filename crocodylia.zip
time datasets download gene symbol apob --include cds --ortholog testudines --filename testudines.zip
time datasets download gene symbol apob --include cds --ortholog squamata --filename squamata.zip
unzip crocodylia.zip -d crocodylia
unzip testudines.zip -d testudines
unzip squamata.zip -d squamata

cd ~/bird_db1/aswin/APOBEC1/function/APOB/other_vertebrates/squamata/ncbi_dataset/data
awk '/^>/{
    acc=$1; sub(/^>/,"",acc); sub(/:.*/,"",acc)
    gene=$2
    match($0,/\[organism=([^]]+)\]/,o); org=o[1]
    match($0,/\[GeneID=([0-9]+)\]/,g); gid=g[1]
    split(org,a," ")
    print ">"gene"_"acc"_"org"_"gid"_squamata"
    next} {print}' cds.fna > squamata.fa
time mafft --maxiterate 1000 --localpair --reorder --thread 4 squamata.fa > squamata.aln 
#sequences to remove:
#APOB_XM_014055581.1_Thamnophis sirtalis_106540476_sirtalis_squamata: only 1 exon & gaps preent
myfasta -vpp squamata.fa APOB_XM_014055581 > squamata_filtered.fa
awk -F"_" '/^>/ {print ">S_"$4; next} {print}' squamata_filtered.fa | sed '/^>/ s/ /_/g' | myfasta -comb > squamata_filtered_renamed.fa

cd ~/bird_db1/aswin/APOBEC1/function/APOB/other_vertebrates/crocodylia/ncbi_dataset/data
awk '/^>/{
    acc=$1; sub(/^>/,"",acc); sub(/:.*/,"",acc)
    gene=$2
    match($0,/\[organism=([^]]+)\]/,o); org=o[1]
    match($0,/\[GeneID=([0-9]+)\]/,g); gid=g[1]
    split(org,a," ")
    print ">"gene"_"acc"_"org"_"gid"_crocodylia"
    next} {print}' cds.fna > crocodylia.fa
time mafft --maxiterate 1000 --localpair --reorder --thread 4 crocodylia.fa > crocodylia.aln 
#sequence to remove:
#APOB_XM_006034385.3_Alligator sinensis_102383976_sinensis_crocodylia : only 9 exons & gaps present & apoB is at assembly end
myfasta -vpp crocodylia.fa APOB_XM_00603438 > crocodylia_filtered.fa
awk -F"_" '/^>/ {print ">C_"$4; next} {print}' crocodylia_filtered.fa | sed '/^>/ s/ /_/g' | myfasta -comb > crocodylia_filtered_renamed.fa

cd ~/bird_db1/aswin/APOBEC1/function/APOB/other_vertebrates/testudines/ncbi_dataset/data
awk '/^>/{
    acc=$1; sub(/^>/,"",acc); sub(/:.*/,"",acc)
    gene=$2
    match($0,/\[organism=([^]]+)\]/,o); org=o[1]
    match($0,/\[GeneID=([0-9]+)\]/,g); gid=g[1]
    split(org,a," ")
    print ">"gene"_"acc"_"org"_"gid"_testudine"
    next} {print}' cds.fna > testudine.fa
time mafft --maxiterate 1000 --localpair --reorder --thread 4 testudine.fa > testudine.aln 
#sequence to remove:
#APOB_XM_006034385.3_Alligator sinensis_102383976_sinensis_testudine : only 9 exons & gaps present & apoB is at assembly end
awk -F"_" '/^>/ {print ">T_"$4; next} {print}' testudine.fa | sed '/^>/ s/ /_/g' | myfasta -comb > testudine_renamed.fa
time mafft --maxiterate 1000 --localpair --reorder --thread 4 testudine_renamed.fa > testudine_renamed.aln 

#Collect all sequence for final alignment
cd ~/bird_db1/aswin/APOBEC1/function/APOB/other_vertebrates
cp ../APOB_representative_mammals_one_per_gene_intact_apobec1_validated_birds_species_manually_checked_APOB_edited.fa .
find . -name "*renamed.fa" | xargs -n1 sh -c 'cp $0 .'
cat APOB_representative_mammals_one_per_gene_intact_apobec1_validated_birds_species_manually_checked_APOB_edited.fa testudine_renamed.fa squamata_filtered_renamed.fa crocodylia_filtered_renamed.fa > all_vertebrates_apob.fa

#Label using jalview alignment
#Load alignment -> change formats: right align sequence id, Show nonconserved -> Set human as reference -> locate motifs : e.g. mooring sequnces
	# -> select the mooring sequence of human & create sequence feature -> select all species aligning in this mooring sequence & create new group 
	# -> group colour based on percentage identity -> select the whole animal group sequence ids & create a new group, add names, e.g. Mammals.

M_Equus_caballus
B_Corvus_moneduloides

B_Limosa_lapponica_baueri

B_Buceros_rhinoceros_silvestris
B_Opisthocomus_hoazin

M_Canis_lupus_familiaris
M_Ursus_maritimus
M_Neovison_vison
M_Felis_catus
M_Panthera_leo
M_Manis_pentadactyla
M_Trichechus_manatus_latirostris
M_Elephas_maximus_indicus
S_Eublepharis_macularius

M_Erinaceus_europaeus
B_Opisthocomus_hoazin



##########################################################################################################################################################################################################################################################################################################

cd ~/bird_db1/aswin/database_details

#species & order count : filtered
grep -f <(less apobec1_supp | sed 's/ /_/g' | awk '{print$1}') ../taxonomy/orders_all_birds -z | awk '$2' | awk '{a[$1]++; b[$2]++} {print length(a),$1,length(b),$2}' | colnum.sh 

#species & order count : before filtering
grep -f <(less apobec1_supp | sed 's/ /_/g' | awk '{print$1}') ../taxonomy/orders_all_birds -z | awk '$2' | awk '{a[$1]++; b[$2]++} {print length(a),$1,length(b),$2}' | colnum.sh 



#Get bed file for custom track in NCBI genome data viewer
sed 's/,/\t/g' Dromaius_novaehollandiae_A1_like_local_blastn.csv | awk '$10' | awk '{if($9>$10) print$2,$10,$9,$1,"1","-"; else print$2,$9,$10,$1,"1","+"}' OFS="\t" 










#manually choose updated genomes, with higher N50 & coverage
species_with_ref_assembly_having_lower_contign50

#87 species
reference_genomes_with_lower_scaffoldN50_contigN50

#compare contigN50 & scaffoldN50 of continous & reference genomes: 87 genomes
compare_rep_genomes

#genome metadata of all versions
quick_compare_genomes

datasets_fetched_birds_genome_metadata_coloured
compare_genomes

eutilities_fetched_birds_genome_metadata
genome_metadata



#with annotation
Acanthisitta_chloris 

#later removed
#Picoides_pubescens - missing data - The genes is at the end of chromosome' the assembly very poor to recover
#Tyto_alba - missing data - exclude this species as it is almost(85%) intact & if it is completely inact it would be very time consuming to prove as it requires more data & even if it has a stop codon it would be in the >85% of the gene length, which is not a reliable way to establish gene loss
#Upupa epops : remove upupa epops since the loss can be due to incomplete assembly (gene is near the end of assembly)



egrep "pube|ScaffoldN50" eutilities_fetched_birds_genome_metadata | colnum.sh
grep -f <(awk '{print$NF}' apobec1_supp) eutilities_fetched_birds_genome_metadata --color=always -z | sort -k1,1 -k10,10nr | sed '1i Input_name Fetched_name Accession Last_Release_Accession Annotation Rep-status Assembly_Status ContigN50 ScaffoldN50 Coverage Submission_Date Last_Update_Date' | colnum.sh


cd ~/soft_links/Cuculus_canorus/aswin/APOBEC1/2nd_gblast/synteny
find ../../../../genome/updated_genome/curated_exon_wise_A1_like_birds/ -name "best_hits" | xargs -n1 sh -c 'grep -v "Query_length" $0' | sed '1i Query Subject Query_length Alignment_length Q_start Q_end S_start S_end E_value Bit_score Raw_score %_Query_covered_per_sub %_Query_covered_per_hsp %_ident Matches Mismatches Gaps Strand' | sort -k9,9g | awk '{print$2,$7,$8,$1,"1",$18}' | sed 's/minus/-/g' | sed 's/plus/+/g' | awk '$1="NC_071401.1"' | grep -v "S_start" | sed 's/ /\t/g' > curated_exon_wise_A1_like_birds.bed
find ../../../../genome/updated_genome/curated_exon_wise_A1_like_birds/ -name "best_hits" | xargs -n1 sh -c 'grep -v "Query_length" $0' | sed '1i Query Subject Query_length Alignment_length Q_start Q_end S_start S_end E_value Bit_score Raw_score %_Query_covered_per_sub %_Query_covered_per_hsp %_ident Matches Mismatches Gaps Strand' | sort -k9,9g | awk '{print$2,$7,$8,$1,"1",$18}' | sed 's/minus/-/g' | sed 's/plus/+/g' | awk '$1="NC_071401.1"' | grep -v "S_start" | awk '{if($2>$3) print$1,$3,$2,$4,$5,$6; else print $0}' | sed 's/ /\t/g' > curated_exon_wise_A1_like_birds.bed





































