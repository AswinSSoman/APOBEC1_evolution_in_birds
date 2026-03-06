

#Download RNA seq data

cd /media/ashutosh/disk3
time /media/ashutosh/disk3/RNA_seq/sratoolkit.3.3.0-ubuntu64/bin/prefetch  --max-size 100000000 SRR10852845







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
    split(org,a," "); last=a[length(a)]
    print ">"gene"_"acc"_"org"_"gid"_"last"_squamata"
    next} {print}' cds.fna > squamata.fa
time mafft --maxiterate 1000 --localpair --reorder --thread 8 squamata.fa > squamata.aln 


XM_014055581.1:156-8057 APOB [organism=Thamnophis sirtalis] [GeneID=106540476] [region=cds]
XM_054970013.1:53-13546 APOB [organism=Eublepharis macularius] [GeneID=129323471] [region=cds]
XM_070732802.1:59-13540 APOB [organism=Erythrolamprus reginae] [GeneID=139156769] [region=cds]

exon_wise -datasets -tid XM_006034385.3
exon_wise -datasets -tid XM_059729267.1
exon_wise -datasets -tid XM_019520673.1
exon_wise -datasets -tid XM_019537339.1





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





































