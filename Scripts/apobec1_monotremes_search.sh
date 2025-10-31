##################################################################################################################################################################################################################################################################################################
#                                                                                                                          APOBEC1 SEARCH IN MAMMALS
##################################################################################################################################################################################################################################################################################################


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mkdir -p ~/mammals_db/monotremes/
cd ~/mammals_db/monotremes/

mkdir -p Ornithorhynchus_anatinus/genome  Tachyglossus_aculeatus/genome

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download genomes of monotremes

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Run cgblast

cd ~/mammals_db/monotremes

time while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}' | sed 's/^/>/g'`
echo $sn
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`

#A1 as query 
mkdir $i1/all_APOBEC1_from_birds/
cd $i1/all_APOBEC1_from_birds/
cp /home/neo/fishes_db/A1_from_all_birds.fa .
time gblast_short ../genome/GC*.fna A1_from_all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > ../APOBEC1_all_birds.bed

#A1-like as query
mkdir ../all_APOBEC1_like_from_birds/
cd ../all_APOBEC1_like_from_birds/
cp /home/neo/bird_db1/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/cgblast/A1_like_from__all_birds.fa .
time gblast_short ../genome/GC*.fna A1_like_from__all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > ../APOBEC1_like_all_birds.bed

unset sn i1
cd ~/mammals_db/monotremes
done < <(find . -name "*.fna") 

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Cgblast using Mammals A1 as query

cd ~/mammals_db/monotremes
awk -v RS=">" '/^M_/ {print">"$0}'  ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cds_alignment/birds_mammals_filtered_8.fa | awk NF > APOBEC1_mammals_cds.fa

time while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}' | sed 's/^/>/g'`
echo $sn
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
mkdir $i1/all_APOBEC1_from_mammals/
cd $i1/all_APOBEC1_from_mammals/
cp /home/neo/mammals_db/monotremes/APOBEC1_mammals_cds.fa .
time gblast_short ../genome/GC*.fna APOBEC1_mammals_cds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > ../APOBEC1_all_mammals.bed
unset sn i1
cd ~/mammals_db/monotremes
done < <(find . -name "*.fna") 

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Show best hits
find . -name "best_hits" | xargs -n1 sh -c 'echo $0 | cut -f2,3 -d "/";sed "1d" $0 | sed "s/^/- /g"' \
	| sed '1i Group/Species/Query Query Subject Query_length Alignment_length Q_start Q_end S_start S_end E_value Bit_score Raw_score %_Query_covered_per_sub %_Query_covered_per_hsp %_ident Matches Mismatches Gaps Strand' | column -t > monotreme_best_hits
#Show Syntenic gene locations
find . -name "*.gff" | xargs -n1 sh -c 'echo $0 | cut -f2,3 -d "/"; egrep -i "nanog|aicda|mfap5|slc2a3" $0' > monotreme_syntenic_genes_gff
egrep -i "nanog|aicda|mfap5|slc2a3" monotreme_syntenic_genes_gff --color=always -z | less -RS

find . -name "Synteny" | xargs -n1 sh -c 'echo $0 | cut -f2,3 -d "/"; cat $0 | sed "s/^/- /g"' | column -t > monotreme_cgblast_synteny


mkdir -p ~/mammals_db/monotremes/Tachyglossus_aculeatus/APOBEC1_Homo_sapiens
cd ~/mammals_db/monotremes/Tachyglossus_aculeatus/APOBEC1_Homo_sapiens
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/mammals/APOBEC1_Homo_sapiens_NM_001644.5.fa .
time gblast_short ../genome/GCF_015852505.1_mTacAcu1.pri_genomic.fna APOBEC1_Homo_sapiens_NM_001644.5.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank ../genome/GCF_015852505.1_mTacAcu1.pri_genomic.gff

mkdir -p ~/mammals_db/monotremes/Ornithorhynchus_anatinus/APOBEC1_Homo_sapiens
cd ~/mammals_db/monotremes/Ornithorhynchus_anatinus/APOBEC1_Homo_sapiens
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/mammals/APOBEC1_Homo_sapiens_NM_001644.5.fa .








