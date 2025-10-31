##################################################################################################################################################################################################################################################################################################
#                                                                                                                          APOBEC1 SEARCH IN FISHES
##################################################################################################################################################################################################################################################################################################

#Mainly 2 types of Bony fishes:
	#- Ray-finned 
	#- Lobe-finned

#Lobe-finned is very small group from which tetrpods emerged
#There are only 3 genomes available in this group in which 1 is coelacanths & other 2 are lungfishes

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mkdir -p ~/fishes_db/lobe_finned_fishes
mkdir ~/fishes_db/ray_finned_fishes

cd ~/fishes_db/lobe_finned_fishes
mkdir -p Protopterus_annectens/genome Neoceratodus_forsteri/genome Latimeria_chalumnae/genome

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download genomes of Lobe finned fishes

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get query

cd ~/fishes_db
for i in `ls ~/bird_db1/aswin/APOBEC1/validated_sequences/*.fa`
do
in=`echo $i | awk -F "/" '{print$NF}' | sed 's/APOBEC1_//g' | sed 's/\.fa//g'`
cat $i | sed "/>/ s/$/_$in/g"
unset in
done > A1_from_all_birds.fa

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run cgblast

cd ~/fishes_db/lobe_finned_fishes

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
cd ~/fishes_db/lobe_finned_fishes
done < <(find . -name "*.fna") 

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Show best hits
find . -name "best_hits" | xargs -n1 sh -c 'echo $0 | cut -f2,3 -d "/";sed "1d" $0 | sed "s/^/- /g"' \
	| sed '1i Group/Species/Query Query Subject Query_length Alignment_length Q_start Q_end S_start S_end E_value Bit_score Raw_score %_Query_covered_per_sub %_Query_covered_per_hsp %_ident Matches Mismatches Gaps Strand' | column -t > fishes_best_hits
#Show Syntenic gene locations
find . -name "*.gff" | xargs -n1 sh -c 'echo $0 | cut -f2,3 -d "/"; egrep -i "nanog|aicda|mfap5|slc2a3" $0' > fishes_syntenic_genes_gff
egrep -i "nanog|aicda|mfap5|slc2a3" fishes_syntenic_genes_gff --color=always -z | less -RS

find . -name "Synteny" | xargs -n1 sh -c 'echo $0 | cut -f2,3 -d "/"; cat $0 | sed "s/^/- /g"' | column -t > fishes_cgblast_synteny


