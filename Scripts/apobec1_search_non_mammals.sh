############################################################################################################################################################################################################################################################################################################
#FIND APOBEC1 IN NON-MAMMALIAN VERTEBRATES
############################################################################################################################################################################################################################################################################################################

mkdir /home/neo/non_mammals_except_birds_db1
cd /home/neo/non_mammals_except_birds_db1

#============================================================================================================================================================================================================================================================================================================
#Manually choosing best genomes from each major taxonomic group (observable in NCBI legacy taxonomy page) based on:
	#genome quality: scaffold N50, contig N50, genome size
	#Annotation availability: number of genes annotated
	#Published date
	#Divergence time based on time calibrated phylogenetic tree

############################################################################################################################################################################################################################################################################################################
#Download assemblies and the QC values

#get genomes
esearch -db assembly -query "Testudines[ORGN] AND Representative [RefSeq]" | esummary > g_raw

#Meta headings
#check if meta headings are same for all entries
grep Meta g_raw | tr "<" "\n" | grep "^Stat" | awk -F "=| |>" '{if($2=="") print$0; else print$3}' | tr -d '"' | awk -v RS="Stats>" '{$1=$1}1' | less
grep Meta g_raw | head -1 | tr "<" "\n" | grep "^Stat" | awk '$3' | cut -f2 -d "=" | cut -f1 -d " " | tr -d '"' | paste -s -d " "

#Print summary
cat g_raw | xtract -pattern DocumentSummary -element AssemblyAccession,Organism,AssemblyStatus,AssemblyStatusSort,Coverage,AsmReleaseDate_GenBank,Stat | sed 's/ /_/g' | sed '1i Accession organism chr_level chr_status coverage date alt_loci_count chromosome_count contig_count contig_l50 contig_n50 non_chromosome_replicon_count replicon_count scaffold_count_all scaffold_count_placed scaffold_count_unlocalized scaffold_count_unplaced scaffold_l50 scaffold_n50 total_length ungapped_length' | awk '{gsub("/.*","",$6)}1' | sed 's/chromosome/chr/g' | sed 's/scaffold/scfld/g' | awk 'NR==1{for(x=1;x<=NF;x++)if($x!="non_chr_replicon_count" && $x!="alt_loci_count" && $x!="chr_status" && $x!="scfld_count_unlocalized" && $x!="scfld_count_placed" && $x!="") l[x]++;}{for(i=1;i<=NF;i++)if(i in l)printf (i==NF)?$i"":$i" ";printf "\n"}' | column -t > genomes

#get taxonomy lineage of downloaded assembly
esearch -db assembly -query "Testudines[ORGN] AND Representative [RefSeq]" | elink -target taxonomy | efetch -format xml > l_raw
#cat l_raw | xtract -pattern Taxon -element  ScientificName  Lineage  | sed 's/ /_/g' | sed 's/;_/ /g' | sed 's/\t_/ /g' | tr "\t" " " | awk '{for(i=1; i<=NF;i++) {if(FNR==1) {header[i]= $i;} else {if(i in data) {if(data[i]!=$i) {data[i]="different";}} else {data[i]=$i;}}} if(FNR==1) {for(i=1;i <= NF; i++) {if (i != NF) {printf"%s\t",header[i];} else {printf "%s\n", header[i];}}} else {for(i=1;i<=NF;i++) {if(data[i]=="different") {if(i!=NF) {printf "%s\t",$i;} else {printf"%s\n",$i;}}}}}' > lineages
cat l_raw | xtract -pattern Taxon -element  ScientificName  Lineage | sed 's/cellular organisms; .*Testudines;//g' | sed 's/ /_/g' | sed 's/;_/ /g' | sed 's/\t_/ /g' | sed 's/[ ]\+/ /g' | sort -k2,2 -k3,3 -k4,4  | column -t > lineages

#Compare genome & lineaege data 
while read i
do
i1=`echo $i | awk '{print$2}' | sed 's/_(.*//g'`
i2=`grep "$i1" lineages`
echo $i $i2
unset i1
done < genomes | sed 's/[ ]\+/ /g' | sort -k18,18 -k19,19 -k20,20 -k21,21 -k22,22  | column -t > compare

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#get genomes
esearch -db assembly -query "Squamata [ORGN] AND Representative [RefSeq]" | esummary > g_raw

#Meta headings
#check if meta headings are same for all entries
grep Meta g_raw | tr "<" "\n" | grep "^Stat" | awk -F "=| |>" '{if($2=="") print$0; else print$3}' | tr -d '"' | awk -v RS="Stats>" '{$1=$1}1' | less
grep Meta g_raw | head -1 | tr "<" "\n" | grep "^Stat" | awk '$3' | cut -f2 -d "=" | cut -f1 -d " " | tr -d '"' | paste -s -d " "

#Print summary
cat g_raw | xtract -pattern DocumentSummary -element AssemblyAccession,Organism,AssemblyStatus,AssemblyStatusSort,Coverage,AsmReleaseDate_GenBank,Stat | sed 's/ /_/g' | sed '1i Accession organism chr_level chr_status coverage date alt_loci_count chromosome_count contig_count contig_l50 contig_n50 non_chromosome_replicon_count replicon_count scaffold_count_all scaffold_count_placed scaffold_count_unlocalized scaffold_count_unplaced scaffold_l50 scaffold_n50 total_length ungapped_length' | awk '{gsub("/.*","",$6)}1' | sed 's/chromosome/chr/g' | sed 's/scaffold/scfld/g' | awk 'NR==1{for(x=1;x<=NF;x++)if($x!="non_chr_replicon_count" && $x!="alt_loci_count" && $x!="chr_status" && $x!="scfld_count_unlocalized" && $x!="scfld_count_placed" && $x!="") l[x]++;}{for(i=1;i<=NF;i++)if(i in l)printf (i==NF)?$i"":$i" ";printf "\n"}' | column -t > genomes

#get taxonomy lineage of downloaded assembly
esearch -db assembly -query "Squamata [ORGN] AND Representative [RefSeq]" | elink -target taxonomy | efetch -format xml > l_raw
#cat l_raw | xtract -pattern Taxon -element  ScientificName  Lineage  | sed 's/ /_/g' | sed 's/;_/ /g' | sed 's/\t_/ /g' | tr "\t" " " | awk '{for(i=1; i<=NF;i++) {if(FNR==1) {header[i]= $i;} else {if(i in data) {if(data[i]!=$i) {data[i]="different";}} else {data[i]=$i;}}} if(FNR==1) {for(i=1;i <= NF; i++) {if (i != NF) {printf"%s\t",header[i];} else {printf "%s\n", header[i];}}} else {for(i=1;i<=NF;i++) {if(data[i]=="different") {if(i!=NF) {printf "%s\t",$i;} else {printf"%s\n",$i;}}}}}' > lineages
cat l_raw | xtract -pattern Taxon -element  ScientificName  Lineage | sed 's/cellular organisms; .*Squamata;//g' | sed 's/ /_/g' | sed 's/;_/ /g' | sed 's/\t_/ /g' | sed 's/[ ]\+/ /g' | sort -k2,2 -k3,3 -k4,4  | column -t > lineages

#Compare genome & lineaege data 
while read i
do
i1=`echo $i | awk '{print$2}' | sed 's/_(.*//g'`
i2=`grep "$i1" lineages`
echo $i $i2
unset i1
done < genomes | sed 's/[ ]\+/ /g' | sort -k18,18 -k19,19 -k20,20 -k21,21 -k22,22  | column -t > compare

############################################################################################################################################################################################################################################################################################################

#Download phylogenetic tree from time tree website

############################################################################################################################################################################################################################################################################################################
#Genome blast

#Create database
cd /home/neo/non_mammals_except_birds_db1

while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}' | sed 's/^/>/g'`
echo $sn
i1=`echo $i | awk -F "/" '!($NF="")' OFS="/"`
cd $i1
g=`ls GC*.fna`
gp=`echo $g | sed 's/.fna//g'`
#create blast database
makeblastdb -in $g -out $g -dbtype nucl
#index genome
samtools faidx $g

#if annotation is present
if [[ -e genomic.gff ]]
then
#rename GFF
mv genomic.gff $gp".gff"
#convert gff to gtf
gffread -E $gp".gff" -T -o $gp".gtf"
#sort gtf
sort -k1,1 -k4,4n $gp".gtf" > $gp"_sorted.gtf"
#index gtf
/home/neo/programmes/IGV_2.12.3/igvtools index $gp"_sorted.gtf"
#convert gff to bed
sortBed -i $gp".gff" | gff2bed > $gp".gff.bed"
else :
fi
unset sn i1 g gp
cd /home/neo/non_mammals_except_birds_db1
done < <(find . -name "*.fna") 

############################################################################################################################################################################################################################################################################################################
#Run Comprehensive blast

############################################################################################################################################################################################################################################################################################################
#Query: APOBEC1_Anas_platyrhynchos

while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}' | sed 's/^/>/g'`
gr=`grep "$sn" taxonomy | awk '{print$(NF-1)}'`
echo $sn $gr
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
cd $i1
#This was mistake that the folder is named "Homo sapiens" but the query used is APOBEC1_Anas_platyrhynchos
mkdir APOBEC1_Homo_sapiens
cd APOBEC1_Homo_sapiens
g=`ls ../genome/GC*.fna`
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/APOBEC1_Anas_platyrhynchos.fa .
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/APOBEC1_Anas_platyrhynchos_flanking.fa .
#Run genome blast
if [ -n "$(ls ../genome/GC*.gff)" ]
then
gff=`ls ../genome/GC*.gff`
time gblast_short $g APOBEC1_Anas_platyrhynchos.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank $gff
else 
time gblast_short $g APOBEC1_Anas_platyrhynchos.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank
fi
unset sn gr i1 g gff
cd /home/neo/non_mammals_except_birds_db1
done < <(find . -name "*.fna") 

#summary
while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}'`
gr=`grep "$sn" taxonomy | awk '{print$(NF-1)}'`
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
#columns: exon (%_Query_covered_per_hsp, Alignment_length, E_value)
cat $i1/APOBEC1_Homo_sapiens/best_hits | grep -v "Query" | awk '{print$1"("$13","$4","$9")"}' | paste -s -d " " | sed "s/^/$sn $gr /g"
unset sn gr i1
cd /home/neo/non_mammals_except_birds_db1
done < <(find . -name "*.fna") | awk '/^>/ {print ++i". "$0;next} {print$0}' | column -t > summary_human

############################################################################################################################################################################################################################################################################################################
#Query: APOBEC1a_Anolis_carolinensis isoform 1

cd /home/neo/non_mammals_except_birds_db1

time while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}' | sed 's/^/>/g'`
gr=`grep "$sn" taxonomy | awk '{print$(NF-1)}'`
echo $sn $gr
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
cd $i1
mkdir APOBEC1a_Anolis_carolinensis
cd APOBEC1a_Anolis_carolinensis
g=`ls ../genome/GC*.fna`
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/APOBEC1a_Anolis_carolinensis.fa .
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/APOBEC1a_Anolis_carolinensis_flanking.fa .
#Run genome blast
if [ -n "$(ls ../genome/GC*.gff)" ]
then
gff=`ls ../genome/GC*.gff`
time gblast_short $g APOBEC1a_Anolis_carolinensis.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank $gff
else 
time gblast_short $g APOBEC1a_Anolis_carolinensis.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank
fi
unset sn gr i1 g gff
cd /home/neo/non_mammals_except_birds_db1
done < <(find . -name "*.fna") 

#summary
while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}'`
gr=`grep "$sn" taxonomy | awk '{print$(NF-1)}'`
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
#columns: exon (%_Query_covered_per_hsp, Alignment_length, E_value)
cat $i1/APOBEC1a_Anolis_carolinensis/best_hits | grep -v "Query" | awk '{print$1"("$13","$4","$9")"}' | paste -s -d " " | sed "s/^/$sn $gr /g"
unset sn gr i1
cd /home/neo/non_mammals_except_birds_db1
done < <(find . -name "*.fna") | awk '/^>/ {print ++i". "$0;next} {print$0}' | column -t > summary_lizard1

############################################################################################################################################################################################################################################################################################################
#Query: APOBEC1a_Anolis_carolinensis isoform 2

time while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}' | sed 's/^/>/g'`
gr=`grep "$sn" taxonomy | awk '{print$(NF-1)}'`
echo $sn $gr
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
cd $i1
mkdir APOBEC1_Anolis_carolinensis
cd APOBEC1_Anolis_carolinensis
g=`ls ../genome/GC*.fna`
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/APOBEC1_Anolis_carolinensis.fa .
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/APOBEC1_Anolis_carolinensis_flanking.fa .
#Run genome blast
if [ -n "$(ls ../genome/GC*.gff)" ]
then
gff=`ls ../genome/GC*.gff`
time gblast_short $g APOBEC1_Anolis_carolinensis.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank $gff
else 
time gblast_short $g APOBEC1_Anolis_carolinensis.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank
fi
unset sn gr i1 g gff
cd /home/neo/non_mammals_except_birds_db1
done < <(find . -name "*.fna") 

#summary
while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}'`
gr=`grep "$sn" taxonomy | awk '{print$(NF-1)}'`
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
#columns: exon (%_Query_covered_per_hsp, Alignment_length, E_value)
cat $i1/APOBEC1_Anolis_carolinensis/best_hits | grep -v "Query" | awk '{print$1"("$13","$4","$9")"}' | paste -s -d " " | sed "s/^/$sn $gr /g"
unset sn gr i1
cd /home/neo/non_mammals_except_birds_db1
done < <(find . -name "*.fna") | awk '/^>/ {print ++i". "$0;next} {print$0}' | column -t > summary_lizard2

####################################################################################################################################################################################################################################################################################################################
#Query: A1-like ostrich

cd /home/neo/non_mammals_except_birds_db1

#cgblast (23m26.232s)
time while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}' | sed 's/^/>/g'`
gr=`grep "$sn" orders | awk '{print$NF}'`
echo $sn $gr
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
cd $i1
mkdir APOBEC1_like_Struthio_camelus
cd APOBEC1_like_Struthio_camelus
g=`ls ../genome/GC*.fna`
cp ~/soft_links/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/APOBEC1_like_Struthio_camelus_1_exonwise.fa .
cp ~/soft_links/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/APOBEC1_like_Struthio_camelus_1_exonwise_flanking.fa .
#Run genome blast
if [ -n "$(ls ../genome/GC*.gff)" ]
then
gff=`ls ../genome/GC*.gff`
time gblast_short $g APOBEC1_like_Struthio_camelus_1_exonwise.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank $gff
else 
time gblast_short $g APOBEC1_like_Struthio_camelus_1_exonwise.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank
fi
unset sn gr i1 g gff
cd /home/neo/non_mammals_except_birds_db1
done < <(find . -name "*.fna") 

#summary
while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}'`
gr=`grep "$sn" orders | awk '{print$NF}'`
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
#columns: exon (%_Query_covered_per_hsp, Alignment_length, E_value)
cat $i1/APOBEC1_like_Struthio_camelus/best_hits | grep -v "Query" | awk '{print$1"("$13","$4","$9")"}' | paste -s -d " " | sed "s/^/$sn $gr /g"
unset sn gr i1
cd /home/neo/non_mammals_except_birds_db1
done < <(find . -name "*.fna") | awk '/^>/ {print ++i". "$0;next} {print$0}' | column -t > summary_A1_like_ostrich

#APOBEC1 like search gave unique hit in all testudines & only alligator but not in other crocodylia

############################################################################################################################################################################################################################################################################################################
#Query: all A1-like exon-wise CDS from all birds

#Using all queries available is much better than a single query, hence run cgblast again from A1 & A1-like

cd ~/non_mammals_except_birds_db1

time while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}' | sed 's/^/>/g'`
echo $sn
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
mkdir $i1/all_APOBEC1_like_from_birds/
cd $i1/all_APOBEC1_like_from_birds/
cp /home/neo/bird_db1/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/cgblast/A1_like_from__all_birds.fa .
time gblast_short ../genome/GC*.fna A1_like_from__all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_like_all_birds.bed
unset sn i1
cd ~/non_mammals_except_birds_db1
done < <(find . -name "*.fna") 

############################################################################################################################################################################################################################################################################################################
#Query: all validated A1 from all birds

cd ~/non_mammals_except_birds_db1

#273m34.064s
time while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}' | sed 's/^/>/g'`
echo $sn
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
mkdir $i1/all_APOBEC1_from_birds/
cd $i1/all_APOBEC1_from_birds/
cp /home/neo/fishes_db/A1_from_all_birds.fa .
time gblast_short ../genome/GC*.fna A1_from_all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_all_birds.bed
unset sn i1
cd ~/non_mammals_except_birds_db1
done < <(find . -name "*.fna") 

############################################################################################################################################################################################################################################################################################################
#Query: A1 CDS from all mammals

cd ~/non_mammals_except_birds_db1

#345m16.914s
time while read i
do
sn=`echo $i | awk -F "/" '{print$(NF-2)}' | sed 's/^/>/g'`
echo $sn
i1=`echo $i | awk -F "/" '{$NF=$(NF-1)=""; print$0}' OFS="/" | sed 's!//!/!g'`
mkdir $i1/all_APOBEC1_from_mammals/
cd $i1/all_APOBEC1_from_mammals/
cp /home/neo/mammals_db/monotremes/APOBEC1_mammals_cds.fa .
time gblast_short ../genome/GC*.fna APOBEC1_mammals_cds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_all_mammals.bed
unset sn i1
cd ~/non_mammals_except_birds_db1
done < <(find . -name "*.fna") 

############################################################################################################################################################################################################################################################################################################
#Check in more non-mammals

#Some groups had only single representative species & these species shows fragmented synteny or very different synteny w.r.t others
#This need to be confirmed by checking more species from the respective group

#Download genomes
cd ~/non_mammals_except_birds_db1/amniota/Squamata/
for i in `readlink -f Gekkos/* Laterata/* Iguania/* Scinciformata/* Snakes/* | egrep -v "Eublepharis_macularius|Aspidoscelis_tigris_stejnegeri|Podarcis_raffonei|Anolis_carolinensis|Pogona_vitticeps|Sceloporus_undulatus|Hemicordylus_capensis|Ahaetulla_prasina|Naja_naja|Protobothrops_mucrosquamatus|Python_bivittatus"`
do
echo ">"$i
cd $i
#A1 birds
mkdir all_APOBEC1_from_birds/
cd all_APOBEC1_from_birds/
cp /home/neo/fishes_db/A1_from_all_birds.fa .
time gblast_short ../genome/GC*.fna A1_from_all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank ../genome/GC*.gff
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_all_birds.bed
#A1 mammals
cd $i
mkdir all_APOBEC1_from_mammals/
cd all_APOBEC1_from_mammals/
cp /home/neo/mammals_db/monotremes/APOBEC1_mammals_cds.fa .
time gblast_short ../genome/GC*.fna APOBEC1_mammals_cds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank ../genome/GC*.gff
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_all_mammals.bed
#A1 like birds
cd $i
mkdir all_APOBEC1_like_from_birds/
cd all_APOBEC1_like_from_birds/
cp /home/neo/bird_db1/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/cgblast/A1_like_from__all_birds.fa .
time gblast_short ../genome/GC*.fna A1_like_from__all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank ../genome/GC*.gff
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_like_all_birds.bed
cd ~/non_mammals_except_birds_db1/amniota/Squamata/
done

#Show summary of for best hits

#For amphibia
cd ~/non_mammals_except_birds_db1/amphibia

#Show best hits
find . -name "best_hits" | xargs -n1 sh -c 'echo $0 | sed "s/.best_hits$//g" ; sed "1d" $0 | sed "s/^/- /g"' \
	| sed '1i Group/Species/Query Query Subject Query_length Alignment_length Q_start Q_end S_start S_end E_value Bit_score Raw_score %_Query_covered_per_sub %_Query_covered_per_hsp %_ident Matches Mismatches Gaps Strand' | column -t > amphibia_best_hits
#Show Syntenic gene locations
find . -name "*.gff" | xargs -n1 sh -c 'echo $0 | cut -f2,3 -d "/"; egrep -i "nanog|aicda|mfap5|slc2a3" $0' > amphibia_syntenic_genes_gff
egrep -i "nanog|aicda|mfap5|slc2a3" amphibia_syntenic_genes_gff --color=always -z | less -RS

#For amniota
cd ~/non_mammals_except_birds_db1/amniota

#Show best hits
find . -name "best_hits" | xargs -n1 sh -c 'echo $0 | sed "s/.best_hits$//g" ; sed "1d" $0 | sed "s/^/- /g"' \
	| sed '1i Group/Species/Query Query Subject Query_length Alignment_length Q_start Q_end S_start S_end E_value Bit_score Raw_score %_Query_covered_per_sub %_Query_covered_per_hsp %_ident Matches Mismatches Gaps Strand' | column -t > amniota_best_hits
#Show Syntenic gene locations
find . -name "*.gff" | xargs -n1 sh -c 'echo $0 | cut -f2,3,4 -d "/"; egrep -i "apobec1|a2ml|nanog|aicda|rimklb|mfap5|slc2a3|foxj2|dnm2|phc1|Ovostatin|m6pr|timm29|yipf2|carm1" $0' | awk '$3=="gene" || $2==""' > amniota_syntenic_genes_gff
egrep -i "nanog|aicda|mfap5|slc2a3" amniota_syntenic_genes_gff --color=always -z | less -RS
find . -name "Synteny" | xargs -n1 sh -c 'echo $0 | cut -f2,3 -d "/"; cat $0 | sed "s/^/- /g"' | column -t > amniota_cgblast_synteny

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Add one more amphibian genome

cd ~/non_mammals_except_birds_db1/amphibia/Salamanders/Ambystoma_mexicanum/

#A1 birds
mkdir all_APOBEC1_from_birds/
cd all_APOBEC1_from_birds/
cp /home/neo/fishes_db/A1_from_all_birds.fa .
time gblast_short ../genome/GC*.fna A1_from_all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_all_birds.bed

#A1 mammals
cd ~/non_mammals_except_birds_db1/amphibia/Salamanders/Ambystoma_mexicanum/
mkdir all_APOBEC1_from_mammals/
cd all_APOBEC1_from_mammals/
cp /home/neo/mammals_db/monotremes/APOBEC1_mammals_cds.fa .
time gblast_short ../genome/GC*.fna APOBEC1_mammals_cds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_all_mammals.bed

#A1 like birds
cd ~/non_mammals_except_birds_db1/amphibia/Salamanders/Ambystoma_mexicanum/
mkdir all_APOBEC1_like_from_birds/
cd all_APOBEC1_like_from_birds/
cp /home/neo/bird_db1/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/cgblast/A1_like_from__all_birds.fa .
time gblast_short ../genome/GC*.fna A1_like_from__all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_like_all_birds.bed

#A1 from publications
cd ~/non_mammals_except_birds_db1/amphibia/Salamanders/Ambystoma_mexicanum/
mkdir A1_from_publication/
cd A1_from_publication/
cp ~/bird_db1/aswin/APOBEC1/outgroup/Diversification_of_AID_APOBEC_like_deaminases_in_metazoa_2018_March/APOBEC1.aa A1_from_publication.aa
time tblastn -task tblastn -query A1_from_publication.aa -db ../genome/GC*.fna -evalue 1e-5 -num_threads 4 \
-outfmt "6 qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps qseq sseq" \
| sed '1i \ Query\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tQuery_sequence\tSubject_sequence\n' > A1_from_publication_tblastn.out
time tblastn -task tblastn -query A1_from_publication.aa -db ../genome/GC*.fna -evalue 0.005 -num_threads 4 -outfmt 3 -line_length 180 > A1_from_publication_tblastn.outfmt3

mkdir ~/non_mammals_except_birds_db1/amphibia/Salamanders/Ambystoma_mexicanum/Query_from_publications

time tblastn -task tblastn -query APOBE1_Ambystoma_mexicanum.fa -db ../genome/GC*.fna -evalue 1e-5 -num_threads 4 \
-outfmt "6 qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps qseq sseq" \
| sed '1i \ Query\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tQuery_sequence\tSubject_sequence\n' > tblastn.out

#Perform online blast in ncbi

cd ~/non_mammals_except_birds_db1/amphibia/ncbi_blast

#A1 from birds, mammals & non-mammals against Caudata nr/nt database
#against Ambystoma_mexicanum SRA
awk '$8=="GENOMIC" && $7=="WGS" && $5=="ILLUMINA"' sra_accessions | sort -k10,10nr | less

############################################################################################################################################################################################################################################################################################################
#Merge all summaries

cd ~/non_mammals_except_birds_db1
for i in `ls summary_*`
do
in=`echo $i | sed 's/summary_//g'`
cat $i | sed "s/^/$in /g"
done | sed 's/[ ]\+/ /g' | sort -k2,2 | column -t > merged_A1_summary

#Manually add "" if the hit is significant i.e., if the alignment length >90bp/100bp and with very low e-value like e^04 

#To view in custom data in NCBI genome viewer

sed 's/[ ]\+/ /g' merged_A1_summary | awk '$NF=="✅"' | grep "A1_like_ostrich" | sort -k1,1 -k3,3 -k4,4 -V| column -t 
sed 's/[ ]\+/ /g' merged_A1_summary | awk '$NF=="✅"' | grep -v "A1_like_ostrich" | sort -k2,2 -k3,3 -k4,4 -V| column -t

for i in `find amniota/Crocodylia/Alligator_sinensis/ -name "test.out.bed"`; do q=`echo $i | awk -F "/" '{print$(NF-2)}'`; cat $i | awk -v n="$q" '{print$1,$2,$3,$4"_"n}' OFS="\t" | grep -v "like" | sed 's/APOBEC1_//g'; done | sed 's/Anolis_carolinensis/lizard/g' | sed 's/APOBEC1a_//g' | sed 's/Homo_sapiens/human/g'
for i in `find amniota/Crocodylia/Alligator_sinensis/ -name "test.out.bed"`; do q=`echo $i | awk -F "/" '{print$(NF-2)}'`; cat $i | awk -v n="$q" '{print$1,$2,$3,$4"_"n}' OFS="\t" | grep "like" | sed 's/APOBEC1_//g'; done | sed 's/Struthio_camelus/ostrich/g' | sed 's/_like//g' 


####################################################################################################################################################################################################################################################################################################################
#DRAFT
####################################################################################################################################################################################################################################################################################################################

#Check latest version of genome
find . -name "*.fna" | awk -F "/" '{print$(NF-2),$NF}' > genomes
awk '{print$1}' genomes | xargs -n1 sh -c 'esearch -db genome -query "$0 [ORGN]" | esummary | xtract -pattern DocumentSummary -element Organism_Name Assembly_Accession' | sed 's/ /_/g' > updated_genome

while read i
do
i1=`echo $i | awk '{print$1}'`
i2=`echo $i | awk '{print$2}'`
j2=`echo $i2 | sed 's/GC[AF]_//g'`
i3=`grep "$i1" genomes | awk '{print$2}' | cut -f1,2 -d "_"`
if [[ "$i3" == "" ]]; then i3="-"; else :; fi
j3=`echo $i3 | sed 's/GC[AF]_//g'`
echo $i1 $i2 $j2 $i3 $j3
done < updated_genome | awk '{if($3==$5) print$0, "same"; else print $0, "diff"}' | column -t > compare_genomes

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Only 2 species had updated genome : Anolis_carolinensis & Pleurodeles_waltl (other 2 species (Protobothrops_mucrosquamatus, Python_bivittatus) number difference is just between GCF & GCA)

# Anolis_carolinensis

cd ~/non_mammals_except_birds_db1/amniota/Squamata/Iguania/Anolis_carolinensis
mkdir old_results
mv APOBEC1_Homo_sapiens APOBEC1a_Anolis_carolinensis APOBEC1_Anolis_carolinensis APOBEC1_like_Struthio_camelus all_APOBEC1_like_from_birds all_APOBEC1_from_birds old_results/
mkdir APOBEC1_Homo_sapiens APOBEC1a_Anolis_carolinensis APOBEC1_Anolis_carolinensis APOBEC1_like_Struthio_camelus all_APOBEC1_like_from_birds all_APOBEC1_from_birds
#
cd APOBEC1a_Anolis_carolinensis
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1a_Anolis_carolinensis.fa .
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1a_Anolis_carolinensis_flanking.fa .
time gblast_short ../genome/GC*.fna APOBEC1a_Anolis_carolinensis.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank ../genome/GC*.gff
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1a_Anolis_carolinensis.bed
#
cd ../APOBEC1_Anolis_carolinensis
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1_Anolis_carolinensis.fa .
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1_Anolis_carolinensis_flanking.fa .
time gblast_short ../genome/GC*.fna APOBEC1_Anolis_carolinensis.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank ../genome/GC*.gff
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_Anolis_carolinensis.bed
#
cd ../all_APOBEC1_like_from_birds
cp /home/neo/bird_db1/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/cgblast/A1_like_from__all_birds.fa .
time gblast_short ../genome/GC*.fna A1_like_from__all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_like_all_birds.bed
#
cd ../all_APOBEC1_from_birds
cp /home/neo/fishes_db/A1_from_all_birds.fa .
time gblast_short ../genome/GC*.fna A1_from_all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_all_birds.bed

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Pleurodeles_waltl

cd ~/non_mammals_except_birds_db1/amphibia/Salamanders/Pleurodeles_waltl
mkdir old_results
mv APOBEC1_Homo_sapiens APOBEC1a_Anolis_carolinensis APOBEC1_Anolis_carolinensis APOBEC1_like_Struthio_camelus all_APOBEC1_like_from_birds all_APOBEC1_from_birds old_results/
mkdir APOBEC1_Homo_sapiens APOBEC1a_Anolis_carolinensis APOBEC1_Anolis_carolinensis APOBEC1_like_Struthio_camelus all_APOBEC1_like_from_birds all_APOBEC1_from_birds
#
cd APOBEC1a_Anolis_carolinensis
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1a_Anolis_carolinensis.fa .
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1a_Anolis_carolinensis_flanking.fa .
time gblast_short ../genome/GC*.fna APOBEC1a_Anolis_carolinensis.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank ../genome/GC*.gff
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1a_Anolis_carolinensis.bed
#
cd ../APOBEC1_Anolis_carolinensis
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1_Anolis_carolinensis.fa .
cp /home/neo/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1_Anolis_carolinensis_flanking.fa .
time gblast_short ../genome/GC*.fna APOBEC1_Anolis_carolinensis.fa -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank ../genome/GC*.gff
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_Anolis_carolinensis.bed
#
cd ../all_APOBEC1_like_from_birds
cp /home/neo/bird_db1/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/cgblast/A1_like_from__all_birds.fa .
time gblast_short ../genome/GC*.fna A1_like_from__all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_like_all_birds.bed
#
cd ../all_APOBEC1_from_birds
cp /home/neo/fishes_db/A1_from_all_birds.fa .
time gblast_short ../genome/GC*.fna A1_from_all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_all_birds.bed

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Align A1 anolis carolensis with A1 birds
mafft --auto --reorder --quiet \
<(cat <(ls ~/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1_Anolis_carolinensis.fa ~/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1a_Anolis_carolinensis.fa | xargs -n1 sh -c 'echo $0 | cut -f9 -d "/" | sed "s/^/>/g"; grep -v ">" $0 | paste -s -d ""') \
<(ls ~/bird_db1/aswin/APOBEC1/validated_sequences/*.fa | xargs -n1 sh -c 'echo $0 | cut -f8 -d "/" | sed "s/^/>/g"; grep -v ">" $0 | paste -s -d ""')) > APOBEC1_Anolis_carolinensis_with_A1_birds_cds.aln

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Main types of APOBEC1 & A1-like queries

cp ~/fishes_db/A1_from_all_birds.fa ~/bird_db1/aswin/APOBEC1/outgroup/
cp /home/neo/bird_db1/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/cgblast/A1_like_from__all_birds.fa ~/bird_db1/aswin/APOBEC1/outgroup/A1_like_from_all_birds.fa
cp ~/mammals_db/monotremes/APOBEC1_mammals_cds.fa ~/bird_db1/aswin/APOBEC1/outgroup/A1_from_all_mammals.fa
ls /home/neo/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1_Anolis_carolinensis.fa
ls /home/neo/bird_db1/aswin/APOBEC1/outgroup/non_mammals/APOBEC1a_Anolis_carolinensis.fa


