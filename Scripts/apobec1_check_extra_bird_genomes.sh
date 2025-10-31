##################################################################################################################################################################################################################################################################################################
#															cgblast on extra birds 
##################################################################################################################################################################################################################################################################################################

#Some birds not included in the initial representatives, but are reported in previous literature

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#PAPER: A Long-Running Arms Race between APOBEC1 Genes and Retroviruses in Tetrapods: 2023 Jan

#This paper also reported 5 birds to have mutliple copies of apobec1, but this may be the A1-like gene next to A1, atleast this is true in Corvus cornix
	#Corvus cornix            - Passeriformes > Corvidae  - A1-like & A1 present like other birds no more copies
	#Balearica regulorum      - Gruiformes
	#Pelecanus crispus        - Pelecaniformes
	#Spheniscus magellanicus  - Sphenisciformes           -
	#Antrostomus carolinensis - Caprimulgiformes

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Balearica_regulorum

mkdir ~/extra_bird_genomes/Balearica_regulorum/APOBEC1
cd ~/extra_bird_genomes/Balearica_regulorum/APOBEC1
cp ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Atlantisia_rogersi.fa .
gblast_short ../genome/GCF_000709895.1_ASM70989v1_genomic.fna APOBEC1_Atlantisia_rogersi.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank ../genome/GCF_000709895.1_ASM70989v1_genomic.gff

#Run comprehensive genome blast
mkdir ~/extra_bird_genomes/Balearica_regulorum/APOBEC1_like/
cd ~/extra_bird_genomes/Balearica_regulorum/APOBEC1_like/
cp ~/soft_links/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/APOBEC1_like_Struthio_camelus_1_exonwise.fa .
gblast_short ../genome/GCF_000709895.1_ASM70989v1_genomic.fna APOBEC1_like_Struthio_camelus_1_exonwise.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank ../genome/GCF_000709895.1_ASM70989v1_genomic.gff

mkdir ~/extra_bird_genomes/Balearica_regulorum/all_APOBEC1_like_from_birds
cd ~/extra_bird_genomes/Balearica_regulorum/all_APOBEC1_like_from_birds
cp ~/extra_bird_genomes/Motacilla_alba/all_APOBEC1_like_from_birds/A1_like_from__all_birds.fa .
gblast_short ../genome/GCF_000709895.1_ASM70989v1_genomic.fna A1_like_from__all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank ../genome/GCF_000709895.1_ASM70989v1_genomic.gff

#IGV snapshot
cd ~/extra_bird_genomes/Balearica_regulorum/
cp all_APOBEC1_like_from_birds/blast_qc/test.out.bed APOBEC1_like_all_birds.bed
cp APOBEC1_like/blast_qc/test.out.bed APOBEC1_like.bed
cp APOBEC1/blast_qc/test.out.bed APOBEC1.bed
my_igv genome/GCF_000709895.1_ASM70989v1_genomic.fna genome/GCF_000709895.1_ASM70989v1_genomic.gff,APOBEC1.bed,APOBEC1_like.bed,APOBEC1_like_all_birds.bed

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Spheniscus_magellanicus

mkdir ~/extra_bird_genomes/Spheniscus_magellanicus/APOBEC1/
cd ~/extra_bird_genomes/Spheniscus_magellanicus/APOBEC1/
cp ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Aptenodytes_forsteri.fa .
gblast_short ../genome/GCA_010076225.1_BGI_Smag.V1_genomic.fna APOBEC1_Aptenodytes_forsteri.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank ../genome/GCA_010076225.1_BGI_Smag.V1_genomic.gff

#Run comprehensive genome blast
mkdir ~/extra_bird_genomes/Spheniscus_magellanicus/APOBEC1_like/
cd ~/extra_bird_genomes/Spheniscus_magellanicus/APOBEC1_like/
cp ~/soft_links/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/APOBEC1_like_Struthio_camelus_1_exonwise.fa .
gblast_short ../genome/GCA_010076225.1_BGI_Smag.V1_genomic.fna APOBEC1_like_Struthio_camelus_1_exonwise.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank ../genome/GCA_010076225.1_BGI_Smag.V1_genomic.gff

mkdir ~/extra_bird_genomes/Spheniscus_magellanicus/all_APOBEC1_like_from_birds
cd ~/extra_bird_genomes/Spheniscus_magellanicus/all_APOBEC1_like_from_birds
cp ~/extra_bird_genomes/Motacilla_alba/all_APOBEC1_like_from_birds/A1_like_from__all_birds.fa .
gblast_short ../genome/GCA_010076225.1_BGI_Smag.V1_genomic.fna A1_like_from__all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank ../genome/GCA_010076225.1_BGI_Smag.V1_genomic.gff

#IGV snapshot
cd ~/extra_bird_genomes/Spheniscus_magellanicus
cp all_APOBEC1_like_from_birds/blast_qc/test.out.bed APOBEC1_like_all_birds.bed
cp APOBEC1_like/blast_qc/test.out.bed APOBEC1_like.bed
cp APOBEC1/blast_qc/test.out.bed APOBEC1.bed
my_igv genome/GCA_010076225.1_BGI_Smag.V1_genomic.fna genome/GCA_010076225.1_BGI_Smag.V1_genomic.gff,APOBEC1.bed,APOBEC1_like.bed,APOBEC1_like_all_birds.bed

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Antrostomus_carolinensis

mkdir ~/extra_bird_genomes/Antrostomus_carolinensis/APOBEC1
cd ~/extra_bird_genomes/Antrostomus_carolinensis/APOBEC1
cp ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Caprimulgus_europaeus.fa .
gblast_short ../genome/GCF_000700745.2_ASM70074v2_genomic.fna APOBEC1_Caprimulgus_europaeus.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank ../genome/GCF_000700745.2_ASM70074v2_genomic.gff

#Run comprehensive genome blast
mkdir ~/extra_bird_genomes/Antrostomus_carolinensis/APOBEC1_like/
cd ~/extra_bird_genomes/Antrostomus_carolinensis/APOBEC1_like/
cp ~/soft_links/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/APOBEC1_like_Struthio_camelus_1_exonwise.fa .
gblast_short ../genome/GCF_000700745.2_ASM70074v2_genomic.fna APOBEC1_like_Struthio_camelus_1_exonwise.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank ../genome/GCF_000700745.2_ASM70074v2_genomic.gff

mkdir ~/extra_bird_genomes/Spheniscus_magellanicus/all_APOBEC1_like_from_birds
cd ~/extra_bird_genomes/Antrostomus_carolinensis/all_APOBEC1_like_from_birds
cp ~/extra_bird_genomes/Motacilla_alba/all_APOBEC1_like_from_birds/A1_like_from__all_birds.fa .
gblast_short ../genome/GCF_000700745.2_ASM70074v2_genomic.fna A1_like_from__all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank ../genome/GCF_000700745.2_ASM70074v2_genomic.gff

#IGV snapshot
cd ~/extra_bird_genomes/Antrostomus_carolinensis
cp all_APOBEC1_like_from_birds/blast_qc/test.out.bed APOBEC1_like_all_birds.bed
cp APOBEC1_like/blast_qc/test.out.bed APOBEC1_like.bed
cp APOBEC1/blast_qc/test.out.bed APOBEC1.bed
my_igv genome/GCF_000700745.2_ASM70074v2_genomic.fna genome/GCF_000700745.2_ASM70074v2_genomic.gff,APOBEC1.bed,APOBEC1_like.bed,APOBEC1_like_all_birds.bed

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Taeniopygia_guttata

cd ~/soft_links/Taeniopygia_guttata/aswin/APOBEC1/all_APOBEC1_like_from_birds
cp ~/extra_bird_genomes/Motacilla_alba/all_APOBEC1_like_from_birds/A1_like_from__all_birds.fa .
time gblast_short ../../../genome/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna A1_like_from__all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank ../../../genome/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.gff

cd ~/soft_links/Taeniopygia_guttata/aswin/APOBEC1/all_APOBEC1_like_from_birds
awk '!($NF="")' OFS="\t" ../test.out.bed | sed 's/\t$//g' > APOBEC1.bed
awk '!($NF="")' OFS="\t" ~/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/cgblast/APOBEC1_like_Struthio_camelus_1_as_query/Taeniopygia_guttata/test.out.bed | sed 's/\t$//g' > APOBEC1_like.bed
awk '!($NF="")' OFS="\t" blast_qc/test.out.bed | sed 's/\t$//g' > APOBEC1_like_all_birds.bed

#Exon_4 is still lacking
cd ~/soft_links/Taeniopygia_guttata/aswin/APOBEC1/all_APOBEC1_like_from_birds/exon_4
#extract region between exon_3 of A1-like and exon_1 of A1 from ncbi
#Search in this narrow region
for i in `ls ~/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/annotated_A1_like_sequences/A1_like_*.fa | grep -v "cds"`
do
in=`echo $i | awk -F "/" '{print$NF}' | sed 's/\.fa//g'`
cat $i | sed "/^>/ s/$/&_$in/g"
done > A1_like_from__all_birds.fa
sed 's/_A1_like.*//g' A1_like_from__all_birds.fa -i
time gblast_short NC_044998.1_28396845..28406198.fa A1_like_from__all_birds.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#This paper also tells these 2 birds have 2 Z domains in their A1 gene
	#Hirundo rustica - Passeriformes > Sylvioidea
	#Motacilla alba  - Passeriformes > Passeroidea > Motacillidae

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#A1 & A1-like: single or seperate genes?

#In 3 species - A1 & A1-like is annaotated as a single gene, but based on observations it suggests that these are likely 2 seperate genes
	#1 Taeniopygia_guttata:
		#- No gaps in assembly
		#- Single isoform
			#- Total 9 exons
			#- Exons 1 to 3 is similar to A1-like
			#- Exons 5 to 9 is similar to A1
			#- Exon 4 doesn't show similarity to A1/A1-like based on blastn search
		#- Intron spanning reads doesn't cover whole gene, there is a gap between exon 3 and 5
		#- Both RNA-seq exon coverage reads and intron spanning reads show that there is an exon after exon 3 and before exon 4
			#- The exon 4 annotation doesn't show any exon/intron read support suggesting that the current exon 4 might be a misannotation
			#- The more accurate exon 4 annotation is ~7kb upstream to current exon 4
	#2 Anser_cygnoides:
		#- Have assembly gaps within this A1 + A1-like annotation
		#- Intron spanning reads are spanning from start to end of gene
	#3 Motacilla alba
		#- No gaps in assembly
		#- Single isoform
			#- Total 7 exons
			#- Exons 1 to 4 is similar to A1-like
			#- Exons 5 to 7 is similar to A1
		#- Intron spanning reads doesn't cover whole gene, there is a gap between exon 4 and 5
		#- Based on exon/intron read support the current annotated exons are likely accurate, but might be missing an extra exon between exon 4 and 5.
			#- This exon is exon 1 of A1 and is detected when exon 1 from all validated birds (with flanking sequence) between exon 4 & 5 
		#- No splice site after exon 4, but a STOP codon is present suggesting this might be the end of the gene 

#
##################################################################################################################################################################################################################################################################################################
#Birds who lost A1 independently, prefferrebly must have a strong evidence to support the loss such as:
	#1. More birds in the same brach having loss even if the mutations are not shared. 
 	#2. Multi platform evidence: DNA (genome, SRA, Long read), RNA

cd ~/bird_db1/aswin/APOBEC1

#Highlight all the birds with A1 loss & identify remotely located ones
grep -if total_validated_losses_Upupa_epops_removed <(cat main_figures/tree/species_in_order_of_loss_events_tree_expanded_upupa_epops_removed | xargs -n1 sh -c 'grep "$0" ../taxonomy/orders_all_birds') -z

#Updated genomes:
 #Columba livia have intact A1 but had an updated genome which shows 2 copies of SLC2A3, NANOG, NANOG, APOBEC1, AICDA, MFAP5, RIMKLB in 2 seperat unplaced scaffolds.
  #links to view:
   # https://www.ncbi.nlm.nih.gov/projects/sviewer/?id=2724854094&tkey=q4HJrKSjr6ahqKWqjp-6ppzUxsnh2tri9snW_N9I_x17a2pCkLlU0RgKHvgEzHf8BJg&assm_context=GCF_036013475.1&mk=98922:104471|APOBEC1|blue,81035:88323|APOBEC1|red,110491:117573|AICDA|green,60151:63581|NANOG|993366,47776:54675|NANOG|993366,13276:25331|SLC2A3|FF00FF,155401:217425|RIMKLB|800000,128401:139706|MFAP5|003300&v=1501:225488&c=00ff00&select=null&slim=0
   # https://www.ncbi.nlm.nih.gov/projects/sviewer/?id=2724854093&tkey=2PK639fQ3NXS29bZ_ezJ1e-ntbqSqamRhbqlj6w7jG4IWRmxrhJaetrEsyapEtoiqUY&assm_context=GCF_036013475.1&mk=751771:758782|AICDA|FF00FF,764729:770337|APOBEC1|FF0000,805638:810125|NANOG-like|blue,814441:821431|NANOG-like|red,633609:713548|RIMKLB|green,843559:875815|SLC2A3|003300&v=628642:889576&c=ff6600&select=gi|2724854093-000badea-000bae37-0200-cc62bce3-5d97c85c;&slim=0

#extra_birds_closer_to_Pterocles_gutturalis
#Identify the bigger clade in which Pterocles_gutturalis is a part of such that atleast 5-6 other bird genomes are available with atleast scaffold level & good scaffold N50 preferrably with an annotation
#Choose 5-6 genomes including the bird with loss with good scaffold N50 & Run ncbi nr blastn with following inputs & parameters:
	#1. Query: validated APOBEC1 gene sequence with flanking region of all birds in this clade
 	#2. Subject: nr/nt database specifying 5-6 birds 
  	#3. Parameter: Somewhat similar sequence
#Download output in xml format

mv HHKC2G1B016-Alignment.xml A1_flanking_Columba_livia_Mesitornis_unicolor_against_Columbimorphea_nr_blast.xml

while read i
do
j=$(cat A1_flanking_Columba_livia_Mesitornis_unicolor_against_Columbimorphea_nr_blast.xml | xtract -pattern Iteration -def "-" -tab "\n" -sep "\t" -element Iteration_query-def,Iteration_query-len,Hit_def,Hsp_score,Hsp_evalue,Hsp_align-len \
| sed 's/ /_/g' | grep "$i" | awk '{print$1}' | awk -F "_" '{print$NF}' | sort -u)
echo $i $j
unset j
done < <(awk '{print$1"_"$2}' extra_birds_closer_to_Pterocles_gutturalis) | column -t

