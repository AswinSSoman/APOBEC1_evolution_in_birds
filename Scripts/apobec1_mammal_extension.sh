##################################################################################################################################################################################################################################################################################################################
#COMPARE THE MAMMALIAN & AVIAN APOBEC1
##################################################################################################################################################################################################################################################################################################################

#It is already known that APOBEC1 of mammals is longer than non-mammalian groups like lizard & zebra finch.
#Here we are trying to find what is the difference b/w mammalian & avian APOBEC1, where is the extension, how mammals gt the extension?

##################################################################################################################################################################################################################################################################################################################
#WORKFLOW for choosing representative mammalian species
										   _______________
												|
							All_annotated including isoforms(309)	|
								      |				|
								      v				|
								NCBI_curated (266)		|
								      |				|
								      v				|
								Complete_ORF (260)		|
								      |				|--> Amino acid sequence of APOBEC1 of mammals being filtered
								      v				|
								Length_filtered (259)		|
								      |				|
								      v				|
								Hybrids_removed (258)		|
								      |				|
								      v				|
						       Duplicates i.e. isoforms removed (185)	|
								      |		   _______________|
								      v	
								RNA-seq for outlier (Not )
								      |
								      v
					        _-----------------------------------------------_                         
					        |                                               |
Alignment manual filtering			Amino acid alignment (156)    		   CDS alignment (150)
  (alignment outliers)			        |                                               |
					        v                                               v
										CDS alignment 2nd round (108)
											  |
											  v
									Time tree filtering (32) 


Time tree filtering	;
  - based on - 1 species from the group that diverged atleast before 35 mya 
  	    - species with highest no of entres records in NCBI for - No of genes
							      - SRA experiments
	   					               - GEO datasets
  - advantages :
    - ensure the chosen representative species of mammals cover a broad taxonomic range within a single group have enough phylogenetic distance to 
    - Consider the diversity within the group & look for a wide range of species that represent different taxonomic families & phylogenetic relationships
																				    
##################################################################################################################################################################################################################################################################################################################
#Download APOBEC1 data of mammals

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download APOBEC1 metadata of mammals

cd~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension
esearch -db gene -query "mammals[ORGN] AND APOBEC1[Gene]" | efetch -format docsum  > apobec1_ncbi_mammals_esummary_raw
cat apobec1_ncbi_mammals_esummary_raw | xtract -pattern DocumentSummary -tab "\t" -def "N/A" -element Name ScientificName CommonName Description Status Chromosome OtherAliases ExonCount \
| sed 's/ /_/g' | sed '1i Gene Scientific_name Common_name Description Status Chromosome Othe_aliases Exon_count' | sed 's/\(DISCONTINUED\)[^\t]\+/\1/g' | column -t > APOBEC1_ncbi_mammals_esummary_table

#Filter only curated genes & fetch exon-count & protein length info of them
cat APOBEC1_ncbi_mammals_esummary_table | awk '$1=="APOBEC1" && $6!~"DISCONTINUED"' | awk '{print$2}' | xargs -n1 sh -c 'echo ">"$0; esearch -db gene -query "$0[ORGN] AND APOBEC1[GENE]" | efetch -format gene_table' > apobec1_ncbi_mammals_basic_info_raw
egrep "^>|length" apobec1_ncbi_mammals_basic_info_raw | egrep "^>|protein" | awk 'BEGIN{cmd="sort -V | head -1"} /^>/{close(cmd); print; next} {print | cmd}' | sed 's/^.*\([0-9]\+\) coding  exon.*AA length: \([0-9]\+\)/\1 \2/g' | awk -v RS=">" '{$1=$1}1' | column -t > APOBEC1_ncbi_mammals_ortholog_info

#Supporting evidence 
awk '{print$2}' APOBEC1_ncbi_mammals_esummary_table | xargs -n2 sh -c 'echo ">"$0; esearch -db nuccore -query "APOBEC1[GENE] AND $0[ORGN]" | efetch -format gb -mode xml' > apobec1_rna_data_raw
cat apobec1_rna_data_raw | egrep "^>|Supporting" | tr "<" "\n" | sed 's/.*Supporting/Supporting/g' | grep -v "/GBQualifier_value" | awk NF | sed 's/ /_/g' | sed '/>/! s/^/- /g' | column -t > supporting_evidence
cat raw_info/apobec1_rna_data_raw | egrep "^>|Supporting|GBSeq_accession-version>XM_|GBSeq_accession-version>NM_|GBSeq_accession-version>XR_" | sed 's/^>/@/g' | sed 's/>/\n/g' | sed 's/</\n/g' | egrep -v "GBSeq_accession-version|GBQualifier_value" | awk NF | sed 's/ /_/g' | sed 's/.*Supporting/Supporting/g' | sed '/^@/! s/^/- /g' | awk '/@/ {print++i"."$0;next}{print$0}' | sed 's/@//g' | column -t | less
awk '/100%/{print$1}' RS=">" supporting_evidence | awk NF > mammals_with_100%_rna_cov
cat apobec1_rna_data_raw | xtract -pattern GBSet -element GBSeq_primary-accession GBSeq_length GBSeq_organism | less

#Download gene info from human as a reference
esearch -db gene -query "Homo_sapiens[ORGN] AND APOBEC1[GENE]" | efetch -format gene_table | egrep "^>|length" | egrep "^>|protein" | awk 'BEGIN{cmd="sort -V | head -1"} /^>/{close(cmd); print; next} {print | cmd}' | sed 's/^.*\([0-9]\+\) coding  exon.*AA length: \([0-9]\+\)/\1 \2/g' | awk -v RS=">" '{$1=$1}1' | sed 's/^/Homo_sapiens /g' | column -t > human_ncbi_apobec1

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download APOBEC1 sequence of mammals

#Download CDS using edirect
#From browser
sed 's/>.*PREDICTED: />/g' APOBEC1_refseq_transcript.fasta | sed 's/ apolipo.*//g' | sed 's/NM_[^ ]\+ //g' | sed 's/ /_/g' | sed 's/\r//g' | awk NF > APOBEC1_mammals.fa

#Something is wrong with this sequence
#From command-line
datasets download ortholog symbol APOBEC1 --exclude-gene --exclude-rna --include-cds
unzip ncbi_dataset.zip
cut -f2 -d "[" ncbi_dataset/data/cds.fna | sed 's/organism=/>/g' | sed 's/]//g' > datasets_apobec1.fa

#Sequnces from previous literature 
#From https://ftp.ncbi.nlm.nih.gov/pub/aravind/AAD/apobec.html#NAD1
sed 's/^>[^_]\+_/>/g' d > krishan_apobec1.aa

##################################################################################################################################################################################################################################################################################################################
#Choose representative mammals

#Approach : Download all annotations of the APOBEC1 gene including all isoforms & accession versions -> Remove species based on QC -> Find mammalian APOBEC1's best showing the extension of APOBEC1 w.r.t birds with curation  

#Download protein fasta from NCBI website
#NCBI -> APOBEC1 in gene -> Download all refseq protein sequences per gene
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cluster
sed 's/\r//g' refseq_protein.fasta | awk NF | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' | sed 's/ .*\[/[/g' | sed 's/ /_/g' | sed 's/>/>M_/g' > APOBEC1_mammals.aa
cat APOBEC1_mammals.aa <(sed '/>/ s/>/>B_/g' ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.aa ) > birds_mammals_APOBEC1.aa
sed 's/\*/X/g' birds_mammals_APOBEC1.aa -i
#Cluster all the sequnces using clans to quick check how much variability is observed within mammals & across birds & mammals

#=================================================================================================================================================================================================================================================================================================================
#Remove mammalian species based existing metadata 

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get Supporting evidence :
 
 #Types of supporting evidences: Similarity to 
   #- ESTs : Number of Expressed Sequence TAGs
   #- long_SRA_reads
   #- mRNAs
   #- Proteins
   #- RNAseq : % coverage of the annotated genomic feature by RNAseq alignments
           # : samples with support for all annotated introns
  #note : NM accessions are validated supporting evidences are not given in the gb format

#get transcript ID from protein ID
echo -n > transcript_id_raw
for i in `cat ../raw_info/refseq_protein.fasta | grep ">" | awk '{print$1}' | tr -d ">"`
do
echo $i
echo "@@"$i >> transcript_id_raw
#direct efetch
efetch -db protein -id "$i" -format gb -mode xml >> transcript_id_raw
done

#Fetch metadata from transcript ID 10m32.122s
echo -n > transcript_rna_raw
while read i
do
i1=`echo $i | awk '{print$1}'`
i2=`echo $i | awk '{print$2}'`
echo $i1 $i2
echo "@@"$i1 $i2 >> transcript_rna_raw
efetch -db nuccore -id "$i2" -format gb -mode xml >> transcript_rna_raw
done < <(cat transcript_id_raw | egrep "^@@|GBSeq_source-db" | sed 's/>/\n/g' | sed 's/</\n/g' | grep -v "GBSeq_source-db" | sed 's/REFSEQ: accession//g' | awk NF | sed 's/@@//g' | paste -d " " - -)

#Create summary table
while read i
do 
i1=`echo -e "$i" | tr "\t" "\n" | head -1 | cut -f1,2 -d "_"`
i2=`echo -e "$i" | tr "\t" "\n" | grep "GBSeq_organism" | sort -u | sed 's/>/\n/g' | sed 's/</\n/g' | grep -v "GBSeq_organism" | sed 's/^_$//g' | awk NF`
i3=`cat transcript_id_raw | egrep "^@@|GBSeq_source-db" | sed 's/>/\n/g' | sed 's/</\n/g' | grep -v "GBSeq_source-db" | sed 's/REFSEQ: accession//g' | awk NF | sed 's/@@//g' | paste -d " " - - | grep "$i1" | awk '{print$1}'`
i4=`echo -e "$i" | tr "\t" "\n" | grep "Supporting" | sed 's/>/\n/g' | sed 's/</\n/g' | egrep -v "GBSeq_accession-version|GBQualifier_value" | sed 's/.*Supporting_evidence_includes_similarity_to:_//g' | sed 's/^_$//g' | awk NF`
i5=`echo -e "$i" | tr "\t" "\n" | grep "GBSeq_taxonomy" | sed 's/>/\n/g' | sed 's/</\n/g' | grep -v "GBSeq_taxonomy" | awk NF | sed 's/.*Mammalia;//g' | awk -F ";" '!($NF="")' \
| sed 's/.*Afrotheria.*/Afrotheria/g' | sed 's/.*Rodentia.*/Rodentia/g' | sed 's/.*Lagomorpha.*/Lagomorpha/g' | sed 's/.*Primates.*/Primates/g' | sed 's/.*Artiodactyla.*/Artiodactyla/g' | sed 's/.*Carnivora.*/Carnivora/g' | sed 's/.*Chiroptera.*/Chiroptera/g' \
| sed 's/.*Perissodactyla.*/Perissodactyla/g' | sed 's/.*Eulipotyphla.*/Eulipotyphla/g' | sed 's/.*Xenarthra.*/Xenarthra/g' | sed 's/.*Metatheria.*/Metatheria/g' | sed 's/.*Monotremata.*/Monotremata/g' | sed 's/.*Pholidota.*/Pholidota/g' | sed 's/.*Scandentia.*/Scandentia/g' \
| sed 's/.*Dermoptera.*/Dermoptera/g'`
echo $i2 $i5 $i3 $i1 $i4
unset i1 i2 i3 i4 i5
done < <(cat transcript_rna_raw | sed 's/[ ]\+/_/g' | awk -v RS="@@" '{$1=$1}1' OFS="\t" | awk NF) | awk NF | sort -k1,1 | awk '{if($1==p) print"-","-",$3,$4,$5; else print$0; p=$1}' | column -t > supporting_evidence

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get basic status about protein (find if some proteins are replaced)

echo -n > protein_status
for i in `cat ../raw_info/refseq_protein.fasta | grep ">" | awk '{print$1}' | tr -d ">"`
do
echo $i
echo "@@"$i >> protein_status
#direct efetch
efetch -db protein -id "$i" -format docsum >> protein_status
done

#Sequnces that have been removed as a result of standard genome annotation processing
awk -v RS="@@" '/removed/{print$1}' protein_status > ncbi_removed_proteins

#Sequnces that are replaced by updated 
egrep "@@[^ ]*|[^;]*replaced[^;]*" transcript_id_raw -o | sed 's/.*replaced by /✅/g' | sed 's/.*replaced /❌/g' | sed 's/ /_/g' | awk -v RS="@@" '{$1=$1}1'| awk '$2!=""' | grep "❌[^ ]\+" -o | sed 's/,_/\n/g' | tr -d "❌" | sort -u | sed 's/\.$//g' > tmp
grep -if tmp <(grep "@@" transcript_id_raw) | tr -d "@@" > ncbi_replaced_proteins
grep -if ncbi_replaced_proteins <(egrep "@@[^ ]*|[^;]*replaced[^;]*" transcript_id_raw -o | sed 's/.*replaced by /✅/g' | sed 's/.*replaced /❌/g' | sed 's/ /_/g' | awk -v RS="@@" '{$1=$1}1'| awk '$2!=""' | column -t) -z

#Remove proteins with "XP_"(non-curated) if the same species have "NP_"(curated) accession also available 
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cluster
grep -if <(grep ">" birds_mammals_APOBEC1.aa | grep "NP_" | awk -F "[" '{print$2}' | tr -d "]" | sort -u) <(grep ">" birds_mammals_APOBEC1.aa) | grep "XP_" | cut -f1 -d "[" | sed 's/>M_//g' > ncbi_non_curated_proteins

cp ../supporting_info/ncbi_removed_proteins ../supporting_info/ncbi_replaced_proteins .
grep -ivf <(cat ncbi_removed_proteins ncbi_replaced_proteins ncbi_non_curated_proteins | sort -u) <(grep ">" birds_mammals_APOBEC1.aa) | tr -d ">" > ncbi_selected_proteins
seqtk subseq birds_mammals_APOBEC1.aa ncbi_selected_proteins > birds_mammals_APOBEC1_filtered_1.aa

#=================================================================================================================================================================================================================================================================================================================
#Remove species based on basic features

#based on ORF
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cluster
sed '/>/! s/.X./-#-/g' birds_mammals_APOBEC1_filtered_1.aa | grep "\-#\-" -B1 | grep ">" | tr -d ">" > incomplete_orf_proteins
seqtk subseq birds_mammals_APOBEC1_filtered_1.aa <(grep -ivf <(cut -f1 -d "[" incomplete_orf_proteins) <(grep ">" birds_mammals_APOBEC1_filtered_1.aa) | tr -d ">") > birds_mammals_APOBEC1_filtered_2.aa

#based on protein length
awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' OFS="\t" birds_mammals_APOBEC1_filtered_2.aa | awk NF | sort -k2 -n > protein_length
termgraph.py protein_length --histogram --bins 30 --format '{:.0f}' --color red | sed 's/^0$/▏ 0/g' | awk NF | paste -d " " - - 

#Remove extreme outliers manually before using IQR method
awk '{print $1,$2,p,q; p=$1; q=$2}' protein_length | awk '$3!=""' | awk '{print$0,$2-$4}' | column -t | awk '$NF>50' 
grep -v "M_NP_001396400.1" protein_length > protein_length_filtered_1
Rscript remove_outliers.r protein_length_filtered_1 F | sed '1,4d' | awk '{print$2}' > protein_length_filtered_2
seqtk subseq birds_mammals_APOBEC1_filtered_2.aa <(awk '{print$1}' protein_length_filtered_2) > birds_mammals_APOBEC1_filtered_3.aa

#Remove hybrids
grep ">" birds_mammals_APOBEC1_filtered_3.aa | grep "_X_" -i > hybrid_proteins
seqtk subseq birds_mammals_APOBEC1_filtered_3.aa <(grep -ivf <(cut -f1 -d "[" hybrid_proteins) <(grep ">" birds_mammals_APOBEC1_filtered_3.aa) | tr -d ">") > birds_mammals_APOBEC1_filtered_4.aa

#=================================================================================================================================================================================================================================================================================================================
#Remove species based alignments

#Manual inspection of birds & mammals APOBEC1 protein alignment alignment along with their supporting info

#Add suporting evidence info to the protein IDs
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/supporting_info
while read i
do 
i1=`echo -e "$i" | tr "\t" "\n" | head -1 | cut -f1,2 -d "_"`
i2=`echo -e "$i" | tr "\t" "\n" | grep "GBSeq_organism" | sort -u | sed 's/>/\n/g' | sed 's/</\n/g' | grep -v "GBSeq_organism" | sed 's/^_$//g' | awk NF`
i3=`cat transcript_id_raw | egrep "^@@|GBSeq_source-db" | sed 's/>/\n/g' | sed 's/</\n/g' | grep -v "GBSeq_source-db" | sed 's/REFSEQ: accession//g' | awk NF | sed 's/@@//g' | paste -d " " - - | grep "$i1" | awk '{print$1}'`
i4=`echo -e "$i" | tr "\t" "\n" | grep "Supporting" | sed 's/>/\n/g' | sed 's/</\n/g' | egrep -v "GBSeq_accession-version|GBQualifier_value" | sed 's/.*Supporting_evidence_includes_similarity_to:_//g' | sed 's/^_$//g' | awk NF | sed 's/coverage_of_the_annotated_genomic_feature_by_RNAseq_alignments/RNAseq/g' | sed 's/_samples_with_support_for_all_annotated_introns/_Intron_samples/g' | sed 's/_sample_with_support_for_all_annotated_introns/_Intron_samples/g' | sed 's/_and_//g' | sed 's/_including_//g' | sed 's/,_/,/g' | sed 's/_,/,/g'`
i5=`echo -e "$i" | tr "\t" "\n" | grep "GBSeq_taxonomy" | sed 's/>/\n/g' | sed 's/</\n/g' | grep -v "GBSeq_taxonomy" | awk NF | sed 's/.*Mammalia;//g' | awk -F ";" '!($NF="")' \
| sed 's/.*Afrotheria.*/Afrotheria/g' | sed 's/.*Rodentia.*/Rodentia/g' | sed 's/.*Lagomorpha.*/Lagomorpha/g' | sed 's/.*Primates.*/Primates/g' | sed 's/.*Artiodactyla.*/Artiodactyla/g' | sed 's/.*Carnivora.*/Carnivora/g' | sed 's/.*Chiroptera.*/Chiroptera/g' \
| sed 's/.*Perissodactyla.*/Perissodactyla/g' | sed 's/.*Eulipotyphla.*/Eulipotyphla/g' | sed 's/.*Xenarthra.*/Xenarthra/g' | sed 's/.*Metatheria.*/Metatheria/g' | sed 's/.*Monotremata.*/Monotremata/g' | sed 's/.*Pholidota.*/Pholidota/g' | sed 's/.*Scandentia.*/Scandentia/g' \
| sed 's/.*Dermoptera.*/Dermoptera/g'`
echo $i2 $i5 $i3 $i4
unset i1 i2 i3 i4
done < <(cat transcript_rna_raw | sed 's/[ ]\+/_/g' | awk -v RS="@@" '{$1=$1}1' OFS="\t" | awk NF) | awk NF | sort -k1,1 > name_with_supporting_evidence

#Rename the protein IDs to protein IDs with supporting info 
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cluster
cp birds_mammals_APOBEC1_filtered_4.aa birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa
for i in `grep ">" birds_mammals_APOBEC1_filtered_4.aa | grep "M_" | cut -f1 -d "[" | sed 's/>M_//g'`
do
j=`grep "$i" /home/neo/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/supporting_info/name_with_supporting_evidence | sed 's/[ ]\+/-/g' | awk '{if($0~"RNA" || $0~"NP_" || $0~"EST") print$0; else print"No_RNA-------------"$0}'`
echo $i $j
sed "/>/ s/.*$i.*/>$j/g" birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa -i
unset
done

#Species with multiple isoforms/seq IDs
grep ">" birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa | grep -v "B_" | sed 's/No_RNA-------------//g' | cut -f1 -d "-" | tr -d ">" | sort | uniq -c | awk '$1>1'

for i in `grep ">" birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa | grep -v "B_" | sed 's/No_RNA-------------//g' | cut -f1 -d "-" | tr -d ">" | sort | uniq -c | awk '$1>1' | awk '{print$2}'`
do
j=`grep "$i" birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa | sed 's/No_RNA-------------//g' | cut -f3 -d "-"`
echo $i $j
unset j
done | column -t > species_with_multiple_sequences

aliview /home/neo/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cluster/birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa

#=================================================================================================================================================================================================================================================================================================================
#Remove duplicates based on : Pairwise alignments all sequences from a each species with curated APOBEC1 mammal consensus -> keep the sequnce with highest alignment score -> remove all others 

#Create consensus of birds APOBEC1
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cluster/remove_duplicate_sequences
cp ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.aa birds_APOBEC1.aa
cp ../birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa .
mafft --auto --reorder --quiet birds_APOBEC1.aa > birds_APOBEC1_mafft.aln
em_cons birds_APOBEC1_mafft.aln -auto -stdout -plurality 0 -identity 0 | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' | sed 's/>.*/>birds_APOBEC1_consensus/g' > birds_APOBEC1_consensus.aa

#Create a consensnus sequnce from curated (NP) annotations of mammals 
mafft --auto --reorder --quiet <(seqtk subseq birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa <(grep "NP_" birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa | tr -d ">")) > mammal_NP_seq.aln
~/programmes/OD-Seq/OD-seq -i mammal_NP_seq.aln -r mammal_NP.out -o mammal_NP_seq.out
mafft --auto --reorder --quiet <(seqtk subseq birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa <(grep -ivf <(grep ">" mammal_NP_seq.out) <(grep ">" mammal_NP_seq.aln) | tr -d ">")) > mammal_NP_seq_filtered.aln
em_cons mammal_NP_seq_filtered.aln -auto -stdout -plurality 0 -identity 0 | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' | sed 's/>.*/>mammal_NP_seq_consensus/g' > mammal_NP_seq_consensus.aa

#Pairwise align each proteins duplicate sequences from a species with birs APOBEC1 conensus to find the species most closest to curated mammal annotation
while read i
do
i1=`echo $i | awk '{print$1}'`
echo ">"$i1
for j in `echo $i | awk '!($1="")'`
do
echo "  - "$j
awk -v RS=">" -v s="$j" '$0~s {print">"$0}' birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa | grep -v "^>$" | awk NF > tmp_seq
needle tmp_seq mammal_NP_seq_consensus.aa -auto -stdout -awidth 280 > "$j"_with_mammal_consensus.aln
rm tmp_seq
done
unset i1 j
done < ../species_with_multiple_sequences

#Summary table of pairwise alignment scores
while read i
do
i1=`echo $i | awk '{print$1}'`
k=$(for j in `echo $i | awk '!($1="")'`
do
j1=`cat "$j"_with_mammal_consensus.aln | egrep "Identity|Similarity|Gaps|Score" | cut -f2 -d ":" | awk '{print$NF}' | tr -d ")(%" | sed 's/\.0\b//g' | paste -s -d ","`
j2=`cat "$j"_with_mammal_consensus.aln | egrep "Length|Score" | cut -f2 -d ":" | awk '{print$NF}' | paste -s -d " " | awk '{print$2/$1}' | awk '{$NF+=0}1' CONVFMT="%.2f"`
echo $j"("$j1")" $j2
unset j1
done | sort -k2 -nr | paste -s -d " ")
echo $i1 $k
unset i1 j j1 j2 k
done < ../species_with_multiple_sequences | column -t > pairwise_alignment_species_all_scores 

sed 's/([^)]\+)//g' pairwise_alignment_species_all_scores | awk '{print$1, $2"("$3")", $4"("$5")", $6"("$7")", $8"("$9")", $10"("$11")", $12"("$13")"}' | sed 's/()//g' | column -t > pairwise_alignment_species_scores
#Distribution of alignment scores
termgraph.py <(awk '{print$3,$5,$7,$9,$11,$13}' pairwise_alignment_species_all_scores | tr " " "\n" | awk NF | nl) --histogram --bins 30 --format '{:.0f}' --color red | sed 's/^0$/▏ 0/g' | awk NF | paste -d " " - -

#Closest species
awk '{print$1,$2}' pairwise_alignment_species_scores | cut -f1 -d "(" | column -t > pairwise_alignment_closest_species

mkdir pairwise_alignments
mv *_with_mammal_consensus.aln pairwise_alignments

cat pairwise_alignments/*_with_mammal_consensus.aln | grep -v "^#" | sed -z 's/\n\n\n\n/\n\n/g' | less

#=================================================================================================================================================================================================================================================================================================================
#Remove duplicates based on : Multiple sequnce alignment of all sequences from a each species with APOBEC1 mammal consensus -> Build phylogeny -> keep the closest species with lowest distance to the mammal consensus

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cluster/remove_duplicate_sequences
mkdir msa

#Run MSA for every species (1m24.541s)
while read i
do
i1=`echo $i | awk '{print$1}'`
echo ">"$i1
i2=`echo $i | awk '!($1="")'`
mkdir -p msa/"$i1"
cd msa/"$i1"
seqtk subseq ../../birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa <(grep -if <(echo $i2 | tr " " "\n") ../../birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa | tr -d ">") > "$i1"_all_seq.aa
#MSA
mafft --auto --reorder --quiet <(cat ../../mammal_NP_seq_consensus.aa "$i1"_all_seq.aa) > "$i1"_msa.aln
#Build tree
iqtree -s "$i1"_msa.aln -quiet
#Find closest species based on a custom R script
Rscript ../../find_clostest_species.r "$i1"_msa.aln.treefile
unset i1 i2
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cluster/remove_duplicate_sequences
done < ../species_with_multiple_sequences

#Alignment scores
while read i
do
i1=`echo $i | awk '{print$1}'`
i2=`awk '$1=="mammal_NP_seq_consensus"' msa/"$i1"/all_possible_pairwise_distances | sort -k3 -n | sed 's/No_RNA-------------//g' | sed 's/-/ /g' | awk -F " " '{print$4,$NF}' | sort -k2 -n | awk '{print$1"("$2")"}' | paste -s -d " "`
echo $i1 $i2
unset i1 i2
done < ../species_with_multiple_sequences | column -t > msa_species_scores

#Closest species
while read i
do
i1=`echo $i | awk '{print$1}'`
i2=`cat msa/$i1/"$i1"_msa.aln.treefile_closest_species_pairs | tr -d '"' | awk '$1=="mammal_NP_seq_consensus"' | awk '{print$2}' | sed 's/No_RNA-------------//g' | sed 's/-/ /g' | awk '{print$3}'`
echo $i1 $i2
unset i1 i2
done < ../species_with_multiple_sequences | column -t > msa_closest_species

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Check if pairwise and MSA alignment gave same closest species

join -1 1 -2 1 pairwise_alignment_closest_species msa_closest_species | column -t > pairwise_and_msa_closest_species
grep -if <(awk '{print$2}' msa_closest_species) pairwise_and_msa_closest_species

#If closest species given by pairwise and MSA alignment are different then manually inspect them 
grep -if <(awk '{print$2}' msa_closest_species) pairwise_alignment_closest_species -v | awk '{print$2}' | xargs -n1 bash -c 'grep "$0" pairwise_alignment_species_scores <(cut -f2- -d " " msa_species_scores) --color=always -h | sed "s/^/|/g" | paste -d " " - -' > compare_alignment_scores
sed 's/^|//g' compare_alignment_scores | awk '{if(NF=="10") {gsub(/\|/,"- - | "); print$0" - -"} else if(NF=="8") {gsub(/\|/,"- - - | "); print$0" - - -"} else if(NF=="6") {gsub(/\|/,"- - - - | "); print$0" - - - -"} else print$0}' | column -t

i=Callithrix_jacchus
grep "$i" pairwise_alignment_species_all_scores
grep "$i" pairwise_alignment_species_all_scores | sed 's/([^)]\+)//g' | awk '!($1="")' | awk '{line=""; for (i=1;i<=NF;i+=2) line=line (" " $i); print line;}' | xargs -n1 bash -c 'grep "$0" <(readlink -f pairwise_alignments/*)' | xargs -n1 sh -c 'echo "$0" | awk -F "/" "{print\$NF}" | cut -f1,2 -d "_" | sed "s/^/>/g"; egrep "Identity|Similarity|Gaps|Score|^mammal_NP_seq" $0 -B2 | sed "s/^#/ /g"' | sed 's/^>/\n>/g' | GREP_COLORS="mt=07;33" grep ">XP.*\|$"
grep "^+" msa/"$i"/"$i"_msa.aln.iqtree -C1
awk '$1~"mammal"' msa/"$i"/all_possible_pairwise_distances | sort -k3 -n | column -t
alen msa/"$i"/"$i"_msa.aln

#Manual inspection shows species chosen by pairwise alignment score is correct than msa
grep -ivf <(awk '{print$2}' pairwise_alignment_closest_species) <(awk '!($1="")' ../species_with_multiple_sequences | tr " " "\n" | awk NF) > seq_to_remove

#Rempve duplicate sequences
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cluster
seqtk subseq birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa <(grep -ivf remove_duplicate_sequences/seq_to_remove <(grep ">" birds_mammals_APOBEC1_filtered_4_with_supporting_info.aa) | tr -d ">") > birds_mammals_APOBEC1_filtered_5_with_supporting_info.aa

#Visualize filtered sequnce alignments based on group
grep ">" birds_mammals_APOBEC1_filtered_5_with_supporting_info.aa | grep -v ">B_" | sed 's/No_RNA-------------//g' | sed 's/^>//g' | cut -f1,2 -d "-" > mammal_names_with_group_names
sed '/>B_/! s/No_RNA-------------//g' birds_mammals_APOBEC1_filtered_5_with_supporting_info.aa | awk -F "-" '{if($1~"^>" && $1!~"^>B_") print">"$2; else if($1~"^>B_") print">Birds"; else print$0}' | awk '(/^>/ && s[$0]++){$0=$0"_"s[$0]}1;' > birds_mammals_APOBEC1_filtered_5_group_names.aa

#Clustering using CLANS
java -Xmx4G -jar ~/softwares/clans.jar -load birds_mammals_APOBEC1_filtered_5_group_names.clans

#=================================================================================================================================================================================================================================================================================================================
#Manual alignment inspection for outliers

aliview birds_mammals_APOBEC1_filtered_5_with_supporting_info.aa

#aliview extreme outliers
```Dipodomys_ordii-Rodentia-XP_012879117.1-11_Proteins,100%_RNAseq,3_Intron_samples
Ochotona_curzoniae-Lagomorpha-XP_040834788.1-15_Proteins,100%_RNAseq,8_Intron_samples
No_RNA-------------Acomys_russatus-Rodentia-XP_051011606.1-3_Proteins
No_RNA-------------Peromyscus_californicus_insignis-Rodentia-XP_052577525.1-3_Proteins
No_RNA-------------Antechinus_flavipes-Metatheria-XP_051857722.1-3_Proteins
Sarcophilus_harrisii-Metatheria-XP_031796053.1-9_Proteins,74%_RNAseq
Tupaia_chinensis-Scandentia-XP_006146678.2-8_Proteins,100%_RNAseq,5_Intron_samples
Bos_indicus-Artiodactyla-XP_019815505.1-1_EST,2_Proteins,90%_RNAseq
Ursus_maritimus-Carnivora-XP_008695735.1-3_Proteins,86%_RNAseq
Leptonychotes_weddellii-Carnivora-XP_030889346.1-2_mRNAs,10_Proteins,28%_RNAseq
Equus_caballus-Perissodactyla-XP_014596027.2-7_Proteins,92%_RNAseq
Colobus_angolensis_palliatus-Primates-XP_011807988.1-13_mRNAs,58_ESTs,8_Proteins
Macaca_mulatta-Primates-XP_028685068.1-20_ESTs,3_Proteins,93%_RNAseq
Papio_anubis-Primates-XP_021777683.1-20_ESTs,3_Proteins,79%_RNAseq
Gorilla_gorilla_gorilla-Primates-XP_030856728.2-1_mRNA,56_ESTs,2_Proteins
Propithecus_coquereli-Primates-XP_012493246.1-1_mRNA,10_Proteins
Galeopterus_variegatus-Dermoptera-XP_008580026.1-7_Proteins,72%_RNAseq
Puma_yagouaroundi-Carnivora-XP_040318567.1-8_Proteins,87%_RNAseq
No_RNA-------------Chrysochloris_asiatica-Afrotheria-XP_006875523.1-8_Proteins
Prionailurus_bengalensis-Carnivora-XP_043420179.1-4_Proteins,92%_RNAseq
No_RNA-------------Jaculus_jaculus-Rodentia-XP_044993035.1-15_Proteins
No_RNA-------------Elephantulus_edwardii-Afrotheria-XP_006901828.1-11_Proteins
Dipodomys_spectabilis-Rodentia-XP_042557467.1-8_Proteins,100%_RNAseq,3_Intron_samples
Fukomys_damarensis-Rodentia-XP_010606145.1-8_Proteins,95%_RNAseq
Octodon_degus-Rodentia-XP_023563190.1-7_Proteins,93%_RNAseq,1_Intron_samples
Vicugna_pacos-Artiodactyla-XP_006218458.1-9_Proteins,98%_RNAseq,3_Intron_samples
Hipposideros_armiger-Chiroptera-XP_019481704.1-10_mRNAs,8_Proteins,18%_RNAseq
Felis_catus-Carnivora-XP_044918591.1-5_Proteins,92%_RNAseq
No_RNA-------------Panthera_uncia-Carnivora-XP_049482614.1-5_Proteins```

#Aliview dispensable outliers
No_RNA-------------Sorex_araneus-Eulipotyphla-XP_054999584.1-13_Proteins
Myodes_glareolus-Rodentia-XP_048300253.1-11_Proteins,100%_RNAseq,72_Intron_samples

#Remove the aliview outliers
seqtk subseq birds_mammals_APOBEC1_filtered_5_with_supporting_info.aa <(grep -ivf aliview_outliers <(grep ">" birds_mammals_APOBEC1_filtered_5_with_supporting_info.aa) | tr -d ">") > birds_mammals_APOBEC1_filtered_5_with_supporting_info_aliview_filtered.aa

##################################################################################################################################################################################################################################################################################################################
#CDS alignment of filtered protein

#=================================================================================================================================================================================================================================================================================================================
#Download CDS

cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cds_alignment
datasets download gene symbol apobec1 --include cds --ortholog all
unzip ncbi_dataset.zip
sed 's/:.*\[organism=/_/g' ncbi_dataset/data/cds.fna | sed 's/].*//g' | sed 's/ /_/g' | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' | less > APOBEC1_all_organisms.fa
awk '/^-/{$0=(x)substr($0,length(x)+1)}{x=$1}1' ../supporting_info/supporting_evidence | awk '{print$1,$3,$4}' | column -t > proteins_and_associated_transcript_ids
grep ">" APOBEC1_all_organisms.fa | cut -f1,2 -d "_" | tr -d ">" > transcript_ids
grep ">" ../cluster/birds_mammals_APOBEC1_filtered_5_with_supporting_info.aa | grep -v ">B_" | sed 's/No_RNA-------------//g' | cut -f3 -d "-" > filtered_5_proteins
seqtk subseq APOBEC1_all_organisms.fa <(grep -if <(grep -if filtered_5_proteins proteins_and_associated_transcript_ids | awk '{print$3}') APOBEC1_all_organisms.fa | tr -d ">") > filtered_5_transcripts.fa
cat <(sed 's/>/>B_/g' ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.fa) <(awk -F "_" '/^>/ {$1=$2=""; print">M_",$0;next} {print$0}' filtered_5_transcripts.fa | sed 's/ /_/g' | sed 's/[_]\+/_/g') > birds_mammals_filtered_5.fa

#=================================================================================================================================================================================================================================================================================================================
#Manual inspection of alignment for outliers

aliview birds_mammals_filtered_5.fa

#Descripencies observed from manual inspection

	#Extra extensions: M_Ornithorhynchus_anatinus, M_Tachyglossus_aculeatus
	#No extensions: M_Ochotona_curzoniae, M_Acomys_russatus, M_Galeopterus_variegatus, M_Leptonychotes_weddellii
	#Extra starting: M_Sarcophilus_harrisii, M_Dipodomys_ordii, M_Tupaia_chinensis, M_Puma_yagouaroundi
	#Gap @ starting : M_Equus_caballus, M_Macaca_mulatta, M_Papio_anubis, M_Colobus_angolensis_palliatus ,M_Gorilla_gorilla_gorilla 
	#Gap inbetween: M_Fukomys_damarensis, M_Antechinus_flavipes, M_Peromyscus_californicus_insignis, M_Bos_indicus
	#Extension inbetween: M_Elephantulus_edwardii, M_Talpa_occidentalis, M_Myodes_glareolus, M_Ornithorhynchus_anatinus, M_Tachyglossus_aculeatus, M_Dipodomys_spectabilis, M_Ursus_maritimus, M_Felis_catus, M_Prionailurus_bengalensis, M_Panthera_uncia, M_Vicugna_pacos, M_Octodon_degus ,M_Hipposideros_armiger
	#birds to remove - B_Malurus_cyaneus_samueli
	#short extra inbetween : M_Phodopus_roborovskii, M_Chrysochloris_asiatica

cat inspected_aliview_data | cut -f2 -d ":" | tr -d " " | tr "," "\n" > aliview_outliers
seqtk subseq birds_mammals_filtered_5.fa <(grep -ivf aliview_outliers <(grep ">" birds_mammals_filtered_5.fa) | tr -d ">") > birds_mammals_filtered_5_aliview_filtered.fa

#=================================================================================================================================================================================================================================================================================================================
#A higher level of filtering : more stringent similarity criteria + minimum number of species per group

cp birds_mammals_filtered_5_aliview_filtered.fa birds_mammals_filtered_5_aliview_filtered_with_group_names.fa
while read i
do
j=`grep "$i" mammal_names_with_group_names | sed 's/-/(/g' | sed 's/$/)/g'`
#echo $i $j
sed "/$i/ s/.*/>M_$j/g" birds_mammals_filtered_5_aliview_filtered_with_group_names.fa -i
unset j
done < <(grep ">M_" birds_mammals_filtered_5_aliview_filtered_with_group_names.fa | sed 's/>M_//g')

#check common names
for i in `grep -v "#" higher_filtering | cut -f1 -d "(" | cut -f2- -d "_"`
do
echo $i
echo @@"$i" >> higher_filtering_common_names_raw
esearch -db taxonomy -query "$i" | esummary >> higher_filtering_common_names_raw
done
egrep "^@@|CommonName" higher_filtering_common_names_raw | sed 's/ /_/g' | sed 's/CommonName//g' | tr -d "></" | paste -d " " - - | sed 's/@@//g' | column -t > higher_filtering_common_names
#manually add names of species
awk '$2==""' higher_filtering_common_names

grep -if <(grep -if <(grep -v "^B_" higher_filtering | cut -f1 -d "(" | cut -f2- -d "_" | grep -v "#") <(grep ">" birds_mammals_filtered_5_aliview_filtered.fa | cut -f2- -d "_") -v) birds_mammals_filtered_5_aliview_filtered_with_group_names.fa | tr -d ">" > higher_filtering_species_to_keep
seqtk subseq birds_mammals_filtered_5_aliview_filtered_with_group_names.fa higher_filtering_species_to_keep > birds_mammals_filtered_6.fa

#=================================================================================================================================================================================================================================================================================================================
#Filtering based on taxonomic redundancy: Max 5 (a specific number) number of species from taxonomic order

grep ">" birds_mammals_filtered_5_aliview_filtered_with_group_names.fa | grep "M_" | cut -f2 -d "(" | tr -d ")" | sort | uniq -c | sort -k1 -nr | awk '{print$2,$1}' | column -t > birds_mammals_filtered_5_aliview_filtered_with_group_names_and_no_of_species
grep ">" birds_mammals_filtered_6.fa | grep "M_" | cut -f2 -d "(" | tr -d ")" | sort | uniq -c | sort -k1 -nr | awk '{print$2,$1}' | column -t > birds_mammals_filtered_6_group_and_no_of_species

#For groups having number of species more than 5 : filter the group again based time tree divergence time
grep ">" birds_mammals_filtered_6.fa | cut -f1 -d "(" | cut -f2- -d "_" > birds_mammals_filtered_6_names
grep ">M_" birds_mammals_filtered_6.fa | cut -f1 -d "(" | cut -f2- -d "_" > birds_mammals_filtered_6_mammals

#get tree order for 108 mammals (birds_mammals_filtered_6_mammals) from time tree website
#nnao mammals_tree_order
sed 's/ /_/g' mammals_tree_order | tr -d "*" | xargs -n1 sh -c 'grep "$0" mammal_names_with_group_names' | sed 's/-/(/g' | sed 's/$/)/g' > birds_mammals_filtered_6_mammals_tree_order
sed 's/ /_/g' mammals_tree_order | tr -d "*" | xargs -n1 bash -c 'grep "$0" <(grep "No_RNA" ../cluster/birds_mammals_APOBEC1_filtered_5_with_supporting_info.aa)' | sed 's/>No_RNA-------------//g' | cut -f1 -d "-" > mammals_tree_order_with_no_rna
for i in `cat birds_mammals_filtered_6_mammals_tree_order`
do
i1=`echo $i | cut -f1 -d "("`
j=`grep "$i1" mammals_tree_order_with_no_rna`
if [[ $j == "" ]]
then
echo $i "-"
else echo $i "No_RNA"
fi
unset i1 j
done 
#Compare each group & take max 5 from each group, choose based on divergence times & the amount of data(entrez records like : gene, SRA, GEO, protein) available
#make a file with final chosen mammals
cat <(awk -v RS=">" '/B_/{print">"$0}' birds_mammals_filtered_6.fa | awk NF) <(seqtk subseq birds_mammals_filtered_6.fa <(awk '{print"M_"$1}' time_tree_ncbi_filtered_mammals)) > birds_mammals_filtered_7.fa

#Check for intact active site 
transeq birds_mammals_filtered_7.fa -auto -stdout | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' | sed 's/\*$/X/g' | sed 's/_1$//g' | egrep "H.E|PC..C" -z --color=always | awk -v RS=">" '{$1=$1}1' | awk NF | column -t

#Filter birds
mafft --auto --maxiterate 1000 --reorder --quiet birds_mammals_filtered_7.fa > birds_mammals_filtered_7.aln
alen birds_mammals_filtered_7.aln
#Because of uncertain exon boundary of Pelecanus crispus this bird is creating a gap in alignment, hence remove it
awk -v RS=">" '!/B_Pelecanus_crispus/ {print">"$0}' birds_mammals_filtered_7.fa | awk NF | grep -v "^>$" > birds_mammals_filtered_8.fa

grep ">" birds_mammals_filtered_8.fa | cut -f2- -d "_" | cut -f1 																														-d "(" > birds_mammals_filtered_8_names
##################################################################################################################################################################################################################################################################################################################
#Final alignment: protein & CDS

#MSA
mafft --auto --maxiterate 1000 --reorder --quiet birds_mammals_filtered_8.fa > birds_mammals_filtered_8.aln
#Codon alignment
guidance=/home/neo/guidance.v2.02/www/Guidance/guidance.pl
mkdir codon_alignment
time perl "$guidance" --program GUIDANCE --seqFile birds_mammals_filtered_8.fa --msaProgram MAFFT --seqType codon --outDir codon_alignment --genCode 1 --bootstraps 100
cp codon_alignment/MSA.MAFFT.aln.With_Names birds_mammals_filtered_8_codon.aln

#Amino acid alignment
transeq birds_mammals_filtered_8.fa -auto -stdout | sed 's/\*/X/g' | sed 's/_1$//g' | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > birds_mammals_filtered_8.aa
mafft --auto --maxiterate 1000 --reorder --quiet birds_mammals_filtered_8.aa > birds_mammals_filtered_8_aa.aln

#DNA alignment based on amino acid translation
java -jar ~/programmes/macse_v2.07.jar -prog alignSequences -seq birds_mammals_filtered_8.fa

#=================================================================================================================================================================================================================================================================================================================
#Visualize exntension 
 
#Based on 
	#DNA alignment using MAFFT : 
		#- The bird sequence stop w.r.t mammals at : 590-600
		#- The gap is from (590-600) till 720 
		#- Difference : 120-130 nt
	#Codon alignement from Guidance based on MAFFT
		#- The bird sequence stop w.r.t mammals at : 590-603
		#- The gaps : (590-603) till 720
		#- Difference : 120-127 nt
	#Amino acid alignement using MAFFT
		#- The bird sequence stop w.r.t mammals at : 196-200
		#- The gaps : (196-200) till 240
		#- Difference : 40-44 aa

##################################################################################################################################################################################################################################################################################################################
#Find where the mamaml extension came from

#2 obvious possibilities:
  #1. Birds lost the extra sequence at the end of APOBEC1 gene, while mammals retained it
  #2. Mammals gained extension:
	#2.1. that is present in birds but outside CDS
	#2.2. that is not present in birds

#=================================================================================================================================================================================================================================================================================================================
#Alignment of birds downstream with mammals: 

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Birds are lacking a C-terminal region w.r.t mammals, hence it is possible that the birds shorter sequence is due to the aquisition of a premature stop codon in birds or extension of C-terminal in mammals via a downstream stop codon
#In both ways the downstream of APOBEC1 stop in birds should be similar to mammalian APOBEC1 extension

cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/birds_downstream
cat ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.fa | grep ">" | tr -d ">" | cut -f1,2 -d "_" | xargs -n1 bash -c 'paste <(grep "$0" /home/neo/bird_db1/aswin/APOBEC1/*queries -l) <(grep "$0" ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.fa | tr -d ">")' | egrep -v "total_validated_queries|qced_queries" > validated_birds
grep -r -f <(cat ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.fa | grep ">" | tr -d ">") /home/neo/bird_db1/aswin/APOBEC1/*queries | grep -v "total_validated_queries" > pipeline_validated
grep -ivf <(cut -f2 -d ":" pipeline_validated) <(grep ">B_" ../cds_alignment/birds_mammals_filtered_5_aliview_filtered.fa) | sed 's/>B_//g' | sed 's!^!~/bird_db1/aswin/APOBEC1/total_validated_queries:!g' > manual_validated

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Extract birds APOBEC1 CDS sequnces but extend the last exon approximately till the C terminal extnsion length of mammalian APOBEC1

while read i
do
i1=`echo $i | awk '{print$NF}'`
i2=`echo $i | awk -F "/" '{print$7}' | cut -f1 -d " "`
echo ">"$i1 $i2
awk -v RS=">" '!/exon_5\>/ {print">"$0}' ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_"$i1".fa | grep -v "^>$" | awk NF | grep -v ">" | paste -s -d ""
awk -v RS=">" '/exon_5\>/ {print">"$0}' ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_"$i1".fa | grep -v "^>$" | awk NF > exon_5
genome=`ls ~/soft_links/"$i1"/genome/GC*.fna`
if [[ "$i2" == "2nd_gblast_validated_queries" ]]
then 
bed=`cat ~/soft_links/"$i1"/aswin/APOBEC1/2nd_gblast/test.out.bed | grep "exon_5" | sed 's/plus/+/g' | sed 's/minus/-/g' | awk '{if($NF=="+") print$1,$2,$3+250,$4,"1",$5; else if($NF=="-") print$1,$2-250,$3,$4,"1",$5}' OFS="\t"`
bedtools getfasta -fi "$genome" -bed <(echo -e "$bed") -s | grep -v ">" | sed '1i >exon_5_ext' > exon_5_ext
elif [[ "$i2" == "exonerate_validated_queries" || "$i2" == "exoneredit_validated_queries" || "$i2" == "gedit_validated_queries" || "$i2" == "spaln_validated_queries" || "$i2" == "annotated_queries" ]]
then
bed=`cat ~/soft_links/"$i1"/aswin/APOBEC1/test.out.bed | grep "exon_5" | sed 's/plus/+/g' | sed 's/minus/-/g' | awk '{if($NF=="+") print$1,$2,$3+250,$4,"1",$5; else if($NF=="-") print$1,$2-250,$3,$4,"1",$5}' OFS="\t"`
bedtools getfasta -fi "$genome" -bed <(echo -e "$bed") -s | grep -v ">" | sed '1i >exon_5_ext' > exon_5_ext
elif [[ "$i2" == "toga_validated_queries" ]]
then
bed12ToBed6 -i <(cat ~/bird_db1/aswin/APOBEC1/TOGA/"$i1"/toga_output/for_vizualization/query_annotation.bed | grep "rna-XM_027449866.2") > bed_tmp
if [[ $(awk '{print$NF}' bed_tmp | sort -u) == "-" ]]
then
bed=`sort -k2 bed_tmp -nr | awk '{print$1,$2,$3,"exon_"++i,"1",$6}' OFS="\t" | grep "exon_5" | awk '{if($NF=="+") print$1,$2,$3+250,$4,$5,$6; else if($NF=="-") print$1,$2-250,$3,$4,$5,$6}' OFS="\t"`
else
bed=`awk '{print$1,$2,$3,"exon_"++i,"1",$6}' OFS="\t" bed_tmp | grep "exon_5" | awk '{if($NF=="+") print$1,$2,$3+250,$4,$5,$6; else if($NF=="-") print$1,$2-250,$3,$4,$5,$6}' OFS="\t"`
fi
bedtools getfasta -fi "$genome" -bed <(echo -e "$bed") -s | grep -v ">" | sed '1i >exon_5_ext' > exon_5_ext
rm bed_tmp
fi
needle exon_5 exon_5_ext -auto -stdout -awidth 100000 > align_tmp
#cat align_tmp | grep "^exon_5_ext" -B2
start=`cat align_tmp | grep "^exon_5_ext" -B2 | grep "|" | grep "|" -aob | head -1 | cut -f1 -d ":" | awk '{print$1+1}'`
cat align_tmp | grep "^exon_5_ext" -B2 | cut -c $start- | awk 'END{print$1}'
unset i1 i2 bed genome start end
rm exon_5 exon_5_ext align_tmp
done < validated_birds | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > APOBEC1_birds_extended.fa

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Align and filter the species based on manual inspection 
cat <(sed 's/^>/>B_/g' APOBEC1_birds_extended.fa) <(awk -v RS=">" '/^M_/ {print">"$0}' ../cds_alignment/birds_mammals_filtered_8.fa | awk NF) > APOBEC1_birds_extended_mammals_filtered_8.fa

#Check alignment 
transeq APOBEC1_birds_extended_mammals_filtered_8.fa -auto -stdout | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' | grep "*" -z
mafft --auto --quiet --reorder <(transeq APOBEC1_birds_extended_mammals_filtered_8.fa -auto -stdout) > APOBEC1_birds_extended_mammals_filtered_8.aln

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Check conservation scores
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/birds_downstream/conservation
cp ../APOBEC1_birds_extended_mammals_filtered_8.fa .
mafft --auto --maxiterate 1000 --quiet --clustalout --reorder <(transeq APOBEC1_birds_extended_mammals_filtered_8.fa -auto -stdout) > APOBEC1_birds_extended_mammals_filtered_8.aln

#Since even adding 250bp region downsrteam doesn't show any reliable alignment its better to see chromosome alignments

##################################################################################################################################################################################################################################################################################################################
#Chromosome alignments

#Before creating chromsome alignments, check the alignment already available in UCSC browser

cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/chromosome_alignments


##################################################################################################################################################################################################################################################################################################################
#Check Protein features

cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/protein_features

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Hydrophobicity

cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/protein_features/hydrophobicity
cp ../../cds_alignment/birds_mammals_filtered_8.aa .
sed 's/X$//g' birds_mammals_filtered_8.aa -i

#Hydrophobicity index
while read i
do
i1=`echo $i | awk '{print$1}'`
i2=`echo $i | awk '{print$2}'`
j=`Rscript hydrophobicity_index.r $i2 | awk '{print$2}'`
echo $i1 $j
unset i1 i2 j
done < <(awk -v RS=">" '{$1=$1}1' birds_mammals_filtered_8.aa | awk NF | sed 's/X$//g') | column -t > hydrophobicity_index_birds_mammals_filtered_8_aa

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Hydrophobicity plot
Rscript plot_hydrophobicity.r


##################################################################################################################################################################################################################################################################################################################
#DOMAIN SEARCH : to find what all domains are present in the final APOBEC1 sequences & distinguish the difference is domains of birds & mammals especially at the C-terminal

#Availbale tools

#Usin protein sequence as query:
  #Motif finder: multiple sequences
  #Motif Scan: one sequence at a time
  #Motif search
  #SMART: 
  #NCBI Batch CD: multiple sequences

#Works on both DNA & protein
  #NCBI CDD - but one sequence at a time
  #ThreaDom online- but one sequence at a time
  #InterPro/InterProScan

cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/domain_search
cp ../cds_alignment/birds_mammals_filtered_8.fa .

#Using NCBI Batch CD-search
#In concise results : 
	#- Birds   : APOBEC4-like
	#- Mammals : NAD1 like

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Using moitif finder : 

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Using Interpro : Image download can done for only one species at a time

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Using PROSITE : works on only 10 sequences at a time

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Using SMART : works on one sequence at a time
#Only tool here that identified the Leucine Rich Repeat at the C terminal of mammals 


sed 's/No_RNA-------------//g' birds_mammals_APOBEC1_filtered_5_with_supporting_info.aa | awk -F "-" '{print$1}' | sed '/^>B_/! s/>/>M_/g'


##################################################################################################################################################################################################################################################################################################################
#Literature
##################################################################################################################################################################################################################################################################################################################

#-----------------------------------------------------------------------------
#C-terminal extension of APOBEC1 w.r.t to other APOBEC members

'APOBEC1 cytosine deaminase activity on single-stranded DNA is suppressed by replication protein A (2020 Dec)'
#A1 also contains an ∼40 amino acid (a.a.) C-terminal domain extension that is not present in any other family members, except A2, which has no known deaminase activity

'The AID/APOBEC family of nucleic acid mutators (2008)'
#Like AID, APOBEC1 acts in the nucleus and shuttles between cytoplasm and nucleus by virtue of an amino-terminal NLS and a carboxy-terminal NES
#The main difference between APOBEC1 and the other AID/APOBEC genes is an extended coding sequence at its 3' end, whose significance has yet to be understood. (later discovered)

'The structure of APOBEC1 and insights into its RNA and DNA substrate selectivity'
#The structure of an individual subunit reveals that the N-terminal APO1 deaminase domain (residues 15–187) has the typical core fold ofan APOBEC deaminase domain, and its C-terminal domain extension has a novel fold that is unique to APO1.

#-----------------------------------------------------------------------------
#C-terminal extension of some mammalian APOBEC1 w.r.t other mammalian APOBEC1

'C→U editing of apolipoprotein B mRNA in marsupials: identification and characterisation of APOBEC–1 from the American opossum Monodelphus domestica'
#short C-terminal extension is found in humans, rabbits(23) and opossum, but not in rodents. 
#In addition, in common with humans and rabbits, marsupial APOBEC-1 has an extended C-terminus seven amino acids longer than that in the rodents, suggesting that this feature was present in the ancestral enzyme. 
#We show here that substitution of the C-terminus of the rodent enzyme for that of the opossum enzyme makes no difference to catalytic activity in vitro , suggesting that this C-terminal extension has little functional significance.

'The antiretroviral potency of APOBEC1 deaminase from small animal species'
#Phlyogenetic tree analyses revealed that the rabbit A1 gene is related to primate A1 genes, while A1s from rodents form a single separate cluster

#-----------------------------------------------------------------------------
#Importance of C-terminal

'Flow-cytometric visualization of C>U mRNA editing reveals the dynamics of the process in live cells'
#Considering the importance of the C-terminal domain of mammalian APOBEC1 for RNA editing,51,53 we tested a fusion protein in which we stitched the rat C-terminal domain to the reptilian homolog. Also in this case there was no evidence of editing

#-----------------------------------------------------------------------------

#Example for literatures using MSA to show extension
'Ô|Identification of GRY-RBP as an Apolipoprotein B RNA-binding Protein That Interacts with Both Apobec-1 and Apobec-1 Complementation Factor to Modulate C to U Editing*' #https://www.jbc.org/article/S0021-9258(19)34301-7/pdf



















##################################################################################################################################################################################################################################################################################################################
#DRAFT SCRIPTS
##################################################################################################################################################################################################################################################################################################################

#based on protein length
awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' OFS="\t" birds_mammals_APOBEC1_filtered_2.aa | awk NF | sort -k2 -n > protein_length
termgraph.py protein_length --histogram --bins 30 --format '{:.0f}' --color red | sed 's/^0$/▏ 0/g' | awk NF | paste -d " " - - 
egrep -v "M_Mus_musculus_3|M_Jaculus_jaculus_2" protein_length > protein_length_filtered_1
Rscript remove_outliers.r protein_length_filtered_1 F | sed '1,4d' | awk '{print$2}' > protein_length_filtered_2
seqtk subseq birds_mammals_APOBEC1_filtered_1.aa <(awk '{print$1}' protein_length_filtered_2) > birds_mammals_APOBEC1_filtered_2.aa

#Remove hybrids
grep ">" birds_mammals_APOBEC1_filtered_2.aa | grep "_X_" -i
awk -v RS=">" '!/M_Bos_indicus_x_Bos_taurus/ {print">"$0}' birds_mammals_APOBEC1_filtered_2.aa | grep -v "^>$" | awk NF > birds_mammals_APOBEC1_filtered_3.aa

#Based on alignment
muscle -in birds_mammals_APOBEC1_filtered_3.aa -out birds_mammals_APOBEC1_filtered_3_muscle.aln
mafft --auto --reorder --quiet birds_mammals_APOBEC1_filtered_3.aa > birds_mammals_APOBEC1_filtered_3_mafft.aln
~/programmes/OD-Seq/OD-seq -i birds_mammals_APOBEC1_filtered_3_muscle.aln -r muscle_qc.out -o birds_mammals_APOBEC1_filtered_3_muscle_outlier.out
~/programmes/OD-Seq/OD-seq -i birds_mammals_APOBEC1_filtered_3_mafft.aln -r mafft_qc.out -o birds_mammals_APOBEC1_filtered_3_mafft_outlier.out
cat birds_mammals_APOBEC1_filtered_3_muscle_outlier.out birds_mammals_APOBEC1_filtered_3_mafft_outlier.out | grep ">" | sort -u | tr -d ">" | paste -s -d "|" | sed 's!|!\\>|!g'
awk -v RS=">" '!/M_Callithrix_jacchus\>|M_Carlito_syrichta\>|M_Ictidomys_tridecemlineatus\>|M_Lutra_lutra\>|M_Meles_meles\>|M_Mustela_putorius_furo\>|M_Ornithorhynchus_anatinus_2\>|M_Tachyglossus_aculeatus_3/ {print">"$0}' birds_mammals_APOBEC1_filtered_3.aa | grep -v "^>$" | awk NF > birds_mammals_APOBEC1_filtered_4.aa 
#2nd iteration (mafft aligner gave more outliers in prebious iteration)
mafft --auto --reorder --quiet birds_mammals_APOBEC1_filtered_4.aa > birds_mammals_APOBEC1_filtered_4_mafft.aln
~/programmes/OD-Seq/OD-seq -i birds_mammals_APOBEC1_filtered_4_mafft.aln -r mafft_qc_2.out -o birds_mammals_APOBEC1_filtered_4_mafft_outlier.out
cat birds_mammals_APOBEC1_filtered_4_mafft_outlier.out | grep ">" | sort -u | tr -d ">" | paste -s -d "|" | sed 's!|!\\>|!g'
awk -v RS=">" '!/M_Aotus_nancymaae\>|M_Dipodomys_spectabilis_3\>|M_Echinops_telfairi_2\>|M_Mesocricetus_auratus\>|M_Monodelphis_domestica\>|M_Ornithorhynchus_anatinus\>|M_Puma_yagouaroundi\>|M_Sapajus_apella_2\>|M_Sorex_araneus_2\>|M_Tachyglossus_aculeatus\>|M_Tachyglossus_aculeatus_2/ {print">"$0}' birds_mammals_APOBEC1_filtered_4.aa | grep -v "^>$" | awk NF > birds_mammals_APOBEC1_filtered_5.aa 
#3rd iteration (mafft aligner gave more outliers in prebious iteration)
mafft --auto --reorder --quiet birds_mammals_APOBEC1_filtered_5.aa > birds_mammals_APOBEC1_filtered_5_mafft.aln
~/programmes/OD-Seq/OD-seq -i birds_mammals_APOBEC1_filtered_5_mafft.aln -r mafft_qc_3.out -o birds_mammals_APOBEC1_filtered_5_mafft_outlier.out
cat birds_mammals_APOBEC1_filtered_5_mafft_outlier.out | grep ">" | sort -u | tr -d ">" | paste -s -d "|" | sed 's!|!\\>|!g'
#This iteration didn't create any outlier - i.e. the all protein sequences are likely to be robust  

#manual alignment
cp birds_mammals_APOBEC1_filtered_5.aa group_named.aa
for i in `awk '{print$2}' ../taxonomy/mammal_taxonomy | sort -u`
do
j=`grep "$i" ../taxonomy/mammal_taxonomy | awk '{print$NF}' | sort -u`
echo $i $j
sed "s/$i/$j/g" group_named.aa -i
done 
sed -e 's/M_//g' -e 's/_[0-9]\+$//g' -e 's/B_.*/Birds/g' group_named.aa -i
awk -iinplace '(/^>/ && s[$0]++){$0=$0"_"s[$0]}1;' group_named.aa

aliview birds_mammals_APOBEC1_filtered_5.aa

#phylogeny
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/taxonomy
cat ../raw_info/refseq_protein.fasta | grep ">" | awk '{print$1}' | tr -d ">" | xargs -n1 sh -c 'echo "@@"$0; efetch -db protein -id "$0" -format xml' >> mammal_taxonomy_raw
awk '/@@/ {p=$0} {a[p]++} END{for(i in a) print i"\n"a[i]-1}' mammal_taxonomy_raw | paste -d " " - - | column -t | less
#manually create the list of major mammal groups
cat mammal_taxonomy_raw | egrep "^@@|Org-ref_taxname|Org-ref_common|OrgName_lineage" | sed 's/>/\n/g' | sed 's/</\n/g' | awk NF | egrep -v "Org-ref_taxname|Org-ref_common|OrgName_lineage" | sed 's/.*Mammalia;//g' > tmp
for i in `cat mammal_groups`
do
sed "s/.*$i.*/$i/g" tmp -i
done
sed 's/ /_/g' tmp -i
awk -v RS="@@" '{$1=$1}1' tmp | column -t > mammal_taxonomy 

#RNA-seq supporting data 
cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/supporting_info
echo -n > mammal_rna_raw
for i in `cat ../raw_info/refseq_protein.fasta | grep ">" | awk '{print$1}' | tr -d ">"`
do
echo $i
echo "@@"$i >> mammal_rna_raw
#esearch using elink
esearch -db protein -query "$i" | elink -target nuccore | efetch -format gb -mode xml >> mammal_rna_raw
done

while read i
do 
#echo -e "$i" | tr "\t" "\n"
i1=`echo -e "$i" | tr "\t" "\n" | head -1`
i2=`echo -e "$i" | tr "\t" "\n" | grep "GBSeq_organism" | sort -u | sed 's/>/\n/g' | sed 's/</\n/g' | grep -v "GBSeq_organism" | sed 's/^_$//g' | awk NF`
echo $i1 $i2
echo -e "$i" | tr "\t" "\n" | egrep "GBSeq_accession-version>XM_|GBSeq_accession-version>NM_|GBSeq_accession-version>XR_|Supporting" | sed 's/>/\n/g' | sed 's/</\n/g' | egrep -v "GBSeq_accession-version|GBQualifier_value" | sed 's/.*Supporting_evidence_includes_similarity_to:_//g' | sed 's/^_$//g' | awk NF | sed 's/^/- - - /g' 
#echo $i1 $i2
done < <(cat mammal_rna_raw | sed 's/[ ]\+/_/g' | awk -v RS="@@" '{$1=$1}1' OFS="\t" ) | awk NF | column -t > supporting_evidence

#Domain info
cat birds_mammals_APOBEC1_hitdata.txt | sed '1,6d' | sed 's/ /_/g' | column -t | awk '$2=="specific"' | awk '{if($1==p) print"-",$2,$3,$4,$5,$6,$7; else print$0; p=$1}' | column -t 

grep "^M_" protein_length | awk '{print$1}' | cut -f2- -d "_" | sed 's/_[0-9]\+//g' | sort -u | xargs -n1 bash -c 'paste <(echo ">"$0) <(esearch -db taxonomy -query "$0" | efetch -format xml | xtract -pattern TaxaSet -element Lineage)' >> Mammal_taxonomy

awk '{print$1}' birds_mammals_APOBEC1_1_length | cut -f2- -d "_" | sed 's/_[0-9]\+//g' | sort -u > unique_mammals
grep "^M_" aliview_filtered | cut -f2- -d "_" | sed 's/_[0-9]\+//g' | sort -u | xargs -n1 bash -c 'paste <(echo ">"$0) <(esearch -db taxonomy -query "$0" | efetch -format xml | xtract -pattern TaxaSet -element Lineage)' >> Mammal_taxonomy

Antechinus
Artiodactyla; Even-toed ungulates
Carnivora;
Chiroptera; Bats

Dasypodidae; Armadillos
Didelphinae; marsupials
Dromiciops marsupials
Elephantidae; asian elephant
Folivora; Xenarthra (placental mammals)

Glires;- Lagomorpha
       - Rodentia

Orycteropodidae; afrotheria
Perissodactyla; Odd-toed ungulates
Phascolarctos marsupials
Pholidota; pangolins
Primates;
Trichosurus marsupials
Vombatus marsupials

cat datasets_apobec1.fa ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.fa | sed 's/ /_/g' | sed 's/_$//g' > birds_mammals_APOBEC1.fa
awk -iinplace '/^>/ && a[$0]++ {$1=$1"_"a[$0]}1' birds_mammals_APOBEC1.fa
awk -F"\t" '{print$4"\t"$6}' raw_info/ncbi_dataset/data/data_table.tsv | sed 's/ /_/g' | sort -u | grep -v "_x_" | xargs -n2 sh -c 'echo ">$0"; efetch -db taxonomy -id "$1" -format xml' > raw_info/taxonomy

seqtk subseq APOBEC1_mammals.fa <(cat mammals_with_100%_rna_cov) > mammals_with_100%_rna_cov.fa

cat APOBEC1_mammals.fa ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.fa > birds_mammals_APOBEC1.fa


