##########################################################################################################################################################################################################################################################################################################################
#                                                                                                                                          Phylogenetic analysis of A1-like 
##########################################################################################################################################################################################################################################################################################################################

#Objective:
	#To understand the relationship b/w A1-like to other members of whole AID/APOBEC family
	#To find when or in which branch A1-like originated and which clades has this gene
	#To decipher the evolution of all A1's from lungfishes to mammals, sauropsides & amphibians

#Major steps:
	#1. Collection of CDS sequences from ncbi using datasets
	#2. Fetch taxonomic data (to use for filtering)
	#3.1. Filtering other members based on: pairwise similarity with refseq, supporting evidence (RNA, protein), taxonomy
	#3.2. Filtering A1 mammals based on: pairwise similarity with refseq, supporting evidence (RNA, protein), median length
	#4. Collect A1 from Sauropsids: from literature, ensembl (filter: manual annotation inspection, complete ORF), orthoDB database (birds excluded), ncbi (synteny + supporting evidence)
	#5. Collect A1-like locally from bird genomes
	#6. Manually curate A1-like sequences: synteny, A1 & A1-like overlap, length variatoins
	#7. Prepare all sequence (remove duplicates + complete ORF + anonymous bases + consecutive stop condons at end) & Codon alignment (Guidance + mafft) & Phylogeny (iqtree2)
	#8. Add (lunfish) & remove (some A3c mammals, A1 squamata) some sequences & build phylogeny again
	#9. Add more species:
			#9.1. Marsupials (filter: one isoform per species i,e, highest pairwise similarity with Ornithorhynchus_anatinus | 106 to 7)
			#9.2. Monotremes (filter: MSA one isoform per species | 5 to 2)
			#9.3. Amphibians (filter annotation & RNA expression)
	#10. Filtering again based on weird phylogeny placement (removed 2 species) & ran phylogeny 3 times:
			#codon alignment + 1000 bootstrap
			#codon + 5000 bootstrap
			#nt alignment + 5000 bootstrap)
	#11. Manually inspect: clade-wise alignment, annotation, genome gaps, repeats, splice sites
	#12. Run phylogeny based on 3 inputs:
			#codon alignment 
			#nt alignment
			#Trim nt alignment based on 4 criterias

#Filtering:
	
						AID_APOBEC_cds.fa (422)
#							    | remove duplicates
#							    v
					AID_APOBEC_cds_unique.fa (411)
#							    | reomve incomplete ORF
#							    v
				AID_APOBEC_cds_unique_complete_orfs.fa (395)
#							    | added lungfish A1, removed A3c mammals, A1 squamata
#							    v
			AID_APOBEC_cds_unique_complete_orfs_refined (380)
#							    | add marsupials, monotremes & amphibians, remove single existing marsupial (Pleurodeles_waltl)
#							    v
		AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians (390)
#							    | remove weird phylogenetic placement (ncbi_A1_mammals_XM_004869517)
#							    v
	AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined.fa (389)
#							    | Manually inpsect cladewise MSA + annotation quality (genome gaps) + (Ran phylogeny 3 times)
#							    v
AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa (379)


##########################################################################################################################################################################################################################################################################################################################
#1. Collection of CDS sequences

#1.1. Download from database : NCBI via commandline

cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download

#Download all 
  datasets download gene symbol apobec1 --ortholog all --include gene,rna,protein,cds,5p-utr,3p-utr,product-report --filename apobec1.zip
  datasets download gene symbol apobec1-like --ortholog all --include gene,rna,protein,cds,5p-utr,3p-utr,product-report --filename apobec1-like.zip
  datasets download gene symbol apobec2 --ortholog all --include gene,rna,protein,cds,5p-utr,3p-utr,product-report --filename apobec2.zip
  datasets download gene symbol apobec3h --ortholog all --include gene,rna,protein,cds,5p-utr,3p-utr,product-report --filename apobec3h.zip
  datasets download gene symbol apobec3c --ortholog all --include gene,rna,protein,cds,5p-utr,3p-utr,product-report --filename apobec3c.zip
  datasets download gene symbol apobec3a --ortholog all --include gene,rna,protein,cds,5p-utr,3p-utr,product-report --filename apobec3a.zip
  datasets download gene symbol apobec4 --ortholog all --include gene,rna,protein,cds,5p-utr,3p-utr,product-report --filename apobec4.zip
  datasets download gene symbol aid --ortholog all --include gene,rna,protein,cds,5p-utr,3p-utr,product-report --filename aid.zip

#1.1.2. Extract data
  unzip apobec1.zip -d apobec1
  unzip apobec2.zip -d apobec2
  unzip apobec3h.zip -d apobec3h
  unzip apobec3c.zip -d apobec3c
  unzip apobec3a.zip -d apobec3a
  unzip apobec4.zip -d apobec4
  unzip aid.zip -d aid

##1.2. Collect bird A1 CDS from manual vaildation
	mkdir /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manually_validated
	scp neo@172.30.1.174:~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_*.fa .
	ls APOBEC1_*.fa | xargs -n1 sh -c 'echo $0 | cut -f2- -d "_" | sed "s/\.fa//g" | sed "s/^/>validated_A1_Birds_/g"; grep -v ">" $0 | paste -s -d ""' | sed '/^>/! s/[a-z]/\U&/g' > validated_A1_Birds.fa

##########################################################################################################################################################################################################################################################################################################################
#2. Fetch taxonomic data (to use for filtering)

#Download vertebrates without restricting rank "i,e, it can include all levels" because some orthologs are from subspecies or some other sub ranks
	datasets download taxonomy taxon "vertebrates" --children --debug --filename vertebrates_children_wo_rank 2> debug_info
	unzip vertebrates_children_wo_rank.zip -d vertebrates_children_wo_rank
	cat vertebrates_children_wo_rank/ncbi_dataset/data/taxonomy_summary.tsv | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | column -t > all_vertebrates

#NOTE: the last column os the summary says if the scientific name is formal or informal which mainly means that the species is unclassified, hence exclude them 
#Check number of species per group
	awk '$NF=="TRUE"' all_vertebrates | awk '{print$10}' | sort | uniq -c | sort -k1,1nr

##########################################################################################################################################################################################################################################################################################################################
#3. Filtering:

#Strategy:
  #1. Extract gene ID's with curated transcripts i,e, "NM" accession
  #2. Find the refseq select transcript for each gene ID's
  #3. If the number of curated transcripts are atleast 30 and cover major clades such as mammals, birds, crocodiles, turtles, squamata & fishes, then just directly use this for phylogeny.
  #4. If the number of curated transcripts are less than 10, then use this curated ones as reference to choose more transcripts from other species.
#Relevant info in metadata
  #1. Gene metadata columns to consider:
        #Annotation_Assembly_Accession : Preferrable to have an annotated assembly
        #Annotation_Release_Name       : Preferrable to have an assembly
        #Common_Name                   : Preferrable to sequences which has common names of associated species which excludes the chance of including hybrid animals, artificial sequences etc.
        #Ensembl_GeneIDs               : Preferrable to have same sequence associated with another database, which suggests more reliability
        #Nomenclature_Authority        : Preferrable to have a seperate authority responsible for nomenclature, also suggests more reliability
        #Proteins                      : Preferrable to have genes with single protein, exclude the complexity in choosing an isoform or paralog
        #Symbol                        : Preferrable to exclude variants apobec2a. If the variant is a paralog then synteny needs to be checked
  #2. Gene-product metadata
        #Transcript_CDS_Sequence_Accession: Preference order is Refseq select, NM, XM
          #Refseq select : experimentally validated
          #NM/NP/NR      : curated
          #XM/XP/XR      : non-curated
        #Transcript_Transcript_Length: Preferrable to choose transcript length closer to the length of validated/curated transcript

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.1. Extract refseq select & print relevant metedata

  cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download
  for i in aid apobec1 apobec2 apobec3h apobec3c apobec3a apobec4
  do
  cd $i && echo ">"$i
  #Print metadata excluding columns with all row values are "-"
  #Each row represents gene
  dataformat tsv gene --inputfile ncbi_dataset/data/data_report.jsonl  | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' \
  | awk -v na="-" 'NR==1{for(i=1;i<=NF;i++)keep[i]=1;header=$0;next} {for(i=1;i<=NF;i++)if($i!=na)keep[i]=0;data[NR]=$0} END{n=split(header,hdr,FS);for(i=1;i<=n;i++)if(!keep[i])printf (i<n?"%s\t":"%s\n",hdr[i]);for(row=2;row<=NR;row++){n=split(data[row],fields,FS);for(i=1;i<=n;i++)if(!keep[i])printf (i<n?"%s\t":"%s\n",fields[i])}}' | column -t > gene_metadata
  #Each row represents exon
  dataformat tsv gene-product --inputfile ncbi_dataset/data/product_report.jsonl | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' \
  | awk -v na="-" 'NR==1{for(i=1;i<=NF;i++)keep[i]=1;header=$0;next} {for(i=1;i<=NF;i++)if($i!=na)keep[i]=0;data[NR]=$0} END{n=split(header,hdr,FS);for(i=1;i<=n;i++)if(!keep[i])printf (i<n?"%s\t":"%s\n",hdr[i]);for(row=2;row<=NR;row++){n=split(data[row],fields,FS);for(i=1;i<=n;i++)if(!keep[i])printf (i<n?"%s\t":"%s\n",fields[i])}}' | column -t > gene_product_metadata
  #download refseq select sequences
  mkdir refseq_select
  esearch -db nuccore -query "$i[GENE] AND Refseq_select[filter" | efetch -format fasta_cds_na > refseq_select/refseq_select_cds.fa
  #download refseq select metadata (genbank format)
  esearch -db nuccore -query "$i[GENE] AND Refseq_select[filter" | efetch -format gene_table > refseq_select/refseq_select_gba
  cd ../
  done

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.2. Print pairwise sequence similarity metrics b/w all cds and refseq select (13m14.470s)

	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download
	time for i in aid apobec1 apobec2 apobec3h apobec3c apobec3a apobec4
	do
	echo ">"$i
	cd $i/refseq_select
	#header for summary (customize based on number of refseq select transcripts available for each apobec members)
	header=$(grep ">" refseq_select_cds.fa | grep "NM_[0-9.]\+" -o | awk '{print"Refseq_select_"++i,"Length Identity Similarity Gaps Score"}' | paste -s -d " " | sed 's/^/Species GeneID transcript Rank Group Sci_Name_is_formal /g')
	#Gene ID column number
	gidc=$(head -1 ../gene_metadata | awk '{for(i=1;i<=NF;i++) if($i=="NCBI_GeneID") print i}')
	#Taxon ID column number
	tidc=$(head -1 ../gene_metadata | awk '{for(i=1;i<=NF;i++) if($i=="Taxonomic_ID") print i}')
	#Run forloop to get pairwise alignment metrics and taxonomic info for each CDS transcripts
	time while read -r d
	do
	#Species, gene ID, Transcript ID of each CDS
	sgt=$(echo $d | tr ":" " " | egrep "[XN]M_[0-9.]+|organism=[a-zA-Z ]+|GeneID=[0-9]+" -o | cut -f2 -d "=" | tr " " "_" | paste -s -d " " | awk '{print$2,$3,$1}')
	#Gene ID
	gid=$(echo $d | grep "GeneID=[0-9]\+" -o | cut -f2 -d "=")
	#Taxon ID
	tid=$(awk -v a="$gidc" -v b="$gid" -v c="$tidc" '{if($a==b) print $c}' ../gene_metadata | sort -u)
	#Get taxon info from taxonomic data downloaded from NCBI using datasets
	#tax=$(awk -v a="$tid" '{if($2==a && $NF=="TRUE") print$5,$10,$27;  else if($2==a && $NF=="FALSE") print$5,$10,"-"}' /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/taxonomy/all_vertebrates)
	tax=$(awk -v a="$tid" '{if($2==a && $19=="Crocodylia") print$5,"crocodiles",$NF; else if($2==a) print$5,$10,$NF}' /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/taxonomy/all_vertebrates | sed -e 's/\bbats\b/mammals/g' -e 's/\bcarnivores\b/mammals/g' -e 's/\beven-toed_ungulates_\&_whales\b/mammals/g' -e 's/\bodd-toed_ungulates\b/mammals/g' -e 's/\binsectivores\b/mammals/g' -e 's/\bmonotremes\b/mammals/g' -e 's/\bplacentals\b/mammals/g' -e 's/\bprimates\b/mammals/g' -e 's/\brabbits_\&_hares\b/mammals/g' -e 's/\brodents\b/mammals/g' -e 's/\bwhales_\&_dolphins\b/mammals/g' -e 's/\bmarsupials\b/mammals/g' -e 's/caecilians/amphibians/g' -e 's/\bfrogs_\&_toads\b/amphibians/g' -e 's/\bsalamanders\b/amphibians/g' -e 's/\bchimaeras\b/cartilaginous_fishes/g' -e 's/\bsharks_\&_rays\b/cartilaginous_fishes/g' -e 's/\bhagfishes\b/jawless_vertebrates/g' -e 's/\blampreys\b/jawless_vertebrates/g' -e 's/\bhawks_\&_eagles\b/birds/g' -e 's/\blizards_\&_snakes\b/squamates/g')
	if [[ $tax == "" ]]; then tax="-"; else :; fi
	#Run pairwise alignment w.r.t each refseq select transcripts
	for r in $(grep -wf <(grep ">" refseq_select_cds.fa | grep "NM_[0-9.]\+" -o) ../ncbi_dataset/data/cds.fna -o)
	do
	ssm=$(needle <(awk -v RS=">" -v a="$r" '{if($0~a) print">"$0}' ../ncbi_dataset/data/cds.fna | awk NF | /media/aswin/programs/myfasta -comb) <(seqtk subseq ../ncbi_dataset/data/cds.fna <(echo "$d")) --auto --stdout | egrep "Length|Identity|Similarity|Gaps|Score" | awk '{print$3}' | cut -f1 -d "/" | paste -s -d " ")
	echo $r $ssm
	unset ssm
	done | paste -s -d " " | sed "s/^/$sgt $tax /g"
	unset sgt gid tid tax r
	done < <(grep ">" ../ncbi_dataset/data/cds.fna | tr -d ">") | sed "1i $header" | column -t > cds_pairwise_seq_similiarity_with_refseq_select
	unset header gidc tidc d
	cd ../../
	done

find . -name "cds_pairwise_seq_similiarity_with_refseq_select" | xargs -n1 sh -c 'echo $0|cut -f2 -d "/"; awk "\$6==\"-\"" $0 | awk "\$10" | wc -l' | paste -d " " - - | column -t

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.3. Collect supporting info and compile all QC info for filtering, run this code section individually / seperately for each apobec's

cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download

time for family in aid apobec1 apobec2 apobec3h apobec3c apobec3a apobec4
do
echo ">QC:" $family
cd $family/refseq_select
#Combine top 25 transcripts with highest similarity score w.r.t each refseq select transcripts
	awk '$6=="TRUE"' cds_pairwise_seq_similiarity_with_refseq_select > filtered_cds_pairwise_seq_similiarity_with_refseq_select
#Get metadata info of selected transcripts
	gid=$(head -1 ../gene_product_metadata | awk '{for(i=1;i<=NF;i++) if($i=="Transcript_CDS_Sequence_Accession") print i}')
	for i in $(awk '{print$3}' filtered_cds_pairwise_seq_similiarity_with_refseq_select | sed "1i Transcript_CDS_Sequence_Accession")
	do
	awk -v a="$gid" -v b="$i" '{if($a==b) print}' ../gene_product_metadata
	done | column -t > gene_product_metadata_filtered_cds_pairwise_seq_similiarity_with_refseq_select
	unset gid
#Some Quality check
	#Check protein length variantion
	awk '!a[$10]++ && NR>1{print$0}' gene_product_metadata_filtered_cds_pairwise_seq_similiarity_with_refseq_select | awk '{print$28}' | ministat -n
	#Check protein count variation (might indicate an isoform or a paralog)
	awk '!a[$10]++ && NR>1{print$0}' gene_product_metadata_filtered_cds_pairwise_seq_similiarity_with_refseq_select | awk '{print$5}' | sort -nr | uniq -c
	#Check gene type & transcript type Protein coding
	awk '!a[$10]++ && NR>1{print$0}' gene_product_metadata_filtered_cds_pairwise_seq_similiarity_with_refseq_select | awk '{print$4}' | sort -nr | uniq -c
	awk '!a[$10]++ && NR>1{print$0}' gene_product_metadata_filtered_cds_pairwise_seq_similiarity_with_refseq_select | awk '{print$30}' | sort -nr | uniq -c
	#Check symbol
	awk '!a[$10]++ && NR>1{print$0}' gene_product_metadata_filtered_cds_pairwise_seq_similiarity_with_refseq_select | awk '{print$6}' | sort -nr | uniq -c
	awk '!a[$10]++ && NR>1{print$0}' gene_product_metadata_filtered_cds_pairwise_seq_similiarity_with_refseq_select | awk '{print$29}' | sort -nr | uniq -c
#Supporting evidence can also be fetched
	echo "   - fetch supporting data..."
	time awk '{print$3}' filtered_cds_pairwise_seq_similiarity_with_refseq_select | xargs -n1 sh -c 'echo "@@"$0; efetch -db nuccore -id "$0" -format gb -mode xml' > genbank_filtered_cds_pairwise_seq_similiarity_with_refseq_select
#Print summary
	for i in $(awk '{print$3}' filtered_cds_pairwise_seq_similiarity_with_refseq_select)
	do
	#i0=$(awk -v RS="@@" -v a="$i" '{if($0~a) print$0}' genbank_filtered_cds_pairwise_seq_similiarity_with_refseq_select | grep "GBSeq_taxonomy" | cut -f2 -d ">" | cut -f1 -d "<" | egrep "Coelacanthiformes|Dipnomorpha|Chondrichthyes|Cyclostomata|Actinopterygii|Mammalia|Crocodylia|Testudines|Aves|Sphenodontia|Squamata|Amphibia" -iwo | sort -u | sed -e 's/Coelacanthiformes/Coelacanth/g' -e 's/Dipnomorpha/Lungfishes/g')
	i0=$(awk -v RS="@@" -v a="$i" '{if($0~a) print$0}' genbank_filtered_cds_pairwise_seq_similiarity_with_refseq_select | grep "GBSeq_taxonomy" | cut -f2 -d ">" | cut -f1 -d "<" | egrep "Coelacanthiformes|Dipnomorpha|Chondrichthyes|Cyclostomata|Actinopterygii|Mammalia|Crocodylia|Testudines|Aves|Sphenodontia|Squamata|Amphibia" -iwo | sort -u | sed -e 's/Coelacanthiformes/Coelacanthes/g' -e 's/Dipnomorpha/Lungfishes/g' -e 's/Chondrichthyes/Cartilaginous_fishes/g' -e 's/Cyclostomata/Jawless_vertebrates/g' -e 's/Actinopterygii/Ray_finned_fishes/g')
	i1=$(awk -v RS="@@" -v a="$i" '{if($0~a) print$0}' genbank_filtered_cds_pairwise_seq_similiarity_with_refseq_select | grep "GBSeq_accession-version" | cut -f2 -d ">" | cut -f1 -d "<")
	i2=$(awk -v RS="@@" -v a="$i" '{if($0~a) print$0}' genbank_filtered_cds_pairwise_seq_similiarity_with_refseq_select | grep "Supporting" | cut -f2 -d ">" | cut -f1 -d "<" | sed 's/^.*Supporting/Supporting/g' | sed 's/Supporting evidence includes similarity to: //g' | sed 's/ /_/g')
	i3=$(awk -v RS="@@" -v a="$i" '{if($0~a) print$0}' genbank_filtered_cds_pairwise_seq_similiarity_with_refseq_select | grep "Evidence-Data-START" | cut -f2 -d ">" | cut -f1 -d "<" | tr ";" "\n" | awk NF | awk '/##Evidence\-Data\-START##/,/##RefSeq-Attributes-END##/' \
	| grep -v "END##" | sed 's/-START##//g' | sed '/#/! s/$/|/g' | awk -v RS="##" '{$1=$1}1' | awk NF | sed 's/ /=>/1' | sed 's/ /_/g' | paste -s -d "|")
	if [[ "$i2" == "" ]] && [[ "$i3" == "" ]]; then i2="-"; else :; fi
	echo $i $i0 $i1 $i2 "$i3"
	unset i0 i1 i2 i3 a
	done | sed '1i Input_accession Group Fetched_accession Supporting_evidence_includes_similarity_to' | column -t > supporting_evidence_filtered_cds_pairwise_seq_similiarity_with_refseq_select
#Combine sequence similarity & annotation supporting evidence for each CDS transcripts
	while read -r i
	do
	gid=$(echo $i | awk '{print$3}')
	sup=$(awk -v a="$gid" '{if($1==a) print$4}' supporting_evidence_filtered_cds_pairwise_seq_similiarity_with_refseq_select)
	echo $i $sup
	unset gid sup
	done < filtered_cds_pairwise_seq_similiarity_with_refseq_select | sed "1i $(head -1 cds_pairwise_seq_similiarity_with_refseq_select | sed 's/$/ Supporting_evidence/g')" | column -t > sequence_similarity_and_supporting_evidence
#Extract CDS (unfiltered)
	for i in $(awk '{print$3}' filtered_cds_pairwise_seq_similiarity_with_refseq_select)
	do
	awk -v RS=">" -v a="$i" '{if($0~a) print">"$0}' ../ncbi_dataset/data/cds.fna | awk NF
	done | /media/aswin/programs/myfasta -comb > $family".fa"
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download
	done
unset family

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.4. Apply the following filters to all members except apobec1:
	#1. Single cds per each gene ID with highest Pairwise similarity score
	#2. With 100%_coverage_of_the_annotated_genomic_feature_by_RNAseq_alignments
	#3. With_support_for_all_annotated_introns
	#4. With simialrity to Proteins
	#5. Take upto 20 sequences per major taxonomic group

	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download
	time for family in aid apobec1 apobec2 apobec3h apobec3c apobec3a apobec4
	do
	echo ">Filtering:" $family
	cd $family/refseq_select
	#Apply filters
	if [[ $family == "aid" || $family == "apobec2" ]]
	then
	#Since aid & apobec2 have a lot of sequences even after filtering in contrast to other members & these are outgroup/external references & we don't need high number of outgroups as they are just used to act as a control (I ran it once & observed)
	sort -k12,12nr sequence_similarity_and_supporting_evidence | awk '!a[$2]++' | awk '($NF~"100%_coverage_of_the_annotated_genomic_feature_by_RNAseq_alignments" && $NF~"with_support_for_all_annotated_introns" && $NF~"[0-9]+_Proteins") || ($NF~"Evidence-Data=")' | awk 'a[$5]++<10' | awk '{print$3}' > rna_protein_taxon_filtered_transcripts
	else
	sort -k12,12nr sequence_similarity_and_supporting_evidence | awk '!a[$2]++' | awk '($NF~"100%_coverage_of_the_annotated_genomic_feature_by_RNAseq_alignments" && $NF~"with_support_for_all_annotated_introns" && $NF~"[0-9]+_Proteins") || ($NF~"Evidence-Data=")' | awk 'a[$5]++<20' | awk '{print$3}' > rna_protein_taxon_filtered_transcripts
	fi
	for i in $(cat rna_protein_taxon_filtered_transcripts)
	do
	g=$(awk -v a="$i" '{if($3==a) print$5}' sequence_similarity_and_supporting_evidence)
	awk -v RS=">" -v a="$i" '{if($0~a) print">"$0}' ../ncbi_dataset/data/cds.fna | /media/aswin/programs/myfasta -comb | sed "/^>/ s/.*/>$family\_$g\_$i/g"
	unset g
	done | awk NF > rna_protein_taxon_filtered_transcripts.fa
	unset i
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download
	done

#Print summary of filtered sequences 
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download
	for family in aid apobec1 apobec2 apobec3h apobec3c apobec3a apobec4
	do
	cd $family/refseq_select
	i1=$(grep ">" rna_protein_taxon_filtered_transcripts.fa -c)
	i4=$(awk 'NR>1{print$5}' sequence_similarity_and_supporting_evidence | sort | uniq -c | sort -k1,1nr | sed 's/^[ ]\+//g' | tr " " "_")
	i5=$(grep "^>" rna_protein_taxon_filtered_transcripts.fa | cut -f2- -d "_" | sed 's/_[XN]M_.*//g' | sort | uniq -c | sort -k1,1nr | sed 's/^[ ]\+//g' | tr " " "_")
	i6=$(for j in $(echo "$i4"); do j1=$(echo $j | cut -f2- -d "_"); j2=$(grep "$j1" <(echo "$i5")); if [[ "$j2" == "" ]]; then j2="-"; else :; fi
	j3=$(/media/aswin/programs/myfasta -l rna_protein_taxon_filtered_transcripts.fa | grep "$j1" | awk '{print$NF}' | ministat -n 2> /dev/null | awk 'END{$1=""; print$0}' | sed 's/^ //g')
	if [[ "$j3" == "" ]]; then j3=$(/media/aswin/programs/myfasta -l rna_protein_taxon_filtered_transcripts.fa | grep "$j1" | awk '{print$NF}' | /media/aswin/programs/ministat_for_less_than_3_number.sh | tail -1); else :; fi
	j4=$(transeq rna_protein_taxon_filtered_transcripts.fa -trim --auto --stdout | /media/aswin/programs/myfasta -l | grep "$j1" | awk '{print$NF}' | ministat -n 2> /dev/null | awk 'END{$1=""; print$0}' | sed 's/^ //g')
	if [[ "$j4" == "" ]]; then j4=$(transeq rna_protein_taxon_filtered_transcripts.fa -trim --auto --stdout | /media/aswin/programs/myfasta -l | grep "$j1" | awk '{print$NF}' | /media/aswin/programs/ministat_for_less_than_3_number.sh | tail -1); else :; fi
	echo $j $j2 $j3 $j4; done | awk 'NR>1 {print"-",$0;next} {print$0}')
	paste <(echo $family) <(echo "$i6") <(echo $i1)
	unset i1 i4 i5 i6
	cd ../../
	done | sed '1i Member Groups_unfiltered Groups_filtered cds_N cds_Min cds_Max cds_Median cds_Avg cds_Stddev/Diff prot_N prot_Min prot_Max prot_Median prot_Avg prot_Stddev/Diff Total_No_filtered_seq' | column -t | GREP_COLORS="mt=01;07;33" grep "^Member.*\|$"  --color=always > filtering_summary
	
	find . -name "rna_protein_taxon_filtered_transcripts.fa" -o -name "sequence_similarity_and_supporting_evidence" | sort | xargs -n2 bash -c 'paste <(echo $0|cut -f2 -d "/") <(grep ">" $0 -c); paste <(awk "{print\$5}" $1 | sort | uniq -c | sort -k1,1nr | sed "s/^/- - /g") <(grep ">" $0 | cut -f2- -d "_" | sort | sed "s/_[XN]M_.*//g" | sort | uniq -c | sort -k1,1nr)' | sed '1i Gene Total_filtered_number_of_orthologs Total_number_of_orthologs-> Group Filtered_number_of_orthologs-> Group' | column -t

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.5. Collect all filtered sequences excluding A1

	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download
	touch /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/aid_a2_a3_a4_filtered.fa
	for family in aid apobec2 apobec3h apobec3c apobec3a apobec4
	do
	cd $family/refseq_select && echo ">"$family
	cat rna_protein_taxon_filtered_transcripts.fa | sed '/^>/! s/[a-z]/\U&/g' >> /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/aid_a2_a3_a4_filtered.fa
	cd ../../
	done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.6. Seperately collect & filter A1 from mammals 

	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/apobec1/refseq_select
	sort -k12,12nr sequence_similarity_and_supporting_evidence | awk '!a[$2]++' | awk '($NF~"100%_coverage_of_the_annotated_genomic_feature_by_RNAseq_alignments" && $NF~"with_support_for_all_annotated_introns" && $NF~"[0-9]+_Proteins") || ($NF~"Evidence-Data=")' | awk '$5=="mammals" {print$3}' > a1_mammals_transcripts
	for i in $(cat a1_mammals_transcripts)
	do
	awk -v RS=">" -v a="$i" '{if($0~a) print">"$0}' ../ncbi_dataset/data/cds.fna | /media/aswin/programs/myfasta -comb | sed "/^>/ s/.*/>ncbi_A1_mammals_$i/g" | sed '/^>/! s/[a-z]/\U&/g' 
	done | awk NF > a1_mammals_transcripts.fa
	unset i
	#Filter based on median length
	median=$(/media/aswin/programs/myfasta -l a1_mammals_transcripts.fa | awk '{print$NF}' | ministat -n | awk 'END{print$5}')
	seqtk subseq a1_mammals_transcripts.fa <(/media/aswin/programs/myfasta -l a1_mammals_transcripts.fa | awk -v a="$median" '{print$0,$2-a}' | tr -d "-" | sort -k3,3n | head -40 | awk '{print$1}') > A1_mammals_transcripts_filtered.fa
	unset median

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Seperately collect A1 from other non-mammlian tetrapods except birds


##########################################################################################################################################################################################################################################################################################################################
#4. Collect A1 from Sauropsids (Since annotations in other clades other than mammals are not cery reliable)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#4.1. Collect from Literature

	mkdir /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/publications
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/publications
	scp neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/outgroup/reptilian_APOBEC1_from_Wang_paper/reptilian_APOBEC1_from_Wang_paper.txt .
	time awk '{print$1}' reptilian_APOBEC1_from_Wang_paper.txt | xargs -n1 sh -c 'echo "@@"$0; efetch -db protein -id $0 -format docsum' > reptilian_APOBEC1_from_Wang_paper.txt_metadata
	
	scp -r neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/outgroup/Analysis_of_Reptilian_APOBEC1_2010_Dec .
	scp -r neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/outgroup/Diversification_of_AID_APOBEC_like_deaminases_in_metazoa_2018_March .
	scp -r neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/outgroup/A_Long_Running_Arms_Race_between_APOBEC1_Jan_2023 .
	scp -r neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/outgroup/A1_like_from_all_birds.fa .
	scp -r neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/outgroup/A1_from_all_mammals.fa .

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#4.2. Collect from ensembl annotation

cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/ensembl/apobec1_sauropsids
#Manually Download the ortholog table: Go to Ensemble -> Search "APOBEC1" in "human" -> Choose "orthologues" in left menu -> Check "sauropsida" group -> Hover over the right most option in then section "Selected orthologues" -> Click "Download what you see"
sed 's/"//g' A1_sauropsids_ensembl_ortholog_table_raw | sed -e 's/ View Gene Tree//g' -e 's/ Compare Regions (.* View Sequence Alignments//g' | tr "," "\t" | awk -F "\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | sed 's/[ ]\+/ /g' | tr " " "_" | column -t > A1_sauropsids_ensembl_ortholog_table
#Manually make a file named "groups", specifying what taxonomic groups the gene belongs to

#Download sequences
time awk '{print$3}' A1_sauropsids_ensembl_ortholog_table | grep -v "Orthologue" | sed 's/_(/ /g' | awk '{if($1!~"^ENS") print$NF; else print$1}' OFS="\t" | tr -d ")" | xargs -n1 sh -c '/media/aswin/programs/ensembl_exon_wise $0 -c -ir'

#Manually check annotations for Birds & distinguish A1 & A1-like based on synteny (order: NANOG A1-like A1 AICDA) & save it in a file "annotations"
	#Based on manually checked annotations, few birds whose A1 sequence was downloaded but also had A1-like sequence available, and was not downloaded since it was not included in the first ortholog table
	#All species with A1 also had A1-like except for Ostrich & duck
#Manually download few more annotations
	/media/aswin/programs/ensembl_exon_wise ENSSCAG00000015283 -c -ir
	/media/aswin/programs/ensembl_exon_wise ENSABRG00000013573 -c -ir
	/media/aswin/programs/ensembl_exon_wise ENSGFOG00000010992 -c -ir 

#Quick check ORF
ls *.fa | xargs -n1 bash -c 'echo ">"$0; transeq <(grep -v ">" $0) --auto --stdout | grep -v ">" | paste -s -d "" | grep "^M.*\*$" --color=always | grep -v "[A-Z]\*[A-Z]"'

#Collect sequences & create summary table
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/ensembl/apobec1_sauropsids
	rm A1_A1_like_Sauropsids_ensembl.fa
	for i in $(find . -name "*_exon_wise.fa")
	do
	i1=$(echo $i | sed 's/_ENS.*//g' | cut -f2 -d "/")
	i2=$(grep "$i1" groups | awk '{print$2}')
	i3=$(echo $i | grep "ENS[^_]\+" -o)
	i4=$(grep "$i3" annotations | awk '{print$4}')
	if [[ $i4 == "" ]]; then i41="A1_paralog"; else i41=$(echo $i4); fi
	i5=$(transeq <(grep -v ">" $i) --auto --stdout | grep -v ">" | paste -s -d "" | grep "^M.*\*$" --color=always | grep -v "[A-Z]\*[A-Z]")
	if [[ $i5 == "" ]]
	then
	i51="disrupted"
	else
	i51="intact"
	echo ">ensembl_"$i41"_"$i2"_"$i3 >> A1_A1_like_Sauropsids_ensembl.fa
	grep -v ">" $i | paste -s -d "" | sed '/^>/! s/[a-z]/\U&/g' >> A1_A1_like_Sauropsids_ensembl.fa
	fi
	echo $i1 $i3 $i2 $i41 $i51
	unset i1 i2 i3 i4 i41 i5 i51
	done | column -t > ortholog_summary

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#4.3. Collect A1 from orthoDB database  

mkdir /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/orthoDB
cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/orthoDB

#NOTE: This download single sequence per gene ID & does not include all isoforms of a gene
#Method: 
	#Go to NCBI & get Gene-ID of a validated gene such as Human Gene-ID: 339
	#Go to "https://www.orthodb.org" -> Choose search string type is NCBI ID -> Paste 339 -> Click "Get Ortholog Groups" link in the results -> Choose taxonomy level eg: Click tetrapoda -> Then choose Sauropsids -> Click "View CDS fasta" -> Copy sequence & paste in terminal file

#Download using API
#Download seperately for groups: #Crocodylia: 1294634, Testudinata: 2841271 
	curl 'https://data.orthodb.org/current/fasta?id=47241at8457&seqtype=cds&species=1294634' -L -o A1_A1_like_Crocodylia.fa
	sed 's/.*organism_name":"\([^,"]\+\).*pub_gene_id":"\([^,"]\+\).*/>orthodb_A1_paralog_Crocodylia_\1_\2/g' A1_A1_like_Crocodylia.fa | sed '/^>/ s/ /_/g' | awk NF > A1_A1_like_Crocodylia_renamed.fa
	curl 'https://data.orthodb.org/current/fasta?id=47241at8457&seqtype=protein&species=1294634' -L -o A1_A1_like_Crocodylia.aa
	sed 's/.*organism_name":"\([^,"]\+\).*pub_gene_id":"\([^,"]\+\).*/>orthodb_A1_paralog_Crocodylia_\1_\2/g' A1_A1_like_Crocodylia.aa | sed '/^>/ s/ /_/g' | awk NF > A1_A1_like_Crocodylia_renamed.aa
	curl 'https://data.orthodb.org/current/tab?id=47241at8457&species=1294634' -L -o A1_A1_like_Crocodylia_annotations.tsv
	
	curl 'https://data.orthodb.org/current/fasta?id=47241at8457&seqtype=cds&species=2841271' -L -o A1_A1_like_Testudines.fa
	sed 's/.*organism_name":"\([^,"]\+\).*pub_gene_id":"\([^,"]\+\).*/>orthodb_A1_paralog_Testudines_\1_\2/g' A1_A1_like_Testudines.fa | sed '/^>/ s/ /_/g' | awk NF > A1_A1_like_Testudines_renamed.fa
	curl 'https://data.orthodb.org/current/fasta?id=47241at8457&seqtype=protein&species=2841271' -L -o A1_A1_like_Testudines.aa
	sed 's/.*organism_name":"\([^,"]\+\).*pub_gene_id":"\([^,"]\+\).*/>orthodb_A1_paralog_Testudines_\1_\2/g' A1_A1_like_Testudines.aa | sed '/^>/ s/ /_/g' | awk NF > A1_A1_like_Testudines_renamed.aa
	curl 'https://data.orthodb.org/current/tab?id=47241at8457&species=2841271' -L -o A1_A1_like_Testudines_annotations.tsv
	
	curl 'https://data.orthodb.org/current/fasta?id=47241at8457&seqtype=cds&species=8509' -L -o A1_A1_like_Squamata.fa
	sed 's/.*organism_name":"\([^,"]\+\).*pub_gene_id":"\([^,"]\+\).*/>orthodb_A1_paralog_Squamata_\1_\2/g' A1_A1_like_Squamata.fa | sed '/^>/ s/ /_/g' | awk NF > A1_A1_like_Squamata_renamed.fa
	curl 'https://data.orthodb.org/current/fasta?id=47241at8457&seqtype=protein&species=8509' -L -o A1_A1_like_Squamata.aa
	sed 's/.*organism_name":"\([^,"]\+\).*pub_gene_id":"\([^,"]\+\).*/>orthodb_A1_paralog_Squamata_\1_\2/g' A1_A1_like_Squamata.aa | sed '/^>/ s/ /_/g' | awk NF > A1_A1_like_Squamata_renamed.aa
	curl 'https://data.orthodb.org/current/tab?id=47241at8457&species=8509' -L -o A1_A1_like_Squamata_annotations.tsv
	
	curl 'https://data.orthodb.org/current/fasta?id=47241at8457&seqtype=cds&species=8505' -L -o A1_A1_like_Sphenodontia.fa
	curl 'https://data.orthodb.org/current/fasta?id=47241at8457&seqtype=protein&species=8505' -L -o A1_A1_like_Sphenodontia.aa

#Collect sequences
	cat A1_A1_like_Crocodylia_renamed.fa A1_A1_like_Testudines_renamed.fa A1_A1_like_Squamata_renamed.fa > A1_A1_like_Sauropsids_except_Birds_orthodb.fa
	cat A1_A1_like_Crocodylia_renamed.aa A1_A1_like_Testudines_renamed.aa A1_A1_like_Squamata_renamed.aa > A1_A1_like_Sauropsids_except_Birds_orthodb.aa

#NOTE:
	#The orthologs groups have specific ID's eg: "47241at8457" for Sauropsida level which doesn't change irrespective of what input NCBI ID is used to initially search database, which means using human or duck or turtle A1 as the initial, all will give the same orthologs group  
	#The most remote homologs are also shown in this database: eg: For APOBEC1 66 ortholgs from Cnidaria, Arthropoda & Lophotrochozoa whose annotation includes APOBECs, CMP deaminases, & uncharacterized proteins

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#4.4. Collect from NCBI annotation

	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download
	cat /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/apobec1/refseq_select/sequence_similarity_and_supporting_evidence | awk '$5=="birds" {print$1,$2,$3}' | column -t > a1_birds_datasets
	time awk '{print$3}' a1_birds_datasets | xargs -n1 sh -c 'echo "@@"$0; efetch -db nuccore -id $0 -format docsum' > a1_birds_datasets_metadata
	
	#Download from annotation
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/annotated
	time datasets download gene gene-id --inputfile ids --include gene,rna,protein,cds,5p-utr,3p-utr,product-report --filename birds.zip
	time datasets download gene gene-id --ortholog all --inputfile ids --include gene,rna,protein,cds,5p-utr,3p-utr,product-report --filename all_orthologs.zip
	unzip birds.zip -d birds
	unzip all_orthologs.zip -d all_orthologs
	dataformat tsv gene-product --inputfile birds/ncbi_dataset/data/product_report.jsonl | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | column -t > birds/gene_product_metadata
	dataformat tsv gene-product --inputfile all_orthologs/ncbi_dataset/data/product_report.jsonl | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | column -t > all_orthologs/gene_product_metadata

#4.4.1. QC of annotated A1-like from sauropsids

cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/annotated
for x in birds all_orthologs
do
cd $x
echo ">"$x
#Get genbank data
	time grep ">" ncbi_dataset/data/cds.fna | cut -f1 -d ":" | tr -d ">" | xargs -n1 sh -c 'echo "@@"$0; efetch -db nuccore -id "$0" -format gb -mode xml' > genbank
#Print summary
	for i in $(grep ">" ncbi_dataset/data/cds.fna | cut -f1 -d ":" | tr -d ">")
	do
	#i0=$(awk -v RS="@@" -v a="$i" '{if($0~a) print$0}' genbank_filtered_cds_pairwise_seq_similiarity_with_refseq_select | grep "GBSeq_taxonomy" | cut -f2 -d ">" | cut -f1 -d "<" | egrep "Coelacanthiformes|Dipnomorpha|Chondrichthyes|Cyclostomata|Actinopterygii|Mammalia|Crocodylia|Testudines|Aves|Sphenodontia|Squamata|Amphibia" -iwo | sort -u | sed -e 's/Coelacanthiformes/Coelacanth/g' -e 's/Dipnomorpha/Lungfishes/g')
	i0=$(awk -v RS="@@" -v a="$i" '{if($0~a) print$0}' genbank | grep "GBSeq_taxonomy" | cut -f2 -d ">" | cut -f1 -d "<" | egrep "Coelacanthiformes|Dipnomorpha|Chondrichthyes|Cyclostomata|Actinopterygii|Mammalia|Crocodylia|Testudines|Aves|Sphenodontia|Squamata|Amphibia" -iwo | sort -u | sed -e 's/Coelacanthiformes/Coelacanthes/g' -e 's/Dipnomorpha/Lungfishes/g' -e 's/Chondrichthyes/Cartilaginous_fishes/g' -e 's/Cyclostomata/Jawless_vertebrates/g' -e 's/Actinopterygii/Ray_finned_fishes/g')
	i1=$(awk -v RS="@@" -v a="$i" '{if($0~a) print$0}' genbank | grep "GBSeq_accession-version" | cut -f2 -d ">" | cut -f1 -d "<")
	i2=$(awk -v RS="@@" -v a="$i" '{if($0~a) print$0}' genbank | grep "Supporting" | cut -f2 -d ">" | cut -f1 -d "<" | sed 's/^.*Supporting/Supporting/g' | sed 's/Supporting evidence includes similarity to: //g' | sed 's/ /_/g')
	i3=$(awk -v RS="@@" -v a="$i" '{if($0~a) print$0}' genbank | grep "Evidence-Data-START" | cut -f2 -d ">" | cut -f1 -d "<" | tr ";" "\n" | awk NF | awk '/##Evidence\-Data\-START##/,/##RefSeq-Attributes-END##/' \
	| grep -v "END##" | sed 's/-START##//g' | sed '/#/! s/$/|/g' | awk -v RS="##" '{$1=$1}1' | awk NF | sed 's/ /=>/1' | sed 's/ /_/g' | paste -s -d "|")
	if [[ "$i2" == "" ]] && [[ "$i3" == "" ]]; then i2="-"; else :; fi
	i4=$(grep ">" ncbi_dataset/data/cds.fna | grep "$i" | grep "GeneID=[0-9]\+" -o | cut -f2 -d "=")
	i5=$(cat gene_product_metadata | awk 'NR==1{for(i=1;i<=NF;i++) {f[$i]=i}} {print$(f["Transcript_Transcript_Length"]),$(f["Transcript_Protein_Length"]),$(f["Transcript_CDS_Sequence_Accession"]),$(f["Taxonomic_Name"]),$(f["Transcript_Protein_Accession"])}' | grep "$i" | awk '{print$5,$4,$1,$2}' | sort -u)
	i6=$(cat gene_product_metadata | awk 'NR==1{for(i=1;i<=NF;i++) {f[$i]=i}} {print$(f["Transcript_CDS_Sequence_Accession"])}' | grep "$i" | wc -l)
	i7=$(awk -v RS="@@" -v a="$i" '{if($0~a) print$0}' genbank | grep "GBSeq_definition" | cut -f2 -d ">" | cut -f1 -d "<" | tr " " "_")
	echo $i $i4 $i0 $i1 $i5 $i6 $i7 $i2 "$i3"
	unset i0 i1 i2 i3 a i4 i5 i6 i7
	done | sed '1i Input_accession Gene_ID Group Fetched_accession Protein_accession Sci_Name Transcript_length Protein_length No_Exons Definition Supporting_evidence_includes_similarity_to' | column -t > supporting_evidence
#gene metadata
	dataformat tsv gene --inputfile ncbi_dataset/data/data_report.jsonl  | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' > gene_metadata
#gene product metadata
	dataformat tsv gene-product --inputfile ncbi_dataset/data/product_report.jsonl | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | column -t > gene_product_metadata
#Get genome accession associated with gene IDs
	for i in $(awk '$3!="Mammalia"' supporting_evidence | awk '{print$2}' | sort -u | grep -v "Gene_ID")
	do
	i1=$(cat gene_metadata | awk 'NR==1{for(i=1;i<=NF;i++) {f[$i]=i}} {print$(f["NCBI_GeneID"]),$(f["Annotation_Assembly_Accession"])}' | grep "$i")
	i2=$(efetch -db gene -id "$i" -format docsum | grep "AssemblyAccVer" | cut -f2 -d ">" | cut -f1 -d "<" | awk '!a[$o]++')
	echo $i $i1 $i2
	done | column -t > gene_ID_associated_genome_accessions
#Download annotations (5m49.823s)
	mkdir annotations
	cd annotations
	for i in $(awk '{if($3=="No_Assembly") print$4; else print$3}' ../gene_ID_associated_genome_accessions)
	do
	echo ">"$i
	datasets download genome accession "$i" --include gff3 --filename $i".zip"
	unzip $i".zip" -d $i
	done
#Get synteny info
	while read -r i
	do
	i1=$(echo $i | awk '{print$1}')
	i2=$(awk -v a="$i1" '$2==a' ../gene_ID_associated_genome_accessions | awk '{if($3=="No_Assembly") print$4; else print$3}' )
	i3=$(echo $i1 | sed 's/^/GeneID:/g')
	#i4=$(awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -C4 | awk '{print$NF}' | cut -f1 -d ";" | sed 's/ID=gene-//g' | paste -s -d " ")
	#i4=$(awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -C4 | tr " " "_" | awk '{print$NF}' | awk -F ";" '{if($1~"LOC") print$0; else print $1}' | awk -F ";" '{for (i=1;i<=NF;i++) if($i~"ID=gene-" || $i~"description") printf$i; print ""}' | sed 's/ID=gene-//g' | sed -e '/^LOC/ s/description=/(/g' -e '/^LOC/ s/$/)/g' | paste -s -d " ")
	loc=$(awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" | awk '{print$1}')
	ung=$(awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -B4 | awk -v a="$loc" '$1==a' | sed '$d' | wc -l)
	if [[ "$ung" == 0 ]]; then echo -e "-\n-\n-\n-" > kala1
	elif [[ "$ung" == 1 ]]; then awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -B4 | awk -v a="$loc" '$1==a' | sed '$d' | sed '1i -\n-\n-' > kala1
	elif [[ "$ung" == 2 ]]; then awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -B4 | awk -v a="$loc" '$1==a' | sed '$d' | sed '1i -\n-' > kala1
	elif [[ "$ung" == 3 ]]; then awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -B4 | awk -v a="$loc" '$1==a' | sed '$d' | sed '1i -' > kala1
	elif [[ "$ung" == 4 ]]; then awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -B4 | awk -v a="$loc" '$1==a' | sed '$d' > kala1
	else :
	fi
	awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" > kala2
	dng=$(awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -A4 | awk -v a="$loc" '$1==a' | sed '1d' | wc -l)
	if [[ "$dng" == 0 ]]; then echo -e "-\n-\n-\n-" > kala3
	elif [[ "$dng" == 1 ]]; then awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -A4 | awk -v a="$loc" '$1==a' | sed '1d' | sed -e '$a -\n-\n-' > kala3
	elif [[ "$dng" == 2 ]]; then awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -A4 | awk -v a="$loc" '$1==a' | sed '1d' | sed -e '$a -\n-' > kala3
	elif [[ "$dng" == 3 ]]; then awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -A4 | awk -v a="$loc" '$1==a' | sed '1d' | sed -e '$a -' > kala3
	elif [[ "$dng" == 4 ]]; then awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -A4 | awk -v a="$loc" '$1==a' | sed '1d' > kala3
	else :
	fi
	cat kala1 kala2 kala3 | tr " " "_" | awk '{print$NF}' | awk -F ";" '{if($1~"LOC[0-9]+") print$0; else print $1}' | awk -F ";" '{for (i=1;i<=NF;i++) if($i~"ID=gene-" || $i~"description" || (NF==1 && $NF=="-")) printf$i; print ""}' | sed 's/ID=gene-//g' | sed -e '/^LOC[0-9]\+/ s/description=/ /g' > kala
	#awk '$3=="gene" || $3=="pseudogene"' $i2/ncbi_dataset/data/$i2/genomic.gff | grep "$i3" -C4 | awk -v a="$loc" '$1==a' | tr " " "_" | awk '{print$NF}' | awk -F ";" '{if($1~"LOC[0-9]+") print$0; else print $1}' | awk -F ";" '{for (i=1;i<=NF;i++) if($i~"ID=gene-" || $i~"description") printf$i; print ""}' | sed 's/ID=gene-//g' | sed -e '/^LOC[0-9]\+/ s/description=/ /g' > kala
	i4=$(while read j
	do
	j1=$(echo $j | awk '{print NF}')
	if [[ "$j1" < 2 ]] && [[ $(echo $j | grep "^LOC[0-9]\+") == "" ]]
	then :
	else
	j2=$(grep "gene=$j;" $i2/ncbi_dataset/data/$i2/genomic.gff | awk '$3=="exon"' | awk -F ";" '{for (i=1;i<=NF;i++) if($i~"product") printf$i; print ""}' | sed 's/product=//g' | sed 's/ transcript variant X[0-9]\+//g' | sort -u | tr " " "_")
	fi
	echo $j $j2 
	unset j1 j2
	done < kala | awk '{if($2=="") print$0; else print$1"("$2")"}' | paste -s -d " ")
	i5=$(awk -v a="$i1" '$2==a' ../supporting_evidence | awk '{print$3,$6}' | sort | uniq -c | awk '{print$2,$3,$1}')
	echo $i5 $i1 $i4
	unset i1 i2 i3 i4 i5 loc ung dng a
	rm kala1 kala2 kala3 kala
	done < ../gene_ID_associated_genome_accessions > synteny
#Replace LOCID/gene descriptions with Gene symbols if possible (get symbols manually) 
	sed -e 's/[^ ]\+alpha-2-macroglobulin-like_protein_1[^ ]\+/A2ML1/g' -e 's/[^ ]\+NANOG[^ ]\+/NANOG/g' -e 's/[^ ]\+solute_carrier_family_2.*_facilitated_glucose_transporter_member_3[^ ]\+/SLC2A3/g' -e 's/[^ ]\+C->U-editing_enzyme_APOBEC-1-like[^ ]\+/APOBEC1-like/g' \
	-e 's/[^ ]\+microfibrillar-associated_protein_5[^ ]\+/MFAP5/g' -e 's/[^ ]\+single-stranded_DNA_cytosine_deaminase[^ ]\+/AICDA/g' -e 's/uncharacterized_LOC[0-9]\+[^ ]\+/Unchara)/g' -e 's/[^ ]\+C->U-editing_enzyme_APOBEC-1[^ ]\+/APOBEC1/g' -e 's/[^ ]\+taste_receptor_type_2_member_40-like[^ ]\+/TAS2R40_like/g' \
	-e 's/[^ ]\+C-type_lectin_domain_family_2_member_B-like[^ ]\+/CLEC2B_like/g' -e 's/[^ ]\+zinc_finger_protein_239-like%2C[^ ]\+/ZNF239_like/g' -e 's/[^ ]\+trafficking_protein_particle_complex_subunit_8-like[^ ]\+/TRAPPC8_like/g' -e 's/[^ ]\+adhesion_G_protein-coupled_receptor_L1-like[^ ]\+/ADGRL1_like/g' \
	-e 's/[^ ]\+major_urinary_protein-like[^ ]\+/MUP_like/g' -e 's/[^ ]\+prostaglandin-H2_D-isomerase-like[^ ]\+/PTGDS_like/g' -e 's/[^ ]\+bestrophin-3-like[^ ]\+/BEST3_like/g' -e 's/[^ ]\+cerebellin-2-like[^ ]\+/CBLN2_like/g' -e 's/[^ ]\+dmX-like_protein_2[^ ]\+/DMXL2/g' \
	-e 's/[^ ]\+C-type_lectin_domain_family_5_member_A-like[^ ]\+/CLEC5A_like/g' -e 's/[^ ]\+nitric_oxide_synthase%2C_brain-like[^ ]\+/nNOS/g' -e 's/[^ ]\+mitotic_spindle_assembly_checkpoint_protein_MAD1-like[^ ]\+/MAD1L/g' -e 's/[^ ]\+beta-citrylglutamate_synthase_B[^ ]\+/RIMKLB/g' \
	-e 's/[^ ]\+taste_receptor_type_2_member_39-like[^ ]\+/TAS2R39_like/g' -e 's/[^ ]\+putative_DBH-like_monooxygenase_protein_2[^ ]\+/MOXD2P/g' -e 's/[^ ]\+myelin_transcription_factor_1-like_protein[^ ]\+/MYT1L/g' -e 's/[^ ]\+homeobox_protein_vent1B-like[^ ]\+/VENT1/g' \
	-e 's/[^ ]\+alpha-(1%2C3)-fucosyltransferase_10-like[^ ]\+/FUT10_like/g' -e 's/[^ ]\+solute_carrier_family_22_member_1_pseudogene[^ ]\+/SLC2A3_pseudo/g' -e 's/[^ ]\+alpha-N-acetylneuraminide_alpha-2%2C8-sialyltransferase_pseudogene[^ ]\+/ST8SIA1/g' \
	-e 's/[^ ]\+N-acetylated-alpha-linked_acidic_dipeptidase_2-like[^ ]\+/NAALAD2_like/g' -e 's/[^ ]\+putative_N-acetylated-alpha-linked_acidic_dipeptidase[^ ]\+/NAALAD2/g' -e 's/[^ ]\+adenosylhomocysteinase_B-like[^ ]\+/AHCYB_like/g' -e 's/[^ ]\+ephrin_type-A_receptor_6-like[^ ]\+/EPHA6_like/g' \
	-e 's/[^ ]\+glycine-rich_cell_wall_structural_protein_1-like[^ ]\+/GRP-1_like/g' -e 's/[^ ]\+neuroligin-4%2C_X-linked[^ ]\+/NLGN4X/g' -e 's/[^ ]\+zinc_finger_protein_3_homolog[^ ]\+/ZNF3/g' -e 's/[^ ]\+ras-related_protein_Rab-40C[^ ]\+/RAB40C/g' -e 's/[^ ]\+sal-like_protein_3[^ ]\+/SALL3/g' synteny \
	| awk '{if($7=="AICDA" || $8=="AICDA") print$0; else print$1,$2,$3,$4,$13,$12,$11,$10,$9,$8,$7,$6,$5}' | column -t > synteny_renamed
unset i
cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/annotated
done

#Collect all sequences 

	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/annotated/birds
	cat ncbi_dataset/data/cds.fna | /media/aswin/programs/myfasta -comb | sed '/^>/! s/[a-z]/\U&/g' > A1_like_Birds_annotated.fa
	for i in $(grep ">" A1_like_Birds_annotated.fa | cut -f1 -d ":" | tr -d ">")
	do
	i1=$(grep "XM_039562095.1" supporting_evidence | awk '{print$3}' | sed 's/Aves/Birds/g')
	#Add name "A1_like" meaning annotated A1-like gene
	sed "/^>/ s/.*$i.*/>ncbi_A1_like_$i1\_$i/g" A1_like_Birds_annotated.fa -i
	unset i1
	done

#####################################################################################################################################################################################################################################################################################################################
#5. Collect A1-like locally from bird genomes

#NOTE: Collect genes with less complications, don't need to get A1-like from all bird genomes.

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#5.1. Extract exon-wise sequence of annotated A1-like from birds

#Manually check if the A1-like is located between A1 & NANOG based on synteny) 
#Extract exon_wise sequences including annotations with single exon (since it is ok to have an extra query and can capture if such exon structure is what is present in some birds)
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/annotated/birds
	mkdir exon_wise
	time for t in $(awk '{print$11}' gene_product_metadata | egrep -v "\-|Transcript_CDS_Sequence_Accession" |  sort -V | uniq)
	do
	echo ">"$t
	id=$(grep $t supporting_evidence | awk '{print$1,$6,$3}' OFS="_")
	utr5=$(awk -v RS=">" -v a="$t" '$0~a {print">"$0}' ncbi_dataset/data/5p_utr.fna | /media/aswin/programs/myfasta -l | awk '{print$NF}')
	utr3=$(awk -v RS=">" -v a="$t" '$0~a {print">"$0}' ncbi_dataset/data/3p_utr.fna | /media/aswin/programs/myfasta -l | awk '{print$NF}')
	awk 'NR==1{for(i=1;i<=NF;i++) {f[$i]=i}} {print$(f["Transcript_Genomic_Exons_Start"]),$(f["Transcript_Genomic_Exons_Stop"]),$(f["Transcript_CDS_Sequence_Accession"])}' gene_product_metadata | grep "$t" | awk '{print$3,$2-$1+1}' | awk -v a="$utr5" 'NR==1{print$2-a; next} {print$2}' \
	| awk -v a="$utr3" 'NR>1 {print last} {last=$0} END{print$1-a}' | awk '{if(p~"-") print$1+p; else print$1; p=$1}' | grep -v "\-" | awk '{a+=$1} {print$1,a}' | awk '{print$1,p,$2; p=$2+1}' | awk 'NR==1{print$1,1,$2; next} {print$1,$2,$3}' | awk '{print NR,$2,$3}' > kala.bed 
	awk -v RS=">" -v a="$t" '$0~a {print">"$0}' ncbi_dataset/data/cds.fna | /media/aswin/programs/myfasta -comb > kalat.fa
	cat kala.bed | xargs -n3 sh -c 'echo ">exon_"$0; grep -v ">" kalat.fa | cut -c$1-$2' > exon_wise/$id"_exon_wise.fa"
	unset utr5 utr3 id
	rm kala.bed kalat.fa
	done

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#5.2. CGBLAST

#In neo
	cd ~/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/annotated_A1_like_sequences
	scp -r ceglab25@172.30.1.131:/media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/annotated/birds/exon_wise .
	
	mkdir /home/neo/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/cgblast/curated_exon_wise_A1_like_birds
	cd /home/neo/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/cgblast/curated_exon_wise_A1_like_birds

#569m55.945s
	time while read i
	do
	i1=`echo $i | awk '{print$1}'`
	i2=`echo $i | awk '{print$2}'`
	echo ">Subject: "$i1
	mkdir $i2/aswin/APOBEC1/curated_exon_wise_A1_like_birds
	cd $i2/aswin/APOBEC1/curated_exon_wise_A1_like_birds
	gn=$(find ../../../genome/ -maxdepth 1 -name "GC*.fna" | xargs readlink -f)
	gff=$(find ../../../genome/ -maxdepth 1 -name "GC*.gff" | xargs readlink -f 2> /dev/null)
	#Run seperately for each exon_wise sequences 
	time for s in $(find /home/neo/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/annotated_A1_like_sequences/exon_wise -name "*_Aves_exon_wise.fa")
	do
	s1=$(echo $s | awk -F "/" '{print$NF}' | sed 's/\.fa//g')
	s2=$(echo $s | awk -F "/" '{print$NF}')
	echo ">Query: "$s1
	mkdir $s1
	cd $s1
	cp $s .
	#Run gblast
	if [[ -z $gff ]]
	then
	gblast_short $gn $s2 -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank
	else
	gblast_short $gn $s2 -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank $gff
	fi
	unset s1 s2
	cd ../
	done
	unset i1 i2 gn gff s
	cd /home/neo/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/cgblast/curated_exon_wise_A1_like_birds
	done < /home/neo/bird_db1/aswin/database_details/all_genome_paths

#####################################################################################################################################################################################################################################################################################################################
#6. Manually curate A1-like from birds

#6.1. Check for complete ORFs from genome blast

	cd /home/neo/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/cgblast/curated_exon_wise_A1_like_birds
	time while read i
	do
	i1=`echo $i | awk '{print$1}'`
	i2=`echo $i | awk '{print$2}'`
	cd $i2/aswin/APOBEC1/curated_exon_wise_A1_like_birds
	i3=$(find . -maxdepth 3 -name "gblast_edited_consensus_orf.fa" | xargs -n1 bash -c 'paste <(echo $0 | awk -F "/" "{print\$2}") <(grep -v ">" $0 | paste -s -d "" | grep -v ">" | paste -s -d "" | grep "^M.*\*$" --color=always | grep -v "[A-Z]\*[A-Z]")' | awk '$2')
	k=$(for j in $(echo "$i3" | awk '{print$1}')
	do
	j1=$(grep ">" $j/*_Aves_exon_wise.fa -c)
	j2=$(grep ">" $j/gblast_edited_consensus.fa -c)
	j31=$(echo "$i3" | grep "$j" | awk '{print$1}')
	j32=$(echo "$i3" | grep "$j" | awk '{print$2}')
	j4=$(grep -v ">" $j/*_Aves_exon_wise.fa | wc | awk '{print$3-$1}')
	j5=$(grep -v ">" $j/gblast_edited_consensus.fa | wc | awk '{print$3-$1}')
	j6=$(awk '$1!~"dup" {$1=""; print}' $j/Synteny | sort | uniq -c | column -t | sed 's/[ ]\+ /_/g')
	if [[ "$j6" == "" ]]; then j6="-"; else :; fi
	if [ $(echo "$j6" | wc -l) -gt 1 ]; then j61=$(echo "$j6" | awk -F "\t" 'NR>1{print"- - - - - "$0; next} {print$0}'); else j61=$(echo "$j6"); fi
	j7=$(awk 'NR==1{$1=""; sub(/^ /, "")} {last=$0} NR>1{print prev} {prev=last} END{$NF=""; print}' $j/splice_sites | paste -s -d "\t" | xargs -n2 sh -c 'echo $0 $1' | tr "[a-z]" "[A-Z]" | awk '{if ($1 == "GT" && $2 == "AG") print "intact"; else print "disrupted"}' | sort | uniq -c | awk '{print$1"_"$2}')
	if [ $(echo "$j6" | wc -l) -eq 1 ] && [ $(echo "$j7" | wc -l) -gt 1 ]; then j71=$(echo "$j7" | awk -F "\t" 'NR>1{print"- - - - - - "$0; next} {print$0}')
	elif [ $(calc $(echo "$j7" | wc -l) - $(echo "$j6" | wc -l)) -gt 0 ]; then j72=$(echo "$j6" | wc -l); j71=$(echo "$j7" | awk -F "\t" -v a="$j72" 'NR>a{print"- - - - - - "$0; next} {print$0}'); else j71=$(echo "$j7"); fi
	j8=$(needle <(transeq <(grep -v ">" $j/*_Aves_exon_wise.fa) --trim --auto --stdout | sed '/^>/ s/.*/>Query/g') <(transeq <(grep -v ">" $j/gblast_edited_consensus.fa) --trim --auto --stdout | sed '/^>/ s/.*/>Subject/g') --auto --stdout -awidth 1000 | grep "^Query" -A2 | cut -c22- | sed 's/[ ]\+[0-9]\+$//g' | tr " " "_")
	if [ $(echo "$j6" | wc -l) -eq 1 ] && [ $(echo "$j7" | wc -l) -eq 1 ] && [ $(echo "$j8" | wc -l) -gt 1 ]; then j81=$(echo "$j8" | awk -F "\t" 'NR>1{print"- - - - - - - "$0; next} {print$0}')
	elif [ $(calc $(echo "$j8" | wc -l) - $(echo "$j7" | wc -l)) -gt 0 ]; then j82=$(echo "$j7" | wc -l); j81=$(echo "$j8" | awk -F "\t" -v a="$j82" 'NR>a{print"- - - - - - - "$0; next} {print$0}'); else j81=$(echo "$j8"); fi
	#echo $j1 $j2 $j4 $j5 "$j6" "$j7" "$j8" $j3 | awk '{print$5,$1,$2,$3,$4,$6}'
	paste <(echo $j31 $j1 $j2 $j4 $j5) <(echo "$j61") <(echo "$j71") <(echo "$j81") <(echo $j32)
	unset j1 j2 j31 j32 j4 j5 j6 j61 j7 j71 j72 j8 j81 j82
	done)
	echo $i1 "$k" | awk 'NR>1{print"-",$0; next} {print$0}'
	unset i1 i2 i3 j k
	cd /home/neo/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/cgblast/curated_exon_wise_A1_like_birds
	done < /home/neo/bird_db1/aswin/database_details/all_genome_paths | sed 's/_Aves_exon_wise//g' | awk '{print$1,$2,$3,$4,$5,$6,$5-$6,$7,$8,$9,$10}' | sed '1i Species Query Query_Exons Subject_Exons Query_Length Subject_Length Length_diff Synteny Splice_sites Query_Subject_Alignment Subject_ORF' | column -t > species_wo_annotation_cgblast_summary
	#done < <(grep -if <(find /home/neo/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/annotated_A1_like_sequences/exon_wise -name "*_Aves_exon_wise.fa" | awk -F "/" '{print$NF}' | cut -f3,4 -d "_") /home/neo/bird_db1/aswin/database_details/all_genome_paths -v) | sed 's/_Aves_exon_wise//g' | awk '{print$1,$2,$3,$4,$5,$6,$5-$6,$7,$8,$9,$10}' | sed '1i Species Query Query_Exons Subject_Exons Query_Length Subject_Length Length_diff Synteny Splice_sites Query_Subject_Alignment Subject_ORF' | column -t > species_wo_annotation_cgblast_summary

	#Quick view if bed entries overlap
	readlink -f best_hits ../../best_hits | xargs -n1 bash -c 'n=$(echo $0 | awk -F "/" "{print\$(NF-1)}"); awk -v a="$n" "!/S_start/ {print a,\$1,\$2,\$7,\$8}" $0' | sort -k4,4nr  | column -t
	#Quick check alignment between gblast_edited_consensus & annotated A1-like if available
	needle <(transeq <(grep -v ">" gblast_edited_consensus.fa) --trim --auto --stdout | sed 's/^>.*/>gedit/g') <(transeq ../XM_009950364.1_Leptosomus_discolor_Aves_exon_wise/LOC104402027.fa --trim --auto --stdout) --auto --stdout -awidth 250 
	
	#NOTE: 
	#Some subject species with longer A1-like query lacked complete ORF, but 1 / 2 exon containing queries gave complete ORF. 
	#This might represent a shorter isoform in subject, neverthless it has a unique sequence in contrast to A1, because the hits of these short A1-like in subject is located next to A1 & NANOG as expected & has no overlap.
	#Some exceptions: 
		#Phoenicopterus_ruber - gave complete ORF when used with longer 4 exon queries, but last exon of these hits are too distant likely after many genes.
		                     #- hence here as well take the shorter version which doesnot overlap with other genes.
		#Merops_nubicus       - Exon 3 missing, but still formed complete ORF
	 	#Nestor_notabilis     - 3 Queries gave hits, with 1, 2 & 3 exon Queries. Annotation was present for A1-like.
	                         #- Annotation sequence is exact match with 1 exon query hit, but 3 exon query hit had extra 54 amni acids! which is quite long, still this consensus is chosen, hence be cautious!!!
		#Phaethon_lepturus    - annotation present with 2 exons, but 1st exon hit from local blast is a bit more downstream than the annotated exon 1
		                     #- But this exon 1 is very short and when whole protein of annotated mRNA and gblast edited consensus protein were pairwise aligned, only 2 amnio acids were extra in gblast edited  with 3 different amni acid mismatches, rest is 100% identical.
		                     #-This seems to not affect much of the analysis such as clustering and phylogeny

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.2. Collect the manually curated sequences 

	cd /home/neo/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/cgblast/curated_exon_wise_A1_like_birds
	rm A1_like_Birds_curated.fa
	while read i
	do
	i1=`echo $i | awk '{print$1}'`
	i2=`echo $i | awk '{print$2}'`
	cd $i2/aswin/APOBEC1/curated_exon_wise_A1_like_birds
	i3=$(find . -name "curated_consensus.fa")
	if [[ $i3 == "" ]]
	then :
	else
	cat $i3 | grep -v ">" | paste -s -d "" | sed "1i >curated_A1_like_Birds_$i1" | sed '/^>/! s/[a-z]/\U&/g' >> ~/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/cgblast/curated_exon_wise_A1_like_birds/A1_like_Birds_curated.fa
	fi
	echo $i1 $i3
	unset i1 i2 i3
	cd /home/neo/bird_db1/aswin/APOBEC1/APOBEC1_duplicate/cgblast/curated_exon_wise_A1_like_birds
	done < /home/neo/bird_db1/aswin/database_details/all_genome_paths

#NOTE: Calypte_anna ORF is intact but has a stop codon just 1 position before last stop.= i,e, 2 stop consecutive stops, hence remove one since it doesn't affect the sequence much.

scp A1_like_Birds_curated.fa ceglab25@172.30.1.131:/media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/annotated/birds 

#####################################################################################################################################################################################################################################################################################################################
#7. Codon alignment & Phylogeny

#7.1. Prepare sequences

#Combine all sequences
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny
	cat /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/aid_a2_a3_a4_filtered.fa /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manually_validated/validated_A1_Birds.fa \
	/media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/apobec1/refseq_select/A1_mammals_transcripts_filtered.fa /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/ensembl/apobec1_sauropsids/A1_A1_like_Sauropsids_ensembl.fa \
	/media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/orthoDB/A1_A1_like_Sauropsids_except_Birds_orthodb.fa /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/annotated/birds/A1_like_Birds_annotated.fa \
	/media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/annotated/birds/A1_like_Birds_curated.fa > AID_APOBEC_cds.fa

#Delete duplicated sequences based on sequences not ID
	seqkit rmdup -s -i < AID_APOBEC_cds.fa -d removed_sequences_AID_APOBEC_cds.fa -D duplicated_sequences_AID_APOBEC_cds | /media/aswin/programs/myfasta -comb > AID_APOBEC_cds_unique.fa

#Delete sequences which doesn't form complete ORF or not divisible by 3
	time for i in $(grep ">" AID_APOBEC_cds_unique.fa | tr -d ">")
	do
	i1=$(transeq <(seqtk subseq AID_APOBEC_cds_unique.fa <(echo $i)) --auto --stdout | grep -v ">" | paste -s -d "" | tr "[a-z]" "[A-Z]" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]")
	if [[ $i1 == "" ]]; then i2="disrupted"; else i2="intact"; fi
	i3=$(seqtk subseq AID_APOBEC_cds_unique.fa <(echo $i) | /media/aswin/programs/myfasta -l | awk '{print$2,$2/3}')
	i4=$(transeq <(seqtk subseq AID_APOBEC_cds_unique.fa <(echo $i)) --auto --stdout | grep -v ">" | paste -s -d "" | tr "[a-z]" "[A-Z]" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]" | grep "[*]\+\*$")
	if [[ $i4 == "" ]]; then i5="no"; else i5="yes"; fi
	i6=$(seqtk subseq AID_APOBEC_cds_unique.fa <(echo $i) | grep -v ">" | tr -d "ATGCatgc")
	if [[ $i6 == "" ]]; then i7="no"; else i7="yes"; fi
	echo $i $i3 $i5 $i7 $i2
	unset i1 i2 i3 i4 i5 i6 i7
	done | sed '1i Sequence Length Divisible_by_3 Consecutive_stops_at_end Anonymous_bases ORF' | column -t > orf_check_AID_APOBEC_cds

#Filter
	seqtk subseq AID_APOBEC_cds_unique.fa <(awk '$4=="no" && $5=="no" && $6=="intact" {print$1}' orf_check_AID_APOBEC_cds) > AID_APOBEC_cds_unique_complete_orfs.fa

#7.2. codon alignment 11m8.273s
	time perl /home/ceglab25/guidance.v2.02/www/Guidance/guidance.pl --program GUIDANCE --seqFile AID_APOBEC_cds_unique_complete_orfs.fa --msaProgram MAFFT --seqType codon --outDir codon_alignment --proc_num 32
	cp codon_alignment/MSA.MAFFT.aln.With_Names AID_APOBEC_cds_unique_complete_orfs.aln

#7.3. Phylogeny
	nohup time /media/aswin/programs/iqtree-2.3.2-Linux-intel/bin/iqtree2 --ufboot 1000 --mem 80G -T AUTO -ntmax 32 -m MFP -s AID_APOBEC_cds_unique_complete_orfs.aln &

#####################################################################################################################################################################################################################################################################################################################
#8. Add missed species & remove wierdly phylogenetically placed species

#Few species are not included after collection sequences based on annotations available in databases
	#Manually add a list of species to include in phylogeny based on previous knowledge "species_to_add_for_phylogeny.fa"
	#Check ORF before adding
	transeq species_to_add_for_phylogeny.fa --auto --stdout | /media/aswin/programs/myfasta -comb | grep -v ">" | tr "[a-z]" "[A-Z]" | grep "^M.*\*$" --color=always| grep -v "[A-Z]\*[A-Z]"
	cat AID_APOBEC_cds_unique_complete_orfs.fa species_to_add_for_phylogeny.fa > AID_APOBEC_cds_unique_complete_orfs_refined.fa

#After doing phylogeny few species are noted as clade to be removed 
	#Manually create a list "species_to_remove_for_phylogeny"
	awk  -iinplace 'FNR==NR { a[$0]; next } /^>/ { f = !(substr($0, 2) in a) }f' species_to_remove_for_phylogeny AID_APOBEC_cds_unique_complete_orfs_refined.fa

#codon alignment & phylogeny
	time perl /home/ceglab25/guidance.v2.02/www/Guidance/guidance.pl --program GUIDANCE --seqFile AID_APOBEC_cds_unique_complete_orfs_refined.fa --msaProgram MAFFT --seqType codon --outDir codon_alignment_refined --proc_num 32
	cp codon_alignment_refined/MSA.MAFFT.aln.With_Names AID_APOBEC_cds_unique_complete_orfs_refined.aln
	nohup time /media/aswin/programs/iqtree-2.3.2-Linux-intel/bin/iqtree2 --ufboot 1000 --mem 80G -T AUTO -ntmax 32 -m MFP -s AID_APOBEC_cds_unique_complete_orfs_refined.aln &

#####################################################################################################################################################################################################################################################################################################################
#9. Find more A1 annotations from marsupials & monotremes & amphibians (these might not be included in datastes package OR was removed as a part of QC) & build phylogeny

#9.1. Download marsupials
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/marsupials
	#Use Ornithorhynchus_anatinus APOBEC1 gene as query to blast against marsupial nr databases
	#download the hit table of blast
	awk '!/#/ {print$2}' Ornithorhynchus_anatinus_APOBEC1_XM_029054257.2_exons_with_UTR_against_marsupials_nr_database | awk NF | grep "[XN]M_" | sort -u > unique_hit_accessions
	#Download exon-wise sequences of hits (26m28.632s)
	time cat unique_hit_accessions | xargs -n1 sh -c 'echo ">"$0; /media/aswin/programs/exon_wise -datasets -tid "$0" -er -s'
	ls *_exon_wise.fa | xargs -n1 sh -c 'echo $0 | sed "s/_exon_wise.fa//g" | sed "s/^/>/g"; grep -v ">" $0 | paste -s -d ""' > marsupials_all_cds.fa

	#Check sequence similarity with monotreme sequence
	for i in $(ls *_exon_wise.fa)
	do
	i1=$(echo $i | sed 's/_exon_wise.fa//g')
	i2=$(needle <(grep -v ">" $i | paste -s -d "" | sed "1i >$i1") <(grep -v ">" ../monotremes/Ornithorhynchus_anatinus_XM_029054255.2_exon_wise.fa | paste -s -d "" | sed '1i >Ornithorhynchus_anatinus') --auto --stdout | egrep "Length|Identity|Similarity|Gaps|Score" | awk '{print$3}' | cut -f1 -d "/" | paste -s -d " ")
	i3=$(grep -v ">" $i | paste -s -d "" | myfasta -l | awk '{print$NF}')
	i4=$(grep -v ">" ../monotremes/Ornithorhynchus_anatinus_XM_029054255.2_exon_wise.fa | paste -s -d "" | myfasta -l | awk '{print$NF}')
	echo $i1 $i3 "Ornithorhynchus_anatinus" $i4 $i2
	unset i1 i2 i3 i4
	done | sed '1i Subject Sub_length Reference Ref_length Length Identity Similarity Gaps Score' | column -t > marsupials_all_cds_similarity_with_Ornithorhynchus_anatinus
	#Choose one isoform / species
	sed 's/_XM_/ XM_/g' marsupials_all_cds_similarity_with_Ornithorhynchus_anatinus | sed 's/_NM_/ NM_/g' | sed 's/[ ]\+/ /g' | grep -v "Subject" | sort -k1,1V -k10,10nr | awk '!a[$1]++' > marsupials_all_cds_similarity_with_Ornithorhynchus_anatinus_selected
	seqtk subseq marsupials_all_cds.fa <(awk '{print$1,$2}' OFS="_" marsupials_all_cds_similarity_with_Ornithorhynchus_anatinus_selected) > marsupials_all_cds_similarity_with_Ornithorhynchus_anatinus_selected.fa
	mafft --auto --reorder --quiet marsupials_all_cds_similarity_with_Ornithorhynchus_anatinus_selected.fa | alv -
	awk -v RS=">" '!/Sarcophilus_harrisii_XM_031940193.1|Petaurus_breviceps_papuanus_XM_069106037.1/ {print">"$0}' marsupials_all_cds_similarity_with_Ornithorhynchus_anatinus_selected.fa | awk NF | grep -v "^>$" | sed '/^>/ s/>/>ncbi_A1_marsupials_/g' > marsupials_selected_cds.fa
	#Check ORF
	for i in $(grep ">" marsupials_selected_cds.fa | tr -d ">")
	do
	i0=$(seqtk subseq marsupials_selected_cds.fa <(echo $i) | myfasta -l | awk '{print$2/3}')
	i1=$(transeq <(seqtk subseq marsupials_selected_cds.fa <(echo $i)) --auto --stdout | grep -v ">" | paste -s -d "" | tr "[a-z]" "[A-Z]" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]" | grep -v '\(*\)\1\{1,\}')
	if [[ $i1 == "" ]]; then i2="disrupted"; else i2="intact"; fi
	echo $i $i0 $i2
	unset i0 i1 i2
	done | column -t
#One sequence had an error during to false info in the 3' utr length, manually correct the sequence in: Monodelphis_domestica_NM_001032982.2_exon_wise.fa & marsupials_selected_cds.fa

#9.2. Download monotremes
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/monotremes
	time awk '{print$NF}' add_these | xargs -n1 sh -c 'echo ">"$0; /media/aswin/programs/exon_wise -datasets -tid "$0" -er -s'
	ls *_exon_wise.fa | xargs -n1 sh -c 'echo $0 | sed "s/_exon_wise.fa//g" | sed "s/^/>/g"; grep -v ">" $0 | paste -s -d ""' > monotremes_all_cds.fa
	#Align & choose single isoform /species
	mafft --auto --reorder --quiet monotremes_all_cds.fa | alv -
	awk -v RS=">" '/Ornithorhynchus_anatinus_XM_029054255.2|Tachyglossus_aculeatus_XM_038767961.1/ {print">"$0}' monotremes_all_cds.fa | awk NF | sed '/^>/ s/>/>ncbi_A1_monotremes_/g' > monotremes_selected_cds.fa
	for i in $(grep ">" monotremes_selected_cds.fa | tr -d ">")
	do
	i0=$(seqtk subseq monotremes_selected_cds.fa <(echo $i) | myfasta -l | awk '{print$2/3}')
	i1=$(transeq <(seqtk subseq monotremes_selected_cds.fa <(echo $i)) --auto --stdout | grep -v ">" | paste -s -d "" | tr "[a-z]" "[A-Z]" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]" | grep -v '\(*\)\1\{1,\}')
	if [[ $i1 == "" ]]; then i2="disrupted"; else i2="intact"; fi
	echo $i $i0 $i2
	unset i0 i1 i2
	done | column -t

#9.3. Download amphibian A1
	#Previous study report A1 detected in a single amphibian species Pleurodeles_waltl 
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/amphibians
	/media/aswin/programs/exon_wise -datasets -tid XM_069201938.1 -er -s
	#Running this XM_069201938.1 as blastn query in ncbi nr database gave one more amphibian hit (rest hits for exon_2 were in birds)
	/media/aswin/programs/exon_wise -datasets -tid XM_069650576.1 -er -s
	#NOTE: RNA expression is present in both species for A1
	scp neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/relaxed_selection/amniotes_as_bg/12_amniotes_as_bg/amniotes_cds_aliview_filtered2.fa .
	mafft --auto --reorder --quiet <(cat amniotes_cds_aliview_filtered2.fa <(grep -v ">" Pleurodeles_waltl_XM_069201938.1_exon_wise.fa | paste -s -d "" | sed '1i >Pleurodeles_waltl') <(grep -v ">" Ambystoma_mexicanum_XM_069650576.1_exon_wise.fa | paste -s -d "" | sed '1i >Ambystoma_mexicanum') <(grep -v ">" Ambystoma_mexicanum_XM_069650577.1_exon_wise.fa | paste -s -d "" | sed '1i >Ambystoma_mexicanum')) | alv -
	ls Pleurodeles_waltl_XM_069201938.1_exon_wise.fa Ambystoma_mexicanum_XM_069650576.1_exon_wise.fa | xargs -n1 sh -c 'echo $0 | sed "s/_exon_wise.fa//g" | sed "s/^/>ncbi_A1_amphibians_/g"; grep -v ">" $0 | paste -s -d ""' > amphibians_selected_cds.fa
	#check ORF
	for i in $(grep ">" amphibians_selected_cds.fa | tr -d ">")
	do
	i0=$(seqtk subseq amphibians_selected_cds.fa <(echo $i) | myfasta -l | awk '{print$2/3}')
	i1=$(transeq <(seqtk subseq amphibians_selected_cds.fa <(echo $i)) --auto --stdout | grep -v ">" | paste -s -d "" | tr "[a-z]" "[A-Z]" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]" | grep -v '\(*\)\1\{1,\}')
	if [[ $i1 == "" ]]; then i2="disrupted"; else i2="intact"; fi
	echo $i $i0 $i2
	unset i0 i1 i2
	done | column -t

#Remove duplicated & combine all selected sequences
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny
	#There was already a marsupial present in the last refined sequence hence remove it & add it from seperate marsupial selected cds
	awk -v RS=">" '!/ncbi_A1_mammals_NM_001032982.2/ {print">"$0}' AID_APOBEC_cds_unique_complete_orfs_refined.fa | grep -v "^>$" | awk NF > AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa
	cat cds_sequences/manual_download/marsupials/marsupials_selected_cds.fa cds_sequences/manual_download/monotremes/monotremes_selected_cds.fa cds_sequences/manual_download/amphibians/amphibians_selected_cds.fa >> AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#9.4. Run phylogeny again with more added species

#codon alignment
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny
	time perl /home/ceglab25/guidance.v2.02/www/Guidance/guidance.pl --program GUIDANCE --seqFile AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa --msaProgram MAFFT --seqType codon --outDir codon_alignment_refined_monotremes_marsupials_amphibians --proc_num 32
	cp codon_alignment_refined_monotremes_marsupials_amphibians/MSA.MAFFT.aln.With_Names AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.aln

#Phylogeny
	nohup time /media/aswin/programs/iqtree-2.3.2-Linux-intel/bin/iqtree2 --ufboot 1000 --mem 80G -T AUTO -ntmax 32 -m MFP -s AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.aln

#####################################################################################################################################################################################################################################################################################################################
#10. Filtering again based on weird phylogeny placement

#When vizualised the tree in phylogeny some sequences are placed in susposious position
	#ncbi_A1_mammals_XM_004869517.2: Heterocephalus glaber (Naked mole rat)
	#In ncbi genome track this shows huge genome gap within intron, hence remove this sequence

cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny
awk -v RS=">" '!/ncbi_A1_mammals_XM_004869517.2/ {print">"$0}' AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa | awk '!/^>$/ && NF' > AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined.fa

#Rerun with removed mammal (10m47.858s)
	time perl /home/ceglab25/guidance.v2.02/www/Guidance/guidance.pl --program GUIDANCE --seqFile AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined.fa --msaProgram MAFFT --seqType codon --outDir codon_alignment_refined_monotremes_marsupials_amphibians_refined --proc_num 32
	cp codon_alignment_refined_monotremes_marsupials_amphibians_refined/MSA.MAFFT.aln.With_Names AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined.aln
	nohup time /media/aswin/programs/iqtree-2.3.2-Linux-intel/bin/iqtree2 --ufboot 1000 --mem 80G -T AUTO -ntmax 32 -m MFP -s AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined.aln

#Rerun with higher bootstrap to check if results vary (102m57.476s)
	time perl /home/ceglab25/guidance.v2.02/www/Guidance/guidance.pl --program GUIDANCE --seqFile AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined.fa --msaProgram MAFFT --seqType codon --outDir codon_alignment_refined_monotremes_marsupials_amphibians_refined_run_2 --bootstraps 1000 --proc_num 32
	cp codon_alignment_refined_monotremes_marsupials_amphibians_refined_run_2/MSA.MAFFT.aln.With_Names AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_run_2.aln
	#99m55.463s
	time /media/aswin/programs/iqtree-2.3.2-Linux-intel/bin/iqtree2 --ufboot 5000 --mem 80G -T AUTO -ntmax 32 -m MFP -s AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_run_2.aln

#Rerun with higher bootstrap & nucleotide alignment instead of codon alignment to check if results vary (1m43.500s)
	time mafft --maxiterate 1000 --localpair --reorder --thread 32 AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined.fa > AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_nt.aln
	#(41m0.840s)
	time /media/aswin/programs/iqtree-2.3.2-Linux-intel/bin/iqtree2 --ufboot 5000 --mem 80G -T AUTO -ntmax 32 -m MFP -s AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_nt.aln


#NOTE:
	#When one single mammal species was remove & ran phylogeny the lung fish A1 was wierdly placed as a seperate branch diverged before A2, A3, & AID instead of as sister brach along with them. 
	#When checked in the phylogeny log it is noticed that the best model model automatically found was different in these 2 runs
	#Articles shows the most consistent phylogeny results were produced by iq-tree2 with "-m MFP" parameter, which I used in this script.
	#If tree inferred is varying despote of this consistency, it might suggest that the input sequence must be wierd in such a way that removal of a single species dramatically changed the placement of gaps in the alignment.
	#Anyway it is recommended to manually inspect the input sequences & their alignment & remove any wierd sequences

#####################################################################################################################################################################################################################################################################################################################
#11. Final QC via Manual Inspection

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#11.1. Manual inpsect sequences before codon alignment

#Align each major clades seperately & manually inspect alignment for weird sequences especially in terms of length & gaps & remove only the most extremely wierd sequences
	mkdir /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/alignment_manual_inspection
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/alignment_manual_inspection
	#Manually create a file major_clades
	for i in $(cat major_clades)
	do
	echo $i
	time mafft --auto --reorder --quiet <(awk -v RS=">" -v a="$i" '$0~a {print">"$0}' ../AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa | grep -v "^>$" | awk NF) > inspect_"$i".aln
	done

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#11.2. Manually check databases for supporting informations & QC for suspicious sequences/clades

#Steps for ensembl: 
	#- Search gene ID (ENSGEVG00005001225)
	#- Check Genome gap, Synteny, Repeats, GC %, EST/mRNA reads, in "Region in Detail" page.
 		#- Make sure the tracks such as Contigs, Ensembl Genes, EST(normal), CpG islands, All repeats(expanded) & %GC are toggled on save this configuration such that these data are visible in the image.
		#- Use "Share this page" for direct link to the page
	#- Check splice sites in the section "Transcript-based displays" under "Sequence" > "Exons"
	#- Check supporting evidence in the "Supporting evidence" section
 		#- The evidence are shown as Type of evidence (protein or intron) & region supporting the evidence
		#- Use "Share this page" for direct link to the page
#Steps for NCBI:
	#- Seatrch Gene ID (127042962)
	#- Toggle on these tracks: Open "Tracks" > "Configure tracks" > "Available tracks" > Sequence Tiling Path, Scaffolds, G+C content, CpG Islands, Repeat features produced using WindowMasker
	#- Open right most icon & click "Link to view" for getting link to the page.

#Manually create a list of sequences to remove & add based on manual curation based on genome gaps, exonic repeats, conserved synteny, splice sites, RNA expression, domains etc ...

#Testudines: Remove & Add
	mkdir /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/Testudines_QC
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/Testudines_QC
	time cat ../../../Testudines_to_add | xargs -n1 sh -c 'echo ">"$0; /media/aswin/programs/exon_wise -datasets -tid "$0" -er -s'
	#Correct if QC shows errors & combine fasta
	for i in $(ls *exon_wise.fa); do j=$(echo $i | awk -F "_" '{print$(NF-3)"_"$(NF-2)}'); grep -v ">" $i | paste -s -d "" | sed "1i >manualncbi_A1_paralog_Testudines_$j"; unset j; done > Testudines_to_add.fa 
	#Check ORF
	for i in $(ls *exon_wise.fa); do j=$(transeq <(grep -v ">" $i | paste -s -d "") --auto --stdout | grep -v ">" | tr "[a-z]" "[A-Z]" | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]" | grep -v '\(*\)\1\{1,\}'); if [[ -z $j ]]; then k="dis"; else k="int"; fi; echo $i $k; unset j k; done | column -t
	#Check for "N"s gaps
	awk -v RS=">" '{$1=$1}1' Testudines_to_add.fa | awk 'tolower($2)~"n" {print$1}'
	awk -v RS=">" '!/XM_038386549.2|XM_039489155.1/ {print">"$0}' Testudines_to_add.fa | grep -v "^>$" | awk NF > Testudines_to_add_refined.fa 

#Mammals: Remove
	mkdir /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/Mammals_QC
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/Mammals_QC
	grep "A1_mammals" /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa | cut -f4,5 -d "_" > ids
	time cat ids | xargs -n1 sh -c 'echo ">"$0; /media/aswin/programs/exon_wise -datasets -tid "$0" -er -s'

#Marsupials: No change
	seqtk subseq AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa <(grep -i "A1_marsupials" AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa | tr -d ">") | myfasta -l

#Monotremes: Remove & Add
	seqtk subseq AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa <(grep -i "A1_monotremes" AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa | tr -d ">") | myfasta -l
	mkdir /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/Monotremes_QC
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/Monotremes_QC
	cat ../../../Monotremes_to_add > ids
	time cat ids | xargs -n1 sh -c 'echo ">"$0; /media/aswin/programs/exon_wise -datasets -tid "$0" -er -s'
	#Check ORF
	for i in $(ls *exon_wise.fa); do j=$(transeq <(grep -v ">" $i | paste -s -d "") --auto --stdout | grep -v ">" | tr "[a-z]" "[A-Z]" | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]" | grep -v '\(*\)\1\{1,\}'); if [[ -z $j ]]; then k="dis"; else k="int"; fi; echo $i $k; unset j k; done | column -t
	for i in $(ls *exon_wise.fa); do j=$(echo $i | awk -F "_" '{print$(NF-3)"_"$(NF-2)}'); grep -v ">" $i | paste -s -d "" | sed "1i >manualncbi_A1_monotremes_$j"; unset j; done > Monotremes_to_add.fa 
	#Check for "N"s gaps
	awk -v RS=">" '{$1=$1}1' Monotremes_to_add.fa | awk 'tolower($2)~"n" {print$1}'

#Squamates: Remove
	seqtk subseq AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa <(grep -i "_A1_paralog_Squamata_" AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa | tr -d ">") | myfasta -l

#Crocodylia: Add
	seqtk subseq AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa <(grep -i "_A1_paralog_Crocodylia_" AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians.fa | tr -d ">") | myfasta -l
	mkdir /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/Crocodylia_QC
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/Crocodylia_QC
	time cat ids | xargs -n1 sh -c 'echo ">"$0; /media/aswin/programs/exon_wise -datasets -tid "$0" -er -s'
	#Check ORF
	for i in $(ls *exon_wise.fa); do j=$(transeq <(grep -v ">" $i | paste -s -d "") --auto --stdout | grep -v ">" | tr "[a-z]" "[A-Z]" | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]" | grep -v '\(*\)\1\{1,\}'); if [[ -z $j ]]; then k="dis"; else k="int"; fi; echo $i $k; unset j k; done | column -t
	for i in $(ls *exon_wise.fa); do j=$(echo $i | awk -F "_" '{print$(NF-3)"_"$(NF-2)}'); grep -v ">" $i | paste -s -d "" | sed "1i >manualncbi_A1_paralog_Crocodylia_$j"; unset j; done > Crocodylia_to_add.fa 
	#Check for "N"s gaps
	awk -v RS=">" '{$1=$1}1' Crocodylia_to_add.fa | awk 'tolower($2)~"n" {print$1}'

#Remove sequences after QC
	cat Testudines_to_remove Mammals_to_remove Monotremes_to_remove Squamata_to_remove > all_manual_selected_to_remove
	awk 'BEGIN{while((getline<"all_manual_selected_to_remove")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined.fa > AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa
#Add sequences after QC
	cat cds_sequences/manual_download/Testudines_QC/Testudines_to_add_refined.fa cds_sequences/manual_download/Monotremes_QC/Monotremes_to_add.fa cds_sequences/manual_download/Crocodylia_QC/Crocodylia_to_add.fa >> AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa

#Check again for gaps & duplicates
	cat alignment_manual_inspection/major_clades | xargs -n1 bash -c 'paste <(echo $0) <(grep $0 AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa -c)' | column -t	#Check major clades
	awk -v RS=">" '{$1=$1}1' AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa | awk 'tolower($2)~"n" {print$1}' #Check gaps
	seqkit rmdup -s -i < AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa | grep ">" -c	#Check duplicated sequences

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#12, Run Phylogeny

#Using codon alignment 81m25.080s (To run without interrruption: nohup bash -c 'time ./run4.sh' &> run4.out &)
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny
	time perl /home/ceglab25/guidance.v2.02/www/Guidance/guidance.pl --program GUIDANCE --seqFile AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa --msaProgram MAFFT --seqType codon --outDir codon_alignment_refined_monotremes_marsupials_amphibians_refined_manual_qc --bootstraps 1000 --proc_num 32
	cp codon_alignment_refined_monotremes_marsupials_amphibians_refined_manual_qc/MSA.MAFFT.aln.With_Names AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_codon.aln
	#29m30.821s
	time /media/aswin/programs/iqtree-2.3.2-Linux-intel/bin/iqtree2 --ufboot 5000 --mem 80G -T AUTO -ntmax 32 -m MFP -s AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_codon.aln

#Using nucleotide alignment
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny
	#1m39.297s
	time mafft --maxiterate 1000 --localpair --reorder --thread 32 AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa > AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_nt.aln
	#46m46.746s
	time /media/aswin/programs/iqtree-2.3.2-Linux-intel/bin/iqtree2 --ufboot 5000 --mem 80G -T AUTO -ntmax 32 -m MFP -s AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_nt.aln

#Trimmed & Run phylogeny
	#For nucleotide alignment
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny
	mkdir -p trimal/nt_alignment
	time /media/aswin/programs/trimal-1.5.0/source/trimal -in AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_nt.aln -out trimal/nt_alignment/AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_nt_trimal_gappyout.aln -gappyout
	time /media/aswin/programs/trimal-1.5.0/source/trimal -in AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_nt.aln -out trimal/nt_alignment/AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_nt_trimal_strict.aln -strict
	time /media/aswin/programs/trimal-1.5.0/source/trimal -in AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_nt.aln -out trimal/nt_alignment/AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_nt_trimal_strictplus.aln -strictplus
	time /media/aswin/programs/trimal-1.5.0/source/trimal -in AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_nt.aln -out trimal/nt_alignment/AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_nt_trimal_automated1.aln -automated1
	#Run phylogeny (To run without interrruption: nohup bash -c 'time ./run1.sh' &> run1.out &)
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/trimal/nt_alignment
	for i in $(ls | grep "AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc_nt_trimal_.*\.aln")
	do
	echo " - Run Phylogeny for Alignment: "$i 
	time /media/aswin/programs/iqtree-2.3.2-Linux-intel/bin/iqtree2 --ufboot 5000 --mem 80G -T AUTO -ntmax 32 -m MFP -s $i
	done


#NOTE:
	#Based on phylogenies published in literature the following pattern is constistenly observed:
		#AID & A3 are closer
		#A2 is closer to (AID+A3) clade than A1
		#A1 is closer to (AID+A3+A2) than A4
		#Timeline of evolution: AID & A2 -> A4 -> A1 -> A3
	#Iqtree2 makes 2 tree output files: ".contree" & ".treefile"
	#When phylogeny is visualized variations in placements observed:
		#nt contree                       : A1 lungfish placed outside of clade combining AID, A1, A2, A3
		#nt treefile                      : 
		#codon contree                    : 
		#codon treefile                   : A1 lungfish placed outside of clade combining AID, A1, A2, A3
		#nt + trimal gappyout contree     : A1 lungfish placed outside of clade combining AID, A1, A2, A3
		#nt + trimal gappyout treefile    : A1 lungfish placed outside of clade combining AID, A1, A2, A3
		#nt + trimal strict contree       : A1 lungfish placed outside of clade combining AID, A1, A2, A3
		#nt + trimal strict treefile      : 
		#nt + trimal strictplus contree   : A1 lungfish placed outside of clade combining AID, A1, A2, A3
		#nt + trimal strictplus treefile  : A1 lungfish placed outside of clade combining AID, A1, A2, A3
		#nt + trimal automated1 contree   : A1 lungfish placed outside of clade combining AID, A1, A2, A3
		#nt + trimal automated1 treefile  : A1 lungfish placed outside of clade combining AID, A1, A2, A3


		
#####################################################################################################################################################################################################################################################################################################################
#13. Clustering

	mkdir /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_clustering/AID_APOBEC_cds_unique_complete_orfs
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_clustering/AID_APOBEC_cds_unique_complete_orfs
	cp /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/AID_APOBEC_cds_unique_complete_orfs_refined.fa .
	transeq AID_APOBEC_cds_unique_complete_orfs_refined.fa --trim --auto --stdout | sed '/^>/ s/_1$//g' | /media/aswin/programs/myfasta -comb > AID_APOBEC_cds_unique_complete_orfs_refined.aa

#Updated clustering after refining final sequences

	mkdir /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_clustering/AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc
	cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_clustering/AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc
	cp /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa .
	transeq AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.fa --trim --auto --stdout | sed '/^>/ s/_1$//g' | /media/aswin/programs/myfasta -comb > AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.aa
	unzip 1789980.clans.zip
	mv 1789980.clans AID_APOBEC_cds_unique_complete_orfs_refined_monotremes_marsupials_amphibians_refined_manual_qc.clans





#####################################################################################################################################################################################################################################################################################################################3333
#####################################################################################################################################################################################################################################################################################################################
#DRAFT
#####################################################################################################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################################################################################

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Filter Gene-product metadata based on curated "NM" transcript accessions
dataformat tsv gene-product --package apobec2.zip | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | awk '$11~"NM" || $11~"Transcript_CDS_Sequence_Accession"' | column -t | less -SN

#Filter gene metadata based on annotation availability (excluding unnecessary columns such as "Summary_SourceDescription" and all columns with "-" in all rows)
dataformat tsv gene --package apobec2.zip | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | awk '$1!="No_Assembly"' \
| awk 'NR==1{for(x=1;x<=NF;x++)if($x!="Summary_SourceDescription")l[x]++;}{for(i=1;i<=NF;i++)if(i in l)printf (i==NF)?$i"":$i" ";printf "\n"}' \
| awk -v na="-" 'NR==1{for(i=1;i<=NF;i++)keep[i]=1;header=$0;next} {for(i=1;i<=NF;i++)if($i!=na)keep[i]=0;data[NR]=$0} END{n=split(header,hdr,FS);for(i=1;i<=n;i++)if(!keep[i])printf (i<n?"%s\t":"%s\n",hdr[i]);for(row=2;row<=NR;row++){n=split(data[row],fields,FS);for(i=1;i<=n;i++)if(!keep[i])printf (i<n?"%s\t":"%s\n",fields[i])}}' \
| column -t | less -SN

#Filter gene-product metadata based on gene metadata which was filtered (based on sequence with annotation & an associated ensemble ID) in the first place
grep -iwf \
<(dataformat tsv gene --package apobec2.zip | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | awk '$1!="No_Assembly"' | awk 'NR==1{for(x=1;x<=NF;x++)if($x!="Summary_SourceDescription")l[x]++;}{for(i=1;i<=NF;i++)if(i in l)printf (i==NF)?$i"":$i" ";printf "\n"}' | awk -v na="-" 'NR==1{for(i=1;i<=NF;i++)keep[i]=1;header=$0;next} {for(i=1;i<=NF;i++)if($i!=na)keep[i]=0;data[NR]=$0} END{n=split(header,hdr,FS);for(i=1;i<=n;i++)if(!keep[i])printf (i<n?"%s\t":"%s\n",hdr[i]);for(row=2;row<=NR;row++){n=split(data[row],fields,FS);for(i=1;i<=n;i++)if(!keep[i])printf (i<n?"%s\t":"%s\n",fields[i])}}' | awk '$13!="-"' | awk 'NR>1{print$14}' | sort -u) \
<(dataformat tsv gene-product --package apobec2.zip | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1') | column -t | less -SN

#Translate cds into protein
transeq ncbi_dataset/data/cds.fna -trim --auto --stdout | /media/aswin/programs/myfasta -comb | /media/aswin/programs/myfasta -rd > translated_cds.aa

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fetching data

#Download gff files given a file containing assembly accession
time datasets download genome accession --inputfile test --include gff3 --filename gff

#Download all taxonomic levels under mammals; including genus & species
datasets download taxonomy taxon 40674 --children --filename mammals_children.zip

datasets summary taxonomy taxon "Snakes" --as-json-lines --children --rank species > t1
dataformat tsv taxonomy --inputfile t1 --template "tax-summary" | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | column -t | less -SN

#Download whole taxonomic info from ncbi (2 versions available: previous & new)
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip taxdmp.zip -d previous_taxonomy
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
unzip new_taxdump.zip -d new_taxonomy

time for i in Chondrichthyes Actinopterygii Coelacanthiformes lungfishes Mammals Amphibia Testudines Crocodylia Aves Sphenodontia Gekkota Iguania Laterata Scinciformata snakes Squamata vertebrates 
do
echo ">Fetching "$i
datasets summary taxonomy taxon "$i" --as-json-lines --children --rank species > "datasets_"$i".json"
dataformat tsv taxonomy --inputfile "datasets_"$i".json" --template "tax-summary" | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | column -t > $i
done
#Rerun datasets manually if some fetching failed
ls -lh *json | awk '$5=="0"'
#Run jawless vertebrates seperately since it has 2 words
time datasets summary taxonomy taxon "jawless vertebrates" --as-json-lines --children --rank species > datasets_jawless_vertebrates.json
dataformat tsv taxonomy --inputfile datasets_jawless_vertebrates.json --template "tax-summary" | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | column -t > jawless_vertebrates
#Combine seperately downloaded groups
cat jawless_vertebrates Chondrichthyes Actinopterygii Coelacanthiformes lungfishes Mammals Amphibia Testudines Crocodylia Aves Sphenodontia Gekkota Iguania Laterata Scinciformata snakes Squamata | column -t | awk '!a[$0]++' | column -t > combined_vertebrate_groups
#Download vertebrates without restricting rank "i,e, it can include all levels" because some orthologs are from subspecies or some other sub ranks
time datasets summary taxonomy taxon "vertebrates" --as-json-lines --children > datasets_vertebrates_all.json
dataformat tsv taxonomy --inputfile datasets_vertebrates_all.json --template "tax-summary" | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | column -t > all_vertebrates
#Download whole taxonomic data from NCBI ftp
mkdir -p /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/taxonomy/ncbi_taxonomy_ftp/previous_taxonomy
cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/taxonomy/ncbi_taxonomy_ftp/previous_taxonomy
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
unzip taxdmp.zip
unzip taxcat.zip
mkdir -p /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/taxonomy/ncbi_taxonomy_ftp/new_taxonomy
cd /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/taxonomy/ncbi_taxonomy_ftp/new_taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip new_taxdump.zip

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

less mammals_children/ncbi_dataset/data/taxonomy_summary.tsv | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | awk '$5=="SPECIES" && $NF=="TRUE"' | wc -l
grep -iwf <(cat mammals_children/ncbi_dataset/data/taxonomy_summary.tsv | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | awk '$5=="SPECIES"' | awk '{print$3}') <(less ../a1_cluster_1 | grep "organism=[a-zA-Z ]\+" -o | cut -f2 -d "=" | sort -u | tr " " "_") | wc -l
grep -ivwf <(cat mammals_children/ncbi_dataset/data/taxonomy_summary.tsv | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | awk '$5=="SPECIES" && $NF=="TRUE"' | awk '{print$3}') <(less ../a1_cluster_1 | grep "organism=[a-zA-Z ]\+" -o | cut -f2 -d "=" | sort -u | tr " " "_")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Print top 50 transcripts with highest similarity score w.r.t refseq select transcripts
paste <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k9,9nr | awk '{print$1,$3}') <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k15,15nr | awk '{print$1,$3}') <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k21,21nr | awk '{print$1,$3}') | less -SN

#Print common transcripts b/w top 50 transcripts with highest similarity score w.r.t each refseq select transcripts
grep -iwf <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k9,9nr | awk '{print$0}' | awk '!a[$1]++' | head -50) <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k15,15nr | awk '{print$0}' | awk '!a[$1]++' | head -50) --color=always -z | column -t | nl
grep -iwf <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k15,15nr | awk '{print$0}' | awk '!a[$1]++' | head -50) <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k21,21nr | awk '{print$0}' | awk '!a[$1]++' | head -50) --color=always -z | column -t | nl
grep -iwf <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k9,9nr | awk '{print$0}' | awk '!a[$1]++' | head -50) <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k21,21nr | awk '{print$0}' | awk '!a[$1]++' | head -50) --color=always -z | column -t | nl

#Print number of common transcripts b/w top 50 transcripts with highest similarity score w.r.t each refseq select transcripts
grep -iwf <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k9,9nr | awk '{print$1,$3}' | awk '!a[$1]++' | head -50) <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k15,15nr | awk '{print$1,$3}' | awk '!a[$1]++' | head -50) | wc -l
grep -iwf <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k15,15nr | awk '{print$1,$3}' | awk '!a[$1]++' | head -50) <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k21,21nr | awk '{print$1,$3}' | awk '!a[$1]++' | head -50) | wc -l
grep -iwf <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k9,9nr | awk '{print$1,$3}' | awk '!a[$1]++' | head -50) <(less -SN cds_pairwise_seq_similiarity_with_refseq_select | sort -k21,21nr | awk '{print$1,$3}' | awk '!a[$1]++' | head -50) | wc -l

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Sort based of seq similarity
cat cds_pairwise_seq_similiarity_with_refseq_select | sort -k12,12nr | awk '!a[$1]++' | sed -e 's/\bbats\b/mammals/g' -e 's/\bcarnivores\b/mammals/g' -e 's/\beven-toed_ungulates__whales\b/mammals/g' -e 's/\bodd-toed_ungulates\b/mammals/g' -e 's/\binsectivores\b/mammals/g' \
 -e 's/\bmonotremes\b/mammals/g' -e 's/\bplacentals\b/mammals/g' -e 's/\bprimates\b/mammals/g' -e 's/\brabbits__hares\b/mammals/g' -e 's/\brodents\b/mammals/g' -e 's/\bwhales__dolphins\b/mammals/g' -e 's/\bmarsupials\b/mammals/g' -e 's/caecilians/amphibians/g' \
 -e 's/\bfrogs__toads\b/amphibians/g' -e 's/\bsalamanders\b/amphibians/g' -e 's/\bchimaeras\b/cartilaginous_fishes/g' -e 's/\bsharks__rays\b/cartilaginous_fishes/g' -e 's/\bhagfishes\b/jawless_vertebrates/g' -e 's/\blampreys\b/jawless_vertebrates/g' -e 's/\bhawks__eagles\b/birds/g' \
 -e 's/\blizards__snakes\b/squamates/g' | colnum.sh

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Print summary of filtered sequences 

for family in aid apobec1 apobec2 apobec3 apobec4
do
cd $family/refseq_select
i1=$(grep ">" rna_protein_taxon_filtered_transcripts.fa -c)
i2=$(/media/aswin/programs/myfasta -l rna_protein_taxon_filtered_transcripts.fa | awk '{print$NF}' | ministat -n | awk 'END{$1=""; print$0}' | sed 's/^ //g')
i3=$(transeq rna_protein_taxon_filtered_transcripts.fa -trim --auto --stdout | /media/aswin/programs/myfasta -l | awk '{print$NF}' | ministat -n | awk 'END{$1=""; print$0}' | sed 's/^ //g')
i4=$(awk 'NR>1{print$5}' sequence_similarity_and_supporting_evidence | sort | uniq -c | sort -k1,1nr | sed 's/^[ ]\+//g' | tr " " "_")
i5=$(grep "^>" rna_protein_taxon_filtered_transcripts.fa | cut -f2- -d "_" | sed 's/_[XN]M_.*//g' | sort | uniq -c | sort -k1,1nr | sed 's/^[ ]\+//g' | tr " " "_")
i6=$(for j in $(echo "$i4"); do j1=$(echo $j | cut -f2- -d "_"); j2=$(grep "$j1" <(echo "$i5")); if [[ "$j2" == "" ]]; then j2="-"; else :;fi; echo $j $j2; done | awk 'NR>1 {print"-",$0;next} {print$0}')
paste <(echo $family) <(echo "$i6") <(echo $i1 $i2 $i3)
unset i1 i2 i3 i4 i5 i6
cd ../../
done | sed '1i Member Groups_unfiltered Groups_filtered Total_No_Seq cds_N cds_Min cds_Max cds_Median cds_Avg cds_Stddev prot_N prot_Min prot_Max prot_Median prot_Avg prot_Stddev' | column -t


for i in $(cat AID_APOBEC_cds_unique_complete_orfs_refined.aa| grep mammals | awk -F "_" '{print$(NF-1),$NF}' OFS="_" | tr -d " \t")
do
echo ">"$i
awk -v a="$i" '$3==a' cds_sequences/datasets_download/apobec1/refseq_select/sequence_similarity_and_supporting_evidence
done | colnum.sh

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Combined fasta

for i in $(find . -name "a*.fa")
do
n=$(echo $i | cut -f2 -d "/")
cat $i | awk -v a="$n" '{if($1~"^>") print">"a"_"++i; else print$0}'
unset n
done > combined_cds.fa

transeq combined_cds.fa -trim --auto --stdout | /media/aswin/programs/myfasta -comb | sed '/^>/ s/_1$//g' > combined_protein.aa

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#OrhtoDB

#Download protein sequence
curl 'https://data.orthodb.org/current/fasta?id=47241at8457&seqtype=cds' -L -o A1_A1_like_Sauropsids.fa
curl 'https://data.orthodb.org/current/fasta?id=47241at8457&seqtype=protein' -L -o A1_A1_like_Sauropsids.aa
#Download other metadata
curl 'https://data.orthodb.org/current/orthologs?id=47241at8457' -L -o A1_A1_like_Sauropsids_orthologs.dat
curl 'https://data.orthodb.org/current/group?id=47241at8457' -L -o A1_A1_like_Sauropsids_group.dat
curl 'https://data.orthodb.org/current/tab?id=47241at8457' -L -o A1_A1_like_Sauropsids_annotations.tsv
cat A1_A1_like_Sauropsids_annotations.tsv | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | column -t > A1_A1_like_Sauropsids_annotations

cat A1_A1_like_Testudines_orthologs.dat | json2xml | xtract -pattern data -def "-" -sep "\t" -tab "\n" -element name -element gene_id/id -element genes/description -element aas  | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | column -t  | less -S

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

grep ">" /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/aid_a2_a3_a4_filtered.fa | awk -F "_" '{$NF=""; $(NF-1)=""; print}' OFS="_" | sed 's/[_]\+$//g' | sort | uniq -c
grep ">" /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manually_validated/validated_A1_Birds.fa | awk -F "_" '{print$1,$2,$3}' OFS="_" | sed 's/[_]\+$//g' | sort | uniq -c 
grep ">" /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/datasets_download/apobec1/refseq_select/A1_mammals_transcripts_filtered.fa | awk -F "_" '{$NF=""; $(NF-1)=""; print}' OFS="_" | sed 's/[_]\+$//g' | sort | uniq -c 
grep ">" /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/ensembl/apobec1_sauropsids/A1_A1_like_Sauropsids_ensembl.fa | awk -F "_" '{$NF=""; print}' OFS="_" | sed 's/[_]\+$//g' | sort | uniq -c 
grep ">" /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/orthoDB/A1_A1_like_Sauropsids_except_Birds_orthodb.fa | awk -F "_" '{print$1,$2,$3,$4}' OFS="_" | sed 's/[_]\+$//g' | sort | uniq -c 
grep ">" /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/annotated/birds/A1_like_Birds_annotated.fa | awk -F "_" '{$NF=""; $(NF-1)=""; print}' OFS="_" | sed 's/[_]\+$//g' | sort | uniq -c 
grep ">" /media/aswin/gene_loss/APOBEC1/APOBEC1_duplicate_phylogeny/DNA_phylogeny/cds_sequences/manual_download/annotated/birds/A1_like_Birds_curated.fa | awk -F "_" '{$NF=""; $(NF-1)=""; print}' OFS="_" | sed 's/[_]\+$//g' | sort | uniq -c 

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Check species names of all included sequences
awk -F "_" '/^>/ {print$(NF-1),$NF}' OFS="_" AID_APOBEC_cds_unique_complete_orfs_refined.fa > AID_APOBEC_cds_unique_complete_orfs_refined_accessions
time datasets summary gene accession --inputfile AID_APOBEC_cds_unique_complete_orfs_refined_accessions --as-json-lines > AID_APOBEC_cds_unique_complete_orfs_refined_accessions_datasets_summary.jsonl
time efetch -db nuccore -input AID_APOBEC_cds_unique_complete_orfs_refined_accessions -mode xml > AID_APOBEC_cds_unique_complete_orfs_refined_accessions_efetch.xml
time cat AID_APOBEC_cds_unique_complete_orfs_refined_accessions | xargs -n1 sh -c 'echo "##"$9; efetch -db nuccore -id "$0" -format docsum -mode xml' > AID_APOBEC_cds_unique_complete_orfs_refined_accessions_efetch_docsum.xml
cat AID_APOBEC_cds_unique_complete_orfs_refined_accessions_efetch.xml | xtract -patterm Seq-entry -element Textseq-id_accession,Org-ref_taxname | less

strand_for_ss=$(awk 'NR==1{for(i=1;i<=NF;++i) if($i=="Transcript_Genomic_Orientation") {n=i;break}} {print$n}' gene_product_metadata | grep -v "Transcript_Genomic_Orientation" | sort -u)
tgsss=$(awk 'NR==1{for(i=1;i<=NF;++i) if($i=="Transcript_Genomic_Start") {n=i;break}} {print$n}' gene_product_metadata | grep -v "Transcript_Genomic_Start" | head -1)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Delete columns if all row values are "-"; NOTE: can't pipe the data, need to be in a file
awk -v na="-" '
BEGIN{OFS=FS}
NR==FNR && NR>1 {
  for(i=1;i<=NF;i++){if($i!=na){s[i]=1}}
}
NR!=FNR {
  for(l in s){true} 
  for(i in s){if (i!=l){printf "%s"OFS,$i} else {printf "%s\n",$i}}
}
' t t | column -t | less -S


#Classify taxonomic groups
mammals: 
bats
carnivores
even-toed_ungulates_&_whales
odd-toed_ungulates
insectivores
marsupials
monotremes
placentals
primates
rabbits_&_hares
rodents
whales_&_dolphins

amphibians:
caecilians
frogs_&_toads
salamanders

cartilaginous_fishes:
chimaeras
sharks_&_rays

jawless_vertebrates:
hagfishes
lampreys

Actinopterygii (bony_fishes) (lobe-finned_fishes)
Coelacanthes
lunfishes

Birds:
hawks_&_eagles

Squamata:
lizards_&_snakes (is actually all squamates except snakes + Sphenodontia (beaked reptiles))
snakes

turtles:

crocodiles:

