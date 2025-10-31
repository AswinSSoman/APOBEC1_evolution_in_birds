####################################################################################################################################################################################################################################################################################################
#APOBEC1 Manual narrow search for each species search 
####################################################################################################################################################################################################################################################################################################

#Species to QC

#TOGA uncertain losses
  #- Buceros_rhinoceros
  #- Cathartes_aura
  #- Chlamydotis_macqueenii
  #- Gymnogyps_californianus
  #- Numida_meleagris        - Don't need to manually inspect here; because this species is in galliform group which we know the gene is  lost by exon loss & some events 
  #- Phoenicopterus_ruber
  #- Tyto_alba
	
#TOGA intact but lacking complete ORF
  #- Camarhynchus_parvulus
  #- Geospiza_fortis
  #- Limosa_lapponica
  #- Picoides_pubescens

#Species with gaps accourding to TOGA  
  #- Anomalopteryx_didiformis
  #- Pavo_muticus
  #- Struthio_camelus
  #- Tinamus_guttatus


####################################################################################################################################################################################################################################################################################################
#Strategy:

	       #Consensus sequences
	               |
		      |manual inspection
		      v
	#identify disruptions sequentlially in ORF
		      |		
		      |alignments : MAFFT,AliView, RevTrans
		      v
#Compare the disruption containing sequence with :- other data(local SRA,ncbi SRA & Longread) of same species
#                       |                         :- cds data of other species
		      |
		      v
#	Create all possbile variant of ORF

####################################################################################################################################################################################################################################################################################################
#Buceros_rhinoceros
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Narrow search

mkdir ~/bird_db1/Buceros_rhinoceros_silvestris/aswin/APOBEC1/manual_search
cd ~/bird_db1/Buceros_rhinoceros_silvestris/aswin/APOBEC1
exalign2 -msa5 ../APOBEC1_Serinus_canaria.fa ../gblast_edited_consensus.fa ../sblast_edited_consensus.fa ~/bird_db1/aswin/APOBEC1/TOGA/Buceros_rhinoceros/toga_rna-XM_027449866.2_exon_wise_consensus.fa ../extracted_flanking_region.fa > exon_wise_consensus.aln

#Primary & only ORF disruption is exon1 miss

#Visualize toga expected location of missing exon_1 in igv 
cd ~/bird_db1/Buceros_rhinoceros_silvestris/aswin/APOBEC1/manual_search
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Buceros_rhinoceros/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" > toga_exon_metadata.bed
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Buceros_rhinoceros/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t"  > toga_exon_metadata_exp_loc.bed 

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Genome narrow search

#Narrow down the subject genomic region for avoiding spurious hits 
bedtools getfasta -fi ../../../genome/GCF_000710305.1_ASM71030v1_genomic.fna -bed <(echo -e "NW_010389062.1\t1\t29531") > scaffold_start_to_apobec1_exon2.fa
#create database
makeblastdb -in scaffold_start_to_apobec1_exon2.fa -out scaffold_start_to_apobec1_exon2.fa -dbtype nucl

#extract all validated sequences for the exon
ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/* | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | sed "s/.fa/_exon_1/g" | sed "s/^/>/g"; grep -A1 "exon_1" $0 | grep -v ">"' > all_validated_queries_exon_1.fa

#genome blast
blastn -task blastn -db scaffold_start_to_apobec1_exon2.fa -query all_validated_queries_exon_1.fa -outfmt '6 sseqid sstart send qseqid qcovhsp sstrand qlen length qstart qend evalue bitscore sseq' | sed '1i #sseqid sstart send qseqid qcovhsp sstrand qlen aln_length qstart qend evalue bitscore sseq ' | awk '{$1=$1}1' OFS="\t" > exon1_hits_in_narrow_region.bed
#To view in IGV
cat exon1_hits_in_narrow_region.bed | sed 's/:1-29531//g' | awk '{print$1,$2,$3,$4,$5,$6}' OFS="\t" > all_gblast_hits.bed 

#sort hits based on first query_cover then evalue
cat exon1_hits_in_narrow_region.bed | sort -k5nr,5 -k11nr,11  | column -t

#Align subject with all queries
cat all_validated_queries_exon_1.fa <(cat exon1_hits_in_narrow_region.bed | sort -k5nr,5 -k11nr,11 | grep -v "#" | head -1 | awk '{print">subject_exon_1\n"$NF}') > subject_and_all_validated_queries_exon_1.fa
clustalo -i subject_and_all_validated_queries_exon_1.fa -t DNA --outfmt=clu --wrap=280 --full-iter --output-order=tree-order --resno --guidetree-out=tree.dnd --outfile=exon_1_msa.aln
#view tree based on alignment
ptree.sh tree.dnd | grep subject -z
awk '/subject/ {print">"$0}' RS=">" subject_and_all_validated_queries_exon_1.fa | grep -v ">" | sed '1i >exon_1' | awk NF > gblasted_exon_1

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#SRA search for exon_1 with start codon

bedtools getfasta -fi scaffold_start_to_apobec1_exon2.fa -bed <(cat exon1_hits_in_narrow_region.bed | sort -k5nr,5 -k11nr,11 | grep -v "#" | head -1 | awk '{print$1,$2-100,$3+100}' OFS="\t") | grep -v ">" | sed '1i >subject_exon_1_flanking' > subject_exon_1_flanking.fa
#validate in SRA
blastn -task blastn -db ../../../SRR952906_7_8.fa -query subject_exon_1_flanking.fa -evalue 0.05 -outfmt 3 -line_length 280 > exon_1_sblast.outfmt3

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Screen for start codon

#exon1 is recovered but the first 2 nucleotides in the CDS is different i,e, "GCG" instead of "ATG"
cat gblasted_exon_1 ~/bird_db1/aswin/APOBEC1/TOGA/Buceros_rhinoceros/toga_rna-XM_027449866.2_exon_wise_consensus.fa > manual_consensus_without_start.fa

#Manually choose the region containing inframe codons from exon_1
bedtools getfasta -fi ../../../genome/GCF_000710305.1_ASM71030v1_genomic.fna -bed <(echo -e "NW_010389062.1\t0\t28190") > upstream_region_of_exon_1.fa
sed '/>/ s/:.*//g' upstream_region_of_exon_1.fa -i
#FInd all start and stop codons within this region using R package "startstop"
library(startstop)
find_start_and_stop_codons("upstream_region_of_exon_1.fa", start="ATG")


#Manually take a long enough "in-frame" sequence from exon_1 to screen for the nearest in-frame start codon 
#highlight all start codons upstream
sed 's/GCTGCTGTCATTGATCGA.*//g' scaffold_start_to_apobec1_exon2.fa | grep -v ">" | rev | sed 's/[A-Z]\{3\}/& /g' | rev | grep -i "atg"
#Take the nearest start codon & look for any stop codons in between
sed 's/GCTGCTGTCATTGATCGA.*//g' scaffold_start_to_apobec1_exon2.fa | grep -v ">" | rev | sed 's/[A-Z]\{3\}/& /g' | sed 's/ GTA .*/ GTA/g' | rev | egrep "tag|tga|taa" -i

sed 's/GCTGCTGTCATTGATCGA.*//g' scaffold_start_to_apobec1_exon2.fa | grep -v ">" | rev | sed 's/[A-Z]\{3\}/& /g' | sed 's/ GTA .*/ GTA/g' | rev | tr -d " " | sed '1i >exon_1' > manual_consensus_with_start.fa
sed 's/^.*\(GCTGCTGTCATTGATCGA\)/\1/g' manual_consensus_without_start.fa | grep -v ">exon_1" >> manual_consensus_with_start.fa

#Compare 2 consenus
exalign2 manual_consensus_without_start.fa manual_consensus_with_start.fa 280

#Search whether this extended exon_1 shares any similarity with anything in NCBI nr-database
cat <(seqtk subseq manual_consensus_with_start.fa <(echo "exon_1") | sed '/>/ s/.*/>Buceros_extended_exon_1/g') <(cat Buceros_extended_exon_1_nrblast.out | sed -e '/>/ s/^[^ ]\+ />/g' -e '/>/ s/(.*//g' -e '/>/ s/,.*//g' -e 's/[ -]/_/g' -e 's/://g' -e 's/_>/_/g' -e 's/_$//g' -e 's/PREDICTED_/(P)/g') > manual_consensus_with_start_and_nrblast_hits.fa

clustalo -i manual_consensus_with_start_and_nrblast_hits.fa -t DNA --outfmt=clu --wrap=280 --full-iter --output-order=tree-order --resno --outfile=extended_exon_1_nrblast_hits.aln

#The exon_1 is missing Start codon; otherwise complete CDS is recovered filling till exon boundaries with intact splice sites
#An unpstrem in-frame START codon is found 123bp upstream
#The conslusion is intact/uncertain loss

cat manual_consensus_with_start.fa > ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Buceros_rhinoceros.fa
echo "Buceros_rhinoceros" > /home/neo/bird_db1/aswin/APOBEC1/manually_validated_queries

####################################################################################################################################################################################################################################################################################################

#Cathartes_aura

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mkdir /home/neo/bird_db1/Cathartes_aura/aswin/APOBEC1/2nd_gblast/manual_search

cd ~/bird_db1/Cathartes_aura/aswin/APOBEC1/2nd_gblast
exalign2 -msa8 APOBEC1_Mesitornis_unicolor.fa gblast_edited_consensus.fa Query_exon_combined_exonerate.out_exon_wise.fa spaln_exon_wise.fa extracted_flanking_region.fa ../sblast_edited_consensus.fa ../optimal_query_sblast.outfmt3.fa ~/bird_db1/aswin/APOBEC1/TOGA/Cathartes_aura/toga_rna-XM_027449866.2_exon_wise_consensus.fa > manual_search/exon_wise_consensus.aln

needle <(grep -v ">" ~/bird_db1/aswin/APOBEC1/TOGA/Cathartes_aura/toga_rna-XM_027449866.2_exon_wise_consensus.fa) <(grep -v ">" APOBEC1_Mesitornis_unicolor.fa) -auto -stdout -awidth 280 > manual_search/pairwise_query_toga.aln

#Visualize toga expected location of missing exons in igv 
cd ~/bird_db1/Cathartes_aura/aswin/APOBEC1/2nd_gblast/manual_search
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Cathartes_aura/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" > toga_exon_metadata_act_loc.bed
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Cathartes_aura/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t"  > toga_exon_metadata_exp_loc.bed 
awk -iinplace '{gsub(/_/,".",$1)}1' OFS="\t" toga_exon_metadata_act_loc.bed 
awk -iinplace '{gsub(/_/,".",$1)}1' OFS="\t" toga_exon_metadata_exp_loc.bed 

#Problem : when including exon_2 till splice sites this introduces a frameshift in cds by deleting last "G" of "GT" which is supposed to be a part of CDS based on query
#and if we disrupt this splice site and try to add this "G" to CDS then it creates another stop codon at exon4 mid region; this stop codon is due to a substitution not frameshift  
cat /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/Mesitornis_unicolor/aswin/APOBEC1/gblast_edited_translated ~/bird_db1/aswin/APOBEC1/TOGA/Cathartes_aura/toga_rna-XM_027449866.2_consensus_translated | less -SN

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#narrow search for exon 2
ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/* | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | sed "s/.fa/_exon_2/g" | sed "s/^/>/g"; grep -A1 "exon_2" $0 | grep -v ">"' > all_validated_queries_exon_2.fa
bedtools getfasta -fi ../../../../genome/GCA_000699945.1_ASM69994v1_genomic.fna -bed <(grep "exon_2" toga_exon_metadata_exp_loc.bed) | grep -v ">" | sed '1i >subject_exon_2_exp' > exon_2_toga_exp_region.fa
cat all_validated_queries_exon_2.fa exon_2_toga_exp_region.fa > subject_and_all_validated_queries_exon_2.fa
cat all_validated_queries_exon_4.fa exon_4_toga_exp_region.fa > subject_and_all_validated_queries_exon_4.fa
clustalo -i subject_and_all_validated_queries_exon_2.fa -t DNA --outfmt=clu --wrap=280 --full-iter --output-order=tree-order --resno --guidetree-out=tree_exon_2.dnd --outfile=exon_2_msa.aln

#narrow search for exon 4
ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/* | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | sed "s/.fa/_exon_4/g" | sed "s/^/>/g"; grep -A1 "exon_4" $0 | grep -v ">"' > all_validated_queries_exon_4.fa
bedtools getfasta -fi ../../../../genome/GCA_000699945.1_ASM69994v1_genomic.fna -bed <(grep "exon_4" toga_exon_metadata_exp_loc.bed) | grep -v ">" | sed '1i >subject_exon_4_exp' > exon_4_toga_exp_region.fa
cat all_validated_queries_exon_4.fa exon_4_toga_exp_region.fa > subject_and_all_validated_queries_exon_4.fa
clustalo -i subject_and_all_validated_queries_exon_4.fa -t DNA --outfmt=clu --wrap=280 --full-iter --output-order=tree-order --resno --guidetree-out=tree_exon_4.dnd --outfile=exon_4_msa.aln

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Since uncertaininty exists in exon boundary, create all possible variants exons 
awk '!/exon_2|exon_4/ {print">"$0}' RS=">" ~/bird_db1/aswin/APOBEC1/TOGA/Cathartes_aura/toga_rna-XM_027449866.2_exon_wise_consensus.fa | grep -v "^>$" | awk NF > base_exons.fa

#Create variants
counter=0
for i in `grep ">" exon_2_variants.fa | tr -d ">"`
do
for j in `grep ">" exon_4_variants.fa | tr -d ">"`
do
counter=$((counter + 1))
echo " - Consensus variant "$counter "=" $i "+" $j
cat base_exons.fa > consensus_variant_"$counter"_"$i"_"$j"_tmp.fa
seqtk subseq exon_2_variants.fa <(echo $i) | sed 's/_v[0-9]\+//g' >> consensus_variant_"$counter"_"$i"_"$j"_tmp.fa
seqtk subseq exon_4_variants.fa <(echo $j) | sed 's/_v[0-9]\+//g' >> consensus_variant_"$counter"_"$i"_"$j"_tmp.fa
done
unset j
done
#sort exon numbers
for i in `ls consensus_variant_* | sort -V`
do
j=`echo $i | sed 's/_tmp//g'`
echo $j
sed 's/^>/\x00&/' $i | sort -z -V | tr -d '\0'  > $j
rm $i
done

#Check for complete ORF in all possible variants 
query=../APOBEC1_Mesitornis_unicolor.fa
for i in `ls consensus_variant_* | sort -V`
do
n=`transeq <(grep -v ">" $i) -auto -stdout | grep -v ">" | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $n == "" ]] ; then o="❌";else o="✅"; fi
echo ">"$i " ORF - " $o
echo -e " - STOP codons"
transeq <(grep -v ">" $i) -auto -stdout | grep -v ">" | paste -s -d "" | grep "[A-Z]\*[A-Z]" --color=always | sed 's/^/      /g'
echo " - Protein alignment against query"
#needle <(transeq <(grep -v ">" $i) -auto -stdout | grep -v ">" | sed "1i >$i") <(transeq <(grep -v ">" $query) -auto -stdout | grep -v ">" | sed '1i >Query') -auto -stdout -awidth 280 | grep -v "#" | sed '1,2d' | sed 's/^/      /g'
needle <(transeq <(grep -v ">" $i) -auto -stdout | grep -v ">" | sed "1i >$i") <(transeq <(grep -v ">" $query) -auto -stdout | grep -v ">" | sed '1i >Query') -auto -stdout -awidth 280 | awk '/Identity/,/Score/; /^consensus/,/Query/' | sed 's/^/      /g'
echo -e "\n"
unset n o
done > all_possible_orfs

#closest variant
sed '/>/ s/_[a-z]\+//g' consensus_variant_8_exon_2_manual_exon_4_exonerate.fa > final_manual_sequence.fa
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Since none of the consensus variants are making complete ORF 
#The conclusion is loss due to substitution in exon_4

mkdir ~/bird_db1/aswin/APOBEC1/validated_lost_genes
cat final_manual_sequence.fa > ~/bird_db1/aswin/APOBEC1/validated_lost_genes/APOBEC1_Cathartes_aura.fa
echo "Cathartes_aura" > /home/neo/bird_db1/aswin/APOBEC1/manually_validated_losses

####################################################################################################################################################################################################################################################################################################
#Chlamydotis_macqueenii

cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/Chlamydotis_macqueenii/aswin/APOBEC1/2nd_gblast
mkdir manual_search
exalign2 -msa6 APOBEC1_Atlantisia_rogersi.fa gblast_edited_consensus.fa extracted_flanking_region.fa ../sblast_edited_consensus.fa ../optimal_query_sblast.outfmt3.fa ~/bird_db1/aswin/APOBEC1/TOGA/Chlamydotis_macqueenii/toga_APOBEC1_mrna_exon_wise_consensus.fa > manual_search/exon_wise_consensus.aln
needle <(grep -v ">" ~/bird_db1/aswin/APOBEC1/TOGA/Chlamydotis_macqueenii/toga_APOBEC1_mrna_exon_wise_consensus.fa) <(grep -v ">" APOBEC1_Atlantisia_rogersi.fa) -auto -stdout -awidth 280 > manual_search/pairwise_query_toga.aln

#Visualize toga expected location of missing exons in igv 
cd ~/bird_db1/Cathartes_aura/aswin/APOBEC1/2nd_gblast/manual_search
grep "APOBEC1_mrna" ~/bird_db1/aswin/APOBEC1/TOGA/Chlamydotis_macqueenii/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" > toga_exon_metadata_act_loc.bed
grep "APOBEC1_mrna" ~/bird_db1/aswin/APOBEC1/TOGA/Chlamydotis_macqueenii/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t"  > toga_exon_metadata_exp_loc.bed 
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' toga_exon_metadata_act_loc.bed -i
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' toga_exon_metadata_exp_loc.bed -i

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#First disruption is exon_2 missing

#Narrow search
ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/* | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | sed "s/.fa/_exon_2/g" | sed "s/^/>/g"; grep -A1 "exon_2" $0 | grep -v ">"' > all_validated_queries_exon_2.fa

bedtools getfasta -fi ../../../../genome/GCF_000695195.1_ASM69519v1_genomic.fna -bed <(grep "exon_2" toga_exon_metadata_exp_loc.bed) | grep -v ">" | sed '1i >toga_exon_2_exp' > exon_2_toga_exp_region.fa
cat all_validated_queries_exon_2.fa exon_2_toga_exp_region.fa > subject_and_all_validated_queries_exon_2.fa
clustalo -i subject_and_all_validated_queries_exon_2.fa -t DNA --outfmt=clu --wrap=280 --full-iter --output-order=tree-order --resno --guidetree-out=tree_exon_2.dnd --outfile=exon_2_msa.aln

#SRA search : since the exon_2 CDS is recovered but splice site is missing
#Manually extract exon)2 sequence from the MSA
blastn -task blastn -db ../../../../SRR954903_4.fa -query exon_2_with_flanking_from_genome_narrow_search.fa -evalue 0.05 -outfmt 3 -line_length 280 > exon_2_sblast.outfmt3

#Since the exon_2 splice sites are mising even in SRA make exon_2 variants

#2nd most visible problem is in exon_4 hence first try to make an ORF first 1st 4 exons 
cat exon_2_consensus_from_genome_narrow_search.fa <(cat ~/bird_db1/aswin/APOBEC1/TOGA/Chlamydotis_macqueenii/toga_APOBEC1_mrna_exon_wise_consensus.fa | awk '!/exon_5/ {print">"$0}' RS=">" | grep -v "^>$" | awk NF) | sed 's/^>/\x00&/' | sort -z -V | tr -d '\0' > manual_exon_1_to_exon_4_for_msa.fa

ls ~/bird_db1/aswin/APOBEC1/validated_sequences/* | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | sed "s/^/>/g"; grep -v ">" $0 | paste -s -d ""' > all_validated_queries_cds.fa
cat <(grep -v ">" manual_exon_1_to_exon_4_for_msa.fa | sed '1i >subject_exn_1_4') all_validated_queries_cds.fa > subject_and_all_validated_queries_cds.fa
clustalo -i subject_and_all_validated_queries_cds.fa -t DNA --outfmt=clu --wrap=280 --full-iter --output-order=tree-order --resno --guidetree-out=tree_cds.dnd --outfile=cds_msa.aln

#looks like the 2st 2 exons are fine but indels in just exon_3 is creating problem
ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/* | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | sed "s/.fa/_exon_3/g" | sed "s/^/>/g"; grep -A1 "exon_3" $0 | grep -v ">"' > all_validated_queries_exon_3.fa
bedtools getfasta -fi ../../../../genome/GCF_000695195.1_ASM69519v1_genomic.fna -bed <(grep "exon_3" toga_exon_metadata_exp_loc.bed) | grep -v ">" | sed '1i >toga_exon_3_exp' > exon_3_toga_exp_region.fa
cat all_validated_queries_exon_3.fa exon_3_toga_exp_region.fa > subject_and_all_validated_queries_exon_3.fa
clustalo -i subject_and_all_validated_queries_exon_3.fa -t DNA --outfmt=clu --wrap=280 --full-iter --output-order=tree-order --resno --guidetree-out=tree_exon_3.dnd --outfile=exon_3_msa.aln

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Compare codon and proteins of clostest query in cds msa and manual consensus 
transeq <(grep -v ">" manual_exon_1_to_exon_4_for_msa.fa | paste -s -d "" | sed "1i \>$j") manual_exon_1_to_exon_4_for_msa_orf.fa
cat <(grep -v ">" manual_exon_1_to_exon_4_for_msa.fa | paste -s -d "" | sed 's/.\{3\}/& /g') <(grep -v ">" manual_exon_1_to_exon_4_for_msa_orf.fa | paste -s -d "" | sed 's/.\{1\}/&   /g' | sed 's/^/ /g') > manual_exon_1_to_exon_4_for_msa_translated

#Exon_wise alignemnt with more closer query since the indels in exon_3 is purely based on query 
cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/Chlamydotis_macqueenii/aswin/APOBEC1/2nd_gblast
exalign2 -msa9 ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Alectura_lathami.fa ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Atlantisia_rogersi.fa ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Melopsittacus_undulatus.fa  gblast_edited_consensus.fa extracted_flanking_region.fa ../sblast_edited_consensus.fa ../optimal_query_sblast.outfmt3.fa ~/bird_db1/aswin/APOBEC1/TOGA/Chlamydotis_macqueenii/toga_APOBEC1_mrna_exon_wise_consensus.fa manual_search/manual_exon_1_to_exon_4_for_msa.fa > manual_search/exon_wise_consensus_with_closer_query.aln
cd manual_search

cat <(grep -v ">" ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Alectura_lathami.fa | paste -s -d "" |sed '1i >Alectura' ) <(grep -v ">" ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Atlantisia_rogersi.fa | paste -s -d "" |sed '1i >Atlantisia') <(grep -v ">" ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Melopsittacus_undulatus.fa|paste -s -d ""|sed '1i >Melopsittacus')  <(grep -v ">" gblast_edited_consensus.fa|paste -s -d ""|sed '1i >gedit')  <(grep -v ">" ../sblast_edited_consensus.fa|paste -s -d "" | sed '1i >sedit') <(grep -v ">" ~/bird_db1/aswin/APOBEC1/TOGA/Chlamydotis_macqueenii/toga_APOBEC1_mrna_exon_wise_consensus.fa|paste -s -d "" | sed '1i >toga') <(grep -v ">" manual_search/manual_exon_1_to_exon_4_for_msa.fa| paste -s -d "" |sed '1i >manual') > manual_search/cds_msa_for_aliview.fa

cat cds_for_aliview_msa.aln | GREP_COLORS="mt=01;32" grep "ATGGCTTCTGTCACTGATAGAG" -iz --color=always | GREP_COLORS="mt=01;33" grep "GTCTTGCACCCAAAGAGGACAACCGAGG" -iz --color=always | GREP_COLORS="mt=01;35" grep "GTGGAAGACACAGCAAAatgacttcaaaataaattatttgcctgGCCAGCATCCAAGTATAT" -iz --color=always | GREP_COLORS="mt=01;35" grep "TGTATGAAGTCAGGTGGAGCAGAGGAACCTGGAGGAGAATCTGGAGGACCTGATGCTCAAACAATTCTACCCAACATGCTGAAGTCAGCTGCTTGGAAAATTGTTTCAAGACCGTGCCATCAGTTTCTTGCTCCATCACCTGGGTCCTACCTACTACCCACTGT" -iz --color=always | GREP_COLORS="mt=01;35" grep "GGGAAATGCTCCAG" -iz --color=always | GREP_COLORS="mt=01;35" grep "AATTCTAGAGTTCCTGAGGGTACATCCCAGTGTGACCTTGGAAATATGTGCAGCCAAGATGTTCAAGCACCTGGATATCCATAACCAGCAAGGTCTCAGGAACCTGGCAATGAATGGAGTCATTATACGTATCATGAATCTTGCAG" -iz --color=always | GREP_COLORS="mt=01;36" grep "ATTACAGTTACTGGTGGAAAAGATTTGTTGCATCCCAACATGGAGAAGATGATTATTGGCCTGACAGCTTCGCTTCACACATCTTTCTGAATTTGATAGAGCTTTGTCATATCCTTT" -iz --color=always | GREP_COLORS="mt=01;36" grep " CA" -iz --color=always | GREP_COLORS="mt=01;31" grep "GGGCTTCCTCCATTGCTGGCAAACTTTGA" -iz --color=always

#Since none of the consensus variants are making complete ORF 
#The conclusion is loss due to indels in exon_3
cat exon_2_consensus_from_genome_narrow_search.fa <(cat ~/bird_db1/aswin/APOBEC1/TOGA/Chlamydotis_macqueenii/toga_APOBEC1_mrna_exon_wise_consensus.fa) | sed 's/^>/\x00&/' | sort -z -V | tr -d '\0' | sed 's/ATACTAA//g'> final_manual_consensus.fa
cat final_manual_sequence.fa > ~/bird_db1/aswin/APOBEC1/validated_lost_genes/APOBEC1_Chlamydotis_macqueenii.fa
echo "Chlamydotis_macqueenii" >> /home/neo/bird_db1/aswin/APOBEC1/manually_validated_losses


################################################################################################################################################################################################################################################################################################################
#Gymnogyps_californianus

cd ~/bird_db1/Gymnogyps_californianus/aswin/APOBEC1/2nd_gblast
mkdir manual_search
exalign2 -msa8 APOBEC1_Mesitornis_unicolor.fa gblast_edited_consensus.fa Query_exon_combined_exonerate.out_exon_wise.fa spaln_exon_wise.fa extracted_flanking_region.fa ../sblast_edited_consensus.fa ../optimal_query_sblast.outfmt3.fa ~/bird_db1/aswin/APOBEC1/TOGA/Gymnogyps_californianus/toga_rna-XM_027449866.2_exon_wise_consensus.fa > manual_search/exon_wise_consensus.aln
needle <(grep -v ">" ~/bird_db1/aswin/APOBEC1/TOGA/Gymnogyps_californianus/toga_rna-XM_027449866.2_exon_wise_consensus.fa) <(grep -v ">" APOBEC1_Mesitornis_unicolor.fa) -auto -stdout -awidth 280  > manual_search/pairwise_query_toga.aln

#Very similar to Cathartes aura; a single del exon_2 makes the tools to lose one base to fit the splice site criteria
#If we relax the splice cirteria and align then a stop codon arises in the exon_4 exactly like Cathartes aura 
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Since uncertaininty exists in exon boundary, create all possible variants exons 
awk '!/exon_2|exon_4/ {print">"$0}' RS=">" /home/neo/bird_db1/aswin/APOBEC1/TOGA/Gymnogyps_californianus/toga_rna-XM_027449866.2_exon_wise_consensus.fa | grep -v "^>$" | awk NF > base_exons.fa

#Create variants
counter=0
for i in `grep ">" exon_2_variants.fa | tr -d ">"`
do
for j in `grep ">" exon_4_variants.fa | tr -d ">"`
do
counter=$((counter + 1))
echo " - Consensus variant "$counter "=" $i "+" $j
cat base_exons.fa > consensus_variant_"$counter"_"$i"_"$j"_tmp.fa
seqtk subseq exon_2_variants.fa <(echo $i) | sed 's/_v[0-9]\+//g' >> consensus_variant_"$counter"_"$i"_"$j"_tmp.fa
seqtk subseq exon_4_variants.fa <(echo $j) | sed 's/_v[0-9]\+//g' >> consensus_variant_"$counter"_"$i"_"$j"_tmp.fa
done
unset j
done
#sort exon numbers
for i in `ls consensus_variant_* | sort -V`
do
j=`echo $i | sed 's/_tmp//g'`
echo $j
sed 's/^>/\x00&/' $i | sort -z -V | tr -d '\0'  > $j
rm $i
done

#Check for complete ORF in all possible variants 
query=../APOBEC1_Mesitornis_unicolor.fa
for i in `ls consensus_variant_* | sort -V`
do
n=`transeq <(grep -v ">" $i) -auto -stdout | grep -v ">" | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $n == "" ]] ; then o="❌";else o="✅"; fi
echo ">"$i " ORF - " $o
echo -e " - STOP codons"
transeq <(grep -v ">" $i) -auto -stdout | grep -v ">" | paste -s -d "" | grep "[A-Z]\*[A-Z]" --color=always | sed 's/^/      /g'
echo " - Protein alignment against query"
#needle <(transeq <(grep -v ">" $i) -auto -stdout | grep -v ">" | sed "1i >$i") <(transeq <(grep -v ">" $query) -auto -stdout | grep -v ">" | sed '1i >Query') -auto -stdout -awidth 280 | grep -v "#" | sed '1,2d' | sed 's/^/      /g'
needle <(transeq <(grep -v ">" $i) -auto -stdout | grep -v ">" | sed "1i >$i") <(transeq <(grep -v ">" $query) -auto -stdout | grep -v ">" | sed '1i >Query') -auto -stdout -awidth 280 | awk '/Identity/,/Score/; /^consensus/,/Query/' | sed 's/^/      /g'
echo -e "\n"
unset n o
done > all_possible_orfs

#closest variant
sed '/>/ s/_[a-z]\+//g' consensus_variant_6_exon_2_manual_exon_4_exonerate.fa > final_manual_sequence.fa
consensus_variant_6_exon_2_manual_exon_4_exonerate.fa

#Hence  directly conclude that it is loss due to either deletion in exon_1 with stop codon at exon_4 or stop codon in exon_3 (if splice site is considered)

cat final_manual_sequence.fa > ~/bird_db1/aswin/APOBEC1/validated_lost_genes/APOBEC1_Gymnogyps_californianus.fa
echo "Gymnogyps_californianus" >> /home/neo/bird_db1/aswin/APOBEC1/manually_validated_losses

################################################################################################################################################################################################################################################################################################################
#Phoenicopterus_ruber

cd ~/bird_db1/Phoenicopterus_ruber/aswin/APOBEC1/2nd_gblast
mkdir manual_search
exalign2 -msa8 APOBEC1_Podilymbus_podiceps.fa gblast_edited_consensus.fa Query_exon_combined_exonerate.out_exon_wise.fa spaln_exon_wise.fa extracted_flanking_region.fa ../sblast_edited_consensus.fa ../optimal_query_sblast.outfmt3.fa ~/bird_db1/aswin/APOBEC1/TOGA/Phoenicopterus_ruber/toga_rna-XM_027449866.2_exon_wise_consensus.fa > manual_search/exon_wise_consensus.aln
needle <(grep -v ">" ~/bird_db1/aswin/APOBEC1/TOGA/Phoenicopterus_ruber/toga_rna-XM_027449866.2_exon_wise_consensus.fa) <(grep -v ">" APOBEC1_Podilymbus_podiceps.fa) -auto -stdout -awidth 280  > manual_search/pairwise_query_toga.aln

#There are a lof of differences between genome & sra : some of these sra are frameshifting variants : So create all posiible varaints with and without frameshifts

#To decipher which one is the best/closest to query:  we need to align amino acid guided dna alignment of all the consensus including query
cat <(grep -v ">" APOBEC1_Podilymbus_podiceps.fa|sed '1i >Query') <(grep -v ">" gblast_edited_consensus.fa| sed '1i >gedit') <(grep -v ">" ../sblast_edited_consensus.fa|sed '1i >sedit') <(grep -v ">" ~/bird_db1/aswin/APOBEC1/TOGA/Phoenicopterus_ruber/toga_rna-XM_027449866.2_exon_wise_consensus.fa| sed '1i >toga') > manual_search/revtrans_fasta.fa

#Create variants of sra edited consensus
#distinguish actual frameshiting variants from artifacts of blast : so comapare sblast query with artifacts & sblast hits with retained artifacts 
sed 's/-/N/g' ../../optimal_query_sblast.outfmt3_processed_data/optimal_query_sblast.outfmt3.fa_with_dashes > sblast_hit_with_N.fa
cat ../../optimal_query_sblast.outfmt3 | grep "Query_" | awk '{print">"$1"\n"$3}' | sed 's/Query/exon/g' | sed 's/-/N/g' | awk '/>/{id=$0} !/>/{seq[id]=seq[id] $0} END{for(id in seq) print id "\n" seq[id]}' > sblast_query_with_N.fa 
#Identify the position of frameshifting variant, extract 3 bases flanking to the variant to use as string to locate the region in msa  
p1=`ls ../../optimal_query_sblast.outfmt3_processed_data/input_?.txt | sort -V | xargs -I {} head -1 {} | paste -s -d "" | cut -c 244-250 | sed 's/-/N/g'`
#Locate the region in msa & seperately create consensus sequnces with all possible combinations of variants manually : replace only the exons with variants in different consensus sequences other exons will remain same
exalign2 sblast_query_with_N.fa sblast_hit_with_N.fa 280 | grep $p1 -z
#Repeat this for other variants
p2=`ls ../../optimal_query_sblast.outfmt3_processed_data/input_?.txt | sort -V | xargs -I {} head -1 {} | paste -s -d "" | cut -c 187-193 | sed 's/-/N/g'`
exalign2 sblast_query_with_N.fa sblast_hit_with_N.fa 280 | grep $p2 -z
#Append these sequences to revtrans_fasta.fa

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Visulaize DNA & amino acid alignments at the same time using Aliview
#NOTE : Aliview coloring is based on translation after alignment not translation of each sequence ; So a codon grouped and coloured as stop codon in Aliview may not be a stop codon in that sequence; but aliview alignment added dashes to    
				#Aliview
				   |
				   V
				#align 
				   |
				   V
	#identify obvious mistakes & edit them in sequences (eg: here splice site 1 was corrected w.r.t query)
	#find which variants are most similar to query & maintain ORF (but don't edit the consensus sequence)

#Perform amino acid guided dna alignment using revtrans server and save it in the output file "revtrans_fasta.aln"
#View the ouput 
for i in `grep "^>" revtrans_fasta.aln | tr -d ">"`; do; echo ">"$i
sed 's/^>/@/g' revtrans_fasta.aln | awk -v n="$i" '$0~n {print">"$0}' RS="@" | awk '/^$/{n=n RS}; /./{printf "%s",n; n=""; print}' | egrep -v "^>|^Reading" | cut -c 4- | sed 's/ [0-9]\+$//g' | tr " " "#" | awk NF | sed 's/#$//g' | awk -v m=90 '{printf("%-90s\n", $0)}' | tr " " "#" | paste -d " " - - - | awk '{for(i=1;i<=NF;i++) {a[NR,i]=$i}} NF>p{p=NF} END{for(j=1;j<=p;j++) {str=a[1,j]; for(i=2; i<=NR; i++){str=str" "a[i,j];}print str}}' | tr -d " " | tr "#" " "
echo -e "\n"; done > revtrans_fasta_linear.aln
cat <(cut -c -310 tmp) <(cut -c 311-619 tmp) | GREP_COLORS="mt=07;33" grep "[A-Z]  [A-Z]" -z
#Normal align
clustalo -i revtrans_fasta.fa -t DNA --outfmt=clu --wrap=280 --full-iter --output-order=tree-order --resno
mafft --maxiterate 1000 --localpair --clustalout --reorder --quiet --anysymbol  --op 1.53 --ep 0.123 revtrans_fasta.fa

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#pairwise align each consensus sequences with Query

for i in `grep ">" revtrans_fasta.fa | tr -d ">" | grep -v "Query"`
do 
n=`transeq <(seqtk subseq revtrans_fasta.fa <(echo $i)) -auto -stdout | grep -v ">" | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $n == "" ]] ; then o="❌";else o="✅"; fi
echo -e ">"$i "ORF - " $o "\n"
echo " - STOP codons :"
transeq <(seqtk subseq revtrans_fasta.fa <(echo $i)) -auto -stdout | grep -v ">" | paste -s -d "" | grep "[A-Z]\*[A-Z]" --color=always | sed 's/^/    /g'
echo ""
echo " - Protein alignment against query"
needle <(transeq <(seqtk subseq revtrans_fasta.fa <(echo "Query")) -auto -stdout) <(transeq <(seqtk subseq revtrans_fasta.fa <(echo $i)) -auto -stdout) -auto -stdout -awidth 280 | grep -v "#" | awk NF | sed 's/^/    /g'
echo -e "\n"
unset n o
done > all_possible_orfs

#The conclusion is gene is intact
awk '/sedit_197_247G_exon_1_corrected/ {print">"$0}' RS=">" revtrans_fasta.fa > manual_consensus.fa
#add exon names as well

cat manual_consensus.fa > ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Phoenicopterus_ruber.fa
echo "Phoenicopterus_ruber" >> /home/neo/bird_db1/aswin/APOBEC1/manually_validated_queries

################################################################################################################################################################################################################################################################################################################
# Tyto_alba

cd ~/bird_db1/Tyto_alba/aswin/APOBEC1/2nd_gblast
mkdir manual_search
exalign2 -msa8 APOBEC1_Galbula_dea.fa gblast_edited_consensus.fa Query_exon_combined_exonerate.out_exon_wise.fa spaln_exon_wise.fa extracted_flanking_region.fa ../sblast_edited_consensus.fa ../optimal_query_sblast.outfmt3.fa ~/bird_db1/aswin/APOBEC1/TOGA/Tyto_alba/toga_rna-XM_027449866.2_exon_wise_consensus.fa > manual_search/exon_wise_consensus.aln
needle <(grep -v ">" ~/bird_db1/aswin/APOBEC1/TOGA/Phoenicopterus_ruber/toga_rna-XM_027449866.2_exon_wise_consensus.fa) <(grep -v ">" APOBEC1_Podilymbus_podiceps.fa) -auto -stdout -awidth 280  > manual_search/pairwise_query_toga.aln


ls APOBEC1_Galbula_dea.fa gblast_edited_consensus.fa Query_exon_combined_exonerate.out_exon_wise.fa spaln_exon_wise.fa ../sblast_edited_consensus.fa ~/bird_db1/aswin/APOBEC1/TOGA/Tyto_alba/toga_rna-XM_027449866.2_exon_wise_consensus.fa | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | cut -f1,2 -d "_" | sed "s/Query_exon/Exonerate/g" | sed "s/_rna-XM//g" | sed "s/^/>/g"; grep -v ">" $0 | paste -s -d ""' > manual_search/all_consensus.fa

for i in `grep ">" all_consensus.fa | tr -d ">" | grep -v "APOBEC1_Galbula"`
do 
n=`transeq <(seqtk subseq all_consensus.fa <(echo $i)) -auto -stdout | grep -v ">" | paste -s -d "" | grep "^M.*\*$" | grep -v "[A-Z]\*[A-Z]"`
if [[ $n == "" ]] ; then o="❌";else o="✅"; fi
echo -e ">"$i "ORF - " $o "\n"
echo " - STOP codons :"
transeq <(seqtk subseq all_consensus.fa <(echo $i)) -auto -stdout | grep -v ">" | paste -s -d "" | grep "[A-Z]\*[A-Z]" --color=always | sed 's/^/    /g'
echo ""
echo " - Protein alignment against query"
needle <(transeq <(seqtk subseq all_consensus.fa <(echo "APOBEC1_Galbula")) -auto -stdout) <(transeq <(seqtk subseq all_consensus.fa <(echo $i)) -auto -stdout) -auto -stdout -awidth 280 | grep -v "#" | awk NF | sed 's/^/    /g'
echo -e "\n"
unset n o
done

#Exon 1-3 is fine, but exon_4 is incomplete in 3' end and SRA reads lacks coverage in this region

ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/* | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | sed "s/.fa/_exon_4/g" | sed "s/^/>/g"; grep -A1 "exon_4" $0 | grep -v ">"' > all_validated_queries_exon_4.fa

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Check RNA

scp ../Query_exon_combined_exonerate.out_exon_wise.fa sagar@172.30.1.172:/media/sagar/disk4/RNAseq_store/Birds/Tyto_alba_alba/aswin/APOBEC1/
blastn -task blastn -db ../../GCA_000687205.1_ASM68720v1_genomic.fna -query Query_exon_combined_exonerate.out_exon_wise.fa -evalue 0.01 -outfmt "6 sseqid sstart send qseqid sstrand" | awk '!a[$4]++' | awk '{if($3>$2) print$1,$2,$3,$4,"1","+"; else print$1,$3,$2,$4,"-"}' OFS="\t" > Tyto_alba_exonerate_APOBEC1_exons.bed

paste <(awk '{print$1}' Tyto_alba_exonerate_APOBEC1_exons.bed | sort -u) <(awk '{print$2,$3}' OFS="\n" Tyto_alba_exonerate_APOBEC1_exons.bed | sort -n | sed -n '1p;$p' | paste -s -d "\t") <(echo "APOBEC1") <(awk '{print$NF}' Tyto_alba_exonerate_APOBEC1_exons.bed | sort -u) -d "\t" > Tyto_alba_exonerate_APOBEC1_gene.bed

#APOBEC1
cd /media/sagar/disk4/RNAseq_store/Birds/Tyto_alba_alba/aswin/APOBEC1
c=`awk '{print$1}' Tyto_alba_exonerate_APOBEC1_gene.bed`
p=`awk '{print$2,$3}' Tyto_alba_exonerate_APOBEC1_gene.bed | tr " " "-"`
/media/sagar/disk4/RNAseq_store/Birds/Tyto_alba_alba/
for i in `find /media/sagar/disk4/RNAseq_store/Birds/Tyto_alba_alba -maxdepth 1 -name "*.bam" -not -name "*subset*"`
do
j=`echo $i | awk -F '/' '{print$NF}' | sed 's/^\([SED]RR[0-9]\+\).*/\1/g'`
k=`~/programs/samtools-1.11/samtools coverage -r $(echo -e "$c:$p") $i | grep -v "#"`
echo $j $k
unset j k
done | sed '1i SRR_id Tissue #rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq' | column -t 
unset c p i

~/aswin/programs/bamcovplot -bam ../../ERR3394421Aligned.sortedByCoord.out.bam -b Tyto_alba_exonerate_APOBEC1_exons.bed:APOBEC1 -ef -isa -id ERR3394421:Testis -o APOBEC1_Tyto_alba
#No RNA expn 

#None of the databases tried recovered the 3' end of exon_4.
#Since it is possible to recover the sequence if we add more data, however iyt is very time consuming & even if it has a stop codon it would be in the >85% of the gene length, which is not a reliable way to establish gene loss
echo "Tyto_alba" >> species_to_exclude

################################################################################################################################################################################################################################################################################################################

#Camarhynchus_parvulus

cd ~/bird_db1/Camarhynchus_parvulus/aswin/APOBEC1/2nd_gblast
mkdir manual_search

exalign2 -msa8 APOBEC1_Zonotrichia_albicollis.fa gblast_edited_consensus.fa Query_exon_combined_exonerate.out_exon_wise.fa spaln_exon_wise.fa ~/bird_db1/aswin/APOBEC1/TOGA/Camarhynchus_parvulus/toga_rna-XM_027449866.2_exon_wise_consensus.fa extracted_flanking_region.fa ../sblast_edited_consensus.fa ../optimal_query_sblast.outfmt3.fa > manual_search/exon_wise_consensus.aln

#Problem - START codon is missing
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#NOTE : The gene is annotated in genomic anntation as LOC_ID & it has an extra exon before the exon_1 we found in our study; however this annotation doesn't have supporting evidence from RNA-seq; it only has similarity to 6 proteins
#SRA local - coverage in exon_1 is above minimu thresold; but not very huigh
#RNA local - No RNA-seq expn found
#SRA ncbi  -

#Screen for start codon

#exon1 is recovered but the first 2 nucleotides in the CDS is different i,e, "AGG" instead of "ATG"
cat Query_exon_combined_exonerate.out_exon_wise.fa > manual_search/manual_consensus_without_start.fa

#Manually choose the region containing inframe codons from exon_1
bedtools getfasta -fi ../../../genome/GCF_901933205.1_STF_HiC_genomic.fna -bed <(blastn -task blastn -db ../../../genome/GCF_901933205.1_STF_HiC_genomic.fna -query <(grep "exon_1" -A1 manual_search/manual_consensus_without_start.fa) -evalue 0.01 -outfmt '6 sseqid sstart send qseqid qcovhsp sstrand' | awk '{if($NF=="minus") print$1,$3,$2+500,$4,$5,"-"; else print$1,$2-500,$3,$4,$5,"+"}' OFS="\t") -s | sed '/>/ s/.*/>exon_1_upstream_500/g' > manual_search/exon_1_upstream_500.fa
#Find the in-frame sequence of exon_1 to screen for in-frame start codon upstream
cd manual_search/
inframe_exon_1=`grep "exon_1" manual_consensus_without_start.fa -A1 | grep -v ">" | sed 's/.\{3\}/& /g' | awk '{for(i=1;i<=NF;i++) if(length($i)==3) print$i}' | paste -s -d " "`
#Visualize all the start (including different varaints) and stop codons upstream f exon_1
grep -v ">" exon_1_upstream_500.fa | rev | sed 's/.\{3\}/& /g' | rev | GREP_COLORS="mt=44" grep -i "$inframe_exon_1" --color=always | GREP_COLORS="mt=42" egrep -i "ttg|ctg|atg" --color=always | GREP_COLORS="mt=41" egrep -i "tag|tga|taa" --color=always
unset inframe_exon_1
#Manually creat consenus sequence with the nearest start codon

#Conslusion : Lost
#Start codon was missing; the nearst in-frame start codon includes a stop codon as well

cat manual_consensus_with_upstream_start.fa > ~/bird_db1/aswin/APOBEC1/validated_lost_genes/APOBEC1_Camarhynchus_parvulus.fa
echo "Camarhynchus_parvulus" >> /home/neo/bird_db1/aswin/APOBEC1/manually_validated_losses

################################################################################################################################################################################################################################################################################################################

#Geospiza_fortis

cd ~/bird_db1/Geospiza_fortis/aswin/APOBEC1/2nd_gblast/
mkdir manual_search

exalign2 -msa8 APOBEC1_Zonotrichia_albicollis.fa gblast_edited_consensus.fa Query_exon_combined_exonerate.out_exon_wise.fa spaln_exon_wise.fa ~/bird_db1/aswin/APOBEC1/TOGA/Geospiza_fortis/toga_rna-XM_027449866.2_exon_wise_consensus.fa extracted_flanking_region.fa ../sblast_edited_consensus.fa ../optimal_query_sblast.outfmt3.fa > manual_search/exon_wise_consensus.aln

#Problem - START codon is missing
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#NOTE : The gene is annotated in genomic anntation as LOC_ID & the gene has only 3 exons & the start codon is present in exon_3; however this annotation doesn't have supporting evidence from RNA-seq; it only has similarity to 6 proteins
#SRA local - Coverage is high
#RNA local - before screening for start codon - NO RNA expn found
#SRA ncbi  -

#Screen for start codon

#exon1 is recovered but the first 2 nucleotides in the CDS is different i,e, "AGG" instead of "ATG"
cat Query_exon_combined_exonerate.out_exon_wise.fa > manual_search/manual_consensus_without_start.fa
#Manually choose the region containing inframe codons from exon_1
bedtools getfasta -fi ../../../genome/GCF_000277835.1_GeoFor_1.0_genomic.fna -bed <(blastn -task blastn -db ../../../genome/GCF_000277835.1_GeoFor_1.0_genomic.fna -query <(grep "exon_1" -A1 manual_search/manual_consensus_without_start.fa) -evalue 0.01 -outfmt '6 sseqid sstart send qseqid qcovhsp sstrand' | awk '{if($NF=="minus") print$1,$3,$2+500,$4,$5,"-"; else print$1,$2-500,$3,$4,$5,"+"}' OFS="\t") -s | sed '/>/ s/.*/>exon_1_upstream_500/g' > manual_search/exon_1_upstream_500.fa
#Find the in-frame sequence of exon_1 to screen for in-frame start codon upstream
cd manual_search/
inframe_exon_1=`grep "exon_1" manual_consensus_without_start.fa -A1 | grep -v ">" | sed 's/.\{3\}/& /g' | awk '{for(i=1;i<=NF;i++) if(length($i)==3) print$i}' | paste -s -d " "`
#Visualize all the start (including different varaints) and stop codons upstream of exon_1
grep -v ">" exon_1_upstream_500.fa | rev | sed 's/.\{3\}/& /g' | rev | GREP_COLORS="mt=44" grep -i "$inframe_exon_1" --color=always | GREP_COLORS="mt=42" egrep -i "ttg|ctg|atg" --color=always | GREP_COLORS="mt=41" egrep -i "tag|tga|taa" --color=always
unset inframe_exon_1
#Manually create consensus sequence with the nearest start codon

#Conslusion : Lost
#Start codon was missing; the nearst in-frame start codon includes a stop codon as well

cat manual_consensus_with_upstream_start.fa > ~/bird_db1/aswin/APOBEC1/validated_lost_genes/APOBEC1_Geospiza_fortis.fa
echo "Geospiza_fortis" >> /home/neo/bird_db1/aswin/APOBEC1/manually_validated_losses

################################################################################################################################################################################################################################################################################################################

#Limosa_lapponica_baueri 

cd ~/bird_db1/Limosa_lapponica_baueri/aswin/APOBEC1/2nd_gblast
mkdir manual_search
exalign2 -msa8 APOBEC1_Calidris_pugnax.fa gblast_edited_consensus.fa Query_exon_combined_exonerate.out_exon_wise.fa spaln_exon_wise.fa ~/bird_db1/aswin/APOBEC1/TOGA/Limosa_lapponica/toga_rna-XM_027449866.2_exon_wise_consensus.fa extracted_flanking_region.fa ../sblast_edited_consensus.fa ../optimal_query_sblast.outfmt3.fa > manual_search/exon_wise_consensus.aln

#Expected exon boundaries from toga
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Limosa_lapponica/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > manual_search/toga_exon_metadata.bed
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Limosa_lapponica/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > manual_search/toga_exon_metadata_exp_loc.bed 

#All SRA data from NCBI
sra-meta -name Limosa_lapponica > manual_search/sra_ncbi

#Problem - STOP codon is missing
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#SRA - good support without variants 
#RNA - No RNA-seq data

#Screen for stop codon

cat Query_exon_combined_exonerate.out_exon_wise.fa > manual_search/manual_consensus_without_stop.fa
#Manually choose the region containing inframe codons from exon_5
bedtools getfasta -fi ../../../genome/GCA_002844005.1_Godwit_v1_genomic.fna -bed <(blastn -task blastn -db ../../../genome/GCA_002844005.1_Godwit_v1_genomic.fna -query <(grep "exon_5" -A1 manual_search/manual_consensus_without_stop.fa) -evalue 0.01 -outfmt '6 sseqid sstart send qseqid qcovhsp sstrand' | awk '{if($NF=="minus") print$1,$3-500,$2+1,$4,$5,"-"; else print$1,$2-1,$3+500,$4,$5,"+"}' OFS="\t") | sed '/>/ s/.*/>exon_5_downstream_500/g' > manual_search/exon_5_downstream_500.fa
inframe_exon_5=`grep "exon_5" manual_consensus_without_stop.fa -A1 | grep -v ">" | sed 's/.\{3\}/& /g' | awk '{for(i=1;i<=NF;i++) if(length($i)==3) print$i}' | paste -s -d " "`
#Visualize all the start (including different varaints) and stop codons downstream of exon_5
grep -v ">" exon_5_downstream_500.fa | sed 's/.\{3\}/& /g' | GREP_COLORS="mt=44" grep -i "$inframe_exon_5" --color=always | GREP_COLORS="mt=41" egrep -i "tag|tga|taa" --color=always
unset inframe_exon_5

#Creat consensus with nearest stop codon
awk '!/exon_5/ {print">"$0}' RS=">" manual_consensus_without_stop.fa | grep -v "^>$" | awk NF > manual_consensus_with_downstream_stop.fa 
grep -v ">" exon_5_downstream_500.fa | sed 's/.\{3\}/& /g' | sed 's/TAG.*/TAG/g' | tr -d " " | sed '1i >exon_5' >> manual_consensus_with_downstream_stop.fa 

#Conslusion : Intact
#Stop codon was missing; the nearest in-frame stop codon is 66bp downstream to the expected position

cat manual_consensus_with_downstream_stop.fa > ~/bird_db1/aswin/APOBEC1/validated_sequences/APOBEC1_Limosa_lapponica.fa
echo "Limosa_lapponica" >> ~/bird_db1/aswin/APOBEC1/manually_validated_queries

################################################################################################################################################################################################################################################################################################################

#Picoides_pubescens

cd ~/bird_db1/Picoides_pubescens/aswin/APOBEC1/2nd_gblast
mkdir manual_search
exalign2 -msa6 APOBEC1_Galbula_dea.fa gblast_edited_consensus.fa ~/bird_db1/aswin/APOBEC1/TOGA/Picoides_pubescens/toga_rna-XM_027449866.2_exon_wise_consensus.fa extracted_flanking_region.fa ../sblast_edited_consensus.fa ../optimal_query_sblast.outfmt3.fa > manual_search/exon_wise_consensus.aln

#Expected exon boundaries from toga
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Picoides_pubescens/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > manual_search/toga_exon_metadata.bed
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Picoides_pubescens/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > manual_search/toga_exon_metadata_exp_loc.bed 

#Narrower search
cd manual_search/
c=`awk '{print$1}' toga_exon_metadata_exp_loc.bed | sort -u`
p=`awk '{print$2,$3}' OFS="\n" toga_exon_metadata_exp_loc.bed | sort -n | sed -n '1p;$p' | paste -s -d "\t"`
bedtools getfasta -fi ../../../../genome/GCF_000699005.1_ASM69900v1_genomic.fna -bed <(echo -e "$c\t$p\tapobec1\t1\t-") -s | sed '/>/ s/.*/>scaffold_containing_apobec1/g' > scaffold_containing_apobec1.fa
makeblastdb -in scaffold_containing_apobec1.fa -out scaffold_containing_apobec1.fa -dbtype nucl
ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/* | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | sed "s/.fa/_exon_1/g" | sed "s/^/>/g"; grep -A1 "exon_1" $0 | grep -v ">"' | awk '!/Buceros/ {print">"$0}' RS=">" | grep -v "^>$" | awk NF > all_validated_queries_exon_1.fa

#blast all exon_1 queries against narrow region
blastn -task blastn -db scaffold_containing_apobec1.fa -query all_validated_queries_exon_1.fa -evalue 0.05 -outfmt '6 sseqid sstart send qseqid qcovhsp sstrand qlen length qstart qend evalue bitscore sseq' | sed '1i #sseqid sstart send qseqid qcovhsp sstrand qlen aln_length qstart qend evalue bitscore sseq ' | awk '{$1=$1}1' OFS="\t" | column -t | grep plus | sort -k5nr,5 -k11nr,11

#MSA all the exon_1
#Concatanate all the validated exon_1 queries and the be hits from narrow region search
bedtools getfasta -fi scaffold_containing_apobec1.fa -bed <(echo -e "scaffold_containing_apobec1\t3999\t4018") | grep -v ">" | sed '1i >PP_3999-4018' > gblasted_exon_1_and_all_validated_queries_exon_1.fa 
bedtools getfasta -fi scaffold_containing_apobec1.fa -bed <(echo -e "scaffold_containing_apobec1\t5880\t5900") | grep -v ">" | sed '1i >PP_5880_5900' >> gblasted_exon_1_and_all_validated_queries_exon_1.fa 
cat all_validated_queries_exon_1.fa >> gblasted_exon_1_and_all_validated_queries_exon_1.fa 
clustalo -i gblasted_exon_1_and_all_validated_queries_exon_1.fa -t DNA --outfmt=clu --wrap=280 --full-iter --output-order=tree-order --resno
mafft --maxiterate 1000 --localpair --clustalout --reorder --quiet --anysymbol  --op 1.53 --ep 0.123 gblasted_exon_1_and_all_validated_queries_exon_1.fa

#blast all exon_1 queries against the whole genome
blastn -task blastn-short -db ../../../../genome/GCF_000699005.1_ASM69900v1_genomic.fna -query gblasted_exon_1_and_all_validated_queries_exon_1.fa -outfmt '6 sseqid sstart send qseqid qcovhsp sstrand qlen length qstart qend evalue bitscore sseq' | sed '1i #sseqid sstart send qseqid qcovhsp sstrand qlen aln_length qstart qend evalue bitscore sseq ' | awk '{$1=$1}1' OFS="\t" | column -t | grep minus | grep "NW_009666054.1" | sort -k5nr,5 -k11nr,11 

sra-meta -name Picoides_pubescens > sra_data

#Colnclusion : missing data
#The terminal exons like exon_1 is very near assembly; although ther is 6.9kb region after exon_2, this region might suffer from poor assembly quality dince it is towards the end 
#The phylogeny also doesn't show any shared loss with nearest species; hence it is better to exclude this species
echo "Picoides_pubescens" >> species_to_exclude

################################################################################################################################################################################################################################################################################################################
#Species with gaps

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Anomalopteryx_didiformis

#The synteny of APOBEC1 is verified from toga chain alignment
#only exon 4 is found but it's incomplete
grep Anomal ~/bird_db1/aswin/APOBEC1/pipeline_summary_row_wise -A3
#Exon 4 contains N's in genome & no hits in SRA 
cat ~/bird_db1/aswin/APOBEC1/TOGA/Anomalopteryx_didiformis/All_consensus_msa.aln
#Expected exon boundaries from toga
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Anomalopteryx_didiformis/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > manual_search/toga_exon_metadata.bed
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Anomalopteryx_didiformis/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > manual_search/toga_exon_metadata_exp_loc.bed 
#Unlike other palaeognathes the exon is not jut missing, but the expected region where we could find the gene contains N's

#Since none of the exons are intact enough to show any ORF disruption, there is a chance that the gene might be intact and couldnot be sequenced
#This bird is an extintc bird as well hence it s not feasible to sequnce it again; hence exclude this species
echo "Anomalopteryx_didiformis" >> species_to_exclude 

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Struthio_camelus

mkdir ~/bird_db1/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/manual_search
cd ~/bird_db1/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/manual_search

#Compare results from pipelines
grep Struthio ~/bird_db1/aswin/APOBEC1/pipeline_summary_row_wise -A3 > pipeline_exon_presence
cat ~/bird_db1/aswin/APOBEC1/TOGA/Struthio_camelus/All_consensus_msa.aln > pipeline_consensus.aln
exmalign -msa4 ../APOBEC1_Calidris_pugnax.fa ../gblast_edited_consensus.fa ../../sblast_edited_consensus.fa ~/bird_db1/aswin/APOBEC1/TOGA/Struthio_camelus/toga_rna-XM_027449866.2_exon_wise_consensus.fa > pipeline_consensus_mafft.aln

#Annotation from TOGA

cp ~/bird_db1/aswin/APOBEC1/TOGA/Struthio_camelus/toga_output/query_* .
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' query_annotation.bed -i
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' query_gene_spans.bed -i

#Expected exon boundaries from toga
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Struthio_camelus/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > toga_exon_metadata.bed
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Struthio_camelus/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > toga_exon_metadata_exp_loc.bed 

#The exon_4 & 5 are completey intact - but other exons are completely degraded including the longest exon_3 which has the conserved motifs
#The N's observed are very short and most probabily occur in intron rather than exon; hence this N's shouldn't be a reason for the inability to recover the exons

#To recover whatever exons that are left : msa consensus varaints 
#the exon 5 has 2 different consensus from different tols. To find which is correct : 
    #1. Check splice sites
    #2. perform an amino acid guided DNA alignment of just exon_4 & 5; NOTE: the alignment should start from a codon not breaking a codon

#Check for inactivating mutations within single intact exon in all possible frames of 5'-3' forward direction
transeq <(cat ~/bird_db1/aswin/APOBEC1/TOGA/Struthio_camelus/toga_rna-XM_027449866.2_exon_wise_consensus.fa | grep -v ">") -auto -stdout -frame F | grep "*" -z --color=always > all_possible_translations_of_intact_exon4_5
#Even if we miss the other exons due to some sequencing error; the inframe inactivating mutations in exon_4 supports the gene must be lost

echo "Struthio_camelus" >> ~/bird_db1/aswin/APOBEC1/manually_validated_losses 
cat ~/bird_db1/aswin/APOBEC1/TOGA/Struthio_camelus/toga_rna-XM_027449866.2_exon_wise_consensus.fa > ~/bird_db1/aswin/APOBEC1/validated_lost_genes/APOBEC1_Struthio_camelus.fa

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Tinamus_guttatus

mkdir ~/bird_db1/Tinamus_guttatus/aswin/APOBEC1/2nd_gblast/manual_search
cd ~/bird_db1/Tinamus_guttatus/aswin/APOBEC1/2nd_gblast/manual_search

#Compare results from pipelines
grep Tinamus ~/bird_db1/aswin/APOBEC1/pipeline_summary_row_wise -A3 > pipeline_exon_presence
cat ~/bird_db1/aswin/APOBEC1/TOGA/Tinamus_guttatus/All_consensus_msa.aln > pipeline_consensus.aln
exmalign2 -msa4 ../APOBEC1_Zonotrichia_albicollis.fa ../gblast_edited_consensus.fa ../../sblast_edited_consensus.fa ~/bird_db1/aswin/APOBEC1/TOGA/Tinamus_guttatus/toga_rna-XM_027449866.2_exon_wise_consensus.fa > pipeline_consensus_mafft.aln

#Annotation from TOGA
cp ~/bird_db1/aswin/APOBEC1/TOGA/Tinamus_guttatus/toga_output/query_* .
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' query_annotation.bed -i
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' query_gene_spans.bed -i

#Expected exon boundaries from toga
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Tinamus_guttatus/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > toga_exon_metadata.bed
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Tinamus_guttatus/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > toga_exon_metadata_exp_loc.bed 
#The gene is towards the end of assembly; Although there is enough length afer the authentic exon_4 hit to find exon_1 to 3 the assembly coverage must have been poor to pick up

#Only exon_4 is initact in genome in the syntenic region; exon_5 is giving hit in another short scaffold lacking other annotations.
#SRA gives hits for exon_1,2,3 with very poor coverage and exon_4 wth high coverage; Exon_5 didn't give any hits; it might be true that exon_1,2,3 are present in genome but the SRA coverage was not enough to assemble in the genome
#To prove the gene is lost we can use the only intact exon_4's all possible ORFs in forward direction
transeq <(grep "exon_4" -A1 ../gblast_edited_consensus.fa) -auto -stdout -frame F | grep "*" -z --color=always > all_possible_translations_of_intact_exon_4
transeq <(grep "exon_4" -A1 ~/bird_db1/aswin/APOBEC1/TOGA/Tinamus_guttatus/toga_rna-XM_027449866.2_exon_wise_consensus.fa) -auto -stdout -frame F | grep "*" -z --color=always  >> all_possible_translations_of_intact_exon_4

#There is STOP codon in all possible frames of exon_4; hence it is lost

echo "Tinamus_guttatus" >> ~/bird_db1/aswin/APOBEC1/manually_validated_losses 
grep "exon_4" -A1 ../gblast_edited_consensus.fa > final_consensus.fa
cp final_consensus.fa ~/bird_db1/aswin/APOBEC1/validated_lost_genes/APOBEC1_Tinamus_guttatus.fa

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Pavo_muticus

mkdir /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Pavo_muticus/aswin/APOBEC1/2nd_gblast/manual_search
cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Pavo_muticus/aswin/APOBEC1/2nd_gblast/manual_search

#Compare results from pipelines
grep Pavo_muticus ~/bird_db1/aswin/APOBEC1/pipeline_summary_row_wise -A3 > pipeline_exon_presence
cat ~/bird_db1/aswin/APOBEC1/TOGA/Pavo_muticus/All_consensus_msa.aln > pipeline_consensus.aln
exmalign2 -msa6 ../APOBEC1_Aythya_fuligula.fa ../gblast_edited_consensus.fa ~/bird_db1/aswin/APOBEC1/TOGA/Pavo_muticus/toga_rna-XM_027449866.2_exon_wise_consensus.fa ../extracted_flanking_region.fa ../../sblast_edited_consensus.fa ../../optimal_query_sblast.outfmt3.fa > pipeline_consensus_mafft.aln

#Annotation from TOGA
cp ~/bird_db1/aswin/APOBEC1/TOGA/Pavo_muticus/toga_output/query_* .
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' query_annotation.bed -i
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' query_gene_spans.bed -i

#Expected exon boundaries from toga
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Pavo_muticus/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > toga_exon_metadata.bed
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Pavo_muticus/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > toga_exon_metadata_exp_loc.bed 

#The gene is in right syntenic location, but location of recovered exon_2,3,4 are same between pipelines, but exon_5 detected only from gblast is in wrong location
#The real assembly gaps are much downstream of expected location of exon_5
 
#The exon_3 boundary is not accurately found hence we can't look for inactivating mutations; but the exon_4 is completely intact

#-----------------------------------------------------
#Verify inactivating mutations in exon_4

transeq <(grep "exon_4" ~/bird_db1/aswin/APOBEC1/TOGA/Pavo_muticus/toga_rna-XM_027449866.2_exon_wise_consensus.fa -A1)  -auto -stdout -frame F | grep "*" -z --color=always
#One frame out of 3 forms an intact exon without any inactivating mutations
#Protein MSA different frames of exon_4 with protein of a validated query : use muscle since it it better with proteins
muscle -in <(cat <( transeq <(grep "exon_4" ~/bird_db1/aswin/APOBEC1/TOGA/Pavo_muticus/toga_rna-XM_027449866.2_exon_wise_consensus.fa -A1) -clean -auto -stdout -frame F) <(cat ~/bird_db1/aswin/APOBEC1/Domain_search/total_validated_queries/total_validated_queries.aa | grep Aythy -A1)) -clw 2>/dev/null
#The frame that doesn't have any inactivating mutations is the right frame

#hence we need to check why TOGA placed an inact mutation there by making ORF using all intact exons(2,3,4 uncluding variants)
egrep "exon_2|exon_4" ~/bird_db1/aswin/APOBEC1/TOGA/Pavo_muticus/toga_rna-XM_027449866.2_exon_wise_consensus.fa -A1 | grep -v "\-\-" > base_exons.fa
grep "exon_3" ../gblast_edited_consensus.fa -A1 | sed '/>/ s/.*/>exon_3_gedit/g' > exon_3_variants.fa
grep "exon_3" ../../sblast_edited_consensus.fa -A1 | sed '/>/ s/.*/>exon_3_sedit/g' >> exon_3_variants.fa
grep "exon_3" ~/bird_db1/aswin/APOBEC1/TOGA/Pavo_muticus/toga_rna-XM_027449866.2_exon_wise_consensus.fa -A1 | sed '/>/ s/.*/>exon_3_toga/g' >> exon_3_variants.fa

#Create variants
counter=0
for j in `grep ">" exon_3_variants.fa | tr -d ">"`
do
counter=$((counter + 1))
echo " - Consensus variant "$counter "=" $j
cat base_exons.fa > consensus_variant_"$counter"_"$j"_tmp.fa
seqtk subseq exon_3_variants.fa <(echo $j) | sed 's/_v[0-9]\+//g' >> consensus_variant_"$counter"_"$j"_tmp.fa
done

#sort exon numbers
for i in `ls consensus_variant_* | sort -V`
do
j=`echo $i | sed 's/_tmp//g'`
echo $j
sed 's/^>/\x00&/' $i | sort -z -V | tr -d '\0'  > $j
rm $i
done

#Choose the correct consensus variant by manual inspection of a iterative amino acid guided alignments
ls ../APOBEC1_Aythya_fuligula.fa consensus_variant_* | xargs -n1 sh -c 'echo $0|sed "s/.fa//g"|sed "s/^/>/g"|sed "s/..\///g" | sed "s/consensus_variant_//g"; grep -v ">" $0 | paste -s -d ""' > all_consensus.fa
#Use aliview and compare nucleotide translated alignment alignment between each consensus with query; when a consensus sequence variant is identified as inaccurate delete that sequence and realign everything again. Repeat this process till find the best match  
#closest variant
cp consensus_variant_2_exon_3_sedit.fa final_consensus.fa

#Since complete exons are missing, best and easiest  way to check gene loss is to look for atleast one completely intact exon possessing atleast one stop codon inside in all 3 forward frames
transeq consensus_variant_2_exon_3_sedit.fa -frame F -auto -stdout | awk -v RS=">" '{$1=$1}1' | sed 's/ //2g' | awk NF | awk '{print$1,gsub(/*/,"*",$2),$2}' | column -t | grep "*" -z
#Exon_3 contains stop codon s in all 3 frames of forward direction 
#hence the gene must be lost

echo "Pavo_muticus" >> ~/bird_db1/aswin/APOBEC1/manually_validated_losses
cp final_consensus.fa ~/bird_db1/aswin/APOBEC1/validated_lost_genes/APOBEC1_Pavo_muticus.fa

################################################################################################################################################################################################################################################################################################################
#Species with missing exons (exon class is N/A)

#Penelope_pileata

mkdir /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Penelope_pileata/aswin/APOBEC1/2nd_gblast/manual_search
cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Penelope_pileata/aswin/APOBEC1/2nd_gblast/manual_search

#Compare results from pipelines
grep Penelope_pileata ~/bird_db1/aswin/APOBEC1/pipeline_summary_row_wise -A3 > pipeline_exon_presence
cat ~/bird_db1/aswin/APOBEC1/TOGA/Penelope_pileata/All_consensus_msa.aln > pipeline_consensus.aln
exmalign2 -msa4 ../APOBEC1_Zonotrichia_albicollis.fa ../gblast_edited_consensus.fa ../../sblast_edited_consensus.fa ~/bird_db1/aswin/APOBEC1/TOGA/Penelope_pileata/toga_rna-XM_027449866.2_exon_wise_consensus.fa > pipeline_consensus_mafft.aln

#Annotation from TOGA
cp ~/bird_db1/aswin/APOBEC1/TOGA/Penelope_pileata/toga_output/query_* .
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' query_annotation.bed -i
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' query_gene_spans.bed -i

#Expected exon boundaries from toga
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Penelope_pileata/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > toga_exon_metadata.bed
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Penelope_pileata/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > toga_exon_metadata_exp_loc.bed 
sed -i '/N\/A/d' toga_exon_metadata_exp_loc.bed

#Extract region where exon_5 is expected
bedtools getfasta -fi ../../../../genome/GCA_013396635.1_ASM1339663v1_genomic.fna -bed <(echo -e "WBMW01004736.1\t39339\t41352")  | grep -v ">" | sed '1i >exon_5_exp_region' > exon_5_exp_region.fa
ls /home/neo/bird_db1/aswin/APOBEC1/validated_sequences/* | grep -v "Limosa" | xargs -n1 sh -c 'echo $0 | awk -F "/" "{print\$NF}" | sed "s/.fa/_exon_5/g" | sed "s/^/>/g"; grep -A1 "exon_5" $0 | grep -v ">"' > all_validated_queries_exon_5.fa
cat exon_5_exp_region.fa all_validated_queries_exon_5.fa > Penelope_exon_5_exp_region_and_all_validated_queries_exon_5.fa
#Best way to find short and diverged sequenes is using a hmm profile bed search  
mafft --maxiterate 1000 --auto --reorder --quiet all_validated_queries_exon_5.fa > all_validated_queries_exon_5.aln
hmmbuild --dna all_validated_queries_exon_5.hmm all_validated_queries_exon_5.aln
hmmsearch --max all_validated_queries_exon_5.hmm exon_5_exp_region.fa
hmmpress all_validated_queries_exon_5.hmm
hmmscan --max all_validated_queries_exon_5.hmm exon_5_exp_region.fa

#There are no reliable hits for exon_5. Since this exon contributes to the last ~5.6% of the gene. Even if this exon is truly lost ~95% of the rest of the gene stays intact
#But the inframe stop codons in previous exons can prove the loss without finding exon_5 given that the previous exons are devoid of assembly and alignment ambiguities
#check for assembly quality using SRA (coverage must not be lower than 10 bases and no variants/SNP's)
awk '$2<10' ../../optimal_query_sblast.outfmt3_processed_data/variants
#check ORF from 1st 4 exons
transeq ~/bird_db1/aswin/APOBEC1/TOGA/Penelope_pileata/toga_rna-XM_027449866.2_consensus.fa -auto -stdout | grep "*" -z
#NOTE: Placement of stop codons are different in TOGA and direct translation using transeq. But in either way there is a stop codon in exon_3
#Best aligner for acounting framshifts is MACSE
java -jar ~/programmes/macse_v2.06.jar -prog alignSequences -seq /home/neo/soft_links/Penelope_pileata/aswin/APOBEC1/2nd_gblast/manual_search/query_and_final_consensus.fa
#To exclude the possibility of having a frameshift in exon_3, check for ORF in all 3 frames in forward direction
transeq -frame F <(grep "exon_3" ~/bird_db1/aswin/APOBEC1/TOGA/Penelope_pileata/toga_rna-XM_027449866.2_exon_wise_consensus.fa -A1) -auto -stdout | awk -v RS=">" '{$1=$1}1' | sed 's/ //2g' | awk NF | awk '{print$1,gsub(/*/,"*",$2),$2}' | column -t | grep "*" -z
#there is a stop codon in all 3 frames of exon_3  

#Conclusion : the gene is lost
echo "Penelope_pileata" >> ~/bird_db1/aswin/APOBEC1/manually_validated_losses
cat ~/bird_db1/aswin/APOBEC1/TOGA/Penelope_pileata/toga_rna-XM_027449866.2_exon_wise_consensus.fa > ~/bird_db1/aswin/APOBEC1/validated_lost_genes/APOBEC1_Penelope_pileata.fa

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Meleagris_gallopavo

mkdir /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Meleagris_gallopavo/aswin/APOBEC1/2nd_gblast/manual_search
cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Meleagris_gallopavo/aswin/APOBEC1/2nd_gblast/manual_search

#Compare results from pipelines
grep Meleagris_gallopavo ~/bird_db1/aswin/APOBEC1/pipeline_summary_row_wise -A3 > pipeline_exon_presence
cat ~/bird_db1/aswin/APOBEC1/TOGA/Meleagris_gallopavo/All_consensus_msa.aln > pipeline_consensus.aln
exmalign2 -msa4 ../APOBEC1_Aythya_fuligula.fa ../gblast_edited_consensus.fa ~/bird_db1/aswin/APOBEC1/TOGA/Meleagris_gallopavo/toga_rna-XM_027449866.2_exon_wise_consensus.fa ../extracted_flanking_region.fa > pipeline_consensus_mafft.aln

#Annotation from TOGA
cp ~/bird_db1/aswin/APOBEC1/TOGA/Meleagris_gallopavo/toga_output/query_* .
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' query_annotation.bed -i
sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' query_gene_spans.bed -i

#Expected exon boundaries from toga
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Meleagris_gallopavo/toga_output/temp/exons_meta_data.tsv | awk '{print$4,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{$1=$1}1' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > toga_exon_metadata.bed
grep "rna-XM_027449866.2" ~/bird_db1/aswin/APOBEC1/TOGA/Meleagris_gallopavo/toga_output/temp/exons_meta_data.tsv | awk '{print$4,$5,"exon_"++i,"1","+"}' | tr ":-" " " | awk '{print$1,$5,$6,$7,$8,$9}' OFS="\t" | sed 's/^\(.*[0-9]\)_\(.*\)/\1.\2/g' > toga_exon_metadata_exp_loc.bed 
sed -i '/N\/A/d' toga_exon_metadata_exp_loc.bed

#From TOGA : exon_2 is very unreliable hit and exon_3 start is also missing, but TOGA shows inactivating mutations in exon_3. These mutations might have been detected based CESAR codon alignment
#However SRA shows zero support for exon_3 start that might suggest that the assembly is poor at this region 
#Manual inspection shows that out of the 3 inactivating mutations detected by TOGA only one mutation (2bp deletion in exon_3) have SRA support other have zero support

#NOTE : the uncertainity in genome remains
#In this uncertainity the best way to find loss is checking all 3 reading frames of exon_3 
transeq -frame F <(grep "exon_3" ~/bird_db1/aswin/APOBEC1/TOGA/Meleagris_gallopavo/toga_rna-XM_027449866.2_exon_wise_consensus.fa -A1) -auto -stdout | awk -v RS=">" '{$1=$1}1' | sed 's/ //2g' | awk NF | awk '{print$1,gsub(/*/,"*",$2),$2}' | column -t | grep "*" -z

echo "Meleagris_gallopavo" >> manually_validated_losses
awk '/exon_3|exon_4/ {print">"$0}' RS=">" ~/bird_db1/aswin/APOBEC1/TOGA/Meleagris_gallopavo/toga_rna-XM_027449866.2_exon_wise_consensus.fa | awk NF > ~/bird_db1/aswin/APOBEC1/validated_lost_genes/APOBEC1_Meleagris_gallopavo.fa

################################################################################################################################################################################################################################################################################################################
#DRAFT

counter=31
for i in `grep -v ">" ~/bird_db1/aswin/APOBEC1/TOGA/Cathartes_aura/toga_rna-XM_027449866.2_exon_wise_consensus.fa`
do 
counter=$((counter + 1))
echo "GREP_COLORS=\"mt=01;"$counter"\" grep -i "$i" --color=always -z"
done | paste -s -d "|" | sed 's!^!cat manual_search/pairwise_query_toga.aln | !g'








