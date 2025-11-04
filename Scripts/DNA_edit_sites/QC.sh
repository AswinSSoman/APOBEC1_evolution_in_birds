
###########################################################################################################################################################################################################################################################################################################
#Quality check of genomic data

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#extarct latest genome accession for the each bird species

time cat all_bird_genomes_used | cut -f1,2 -d "_" | xargs -n1 sh -c 'esearch -db genome -query "$0 [ORGN]" | esummary | xtract -pattern DocumentSummary -element Organism_Name Assembly_Accession' > results/QC/latest_version_of_all_bird_genomes_used
sed -e 's/ /_/g' -e 's/Chloebia_gouldiae/Erythrura_gouldiae/g' -e 's/Cyanoderma_ruficeps/Stachyris_ruficeps/g' results/QC/latest_version_of_all_bird_genomes_used -i

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Compare local genome version with latest version, to check if a significantly better genome is available

while read -r i
do
cd $i
#in=`echo $i | cut -f1,2 -d "_"`
if [[ -d ucsc ]]; then cd ucsc; genome=$(ls GC*.fa | cut -f1,2 -d "_" | egrep -v "rmsk|repeatModeler" | sed 's/\.fa//g'); rmsk=$(ls *.repeatMasker.out); else genome=$(ls GC*.fna | cut -f1,2 -d "_"); cd repeat_masker_using_"$i"_library; rmsk=$(ls *.fna.out); fi
latest_genome=$(grep -if <(echo "$i" | cut -f1,2 -d "_") <(sed 's/ /_/g' /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/QC/latest_version_of_all_bird_genomes_used) | awk '{print$2}')
if [[ "$rmsk" == *"$genome"* ]]; then rg="same"; else rg="diff"; fi
if [[ $(echo $latest_genome | cut -f2 -d "_") == $(echo $genome | cut -f2 -d "_") ]]; then lg="same"; else lg="diff"; fi
echo $i $genome $rmsk $latest_genome $rg $lg
unset genome rmsk latest_genome rg lg
cd $base_dir
done < all_bird_genomes_used | sed '1i Species Genome(G) RepeatMasker(R) Latest_genome(lG) compare_GR compare_GlG' > results/QC/genome_version_compare_local_vs_latest_local_vs_latest


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fetch genome metadata is the genome versions differ (If genome versions are same skip) 

time for i in $(awk '/diff/ {print$1}' results/QC/genome_version_compare_local_vs_latest)
do
echo ">"$i
esearch -db assembly -query "$i [ORGN]" | esummary | xtract -pattern DocumentSummary -def "-" \
-element SpeciesName Sub_value UCSCName AssemblyName AssemblyDescription RefSeq_category AssemblyAccession Genbank RefSeq LastMajorReleaseAccession LatestAccession RefSeqAnnotationRelease AssemblyStatus Coverage AsmReleaseDate_GenBank LastUpdateDate ContigN50 ScaffoldN50 \
| sed '1i Species_Name\tBreed\tUCSC_Name\tAssembly_Name\tDescription\tCategory\tAssembly_Accession\tGenbank\tRefSeq\tLast_Major_Release\tLatest_Accession\tAnnotation\tAssembly_level\tCoverage\tRelease_Date\tlast_update\tContig_N50\tScaffold_N50\n' | sed 's/ /_/g'
done > results/QC/check_metadata_of_species_with_diff_genome_versions

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Check which genome version is better in quality (based on coverage, scaffoldN50, ContigN50)

while read -r i
do
i2=$(echo $i | awk '{print$2}')
i4=$(echo $i | awk '{print$4}')
echo ">"$i | awk '{$3=""; print$0}' | GREP_COLORS="mt=01;33" grep "$i2" --color=always | GREP_COLORS="mt=01;36" grep "$i4" -a --color=always
echo -e
GREP_COLORS="mt=01;33" grep "$i2" results/QC/check_metadata_of_species_with_diff_genome_versions -az --color=always | GREP_COLORS="mt=01;36" grep "$i4" -az --color=always | awk '{$1="-"; print$0}' | egrep "$i2|$i4" -a 
#echo ---$a{1..19}
#echo ---$a{1..19}
unset i2 i4
done < <(grep "diff" results/QC/genome_version_compare_local_vs_latest) | sed "1i $(sed -n '2p' results/QC/check_metadata_of_species_with_diff_genome_versions)" | sed '1i Species_Name Local_Genome Latest_Genome' | column -t | sed "s/^>/$(printf '_%.s' {1..380})\n/g" > results/QC/manual_check_4_better_genome_version

#########-----------------------#########-----------------------#########-----------------------##########-----------------------##########-----------------------##########-----------------------##########-----------------------##########-----------------------##########-----------------------#
#Observations:

#These genomes are found to have a better version in database than what was locally dowloaded 
  #1. Colius_striatus - higher sN50, similar cN50 & coverage
  #2. Passer_domesticus - higher sN50 & cN50, lower coverage
  #3. Chrysolophus_pictus - higher sN50 & cN50, lower coverage
  #4. Pavo_cristatus - higher sN50, cN50 & coverage

#########-----------------------#########-----------------------#########-----------------------##########-----------------------##########-----------------------##########-----------------------##########-----------------------##########-----------------------##########-----------------------#
#Check LTR content & number of edit sites after each filtering steps. the final edit site count & relax selection results for each species

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method

#RNA expression data
cat RNA_expression | sed 's/ /_/g' | awk '{if($0~"Data_not_checked") print$0,"-","-"; else print$1,$(NF-2),$(NF-1)}' | column -t | less -SN

#Relax selection results
scp neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/relaxed_selection/amniotes_as_bg/12_amniotes_as_bg/one_birds_as_fg_all_amniotes_as_bg/relax_all_models/relax_summary relax_all_models_summary
scp neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/relaxed_selection/amniotes_as_bg/12_amniotes_as_bg/one_birds_as_fg_all_amniotes_as_bg/relax_all_models_srv/relax_summary relax_all_models_srv_summary
scp neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/relaxed_selection/amniotes_as_bg/12_amniotes_as_bg/one_birds_as_fg_all_amniotes_as_bg/relax_summary relax_minimal_models_summary
scp neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/relaxed_selection/amniotes_as_bg/12_amniotes_as_bg/one_birds_as_fg_all_amniotes_as_bg/relax_minimal_model_srv_multiple_hits/relax_summary relax_minimal_models_srv_summary

while read i
do
i1=$(echo $i | awk '{print$1}')  #species name
i2=$(echo $i | awk '{print$2}')  #gene loss info
i3=$(grep $i1 /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/QC/size_of_each_repeatmasker_identified_transposon_type_per_genome | awk '{print$5}')  #LTR size
i4=$(awk -v c=1 'BEGIN {FS=OFS=" "} $c ~ /\-/ {$c = p} {p = $c} 1' glimpse_output_of_knisbacher_script_step_5_create_filtered_tracks_GA_81_birds | grep $i1 | awk '{for (i=1; i<=NF; i++) sum[i]+=$i} 1; END {for (i=1; i<=NF; i++) printf sum[i] " "; print ""}' | awk 'END{print$3,$4,$5,$6,$7}')  #cluster count after each clusters
i5=$(echo $i | awk '{$1=""; $2=""; print$0}')  #GA & CT count
i6=$(echo $i | awk '{if($3=="0" && $4=="0") print"- -"; else if($3==0 && $4>0) print"0 100"; else if($4==0 && $3>0) print"100 0"; else if($3>0 && $4>0) print($3/($3+$4))*100,($4/($3+$4))*100}')  #GE & CT edit %
j1=$(echo $i | awk '{if($3=="0" && $4=="0") print"-"; else if($4==0 && $3>0) print"-"; else print$3/$4}')  #GA/CT
j2=$(echo $i | awk -v a="$i3" '{if(a=="0") print"-"; else print$3/a}')  #GA/LTR bp
if [[ -z $j2 ]]; then j2="-"; else :; fi
j2a=$(echo $i | awk -v a="$i3" '{print$3/a}' | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}')
#j3=$(echo "l("$j2a")/l(10)" | bc -l)  #Log10(GA/LTR bp) in decimals
j3=$(echo "l("$j2a")/l(10)" | bc -l | xargs -n1 sh -c 'printf "%.3e" $0')  #Log10(GA/LTR bp) in scietific notation & rounded hence not very accurate
if [[ -z $j3 ]]; then j3="-"; else :; fi
i7=$(grep $i1 relax_all_models_summary | awk '{print$4}')
i8=$(grep $i1 relax_all_models_srv_summary | awk '{print$4}')
i9=$(grep $i1 relax_minimal_models_summary | awk '{print$4}')
i10=$(grep $i1 relax_minimal_models_srv_summary | awk '{print$4}')
if [[ -z $i7 ]]; then i7="-"; else :; fi
if [[ -z $i8 ]]; then i8="-"; else :; fi
if [[ -z $i9 ]]; then i9="-"; else :; fi
if [[ -z $i10 ]]; then i10="-"; else :; fi
echo $i1 $i2 $i3 $i4 $i5 $i6 $j1 $j2 $j3 $i7 $i8 $i9 $i10
unset i1 i2 i3 i4 i5 i6 j1 j2 j2a j3 i7 i8 i9 i10
done < <(awk 'NR>1' birds_total_GA_CT_edit_sites_count_with_loss_status) | sed '1i Species Loss_info LTR_size Unfiltered_GA_clusters Pairwise_filtered_GA_clusters Divergence_filtered_GA_clusters Map_to_G_than_A_filtered_GA_clusters Best_sources_filtered_GA_clusters GA_edit_sites CT_edit_sites GA_%_edit CT_%_edit GA/CT GA/LTR_bp log10(GA/LTR_bp) relax_all_models relax_all_models_srv relax_min_models relax_min_models_srv' | column -t > birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax

#min=$(awk 'NR>1{print$3}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax | sort -n | awk '$1>0' | head -1)
#sed -e 's/_clusters//g' -e 's/_filtered//g' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax | awk -v a="$min" '{print$1,$2,$3,$3/a,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' | awk '{$4+=0}1' CONVFMT="%.0f" | awk 'NR==1{$4="LTR_size_norm"; $1=$1}1' | sort -k2,2n -k3,3n | column -t

#########-----------------------#########-----------------------#########-----------------------##########-----------------------##########-----------------------##########-----------------------##########-----------------------##########-----------------------##########-----------------------#
#Normalize Output (GA edit sites) considering negative control (CT edit sites) & Input (Total LTR sequence size)

#Use different normaliztion methods
  #1. GA_CT_Ratio_Length_Adjusted: normalized GA sites relative to CT sites, scaled by the sequence length.
  #2. GA_Log_Normalized: If the distribution of GA edit sites is highly variable, you could use a log transformation approach to reduce skew, normalizing by both CT sites and total LTR size
                       #: The +1 ensures that values with zero GA edit sites can still be handled in the log transformation.
                       #: This approach allows you to control for variance by using both CT and LTR size together as normalizing factors, offering a scaled and normalized output.
  #3. Composite_Score: combines the effects of both the CT sites and LTR sequence size in a single formula.
                    #: uses the geometric mean of CT sites and LTR size to avoid any single variable from dominating the normalization process.
                    #: The geometric mean gives a balanced scaling factor that considers both variables equally.
  #4. Z_Score: if you have a large dataset across many bird species, you could standardize GA edit sites across species by calculating z-scores based on the CT and LTR values
            #: Here, 𝜇 and 𝜎 represent the mean and standard deviation of the adjusted values.
            #: This approach gives a normalized, unit-free score that accounts for variability across species.

awk 'BEGIN { OFS="\t" } 
NR == 1 {
    # Print header line with new column names
    print $0, "GA_CT_Ratio_Length_Adjusted", "GA_Log_Normalized", "Composite_Score"
    next
}
NR > 1 {
    # Read in values
    LTR_Size = $3
    GA_Edit_Sites = $9
    CT_Edit_Sites = $10
    # 1. GA_CT_Ratio_Length_Adjusted: ((GA edit sites) / (CT edit sites)) / LTR size
    if (CT_Edit_Sites > 0 && LTR_Size > 0) {
        GA_CT_Ratio_Length_Adjusted = (GA_Edit_Sites / CT_Edit_Sites) / LTR_Size
    } else {
        GA_CT_Ratio_Length_Adjusted = "NA"
    }
    # 2. GA_Log_Normalized: log10((GA edit sites + 1) / (CT edit sites * LTR size))
    if (CT_Edit_Sites > 0 && LTR_Size > 0) {
        GA_Log_Normalized = log((GA_Edit_Sites + 1) / (CT_Edit_Sites * LTR_Size)) / log(10)
    } else {
        GA_Log_Normalized = "NA"
    }
    # 3. Composite_Score: (GA edit sites) / sqrt(CT edit sites * LTR size)
    if (CT_Edit_Sites > 0 && LTR_Size > 0) {
        Composite_Score = GA_Edit_Sites / sqrt(CT_Edit_Sites * LTR_Size)
    } else {
        Composite_Score = "NA"
    }
    # Print row with new calculations (excluding Z-Score)
    print $0, GA_CT_Ratio_Length_Adjusted, GA_Log_Normalized, Composite_Score
}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax > intermediate_table.txt

# Calculate mean of Composite_Score
mean=$(awk 'NR>1 && $22!="NA" {sum+=$22; count++} END{if(count>0) print sum/count}' intermediate_table.txt)

# Calculate standard deviation of Composite_Score
stdev=$(awk -v mean=$mean 'NR>1 && $22!="NA" {sum_sq+=($22-mean)^2; count++} END{if(count>1) print sqrt(sum_sq/(count-1))}' intermediate_table.txt)

# Use mean and stdev to calculate Z_Score and output final table
awk -v mean=$mean -v stdev=$stdev 'BEGIN { OFS="\t" } 
NR == 1 {
    # Print header line with Z_Score column added
    print $0, "Z_Score"
    next
}
NR > 1 {
    # Read in Composite_Score
    Composite_Score = $8
    
    # Calculate Z-Score if Composite_Score is not "NA" and stdev > 0
    if (Composite_Score != "NA" && stdev > 0) {
        Z_Score = (Composite_Score - mean) / stdev
    } else {
        Z_Score = "NA"
    }

    # Print row with calculated Z_Score
    print $0, Z_Score
}' intermediate_table.txt | column -t > birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized

rm intermediate_table.txt birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax
unset mean stdev

#Shorten titles view in single frame
grep -v "^Species" birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized \
| sed "1i - Loss - Unfiltered Pairwise Divergence Map_to_G_than_A Best_sources #GA_edit #CT_edit - - - - log relax_all relax_all relax_min relax_min GA_CT_Ratio GA_Log Composite - \n \
Species info LTR_size GA_clusters filtered filtered filtered filtered sites sites GA_%_edit CT_%_edit GA/CT GA/LTR_bp (GA/LTR_bp) models models_srv models models_srv Length_Adjusted Normalized Score Z_score" | column -t \
| GREP_COLORS="mt=01;32" grep "^-.*\|$" --color=always | GREP_COLORS="mt=01;04;32" grep "^Species.*\|$" --color=always | less -SR

#Print Species with A1 loss having higher GA % edit than CT % edit & Species with intact A1 having lower GA % edit than CT % edit
#NOTE: these species must be having very low number of overall edit site count for both GA & CT edit sites
grep -v "^Species" birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized | sed "1i - Loss - Unfiltered Pairwise Divergence Map_to_G_than_A Best_sources #GA_edit #CT_edit - - - - log relax_all relax_all relax_min relax_min GA_CT_Ratio GA_Log Composite - \n \
Species info LTR_size GA_clusters filtered filtered filtered filtered sites sites GA_%_edit CT_%_edit GA/CT GA/LTR_bp (GA/LTR_bp) models models_srv models models_srv Length_Adjusted Normalized Score Z_score" | column -t | GREP_COLORS="mt=01;32" grep "^-.*\|$" --color=always | GREP_COLORS="mt=01;04;32" grep "^Species.*\|$" --color=always  | awk '{if($2=="Loss" || $2=="info") print$0; else if(($2==0) && ($11>$12)) print$0; else if($2=="1" && $11<$12) print$0}'  | less -SR

#Add bird order info
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method
scp neo@172.30.1.174:~/bird_db1/aswin/taxonomy/orders_all_birds .
sed 's/Erythrura_gouldiae/Chloebia_gouldiae/g' orders_all_birds -i
sed 's/Stachyris_ruficeps/Cyanoderma_ruficeps/g' orders_all_birds -i

for i in $(awk '{print$1}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed | grep -v "Species")
do
i1=$(grep "$i" /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/orders_all_birds | awk '{print$2}')
i2=$(grep "$i" birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed | awk '!($1="")')
echo $i $i1 $i2
unset i1 i2
done | sed '1i Species Order Loss_info LTR_size Unfiltered_GA_clusters Pairwise_filtered_GA_clusters Divergence_filtered_GA_clusters Map_to_G_than_A_filtered_GA_clusters Best_sources_filtered_GA_clusters GA_edit_sites CT_edit_sites GA_%_edit CT_%_edit GA/CT GA/LTR_bp log10(GA/LTR_bp) relax_all_models relax_all_models_srv relax_min_models relax_min_models_srv GA_CT_Ratio_Length_Adjusted GA_Log_Normalized Composite_Score Z_Score' > birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed_with_bird_orders
awk '{print$1,$2,$3,$4,$10}' ../results/knisbacher_method/birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed_with_bird_orders | sort -k5,5nr | column -t | grep "Passeriformes" -z

#NOTE: The GA edit sites are higher in Passeroidea within passeriformes, but not in ther passeriformes

###########################################################################################################################################################################################################################################################################################################
#Print DNA edit sites results with ERV total & edited sizes

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons
for i in $(cat all_bird_genomes_used)
do
#Species name suffix
o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
ksum=$(readlink -f results/knisbacher_method/birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed_with_bird_orders)
i1=$(echo $i | cut -f1,2 -d "_" | sed 's/Stachyris_ruficeps/Cyanoderma_ruficeps/g' | sed 's/Erythrura_gouldiae/Chloebia_gouldiae/g')
#Bird order
order=$(grep $i1 $ksum | awk '{print$2}')
#Gene loss/intact info
gli=$(grep $i1 $ksum | awk '{print$3}')
#GA edit sites
gaes=$(grep $i1 $ksum | awk '{print$10}')
#CT edit sites
ctes=$(grep $i1 $ksum | awk '{print$10}')
#GA log normalized
ganorm=$(grep $i1 $ksum | awk '{print$22}')
#Relax results (min models)
relax=$(grep $i1 $ksum | awk '{print$19}')
#LTR size
ltrsize=$(grep $i1 $ksum | awk '{print$4}')
#Calculate only if LTRs are present in genome
if [[ -d $i/knisbacher/Data/$o/LTR/db ]]
then
cd $i/knisbacher/Data/$o/LTR/db
#ERV length
ervlength=$(for erv in $(find . -name "Len_*.txt" | grep -i "erv")
do
awk '{a=0; for(n=2;n<=NF;n++) a+=$n; print $1, a}' $erv
done | awk '{a+=$2;} END{print a;}')
if [[ -z $ervlength ]]; then ervlength="0"; else :; fi
#Total subfamily names
subfam=$(ls -d */ | tr -d "/" | sed 's/files_//g' | awk NF | paste -s -d",")
if [[ -z $subfam ]]; then subfam="0"; else :; fi
ervsubfam=$(ls -d */ | tr -d "/" | sed 's/files_//g' | awk NF | grep -i "erv" | paste -s -d",")
if [[ -z $ervsubfam ]]; then ervsubfam="0"; else :; fi
cd ../../../../
#Total size of edited targets
bpc_totalsize=$(for bpc in $(find Data/$o/LTR/results -maxdepth 5 -mindepth 5 -path "*/Tracks/*/GA/pairwise_filter*" -name "bestPairsClusters_*.tab")
do
bpc_subfam=$(echo $bpc | cut -f6 -d "/" | sed 's/.*_LTR_//g' | sed 's/_1e\-0.*//g')
totalsize_bpc_subfam=$(awk '{print$7}' $bpc | cut -f2 -d ":" | sed 's/\-$//g' | sed 's/+$//g' | awk -F "-" '{print$2-$1}' | awk '{a+=$1;} END{print a;}')
echo $bpc_subfam"("$totalsize_bpc_subfam")"
unset bpc_subfam totalsize_bpc_subfam
done | paste -s -d ",")
if [[ -z $bpc_totalsize ]]; then bpc_totalsize="0"; else :; fi
combined_bpc_totalsize=$(echo $bpc_totalsize | tr "," "\n" | tr "(" " " | cut -f2 -d " " | tr -d ")" | awk '{a+=$1;} END{print a}')
if [[ -z $combined_bpc_totalsize ]]; then combined_bpc_totalsize="0"; else :; fi
echo $i $order $gli $gaes $ctes $ganorm $relax $ltrsize $ervlength $combined_bpc_totalsize $subfam $ervsubfam $bpc_totalsize
else
echo $i $order $gli $gaes $ctes $ganorm $relax $ltrsize "- - - - -"
fi
unset o i1 ksum order gli gaes ctes ganorm relax ltrsize ervlength subfam ervsubfam bpc_totalsize combined_bpc_totalsize
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons
done | sed '1i Species Order Loss/Intact #GA_edit_sites #CT_edit_sites GA_Log_norm Relaxed_selection LTR_size ERV_size Edited_ERV_size Total_subfam ERV_subfam Edited_ERV_subfam' | column -t > results/knisbacher_method/summary_GA_edit_sites_ERV_size

###########################################################################################################################################################################################################################################################################################################

#Correlation plot

#Check names not matching b/w table & tree
grep -vif  <(cat list_of_species_to_keep_upupa_epops_removed.nwk | tr -d "():[0-9]'.;" | tr "," "\n" | sort) <(awk '{print$1}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized | sort)
grep -wvif  <(cat list_of_species_to_keep_upupa_epops_removed.nwk | tr -d "():[0-9]'.;" | tr "," "\n" | sort) <(awk '{print$1}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized | sort) | awk -F "_" '$3' \
    | awk -F "_" '{print"s/"$0"/"$1"_"$2"/g"}' | sed "s/^/-e '/g" | sed "s/$/'/g" | paste -s -d " " | sed 's/^/sed /g'

#Rename names
sed -e 's/Erythrura_gouldiae/Chloebia_gouldiae/g' -e 's/Stachyris_ruficeps/Cyanoderma_ruficeps/g' -e 's/Aquila_chrysaetos_chrysaetos/Aquila_chrysaetos/g' -e 's/Buceros_rhinoceros_silvestris/Buceros_rhinoceros/g' -e 's/Limosa_lapponica_baueri/Limosa_lapponica/g' -e 's/Lonchura_striata_domestica/Lonchura_striata/g' -e 's/Malurus_cyaneus_samueli/Malurus_cyaneus/g' -e 's/Struthio_camelus_australis/Struthio_camelus/g' list_of_species_to_keep_upupa_epops_removed.nwk > correlation_plot/list_of_species_to_keep_upupa_epops_removed_renamed.nwk

sed -e 's/Erythrura_gouldiae/Chloebia_gouldiae/g' -e 's/Stachyris_ruficeps/Cyanoderma_ruficeps/g' -e 's/Aquila_chrysaetos_chrysaetos/Aquila_chrysaetos/g' -e 's/Buceros_rhinoceros_silvestris/Buceros_rhinoceros/g' -e 's/Limosa_lapponica_baueri/Limosa_lapponica/g' -e 's/Lonchura_striata_domestica/Lonchura_striata/g' -e 's/Malurus_cyaneus_samueli/Malurus_cyaneus/g' -e 's/Struthio_camelus_australis/Struthio_camelus/g' \
birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized | sed -e '/^Species/ s/%/percent/g' -e '/^Species/ s!/!_by_!g' -e '/^Species/ s/log10(/log10_/g' -e '/^Species/ s/)//g' | sed 's/[ ]\+/\t/g' > correlation_plot/birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed

#Print how the normalized GA edit count is distributed for each normalization method, NOTE: Normalization method without extreme outliers & with 2 peaks are preferred 
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/correlation_plot
for i in GA_percent_edit CT_percent_edit GA_by_CT GA_by_LTR_bp GA_CT_Ratio_Length_Adjusted GA_Log_Normalized Composite_Score Z_Score
do
echo " "$i | GREP_COLORS="mt=01;07;33" grep ".*" --color=always
awk -v h="$i" 'NR==1{for(i=1;i<=NF;++i) if($i==h) {n=i;break}} {print$n}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed | egrep -v "$i|NA|^\-$" \
  | ministat -w 150 | grep -v "stdin" | GREP_COLORS="mt=01;34" egrep "^x.*|$|.*Stddev.*|$" --color=always |  sed 's/^/   /g' && printf "\n"
done | less -SR

#Prune tree for each normalization method: each method has different normalized values, & if the normalized count of any species is "NA" then that species should be skipped & pruned correpondingly in tree 
start_time=$(date +%s)
for i in GA_percent_edit CT_percent_edit GA_by_CT GA_by_LTR_bp GA_CT_Ratio_Length_Adjusted GA_Log_Normalized Composite_Score Z_Score
do
{
  echo ">"$i
  awk -v h="$i" 'NR==1{for(i=1;i<=NF;i++) {f[$i]=i}} {print$(f["Species"]),$(f["Loss_info"]),$(f[h]) }' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed | awk '$3!="NA"' | awk '$3!="-"' > $i"_filtered"
  awk '{print$1}' $i"_filtered" | grep -v "^Species" | paste -s -d "," > $i"_species_to_keep"
  Rscript /media/aswin/programs/prune_tree.r list_of_species_to_keep_upupa_epops_removed_renamed.nwk $i"_species_to_keep" $i"_filtered".nwk
  sed "s/)'[0-9]\+'/)/g" $i"_filtered".nwk -i
  Rscript test_normalization_method_correlation.R $i"_filtered" $i"_filtered".nwk $i &
}
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

###########################################################################################################################################################################################################################################################################################################
#Add these column next to tree

sed 's/Erythrura_gouldiae/Chloebia_gouldiae/g' ordered_list_of_species_to_keep_upupa_epops_removed -i
sed 's/Stachyris_ruficeps/Cyanoderma_ruficeps/g' ordered_list_of_species_to_keep_upupa_epops_removed -i
cat <(head -1 birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed) <(cat ordered_list_of_species_to_keep_upupa_epops_removed | cut -f1,2 -d "_" | xargs -n1 sh -c 'grep "$0" birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed') | awk '{print$1,$3,$9,$10,$21,$17}' | awk '{if($NF~"NS") print$1,$2,$3,$4,$5,"Non-significant"; else print$0}'

###########################################################################################################################################################################################################################################################################################################
#DRAFT SCRIPTS
###########################################################################################################################################################################################################################################################################################################

#Check differences in species names b/w tree & data table
sed -e 's/_clusters//g' -e 's/_filtered//g' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax | awk '{print$1,$2,$3,$3/23328,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' | awk '{$4+=0}1' CONVFMT="%.0f" | awk 'NR==1{$4="LTR_size_norm"; $1=$1}1' | sort -k2,2n -k3,3n | column -t
sed -e 's/_clusters//g' -e 's/_filtered//g' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_relax | awk '{print$1,$2,$3,$3/23328,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | awk '{$4+=0}1' CONVFMT="%.0f" | awk 'NR==1{$4="LTR_size_norm"; $1=$1}1' | sort -k2,2n -k3,3n | column -t
awk '{print$1,$2,$3,$3/23328,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings | q | column -t | less -SN

#Print Species with A1 loss having higher GA % edit than CT % edit & Species with intact A1 having lower GA % edit than CT % edit
awk '{if($2=="Loss_info") print$0; else if(($2==0) && ($11>$12)) print$0; else if($2=="1" && $11<$12) print$0}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized | less -S

awk 'NR>1{print$20}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized | grep -v "NA" | ministat -w 180
awk 'NR>1{print$21}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized | grep -v "NA" | ministat -w 180
awk 'NR>1{print$22}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized | grep -v "NA" | ministat -w 180
awk 'NR>1{print$23}' birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized | grep -v "NA" | ministat -w 180

scp list_of_species_to_keep_upupa_epops_removed_renamed.nwk birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed neo@172.30.1.174:~/bird_db1/aswin/APOBEC1/function/DNA_edit/correlation_plot/

Rscript test_normalization_method_correlation.R GA_percent_edit_filtered list_of_species_to_keep_upupa_epops_removed_renamed.nwk GA_percent_edit_filtered.nwk GA_percent_edit 100 4
Rscript test_normalization_method_correlation.R CT_percent_edit_filtered list_of_species_to_keep_upupa_epops_removed_renamed.nwk CT_percent_edit_filtered.nwk CT_percent_edit 10 4
Rscript test_normalization_method_correlation.R GA_by_CT_filtered list_of_species_to_keep_upupa_epops_removed_renamed.nwk GA_by_CT_filtered.nwk GA_by_CT
Rscript test_normalization_method_correlation.R GA_by_LTR_bp_filtered list_of_species_to_keep_upupa_epops_removed_renamed.nwk GA_by_LTR_bp_filtered.nwk GA_by_LTR_bp
Rscript test_normalization_method_correlation.R GA_CT_Ratio_Length_Adjusted_filtered list_of_species_to_keep_upupa_epops_removed_renamed.nwk GA_CT_Ratio_Length_Adjusted_filtered.nwk GA_CT_Ratio_Length_Adjusted
Rscript test_normalization_method_correlation.R GA_Log_Normalized_filtered list_of_species_to_keep_upupa_epops_removed_renamed.nwk GA_Log_Normalized_filtered.nwk GA_Log_Normalized
Rscript test_normalization_method_correlation.R Composite_Score_filtered list_of_species_to_keep_upupa_epops_removed_renamed.nwk Composite_Score_filtered.nwk Composite_Score
Rscript test_normalization_method_correlation.R Z_Score_filtered list_of_species_to_keep_upupa_epops_removed_renamed.nwk Z_Score_filtered.nwk Z_Score


ls -rt | egrep "_filtered$|_filtered.nwk" | paste -d " " - - | awk '{print$1,$2,$2}' | sed 's/_filtered.nwk$//g' | awk '{print"a"++i"=read.table(\"Downloads/APOBEC1/correlation_plot/"$1"\",header=T,row.names=1)""\n""t"++j"=read.tree(\"Downloads/APOBEC1/correlation_plot/"$2"\")""\n""t(table(a"++k"$Loss_info,a"++l"$"$3"))->M"++m}'
ls *_filtered -rt | sed 's/_filtered$//g' | xargs -n1 sh -c 'ls Regression_R_data_* | grep "$0"' | awk '{print"load(\"Downloads/APOBEC1/correlation_plot/"$1"\", c"++i "<-new.env())"}'

#END

