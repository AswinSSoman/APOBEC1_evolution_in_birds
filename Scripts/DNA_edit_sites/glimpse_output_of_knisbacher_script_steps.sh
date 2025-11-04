#!/bin/bash

###########################################################################################################################################################################################################################################################################################################

#Function: Quality Check: Quick check the main parts of output to ensure there are no errors

###########################################################################################################################################################################################################################################################################################################

#Set path variables
base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"
script_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/scripts/Perl_DNAE"

if [[ "$1" == "-h" ]] || [[ "$1" == "" ]]
then

echo -e
echo "Function: Quality Check: Quick check the main parts of output to ensure there are no errors"
echo -e
echo "Usage: ./glimpse_output_of_knisbacher_script_steps.sh <script_name> <mismatch type>"
echo -e
echo "Example runs:"
echo -e
echo "  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_1_prepare_data"
echo "  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_2_find_clusters"
echo "  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_3_create_tracks GA"
echo "  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_3_create_tracks CT"
echo "  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_4_Filter_clusters GA"
echo "  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_4_Filter_clusters CT"
echo "  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_5_create_filtered_tracks GA"
echo "  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_5_create_filtered_tracks CT"
echo "  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_6_Calculate_total_edit_sites GA"
echo "  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_6_Calculate_total_edit_sites CT"
echo -e

###########################################################################################################################################################################################################################################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#glimpse_output_of_knisbacher_script_step_1_prepare_data_81_birds

elif [[ "$1" == "knisbacher_script_step_1_prepare_data" ]]
then

cd $base_dir
while read -r i
  do
  o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
  cd $i/knisbacher
  j1=$(ls -lhrtR "$o"* | egrep "\.interval|\.bed|\.fa" | awk '{print$5}')
  j2=$(head -1 $o".fa")  #header
  j3=$(grep ">" $o".fa" -c)  #sequences
  if [[ -d Data ]]
    then
    j4=$(find Data/"$o"/"$2"/db/ -maxdepth 2 -mindepth 2 -path "*/files*" -type f | sort -V | head -1 | xargs head -1)  #family header
    j5=$(ls Data/"$o"/"$2"/db/*/ -d | wc -l)  #folders
    j6=$(find Data/"$o"/"$2"/db/ -maxdepth 1 -mindepth 1 -path "*/files*" -type d | xargs -n1 bash -c 'paste <(echo $0 | cut -f5 -d "/") <(ls $0 | wc -l)' | awk '{print$2"~("$1")"}' | paste -s -d ",")  #family sequences
    else j4="-"; j5="0"; j6="0"
  fi
  echo $i $o $j1 $j2 $j3 $j4 $j5 $j6
  unset o j1 j2 j3 j4 j5 j6
  cd $base_dir
done < all_bird_genomes_used | sed '1i Species Folder Interval_file Bed_file Fasta_file Fasta_header_sample Num_fasta Family_header_sample Num_subfamilies Num_fasta_per_subfamily_comma_seperated' | column -t > results/knisbacher_method/glimpse_output_of_knisbacher_script_step_1_prepare_data_"$2"_81_birds

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#glimpse_output_of_knisbacher_script_step_2_runClusterFinder_81_birds

elif [[ "$1" == "knisbacher_script_step_2_find_clusters" ]]
then

cd $base_dir
while read -r i
  do
  o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
  cd $i/knisbacher
  if [[ -d Data ]]
    then
    if [[ $(ls Data/"$o"/"$2"/results/GA/) != "" ]]
      then
      #j1=$(find Data/"$o"/"$2"/results/GA/ -name "cluster*.tab" -type f | wc -l)  #clusters made
      j1=$(find Data/"$o"/"$2"/results/ -maxdepth 1 -mindepth 1 | egrep -v "blasts|Tracks" | sort -V | xargs -n1 bash -c 'paste <(ls $0 | wc -l) <(wc -l $0/clusters_*.tab | awk "END{print\$1}")' | paste -s -d "\t" | sed 's/[\t]\+/ /g')  #number of cluster files & alignments/ mismatch type
      j2=$(find Data/"$o"/"$2"/results/GA/ -maxdepth 1 -mindepth 1 -name "cluster*.tab" -type f | xargs -n1 bash -c 'paste <(wc -l $0 | awk "{print\$1}") <(echo $0 | awk -F "/" "{print\$NF}" | cut -f4 -d "_")' | awk '{print$1"~("$2")"}' | paste -s -d ",")  #Lines in cluster files
      else j1="- - - - - - - - - - - -"; j2="- - - - - - - - - - - - -"
    fi
    else j1="- - - - - - - - - - - -"; j2="- - - - - - - - - - - - -"
  fi
  echo $i $j1 $j2
  unset o j1 j2
  cd $base_dir
done < all_bird_genomes_used \
| sed '1i Species AC_clus AC_elem AG_clus AG_elem AT_clus AT_elem CA_clus CA_elem CG_clus CG_elem CT_clus CT_elem GA_clus GA_elem GC_clus GC_elem GT_clus GT_elem TA_clus TA_elem TC_clus TC_elem TG_clus TG_elem Size_of_clusters' \
| column -t > results/knisbacher_method/glimpse_output_of_knisbacher_script_step_2_runClusterFinder_"$2"_81_birds

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#glimpse_output_of_knisbacher_script_step_3_create_tracks_81_birds

elif [[ "$1" == "knisbacher_script_step_3_create_tracks" ]]
then

#Specify mismatch type: GA / CT
if [[ "$3" == "" ]]
then
echo -e
echo " Specify mismatch type Eg: GA, CT"
echo -e
else

cd $base_dir
while read -r i
  do
  o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
  cd $i/knisbacher
  if [[ -d Data ]]
    then
    if [[ -d Data/"$o"/"$2"/results/Tracks ]]
      then
      j1=$(find Data/"$o"/"$2"/results/Tracks/ -maxdepth 4 -mindepth 3 -path "*/$3/raw/*" -type f | wc -l)  #total output files
      j2=$(find Data/"$o"/"$2"/results/Tracks/ -maxdepth 4 -mindepth 3 -path "*/$3/raw/*" -name "clusters*.tab" -type f | wc -l)  #clusters analyzed
      j3=$(find Data/"$o"/"$2"/results/Tracks/ -maxdepth 4 -mindepth 3 -path "*/$3/raw/*" -name "nucListFreqInConsPerPair_*.txt" -type f | wc -l)  #nucListFreqInConsPerPair_*.txt files
      j4=$(find Data/"$o"/"$2"/results/Tracks/ -maxdepth 5 -mindepth 4 -path "*/$3/raw/*" -name "blastToClass*.tab" -type f | wc -l)  #blast output files
      #j5=$(find Data/"$o"/"$2"/results/Tracks/ -maxdepth 5 -mindepth 4 -path "*/$3/raw/*" -name "blastToClass*.tab" -type f | xargs -n1 sh -c 'wc -l $0' | awk '{a+=$1;} END{print a;}')  #blast hits
      else j1="0"; j2="0"; j3="0"; j4="0"
    fi
    else j1="0"; j2="0"; j3="0"; j4="0"
  fi
  echo $i $j1 $j2 $j3 $j4
  unset o j1 j2 j3 j4
  cd $base_dir
done < all_bird_genomes_used | sed "1i Species Total_num_output_files Num_$3_Clusters_analyzed Num_Consensus_maps_made Num_Consensus_blast_outputs" | column -t > results/knisbacher_method/glimpse_output_of_knisbacher_script_step_3_create_tracks_"$2"_"$3"_81_birds

fi

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#glimpse_output_of_knisbacher_script_step_4_Filter_clusters_81_birds

elif [[ "$1" == "knisbacher_script_step_4_Filter_clusters" ]]
then

#Specify mismatch type: GA / CT
if [[ "$3" == "" ]]
then
echo -e
echo " Specify mismatch type Eg: GA, CT"
echo -e
else

cd $base_dir
while read -r i
  do
  o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
  cd $i/knisbacher
  if [[ -d Data ]]
    then
    if [[ -d Data/"$o"/"$2"/results/Tracks ]]
      then
      for track in $(find Data/"$o"/"$2"/results/Tracks/ -maxdepth 1 -mindepth 1 -type d)
        do
	subfam=$(echo $track | awk -F "/" '{print$6}' | awk -F "_" '{print$4}')
	if [[ $(find $track -maxdepth 1 -mindepth 1 -name "${3}" -type d) != "" ]]
	then
        size=$(find $track  -maxdepth 4 -mindepth 3 -path "*/$3/raw/*" -name "clusters_*.tab" | sort -V | xargs -n1 bash -c 'paste <(wc -l $0 | awk "{print\$1}")' | paste -s -d " ")
        echo "-" $subfam $size
        unset subfam size
	else echo "-" $subfam "0 0 0 0"
	unset subfam
	fi
      done | awk -v a="$i" 'NR==1{$1=a}1'
      else subfam="-"; size="0 0 0 0"
      echo $i $subfam $size
      unset subfam size
    fi
    else subfam="-"; size="0 0 0 0"
    echo $i $subfam $size
    unset subfam size
  fi
  unset o track
  cd $base_dir
done < all_bird_genomes_used | sed "1i Species Subfamily Unfiltered_$3_clusters Pairwise_filtered_$3_clusters Divergence_filtered_$3_clusters Map_to_G_than_A_filtered_$3_clusters" | column -t > results/knisbacher_method/glimpse_output_of_knisbacher_script_step_4_Filter_clusters_"$2"_"$3"_81_birds

fi

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------$
#glimpse_output_of_knisbacher_script_step_5_create_filtered_tracks

elif [[ "$1" == "knisbacher_script_step_5_create_filtered_tracks" ]]
then

#Specify mismatch type: GA / CT
if [[ "$3" == "" ]]
then
echo -e
echo " Specify mismatch type Eg: GA, CT"
echo -e
else
a=$(echo $3 | cut -c1)
b=$(echo $3 | cut -c2)

cd $base_dir
while read -r i
  do
  o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
  cd $i/knisbacher
  if [[ -d Data ]]
    then
    if [[ -d Data/"$o"/"$2"/results/Tracks ]]
      then
      if [[ $(find Data/"$o"/"$2"/results/Tracks/ -maxdepth 1 -mindepth 1 -type d) != "" ]]
        then
        for track in $(find Data/"$o"/"$2"/results/Tracks/ -maxdepth 1 -mindepth 1 -type d)
          do
	subfam=$(echo $track | awk -F "/" '{print$6}' | awk -F "_" '{print$4}')
          if [[ $(find $track -maxdepth 1 -mindepth 1 -name "${3}" -type d) != "" ]]
            then
            #motifs=$(find $track -maxdepth 3 -mindepth 3 -path "*/$3/pairwise_filter/*" -name "motifs_*_1e-3_*" | sort -V | xargs -n1 sh -c 'sed -n "2p" $0' | sed 's/: /(/g' | sed 's/$/)/g' | paste -s -d " ")
            if [[ $(find $track -maxdepth 3 -mindepth 3 -path "*/$3/pairwise_filter/*" -name "motifs_${a}_1e-3_*") == "" ]]; then motif1="-"; else
              motif1=$(find $track -maxdepth 3 -mindepth 3 -path "*/$3/pairwise_filter/*" -name "motifs_${a}_1e-3_*" | xargs -n1 sh -c 'grep "^Motifs statistics," -A1 $0 | egrep -v "Motifs statistics,|More statistics:"' | sed 's/: /(/g' | sed 's/$/)/g')
              if [[ $motif1 == "" ]]; then motif1="-"; else :;fi
            fi
            if [[ $(find $track -maxdepth 3 -mindepth 3 -path "*/$3/pairwise_filter/*" -name "motifs_${b}_1e-3_*") == "" ]]; then motif2="-"; else
              motif2=$(find $track -maxdepth 3 -mindepth 3 -path "*/$3/pairwise_filter/*" -name "motifs_${b}_1e-3_*" | xargs -n1 sh -c 'grep "^Motifs statistics," -A1 $0 | egrep -v "Motifs statistics,|More statistics:"' | sed 's/: /(/g' | sed 's/$/)/g')
              if [[ $motif2 == "" ]]; then motif2="-"; else :;fi
            fi
            size1=$(find $track -maxdepth 4 -mindepth 3 -path "*/$3/raw/*" -name "clusters_*.tab" | sort -V | xargs -n1 bash -c 'paste <(wc -l $0 | awk "{print\$1}")' | paste -s -d " ")
            if [[ $(find $track -maxdepth 3 -mindepth 3 -path "*/$3/pairwise_filter/*" -name "bestPairsClusters_*.tab") == "" ]]; then size2="0"; else
              size2=$(find $track -maxdepth 3 -mindepth 3 -path "*/$3/pairwise_filter/*" -name "bestPairsClusters_*.tab" | xargs -n1 bash -c 'paste <(wc -l $0 | awk "{print\$1}")' | paste -s -d " ")
            fi
            echo "-" $subfam $size1 $size2 $motif1 $motif2
            unset subfam size1 size2 motif1 motif2
            else
            echo "-" $subfam "0 0 0 0 0 - -"
	    unset subfam
          fi
        done | awk -v a="$i" 'NR==1{$1=a}1'
        else subfam="-"; size="0 0 0 0 0 - -"
        echo $i $subfam $size
        unset subfam size
      fi
      else subfam="-"; size="0 0 0 0 0 - -"
      echo $i $subfam $size
      unset subfam size
    fi
    else subfam="-"; size="0 0 0 0 0 - -"
    echo $i $subfam $size
    unset subfam size
  fi
  unset o track
  cd $base_dir
done < all_bird_genomes_used | sed "1i Species Subfamily Unfiltered_$3_clusters Pairwise_filtered_$3_clusters Divergence_filtered_$3_clusters Map_to_${a}_than_${b}_filtered_$3_clusters Best_sources_filtered_$3_clusters ${a}_Motifs(Num) ${b}_Motifs(Num)" | column -t > results/knisbacher_method/glimpse_output_of_knisbacher_script_step_5_create_filtered_tracks_"$2"_"$3"_81_birds

unset a b
fi

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------$
#glimpse_output_of_knisbacher_script_step_6_Calculate_total_edit_sites

elif [[ "$1" == "knisbacher_script_step_6_Calculate_total_edit_sites" ]]
then

#Specify mismatch type: GA / CT
if [[ "$3" == "" ]]
then
echo -e
echo " Specify mismatch type Eg: GA, CT"
echo -e
else

cd $base_dir
while read -r i
        do
        o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
        cd $i/knisbacher
        if [[ -d Data ]]
                then
                if [[ -d Data/"$o"/"$2"/results/Total_edit_sites ]] && ls Data/"$o"/"$2"/results/Total_edit_sites/"$3"_TotalEditSites_*_fam_*.txt &> /dev/null
                        then
                        j1=$(awk '{a+=$NF} END{print a;}' Data/"$o"/"$2"/results/Total_edit_sites/"$3"_TotalEditSites_*_fam_*.txt)
                        j2=$(awk '{print$4"("$5,$6,$7")"}' OFS="," Data/"$o"/"$2"/results/Total_edit_sites/"$3"_TotalEditSites_*_fam_*.txt  | paste -s -d " ")
                        else j1="0"; j2="-"
                fi
                else j1="0"; j2="-"
        fi
        echo $i $j1 $j2
        unset j1 j2
        cd $base_dir
done < all_bird_genomes_used | sed "1i Species Total_${3}_edit_sites Family(${3}_count,2nd_highest_count,Corrected_${3}_count)" | column -t > results/knisbacher_method/glimpse_output_of_knisbacher_script_step_6_Calculate_total_"$2"_"$3"_edit_sites_81_birds

fi

##################################################################################################################################################################################################################################################################################$

else :
fi

##################################################################################################################################################################################################################################################################################$

#Draft
#find Data/"$o"/LTR/results/ -maxdepth 1 -mindepth 1 | egrep -v "blasts|Tracks" | sort -V | xargs -n1 bash -c 'paste <(echo $0 | awk -F "/" "{print\$NF \"_clusters\",\$NF\"_elements\"}")' | paste -s -d " "
#find Data/"$o"/LTR/results/ -maxdepth 1 -mindepth 1 | egrep -v "blasts|Tracks" | sort -V | xargs -n1 bash -c 'paste <(echo $0 | awk -F "/" "{print\$NF \"_clus\",\$NF\"_elem\"}")' | paste -s -d " "

