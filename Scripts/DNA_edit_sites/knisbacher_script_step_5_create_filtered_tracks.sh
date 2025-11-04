#!/bin/bash

##################################################################################################################################################################################################################################################################################

#Function: To apply best sources filter to clusters & find motif.

#Concept: Apply more filters & sequence motif context for accurate identification of APOBEC-associated mismatches
	  #1. Best sources filter:
		#The screen for clustered G-to-A mutations reduces false positives but underestimates the total number of edited sites per retroelement.
		#To get a complete count of edited sites, all G-to-A mismatches in the alignment to the source element were counted.
		#The second most common mismatch was subtracted to account for random mutations after insertion.
		#The most likely source element was determined by finding elements with signs of DNA editing and selecting the one with the highest BLAST bitscore.
	  #2. Find motifs:
		#APOBECs have been shown to preferentially edit cytidines in specific sequence contexts.
		#Eg: motif most preferred by:
			# mouse APOBEC3 is GxA (G is edited; x is any nucleotide)
			# human APOBEC3G and APOBEC3F are GG and GA
		#Crossmatching known motifs with motifs found around edit sites suggests that serve as an additional layer of confidence.
		#Finding consistent motifs validate both the editing sites identified and the motifs themselves.
		#We could infer APOBEC preferences in novel genomes as well and supply a broader understanding of APOBEC editing preferences.

	#The main script used is written to parses filtered cluster files and creates several UCSC-tracks and analysis-output files
	#More info about the steps is explained in the script createTrackFiles2.pl

#Input:
	#Cluster files in the directory: Tracks/tracks.*/GA/raw/filteredClusts
	  #Cluster file is not exclusively given but automatically chosen using the paramaters values such as: organism, class, family, pval, th & filter
	#Script used: createTrackFiles2.pl

##################################################################################################################################################################################################################################################################################

#Usage: ./knisbacher_script_step_5_create_filtered_tracks.sh mm

#Options:
        #1. family (Repeat type): available options are: LTR, DNA, LINE, SINE
        #1. mm (mismatch type) : available options are: GA, CT

#Example runs:
        #nohup ./knisbacher_script_step_5_create_filtered_tracks.sh LTR GA &
        #nohup ./knisbacher_script_step_5_create_filtered_tracks.sh DNA GA &
        #nohup ./knisbacher_script_step_5_create_filtered_tracks.sh LTR CT &

##################################################################################################################################################################################################################################################################################

#Set path variables
base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"
script_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/scripts/Perl_DNAE"

#Root directory
cd $base_dir

#Set perl version
perlbrew use perl-5.16.3 &> /dev/null

#Project
echo -e "\nESTIMATE APOBEC INDUCED DNA EDITED SITES IN AVIAN LTR RETROTRANSPOSONS\n\n - Create Tracks for Filtered $2 Clusters in $1:\n" > stdout_knisbacher_script_step_5_create_tracks_for_filtered_"$1"_"$2"_clusters_81_birds

###########################################################################################################################################################################################################################################$

#Run createTrackFiles2.pl for each families of each birds
{
  time while read -r i
  do
  o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
  cd $i/knisbacher
  echo ">"$i "(" $o ")"
  echo -e "\nESTIMATE APOBEC INDUCED DNA EDITED SITES IN AVIAN LTR RETROTRANSPOSONS\n\n - Create Tracks for Filtered $2 Clusters in $1:\n" > stdout_knisbacher_script_step_5_create_tracks_for_filtered_"$1"_"$2"_clusters
  echo -e "  - Run createTrackFiles2.pl..." >> stdout_knisbacher_script_step_5_create_tracks_for_filtered_"$1"_"$2"_clusters

  #Parses clusters & apply best sources filter & find motifs
  {
    if [[ -d Data ]]
      then
      if [[ -d Data/"$o"/"$1"/results/Tracks ]]
        then
        if [[ $(find Data/"$o"/"$1"/results/Tracks/ -maxdepth 5 -mindepth 4 -path "*/$2/raw/filteredClusts/*" -name "clusters_*_cleanAlign_targetMoreDiv_moreMapPre.tab") != "" ]]
          then
          echo "   - Apply best sources filter:"
          time for filter in $(find Data/"$o"/"$1"/results/Tracks/ -maxdepth 5 -mindepth 4 -path "*/$2/raw/filteredClusts/*" -name "clusters_*_cleanAlign_targetMoreDiv_moreMapPre.tab")
            do
            family=$(echo $filter | awk -F "/" '{print$NF}' | cut -f4 -d "_")
            echo "    - "$family $filter
            time perl $script_dir/analysis_scripts/createTrackFiles2.pl \
              --dataDir Data \
              --organism $o \
              --class $1 \
              --family $family \
              --filter "clusters|1|$filter|pairwise_filter" \
              --pval 0 \
              --th 5 \
              --cores 32 \
              --mm $2 \
              --GEPIC_consRoot /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/RepBase/Libraries \
              --gepic_consall /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/RepBase/Libraries/RMRBSeqs.fa &>> stdout_knisbacher_script_step_5_create_tracks_for_filtered_"$1"_"$2"_clusters
            unset family
          done
          else
          echo -e "   - No filtered clusters to apply best sources filter! \n"
        fi
        else
        echo -e "   - No filtered clusters to apply best sources filter! \n"
      fi
      else
      echo -e "   - No filtered clusters to apply best sources filter! \n"
    fi
  } &>> stdout_knisbacher_script_step_5_create_tracks_for_filtered_"$1"_"$2"_clusters

  unset o filter
  cd $base_dir
  done < all_bird_genomes_used
} &>> stdout_knisbacher_script_step_5_create_tracks_for_filtered_"$1"_"$2"_clusters_81_birds


#END#
###########################################################################################################################################################################################################################################$

#Parameters to tune with createTrackFiles2.pl
#--skip_pos_in_cons 1  #Don't need to do consensus mapping, consensus filter already applied in knisbacher_script_step_3_create_tracks.sh
