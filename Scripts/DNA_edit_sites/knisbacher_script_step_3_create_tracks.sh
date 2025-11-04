#!/bin/bash

##################################################################################################################################################################################################################################################################################

#Function: To create pairwise alignments between edited element & the closest repbase consensus sequence & create a frequency data for source/target sequence for each mismatch pair.

#Concept: The consensus filter uses the consensus sequence of each element’s subfamily to identify which of the two aligned sequences is most probably the ancestral one.
	  #This information is crucial to ascertain that the direction of the mutations is G-to-A and not the opposite.
	  #The maion script used is written to parses cluster files and creates several UCSC-tracks and analysis-output files
	  #More info about the steps is explained in the script createTrackFiles2.pl

#Input:
	#Cluster files in sepcific directory heirarchy (made when running sortGenome.pl & runClusterFinder.pl)
	  #Cluster file is not exclusively given but automatically chosen using the paramaters values such as: organism, class, family, pval, th
	#Repbase consensus sequences: single multifasta & individual fasta per sequences.
	#Script used: createTrackFiles2.pl

##################################################################################################################################################################################################################################################################################

#Usage: ./knisbacher_script_step_3_create_tracks.sh mm

#Options:
	#1. family (Repeat type): available options are: LTR, DNA, LINE, SINE
        #2. mm (mismatch type) : available options are: GA, CT

#Example runs:
        #nohup ./knisbacher_script_step_3_create_tracks.sh LTR GA &
	#nohup ./knisbacher_script_step_3_create_tracks.sh DNA GA &
        #nohup ./knisbacher_script_step_3_create_tracks.sh LTR CT &

##################################################################################################################################################################################################################################################################################

#Set path variables
base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"
script_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/scripts/Perl_DNAE"

#Root directory
cd $base_dir

#Set perl version
perlbrew use perl-5.16.3 &> /dev/null

#Project
echo -e "\nESTIMATE APOBEC INDUCED DNA EDITED SITES IN AVIAN LTR RETROTRANSPOSONS\n\n - Create Tracks for $2 clusters in $1:\n" > stdout_knisbacher_script_step_3_create_tracks_for_"$1"_"$2"_clusters_81_birds

###########################################################################################################################################################################################################################################$

#Run createTrackFiles2.pl for each families of each birds
{
  time while read -r i
  do
  o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
  cd $i/knisbacher
  echo ">"$i "(" $o ")"
  echo -e "\nESTIMATE APOBEC INDUCED DNA EDITED SITES IN AVIAN LTR RETROTRANSPOSONS\n\n - Create Tracks for $2 clusters in $1:\n" > stdout_knisbacher_script_step_3_create_tracks_for_"$1"_"$2"_clusters
  echo -e "  - Run createTrackFiles2.pl..." >> stdout_knisbacher_script_step_3_create_tracks_for_"$1"_"$2"_clusters

  #Create pairwise alignments b/w edited element & consensus sequence (No need to find motifs)
  {
    time for family in $(ls Data/"$o"/"$1"/results/blasts/)
    do
    echo "   - "$family
    time perl $script_dir/analysis_scripts/createTrackFiles2.pl \
	--dataDir Data \
	--organism $o \
	--class $1 \
	--family $family \
	--pval 0 \
	--th 5 \
	--cores 32 \
	--mm $2 \
	--no_find_motifs 1 \
	--GEPIC_consRoot /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/RepBase/Libraries \
	--gepic_consall /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/RepBase/Libraries/RMRBSeqs.fa &>> stdout_knisbacher_script_step_3_create_tracks_for_"$1"_"$2"_clusters
    done
    } &>> stdout_knisbacher_script_step_3_create_tracks_for_"$1"_"$2"_clusters

  unset o family
  cd $base_dir
  done < all_bird_genomes_used
} &>> stdout_knisbacher_script_step_3_create_tracks_for_"$1"_"$2"_clusters_81_birds

#END#
###########################################################################################################################################################################################################################################$
