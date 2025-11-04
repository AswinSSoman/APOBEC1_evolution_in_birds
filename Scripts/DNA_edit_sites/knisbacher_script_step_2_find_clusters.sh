#!/bin/bash

##################################################################################################################################################################################################################################################################################

#Function: Create pairwise alignments using BLAST only between same subfamilies & finds clusters in the BLAST output. Then parses the cluster, for each Name (subfamily) of sequences of the organism.
           #More info about the steps is explained in the script runClusterFinder.pl

#Concept:
	#DNA editing of retroelements by APOBECs introduces G-to-A mutations specifically in the retroelements’ sense strand.
	#The mutations are nonuniformly distributed throughout the sequence and tend to be found in clusters.
	#In this approach, assume that edited elements have unedited “ancestral” elements in the genome.
	#Thus, an edited element should be very similar to one or more ancestral or “source” elements, except for the editing sites, where the source will contain guanosines and the edited element adenosines.
	#To detect clustered mutations, align pairs of LTRs of the same subfamily in each genome using BLAST.
	#In order to create alignments only between closely related retroelements, set stringent criterias in blastn search. E
		#Eg: E-value 1E-50, 250 top alugnments, plus strand, disable sequence filtering

#Input:
	#Repeat sequence fasta file in sepcific directory heirarchy (Data > org > class > family > (2 folders: db & results) -> made when running sortGenome.pl)
	#P-value higher & lower threshold: Control loose to stringent
	#Cluster length higher & lower threshold: Minimum number of consecutive mismatches in pairwise alignments to define a cluster eg: 5
	#blastn parameters
	#Mismatch type: eg: GA, CT or all 12 types
	#Script used: runClusterFinder.pl

##################################################################################################################################################################################################################################################################################

#Usage: ./knisbacher_script_step_2_find_clusters.sh Repeat_family

#Options: To fine tune the steps, change parameters used in the script `runClusterFinder.pl`.

#Options:
        #1. family (Repeat type): available options are: LTR, DNA, LINE, SINE

#Example run:
	#nohup ./knisbacher_script_step_2_find_clusters.sh LTR &
	#nohup ./knisbacher_script_step_2_find_clusters.sh DNA &

##################################################################################################################################################################################################################################################################################

#Set path variables
base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"
script_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/scripts/Perl_DNAE"

#Root directory
cd $base_dir

#Set perl version
perlbrew use perl-5.16.3 &> /dev/null

#Project
echo -e "\nESTIMATE APOBEC INDUCED DNA EDITED SITES IN AVIAN LTR RETROTRANSPOSONS\n\n - Find clusters in $1:\n" > stdout_knisbacher_script_step_2_runClusterFinder_"$1"_81_birds

###########################################################################################################################################################################################################################################

#Run runClusterFinder.pl for each bird
{
  time while read -r i
  do
  o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
  cd $i/knisbacher
  echo ">"$i "(" $o ")"
  echo -e "\nESTIMATE APOBEC INDUCED DNA EDITED SITES IN AVIAN LTR RETROTRANSPOSONS\n\n - Find clusters in $1:\n" > stdout_knisbacher_script_step_2_runClusterFinder_"$1"
  echo -e "   - Run runClusterFinder.pl..." >> stdout_knisbacher_script_step_2_runClusterFinder_"$1"

  {
    time perl $script_dir/runClusterFinder.pl \
	--dataDir Data \
	--organism $o \
	--pval_h 0 \
	--pval_l 0 \
	--th_l 5 \
	--th_h 5 \
	--allmms 0 \
	--directional 1 \
	--blastEvalue "1e-50" \
	--blastargs "-num_alignments 250" \
	--blastn_threads 32 \
	--num_parallel_analyze_blast 32 \
	--classes $1
  } &>> stdout_knisbacher_script_step_2_runClusterFinder_"$1"

  unset o
  cd $base_dir
  done < all_bird_genomes_used
} &>> stdout_knisbacher_script_step_2_runClusterFinder_"$1"_81_birds


#END#
###########################################################################################################################################################################################################################################$
