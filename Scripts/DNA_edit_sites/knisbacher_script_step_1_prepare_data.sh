#!/bin/bash

##################################################################################################################################################################################################################################################################################

#Project
echo -e "\nESTIMATE APOBEC INDUCED DNA EDITED SITES IN AVIAN LTR RETROTRANSPOSONS\n\n - Prepare Data:\n"

#Function: Prepare Data to run downstream analysis of detecting & counting the number APOBEC associated GA edit sites in avian LTR retrotransposons

#Concept: #Hypermutation of retroelement DNA by APOBECs creates an abundance of G-to-A mutations that tend to be in clusters and not uniformly distributed throughout the targeted element.
	  #Here a sequence alignment based approach was used to identify such clusters in LTR retrotransposons of genomes.
	  #Hypotheses: Every edited element in the genome has a closely related non-edited ‘parental’ element, such that their alignment will contain identical regions containing a series of unidirectional G-to-A mutations.
	  #Thus, for every given genome, pairwise alignments of retroelements of the same repeat subfamily (based on annotation in RepeatMasker output) were generated (Further screened for clusters of G>A mutations in downstream steps).

#Runs mainly 3 steps:
	#1. Run rmskOutToInterval.pl: Converts the repeatmasker (rmsk) file to the Interval file formats
	#2. Run Genome2Fasta.pl: Extract fasta from genome
	#3. Run sortGenome.pl: Sorts genome-wide repeats to files by classes, families and names. Also creates a directory, subdirectory and file for each Class, Family and Name respectively

#Input:
	#1. A file containing list of path to directories of each bird genomes. This directory must contain genome and repeatmasker file.

	   #Eg: file: all_bird_genomes_used

		#>head -3 all_bird_genomes_used
		  # Apteryx_rowi
		  # Aquila_chrysaetos_chrysaetos
		  # Atlantisia_rogersi

		#>ls Apteryx_rowi/ucsc/
		  # rw-rw-r-- 1 ceglab25 ceglab25  66M Dec 30  2022 GCF_003343035.1.repeatMasker.out
		  # -rw-rw-r-- 1 ceglab25 ceglab25 1.2G Dec 30  2022 GCF_003343035.1.fa

##################################################################################################################################################################################################################################################################################

#Usage: ./knisbacher_script_step_1_prepare_data.sh

#Options: This script is written in a way to just execute without giving any paramaters!
          #To fine tune the steps, change parameters used in the script `sortGenome.pl` or replace `Genome2Fasta.pl` with some other ways to extrat repeat sequences.

#Example run:
        #nohup ./knisbacher_script_step_1_prepare_data.sh &

##################################################################################################################################################################################################################################################################################

#Set input file
if [[ -z $1 ]]; then
    input_file="all_bird_genomes_used"   # Default input file if no argument is provided
else
    input_file="$1"                      # Set input file to the provided argument
fi

#Set path variables
base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"
script_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/scripts/Perl_DNAE"

#Root directory
cd $base_dir

#Set perl version
perlbrew use perl-5.16.3 &> /dev/null

##################################################################################################################################################################################################################################################################################

#Prepare data for 81 bird genome
{
	time while read -r i
	do

	# Organism suffix
	  o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
	  echo -e ">Org: "$i "("$o")"

	# Create folder
	  cd $i
	  mkdir knisbacher

	# Get genome file & repeatmasker file
	  if [[ -d ucsc ]]
	  then
	  cd ucsc
	  genome=$(readlink -f GC*.fa | egrep -v "rmsk|repeatModeler")
	  rmsk=$(ls *.repeatMasker.out)
	  else
	  genome=$(readlink -f GC*.fna)
	  cd repeat_masker_using_"$i"_library
	  rmsk=$(ls *.fna.out)
	  fi

	# Check if accession versions are same between genome & repeatmasker table
	  ga=$(echo $genome | awk -F "/" '{print$NF}' | sed 's/\(^.*\.[0-9]\).*/\1/g')
	  ra=$(echo $rmsk | awk -F "/" '{print$NF}' | sed 's/\(^.*\.[0-9]\).*/\1/g')
	  if [[ "$ga" == "$ra" ]]
	  then
	  echo -e "  - Inputs:\n   |- Genome: $genome\n   |- RepeatMasker file: $rmsk\n"
	  cp $rmsk ../knisbacher/"$o"_RepeatMasker.out
	  cd ../knisbacher

	echo -e "\nESTIMATE APOBEC INDUCED DNA EDITED SITES IN AVIAN LTR RETROTRANSPOSONS\n\n - Prepare Data:\n" &> stdout_knisbacher_script_step_1_prepare_data

	# Converts the repeatmasker (rmsk) file to the Interval file formats
	  echo -e "\n   - Run rmskOutToInterval.pl..." &>> stdout_knisbacher_script_step_1_prepare_data
	  { time perl $script_dir/Tools/rmskOutToInterval.pl "$o"_RepeatMasker.out $o 1 ; } &>> stdout_knisbacher_script_step_1_prepare_data

	# Extract fasta from genome
	  echo -e "\n   - Run Genome2Fasta.pl..." &>> stdout_knisbacher_script_step_1_prepare_data
	  { time perl $script_dir/Tools/Genome2Fasta.pl $genome "$o".interval $o".fa" 1 ; } &>> stdout_knisbacher_script_step_1_prepare_data

	# Rename fasta headers: Need to add assembly name at head of each defline & replaces spacers with underscores (0m18.602s)
	  sed -e "/^>/ s/>[^_]\+_/>$o\_/g" -e '/^>/ s/ /_/g' $o".fa" -i

	# Sorts genome-wide repeats to file by classes, families and names. Also creates a directory, subdirectory and file for each Class, Family and Name respectively
	  echo -e "\n   - Run sortGenome.pl..." &>> stdout_knisbacher_script_step_1_prepare_data
	  { time perl $script_dir/sortGenome.pl --interval "$o".interval --fasta $o".fa" --dataDir Data --org $o --lc --classes LTR ; } &>> stdout_knisbacher_script_step_1_prepare_data
          { time perl $script_dir/sortGenome.pl --interval "$o".interval --fasta $o".fa" --dataDir Data --org $o --lc --classes DNA ; } &>> stdout_knisbacher_script_step_1_prepare_data
          { time perl $script_dir/sortGenome.pl --interval "$o".interval --fasta $o".fa" --dataDir Data --org $o --lc --classes LINE ; } &>> stdout_knisbacher_script_step_1_prepare_data
	  echo -e "\nDONE.....\n" &>> stdout_knisbacher_script_step_1_prepare_data

	  echo -e "DONE.....\n"
	  else
	  echo -e " Bad input: Accession versions are different between genome & RepeatMasker file :\n  - Genome Accession: $ga \n  - RepeatMasker accession: $ra\n"
	  fi

	#Unset temporary variables
	unset o genome rmsk ga ra
	cd $base_dir
	done < $input_file
}

#END#
##################################################################################################################################################################################################################################################################################
