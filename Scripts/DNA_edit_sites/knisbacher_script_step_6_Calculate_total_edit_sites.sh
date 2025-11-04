#!/bin/bash

##################################################################################################################################################################################################################################################################################

#Function: Calculate background mutation rate (since insertion of repeat in genome) & total edit sites

#Concept:
        #The Screen for APOBEC associated edit sites in genomic region is done by identifying clusters of G-to-A mutations with no other mismatches allowed.
        #This approach reduces false positives but also underestimates the number of edited sites, due to random mutations that accumulate within the edited sequence, preventing detection of many clusters and individual editing sites.
        #This is evident in the reduced number of editing sites found after applying the Consensus filter and Best source selection.
        #To get a more precise estimate, analyze the alignment between edited element & it's `best sources` in the following way:
                #Count the total number of mismatches (mismatch of interest: GA) not just from the cluster include even individual mismatches
                #Count the 2nd highest occuring mismatch
                #Subtract the 2nd highest occuring mismatch from total count of mimsatch of interest
        #This estimates the background mutation rate caused by random mutagenesis since the edited element was inserted.
        #This can increase the total number of detected editing sites by ~ 50%.

#Input:
        #`Best sources` filtered clusters named `bestPairsClusters_*.tab`.

##################################################################################################################################################################################################################################################################################

#Usage: ./knisbacher_script_step_6_Calculate_total_edit_sites.sh mm

#Options:
        #1. family (Repeat type): available options are: LTR, DNA, LINE, SINE
        #1. mm (mismatch type) : available options are: GA, CT

#Example runs:
        #nohup ./knisbacher_script_step_6_Calculate_total_edit_sites.sh LTR GA &
        #nohup ./knisbacher_script_step_6_Calculate_total_edit_sites.sh DNA GA &
        #nohup ./knisbacher_script_step_6_Calculate_total_edit_sites.sh LTR CT &

##################################################################################################################################################################################################################################################################################

#Set path variables
base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"
#script_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/scripts/Perl_DNAE"

#Root directory
cd $base_dir

#Set perl version
#perlbrew use perl-5.16.3

#Project
echo -e "\nESTIMATE APOBEC INDUCED DNA EDITED SITES IN AVIAN LTR RETROTRANSPOSONS\n\n - Calculate Total Edit Sites for Filtered Best Pair $2 Clusters in $1:\n" > stdout_knisbacher_script_step_6_Calculate_total_edit_sites_for_filtered_best_pairs_"$1"_"$2"_clusters_81_birds

###########################################################################################################################################################################################################################################$

#Mismatch type
mm=$2

#Column to extract from cluster table
num=$(echo "GA CT GC GT CA TA AG TC CG TG AC AT" | awk -v m="$mm" '{for(i=1;i<=NF;i++) if($i == m) print i}')

###########################################################################################################################################################################################################################################$
#Calculate for each families of each birds

{
	time while read -r i
	do
	o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
	cd $i/knisbacher
	echo ">"$i "(" $o ")"
	echo -e "\nESTIMATE APOBEC INDUCED DNA EDITED SITES IN AVIAN LTR RETROTRANSPOSONS\n\n - Calculate Total Edit Sites for Filtered Best Pair $2 Clusters in $1:\n" > stdout_knisbacher_script_step_6_Calculate_total_edit_sites_for_filtered_best_pairs_"$1"_"$2"_clusters

	#Calculate count of mismatch of interest (would have the highest count), 2nd most common & their difference
	{
		if [[ -d Data ]]
			then
			if [[ $(find Data/"$o"/"$1"/results -maxdepth 5 -mindepth 5 -path "*/Tracks/*/${mm}/pairwise_filter*" -name "bestPairsClusters_*.tab") != "" ]]
				then

				#Create rseperate folder for total edit sites
				if [[ -d Data/"$o"/"$1"/results/Total_edit_sites ]]
					then :
					else
					mkdir -p Data/"$o"/"$1"/results/Total_edit_sites
				fi

				#Delete if the count file already exists. Be cautious of what is being deleted!
				if ls Data/"$o"/"$1"/results/Total_edit_sites/"$mm"_TotalEditSites_*.txt 1> /dev/null 2>&1
					then
					rm Data/"$o"/"$1"/results/Total_edit_sites/"$mm"_TotalEditSites_*.txt
					else :
				fi
				time for bestpairs in $(find Data/"$o"/"$1"/results -maxdepth 5 -mindepth 5 -path "*/Tracks/*/${mm}/pairwise_filter*" -name "bestPairsClusters_*.tab")
					do
					#Create suffixes / files to save counts for subfamily & family
					sf=$(echo $bestpairs | awk -F "/" '{print$NF}' | awk -v m="$mm" 'BEGIN{FS=OFS="_"} {print m,"TotalEditSites",$2,$3,$4,"subfam",$5".txt"}')
					f=$(echo $bestpairs | awk -F "/" '{print$NF}' | awk -v m="$mm" 'BEGIN{FS=OFS="_"} {print m,"TotalEditSites",$2,$3,"fam",$5".txt"}')
					#touch Data/"$o"/LTR/results/Total_edit_sites/"$f"
					echo "  - " $f

					while read -r bp
						do
						s=$(echo $bp | awk '{print$1,$2,$3,$4,$5}')
						c=$(echo $bp | awk '{print$NF}')
						i1=$(echo $c | cut -f $num -d "|")
						i2=$(echo $c | tr "|" "\n" | sed "$num d" | sort -nr | head -1)
						i3=$(calc $i1 - $i2)
						echo $s $i1 $i2 $i3
						unset s c i1 i2 i3
					done < $bestpairs > Data/"$o"/"$1"/results/Total_edit_sites/"$sf"
					awk '{a+=$6; b+=$7; c+=$8;} END{print$1,$2,$3,$4,a,b,c}' Data/"$o"/"$1"/results/Total_edit_sites/"$sf" >> Data/"$o"/"$1"/results/Total_edit_sites/"$f"
					unset f bp
				done
				else "  - No Best Pair Clusters!"
			fi
			else "  - No Data!"
		fi
	} &>> stdout_knisbacher_script_step_6_Calculate_total_edit_sites_for_filtered_best_pairs_"$1"_"$2"_clusters

	unset o bestpairs
	cd $base_dir
	done < all_bird_genomes_used
} &>> stdout_knisbacher_script_step_6_Calculate_total_edit_sites_for_filtered_best_pairs_"$1"_"$2"_clusters_81_birds

unset mm num

#END#
###########################################################################################################################################################################################################################################$
