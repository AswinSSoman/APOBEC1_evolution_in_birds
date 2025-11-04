###########################################################################################################################################################################################################################################################################################################
#                                                                                                               ESTIMATE APOBEC INDUCED DNA EDITED SITES IN AVIAN LTR RETROTRANSPOSONS
###########################################################################################################################################################################################################################################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. DOWNLOAD DATA

#__________________________________________________________________________________________________________________________________________________________
#1.1. Download genome & repeatmasker tables

  #Genomic & repeatmasker data is already downloaded for 81 birds & is located in:
  base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"
  cd $base_dir
  #Better abbotation of bird repeats have been done by several studies including works done by Alexander Suh [https://scholar.google.de/citations?user=-J4CG5QAAAAJ&hl=de]
  
  #Birds of interest that is lacking in this: Parus major
  
  #Stats & QC the genome & repeatmasker table
  ./QC_genomic_data.sh

#__________________________________________________________________________________________________________________________________________________________
#1.2. Download scripts

  git clone https://github.com/binknis/DNA-editing.git scripts  #DNA editing script
  #git clone https://github.com/adadiehl/repeatMaskerUtils.git  #repeatMaskerUtils
  script_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/scripts/Perl_DNAE"
  
  #Modifications to be done inside scripts
  # set perl version
  perlbrew use perl-5.16.3

#__________________________________________________________________________________________________________________________________________________________
#1.2.1. Modifications to be done inside scripts

  #------------------------------------------------------------------------------------
    #Change hardcoded perl version
    find scripts/ -type f -exec sed -i 's/perl516/perl/g' {} +
    
  #------------------------------------------------------------------------------------
    #Irrepspective of what parameters are used modify the following in analysisSubs.pm       
    #uncomment line in 961 & 962
    #comment line 963 & 964
    #In line 962: replace -num_threads 8 to -num_threads 32
  
  #------------------------------------------------------------------------------------
    #Modification in createTrackFiles2.pl
    #Comment line 207 in 
    #symlink($cluster_file, $tabular_file);
    #insert the following line after line 206
     copy($cluster_file, $tabular_file) or die "Copy failed: $!";
    #insert the folowing after line 39
    use File::Copy; 
  
  #------------------------------------------------------------------------------------
    #Modification in DNA-editing/Perl_DNAE/analysis_scripts/getFilteredClusterFile.p
    #Replace in line 87:  "nucListFreqPerPair_" with "nucListFreqInConsPerPair_"

#__________________________________________________________________________________________________________________________________________________________
#1.3. Download repbase consensus sequence from github (1m22.319s)

  git clone https://github.com/yjx1217/RMRB.git RepBase
  cd RepBase
  tar -xvzf RepBaseRepeatMaskerEdition-20181026.tar.gz
  cd Libraries
  perl /media/aswin/programs/RepeatMasker/util/buildRMLibFromEMBL.pl RMRBSeqs.embl > RMRBSeqs.lib
  /media/aswin/programs/repeatMaskerUtils/extractAncestral.py RMRBSeqs.embl > RMRBSeqs.fa
  
  #split multifasta consensus to individual fasta files, this is required for a downstream analysis
  mkdir perSeq
  time while read line
  do
  if [[ ${line:0:1} == '>' ]]
  then
  outfile=${line#>}.fa
  echo $line > perSeq/$outfile
  else
  echo $line >> perSeq/$outfile
  fi
  done < RMRBSeqs.fa

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#2. PREPARE DATA FOR ANALYSIS

#__________________________________________________________________________________________________________________________________________________________
#Function to ensure max number of jobs are running in the background at a time
  
  check_jobs() {
    while [ "$counter" -ge "$max_jobs" ]; do
      # Continue to wait for any job to finish and decrement the counter for each finished job
      while wait -n; do
        counter=$((counter - 1))
        # Break out of the inner loop if the counter is below the max_jobs threshold
        [ "$counter" -lt "$max_jobs" ] && break
      done
    done
  }

#__________________________________________________________________________________________________________________________________________________________
#2.1. Sort & extract repeats based on class, family & subfamilies

  #Set parameters to run parallel runs
  start_time=$(date +%s) && max_jobs=32 && counter=0
  
  #Run script knisbacher_script_step_1_prepare_data.sh for 32 birds in parallel background jobs
  for species in $(cat all_bird_genomes_used)
    do
    echo $species > "split_"$species
    nohup ./knisbacher_script_step_1_prepare_data.sh "split_"$species > "stdout_knisbacher_script_step_1_prepare_data_81_birds_split_"$species 2>&1 &
    counter=$((counter + 1))
    echo ">Species: "$species "| Background job: " $counter
    # Check if the number of running jobs has reached max_jobs, if so, wait for one or more to finish
    check_jobs
  done
  wait
  end_time=$(date +%s) && elapsed_time=$((end_time - start_time))  #total time taken

  echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e > stdout_knisbacher_script_step_1_prepare_data_81_birds
  cat stdout_knisbacher_script_step_1_prepare_data_81_birds_split_* >> stdout_knisbacher_script_step_1_prepare_data_81_birds
  unset start_time end_time elapsed_time max_jobs counter

#QC: Quick check the main parts of output to ensure there are no errors
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_1_prepare_data LTR
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_1_prepare_data DNA
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_1_prepare_data LINE
  
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. Main analysis

#__________________________________________________________________________________________________________________________________________________________
# 3.1. Runs BLAST, finds clusters in the BLAST output and parses the cluster, for each Name LTR subfamily of sequences of the organism.

  nohup ./knisbacher_script_step_2_find_clusters.sh LTR &
  nohup ./knisbacher_script_step_2_find_clusters.sh DNA &  #Negative control
  nohup ./knisbacher_script_step_2_find_clusters.sh LINE &  #control

#QC: Quick check the main parts of output to ensure there are no errors
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_2_find_clusters LTR
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_2_find_clusters DNA
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_2_find_clusters LINE

  
  #grep "clusters.*fields" -B3 scripts/Perl_DNAE/analysis_scripts/getFilteredClusterFile.pl | grep "my \$[A-Za-z_]\+" -o | awk '{print$2}' | tr -d "$" | paste -s -d " "  #Print fields of cluster file

#__________________________________________________________________________________________________________________________________________________________
#3.2. Parses cluster files and creates several UCSC-tracks and analysis-output files

  nohup ./knisbacher_script_step_3_create_tracks.sh LTR GA &
  nohup ./knisbacher_script_step_3_create_tracks.sh DNA GA &
  #Negative controls: To confirm that the G>A hypermutants identified were the product of APOBEC editing and not of random mutagenesis or computational or sequencing artefacts.
                      #The strand-biased nature of APOBEC DNA editing (with regard to the retroelement’s ORFs) enables to use C-to-T, the event complementary to G-to-A, as control.
  nohup ./knisbacher_script_step_3_create_tracks.sh LTR CT &

  #Run for all mismatches to look for more patterns & ensure the G-A enrichment pattern is unique
  for mismatch in GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  nohup ./knisbacher_script_step_3_create_tracks.sh LTR $mismatch &
  done

#QC: Quick check the main parts of output to ensure there are no errors
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_3_create_tracks LTR GA
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_3_create_tracks DNA GA
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_3_create_tracks LTR CT

  for mismatch in GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_3_create_tracks LTR $mismatch &
  done

#__________________________________________________________________________________________________________________________________________________________
#3.3. Apply pairwise & consensus filters on all cluster files

  nohup ./knisbacher_script_step_4_Filter_clusters.sh LTR GA &
  nohup ./knisbacher_script_step_4_Filter_clusters.sh DNA GA &
  nohup ./knisbacher_script_step_4_Filter_clusters.sh LTR CT &  #Control

  for mismatch in GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  nohup ./knisbacher_script_step_4_Filter_clusters.sh LTR $mismatch &
  done

#QC: Quick check the main parts of output to ensure there are no errors
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_4_Filter_clusters LTR GA
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_4_Filter_clusters DNA GA
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_4_Filter_clusters LTR CT

  for mismatch in GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_4_Filter_clusters LTR $mismatch &
  done
  
#__________________________________________________________________________________________________________________________________________________________
#3.4. Apply best sources filter & find motifs & parses filtered cluster files and creates several UCSC-tracks and analysis-output files.

  nohup ./knisbacher_script_step_5_create_filtered_tracks.sh LTR GA &
  nohup ./knisbacher_script_step_5_create_filtered_tracks.sh DNA GA &
  nohup ./knisbacher_script_step_5_create_filtered_tracks.sh LTR CT &  #Control

  for mismatch in GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  nohup ./knisbacher_script_step_5_create_filtered_tracks.sh LTR $mismatch &
  done
  
#QC: Quick check the main parts of output to ensure there are no errors
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_5_create_filtered_tracks LTR GA
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_5_create_filtered_tracks DNA GA
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_5_create_filtered_tracks LTR CT

  for mismatch in GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_5_create_filtered_tracks LTR $mismatch &
  done

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.5. Calculate background mutation rate (since insertion of repeat in genome) & total edit sites

  ./knisbacher_script_step_6_Calculate_total_edit_sites.sh LTR GA
  ./knisbacher_script_step_6_Calculate_total_edit_sites.sh DNA GA
  ./knisbacher_script_step_6_Calculate_total_edit_sites.sh LTR CT

  for mismatch in GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  nohup ./knisbacher_script_step_6_Calculate_total_edit_sites.sh LTR $mismatch &
  done

#QC: Quick check the main parts of output to ensure there are no errors
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_6_Calculate_total_edit_sites LTR GA
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_6_Calculate_total_edit_sites DNA GA
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_6_Calculate_total_edit_sites LTR CT

  for mismatch in GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  ./glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_6_Calculate_total_edit_sites LTR $mismatch &
  done

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

mkdir /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/stdout

#Move outputs used as control to a seperate folder
mkdir -p /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/control/LTR_other_mismatches
ls glimpse_output_of_knisbacher_script_step_*_81_birds | egrep -v "GA|CT|step_1|step_2" | xargs -n1 sh -c 'mv $0 /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/control/LTR_other_mismatches'

mkdir -p /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/control/Bird_DNA_transposons


mkdir -p /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/control/LINE

##############################################################################################################################################################################################################################################################################################################
#View & compare step 6 results

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method

for final in glimpse_output_of_knisbacher_script_step_6_Calculate_total_*_edit_sites_81_birds
do
j0=$(echo $final | grep "[A-Z][A-Z]" -o)
j1=$(grep -v "Total_.*_edit_sites" $final | awk '{a+=$2} END{print a}')  #Cumulative sum of total edit sites from all species
j2=$(grep -v "Total_.*_edit_sites" $final | awk '$2>0{print$2}' | ministat -n 2> /dev/null | awk 'END{$1=""; print$0}') #Stats like N, min, max, med, avg, stddev
if [[ -z $j2 ]]; then j2="0 0 0 0 0 0"; else :; fi
j3=$(grep -v "Total_.*_edit_sites" $final | awk '$3!="-"{$1=$2=""; print$0}' OFS="\n" | awk NF | cut -f1 -d "(" | sort | uniq -c | awk '{print$2,$1}' | sort -k2,2nr | awk 'NR==1{print$0} NR>1{print". - - - - - - - "$0}')  #Subfamily enriched: no of species with giving edit sites in each subfamily
j4=$(grep -v "Total_.*_edit_sites" $final | awk '$3!="-"{$1=$2=""; print$0}' OFS="\n" | awk NF | tr "()," " " | awk '{print$1,$NF}' | awk '{sum[$1] += $2} END {for (i in sum) print i, sum[i]}' | sort -k2,2nr) #Cumulative count of edit sites each subfamily 
l=$(echo -e "$j4" | wc -l)  #total no of subfamily found considering all species
j5=$(grep -v "Total_.*_edit_sites" $final | awk '$2>999 {print$1,$2}' | sort -k2,2nr | awk -v a="$l" 'NR<=a{print$0} NR> a {print". - - - - - - - - - - - "$0}')  #Spiecies with edit count more than 1000
paste <(echo $j0 $j1 $j2) <(echo -e "$j3") <(echo -e "$j4") <(echo -e "$j5")
unset j0 j1 j2 j3 j4 j5
done | sed '1i Mismatch Total_sites N Min Max Median Avg Stddev subfam- species_enriched Subfam- total_edit_sites Species_with_high_count- Count' | column -t | GREP_COLORS="mt=01;04;33" grep "Mismatch.*\|$" --color=always > compare_total_edit_sites_of_all_mismatches

paste glimpse_output_of_knisbacher_script_step_5_create_filtered_tracks_GA_81_birds glimpse_output_of_knisbacher_script_step_5_create_filtered_tracks_AG_81_birds | sed 's/_clusters//g' | sed 's/filtered//g' | column -t | less -SN
#CT edit sites are less, but they motifs are exact reverse complement of GA: it is possible that many of the C>T mutations are actually caused by APOBEC-mediated editing
paste glimpse_output_of_knisbacher_script_step_5_create_filtered_tracks_GA_81_birds glimpse_output_of_knisbacher_script_step_5_create_filtered_tracks_CT_81_birds | sed 's/_clusters//g' | sed 's/filtered//g' | column -t | less -SN

##############################################################################################################################################################################################################################################################################################################
#DRAFT
##############################################################################################################################################################################################################################################################################################################

#Print folder tree heirarchy of each bird main folder
cd $base_dir
while read -r i; do
cd "$i"/knisbacher
echo ">"$i && pwd | sed 's/^/  /g' | GREP_COLORS="mt=01;34" grep ".*" --color=always
tree -d | grep -v "^.$" | awk NF | sed 's/^/    /g' | GREP_COLORS="mt=01;36" grep ".*" --color=always
echo -e
cd $base_dir
done < all_bird_genomes_used | awk '/^>/{print ++i,$0;next} {print$0}' | less -SR

#Print all major outputs from running each steps
paste /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/QC/glimpse_output_of_knisbacher_script_step_* | column -t | less -S

#Print Clusters that became empty after 2nd consensus filtering
awk '{if ($1 == "-") $1 = prev; else prev = $1}1' results/QC/glimpse_output_of_knisbacher_script_step_4_Filter_clusters_81_birds | awk '$NF=="0" && $3!="0"' | column -t

for track in $(find Data/"$o"/LTR/results/Tracks/ -maxdepth 1 -mindepth 1 -type d)
do
subfam=$(echo $track | awk -F "/" '{print$6}' | awk -F "_" '{print$4}')
size=$(find $track -path "*/GA/raw/filteredClusts/*" -name "clusters_*.tab" | xargs -n1 bash -c 'paste <(ls -lh $0 | awk "{print\$5}")' | paste -s -d " ")
echo "-" $track $subfam $size
unset subfam size
done | awk -v a="$i" 'NR==1{$1=a}1'

#Look for error messages in standard output & standard error from each step
cd $base_dir
while read -r i; do
cd "$i"/knisbacher
echo ">"$i | sed 's/^/  /g' | GREP_COLORS="mt=01;34" grep ".*" --color=always
cat stdout_knisbacher_script_step_1_prepare_data
cat stdout_knisbacher_script_step_2_runClusterFinder
cat stdout_knisbacher_script_step_3_create_tracks
cat stdout_knisbacher_script_step_3_create_tracks_for_CT_clusters
cat stdout_knisbacher_script_step_4_Filter_clusters
stdout_knisbacher_script_step_4_Filter_clusters_of_CT
cat stdout_knisbacher_script_step_5_create_filtered_tracks 
cd $base_dir
done < all_bird_genomes_used | awk NF \
| egrep -v "^ESTIMATE| - Create Filtered Tracks:| - Run createTrackFiles2.pl| - Apply best sources filter:| - No filtered clusters to apply best sources filter|^real|^user|^sys|    - ERV" \
| egrep -v " Prepare Data| Run rmskOutToInterval.pl| Run Genome2Fasta.pl| Run sortGenome.pl|DONE" \
| egrep -v " - Find clusters:|   - Run runCluste|Building a new DB,|New DB name: |New DB title:  Dat|Sequence type: Nuc|Keep MBits: T|Maximum file size:|Adding sequences" | less -SR

#Print folders
cd $base_dir
while read  i; do
o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
if [[ -d "$i"/knisbacher/Data/"$o" ]]; then cd "$i"/knisbacher/Data/"$o"/; j=$(ls); else j="-"; fi
echo $i $j
unset o j
cd $base_dir
done < all_bird_genomes_used | column -t | less -SN

#Delete DNA transposons & LINE results
cd $base_dir
while read  i; do
echo ">"$i
o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
if [[ -d "$i"/knisbacher/Data/"$o" ]]; then cd "$i"/knisbacher/Data/"$o"/; rm -rv DNA LINE; else :; fi
unset o
cd $base_dir
done < all_bird_genomes_used 

nohup ./dna_transposons_as_control_knisbacher_script_step_1_prepare_data.sh > stdout_knisbacher_script_step_1_prepare_data_dna_transposons_as_control_81_birds 2>&1 &
nohup ./line_as_control_knisbacher_script_step_1_prepare_data.sh > stdout_knisbacher_script_step_1_prepare_data_line_as_control_81_birds 2>&1 &

nohup ./knisbacher_script_step_2_find_clusters.sh DNA &

#print strand info
cd $base_dir
while read  i; do
o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
if [[ -d "$i"/knisbacher/Data/"$o"/LTR/results/blasts ]]; then
j=$(find "$i"/knisbacher/Data/"$o"/LTR/results/blasts/ -name "*.gz" | xargs -n1 sh -c 'zcat $0 | grep "Strand=" | sort | uniq -c | sed "s/Strand=//g" | sed "s/Plus/+/g"' | awk '{sum[$2]+=$1; if(user[$2]=="") {user[$2]=$3}} END{for(idx in sum) {print sum[idx],idx,user[idx]}}')
else j="0"
fi
echo $i $j
unset o j
cd $base_dir
done < <(head -5 all_bird_genomes_used) | column -t | less -SN

