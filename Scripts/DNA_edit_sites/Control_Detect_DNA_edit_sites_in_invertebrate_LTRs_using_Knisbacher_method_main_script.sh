###########################################################################################################################################################################################################################################################################################################
#                                                                                                               NEGATIVE CONTROL: ESTIMATE APOBEC INDUCED DNA EDITED SITES IN INVERTEBRATE LTR RETROTRANSPOSONS
###########################################################################################################################################################################################################################################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. DOWNLOAD DATA

#Main folder
base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/invertebrates"
cd $base_dir

#Use the same modified script used for birds
script_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/scripts/Perl_DNAE"

# set perl version
perlbrew use perl-5.16.3


#__________________________________________________________________________________________________________________________________________________________
#1.1. Manually download genomes & respective repeatmasker tables from UCSC

#Some invertebrates had genome sequences and RepeatMasker output files in seperate individual files for each chromosomes
#Concatenate these into single file
cd $base_dir
for i in $(find . -name "*Out.tar.gz" | cut -f2 -d "/")
do
echo $i
cd $i
mkdir ucsc
#Combine genome fasta
o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
ls *.fa | sort -V | xargs -n1 sh -c 'cat $0' > ucsc/$o".fa"
#Combine rRepeatMasker output files
find . -name "*.fa.out" | grep -v "\./ucsc\/.*\.fa\.out" | sort -V | xargs -n1 sh -c 'cat $0' > ucsc/$o".fa.out"
awk -iinplace '!a[$0]++' ucsc/$o".fa.out"
unset o
cd ../
done

#Save all loist of species in a file
ls -d */ | egrep -v "all|knisbacher" | tr -d "/" > all_invertebrate_genomes 

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
  cd $base_dir
  start_time=$(date +%s) && max_jobs=32 && counter=0

  #Run script knisbacher_script_step_1_prepare_data.sh for 32 birds in parallel background jobs
  counter=0
  for species in $(cat all_invertebrate_genomes)
    do
    {
    echo $species > "split_"$species
    nohup ./invertebrates_knisbacher_script_step_1_prepare_data.sh "split_"$species > "stdout_invertebrates_knisbacher_script_step_1_prepare_data_"$species 2>&1 &
    counter=$((counter + 1))
    echo ">Species: "$species "| Background job: " $counter
    }
    # Check if the number of running jobs has reached max_jobs, if so, wait for one or more to finish
    check_jobs
  done
  wait
  end_time=$(date +%s) && elapsed_time=$((end_time - start_time))  #total time taken

  echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' | sed '$a\\n' > stdout_16_invertebrates_knisbacher_script_step_1_prepare_data
  cat stdout_invertebrates_knisbacher_script_step_1_prepare_data_* >> stdout_16_invertebrates_knisbacher_script_step_1_prepare_data
  rm split_* stdout_invertebrates_knisbacher_script_step_1_prepare_data_* 
  unset species start_time end_time elapsed_time max_jobs counter

#QC: Quick check the main parts of output to ensure there are no errors
  mkdir -p /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/invertebrates/results/summary/LTR_other_mismatches
  ./invertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_1_prepare_data LTR results/summary

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. Main analysis

#__________________________________________________________________________________________________________________________________________________________
# 3.1. Runs BLAST, finds clusters in the BLAST output and parses the cluster, for each Name LTR subfamily of sequences of the organism.

  nohup ./invertebrates_knisbacher_script_step_2_find_clusters.sh LTR &

#QC: Quick check the main parts of output to ensure there are no errors
  ./invertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_2_find_clusters LTR results/summary

#__________________________________________________________________________________________________________________________________________________________
#3.2. Parses cluster files and creates several UCSC-tracks and analysis-output files

  for mismatch in GA CT GC GT CA TA AG TC CG TG AC AT
  do
  {
  echo ">"$mismatch
  nohup ./invertebrates_knisbacher_script_step_3_create_tracks.sh LTR $mismatch &
  }
  done

  #QC: Quick check the main parts of output to ensure there are no errors
  ./invertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_3_create_tracks LTR GA results/summary
  for mismatch in CT GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  ./invertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_3_create_tracks LTR $mismatch results/summary/LTR_other_mismatches
  done

#__________________________________________________________________________________________________________________________________________________________
#3.3. Apply pairwise & consensus filters on all cluster files

  for mismatch in GA CT GC GT CA TA AG TC CG TG AC AT
  do
  {
  echo ">"$mismatch
  nohup ./invertebrates_knisbacher_script_step_4_Filter_clusters.sh LTR $mismatch &
  }
  done

 #QC: Quick check the main parts of output to ensure there are no errors
  ./invertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_4_Filter_clusters LTR GA results/summary
  for mismatch in CT GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  ./invertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_4_Filter_clusters LTR $mismatch results/summary/LTR_other_mismatches
  done

#__________________________________________________________________________________________________________________________________________________________
#3.4. Apply best sources filter & find motifs & parses filtered cluster files and creates several UCSC-tracks and analysis-output files.

  for mismatch in GA CT GC GT CA TA AG TC CG TG AC AT
  do
  {
  echo ">"$mismatch
  nohup ./invertebrates_knisbacher_script_step_5_create_filtered_tracks.sh LTR $mismatch &
  }
  done
  
#QC: Quick check the main parts of output to ensure there are no errors
  ./invertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_5_create_filtered_tracks LTR GA results/summary
  for mismatch in CT GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  ./invertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_5_create_filtered_tracks LTR $mismatch results/summary/LTR_other_mismatches
  done

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.5. Calculate background mutation rate (since insertion of repeat in genome) & total edit sites

  for mismatch in GA CT GC GT CA TA AG TC CG TG AC AT
  do
  {
  echo ">"$mismatch
  nohup ./invertebrates_knisbacher_script_step_6_Calculate_total_edit_sites.sh LTR $mismatch &
  }
  done

#QC: Quick check the main parts of output to ensure there are no errors
  ./invertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_6_Calculate_total_edit_sites LTR GA results/summary
  for mismatch in CT GC GT CA TA AG TC CG TG AC AT
  do
  echo ">"$mismatch
  ./invertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_6_Calculate_total_edit_sites LTR $mismatch results/summary/LTR_other_mismatches
  done

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mkdir /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/stdout

#Move outputs used as control to a seperate folder
mkdir -p /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/control/LTR_other_mismatches
ls glimpse_output_of_knisbacher_script_step_*_81_birds | egrep -v "GA|CT|step_1|step_2" | xargs -n1 sh -c 'mv $0 /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/control/LTR_other_mismatches'
############################################################################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################################################################
#DRAFT SCRIPTS

#Collect all inputs in a single folder called "ucsc"
for x in $(for i in `cat all_invertebrate_genomes`; do if [ -d $i/ucsc ]; then :; else echo $i; fi; done)
do
echo $x
cd $x
mkdir ucsc
#Combine genome fasta
o=$(echo $x | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
ls *.fa | sort -V | xargs -n1 sh -c 'cat $0' > ucsc/$o".fa"
#Combine rRepeatMasker output files
find . -name "*.fa.out" | grep -v "\./ucsc\/.*\.fa\.out" | sort -V | xargs -n1 sh -c 'cat $0' > ucsc/$o".fa.out"
awk -iinplace '!a[$0]++' ucsc/$o".fa.out"
unset o
cd ../
done

#Print input folder
for i in `cat all_invertebrate_genomes`
do
j=`ls $i | grep ucsc`
echo $i $j
unset j
done | column -t | nl

#Print main results folder
for i in `cat all_invertebrate_genomes`
do
j=`ls $i | grep knisbacher`
if [[ -z $j ]]
then
j="-"
else :
fi
echo $i $j
unset j
done | column -t | nl

#Print folder heirarchy of main result folder
cd $base_dir
while read -r i; do
cd "$i"/knisbacher
echo ">"$i && pwd | sed 's/^/  /g' | GREP_COLORS="mt=01;34" grep ".*" --color=always
tree -d | grep -v "^.$" | awk NF | sed 's/^/    /g' | GREP_COLORS="mt=01;36" grep ".*" --color=always
echo -e
cd $base_dir
done < all_invertebrate_genomes | awk '/^>/{print ++i,$0;next} {print$0}' | less -SR

#Delete main results folder
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/invertebrates
for i in `cat all_invertebrate_genomes`
do
j=`ls $i | grep knisbacher`
if [[ -z $j ]]
then
j="-"
else
cd $i
rm -r knisbacher
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/invertebrates
fi
echo $i $j
unset j
done | column -t | nl

#Print time taken for running each steps 
ls stdout_16_invertebrates_knisbacher_script_step_4_Filter_clusters_of_LTR_* | xargs -n1 bash -c 'paste <(echo $0) <(tail -3 $0 | head -1 | awk "{print\$NF}")' | sort -k2,2V


