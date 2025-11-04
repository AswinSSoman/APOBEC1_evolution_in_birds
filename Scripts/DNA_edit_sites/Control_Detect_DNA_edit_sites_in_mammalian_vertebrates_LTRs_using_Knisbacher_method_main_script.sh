###########################################################################################################################################################################################################################################################################################################
#                                                                                                               NEGATIVE CONTROL: ESTIMATE APOBEC INDUCED DNA EDITED SITES IN MAMMALIAN VERTEBRATES LTR RETROTRANSPOSONS
###########################################################################################################################################################################################################################################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. DOWNLOAD DATA

#Main folder
base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/mammalian_vertebrates"
cd $base_dir

#Few examples of how to download repeatmasker output & associated genomes
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz ./
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.out.gz ./
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz ./
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz ./
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/oryCun2/bigZips/oryCun2.fa.gz ./
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/oryCun2/bigZips/oryCun2.fa.out.gz ./
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/ornAna2/bigZips/ornAna2.fa.gz ./
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/ornAna2/bigZips/ornAna2.fa.out.gz ./
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/monDom5/bigZips/chromFa.tar.gz ./
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/monDom5/bigZips/chromOut.tar.gz ./
find . -name "*chr*.fa" | xargs -n1 sh -c 'cat $0' > Mondom.fa
find . -name "*chr*.fa.out" | xargs -n1 sh -c 'cat $0' > Mondom.fa.out

ls -d */ | egrep -v "^all|knisbacher" | tr -d "/" > all_mammalian_vertebrates_genomes 

#Steps to download :
#Go to UCSC Download section (https://hgdownload.soe.ucsc.edu/downloads.html) > Search Species > Click "Genome sequence files and select annotations (2bit, GTF, GC-content, etc)" > Locate & copy the links of genome & repeatmasker output file (remove "https" at the starting of the link) > Download in commandline using rsync
#Eg: Alligator_mississippiensis
#rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/allMis1/bigZips/allMis1.fa.gz ./
#rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/allMis1/bigZips/allMis1.fa.out.gz ./
#gzip -d allMis1.fa.*
#for i in `ls -d */`; do cd $i; echo $i; mkdir ucsc; mv *.fa *.out ucsc/; cd ../; done

#Use the same modified script used for birds
script_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/scripts/Perl_DNAE"

# set perl version
perlbrew use perl-5.16.3

#Transfer execution ready scripts & add group specific modifications
cp ../invertebrates/*.sh .
rename 's/invertebrates/mammalian_vertebrates/g' *.sh
sed 's/invertebrates/mammalian_vertebrates/g' mammalian_vertebrates_*.sh -i
sed 's/all_invertebrate_genomes/all_mammalian_vertebrates_genomes/g' mammalian_vertebrates_*.sh -i
#sed 's!/mammalian_vertebrates!/mammalian_vertebrates!g' *.sh -i
sed 's/stdout_16_mammalian_vertebrates/stdout_14_mammalian_vertebrates/g' mammalian_vertebrates_*.sh -i
sed 's/INVERTEBRATE/MAMMALIAN VERTEBRATE/g' mammalian_vertebrates_*.sh -i

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
  for species in $(cat all_mammalian_vertebrates_genomes)
    do
    {
    echo $species > "split_"$species
    nohup ./mammalian_vertebrates_knisbacher_script_step_1_prepare_data.sh "split_"$species > "stdout_mammalian_vertebrates_knisbacher_script_step_1_prepare_data_"$species 2>&1 &
    counter=$((counter + 1))
    echo ">Species: "$species "| Background job: " $counter
    }
    # Check if the number of running jobs has reached max_jobs, if so, wait for one or more to finish
    check_jobs
  done
  wait
  end_time=$(date +%s) && elapsed_time=$((end_time - start_time))  #total time taken

  echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' | sed '$a\\n' > stdout_9_mammalian_vertebrates_knisbacher_script_step_1_prepare_data
  cat stdout_mammalian_vertebrates_knisbacher_script_step_1_prepare_data_* >> stdout_9_mammalian_vertebrates_knisbacher_script_step_1_prepare_data
  rm split_* stdout_mammalian_vertebrates_knisbacher_script_step_1_prepare_data_* 
  unset species start_time end_time elapsed_time max_jobs counter

#QC: Quick check the main parts of output to ensure there are no errors
  mkdir -p /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/mammalian_vertebrates/results/summary/LTR_other_mismatches
  ./mammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_1_prepare_data LTR results/summary

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. Main analysis

#__________________________________________________________________________________________________________________________________________________________
# 3.1. Runs BLAST, finds clusters in the BLAST output and parses the cluster, for each Name LTR subfamily of sequences of the organism.

  nohup ./mammalian_vertebrates_knisbacher_script_step_2_find_clusters.sh LTR &

#QC: Quick check the main parts of output to ensure there are no errors
  ./mammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_2_find_clusters LTR results/summary

#__________________________________________________________________________________________________________________________________________________________
#3.2. Parses cluster files and creates several UCSC-tracks and analysis-output files

  for mismatch in GA CT
  do
  {
  echo ">"$mismatch
  nohup ./mammalian_vertebrates_knisbacher_script_step_3_create_tracks.sh LTR $mismatch &
  }
  done

  #QC: Quick check the main parts of output to ensure there are no errors
  for mismatch in GA CT
  do
  echo ">"$mismatch
  ./mammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_3_create_tracks LTR $mismatch results/summary
  done

#__________________________________________________________________________________________________________________________________________________________
#3.3. Apply pairwise & consensus filters on all cluster files

  for mismatch in GA CT 
  do
  {
  echo ">"$mismatch
  nohup ./mammalian_vertebrates_knisbacher_script_step_4_Filter_clusters.sh LTR $mismatch &
  }
  done

 #QC: Quick check the main parts of output to ensure there are no errors
  for mismatch in GA CT
  do
  echo ">"$mismatch
  ./mammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_4_Filter_clusters LTR $mismatch results/summary
  done

#__________________________________________________________________________________________________________________________________________________________
#3.4. Apply best sources filter & find motifs & parses filtered cluster files and creates several UCSC-tracks and analysis-output files.

  for mismatch in GA CT
  do
  {
  echo ">"$mismatch
  nohup ./mammalian_vertebrates_knisbacher_script_step_5_create_filtered_tracks.sh LTR $mismatch &
  }
  done
  
#QC: Quick check the main parts of output to ensure there are no errors
  for mismatch in GA CT
  do
  echo ">"$mismatch
  ./mammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_5_create_filtered_tracks LTR $mismatch results/summary
  done

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.5. Calculate background mutation rate (since insertion of repeat in genome) & total edit sites

  for mismatch in GA CT
  do
  {
  echo ">"$mismatch
  nohup ./mammalian_vertebrates_knisbacher_script_step_6_Calculate_total_edit_sites.sh LTR $mismatch &
  }
  done

#QC: Quick check the main parts of output to ensure there are no errors
  for mismatch in GA CT
  do
  echo ">"$mismatch
  ./mammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_6_Calculate_total_edit_sites LTR $mismatch results/summary
  done


##


