###########################################################################################################################################################################################################################################################################################################
#                                                                                                               NEGATIVE CONTROL: ESTIMATE APOBEC INDUCED DNA EDITED SITES IN NON-MAMMALIAN VERTEBRATES LTR RETROTRANSPOSONS
###########################################################################################################################################################################################################################################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. DOWNLOAD DATA

#Main folder
base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/non_mammalian_vertebrates"
cd $base_dir

ls -d */ | egrep -v "^all|knisbacher" | tr -d "/" > all_nonmammalian_vertebrates_genomes 

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
rename 's/invertebrates/nonmammalian_vertebrates/g' *.sh
sed 's/invertebrates/nonmammalian_vertebrates/g' nonmammalian_vertebrates_*.sh -i
sed 's/all_invertebrate_genomes/all_nonmammalian_vertebrates_genomes/g' nonmammalian_vertebrates_*.sh -i
sed 's!/nonmammalian_vertebrates!/non_mammalian_vertebrates!g' *.sh -i
sed 's/stdout_16_nonmammalian_vertebrates/stdout_10_nonmammalian_vertebrates/g' nonmammalian_vertebrates_*.sh -i
sed 's/INVERTEBRATE/NON-MAMMALIAN VERTEBRATE/g' nonmammalian_vertebrates_*.sh -i
sed 's/81_birds/10_nonmammalian_vertebrates/g' nonmammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh -i

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
  for species in $(cat all_nonmammalian_vertebrates_genomes)
    do
    {
    echo $species > "split_"$species
    nohup ./nonmammalian_vertebrates_knisbacher_script_step_1_prepare_data.sh "split_"$species > "stdout_nonmammalian_vertebrates_knisbacher_script_step_1_prepare_data_"$species 2>&1 &
    counter=$((counter + 1))
    echo ">Species: "$species "| Background job: " $counter
    }
    # Check if the number of running jobs has reached max_jobs, if so, wait for one or more to finish
    check_jobs
  done
  wait
  end_time=$(date +%s) && elapsed_time=$((end_time - start_time))  #total time taken

  echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' | sed '$a\\n' > stdout_10_nonmammalian_vertebrates_knisbacher_script_step_1_prepare_data
  cat stdout_nonmammalian_vertebrates_knisbacher_script_step_1_prepare_data_* >> stdout_10_nonmammalian_vertebrates_knisbacher_script_step_1_prepare_data
  rm split_* stdout_nonmammalian_vertebrates_knisbacher_script_step_1_prepare_data_* 
  unset species start_time end_time elapsed_time max_jobs counter

#QC: Quick check the main parts of output to ensure there are no errors
  mkdir -p /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/non_mammalian_vertebrates/results/summary/LTR_other_mismatches
  ./nonmammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_1_prepare_data LTR results/summary

#Since 2 genomes has no LTR remove them from further analysis
sed -e '/Callorhinchus_milii/d' -e '/Pleurodeles_waltl/d' -e '/Latimeria_chalumnae/d' all_nonmammalian_vertebrates_genomes -i
sed 's/stdout_10_nonmammalian_vertebrates/stdout_7_nonmammalian_vertebrates/g' nonmammalian_vertebrates_*.sh -i

#grep -v "#" common_names | awk '{print$1}' > all_nonmammalian_vertebrates_genomes_filtered
#awk '{print$1}' common_names | tr -d "#" > all_nonmammalian_vertebrates_genomes
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. Main analysis

#__________________________________________________________________________________________________________________________________________________________
# 3.1. Runs BLAST, finds clusters in the BLAST output and parses the cluster, for each Name LTR subfamily of sequences of the organism.

  nohup ./nonmammalian_vertebrates_knisbacher_script_step_2_find_clusters.sh LTR &

#QC: Quick check the main parts of output to ensure there are no errors
  ./nonmammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_2_find_clusters LTR results/summary

#__________________________________________________________________________________________________________________________________________________________
#3.2. Parses cluster files and creates several UCSC-tracks and analysis-output files

  for mismatch in GA CT
  do
  {
  echo ">"$mismatch
  nohup ./nonmammalian_vertebrates_knisbacher_script_step_3_create_tracks.sh LTR $mismatch &
  }
  done

  #QC: Quick check the main parts of output to ensure there are no errors
  for mismatch in GA CT
  do
  echo ">"$mismatch
  ./nonmammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_3_create_tracks LTR $mismatch results/summary
  done

#__________________________________________________________________________________________________________________________________________________________
#3.3. Apply pairwise & consensus filters on all cluster files

  for mismatch in GA CT 
  do
  {
  echo ">"$mismatch
  nohup ./nonmammalian_vertebrates_knisbacher_script_step_4_Filter_clusters.sh LTR $mismatch &
  }
  done

 #QC: Quick check the main parts of output to ensure there are no errors
  for mismatch in GA CT
  do
  echo ">"$mismatch
  ./nonmammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_4_Filter_clusters LTR $mismatch results/summary
  done

#__________________________________________________________________________________________________________________________________________________________
#3.4. Apply best sources filter & find motifs & parses filtered cluster files and creates several UCSC-tracks and analysis-output files.

  for mismatch in GA CT
  do
  {
  echo ">"$mismatch
  nohup ./nonmammalian_vertebrates_knisbacher_script_step_5_create_filtered_tracks.sh LTR $mismatch &
  }
  done
  
#QC: Quick check the main parts of output to ensure there are no errors
  for mismatch in GA CT
  do
  echo ">"$mismatch
  ./nonmammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_5_create_filtered_tracks LTR $mismatch results/summary
  done

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.5. Calculate background mutation rate (since insertion of repeat in genome) & total edit sites

  for mismatch in GA CT
  do
  {
  echo ">"$mismatch
  nohup ./nonmammalian_vertebrates_knisbacher_script_step_6_Calculate_total_edit_sites.sh LTR $mismatch &
  }
  done

#QC: Quick check the main parts of output to ensure there are no errors
  for mismatch in GA CT
  do
  echo ">"$mismatch
  ./nonmammalian_vertebrates_glimpse_output_of_knisbacher_script_steps.sh knisbacher_script_step_6_Calculate_total_edit_sites LTR $mismatch results/summary
  done

###########################################################################################################################################################################################################################################################################################################
#Rerun

#Anolis_carolinensis is a species where we expect DNA editing, but not observed as there were fo ERVs in the repeatmasker output file from UCSC, hence run  it locally with updated genome
#Add 2 more species available in UCSC

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/non_mammalian_vertebrates/Anolis_carolinensis/local
#417m16.685s
time /home/ceglab25/ajs/maker/exe/RepeatMasker/RepeatMasker -pa 32 -engine RMBlast -a -s -gff -species "Anolis carolinensis" GCF_035594765.1_rAnoCar3.1.pri_genomic.fna
cp GCF_035594765.1_rAnoCar3.1.pri_genomic.fna Anocar.fa
cp GCF_035594765.1_rAnoCar3.1.pri_genomic.fna.out Anocar.fa.out

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/non_mammalian_vertebrates/Anolis_carolinensis
mv knisbacher ucsc_repeatmasker_knisbacher

base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/non_mammalian_vertebrates"
cd $base_dir



