##################################################################################################################################################################################################################################################################################################################
#                                                                                                                                     Compare functional impact in apobec1 b/w birds & mammals
##################################################################################################################################################################################################################################################################################################################

mkdir /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison

##################################################################################################################################################################################################################################################################################################################
#Check functional impact of amino acid mutations

#PROVEAN

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Setting up resources
  
 #Download in ceglab8 directly from sourceforge & transfer to ceglab25:
  scp provean-1.1.5.tar.gz ceglab25@172.30.1.131:/media/aswin/programs/
  tar -xvzf provean-1.1.5.tar.gz
  cd /media/aswin/programs/provean-1.1.5
  ./configure
  make
  sudo make install

 #Download protein nr database
  cd /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison
  #Using wget OR aria2c
  #wget -r -nd -A "nr.*.tar.gz" ftp://ftp.ncbi.nlm.nih.gov/blast/db/
  wget -r -nd -A "nr.*.tar.gz" ftp://ftp.ncbi.nlm.nih.gov/blast/db/v4/
  #Total wall clock time: 7h 0m 39s, Downloaded: 123 files, 296G in 6h 58m 13s (12.1 MB/s)
  #151m53.782s
  cat nr*.tar.gz | tar -zxvi -f - -C .

  #other way to download
  #aria2c -x 10 -s 10 -i urls.txt

#Install specific version of version of blast (based on this forum: https://www.biostars.org/p/443248/)
  cd /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/database/older_blast
  wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-x64-linux.tar.gz
  tar -xvzf ncbi-blast-2.4.0+-x64-linux.tar.gz
  #Add this path of blastdbsmd in /usr/local/bin/provean.sh
  
 #Specified paths inside /usr/local/bin/provean.sh
  #BLAST_DB="/media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/database/nr"
  #PSIBLAST="/home/ceglab25/miniconda3/bin/psiblast"
  #CD_HIT="/usr/local/bin/cd-hit"
  #BLASTDBCMD="/media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/database/older_blast/ncbi-blast-2.4.0+/bin/blastdbcmd"

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run provean

cd /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison
#26m20.839s
time provean.sh -q human.aa -v all_human.var


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run for all mutations observed in bird mammal A1 protein alignment

#create mutations list
  #In neo
  cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cds_alignment
  count=0
  while read i
  do
  count=$(calc $count + 1)
  #echo ">"$i $count
  sed '1,3d;241d' compare_residues_43_birds_32_mammals | sed '54d' | sed -n "$count p" | awk '{print$3,$4}' | tr "," " " | tr " " "\n" | cut -f1 -d "(" | grep -v "$i" | sed "s/^/$i$count/g" | sed 's/-/del/g' | tr -d "\t" | grep -v "X"
  done < <(grep Homo -A1 birds_mammals_filtered_8.aa | tail -1 | fold -1 ) > all_bird_mammals_A1_mutations_wrt_human.var

#Run provean
  #In ceglab25
  mkdir /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/all_bird_mammals_A1_mutations_wrt_human
  cd /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/all_bird_mammals_A1_mutations_wrt_human
  scp neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cds_alignment/all_bird_mammals_A1_mutations_wrt_human.var .
  cp ../human.aa .
  split -n l/32 -d -a2 --additional-suffix=.var all_bird_mammals_A1_mutations_wrt_human.var chunk_
  time ls chunk_*.var | parallel -j 32 'provean.sh -q human.aa -v {} > results_{}.out'
  cat results_chunk_*.out > combined_provean_results.txt
  provean.sh -q human.aa -v all_human.var --save_supporting_set human.sss

#Results
  ls results_chunk_*.var.out | xargs -n1 sh -c 'awk "/# VARIATION/,/^$/" $0' | grep -v "# VARIATION" | awk '{if($2<=-2.5) print$0,"Deleterious"; else print$0,"Neutral"}' | column -t > combined_scores
  #Look at distribution of provean scores
  sort -k2,2nr combined_scores | awk '{print$2}' | ministat

scp ../all_bird_mammals_A1_mutations_wrt_human neo@172.30.1.174:~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cds_alignment/provean/

#Provean results on bird-mammals A1 amino acid alignment table
  #In neo
  cd ~/bird_db1/aswin/APOBEC1/main_figures/mammal_extension/cds_alignment
  while read i
  do
  i1=$(echo $i | awk '{print$2}')
  j1=$(awk -v a="$i1" '$2==a' provean/all_bird_mammals_A1_mutations_wrt_human/deleterious_mutants | awk '{print$3}' | awk '{print$1"([0-9]*)"}' | sed 's/^-/\\-/g')
  if [[ -z $j1 ]]
  then
  echo $i
  else
  echo $i | grep "$j1" --color=always
  fi
  unset i1 j1
  done < compare_residues_43_birds_32_mammals | column -t > provean/all_bird_mammals_A1_mutations_wrt_human/compare_residues_43_birds_32_mammals_provean_deleterious_coloured


##################################################################################################################################################################################################################################################################################################################
#Sources
##################################################################################################################################################################################################################################################################################################################

#Tutorail:
  #how to install & run PROVEAN: https://gist.github.com/darencard/d3ff1399f1fa05b35429eb148383aea4

#About PROVEAN SCORING
  #Based on the PROVEAN score, the program reports a predicted functional category, either deleterious or neutral, based on a pre-set threshold.
  #Though it is possible for a mutant protein to receive a higher mean alignment score than the wild type, there is no category for beneficial effects.
  #The default threshold value is 2.5, and variants with scores below this threshold are classified as deleterious.
  #This threshold was chosen to maximize sensitivity (detection) and specificity (accuracy) when assigning functional effects to common versus disease-causing human protein variants (Choi et al. 2012; Choi and Chan 2015).
  #The creators of PROVEAN suggest adjusting the threshold score for defining deleterious alleles depending on the user’s needs.
  #However, in the studies summarized in supplementary table S1, Supplementary Material online, that applied PROVEAN to nonhuman organisms, nearly all used the default cutoff value.
  #Although there have been studies evaluating different methods of functional annotation of mutations, they do not consider the possibility of species-specific thresholds (Kono et al. 2018).
  #How well does this tool, which was developed with data from humans, work in predicting fitness effects in other organisms, when many variants are present?


