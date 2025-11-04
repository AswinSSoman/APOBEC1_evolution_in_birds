###########################################################################################################################################################################################################################################################################################################
#Control for APOBEC associated DNA edit sites in birds LTR retrotransposons
###########################################################################################################################################################################################################################################################################################################

# 3 Negative controls used:
  #1. Estimate C-T edit sites in LTR (Strand specificity)
  #2. Estimate G-A edit sites in DNA transposons
  #3. Estimate G-A edit sites in invertebrate LTRs

#Concept:
  #DNA transposons are seemingly not targeted by APOBECs, thus should not contain strand-biased mutations, hence we searched for all types of clustered mutations in DNA transposons of the 80 genomes to validate this.
  #DNA transposons can serve as additional control. These elements transpose in the genome through a “cut-and-paste” mechanism that does not involve reverse transcription. 
  #As APOBECs specifically edit single-stranded DNA synthesized during reverse transcription, these elements are not expected to contain strand biased clusters of G-to-A mutations, characteristic of APOBEC editing.
  #If the number of elements bearing other types of clusters was greater than those with G>A shows weak signal.
  #This weak, nonspecific and nonstrand-biased signal implies that the few clusters found in DNA transposons are not associated with APOBECs, unlike the clusters found in LTRs.

###########################################################################################################################################################################################################################################################################################################

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons

base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"
script_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/scripts/Perl_DNAE"

#Root directory
cd $base_dir

#Set perl version
perlbrew use perl-5.16.3 &> /dev/null

#__________________________________________________________________________________________________________________________________________________________
#2.1. Sort & extract repeats based on class, family & subfamilies

#Inputs & scripts required to run this are already present in the directory from running the knisbacher method
#Hence just start running from sortgenome.pl

nohup ./dna_transposons_as_control_knisbacher_script_step_1_prepare_data.sh > stdout_knisbacher_script_step_1_prepare_data_dna_transposons_as_control_81_birds 2>&1 &

