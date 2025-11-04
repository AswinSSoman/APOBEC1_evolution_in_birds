###########################################################################################################################################################################################################################################################################################################
#Control for APOBEC associated DNA edit sites in birds LTR retrotransposons
###########################################################################################################################################################################################################################################################################################################

#The clusters of G-to-A mutations identified were strand-biased, appearing predominantly in the retroelements’ sense strand.
#This alone is sufficient to demonstrate that the G>A mutations aren’t sequencing errors, because WGS methods sequence both strands of DNA.

#The knisbacher method is repeated with a small change: change strand info in from "-strand plus" to "-strand minus" `runClusterFinder.pl`
sed 's/\-strand plus/\-strand minus/g' scripts/Perl_DNAE/runClusterFinder.pl -i

#Inputs & scripts required to run this are already present in the directory from running the knisbacher method
#Hence just start running from sortgenome.pl

base_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"
script_dir="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/scripts/Perl_DNAE"

#Root directory
cd $base_dir

#Set perl version
perlbrew use perl-5.16.3 &> /dev/null

nohup ./blast_minus_strand_as_control_knisbacher_script_step_1_prepare_data.sh > stdout_knisbacher_script_step_1_prepare_data_blast_minus_strand_as_control_81_birds 2>&1 &

nohup ./blast_minus_strand_as_control_knisbacher_script_step_2_find_clusters.sh &
