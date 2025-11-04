###########################################################################################################################################################################################################################################################################################################
#                                                                                                              PLOT CORRELATION BETWEEN GA EDIT SITES & APOBEC1 GENE PRESENCE/ABSENCE
###########################################################################################################################################################################################################################################################################################################

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons

#Get corrected total count of GA & CT edit sites
for i in $(awk '{print$1}' results/knisbacher_method/glimpse_output_of_knisbacher_script_step_6_Calculate_total_GA_edit_sites_81_birds)
do
j1=$(grep "$i" results/knisbacher_method/glimpse_output_of_knisbacher_script_step_6_Calculate_total_GA_edit_sites_81_birds | awk '{print$1,$2}')
j2=$(grep "$i" results/knisbacher_method/glimpse_output_of_knisbacher_script_step_6_Calculate_total_CT_edit_sites_81_birds | awk '{print$2}')
echo $j1 $j2
unset j1 j2
done > results/knisbacher_method/birds_total_GA_CT_edit_sites_count

#Convert count to percent
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method
awk 'NR>1{
    if ($2 == 0 && $3 > 0) {
        print $1, 0, 100
    } else if ($3 == 0 && $2 > 0) {
        print $1, 100, 0
    } else if ($2 == 0 && $3 == 0) {
        print $1, 0, 0
    } else {
        print $1, $2/($2+$3)*100, $3/($2+$3)*100
    }
}' birds_total_GA_CT_edit_sites_count > birds_total_GA_CT_edit_sites_percent

#Add gene status (loss/retain) info to the percent table
awk 'NR==FNR {file1[$1]; next} {print $1, ($1 in file1) ? 0 : 1, $2, $3}' /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/species_with_APOBEC1_loss birds_total_GA_CT_edit_sites_percent > birds_total_GA_CT_edit_sites_percent_with_loss_status
awk 'NR==FNR {file1[$1]; next} {print $1, ($1 in file1) ? 0 : 1, $2, $3}' /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/species_with_APOBEC1_loss birds_total_GA_CT_edit_sites_count > birds_total_GA_CT_edit_sites_count_with_loss_status

#Look at species with ZERO GA edit sites but has CT edit sites (if the difference is very low, then it might not have any biological meaning it is negligible)
awk 'NR==1{print$0,"GA-CT";next} NR>1{print$0,$2-$3 | "sort -k4,4 -nr"}' birds_total_GA_CT_edit_sites_count | column -t | less
cat birds_total_GA_CT_edit_sites_percent | awk '$2=="0" && $3!="0"'
cat birds_total_GA_CT_edit_sites_count | awk '$2=="0" && $3!="0"'
#Remove species with zero GA edit sites
awk '$3!="0"' birds_total_GA_CT_edit_sites_percent_with_loss_status > birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered

#Obtain the time calibrated phylogenetic tree
awk '{print$1}' birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered | paste -s -d "," > birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree_names
#scp neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/tree_update/list_of_species_to_keep_upupa_epops_removed.nwk .
#grep -wvif  <(cat list_of_species_to_keep_upupa_epops_removed.nwk | tr -d "():[0-9]'.;" | tr "," "\n" | sort) <(cat birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree_names | tr "," "\n" | sort)
sed  -e 's/Erythrura_gouldiae/Chloebia_gouldiae/g' -e 's/Aquila_chrysaetos_chrysaetos/Aquila_chrysaetos/g' -e 's/Limosa_lapponica_baueri/Limosa_lapponica/g' -e 's/Lonchura_striata_domestica/Lonchura_striata/g' birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree_names -i
sed  -e 's/Erythrura_gouldiae/Chloebia_gouldiae/g' -e 's/Aquila_chrysaetos_chrysaetos/Aquila_chrysaetos/g' -e 's/Limosa_lapponica_baueri/Limosa_lapponica/g' -e 's/Lonchura_striata_domestica/Lonchura_striata/g' birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered -i

#Prune tree keeping only species with non-zero GA edit sites 
Rscript /media/aswin/programs/prune_tree.r list_of_species_to_keep_upupa_epops_removed.nwk birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree_names birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree.nwk
sed "s/)'[0-9]\+'/)/g" birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree.nwk -i

#plot
time Rscript DNA_edit_site_gene_loss_correlation_plot.r birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree.nwk

############################################################################################################################################################################################################################################################################################################
#DRAFT

grep -wvif  <(cat list_of_species_to_keep_upupa_epops_removed.nwk | tr -d "():[0-9]'.;" | tr "," "\n" | sort) <(cat birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree_names | tr "," "\n" | sort)
grep -vif  <(cat list_of_species_to_keep_upupa_epops_removed.nwk | tr -d "():[0-9]'.;" | tr "," "\n" | sort) <(cat birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree_names | tr "," "\n" | sort)
grep -wvif  <(cat list_of_species_to_keep_upupa_epops_removed.nwk | tr -d "():[0-9]'.;" | tr "," "\n" | sort) <(cat birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree_names | tr "," "\n" | sort) | awk -F "_" '$3'
grep -wvif  <(cat list_of_species_to_keep_upupa_epops_removed.nwk | tr -d "():[0-9]'.;" | tr "," "\n" | sort) <(cat birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree_names | tr "," "\n" | sort) | awk -F "_" '$3' \
    | awk -F "_" '{print"s/"$0"/"$1"_"$2"/g"}' | sed "s/^/-e '/g" | sed "s/$/'/g" | paste -s -d " " | sed 's/^/sed /g'

sed -e 's/Erythrura_gouldiae/Chloebia_gouldiae/g' -e 's/Stachyris_ruficeps/Cyanoderma_ruficeps/g' -e 's/Aquila_chrysaetos_chrysaetos/Aquila_chrysaetos/g' -e 's/Buceros_rhinoceros_silvestris/Buceros_rhinoceros/g' -e 's/Limosa_lapponica_baueri/Limosa_lapponica/g' -e 's/Lonchura_striata_domestica/Lonchura_striata/g' -e 's/Malurus_cyaneus_samueli/Malurus_cyaneus/g' -e 's/Struthio_camelus_australis/Struthio_camelus/g' birds_total_GA_CT_edit_sites_percent_with_loss_status_tree_names -i
sed -e 's/Erythrura_gouldiae/Chloebia_gouldiae/g' -e 's/Stachyris_ruficeps/Cyanoderma_ruficeps/g' -e 's/Aquila_chrysaetos_chrysaetos/Aquila_chrysaetos/g' -e 's/Buceros_rhinoceros_silvestris/Buceros_rhinoceros/g' -e 's/Limosa_lapponica_baueri/Limosa_lapponica/g' -e 's/Lonchura_striata_domestica/Lonchura_striata/g' -e 's/Malurus_cyaneus_samueli/Malurus_cyaneus/g' -e 's/Struthio_camelus_australis/Struthio_camelus/g' birds_total_GA_CT_edit_sites_percent_with_loss_status -i

#In mac
cd
scp ceglab25@172.30.1.131:/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered_tree.nwk .
scp ceglab25@172.30.1.131:/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/birds_total_GA_CT_edit_sites_percent_with_loss_status_filtered .         



