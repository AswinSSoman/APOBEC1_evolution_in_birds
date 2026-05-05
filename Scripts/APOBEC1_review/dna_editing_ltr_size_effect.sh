
################################################################################################################################################################################################################################################################################################################

mkdir /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/effect_of_abundance
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/effect_of_abundance
awk '{print$1,$3,$4,$10}' OFS="\t" /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed_with_bird_orders > birds_ga_ltr.out


for r in $(ls ../../Taeniopygia_guttata/knisbacher/Data/Taegut/LTR/db/Len_*.txt | grep -v "Len_.txt")
do
awk '!($1="")' OFS="\n" $r | awk '{a+=$1;} END{print a;}'
done | awk '{a+=$1;} END{print a;}'


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Mammals

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/mammalian_vertebrates

#LTR/ERV size
for i in $(cat all_mammalian_vertebrates_genomes)
do
o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
s=$(find "$i"/knisbacher/Data/"$o"/LTR/db/ -name "Len_*.txt" ! -name "Len_.txt" -exec awk '{for(i=2;i<=NF;i++) sum+=$i} END {print sum}' {} +)
ga=$(grep $i /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/mammalian_vertebrates/results/summary/glimpse_output_of_knisbacher_script_step_6_Calculate_total_LTR_GA_edit_sites_14_mammalian_vertebrates | awk '{print$2}')
echo $i $s $ga
unset o s ga
done > mammals_ga_ltr.out

cp mammals_ga_ltr.out /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/effect_of_abundance/

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Non-mammalian vertebrates

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/non_mammalian_vertebrates
#LTR/ERV size
for i in $(cat all_nonmammalian_vertebrates_genomes_filtered)
do
o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
s=$(find "$i"/knisbacher/Data/"$o"/LTR/db/ -name "Len_*.txt" ! -name "Len_.txt" -exec awk '{for(i=2;i<=NF;i++) sum+=$i} END {print sum}' {} +)
ga=$(grep $i /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/non_mammalian_vertebrates/results/summary/glimpse_output_of_knisbacher_script_step_6_Calculate_total_LTR_GA_edit_sites_10_nonmammalian_vertebrates | awk '{print$2}')
echo $i $s $ga
unset o s ga
done > non_mammalian_vertebrates_ga_ltr.out

cp non_mammalian_vertebrates_ga_ltr.out /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/effect_of_abundance/

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Invertebrates

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/invertebrates
#LTR/ERV size
for i in $(cat all_invertebrate_genomes)
do
o=$(echo $i | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
s=$(find "$i"/knisbacher/Data/"$o"/LTR/db/ -name "Len_*.txt" ! -name "Len_.txt" -exec awk '{for(i=2;i<=NF;i++) sum+=$i} END {print sum}' {} +)
ga=$(grep $i /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/invertebrates/results/summary/glimpse_output_of_knisbacher_script_step_6_Calculate_total_LTR_GA_edit_sites_16_invertebrates | awk '{print$2}')
echo $i $s $ga
unset o s ga
done > invertebrate_genomes_ga_ltr.out

cp invertebrate_genomes_ga_ltr.out /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/effect_of_abundance/

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Compile data

for s in birds_ga_ltr.out invertebrate_genomes_ga_ltr.out mammals_ga_ltr.out non_mammalian_vertebrates_ga_ltr.out
do
n=$(echo $s | sed 's/_ga_ltr.out//g' | sed 's/_genomes//g')
if [[ "$s" == "birds_ga_ltr.out" ]]
then
cat $s | awk '{print$1,$3,$4}' | sed "s/^/$n /g"
else
cat $s | awk '{print$1,$2,$3}' | sed "s/^/$n /g"
fi
unset n
done | sed 's/[ \t]\+/\t/g' > all_ga_ltr.out

awk '$3>0' all_ga_ltr.out > all_ga_ltr_filtered.out

################################################################################################################################################################################################################################################################################################################
#Plot & statistical analysis

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/effect_of_abundance
python3 analyze_ltr_ga_stat.py

