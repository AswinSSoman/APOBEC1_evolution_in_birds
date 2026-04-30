

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/Geospiza_fortis/knisbacher

#TE size
egrep "Species|Geo" /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/QC/size_of_main_TE_per_species
egrep "Species|Geo" /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/QC/size_of_each_repeatmasker_identified_transposon_type_per_genome | column -t
egrep "Species|Geo" /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/QC/size_of_each_repeatmasker_identified_octa_processed_transposon_type_per_genome | column -t

#LTR size
egrep "Species|Geo" /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/QC/size_of_LTR_harvest_identified_LTRs_per_genome | column -t
egrep "Species|Geo" /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/QC/size_of_all_LTRs_identified | column -t
egrep "Species|Geo" /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/supplementary/extracted_TEs | column -t

#summary
egrep "Species|Geo" /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/summary_GA_edit_sites_ERV_size | column -t
colnum.sh /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/birds_total_GA_CT_edit_sites_count_with_loss_status_ltr_size_clusters_after_each_filterings_all_relax_normalized_renamed_with_bird_orders

#normalized log GA
cat /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/correlation_plot/GA_Log_Normalized_filtered | column -t

#Total edit dites
colnum.sh /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/results/knisbacher_method/glimpse_output_of_knisbacher_script_step_6_Calculate_total_GA_edit_sites_81_birds 

#repeatmasker
head /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/Geospiza_fortis/ucsc/GCF_000277835.1.repeatMasker.out

awk '{print $1,NF-1}'/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/Geospiza_fortis/knisbacher/Data/Geofor/LTR/db/Len_ERVL.txt

grep TguLTRL4d ../../../../../ucsc/GCF_000277835.1.repeatMasker.out | awk '$NF!="*"' | awk '{print$10,$7-$6+1}' | less

#Total ERVs used by knisbacher script
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/Geospiza_fortis/knisbacher/Data/Geofor/LTR/db
awk '$11~"ERV"' ../../../../../ucsc/GCF_000277835.1.repeatMasker.out | awk '$NF!="*"' | wc -l
awk '{print $1,NF-1}' Len_ERVL.txt | awk '{a+=$2} END{print a}'





################################################################################################################################################################################################################################################################################################################
#ERV divergence and edit site count table summary

#Geospiza_fortis
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/Geospiza_fortis/knisbacher/Data/Geofor/LTR/results

ga="GA"
ganum=$(echo "GA CT GC GT CA TA AG TC CG TG AC AT" | awk -v m="$ga" '{for(i=1;i<=NF;i++) if($i == m) print i}')
ct="CT"
ctnum=$(echo "GA CT GC GT CA TA AG TC CG TG AC AT" | awk -v m="$ct" '{for(i=1;i<=NF;i++) if($i == m) print i}')

#Get edit count
while read -r bp
do
s=$(echo $bp | awk '{print$1,$2,$3,$4,$5,$6}' | tr ":-" " " | tr -d "+-")
c=$(echo $bp | awk '{print$NF}')
i1=$(echo $c | cut -f $ganum -d "|")
i2=$(echo $c | tr "|" "\n" | sed "$ganum d" | sort -nr | head -1)
i3=$(calc $i1 - $i2)
i4=$(echo $c | cut -f $ctnum -d "|")
echo $s $i1 $i4 $i2 $i3
unset s c i1 i2 i3 i4
done < Tracks/tracks_Geofor_LTR_ERVL_1e-0_5/GA/pairwise_filter/bestPairsClusters_Geofor_LTR_ERVL_1e-0_5.tab > edited_GA_count

#Get repeat divergence
time while read i
do
subfam=$(echo $i | awk '{print$5}')
chr=$(echo $i | awk '{print$6}' | cut -f1 -d ":")
rstart=$(echo $i | awk '{print$6}' | cut -f2 -d ":" | cut -f1 -d "-" | awk '{print$1+1}')
cstart=$(echo $i | awk '{print$6}' | cut -f2 -d ":" | cut -f1 -d "-" | awk '{print$1}')
end=$(echo $i | awk '{print$6}' | cut -f2 -d ":" | cut -f2 -d "-" | tr -d "+")
#echo  $subfam $chr $start $end
p1=$(awk -v chr="$chr" -v start="$rstart" -v end="$end" -v subfam="$subfam" '$10~subfam && $5==chr && $6==start && $7==end' /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/Geospiza_fortis/ucsc/GCF_000277835.1.repeatMasker.out | awk '{print$1,$2,$3,$4}')
p2=$(awk -v chr="$chr" -v start="$cstart" -v end="$end" -v subfam="$subfam" '$5==subfam && $6==chr && $7==start && $8==end' edited_GA_count | awk '{print$1,$2,$3,$4,$5,$6,$7,$8,$8-$7-1,$9,$10,$11,$12}')
echo $p2 $p1
unset chr rstart cstart end subfam p1 p2
done < Tracks/tracks_Geofor_LTR_ERVL_1e-0_5/GA/pairwise_filter/bestPairsClusters_Geofor_LTR_ERVL_1e-0_5.tab \
 | sed '1i mismatch species Repeat_type family subfamily chromosome start end length GA_count CT_count 2nd_highest_count subtrated_GA_count SW_score Per_Div Per_Del Per_Ins' | sed 's/[ \t]\+/\t/g' > edited_repeatmasker_edit_count.out


#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Filter repeatmasker for selected ERVs

awk '$11~"ERV"' ../../../../../ucsc/GCF_000277835.1.repeatMasker.out | awk '$NF!="*"' > erv_repeatmasker.out
awk '{print$2}' erv_repeatmasker.out | statplot.R --biplot divergence_distribution_erv_repeatmasker.out.png

awk '$11~"ERV"' ../../../../../ucsc/GCF_000277835.1.repeatMasker.out | awk '$NF!="*"' | awk '{print$10,$5,$6-1,$7,$2}' > f1.out
awk '{print$5,$6,$7,$8,$NF}' edited_GA_count > f2.out

awk 'NR==FNR {key=$1$2$3$4; val[key]=$5; next} {key=$1$2$3$4; if (key in val) print $0, val[key]; else print $0, 0}' f2.out f1.out \
	| sed '1i Subfamily chr start end percent_divergence GA_edit_count' | sed 's/[ \t]\+/\t/g' > erv_div_ga_edit_count.out

Rscript plot_div_vs_edit.R erv_div_ga_edit_count.out erv_div_ga_edit_count.pdf

################################################################################################################################################################################################################################################################################################################

species=Gallus_gallus
o=$(echo $species | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher
mkdir repeat_insertion_time_DNA_editing_association

#Set mismatch type
ga="GA"
ganum=$(echo "GA CT GC GT CA TA AG TC CG TG AC AT" | awk -v m="$ga" '{for(i=1;i<=NF;i++) if($i == m) print i}')
ct="CT"
ctnum=$(echo "GA CT GC GT CA TA AG TC CG TG AC AT" | awk -v m="$ct" '{for(i=1;i<=NF;i++) if($i == m) print i}')

#Get repeatmasker file
if [[ -d ../ucsc ]]
then
rmsk=$(ls ../ucsc/*.repeatMasker.out)
else
rmsk=$(ls ../repeat_masker_using_"$species"_library/*.fna.out)
fi

#Total ERV count from blastdb of knisbacher pipeline & manual selection of repeatmasker table: These 2 values should match
awk '$11~"ERV"' $rmsk | awk '$NF!="*"' | wc -l
for r in $(ls Data/"$o"/LTR/db/Len_*.txt | grep -v "Len_.txt"); do r1=$(echo $r | awk -F "/" '{print$NF}')
awk '{print $1,NF-1}' $r | awk '{a+=$2} END{print a}' | sed "s/^/$r1 /g"
done |  awk '{a+=$2} END{print a}'

#Filter repeatmasker for selected ERVs
awk '$11~"ERV"' $rmsk | awk '$NF!="*"' > repeat_insertion_time_DNA_editing_association/erv_repeatmasker.out
awk '{print$2}' erv_repeatmasker.out | statplot.R --biplot divergence_distribution_erv_repeatmasker.out.png

#Get edited ERV's GA edit site count & divergence based on coordinates
for bestpairs in $(find Data/"$o"/LTR/results -maxdepth 5 -mindepth 5 -path "*/Tracks/*/${ga}/pairwise_filter*" -name "bestPairsClusters_*.tab")
do
while read -r bp
do
s=$(echo $bp | awk '{print$1,$2,$3,$4,$5,$6}' | awk '{split($NF,a,/[:\-+]/); s=substr($NF,length($NF),1); $(NF)=a[1]; print $0, a[2], a[3], s}')
c=$(echo $bp | awk '{print$NF}')
i1=$(echo $c | cut -f $ganum -d "|")
i2=$(echo $c | tr "|" "\n" | sed "$ganum d" | sort -nr | head -1)
i3=$(calc $i1 - $i2)
i4=$(echo $c | cut -f $ctnum -d "|")
subfam=$(echo $s | awk '{print$5}')
chr=$(echo $s | awk '{print$6}')
start=$(echo $s | awk '{print$7+1}')
end=$(echo $s | awk '{print$8}')
strand=$(echo $s | awk '{print$9}' | awk '{if($1=="+") print$1; else print"C"}')
p1=$(awk -v chr="$chr" -v start="$start" -v end="$end" -v subfam="$subfam" -v strand="$strand" '$10~subfam && $5==chr && $6==start && $7==end && $9==strand' $rmsk | awk '{print$1,$2,$3,$4}')
echo $s $i1 $i4 $i2 $i3 $p1
unset s c i1 i2 i3 i4 subfam chr start end p1 strand
done < $bestpairs
unset bp
done | awk '{print$1,$2,3,4,5,6,7}' > repeat_insertion_time_DNA_editing_association/edited_erv_per_div_GA_edit_count.out

#Get all ERV's GA edit site count & divergence based on coordinates
awk '{print$5,$6,$7+1,$8,$9,$13}' repeat_insertion_time_DNA_editing_association/edited_erv_per_div_GA_edit_count.out | awk '{if($5=="-") $5="C"; print}' > repeat_insertion_time_DNA_editing_association/c1
awk '{print$10,$5,$6,$7,$9,$2}' repeat_insertion_time_DNA_editing_association/erv_repeatmasker.out > repeat_insertion_time_DNA_editing_association/c2

awk 'FNR==NR {
  key[$1 FS $2 FS $3 FS $4 FS $5] = $6
  next}
{  k = $1 FS $2 FS $3 FS $4 FS $5
  print $0, (k in key ? key[k] : 0)
}' repeat_insertion_time_DNA_editing_association/c1 repeat_insertion_time_DNA_editing_association/c2 \
 | sed '1i Subfamily chr start end strand percent_divergence GA_edit_count' | sed 's/[ \t]/\t/g' > repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.out

#Plot
Rscript /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/Geospiza_fortis/knisbacher/Data/Geofor/LTR/results/plot_div_vs_edit.R repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.out repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.pdf

################################################################################################################################################################################################################################################################################################################
#DRAFT FOR DRAFT

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/Geospiza_fortis/knisbacher/Data/Geofor/LTR/results

#Filter repeatmasker for selected ERVs
while read i
do
chr=$(echo $i | cut -f1,2 -d "_" | tr -d ">")
start=$(echo $i | sed 's/\.\..*//g' | awk -F "_" '{print$NF}')
end=$(echo $i | sed 's/.*\.\.//g' | awk -F "_" '{print$1}')
subfam=$(echo $i | sed 's/.*\.\.//g' | awk -F "_" '{print$4}')
awk -v chr="$chr" -v start="$start" -v end="$end" -v subfam="$subfam" '$10~subfam && $5==chr && $6==start && $7==end' ../../../../../ucsc/GCF_000277835.1.repeatMasker.out
unset chr start end subfam
done < <(grep ">" ../../../../../transposons/ERVs.fa) > erv_repeatmasker.out
  
##Filter repeatmasker for selected ERVs (faster)
awk 'FNR==NR && /^>/{
  sub(/^>/,""); split($0,a,"_"); chr=a[1]"_"a[2]
  split($0,b,/\.\./)
  split(b[1],c,"_"); start=c[length(c)]
  split(b[2],d,"_"); key[chr"|"start"|"d[1]"|"d[4]]
  next}
{k=$5"|"$6"|"$7"|"$10
  if (k in key) print}' ../../../../../transposons/ERVs.fa ../../../../../ucsc/GCF_000277835.1.repeatMasker.out > erv_repeatmasker.out

