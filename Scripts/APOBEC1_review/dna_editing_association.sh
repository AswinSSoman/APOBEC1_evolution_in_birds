################################################################################################################################################################################################################################################################################################################
# Plot divergence & time Vs DNA edit sites
################################################################################################################################################################################################################################################################################################################



################################################################################################################################################################################################################################################################################################################

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons

time while read species
do
echo ">"$species

#species specific folder
o=$(echo $species | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher
mkdir repeat_insertion_time_DNA_editing_association

#Set mismatch type
ga="GA"
ganum=$(echo "GA CT GC GT CA TA AG TC CG TG AC AT" | awk -v m="$ga" '{for(i=1;i<=NF;i++) if($i == m) print i}')
ct="CT"
ctnum=$(echo "GA CT GC GT CA TA AG TC CG TG AC AT" | awk -v m="$ct" '{for(i=1;i<=NF;i++) if($i == m) print i}')

#Get repeatmasker file
if [[ -d ../ucsc ]]; then rmsk=$(ls ../ucsc/*.repeatMasker.out); else rmsk=$(ls ../repeat_masker_using_"$species"_library/*.fna.out); fi

#Total ERV count from blastdb of knisbacher pipeline & manual selection of repeatmasker table: These 2 values should match
echo "Number of ERVs:"
awk '$11~"ERV"' $rmsk | awk '$NF!="*"' | wc -l
for r in $(ls Data/"$o"/LTR/db/Len_*.txt | grep -v "Len_.txt"); do r1=$(echo $r | awk -F "/" '{print$NF}')
awk '{print $1,NF-1}' $r | awk '{a+=$2} END{print a}' | sed "s/^/$r1 /g"
done |  awk '{a+=$2} END{print a}'

#Filter repeatmasker for selected ERVs
awk '$11~"ERV"' $rmsk | awk '$NF!="*"' > repeat_insertion_time_DNA_editing_association/erv_repeatmasker.out
awk '{print$2}' repeat_insertion_time_DNA_editing_association/erv_repeatmasker.out | statplot.R repeat_insertion_time_DNA_editing_association/divergence_distribution_erv_repeatmasker.out.png

#Get edited ERV's GA edit site count & divergence based on coordinates
time for bestpairs in $(find Data/"$o"/LTR/results -maxdepth 5 -mindepth 5 -path "*/Tracks/*/${ga}/pairwise_filter*" -name "bestPairsClusters_*.tab")
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
done > repeat_insertion_time_DNA_editing_association/edited_erv_per_div_GA_edit_count.out

#Get all ERV's GA edit site count & divergence based on coordinates
awk '{print$4,$5,$6,$7+1,$8,$9,$13}' repeat_insertion_time_DNA_editing_association/edited_erv_per_div_GA_edit_count.out | awk '{if($6=="-") $6="C"; print}' > repeat_insertion_time_DNA_editing_association/c1
awk '{print$11,$10,$5,$6,$7,$9,$2}' repeat_insertion_time_DNA_editing_association/erv_repeatmasker.out | sed 's|^LTR/||g' > repeat_insertion_time_DNA_editing_association/c2

#combine edit site & divergence
awk '
FNR==NR {
    k = $1 FS $2 FS $3 FS $4 FS $5 FS $6
    key[k] = key[k] ? key[k] SUBSEP $7 : $7
    next}
{k = $1 FS $2 FS $3 FS $4 FS $5 FS $6
    if(k in key){
        n = split(key[k], a, SUBSEP)
        for(i=1; i<=n; i++)
            print $0, a[i]}
    else{print $0, 0}}' repeat_insertion_time_DNA_editing_association/c1 repeat_insertion_time_DNA_editing_association/c2 \
| sed '1i family Subfamily chr start end strand percent_divergence GA_edit_count' | sed 's/[ \t]/\t/g' > repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.out

#Plot
#Rscript /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/plot_repeat_divergence_editing_count.R repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.out repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.png "$species"

unset o ga ct ganum ctnum rmsk bestpairs r r1
done < all_bird_genomes_used

################################################################################################################################################################################################################################################################################################################
#Collect plots in folder

mkdir /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/repeat_insertion_time_DNA_editing_association

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons

while read species
do
echo ">"$species
cp /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.out /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/repeat_insertion_time_DNA_editing_association/"$species"_all_erv_per_div_GA_edit_count.out
cp /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.png /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/repeat_insertion_time_DNA_editing_association/"$species"_all_erv_per_div_GA_edit_count.png
done < all_bird_genomes_used

################################################################################################################################################################################################################################################################################################################
#Add time in plot

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons

#Use human mutation rates:
average mutation rate:  2.5 x 10(-8)

while read species
do
echo ">"$species
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/repeat_insertion_time_DNA_editing_association/
#awk 'BEGIN{OFS="\t"; r=0.5e-9} {print $0, $7/(2*r)}' all_erv_per_div_GA_edit_count.out | awk 'BEGIN{OFS="\t"} NR==1{$9="Time"} {print}' > all_erv_per_div_GA_edit_count_human_time.out
#awk 'BEGIN{OFS="\t"; r=0.5e-9} {print $0, ($7/100)/(2*r)/1e6}' all_erv_per_div_GA_edit_count.out | awk 'BEGIN{OFS="\t"} NR==1{$9="Time"} {print}' > time2
#awk 'BEGIN{OFS="\t"; r=2.5e-8} {print $0, ($7/100)/(2*r)/1e6}' all_erv_per_div_GA_edit_count.out | awk 'BEGIN{OFS="\t"} NR==1{$9="Time"} {print}' > time3
#awk 'BEGIN{OFS="\t"; r=2.2e-9} {print $0, ($7/100)/(2*r)/1e6}' all_erv_per_div_GA_edit_count.out | awk 'BEGIN{OFS="\t"} NR==1{$9="Time"} {print}' > time4
awk 'BEGIN{OFS="\t"; r=0.5e-9} {print $0, ($7/100)/(r)/1e6}' all_erv_per_div_GA_edit_count.out | awk 'BEGIN{OFS="\t"} NR==1{$9="Time"} {print}' > time_all_erv_per_div_GA_edit_count.out
Rscript plot.R time_all_erv_per_div_GA_edit_count.out time_all_erv_per_div_GA_edit_count.png $species
done < all_bird_genomes_used

################################################################################################################################################################################################################################################################################################################

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/
while read species
do
  o=$(echo $species | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')
  cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher
  #find Data/"$o"/LTR/db/ -mindepth 1 -type d | grep -vw "files_" | awk -F "/" '{print$NF}' | sed 's/files_//g'
  ec=$(find Data/"$o"/LTR/db/ -mindepth 1 -type d | grep -vw "files_" | awk -F "/" '{print$NF}' | sed 's/files_//g' | grep ERV -c)
  nec=$(find Data/"$o"/LTR/db/ -mindepth 1 -type d | grep -vw "files_" | awk -F "/" '{print$NF}' | sed 's/files_//g' | grep ERV -vc)
  echo $species $ec $nec
  unset o ec nec
done < all_bird_genomes_used > /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/repeat_insertion_time_DNA_editing_association/number_of_subfams_per_species

for s in $(ls | grep "_all_erv_per_div_GA_edit_count.out")
do
  c1=$(wc -l < $s)
  c2=$(awk '$8>0' $s | grep -v percent_divergence | wc -l)
  c3=$(awk '$8>0' $s | grep -v percent_divergence | awk '{a+=$NF} END{print a}')
  if [[ -z $c3 ]]; then c3="0"; else :; fi
  n=$(echo $s | cut -f1,2 -d "_")
  echo $n $c1 $c3 $c2
  unset c1 c2 c3 n
done > all_species_edit_count

join -1 1 -2 1 <(sort -k1,1 number_of_subfams_per_species) <(sort -k1,1 all_species_edit_count) | sort -k2nr,2 -k3nr,3 | sed '1i Species Num_ERVs Num_Edits Num_uniq_ERV_fams Num_uniq_Non_ERV_fams' > number_edit_sites_ervs_unique_erv_families

################################################################################################################################################################################################################################################################################################################

#Python scripts
apobec_blast_consensus_report_v3="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/knisbacher_blast_filter_visualize/apobec_blast_consensus_report_v3/apobec_blast_consensus_report_v3.py"
apobec_blast_consensus_report_v3_source_analysis_plots="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/knisbacher_blast_filter_visualize/apobec_blast_consensus_report_v3_source_analysis_plots/apobec_blast_consensus_report_v3_source_analysis_plots.py"

#Set species
species="Gallus_gallus"
o=$(echo $species | awk -F "_" '{print substr($1,1,3)""substr($2,1,3)}')

#Create folder
mkdir /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/Data/"$o"/LTR/results/knisbacher_blast_filter_visualize

#Check edit sites from knisbacher pipeline
cat /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/Data/"$o"/LTR/results/Total_edit_sites/GA_TotalEditSites_*_LTR_fam_1e-0.txt | column -t
find /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/Data/"$o"/LTR/results/ -maxdepth 5 -mindepth 5 -path "*/Tracks/*/GA/pairwise_filter*" -name "bestPairsClusters_*.tab" | xargs -n1 sh -c 'cat $0' | awk '{print$1,$2,$3,$4,$5,$6,$8,$9,$10,$13,$15,$16,$17,$18}' | column -t

subfam="GGLTR10D"
consensus_seq=$(find /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/RepBase/Libraries/perSeq/ -name "$subfam.fa" -type f)
blast=$(find /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/Data/"$o"/LTR/results/unzipped_blasts/ -name "*$subfam*" -type f)
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/Data/"$o"/LTR/results/knisbacher_blast_filter_visualize

#Run simple report
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/Data/"$o"/LTR/results/knisbacher_blast_filter_visualize
mkdir DNA_edit_report_"$subfam"/apobec_blast_consensus_report_v3/
cd DNA_edit_report_"$subfam"/apobec_blast_consensus_report_v3/
time python $apobec_blast_consensus_report_v3 \
  --consensus $consensus_seq \
  --blast $blast \
  --out-html "$subfam"_APOBEC_report_v3.html \
  --out-prefix "$subfam"_v3 \
  --mismatch GA \
  --include-controls \
  --cluster-mode consecutive \
  --min-cluster-sites 5 \
  --min-clustered-sites-per-alignment 10 \
  --include-failed-alignment-panels \
  --max-alignment-panels 30 \
  --max-failed-alignment-panels 20

#Run report including source analysis
mkdir -p DNA_edit_report_"$subfam"/apobec_blast_consensus_report_v3_source_analysis_plots/
cd DNA_edit_report_"$subfam"/apobec_blast_consensus_report_v3_source_analysis_plots/
time python3 $apobec_blast_consensus_report_v3_source_analysis_plots \
  --consensus $consensus_seq \
  --blast $blast \
  --out-html "$subfam"_APOBEC_source_analysis_plots.html \
  --out-prefix "$subfam"_source_analysis_plots \
  --mismatch GA \
  --include-controls \
  --cluster-mode consecutive \
  --min-cluster-sites 5 \
  --min-clustered-sites-per-alignment 5 \
  --consensus-map-method kmer \
  --max-alignment-panels 40 \
  --include-failed-alignment-panels \
  --max-matrix-sites 80 \
  --max-matrix-copies 80 \
  --max-distance-copies 40

################################################################################################################################################################################################################################################################################################################
#plot heatmap

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/repeat_insertion_time_DNA_editing_association
python3 plot_ga_edits_vs_divergence_mya_overlay_v4_spacer_labels.py /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/repeat_insertion_time_DNA_editing_association \
  --gene-loss-file all_gene_loss_dates.out \
  --pattern "*_all_erv_per_div_GA_edit_count.out" \
  --bin-width 1 \
  --metric per_1000 \
  --min-elements 50 \
  --date-model jc69 \
  --substitution-rate-per-myr 0.0019 \
  --scale log1p \
  --output ga_edits_loss_events_mya_spacer_labels.png \
  --summary-csv ga_edits_by_divergence_summary.csv \
  --group-csv species_apobec1_group_assignments.csv \
  --event-csv apobec1_loss_event_ranges_projected.csv \
  --display-layout-csv heatmap_display_layout.csv

python plot_ga_edits_vs_divergence_mya_overlay_v4_spacer_labels_font18_spacing_x40_species_clear.py \
  /home/workstation/APOBEC1_evolution_in_birds/Supplementary_files/Gene_loss_dating/DNA_editing_association/repeat_insertion_time_DNA_editing_association \
  --gene-loss-file all_gene_loss_dates.out \
  --output ga_edits_loss_events_mya_font18_spacing.pdf

################################################################################################################################################################################################################################################################################################################
#Script folder & usage for html report of knisbacher DNA editing detetction of ERVs by blast alignments & their QC & analysis by Nagarjun sir

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/knisbacher_blast_filter_visualize
ls  | xargs -n1 sh -c 'echo ">Script folder"; tree -hD $0 | sed "s/^/   /g" ; cat $0/readme | sed "s/^/   /g"' | sed 's/.*directories.*/   Usage:/g'

>Script folder
   apobec_blast_consensus_report_v3
   ├── [ 84K May  2 10:20]  apobec_blast_consensus_report_v3.py
   ├── [1.0M May  2 10:20]  GGLTR10D_APOBEC_report_v3.html
   ├── [ 358 May  2 10:20]  GGLTR10D.fa
   └── [ 417 May  2 11:30]  readme

   Usage:
   python apobec_blast_consensus_report_v3.py \
     --consensus GGLTR10D.fa \
     --blast Seq_GGLTR10D \
     --out-html GGLTR10D_APOBEC_report_v3.html \
     --out-prefix GGLTR10D_v3 \
     --mismatch GA \
     --include-controls \
     --cluster-mode consecutive \
     --min-cluster-sites 5 \
     --min-clustered-sites-per-alignment 10 \
     --include-failed-alignment-panels \
     --max-alignment-panels 30 \
     --max-failed-alignment-panels 20
>Script folder
   apobec_blast_consensus_report_v3_source_analysis
   ├── [103K May  2 11:13]  apobec_blast_consensus_report_v3_source_analysis.py
   ├── [1.6M May  2 11:13]  GGLTR10D_APOBEC_source_analysis.html
   └── [ 450 May  2 11:29]  readme

   Usage:
   python3 -S apobec_blast_consensus_report_v3_source_analysis.py \
     --consensus GGLTR10D.fa \
     --blast Seq_GGLTR10D \
     --out-html GGLTR10D_APOBEC_source_analysis.html \
     --out-prefix GGLTR10D_source_analysis \
     --mismatch GA \
     --include-controls \
     --cluster-mode consecutive \
     --min-cluster-sites 5 \
     --min-clustered-sites-per-alignment 5 \
     --consensus-map-method kmer \
     --max-alignment-panels 40 \
     --include-failed-alignment-panels
>Script folder
   apobec_blast_consensus_report_v3_source_analysis_plots
   ├── [143K May  2 11:27]  apobec_blast_consensus_report_v3_source_analysis_plots.py
   ├── [1.9M May  2 11:27]  GGLTR10D_APOBEC_source_analysis_plots.html
   └── [ 547 May  2 11:28]  readme

   Usage:
   python3 apobec_blast_consensus_report_v3_source_analysis_plots.py \
     --consensus GGLTR10D.fa \
     --blast Seq_GGLTR10D \
     --out-html GGLTR10D_APOBEC_source_analysis_plots.html \
     --out-prefix GGLTR10D_source_analysis_plots \
     --mismatch GA \
     --include-controls \
     --cluster-mode consecutive \
     --min-cluster-sites 5 \
     --min-clustered-sites-per-alignment 5 \
     --consensus-map-method kmer \
     --max-alignment-panels 40 \
     --include-failed-alignment-panels \
     --max-matrix-sites 80 \
     --max-matrix-copies 80 \
     --max-distance-copies 40
