



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

awk 'FNR==NR {key[$1 FS $2 FS $3 FS $4 FS $5 FS $6] = $7
  next}
{k = $1 FS $2 FS $3 FS $4 FS $5 FS $6
  print $0, (k in key ? key[k] : 0)}' repeat_insertion_time_DNA_editing_association/c1 repeat_insertion_time_DNA_editing_association/c2 \
 | sed '1i family Subfamily chr start end strand percent_divergence GA_edit_count' | sed 's/[ \t]/\t/g' > repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.out

#Plot
Rscript /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/plot_repeat_divergence_editing_count.R repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.out repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.png "$species"

unset o ga ct ganum ctnum rmsk bestpairs r r1
done < all_bird_genomes_used

################################################################################################################################################################################################################################################################################################################

mkdir /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/repeat_insertion_time_DNA_editing_association

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons

while read species
do
echo ">"$species
cp /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.out /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/repeat_insertion_time_DNA_editing_association/"$species"_all_erv_per_div_GA_edit_count.out
cp /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/repeat_insertion_time_DNA_editing_association/all_erv_per_div_GA_edit_count.png /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/repeat_insertion_time_DNA_editing_association/"$species"_all_erv_per_div_GA_edit_count.png
done < all_bird_genomes_used

################################################################################################################################################################################################################################################################################################################

cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons

#Use human mutation rates:
average mutation rate:  2.5 x 10(-8)

while read species
do
echo ">"$species
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/"$species"/knisbacher/repeat_insertion_time_DNA_editing_association/
awk 'BEGIN{OFS="\t"; r=0.5e-9} {print $0, $7/(2*r)}' all_erv_per_div_GA_edit_count.out | awk 'BEGIN{OFS="\t"} NR==1{$9="Time"} {print}' > all_erv_per_div_GA_edit_count_human_time.out
awk 'BEGIN{OFS="\t"; r=0.5e-9} {print $0, ($7/100)/(2*r)/1e6}' all_erv_per_div_GA_edit_count.out | awk 'BEGIN{OFS="\t"} NR==1{$9="Time"} {print}' > time2
awk 'BEGIN{OFS="\t"; r=2.5e-8} {print $0, ($7/100)/(2*r)/1e6}' all_erv_per_div_GA_edit_count.out | awk 'BEGIN{OFS="\t"} NR==1{$9="Time"} {print}' > time3
awk 'BEGIN{OFS="\t"; r=2.2e-9} {print $0, ($7/100)/(2*r)/1e6}' all_erv_per_div_GA_edit_count.out | awk 'BEGIN{OFS="\t"} NR==1{$9="Time"} {print}' > time4
awk 'BEGIN{OFS="\t"; r=0.5e-9} {print $0, ($7/100)/(r)/1e6}' all_erv_per_div_GA_edit_count.out | awk 'BEGIN{OFS="\t"} NR==1{$9="Time"} {print}' > time5
Rscript plot.R time5 time5.png Gallus_gallus
done < all_bird_genomes_used


################################################################################################################################################################################################################################################################################################################

apobec_blast_consensus_report_v3_source_analysis_plots="/media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/knisbacher_blast_filter_visualize/apobec_blast_consensus_report_v3_source_analysis_plots/apobec_blast_consensus_report_v3_source_analysis_plots.py"

subfam="GGLTR7A"
consensus_seq=$(find /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/RepBase/Libraries/perSeq/ -name "$subfam.fa" -type f)
blast=$(find /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/Gallus_gallus/knisbacher/Data/Galgal/LTR/results/unzipped_blasts/ -name "*$subfam*" -type f)
mkdir /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/Gallus_gallus/knisbacher/Data/Galgal/LTR/results/knisbacher_blast_filter_visualize
cd /media/aswin/gene_loss/APOBEC1/hypermutation_analyses/identify_retrotransposons/Gallus_gallus/knisbacher/Data/Galgal/LTR/results/knisbacher_blast_filter_visualize

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
