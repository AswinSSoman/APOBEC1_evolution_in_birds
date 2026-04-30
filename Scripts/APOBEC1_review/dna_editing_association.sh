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

