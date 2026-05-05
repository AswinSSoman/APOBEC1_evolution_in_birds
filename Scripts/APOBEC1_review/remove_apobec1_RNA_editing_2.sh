
##Check how output of RNA editing sites varies with different filters in heart tissue

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947395_editing/DnaRna_402150276
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947395_editing/DnaRna_402150276
time python2.7 /media/aswin/programs/REDItools/accessory/AnnotateTable.py -a $rmsk -n rmsk -u -i outTable_402150276_filtered.out -o outTable_402150276_filtered_rmsk.out
time python2.7 /media/aswin/programs/REDItools/accessory/AnnotateTable.py -a $splice -n splice -u -i outTable_402150276_filtered_rmsk.out -o outTable_402150276_filtered_rmsk_splice.out


#Remove unknown as well
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947395_editing/DnaRna_402150276/filtered
awk -F "\t" '$15!~"Simple_repeat" && $15!~"Low_complexity" && $15!~"Unknown"' outTable_402150276_filtered_rmsk_splice.out | awk -F "\t" '$17=="-"' > filtered/outTable_402150276_filtered_rmsk_splice_filtered.out

#NOTE: don't use sup parameter & flag -u

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947395_editing/DnaRna_402150276/filtered
infile="outTable_402150276_filtered_rmsk_splice_filtered.out"
p=$(echo "$infile" | sed 's/\.out//g')

# parameter sets
c_vals=(10 15)
C_vals=(10 15)
f_vals=(0.0 0.1)
F_vals=(0.95)

# flag combinations
flags_list=("" "-e" "-r" "-e -r")

# header file (only once)
echo -e "SRR_ID\tTissue\tcov\tfreq\tflags\tSubstitution\tRead_count\tTotal_reads\tPct_Read\tSite_count\tTotal_sites\tPct_Site" > all_parameters_summary.tsv

for c in "${c_vals[@]}"; do
for C in "${C_vals[@]}"; do
for f in "${f_vals[@]}"; do
for F in "${F_vals[@]}"; do
for flags in "${flags_list[@]}"; do

    # clean flag name for filename
    flag_tag=$(echo "$flags" | sed 's/ /_/g; s/^-//; s/-/_/g')
    [ -z "$flag_tag" ] && flag_tag="none"

    outprefix="${p}_cov${c}.${C}_freq${f}${F}_${flag_tag}"

    echo "Running: c=$c C=$C f=$f F=$F flags=$flags"

    # run REDItools selection
    python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
        -i "$infile" -c "$c" -C "$C" -f "$f" -F "$F" $flags \
        -o "${outprefix}.out"

    # read-based counts
    python2.7 /media/aswin/programs/REDItools/accessory/subCount.py "${outprefix}.out" \
    | sort -t $'\t' -k1,1 \
    | awk 'BEGIN{print "Substitution\tRead_count\tTotal_reads\tPercentage"}1' \
    > "${outprefix}_readcount.out"

    # site-based counts
    python2.7 /media/aswin/programs/REDItools/accessory/subCount2.py "${outprefix}.out" \
    | sort -t $'\t' -k1,1 \
    | awk 'BEGIN{print "Substitution\tSite_count\tTotal_sites\tPercentage"}1' \
    > "${outprefix}_sitecount.out"

    # join results (skip headers before joining)
    tail -n +2 "${outprefix}_readcount.out" > tmp1
    tail -n +2 "${outprefix}_sitecount.out" > tmp2

    join -1 1 -2 1 tmp1 tmp2 > "${outprefix}_joined.out"

    # add metadata columns
    awk -v c="$c" -v f="$f" -v flags="$flag_tag" '
    BEGIN{OFS="\t"}
    {
        print "SRR1947394","dna",c,f,flags,$0
    }' "${outprefix}_joined.out" >> all_parameters_summary.tsv

done
done
done
done
done

# cleanup
rm -f tmp1 tmp2

echo -e "SRR_ID\tTissue\tcov\tfreq\tflags\tSubstitution\tRead_count\tTotal_reads\tPct_Read\tSite_count\tTotal_sites\tPct_Site" > all_parameters_summary_narrow.tsv
find . -name "${p}_cov*_freq*_*joined.out"  -print0 | while IFS= read -r -d '' file
do
    base=$(basename "$file")
    srr=$(readlink -f $file | cut -f10 -d "/" | cut -f1 -d "_")
    cov=$(echo "$base"  | grep -o 'cov[0-9\.]\+' | sed 's/cov//' )
    freq=$(echo "$base" | grep -o 'freq[0-9\.]\+' | sed 's/freq//')
    flags="none"
    [[ "$base" == *_e_* ]] && flags="e"
    [[ "$base" == *_r_* ]] && flags="r"
    [[ "$base" == *_e_r_* ]] && flags="e,r"
    # set NA if empty
    cov=${cov:-NA}
    freq=${freq:-NA}
    # ---- metadata lookup (Tissue etc.) ----
    tissue=$(awk -F "\t" -v r="$srr" '$1==r {print$32,$36,$42,$45}' OFS="\t" /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/ncbi_SRP022901.tsv | tr " " "_" | tr "\t" "\n" | grep -v "^-$" | sort -u)
    # ---- process substitution table ----
    tail -n +2 "$file" | while read -r sub rc tr pr sc ts ps
    do
        echo -e "${srr}\t${tissue}\t${cov}\t${freq}\t${flags}\t${sub}\t${rc}\t${tr}\t${pr}\t${sc}\t${ts}\t${ps}"
    done
    unset base srr cov sup freq flags tissue
done >> all_parameters_summary_narrow.tsv

Rscript /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643/unknown/gemini3.R all_parameters_summary_narrow.tsv all_parameters_summary_narrow.pdf















