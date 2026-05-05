
#Check how output of RNA editing sites varies with different filters

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Apply filters

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643
time python2.7 /media/aswin/programs/REDItools/accessory/AnnotateTable.py -a $rmsk -n rmsk -u -i outTable_154256643_filtered.out -o outTable_154256643_filtered_rmsk.out
time python2.7 /media/aswin/programs/REDItools/accessory/AnnotateTable.py -a $splice -n splice -u -i outTable_154256643_filtered_rmsk.out -o outTable_154256643_filtered_rmsk_splice.out

#Filter Low complexity, simple repeats & splice sites

#Mainly 
Simple_repeat
Low_complexity
DNA
LINE
LTR
SINE
Unknown
rRNA, snRNA, tRNA

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Remove repeats & splice sites

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643

#Remove simple repeat & low complexity
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643/test
awk -F "\t" '$15!~"Simple_repeat" && $15!~"Low_complexity"' outTable_154256643_filtered_rmsk_splice.out | awk -F "\t" '$17=="-"' > test/outTable_154256643_filtered_rmsk_splice_filtered.out

#Remove unknown as well
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643/unknown
awk -F "\t" '$15!~"Simple_repeat" && $15!~"Low_complexity" && $15!~"Unknown"' outTable_154256643_filtered_rmsk_splice.out | awk -F "\t" '$17=="-"' > unknown/outTable_154256643_filtered_rmsk_splice_filtered_unknown.out

mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643/dna
awk -F "\t" '$15!~"Simple_repeat" && $15!~"Low_complexity" && $15!~"Unknown" && $15!~"DNA"' outTable_154256643_filtered_rmsk_splice.out | awk -F "\t" '$17=="-"' > dna/outTable_154256643_filtered_rmsk_splice_filtered_dna.out

mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643/line
awk -F "\t" '$15!~"Simple_repeat" && $15!~"Low_complexity" && $15!~"Unknown" && $15!~"LINE"' outTable_154256643_filtered_rmsk_splice.out | awk -F "\t" '$17=="-"' > line/outTable_154256643_filtered_rmsk_splice_filtered_line.out

mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643/sine
awk -F "\t" '$15!~"Simple_repeat" && $15!~"Low_complexity" && $15!~"Unknown" && $15!~"SINE"' outTable_154256643_filtered_rmsk_splice.out | awk -F "\t" '$17=="-"' > sine/outTable_154256643_filtered_rmsk_splice_filtered_sine.out

mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643/ltr
awk -F "\t" '$15!~"Simple_repeat" && $15!~"Low_complexity" && $15!~"Unknown" && $15!~"LTR"' outTable_154256643_filtered_rmsk_splice.out | awk -F "\t" '$17=="-"' > ltr/outTable_154256643_filtered_rmsk_splice_filtered_ltr.out

mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643/all_repeats
awk -F "\t" '$15!~"Simple_repeat" && $15!~"Low_complexity" && $15!~"Unknown" && $15!~"DNA" && $15!~"LINE" && $15!~"SINE" && $15!~"LTR"' outTable_154256643_filtered_rmsk_splice.out | awk -F "\t" '$17=="-"' > all_repeats/outTable_154256643_filtered_rmsk_splice_filtered_all_repeats.out


./check_different_filtering_to_remove_A1_editing.sh line
./check_different_filtering_to_remove_A1_editing.sh sine
./check_different_filtering_to_remove_A1_editing.sh ltr
./check_different_filtering_to_remove_A1_editing.sh all_repeats

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
back="/media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643/sine"

#Count substitutions (1m8.090s)
cd $back

time while read i
do
    p=$(find . -name "outTable_*" | grep "$i")
    echo ">$p"
    infile="${p}"
    # ---- Frequency variations ----
    for f in 0.0 0.1
    do
        python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
            -i "$infile" -f $f \
            -o "${p}_filtered_freq${f/./}.out"
    done
    # ---- Coverage variations ----
    for cov in 5 10 15
    do
        python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
            -i "$infile" -c $cov -C $cov \
            -o "${p}_filtered_cov${cov}.out"
    done
    # ---- Support variations ----
#    for sup in 5 10 15
#    do
#        python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
#            -i "$infile" -v $sup -V $sup \
#            -o "${p}_filtered_sup${sup}.out"
#    done
    # ---- Combined conditions ----
#    python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
#        -i "$infile" -c 15 -C 15 -v 15 -V 15 -f 0.1 -F 0.95 -e -r \
#        -o "${p}_filtered_cov15_sup15_freq01_e_r_both.out"
#    python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
#        -i "$infile" -c 15 -C 15 -v 15 -V 15 -f 0.1 -F 0.95 \
#        -o "${p}_filtered_cov15_sup15_freq01_both.out"
    #python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
    #    -i "$infile" -c 15 -v 15 -f 0.1 \
    #    -o "${p}_filtered_cov15_sup15_freq01_rna.out"
    # ---- Flags ----
#    python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
#        -i "$infile" -e -r -u \
#        -o "${p}_filtered_all_excl.out"
    for flag in e r
    do
        python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
            -i "$infile" -$flag \
            -o "${p}_filtered_${flag}.out"
    done
done < <(echo 154256643)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cd $back
time find . -name "*_filtered*.out" -print0 | while IFS= read -r -d '' file
do
    base="${file%.out}"
    echo "Processing: $file"
    # ---- Read-level substitution counts ----
    python2.7 /media/aswin/programs/REDItools/accessory/subCount.py "$file" | sort -k1,1 | awk 'BEGIN{print "Substitution\tRead_count\tTotal_reads\tPercentage"}1' > "${base}_all_subs_readcount.out"
    # ---- Site-level substitution counts ----
    python2.7 /media/aswin/programs/REDItools/accessory/subCount2.py "$file" | sort -k1,1 | awk 'BEGIN{print "Substitution\tSite_count\tTotal_sites\tPercentage"}1' > "${base}_all_subs_sitecount.out"
    # ---- Join both tables safely ----
    join -1 1 -2 1 "${base}_all_subs_readcount.out" "${base}_all_subs_sitecount.out" > "${base}_all_subs_count.out"
done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Final table
cd $back
echo -e "SRR_ID\tTissue\tcov\tsup\tfreq\tflags\tSubstitution\tRead_count\tTotal_reads\tPct_Read\tSite_count\tTotal_sites\tPct_Site" > all_parameters_summary_substitutions.tsv
find . -name "*_filtered*_all_subs_count.out" -print0 | \
while IFS= read -r -d '' file
do
    base=$(basename "$file")
    # ---- extract SRR/sample ID ----
    srr=$(readlink -f $file | cut -f10 -d "/" | cut -f1 -d "_")
    # ---- extract parameters from filename ----
    cov=$(echo "$base"  | grep -o 'cov[0-9]\+' | sed 's/cov//' )
#    sup=$(echo "$base"  | grep -o 'sup[0-9]\+' | sed 's/sup//' )
    freq=$(echo "$base" | grep -o 'freq[0-9\.]\+' | sed 's/freq//' | sed 's/0/0\./1')
    # flags
    flags="none"
    [[ "$base" == *all_excl* ]] && flags="e,r,u"
    [[ "$base" == *_e_* ]] && flags="e"
    [[ "$base" == *_r_* ]] && flags="r"
    [[ "$base" == *_u_* ]] && flags="u"
    [[ "$base" == *_e_r_* ]] && flags="e,r"
    # set NA if empty
    cov=${cov:-NA}
    sup=${sup:-NA}
    freq=${freq:-NA}
    # ---- metadata lookup (Tissue etc.) ----
    tissue=$(awk -F "\t" -v r="$srr" '$1==r {print$32,$36,$42,$45}' OFS="\t" /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/ncbi_SRP022901.tsv | tr " " "_" | tr "\t" "\n" | grep -v "^-$" | sort -u)
    # ---- process substitution table ----
    tail -n +2 "$file" | while read -r sub rc tr pr sc ts ps
    do
        echo -e "${srr}\t${tissue}\t${cov}\t${sup}\t${freq}\t${flags}\t${sub}\t${rc}\t${tr}\t${pr}\t${sc}\t${ts}\t${ps}"
    done
    unset base srr cov sup freq flags tissue
done >> all_parameters_summary_substitutions.tsv

Rscript /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/gemini2.R all_parameters_summary_substitutions.tsv all_parameters_summary_substitutions.pdf

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

awk -iinplace 'BEGIN{FS=OFS="\t"} {$2="heart"; print}' all_parameters_summary_substitutions.tsv
awk -iinplace 'BEGIN{FS=OFS="\t"} {$1="SRR1947394"; print}' all_parameters_summary_substitutions.tsv
Rscript /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/gemini2.R all_parameters_summary_substitutions.tsv SRR1947394_editing_without_selectPositions.pdf


python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py -i "$infile" -c $cov -C $cov -o "${p}_filtered_cov${cov}.out"


python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py -i "$infile" -c $cov1 -C $cov2 -f $freq1 -F $freq2 -e -r -o "${p}_filtered_c_${cov1}_C_${cov2}.out"

infile="outTable_154256643_filtered_rmsk_splice_filtered_unknown.out"
p=$(echo $infile | sed 's/\.out//g')
python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py -i "$infile" -c 5 -C 5 -f 0.0 -F 0.95 -e -r -o "${p}_cov5.5_freq0.00.95_e_r.out"
python2.7 /media/aswin/programs/REDItools/accessory/subCount.py "${p}_cov5.5_freq0.00.95_e_r.out" | sort -k1,1 | awk 'BEGIN{print "Substitution\tRead_count\tTotal_reads\tPercentage"}1' > "${p}_cov5.5_freq0.00.95_e_r_all_subs_readcount.out"
python2.7 /media/aswin/programs/REDItools/accessory/subCount2.py "${p}_cov5.5_freq0.00.95_e_r.out" | sort -k1,1 | awk 'BEGIN{print "Substitution\tSite_count\tTotal_sites\tPercentage"}1' > "${p}_cov5.5_freq0.00.95_e_r_all_subs_sitecount.out"
join -1 1 -2 1 "${p}_cov5.5_freq0.00.95_e_r_all_subs_readcount.out" "${p}_cov5.5_freq0.00.95_e_r_all_subs_sitecount.out" > "${p}_cov5.5_freq0.00.95_e_r_all_subs_count.out"
echo -e "SRR_ID\tTissue\tcov\tfreq\tflags\tSubstitution\tRead_count\tTotal_reads\tPct_Read\tSite_count\tTotal_sites\tPct_Site" > all_parameters_cov5.5_freq0.00.95_e_r.tsv
find . -name "${p}_cov5.5_freq0.00.95_e_r_all_subs_count.out" -print0 | while IFS= read -r -d '' file
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
done >> all_parameters_cov5.5_freq0.00.95_e_r.tsv
Rscript gemini3.r all_parameters_cov5.5_freq0.00.95_e_r.tsv all_parameters_cov5.5_freq0.00.95_e_r.pdf

#NOTE: Filtering different Repeat types did not improve output, I expected A1 related sites to reduce, it reduced but not much compared to ADA edit sites
#NOTE: don't use sup parameter & flag -u
#Use cov, freq & flag -e & -r

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Try more combinations with these filters: cov, freq, flags -e -r

infile="outTable_154256643_filtered_rmsk_splice_filtered_unknown.out"
p=$(echo "$infile" | sed 's/\.out//g')

# parameter sets
c_vals=(5 10 15)
C_vals=(5 10 15)
f_vals=(0.0 0.1)
F_vals=(0.1 0.95)

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

    # ---------------------------
    # Extract SRR ID (safer)
    # ---------------------------
    srr=$(readlink -f "$file" | awk -F'/' '{print $10}' | cut -d'_' -f1)

    # ---------------------------
    # Extract parameters
    # ---------------------------
    cov=$(echo "$base"  | grep -oP 'cov\K[0-9]+\.[0-9]+')
    freq=$(echo "$base" | grep -o 'freq[0-9\.]\+' | sed 's/freq//')

    # ---------------------------
    # Extract flags properly
    # ---------------------------
    flags="none"
    [[ "$base" == *_e_r_* ]] && flags="e,r"
    [[ "$base" == *_e_* && "$base" != *_e_r_* ]] && flags="e"
    [[ "$base" == *_r_* && "$base" != *_e_r_* ]] && flags="r"

    # fallback
    cov=${cov:-NA}
    freq=${freq:-NA}

    # ---------------------------
    # Metadata lookup (Tissue etc.)
    # ---------------------------
    tissue=$(awk -F "\t" -v r="$srr" '
        $1==r {print $32,$36,$42,$45}
    ' OFS="\t" /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/ncbi_SRP022901.tsv \
    | tr " " "_" | tr "\t" "\n" | grep -v "^-$" | sort -u | paste -sd "," -)

    # ---------------------------
    # Process substitution table
    # ---------------------------
    tail -n +2 "$file" | while IFS=$'\t' read -r sub rc tr pr sc ts ps
    do
        echo -e "${srr}\t${tissue}\t${cov}\t${freq}\t${flags}\t${sub}\t${rc}\t${tr}\t${pr}\t${sc}\t${ts}\t${ps}"
    done

done >> all_parameters_summary_narrow.tsv

Rscript gemini3.R all_parameters_summary_narrow.tsv all_parameters_summary_narrow.pdf


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

Rscript gemini3.R all_parameters_summary_narrow.tsv all_parameters_summary_narrow.pdf

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643
time python2.7 /media/aswin/programs/REDItools/accessory/AnnotateTable.py -a $rmsk -n rmsk -u -i outTable_154256643_filtered.out -o outTable_154256643_filtered_rmsk.out
time python2.7 /media/aswin/programs/REDItools/accessory/AnnotateTable.py -a $splice -n splice -u -i outTable_154256643_filtered_rmsk.out -o outTable_154256643_filtered_rmsk_splice.out



