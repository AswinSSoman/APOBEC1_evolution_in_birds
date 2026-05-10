################################################################################################################################################################################################################################################
#Comment 10: Line 273, the title is too strong. It is not shown here unequivocally that Apobec1 does not edit RNA in general.
################################################################################################################################################################################################################################################

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download data for Corvus

#Based on : https://doi.org/10.1038/ncomms13195


mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus

#Get metadata from NCBI & ENA browsers
csvformat -T ncbi_SRP022901.csv | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' > ncbi_SRP022901.tsv
csvformat -T ncbi_SRP342045.csv | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' > ncbi_SRP342045.tsv

#Choose SRA Runs to download
awk '{print$1,$2,$4/1000/1000/1000,$6,$7/1024/1024/1024,$15,$18,$19,$22,$23,$27,$33,$32,$42,$45}' ncbi_SRP022901.tsv | sort -k8,8 -k10,10 | egrep -v "Corvus_brachyrhynchos|Corvus_frugilegus|Coloeus_monedula" \
 | sort -k8,8 -k10,10 -k9,9 | grep S_Up_H32 | awk '$2=="WGS"' | awk '{print$1}' > dna_ids
awk '{print$1,$2,$4/1000/1000/1000,$6,$7/1024/1024/1024,$15,$18,$19,$22,$23,$27,$33,$32,$42,$45}' ncbi_SRP022901.tsv | sort -k8,8 -k10,10 | egrep -v "Corvus_brachyrhynchos|Corvus_frugilegus|Coloeus_monedula" \
	| sort -k8,8 -k10,10 -k9,9 | grep S_Up_H32 | awk '$2=="RNA-Seq"' | awk '{print$1}' > rna_ids

#Create aspera links
grep -f dna_ids ena_SRP022901.tsv | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | tr " " "_" | awk '{print$20}' | tr ";" "\n"  | awk '{print"time /home/ceglab25/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab25/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@"$0,"."}' > dna_ascp_links.sh
grep -f rna_ids ena_SRP022901.tsv | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | tr " " "_" | awk '{print$20}' | tr ";" "\n"  | awk '{print"time /home/ceglab25/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab25/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@"$0,"."}' > rna_ascp_links.sh

#Download data using aspera
chmod +x dna_ascp_links.sh rna_ascp_links.sh
#RNA (165m32.231s)
time ./rna_ascp_links.sh
#DNA (624m16.845s)
time ./dna_ascp_links.sh

#Uncompress (129m10.748s)
time for i in $(ls SRR*.fastq.gz); do echo ">"$i; time gzip -d $i; done

#Download genome
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/genome
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/genome
datasets download genome accession GCF_000738735.6 --include genome,gtf,seq-report --dehydrated
unzip ncbi_dataset.zip -d GCF_000738735.6
time datasets rehydrate --directory GCF_000738735.6
mv GCF_000738735.6/ncbi_dataset/data/GCF_000738735.6/GCF_000738735.6_ASM73873v6_genomic.fna .
mv GCF_000738735.6/ncbi_dataset/data/GCF_000738735.6/genomic.gtf GCF_000738735.6_ASM73873v6_genomic.gtf
samtools faidx GCF_000738735.6_ASM73873v6_genomic.fna

#Get repeatmasker data: from https://hgdownload.soe.ucsc.edu/hubs/birds/index.html
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/repeatmasker
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/repeatmasker
rsync -a -P rsync://hgdownload.soe.ucsc.edu/hubs/GCF/000/738/735/GCF_000738735.6/GCF_000738735.6.repeatMasker.out.gz ./
wget https://hgdownload.soe.ucsc.edu/hubs/GCF/000/738/735/GCF_000738735.6/GCF_000738735.6.repeatMasker.version.txt
gzip -d GCF_000738735.6.repeatMasker.out.gz

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC

#Fastqc (158m41.452s)
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus
mkdir fastqc
time for i in $(cat rna_ids dna_ids); do echo ">"$i; fastqc "$i"_1.fastq "$i"_2.fastq; done
mv SRR*fastqc.zip *_fastqc.html fastqc/
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/fastqc
multiqc .

#fastp
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus
#DNA (31m36.734s)
mkdir fastp
time for i in $(cat rna_ids)
do
echo ">" $i
time fastp -i "$i"_1.fastq -I "$i"_2.fastq -o fastp/out_"$i"_1.fastq -O fastp/out_"$i"_2.fastq -q 25 -u 10 -l 50 -y -x -w 32 -h fastp/fastp_"$i".html -j fastp/fastp_"$i".json
done
#DNA (87m5.424s)
time for i in $(cat dna_ids)
do
echo ">" $i
time fastp -i "$i"_1.fastq -I "$i"_2.fastq -o fastp/out_"$i"_1.fastq -O fastp/out_"$i"_2.fastq -w 32 -h fastp/fastp_"$i".html -j fastp/fastp_"$i".json
done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Alignment

#RNA-seq
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/rna
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/rna
#Indexing (7m42.895s)
time /media/aswin/programs/STAR_2.7.11b/Linux_x86_64_static/STAR --runThreadN 32 --runMode genomeGenerate --genomeDir . --genomeFastaFiles ../genome/GCF_000738735.6_ASM73873v6_genomic.fna --limitGenomeGenerateRAM 31000000000
#Mapping using trimmed fastq (69m4.141s)
ulimit -n 65535
time for i in $(cat ../rna_ids)
do
echo ">" $i
time /media/aswin/programs/STAR_2.7.11b/Linux_x86_64_static/STAR --runThreadN 16 --genomeDir . --readFilesIn ../fastp/out_"$i"_1.fastq ../fastp/out_"$i"_2.fastq --outFileNamePrefix "$i"_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outFilterMultimapNmax 1
samtools index "$i"_Aligned.sortedByCoord.out.bam
done

#WGS
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/dna
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/dna
#Indexing ()
time bwa index -a is ../genome/GCF_000738735.6_ASM73873v6_genomic.fna -p GCF_000738735.6
#Mapping using trimmed fastq (380m39.028s)
time for i in $(cat ../dna_ids)
do
echo ">"$i
time bwa mem -t 32 GCF_000738735.6 -Y ../fastp/out_"$i"_1.fastq ../fastp/out_"$i"_2.fastq > "$i".sam
done &> bwa_run.stdout

#Delete unncessary input files after use
time rm /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/fastp/*.fastq
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prepare inputs

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/dna
#Sort bam files ()
start_time=$(date +%s)
for f in $(ls *.sam | sort -V | sed -n '21,30p')
do
(
    filename="${f%%.*}"
    echo $filename
#    samtools view -bS "$f" -o "${filename}.bam"
    samtools view -bS "$f" | samtools sort -@ 6 -m 1G -o "${filename}.sorted.bam"
) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time)) && echo "- " "$j" " : " $elapsed_time "secs"
unset start_time end_time elapsed_time

#Merge bam files (111m44.404s)
time samtools merge -@ 24 dna_merged.bam *.bam

#Sort (239m10.274s)
time samtools sort dna_merged.bam -@ 30 -m 2G -o dna_merged_sorted.bam
#Index (47m56.343s)
samtools index dna_merged_sorted.bam

#Compress & keep raw data (512m16.535s)
#for f in SRR*.fastq; do gzip "$f"; done
#grep -f dna_ids <(find fastp) | grep fastq | xargs rm
grep -f rna_ids <(ls SRR*) | xargs rm
rm SRR*.sam SRR*.bam

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Identify editing sites

#Using reditools 1
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing

#Run REDitools 1 (25244m55.598s)
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus
time for b in $(ls rna/SRR*.bam)
do
p=$(echo $b | cut -f2 -d "/" | cut -f1 -d "_")
echo ">"$b ":" $p
time python2.7 /media/aswin/programs/REDItools/NPscripts/REDItoolDnaRnav13.py -i $b -j dna/dna_merged_sorted.bam -o editing/"$p"_editing -f genome/GCF_000738735.6_ASM73873v6_genomic.fna -t32 -c1,1 -m30,255 -v1 -q30,30 -v1 -e -n0.0 -N0.0 -u -l -p -s2 -g2 -S &> editing/"$p"_run_std.out
done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Filter & plot

/media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/filter_classify_rna_edits_py36_v6.py

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing
ls -d */ | grep -v "old_plots/" | grep -v "plot" | tr -d "/" > rna_folders

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing

while read id
do
  echo ">" $id
  cd $id/DnaRna_*
  in=$(find . -name "outTable_*" | grep -v filtered)
  ac=$(echo $id | cut -f1 -d "_")
  bam=$(find /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/rna/ -name "*$ac*.bam")
python3 /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/filter_classify_rna_edits_py36_v6.py \
  --reditools $in \
  --repeatmasker /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/repeatmasker/GCF_000738735.6.repeatMasker.out \
  --genome /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/genome/GCF_000738735.6_ASM73873v6_genomic.fna \
  --rna-bam $bam \
  --dna-bam /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/dna/dna_merged_sorted.bam \
  --gtf /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/genome/GCF_000738735.6_ASM73873v6_genomic.gtf \
  --out-prefix results/sample \
  --write-failures \
  --use-dna-bam-counts-if-reditools-dna-missing
unset in ac bam
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing
done < <(cat rna_folders | egrep -v "SRR1947394_editing|SRR1947476_editing")

nohup bash -c 'time ./filter_and_plot.sh' &> filter_and_plot.stdout &

mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/final_plot

while read id
do
  echo ">" $id
  cd "$id"/DnaRna_*/results
  ac=$(echo $id | cut -f1 -d "_")
  p1=$(awk -F "\t" -v r="$r" '$1==r {print$32,$36,$42,$45}' OFS="\t" /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/ncbi_SRP022901.tsv | tr " " "_" | tr "\t" "\n" | grep -v "^-$" | sort -u)
  cp sample.event_class_counts.png /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/final_plot/"$ac"_"$p1"_event_class_counts.png
  cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing
  unset ac p1
done < rna_folders







#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Filtering:

#51m3.882s
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus
time for i in $(cat rna_ids)
do
p=$(find editing/ -name "outTable_*" | grep "$i" | grep -v "filtered")
echo ">"$p
#Exclude invariant positions as well as positions not supported by ≥10 WGS reads
time awk 'FS="\t" {if ($8!="-" && $10>=10 && $13=="-") print}' $p > $p"_filtered.out"
#selecting sites with at least five RNAseq reads and a single mismatch:
#time python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py -i $p"_filtered.out" -c 5 -v 1 -f 0.0 -o $p"_filtered.sel1"
#selecting sites with ≥10 RNAseq reads, three mismatches and minimum editing frequency of 0.1:
#time python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py -i $p"_filtered.out" -c 10 -v 3 -f 0.1 -o $p"_filtered.sel2"
unset p
done

#Count substitutions
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing
time for i in $(find . -name "*_filtered.out")
do
p=$(echo $i | sed 's/\.out//g')
python2.7 /media/aswin/programs/REDItools/accessory/subCount.py "$i" | sed '1i Substitution Read_count Total_reads Percentage' > "$p"_all_subs_readcount.out
python2.7 /media/aswin/programs/REDItools/accessory/subCount2.py "$i" | sed '1i Substitution Site_count Total_sites Percentage' > "$p"_all_subs_sitecount.out
join -1 1 -2 1 "$p"_all_subs_readcount.out "$p"_all_subs_sitecount.out | sed 's/[ ]\+/\t/g' > "$p"_all_subs_count.out
unset p
done

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing
for i in $(find . -name "*_filtered.out")
do
p=$(echo $i | sed 's/\.out//g')
r=$(echo $i | cut -f2 -d "/" | cut -f1 -d "_")
p1=$(awk -F "\t" -v r="$r" '$1==r {print$32,$36,$42,$45}' OFS="\t" /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/ncbi_SRP022901.tsv | tr " " "_" | tr "\t" "\n" | grep -v "^-$" | sort -u)
tail -n +2 "$p"_all_subs_count.out | sort -k 2nr | sed "s/^/$r $p1 /g"
unset p r p1
done | sed '1i SRR_ID\tTissue\tSubstitution\tRead_count\tTotal_reads\t%_Read\tSite_count\tTotal_sites\t%_Site' | sed 's/[ ]\+/\t/g' > summary_substitutions.tsv

#Plot
Rscript /media/aswin/gene_loss/APOBEC1/RNA_editing/plot_read_site_count.R summary_substitutions.tsv plot_read_site_count.pdf
Rscript /media/aswin/gene_loss/APOBEC1/RNA_editing/plot_tissue_count.R summary_substitutions.tsv plot_tissue_count.pdf

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Apply different filters & plot

#Apply 2 types of filters:

#nohup bash -c 'time ./filter_sets.sh' &> filter_sets.stdout &
#find editing/ -name "outTable_*_filtered*" -type f | grep -v "_filtered.out" | xargs rm

#Try different parameters (9m49.070s)
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus
time while read i
do
    p=$(find editing/ -name "outTable_*" | grep "$i" | grep -v "filtered")
    echo ">$p"
    infile="${p}_filtered.out"
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
    for sup in 5 10 15
    do
        python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
            -i "$infile" -v $sup -V $sup \
            -o "${p}_filtered_sup${sup}.out"
    done
    # ---- Combined conditions ----
    python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
        -i "$infile" -c 15 -C 15 -v 15 -V 15 -f 0.1 -F 0.95 -e -r \
        -o "${p}_filtered_cov15_sup15_freq01_e_r_both.out"
    python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
        -i "$infile" -c 15 -C 15 -v 15 -V 15 -f 0.1 -F 0.95 \
        -o "${p}_filtered_cov15_sup15_freq01_both.out"
    #python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
    #    -i "$infile" -c 15 -v 15 -f 0.1 \
    #    -o "${p}_filtered_cov15_sup15_freq01_rna.out"
    # ---- Flags ----
    python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
        -i "$infile" -e -r -u \
        -o "${p}_filtered_all_excl.out"
    for flag in e r u
    do
        python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py \
            -i "$infile" -$flag \
            -o "${p}_filtered_${flag}.out"
    done
done < <(grep -f <(ls -d editing/*/ | cut -f1 -d "_" | cut -f2 -d "/") rna_ids)

#Summary of parameters used
echo -e "sample\tcov\tsup\tfreq\tflags\tlines" > editing/parameter_summary.tsv
find editing/ -name "*_filtered*.out" | while read f
do
    # ---- sample name ----
    sample=$(basename "$f" | cut -d "_" -f2)
    # ---- defaults ----
    cov="NA"
    sup="NA"
    freq="NA"
    flags="none"
    name=$(basename "$f")
    # ---- extract coverage ----
    if [[ $name =~ cov([0-9]+) ]]; then
        cov="${BASH_REMATCH[1]}"
    fi
    # ---- extract support ----
    if [[ $name =~ sup([0-9]+) ]]; then
        sup="${BASH_REMATCH[1]}"
    fi
    # ---- extract frequency ----
    if [[ $name =~ freq([0-9]+) ]]; then
        freq="0.${BASH_REMATCH[1]}"
    fi
    if [[ $name =~ freq([0-9]+) ]]; then
    freq="${BASH_REMATCH[1]}"
    freq="${freq:0:-1}.${freq: -1}"
    fi
    # ---- extract flags ----
    if [[ $name =~ all_excl ]]; then
        flags="e,r,u"
    else
        tmp=""
        [[ $name =~ _e\.out ]] && tmp+="e,"
        [[ $name =~ _r\.out ]] && tmp+="r,"
        [[ $name =~ _u\.out ]] && tmp+="u,"
        [[ $name =~ _e_r_both\.out ]] && tmp+="e,r"
        flags=${tmp%,}
        [[ -z "$flags" ]] && flags="none"
    fi
    # ---- count lines (excluding header if needed) ----
    lines=$(wc -l < "$f")
    echo -e "$sample\t$cov\t$sup\t$freq\t$flags\t$lines"
done >> editing/parameter_summary.tsv

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing
mkdir plot_parameter_effects
Rscript plot_parameter_summary.R parameter_summary.tsv plot_parameter_effects/plot_parameter
Rscript plot_parameter_simple.R parameter_summary.tsv plot_parameter_effects/plot_parameter_simple.pdf

#Count substitutions (1m8.090s)
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing
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

#Final table
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing
echo -e "SRR_ID\tTissue\tcov\tsup\tfreq\tflags\tSubstitution\tRead_count\tTotal_reads\tPct_Read\tSite_count\tTotal_sites\tPct_Site" > all_parameters_summary_substitutions.tsv
find . -name "*_filtered*_all_subs_count.out" -print0 | \
while IFS= read -r -d '' file
do
    base=$(basename "$file")
    # ---- extract SRR/sample ID ----
    srr=$(echo $file | cut -f2 -d "/" | cut -f1 -d "_")
    # ---- extract parameters from filename ----
    cov=$(echo "$base"  | grep -o 'cov[0-9]\+' | sed 's/cov//' )
    sup=$(echo "$base"  | grep -o 'sup[0-9]\+' | sed 's/sup//' )
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

#Make sure counts are fine
awk '{print$3,$4,$5,$6}' all_parameters_summary_substitutions.tsv | sort | uniq -c

awk 'BEGIN{FS=OFS="\t"} NR==1 {print; next}
{if($6=="u") $6="only_positions_supported_by_DNA"
    else if($6=="e") $6="exclude_multiple_subs_in_RNA"
    else if($6=="r") $6="exclude_invariant_sites_in_RNA"
    else if($6=="e,r,u") $6="exclude_multiple_subs_in_RNA,exclude_invariant_sites_in_RNA,only_positions_supported_by_DNA"
    print}' all_parameters_summary_substitutions.tsv > all_parameters_summary_substitutions_edited.tsv

Rscript /media/aswin/gene_loss/APOBEC1/RNA_editing/plot_read_site_count.R summary_substitutions.tsv plot_read_site_count.pdf

Rscript plot_dual_axis.R all_parameters_summary_substitutions.tsv plot_SRR_dual.pdf 15 15 1 none
Rscript plot_dual_canonical.R all_parameters_summary_substitutions.tsv plot_dual_canonical.pdf 15 15 1 none
Rscript plot_dual_canonical_col.R all_parameters_summary_substitutions.tsv plot_dual_canonical_col.pdf 15 15 1 none
Rscript plot_dual_Y.R all_parameters_summary_substitutions.tsv plot_dual_Y.pdf 15 15 1 none

Rscript plot_gemini.R all_parameters_summary_substitutions.tsv plot_gemini.pdf 15 15 1 none
Rscript gemini3.R all_parameters_summary_substitutions.tsv gemini3.pdf
Rscript gemini4.R all_parameters_summary_substitutions.tsv gemini4.pdf
Rscript gemini5.R all_parameters_summary_substitutions.tsv gemini5.pdf
Rscript gemini6.R all_parameters_summary_substitutions.tsv gemini6.pdf
#7 is good
Rscript gemini7.R all_parameters_summary_substitutions.tsv gemini7.pdf
Rscript gemini8.R all_parameters_summary_substitutions.tsv gemini8.pdf

#10 is good
Rscript gemini10.R all_parameters_summary_substitutions.tsv gemini10.pdf
Rscript gemini11.R all_parameters_summary_substitutions.tsv gemini11.pdf

awk '!a[$3]++' all_parameters_summary_substitutions.tsv | awk '{print$3}' | paste -s -d " "
awk '!a[$4]++' all_parameters_summary_substitutions.tsv | awk '{print$4}' | paste -s -d " "
awk '!a[$5]++' all_parameters_summary_substitutions.tsv | awk '{print$5}' | paste -s -d " "
awk '!a[$6]++' all_parameters_summary_substitutions.tsv | awk '{print$6}' | paste -s -d " "
flags none u e,r,u r e


cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing
for i in $(find . -name "*_filtered_*_all_subs_count.out")
do
p=$(echo $i | sed 's/\.out//g')
r=$(echo $i | cut -f2 -d "/" | cut -f1 -d "_")
p1=$(awk -F "\t" -v r="$r" '$1==r {print$32,$36,$42,$45}' OFS="\t" /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/ncbi_SRP022901.tsv | tr " " "_" | tr "\t" "\n" | grep -v "^-$" | sort -u)
params=$()
tail -n +2 "$p".out | sort -k 2nr | sed "s/^/$r $p1 $params /g"
unset p r p1 folder
done | sed '1i SRR_ID\tTissue\tcov\tsup\tfreq\tflags\tSubstitution\tRead_count\tTotal_reads\t%_Read\tSite_count\tTotal_sites\t%_Site' | sed 's/[ ]\+/\t/g' > all_parameters_summary_substitutions.tsv

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#prepare splice sites annotations for REDItools

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
#gunzip gencode.v46.annotation.gtf.gz
#perl /media/aswin/programs/GMAP-GSNAP/util/gtf_splicesites.pl gencode.v46.annotation.gtf > ss
#awk -F" " '{split($2,a,":"); split(a[2],b,"."); if (b[1]>b[3]) print a[1],b[3],b[1],toupper(substr($3,1,1)),"-"; else print a[1],b[1],b[3],toupper(substr($3,1,1)),"+"}' ss > gencode.v46.ss.txt

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/genome
perl /media/aswin/programs/GMAP-GSNAP/util/gtf_splicesites.pl GCF_000738735.6_ASM73873v6_genomic.gtf > splice_sites
awk -F" " '{split($2,a,":"); split(a[2],b,"."); if (b[1]>b[3]) print a[1],b[3],b[1],toupper(substr($3,1,1)),"-"; else print a[1],b[1],b[3],toupper(substr($3,1,1)),"+"}' splice_sites \
	| awk '{print$1,"GCF_000738735.6",$4,$2-4,$3+4,".",$5,".","gene_id \"NA\"; transcript_id \"NA\";"}' OFS="\t" > reditools_splice_sites.txt
sort -k1,1 -k4,4n reditools_splice_sites.txt > reditools_splice_sites.sorted.txt
bgzip reditools_splice_sites.sorted.txt
tabix -p gff reditools_splice_sites.sorted.txt.gz
splice=$(readlink -f reditools_splice_sites.sorted.txt.gz)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#prepare RepeatMasker annotations for REDItools

#cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/repeatmasker
#awk 'OFS="\t"{print $6,"rmsk_mm10",$12,$7+1,$8,".",$10,".","gene_id \""$11"\"; transcript_id \""$13"\";"}' rmsk.txt > rmsk.gtf
#sort -k1,1 -k4,4n rmsk.gtf > rmsk.sorted.gtf
#bgzip rmsk.sorted.gtf
#tabix -p gff rmsk.sorted.gtf.gz

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/repeatmasker
sed '1,3d' GCF_000738735.6.repeatMasker.out | awk -v OFS="\t" '{print$5,"GCF_000738735.6",$11,$6+1,$7,".",$9,".","gene_id \""$10"\"; transcript_id \""$11"\";"}' | awk 'BEGIN{FS=OFS="\t"} $7=="C" {$7="-"} {print}' > GCF_000738735.6.repeatMasker.gtf
sort -k1,1 -k4,4n GCF_000738735.6.repeatMasker.gtf > GCF_000738735.6.repeatMasker.sorted.gtf
bgzip GCF_000738735.6.repeatMasker.sorted.gtf
tabix -p gff GCF_000738735.6.repeatMasker.sorted.gtf.gz
rmsk=$(readlink -f GCF_000738735.6.repeatMasker.sorted.gtf.gz)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Apply filters

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643
time python2.7 /media/aswin/programs/REDItools/accessory/AnnotateTable.py -a $rmsk -n rmsk -u -i outTable_154256643_filtered.out -o outTable_154256643_filtered_rmsk.out
time python2.7 /media/aswin/programs/REDItools/accessory/AnnotateTable.py -a $splice -n splice -u -i outTable_154256643_filtered_rmsk.out -o outTable_154256643_filtered_rmsk_splice.out

#Filter Low complexity, simple repeats & splice sites
awk -F "\t" '$15!~"Simple_repeat" && $15!~"Low_complexity"' outTable_154256643_filtered_rmsk_splice.out | awk -F "\t" '$17=="-"' > outTable_154256643_filtered_rmsk_splice_filtered.out

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/SRR1947394_editing/DnaRna_154256643/test
awk -iinplace 'BEGIN{FS=OFS="\t"} {$2="heart"; print}' all_parameters_summary_substitutions.tsv
awk -iinplace 'BEGIN{FS=OFS="\t"} {$1="SRR1947394"; print}' all_parameters_summary_substitutions.tsv
Rscript /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing/gemini2.R all_parameters_summary_substitutions.tsv gemini2.pdf







#Create a first set of positions selecting sites supported by at least five RNAseq reads and a single mismatch
python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py -i outTable_154256643_filtered_rmsk_splice.out -c 5 -v 1 -f 0.0 -o outTable_154256643_filtered_rmsk_splice.out.sel1
#Create a second set of positions selecting sites supported by ≥10 RNAseq reads, three mismatches and minimum editing frequency of 0.1
python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py -i outTable_154256643_filtered_rmsk_splice.out -c 10 -v 3 -f 0.1 -o outTable_154256643_filtered_rmsk_splice.out.sel2

#Select ALU sites from the first set of positions
awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)=="Alu" && $17=="-"&& $8!="-" && $13=="-") print}' parallel_table.txt_all_chr.out.rmsk.snp.sel1 > parallel_table.txt_all_chr.out.rmsk.snp.alu

#Select REP NON ALU sites from the second set of positions, excluding sites in Simple repeats or Low complexity regions
awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)!="Alu" && $15!="-" && $15!="Simple_repeat" && $15!="Low_complexity" && $17=="-" && $8!="-" && $9>=0.1) print}' parallel_table.txt_all_chr.out.rmsk.snp.sel2 > parallel_table.txt_all_chr.out.rmsk.snp.nonalu

#Select NON REP sites from the second set of positions
awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)!="Alu" && $15=="-" && $17=="-" && $8!="-" && $9>=0.1) print}' parallel_table.txt_all_chr.out.rmsk.snp.sel2 > parallel_table.txt_all_chr.out.rmsk.snp.nonrep


awk -F "\t" '{print$15}' outTable_154256643_filtered_rmsk_splice.out | sort | uniq -c

#Mainly
Simple_repeat
Low_complexity
DNA
LINE
LTR
SINE
Unknown
rRNA, snRNA, tRNA


while read -r id
do
(
    echo ">$id"
    cd "$id"/DnaRna_* || exit
    in=$(find . -name "outTable_*" | grep -v filtered | head -n1)
    ac=$(echo "$id" | cut -f1 -d "_")
    bam=$(find /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/rna/ \
          -name "*$ac*.bam" | head -n1)
    [[ -z "$bam" ]] && {
        echo "No BAM found for $id"
        exit
    }
    ) &
    done < <(grep -Ev "SRR1947394_editing|SRR1947476_editing" rna_folders
    )
    wait
