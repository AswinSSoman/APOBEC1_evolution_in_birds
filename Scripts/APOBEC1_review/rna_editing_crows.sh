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

#Using reditools 1 ()
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/editing

#Run REDitools 1 ()
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus
time for b in $(ls rna/SRR*.bam)
do
p=$(echo $b | cut -f2 -d "/" | cut -f1 -d "_")
echo ">"$b ":" $p
time python2.7 /media/aswin/programs/REDItools/NPscripts/REDItoolDnaRnav13.py -i $b -j dna/dna_merged_sorted.bam -o editing/"$p"_editing -f genome/GCF_000738735.6_ASM73873v6_genomic.fna -t32 -c1,1 -m30,255 -v1 -q30,30 -v1 -e -n0.0 -N0.0 -u -l -p -s2 -g2 -S &> editing/"$p"_run_std.out
done

>rna/SRR1928171_Aligned.sortedByCoord.out.bam : SRR1928171:	1307m59.402s
>rna/SRR1947394_Aligned.sortedByCoord.out.bam : SRR1947394: 1314m33.407s

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Filtering:

# 
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus
time for i in $(cat rna_ids)
do
p=$(find editing/ -name "outTable_*" | grep "$i")
echo ">"$p
#Exclude invariant positions as well as positions not supported by ≥10 WGS reads
time awk 'FS="\t" {if ($8!="-" && $10>=10 && $13=="-") print}' $p > $p"_filtered.out"
#selecting sites with at least five RNAseq reads and a single mismatch:
time python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py -i $p"_filtered.out" -c 5 -v 1 -f 0.0 -o $p"_filtered.sel1"
#selecting sites with ≥10 RNAseq reads, three mismatches and minimum editing frequency of 0.1:
time python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py -i $p"_filtered.out" -c 10 -v 3 -f 0.1 -o $p"_filtered.sel2"
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
p1=$(awk -F "\t" -v r="$r" '$1==r {print$1,$32,$36,$42,$45}' OFS="\t" /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus/ncbi_SRP022901.tsv | tr " " "_" | cut -f3- -d "_" | cut -f1 -d "-")
tail -n +2 "$p"_all_subs_count.out | sort -k 2nr | sed "s/^/$r $p1 /g" 
unset p r p1
done | sed '1i SRR_ID\tTissue\tSubstitution\tRead_count\tTotal_reads\t%_Read\tSite_count\tTotal_sites\t%_Site' | sed 's/[ ]\+/\t/g' > summary_substitutions.tsv

Rscript /media/aswin/gene_loss/APOBEC1/RNA_editing/plot_read_site_count.R summary_substitutions.tsv plot_read_site_count.pdf
Rscript /media/aswin/gene_loss/APOBEC1/RNA_editing/plot_tissue_count.R summary_substitutions.tsv plot_tissue_count.pdf






