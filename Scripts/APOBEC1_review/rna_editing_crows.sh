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
for f in $(ls *.sam | sort -V | head -n 20)
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

#Merge bam files ()
time samtools merge -@ 24 -o dna_merged.bam *.bam

#Sort ()
time samtools sort dna_merged.bam -@ 30 -m 3G dna_merged_sorted.bam
#Index ()
samtools index dna_merged_sorted.bam

#Compress & keep raw data (512m16.535s)
#for f in SRR*.fastq; do gzip "$f"; done
#grep -f dna_ids <(find fastp) | grep fastq | xargs rm
grep -f rna_ids <(ls SRR*) | xargs rm

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

















