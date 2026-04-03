################################################################################################################################################################################################################################################
#Comment 10: Line 273, the title is too strong. It is not shown here unequivocally that Apobec1 does not edit RNA in general.
################################################################################################################################################################################################################################################

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download data for Flycatcher

#Based on : https://doi.org/10.1093/gbe/evt114

mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis

#Get metadata from NCBI & ENA browsers
mv PRJEB2984.csv ncbi_ERP001377.csv
csvformat -T ncbi_ERP001377.csv | tr " " "_" > ncbi_ERP001377.tsv

#All female samples are transcriptmic not genomic, hence use male
awk '{print$6,$17,$18,$24,$25,$27,$28,$33,$34,$36}' ncbi_ERP001377.tsv | sort -k7,7 -k10,10 | grep female | colnum.sh
#Male had samples from DNA & RNA, but unclear if they are from same individual.
#Individual information is not present in metadata, multiple sample can be from same individual.
awk '{print$1,$4/1000/1000/1000,$7/1024/1024/1024,$6,$17,$18,$24,$25,$27,$28,$33,$34,$36}' ncbi_ERP001377.tsv | awk '{$2+=0; $3+=0; print}' CONVFMT="%.2f" | sort -k10,10 -k12,12 | grep male -w | colnum.sh

#Choose SRA Runs to download
grep SAMEA1570901 ena_ERP001377.tsv | tr " " "_" | grep TRANSCRIPTOMIC | awk '{print$1}' > rna_ids
grep SAMEA1570901 ena_ERP001377.tsv | tr " " "_" | grep GENOMIC | awk '{print$1}' > dna_ids

#Create aspera links
grep SAMEA1570901 ena_ERP001377.tsv | tr " " "_" | grep TRANSCRIPTOMIC | awk '{print$14}' | tr ";" "\n"  | awk '{print"time /home/ceglab25/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab25/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@"$0,"."}' > rna_ascp_links.sh
grep SAMEA1570901 ena_ERP001377.tsv | tr " " "_" | grep GENOMIC | awk '{print$14}' | tr ";" "\n"  | awk '{print"time /home/ceglab25/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab25/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@"$0,"."}' > dna_ascp_links.sh

#Download data using aspera
chmod +x dna_ascp_links.sh rna_ascp_links.sh
#RNA (87m30.136s)
time ./rna_ascp_links.sh
#DNA (244m27.078s)
time ./dna_ascp_links.sh

#Uncompress (22m50.376s)
time for i in $(ls ERR*.fastq.gz); do echo ">"$i; time gzip -d $i; done

#Download genome
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/genome
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/genome
datasets download genome accession GCF_000247815.1 --include genome,gtf,seq-report --dehydrated
unzip ncbi_dataset.zip -d GCF_000247815.1 
time datasets rehydrate --directory GCF_000247815.1
mv GCF_000247815.1/ncbi_dataset/data/GCF_000247815.1/GCF_000247815.1_FicAlb1.5_genomic.fna .
mv GCF_000247815.1/ncbi_dataset/data/GCF_000247815.1/genomic.gtf GCF_000247815.1_FicAlb1.5_genomic.gtf
samtools faidx GCF_000247815.1_FicAlb1.5_genomic.fna

#One SRA data seems to be corrupted, hence remove it
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis
sed '/ERR117157/d' dna_ids -i

#Get repeatmasker data: from https://hgdownload.soe.ucsc.edu/hubs/birds/index.html
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/repeatmasker
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/repeatmasker
wget https://hgdownload.soe.ucsc.edu/hubs/GCF/000/247/815/GCF_000247815.1/GCF_000247815.1.repeatMasker.out.gz
wget https://hgdownload.soe.ucsc.edu/hubs/GCF/000/247/815/GCF_000247815.1/GCF_000247815.1.repeatMasker.version.txt
gzip -d GCF_000247815.1.repeatMasker.out.gz

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC

#Fastqc (77m2.297s)
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis
mkdir fastqc
time for i in $(cat rna_ids dna_ids); do echo ">"$i; fastqc "$i"_1.fastq "$i"_2.fastq; done
mv ERR*fastqc.zip *_fastqc.html fastqc/
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/fastqc
multiqc .

#fastp (+ 21m42.893s)
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis
#DNA ()
mkdir fastp
time for i in $(cat rna_ids)
do
echo ">" $i
time fastp -i "$i"_1.fastq -I "$i"_2.fastq -o fastp/out_"$i"_1.fastq -O fastp/out_"$i"_2.fastq -q 25 -u 10 -l 50 -y -x -w 32 -h fastp/fastp_"$i".html -j fastp/fastp_"$i".json 
done
#DNA (37m50.847s)
time for i in $(cat dna_ids)
do
echo ">" $i
time fastp -i "$i"_1.fastq -I "$i"_2.fastq -o fastp/out_"$i"_1.fastq -O fastp/out_"$i"_2.fastq -w 32 -h fastp/fastp_"$i".html -j fastp/fastp_"$i".json 
done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Alignment

#RNA-seq
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/rna
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/rna
#Indexing (7m42.146s)
time /media/aswin/programs/STAR_2.7.11b/Linux_x86_64_static/STAR --runThreadN 32 --runMode genomeGenerate --genomeDir . --genomeFastaFiles ../genome/GCF_000247815.1_FicAlb1.5_genomic.fna --limitGenomeGenerateRAM 31000000000
#Mapping using trimmed fastq (5m4.789s)
ulimit -n 65535
time for i in $(cat ../rna_ids ../dna_ids)
do
echo ">" $i
time /media/aswin/programs/STAR_2.7.11b/Linux_x86_64_static/STAR --runThreadN 16 --genomeDir . --readFilesIn ../fastp/out_"$i"_1.fastq ../fastp/out_"$i"_2.fastq --outFileNamePrefix "$i"_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outFilterMultimapNmax 1
samtools index "$i"_Aligned.sortedByCoord.out.bam
done

#WGS
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/dna
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/dna
#Indexing (7m8.397s)
time bwa index -a is ../genome/GCF_000247815.1_FicAlb1.5_genomic.fna -p GCF_000247815.1
#Mapping using trimmed fastq (67m19.119s)
time for i in $(cat ../dna_ids)
do
echo ">"$i
time bwa mem -t 32 GCF_000247815.1 -Y ../fastp/out_"$i"_1.fastq ../fastp/out_"$i"_2.fastq > "$i".sam
done

#NOTE: Got error wth ERR117158: Segmentation fault (core dumped)
rm ERR117158.sam
#Delete unncessary input files after use
rm /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/fastp/*.fastq

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prepare inputs

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/dna
#Sort bam files
start_time=$(date +%s)
for f in *.sam
do
(
    filename="${f%%.*}"
    echo $filename
    samtools view -bS "$f" -o "${filename}.bam"
) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time)) && echo "- " "$j" " : " $elapsed_time "secs"
unset start_time end_time elapsed_time

#Merge bam files ()
time samtools merge -@ 24 -o dna_merged.bam *.bam

#Sort ()
time samtools sort dna_merged.bam -@ 30 -m 3G -o dna_merged_sorted.bam
#Index ()
samtools index dna_merged_sorted.bam

#Compress & keep raw data (512m16.535s)
#for f in SRR*.fastq; do gzip "$f"; done
#grep -f dna_ids <(find fastp) | grep fastq | xargs rm
grep -f rna_ids <(ls SRR*) | xargs rm





#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#DRAFT SCRIPTS

#3 files got error during mapping, hence redownload raw data & rerun bwa mapping
ERR117158|ERR117157|ERR117153
egrep "ERR117158|ERR117153"
time fasterq-dump --split-files ERR117158
time fasterq-dump --split-files ERR117157
time fasterq-dump --split-files ERR117153

#ERR117157 consistently gave error in downloaded fastq, the raw data itself have some quality issues
time for i in $(cat dna_ids | egrep "ERR117158|ERR117153"); do echo ">"$i; fastqc "$i"_1.fastq "$i"_2.fastq; done

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis
sed '/ERR117157/d' dna_ids -i


