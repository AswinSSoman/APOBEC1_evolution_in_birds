################################################################################################################################################################################################################################################
#Comment 10: Line 273, the title is too strong. It is not shown here unequivocally that Apobec1 does not edit RNA in general.
################################################################################################################################################################################################################################################


#Get input data

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Chicken: data is based on: Comprehensive sequencing of the genome and transcriptome of the Xishuangbanna game fowl: https://doi.org/10.1038/s41597-024-04014-4
#Crows: The genomic landscape underlying phenotypic integrity in the face of gene flow in crows : https://doi.org/10.1126/science.1253226

#Get SRA data

#Get metadata from NCBI & ENA browsers
	#scp ENA_SRP459583_metadata.tsv ceglab25@172.28.65.125:/media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/

#For WXS data (RNA) (167m24.871s + 90m)
	grep "RNA\-Seq" ncbi_SRP459583_metadata.csv | sort -k41,41 -t $'\t' | awk -F "," '!a[$41]++ {print$1}' > rna_ids
	grep -f rna_ids ENA_SRP459583_metadata.tsv | awk -F "\t" '{print$21}' | tr ";" "\n" | awk '{print"time /home/ceglab25/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab25/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@"$0,"."}' > rna_ascp_links.sh
	chmod +x rna_ascp_links.sh
	time ./rna_ascp_links.sh

#Download WGS data (DNA) (94m57.447s)
	grep -if <(grep "RNA\-Seq" ENA_SRP459583_metadata.tsv | tr " " "_" | awk '{print$22}') <(cat ENA_SRP459583_metadata.tsv | tr " " "_" | awk '$13=="WGS" && $14=="GENOMIC"') | awk '{print$1}' > dna_ids
	grep -f dna_ids ENA_SRP459583_metadata.tsv | awk -F "\t" '{print$21}' | tr ";" "\n" | awk '{print"time /home/ceglab25/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab25/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@"$0,"."}' > dna_ascp_links.sh
	chmod +x dna_ascp_links.sh
	time ./dna_ascp_links.sh

#Uncompress (81m2.399s)
	time for i in $(ls SRR*.fastq.gz); do echo ">"$i; time gzip -d $i; done

#Download genome
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/genome
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/genome
time datasets download genome accession GCA_041920315.1 --include genome,gtf,seq-report --dehydrated
unzip ncbi_dataset.zip -d GCA_041920315.1 
mv GCA_041920315.1/ncbi_dataset/data/GCA_041920315.1/GCA_041920315.1_ASM4192031v1_genomic.fna .
samtools faidx GCA_041920315.1_ASM4192031v1_genomic.fna

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC

#Fastqc (133m36.377s)
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools
mkdir fastqc
 time for i in $(cat rna_ids dna_ids); do echo ">"$i; fastqc "$i"_1.fastq "$i"_2.fastq; done
mv SRR*fastqc.zip *_fastqc.html fastqc/

#fastp (94m56.508s)
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools
mkdir fastp
time for i in $(cat rna_ids dna_ids)
do
echo ">" $i
time fastp -i "$i"_1.fastq -I "$i"_2.fastq -o fastp/out_"$i"_1.fastq -O fastp/out_"$i"_2.fastq -q 25 -u 10 -l 50 -y -x -w 32 -h fastp/fastp_"$i".html -j fastp/fastp_"$i".json 
done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Alignment

mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/alignment
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/alignment

#STAR mapping (5m26.924s)
time /media/aswin/programs/STAR_2.7.11b/Linux_x86_64_static/STAR --runThreadN 32 --runMode genomeGenerate --genomeDir . --genomeFastaFiles ../genome/GCA_041920315.1_ASM4192031v1_genomic.fna  --limitGenomeGenerateRAM 31000000000
#134m49.188s
ulimit -n 65535
time for i in $(cat ../rna_ids | grep -v "SRR28369625")
do
echo ">" $i
time /media/aswin/programs/STAR_2.7.11b/Linux_x86_64_static/STAR --runThreadN 16 --genomeDir . --readFilesIn ../fastp/out_"$i"_1.fastq ../fastp/out_"$i"_2.fastq --outFileNamePrefix "$i"_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outFilterMultimapNmax 1
samtools index "$i"_Aligned.sortedByCoord.out.bam
done

#BWA mapping
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/dna
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/dna
time bwa index -a is ../genome/GCA_041920315.1_ASM4192031v1_genomic.fna -p GCA_041920315.1
#52m8.335s
time bwa mem -t 32 GCA_041920315.1 -Y ../SRR30595317_1.fastq ../SRR30595317_2.fastq > SRR30595317.sam
#
time bwa mem -t 32 GCA_041920315.1 -Y ../SRR28002323_1.fastq ../SRR28002323_2.fastq > SRR28002323.sam








