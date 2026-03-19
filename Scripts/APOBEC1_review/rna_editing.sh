################################################################################################################################################################################################################################################
#Comment 10: Line 273, the title is too strong. It is not shown here unequivocally that Apobec1 does not edit RNA in general.
################################################################################################################################################################################################################################################


#Get input data

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Chicken: data is based on: Comprehensive sequencing of the genome and transcriptome of the Xishuangbanna game fowl: https://doi.org/10.1038/s41597-024-04014-4

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

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC

#Fastqc (133m36.377s)
 time for i in $(cat rna_ids dna_ids); do echo ">"$i; fastqc "$i"_1.fastq "$i"_2.fastq; done

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools
mkdir fastp
time for i in $(cat rna_ids dna_ids)
do
echo ">" $i
time fastp -i "$i"_1.fastq -I "$i"_2.fastq -o fastp/out_"$i"_1.fastq -O fastp/out_"$i"_2.fastq -q 25 -u 10 -l 50 -y -x -w 32 -h fastp/fastp_"$i".html -j fastp/fastp_"$i".json 
done

time fastp -i SRR28369625_1.fastq -I SRR28369625_2.fastq -o out_SRR28369625_1.fastq -O out_SRR28369625_2.fastq -q 25 -u 10 -l 50 -y -x -w 32 -h 

















