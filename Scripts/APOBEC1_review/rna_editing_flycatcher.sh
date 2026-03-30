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

#Uncompress ()
time for i in $(ls ERR*.fastq.gz); do echo ">"$i; time gzip -d $i; done

#Download genome


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC

#Fastqc (133m36.377s)
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis
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


