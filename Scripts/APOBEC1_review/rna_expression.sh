##########################################################################################################################################################################################################################################################################################################
#Comment 6: RNA-seq 
##########################################################################################################################################################################################################################################################################################################


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#To get RNA seq info from NCBI database:

	#Go to genome browser : https://www.ncbi.nlm.nih.gov/datasets/genome/
	#Search species name (without underscores)
	#Get most updated genome with annotation
	#Go to it's genome data viewer
	#Search syntenic genes in Genome data viewer
	#Zoom out such that A1 & A1-like genes are visible & syntenic genes with robust RNA expression for reference.
	#Change view: show protein fetures, 
	#Add local track data: A1 local blastn, TOGA expected exon hits, A1-like local blastn 
	#Add exon coverage data from RNA-seq track & check the count: half the total sample i.e. 74 items means 37 samples
	#Tracks -> configure tracks -> Expression -> RNA-Seq , samples (74 items) -> check all exon coverafe tracks -> Configure

##########################################################################################################################################################################################################################################################################################################
#Download RNA seq data

##########################################################################################################################################################################################################################################################################################################
#Casuarius casuarius

#Download SRA: unknown sex , unknown tissue type: 1.5Gb
cd /media/ashutosh/disk3/RNA_seq
time /media/ashutosh/disk3/RNA_seq/sratoolkit.3.3.0-ubuntu64/bin/prefetch  --max-size 100000000 SRR10852845

#Download new genome since updated genome has annotation but lacks any genes named APOBEC1, AICDA, NANOG etc. which is strange!
datasets download genome accession GCA_013396415.1 --include genome,gtf,seq-report --dehydrated --filename GCA_013396415.1.zip
unzip GCA_013396415.1.zip -d GCA_013396415.1
datasets rehydrate --directory GCA_013396415.1

##########################################################################################################################################################################################################################################################################################################
#Leptosomus discolor

#Download SRA : unknown sex , blood: 1.1Gb
cd /media/ashutosh/disk3/RNA_seq
time /media/ashutosh/disk3/RNA_seq/sratoolkit.3.3.0-ubuntu64/bin/prefetch  --max-size 100000000 SRR10853056


##########################################################################################################################################################################################################################################################################################################
#Rhea_americana

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Download using aspera 

#Fetch SRA metadata first to get filter samples & extract their ftp links to download using aspera
cd ~/soft_links/Rhea_americana/genome
time sra-meta -s Rhea_americana > Rhea_americana_sra

#There are mainly 2 projects: one project has only ine sample & 2nd one has 24 samples
#Download one manually: 4m28.902s
time /home/neo/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/neo/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR108/033/SRR10852933/SRR10852933_1.fastq.gz .
time /home/neo/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/neo/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR108/033/SRR10852933/SRR10852933_2.fastq.gz .

#Get project number & search in ENA browser & select the column names & download table as tsv file, 
#Check columns, especially "fastq ftp" & "Library Strategy" & "library layout" & "Library Source" & "Experiment Title" 
grep -if <(grep RNA Rhea_americana_sra | awk '{print$1}') ENA_PRJNA1338461_table.tsv | awk -F "\t" '{print$NF}' | awk -F ";" '{print$2"\n"$3}' | awk '{print"time /home/neo/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/neo/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@"$0,"."}' > ascp_links.sh
chmod +x ascp_links.sh
time ./ascp_links.sh

#Uncompress data : 15m35.114s
time for s in $(ls SRR*.gz)
do
sn=$(echo $s | sed 's/\.gz//g' | sed 's/fastq/fa/g')
echo $sn
zcat $s | sed -n '1~4s/^@/>/p;2~4p' | sed '/>/ s!/!_!g' > $sn
unset sn
done

#In Z840
mkdir /media/ashutosh/disk3/RNA_seq/Rhea_americana
#In neo
time scp SRR*.fa workstation@172.30.1.172:/media/ashutosh/disk3/RNA_seq/Rhea_americana/
time scp SRR*.gz workstation@172.30.1.172:/media/ashutosh/disk3/RNA_seq/Rhea_americana
#In Z840
cd /media/ashutosh/disk3/RNA_seq/Rhea_americana
#14m13.799s
time cat SRR*.fa > combined_SRA_Rhea_americana.fa

#In neo
scp -r /home/neo/soft_links/Rhea_americana/aswin/APOBEC1/2nd_gblast/synteny/ workstation@172.30.1.172:/media/ashutosh/disk3/RNA_seq/Rhea_americana/bed_files

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#RNA mapping:

samtools faidx GCA_003343005.1_rheAme1_genomic.fna

#Index the genome from mapping
time /home/workstation/aswin/programs/STAR/source/STAR --runThreadN 64 --runMode genomeGenerate --genomeDir . --genomeFastaFiles GCA_003343005.1_rheAme1_genomic.fna --limitGenomeGenerateRAM 31000000000

read_1="SRR10852933_1.fastq.gz"
read_2="SRR10852933_2.fastq.gz"
id="SRR10852933"
ulimit -n 65535
time /home/workstation/aswin/programs/STAR/source/STAR --runThreadN 16 --genomeDir . --readFilesCommand zcat --readFilesIn $read_1 $read_2 --outFileNamePrefix "${id}_aligned_" --outSAMtype BAM SortedByCoordinate
samtools index SRR10852933_aligned_Aligned.sortedByCoord.out.bam

#STAR Mapping
cd /media/ashutosh/disk3/RNA_seq/Rhea_americana
ls SRR*.fastq.gz | cut -f1 -d "_" | sort -u > ids

time for id in $(cat ids)
do
read1=$(ls SRR*.fastq.gz | grep "$id"_1.fastq.gz)
read2=$(ls SRR*.fastq.gz | grep "$id"_2.fastq.gz)
echo ">" $id $read1 $read2
time /home/workstation/aswin/programs/STAR/source/STAR --runThreadN 16 --genomeDir . --readFilesCommand zcat --readFilesIn $read1 $read2 --outFileNamePrefix "${id}_" --outSAMtype BAM SortedByCoordinate
samtools index "$id"_Aligned.sortedByCoord.out.bam
unset read1 read2
done


