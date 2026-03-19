################################################################################################################################################################################################################################################
#Comment: 
################################################################################################################################################################################################################################################


cd /media/aswin/gene_loss/APOBEC1/RNA_editing/rnabam

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get inputs

#Chicken: data is based on: Comprehensive sequencing of the genome and transcriptome of the Xishuangbanna game fowl: https://doi.org/10.1038/s41597-024-04014-4

#Download genome
datasets download genome accession GCF_000002315.6 --include gff3,rna,cds,protein,genome,seq-report
datasets download genome accession GCF_000002315.6 --include gff3,genome,seq-report
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reference
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reference
datasets download genome accession GCF_000002315.6 --include genome,gtf,seq-report --dehydrated --filename GCF_000002315.6.zip
unzip GCF_000002315.6.zip -d GCF_000002315.6
datasets rehydrate --directory GCF_000002315.6

#Download SRA
sra-meta -id  SRP459583 > sra-meta -id  SRP459583_metadata.tsv

#There are mainly 2 projects: one project has only ine sample & 2nd one has 24 samples
#Download one manually: 4m28.902s
time /home/neo/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/neo/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR108/033/SRR10852933/SRR10852933_1.fastq.gz .
time /home/neo/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/neo/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR108/033/SRR10852933/SRR10852933_2.fastq.gz .

#Get project number & search in ENA browser & select the column names & download table as tsv file, 
#Check columns, especially "fastq ftp" & "Library Strategy" & "library layout" & "Library Source" & "Experiment Title" 
grep -if <(grep RNA Rhea_americana_sra | awk '{print$1}') ENA_PRJNA1338461_table.tsv | awk -F "\t" '{print$NF}' | awk -F ";" '{print$2"\n"$3}' | awk '{print"time /home/neo/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/neo/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@"$0,"."}' > ascp_links.sh
chmod +x ascp_links.sh
time ./ascp_links.sh


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get metadata

#Compare ids of RNA & DNA & enusre the sample is same between them
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/rnabam
ls | grep "SRR.*bam$" | cut -f1 -d "." | cut -f1 -d "_" | sort -u > ids
time sra-meta -id ids > ids_metadata

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/dnabam
ls | grep "SRR.*bam$" | cut -f1 -d "." | cut -f1 -d "_" | sort -u > ids
time sra-meta -id ids > ids_metadata

grep -f <(awk '{print$6}' ids_metadata) ../rnabam/ids_metadata ids_metadata | column -t

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QC inputs

#RNA-seq bam files
#Remove duplicates112m3.896s
time java -jar /home/ceglab25/ajs/smart-phase/Phase/libraries/picard.jar MarkDuplicates I=SRR30595317.coord.sorted.bam O=SRR30595317.coord.sorted.dedup.bam M=marked_dup_metrics.txt
	

################################################################################################################################################################################################################################################
#Identify RNA editing


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Using reditools

#Using reditools2
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/rnabam
time $PYTHON /media/aswin/programs/reditools2.0/src/cineca/reditools.py -f SRR30595317.coord.sorted.dedup.bam -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GCF_000002315.6_GRCg6a_genomic.fna -o rna_table_reditools3.txt 

#Using reditools3 (in ceglab25: 2m36.615s) 

#Rerun reditools3 with default parameters
#Run Reditools3 on RNA bam (379m12.303s)
time for b in $(ls *.bam | egrep -v "SRR30595317|SRR28369619.haplotagged.bam|SRR28369629.MAPQ20_40.bam")
do
echo $b
time python3.10 -m reditools analyze $b -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GRCg6a_clean.fa -o "rna_table_reditools3_"$b".txt" -t 32 -V &> "rna_table_reditools3_"$b".stdout"
done

#Index the reditools RNA output
time python3.10 -m reditools index rna_table_reditools3_SRR28208762_Aligned.sortedByCoord.out.bam.txt -o rna_table_reditools3_SRR28208762_Aligned.sortedByCoord.out.bam.index.txt

#Reditools on DNA bam: 216m31.053s
time python3.10 -m reditools analyze SRR28002323.coord.sorted.bam -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GRCg6a_clean.fa -o "dna_table_reditools3_SRR28002323.txt" -t 32 -V &> "dna_table_reditools3_SRR28002323.stdout"
time python3.10 -m reditools analyze SRR28002323.coord.sorted.bam -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GRCg6a_clean.fa -o "dna_table_reditools3_SRR28002323.txt" -t 64 -s2	-q30,30 -C -N &> "dna_table_reditools3_SRR28002323.stdout"

#Annotate in reditools 3
time python3.10 -m reditools annotate rna_table_reditools3_SRR28208762_Aligned.sortedByCoord.out.bam.txt ../dnabam/dna_table_reditools3_SRR28002323.txt > SRR28208762_SRR28002323_rna_dna.tsv
time python3 -m reditools index SRR28208762_SRR28002323_rna_dna.tsv -o SRR28208762_SRR28002323_rna_dna_index.tsv
#Exclude invariant positions as well as positions not supported by ≥10 WGS reads:
#where $8!="-" requires only variant positions (from column 8 of the output table), $10>=10 selects sites covered by ≥10 WGS reads and $13=="-" considers only WGS homozygous positions.
time awk 'FS="\t" {if ($8!="-" && $10>=10 && $13=="-") print}' SRR28208762_SRR28002323_rna_dna.tsv > SRR28208762_SRR28002323_rna_dna_filtered.out

#Commonly used filters:
#MeanQ = 30
#gMeanQ = 30
#coverage 1

python ../corso_epitrascrittomica/data_reditools/src/REDItools/main/REDItoolDnaRna.py -o /home/student_7/RNAseq -i SRR1319672.bam -f /data/annotations/GRCh37.primary_assembly.genome.fa -t 4 -c 0,1 -m 0,255 -v 1 -q 0,30 -e -n 0.0 -N 0.0 -u -l -p

#Reditools 1 takes a lot of time : 1026m36.015s
#Filters on : 
#-c : Min. read coverage (dna,rna)
#-m : Min. mapping quality score (dna,rna)
#-v : Min. num. of reads supporting the variation [3] for RNA-Seq
#-q : Min. quality score (dna,rna)
#-e : Exclude multi hits in RNA-Seq 
#-n : Min. editing frequency [0.1] for RNA-Seq
#-N : Min. variation frequency [0.1] for DNA-Seq
#-u : Consider mapping quality in RNA-Seq
#-l : Remove substitutions in homopolymeric regions in RNA-Seq
#-p : Use paired concardant reads only in RNA-Seq
#-s : Infer strand (for strand oriented reads)
#-g : Strand inference type 1:maxValue 2:useConfidence
#-S : Strand correction

#-d : Exclude duplicates in RNA-Seq
#-U : Consider mapping quality in DNA-Seq
#-G : Infer strand by GFF annotation (must be GFF and sorted, otherwise use -X)


#Identification of all DNA–RNA variants using Reditools 1
time python2.7 /media/aswin/programs/REDItools/main/REDItoolDnaRna.py -i SRR28208762_Aligned.sortedByCoord.out.bam -j ../dnabam/SRR28002323.coord.sorted.bam -o dna_corrected_rna_editing_SRR28002323_SRR28208762 -f /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GRCg6a_clean.fa -t32 -c1,1 -m30,255 -v1 -q30,30 -e -n0.0 -N0.0 -u -l -p -s2 -g2 -S

#Exclude invariant positions as well as positions not supported by ≥10 WGS reads:
#where $8!="-" requires only variant positions (from column 8 of the output table), $10>=10 selects sites covered by ≥10 WGS reads and $13=="-" considers only WGS homozygous positions.
time awk 'FS="\t" {if ($8!="-" && $10>=10 && $13=="-") print}' dna_corrected_rna_editing_SRR28002323_SRR28208762/DnaRna_297244823/outTable_297244823 > dna_corrected_rna_editing_SRR28002323_SRR28208762/DnaRna_297244823/outTable_297244823_filtered.out

#Create a first set of positions selecting sites supported by at least five RNAseq reads and a single mismatch:
time python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py -i dna_corrected_rna_editing_SRR28002323_SRR28208762/DnaRna_297244823/outTable_297244823_filtered.out -c 5 -v 1 -f 0.0 -o dna_corrected_rna_editing_SRR28002323_SRR28208762/DnaRna_297244823/outTable_297244823_filtered.out.sel1

#Create a second set of positions selecting sites supported by ≥10 RNAseq reads, three mismatches and minimum editing frequency of 0.1:
time python2.7 /media/aswin/programs/REDItools/accessory/selectPositions.py -i dna_corrected_rna_editing_SRR28002323_SRR28208762/DnaRna_297244823/outTable_297244823_filtered.out -c 10 -v 3 -f 0.1 -o dna_corrected_rna_editing_SRR28002323_SRR28208762/DnaRna_297244823/outTable_297244823_filtered.out.sel2



################################################################################################################################################################################################################################################
#Using jacusa1

java -jar /media/aswin/programs/JACUSA2/target/JACUSA2-2.1.15.jar

################################################################################################################################################################################################################################################

#Example runs
reditools.py -f rna.bam -r genome.fa -o rna_table.txt
reditools.py -f dna.bam -r genome.fa -B rna_table.txt -o dna_table.txt
annotate_with_DNA.py -r rna_table.txt -d dna_table.txt

python REDItoolDnaRna.py -i WT1.bam -j KO1.bam -f mm10.fa -c10,10 -s2 -g1 -v2 -V -d -D -o WT1vsKO1


#DRAFT
time $PYTHON /media/aswin/programs/reditools2.0/src/cineca/reditools.py \
-f SRR30595317.coord.sorted.dedup.bam \
-r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GCF_000002315.6_GRCg6a_genomic.fna \
-o rna_table.txt \
-q 30 \
-bq 30 \
-l 10


time python src/cineca/reditools.py \
-f ../rnabam/SRR30595317.coord.sorted.dedup.bam \
-r reference.fa \
-o rna_table.txt \
-q 30 \
-bq 30 \
-l 10

python src/cineca/reditools.py \
-f ../dnabam/merged.coord.sorted.bam \
-r reference.fa \
-B rna_table.txt \
-o dna_table.txt


python src/cineca/annotate_with_DNA.py \
-r rna_table.txt \
-d dna_table.txt

################################################################################################################################################################################################################################################

#INSTALL
conda install -c bioconda -c conda-forge pysam
pip install REDItools3



################################################################################################################################################################################################################################################
#DRAFT SCRIPTS
################################################################################################################################################################################################################################################

#reditools3 ran very fast but produced empty output file with following parameters
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/rnabam
time python3.10 -m reditools analyze SRR30595317.coord.sorted.dedup.bam -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GCF_000002315.6_GRCg6a_genomic.fna -o rna_table_reditools3.txt -q 255 -t 30 -s 2 -C -V

#Gives error
time python3.10 -m reditools analyze SRR30595317.coord.sorted.dedup.bam -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GCF_000002315.6_GRCg6a_genomic.fna -o rna_table_reditools3_run2.txt -t 30 -V

#empty output
time python3.10 -m reditools analyze SRR30595317.coord.sorted.dedup.bam -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GRCg6a_clean.fa -o rna_table_reditools3_run_q255.txt -q 255 -t 30 -V

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/rnabam
# 233m19.541s
time python3.10 -m reditools analyze SRR30595317.coord.sorted.dedup.bam -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GRCg6a_clean.fa -o rna_table_reditools3_run2.txt -t 30 -V

#nonempty output (230m0.046s)
time python3.10 -m reditools analyze SRR30595317.coord.sorted.dedup.bam -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GRCg6a_clean.fa -o rna_table_reditools3_run_q255.txt -t 30 -V -s 2 -C

#Using reditools3 (in Z840)
#cd ~/aswin/RNA_editing/rnabam
#time python3.10 -m reditools analyze SRR30595317.coord.sorted.dedup.bam -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GCF_000002315.6_GRCg6a_genomic.fna -o rna_table_reditools3.txt -q 255 -t 30 -s 2 -C -V

