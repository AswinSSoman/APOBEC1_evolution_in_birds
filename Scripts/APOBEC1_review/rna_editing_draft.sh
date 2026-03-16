################################################################################################################################################################################################################################################
#Comment: 
################################################################################################################################################################################################################################################


cd /media/aswin/gene_loss/APOBEC1/RNA_editing/rnabam

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get inputs

#Download genome
datasets download genome accession GCF_000002315.6 --include gff3,rna,cds,protein,genome,seq-report
datasets download genome accession GCF_000002315.6 --include gff3,genome,seq-report
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reference
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reference
datasets download genome accession GCF_000002315.6 --include genome,gtf,seq-report --dehydrated --filename GCF_000002315.6.zip
unzip GCF_000002315.6.zip -d GCF_000002315.6
datasets rehydrate --directory GCF_000002315.6

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
#Run Reditools3 on RNA bam
time for b in $(ls *.bam | egrep -v "SRR30595317|SRR28369619.haplotagged.bam|SRR28369629.MAPQ20_40.bam")
do
echo $b
time python3.10 -m reditools analyze $b -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GRCg6a_clean.fa -o "rna_table_reditools3_"$b".txt" -t 32 -V &> "rna_table_reditools3_"$b".stdout"
done

#Index the reditools RNA output
time python3.10 -m reditools index rna_table_reditools3_SRR28208762_Aligned.sortedByCoord.out.bam.txt -o rna_table_reditools3_SRR28208762_Aligned.sortedByCoord.out.bam.index.txt

#Reditools on DNA bam: 216m31.053s
time python3.10 -m reditools analyze SRR28002323.coord.sorted.bam -r /media/aswin/gene_loss/APOBEC1/RNA_editing/reference/GRCg6a_clean.fa -o "dna_table_reditools3_SRR28002323.txt" -t 32 -V &> "dna_table_reditools3_SRR28002323.stdout"

time python3.10 -m reditools annotate rna_table_reditools3_SRR28208762_Aligned.sortedByCoord.out.bam.txt ../dnabam/dna_table_reditools3_SRR28002323.txt > SRR28208762_SRR28002323_rna_dna.tsv

3) Detect all potential RNA variants in your input BAM using the REDItoolDnaRNA.py script:
python ../corso_epitrascrittomica/data_reditools/src/REDItools/main/REDItoolDnaRna.py -o /home/student_X/RNAseq -i SRRXXXXXXX.bam -f /data/annotations/GRCh37.primary_assembly.genome.fa -t 4 -c 0,1 -m 0,255 -v 1 -q 0,30 -e -n 0.0 -N 0.0 -u -l -p

python ../corso_epitrascrittomica/data_reditools/src/REDItools/main/REDItoolDnaRna.py -o /home/student_7/RNAseq -i SRR1319672.bam -f /data/annotations/GRCh37.primary_assembly.genome.fa -t 4 -c 0,1 -m 0,255 -v 1 -q 0,30 -e -n 0.0 -N 0.0 -u -l -p



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

