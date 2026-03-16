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

#Guidelines for RNA editing:

	#Preprocessing:
		#You can and should use adapter trimming as you would for a typical NGS workflow. No special considerations need to be taken for REDItools.
	#Mapping:
		#When aligning with tools like STAR mapper: With any aligner, I strongly encourage using whatever parameters are needed to include MD tags in the output. A BAM file with MD tags will be processed many times faster than one without. For STAR specifically, you'll need this: --outSAMattributes NH HI AS nM MD
	#Strand:
		#need to determine the strandedness of your BAM file and supply that using the --strand option.
	#Finding repeats:
		#The find-repeats tool is legacy code. You do not need to run it. If you want to focus your analysis on repetitive elements (or exclude them), I suggest using Repeat Masker data from UCSC.
	#Running command `analyze` command
		#This is my best translation of the REDI1 parameters from the REDInet tutorial:
		
		`	python3 -m reditools analyze -l 1 -m 255 -me 1 -bq 30 -Men 1 -s 2 -C.
		`
		#You may want to use --threads and --window for parallel processing. And if you want to exclude repetitive regions, the -k option can be used with a BED file.
		If your BAM file has MD tags, do not provide a reference with -r. The reference option will override the presence of MD tags and will slow REDItools down significantly. You also do not need to add "True" after the -C option. Simply adding -C should turn on strand correction.

	`python3 -m reditools analyze <BAM file> -r <hg38_reference.fasta.gz> -o <reditools_output_table.txt> -q 255 -s 2 -C True
	`
	#I made the following modifications from the sample code:
		#Removed arguments where the value matches the default in the REDItools 3 help output in the command line
		#Sample code has -m 255 -> realised this maps to -q MIN_READ_QUALITY in REDItools 3
		#Sample code has -Men 1 -> removed this to use REDItools 3 default, which is 4
		#Sample code has -s 2 -> remove this to use REDItools 3 default, which is 0. This means unstranded, which I understand is the default setting
		#Sample code has -C, realised I need to put True here, based on the REDItools 3 help output

#After suitable low-quality read and adapter trimming (see Note 1), the mapping of DNA-Seq reads in fastq format (.fq) onto the human genome assembly hg19 can be carried out using BWA, after indexing the reference sequence (see Note 3), with the following command line:
	`bwa mem hg19.fa DNA_R1.fq DNA_R2.fq > DNA.sam
	`
	#In this case, the input sequences are paired end (DNA_R1.fq and DNA_R2.fq). At the end of the run, BWA will generate a SAM file containing the aligned reads

#RNA-Seq reads can be more accurately aligned using STAR. In this case, the aligner requires an index which should be previously generated, by providing both the reference sequence and known splice junctions (from RefSeq, Gencode, or other specialized databases) (see Note 4):
	`STAR --genomeDir STAR_hg19_index --readFilesIn RNA_R1.fq RNA_R2.fq --outSAMtype BAM SortedByCoordinate
	`
	#Input files may be paired end, as in the command line provided (RNA_R1.fq and RNA_R2.fq). In contrast with BWA, STAR allows to directly choose BAM as output format, with reads sorted by coordinate.

#DNA-Seq reads in SAM format need to be converted in BAM and sorted, using SAMtools:
	`samtools view -b DNA.sam -o DNA.bam
		samtools sort DNA.bam -o DNA_sorted.bam
	`
#Read groups are generally required in both DNA-Seq and RNA-Seq alignments. We recommend the following command line, adapting input and output file names:
	`java -Xmx8g -jar AddOrReplaceReadGroups.jar INPUT=X_sorted.bam OUTPUT=X_sorted_RG.bam VALIDATION_STRINGENCY=SILENT TMP_DIR=tmp CREATE_INDEX=True SORT_ORDER=coordinate RGID=sample RGLB=sample RGPL=illumina RGSM=sample RGPU=sample
	`



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

