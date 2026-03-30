################################################################################################################################################################################################################################################
#Comment 10: Line 273, the title is too strong. It is not shown here unequivocally that Apobec1 does not edit RNA in general.
################################################################################################################################################################################################################################################


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get input data

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools

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

#Compress & keep raw data (512m16.535s)
for f in SRR*.fastq; do gzip "$f"; done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Alignment

mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/rna
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/rna

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

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prepare inputs

#Sort bam files
for f in *.sam
do
(
    filename="${f%%.*}"
    echo $filename
    samtools view -bS "$f" | samtools sort -@ 6 -m 2G -o "${filename}.sorted.bam"
) &
done
wait

#Merge bam files (25m35.751s)
time samtools merge SRR28002323_SRR30595317_merged.bam SRR28002323.sorted.bam SRR30595317.sorted.bam -@ 24

#Sort (35m59.601s)
time samtools sort SRR28002323_SRR30595317_merged.bam -@ 30 -m 3G -o SRR28002323_SRR30595317_merged_sorted.bam
#Index (13m29.160s)
samtools index SRR28002323_SRR30595317_merged_sorted.bam


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Identify editing sites

#Using reditools 1 (231m18.319s)
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/editing

#Run REDitools 1 (2400m50.715s)
time for b in $(ls rna/SRR*_Aligned.sortedByCoord.out.bam)
do
p=$(echo $b | cut -f2 -d "/" | cut -f1 -d "_")
echo ">"$b ":" $p
time python2.7 /media/aswin/programs/REDItools/NPscripts/REDItoolDnaRnav13.py -i $b -j dna/SRR28002323_SRR30595317_merged_sorted.bam -o editing/"$p"_editing -f genome/GCA_041920315.1_ASM4192031v1_genomic.fna -t32 -c1,1 -m30,255 -v1 -q30,30 -v1 -e -n0.0 -N0.0 -u -l -p -s2 -g2 -S &> editing/"$p"_run_std.out
done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Filtering:

# 41m32.908s
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools
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
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/editing
time for i in $(find . -name "*_filtered.out")
do
p=$(echo $i | sed 's/\.out//g')
python2.7 /media/aswin/programs/REDItools/accessory/subCount.py "$i" | sed '1i Substitution Read_count Total_reads Percentage' > "$p"_all_subs_readcount.out
python2.7 /media/aswin/programs/REDItools/accessory/subCount2.py "$i" | sed '1i Substitution Site_count Total_sites Percentage' > "$p"_all_subs_sitecount.out
join -1 1 -2 1 "$p"_all_subs_readcount.out "$p"_all_subs_sitecount.out | sed 's/[ ]\+/\t/g' > "$p"_all_subs_count.out
unset p
done 

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/editing	
for i in $(find . -name "*_filtered.out")
do
p=$(echo $i | sed 's/\.out//g')
r=$(echo $i | cut -f2 -d "/" | cut -f1 -d "_")
p1=$(awk -F "," -v r="$r" '$1==r {print$1,$NF}' OFS="\t" /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/ncbi_SRP459583_metadata.csv | tr " " "_")
tail -n +2 "$p"_all_subs_count.out | sort -k 2nr | sed "s/^/$p1 /g" 
unset p r p1
done | sed '1i SRR_ID\tTissue\tSubstitution\tRead_count\tTotal_reads\t%_Read\tSite_count\tTotal_sites\t%_Site' | sed 's/[ ]\+/\t/g' > summary_substitutions.tsv
sed 's/the_large_intestine/large_intestine/g' summary_substitutions.tsv -i

Rscript plot_read_site_count.R summary_substitutions.tsv plot_read_site_count.pdf
Rscript plot_tissue_count.R summary_substitutions.tsv plot_tissue_count.pdf





#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#TO compare: Run REDitools3

#cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools
#time scp genome/GCA_041920315.1_ASM4192031v1_genomic.fna dna/SRR28002323_SRR30595317_merged_sorted.bam rna/SRR28369623_Aligned.sortedByCoord.out.bam workstation@172.30.1.172:~/aswin/RNA_editing/reditools/
#cd ~/aswin/RNA_editing/reditools
#samtools faidx GCA_041920315.1_ASM4192031v1_genomic.fna
#time samtools index SRR28369623_Aligned.sortedByCoord.out.bam					#0m28.911s
#time samtools index SRR28002323_SRR30595317_merged_sorted.bam					#14m45.851s

time python3.10 -m reditools analyze SRR28369623_Aligned.sortedByCoord.out.bam -r GCA_041920315.1_ASM4192031v1_genomic.fna -o rna_reditools3_SRR28369623.out -t 64 -l 1 -s 2 -C -e -q 30 
time python3.10 -m reditools analyze SRR28002323_SRR30595317_merged_sorted.bam -r GCA_041920315.1_ASM4192031v1_genomic.fna -o dna_reditools3_SRR28369623.out -N -t 64 -l 1 -s 2 -C -e -q 255

#65m39.031s
time python3.10 -m reditools analyze SRR28369623_Aligned.sortedByCoord.out.bam -r GCA_041920315.1_ASM4192031v1_genomic.fna -o rna_reditools3_SRR28369623.out -t 64 -s 2 -C -q 30 
#0m0.468s
time python3.10 -m reditools analyze SRR28002323_SRR30595317_merged_sorted.bam -r GCA_041920315.1_ASM4192031v1_genomic.fna -o dna_reditools3_SRR28369623.out -N -t 64 -s 2 -C -q 255 









