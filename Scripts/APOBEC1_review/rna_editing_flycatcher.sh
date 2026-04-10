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
#Sort bam files (34.63 minutes)
start_time=$(date +%s)
for f in *.sam
do
(
    filename="${f%%.*}"
    echo $filename
    samtools view -bS "$f" | samtools sort -@ 6 -m 3G -o "${filename}.sorted.bam"
) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time)) && echo "- " "$j" " : " $elapsed_time "secs"
unset start_time end_time elapsed_time

#Merge bam files (15m37.526s)
time samtools merge -@ 24 -o dna_merged.bam *.bam

#Sort (13m58.534s)
time samtools sort dna_merged.bam -@ 30 -m 3G -o dna_merged_sorted.bam
#Index (10m15.311s)
time samtools index dna_merged_sorted.bam

#Compress & keep raw data (512m16.535s)
#for f in SRR*.fastq; do gzip "$f"; done
#grep -f dna_ids <(find fastp) | grep fastq | xargs rm
grep -f rna_ids <(ls SRR*) | xargs rm

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Identify editing sites

#Using reditools 1
mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/editing

#Run REDitools 1 (94m11.269s)
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis
time for b in $(ls rna/ERR*.bam)
do
p=$(echo $b | cut -f2 -d "/" | cut -f1 -d "_")
echo ">"$b ":" $p
time python2.7 /media/aswin/programs/REDItools/NPscripts/REDItoolDnaRnav13.py -i $b -j dna/dna_merged_sorted.bam -o editing/"$p"_editing -f genome/GCF_000247815.1_FicAlb1.5_genomic.fna -t32 -c1,1 -m30,255 -v1 -q30,30 -v1 -e -n0.0 -N0.0 -u -l -p -s2 -g2 -S &> editing/"$p"_run_std.out
done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Filtering:

# (5m36.601s)
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis
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
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/editing
time for i in $(find . -name "*_filtered.out")
do
p=$(echo $i | sed 's/\.out//g')
python2.7 /media/aswin/programs/REDItools/accessory/subCount.py "$i" | sed '1i Substitution Read_count Total_reads Percentage' > "$p"_all_subs_readcount.out
python2.7 /media/aswin/programs/REDItools/accessory/subCount2.py "$i" | sed '1i Substitution Site_count Total_sites Percentage' > "$p"_all_subs_sitecount.out
join -1 1 -2 1 "$p"_all_subs_readcount.out "$p"_all_subs_sitecount.out | sed 's/[ ]\+/\t/g' > "$p"_all_subs_count.out
unset p
done 

cd /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/editing
for i in $(find . -name "*_filtered.out")
do
p=$(echo $i | sed 's/\.out//g')
r=$(echo $i | cut -f2 -d "/" | cut -f1 -d "_")
p1=$(awk -F "\t" -v r="$r" '$1==r {print$1,$24}' OFS="\t" /media/aswin/gene_loss/APOBEC1/RNA_editing/Ficedula_albicollis/ncbi_ERP001377.tsv | tr " " "_" | cut -f3- -d "_" | cut -f1 -d "-")
tail -n +2 "$p"_all_subs_count.out | sort -k 2nr | sed "s/^/$r $p1 /g" 
unset p r p1
done | sed '1i SRR_ID\tTissue\tSubstitution\tRead_count\tTotal_reads\t%_Read\tSite_count\tTotal_sites\t%_Site' | sed 's/[ ]\+/\t/g' > summary_substitutions.tsv

Rscript /media/aswin/gene_loss/APOBEC1/RNA_editing/plot_read_site_count.R summary_substitutions.tsv plot_read_site_count.pdf
Rscript /media/aswin/gene_loss/APOBEC1/RNA_editing/plot_tissue_count.R summary_substitutions.tsv plot_tissue_count.pdf





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


