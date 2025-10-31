###########################################################################################################################################################################################################################################################################################################
#APOBEC1 GENE LOSS VALIDATION: LONG READ MAPPING
###########################################################################################################################################################################################################################################################################################################

###########################################################################################################################################################################################################################################################################################################
#PacBio Long read Gallus gallus

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Send exon (final gallus apobec1 consensus sequence) & intron files

ssh -Y neo@172.30.1.174
cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Gallus_gallus/aswin/APOBEC1
scp manual_consensus.fa sagar@172.30.1.172:/media/sagar/disk4/Gallus_gallus_PacBio_bam/Gal6_pacbio/aswin/
cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/Anas_platyrhynchos/aswin/APOBEC1
scp apobec1_duck_intron.fa sagar@172.30.1.172:/media/sagar/disk4/Gallus_gallus_PacBio_bam/Gal6_pacbio/aswin/

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#New merged data (521Gb Gal6_packbio_merged.sorted.bam)

mkdir -p /media/sagar/disk4/Gallus_gallus_PacBio_bam/Gal6_pacbio/aswin
cd /media/sagar/disk4/Gallus_gallus_PacBio_bam/Gal6_pacbio/aswin

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Create bed files for focal gene
#NOTE: confirm the genome & annotation used for blast is exactly equal to the genome in UCSC at base pair level

cp ../Gallus_gallus_gca000002315v5.GRCg6a.dna.toplevel.fa .
makeblastdb -in Gallus_gallus_gca000002315v5.GRCg6a.dna.toplevel.fa -dbtype nucl -out Gallus_gallus_gca000002315v5.GRCg6a.dna.toplevel.fa
#for exon bed
blastn -task blastn -db Gallus_gallus_gca000002315v5.GRCg6a.dna.toplevel.fa -query manual_consensus.fa -evalue 0.001 -outfmt 6 | awk '!a[$1]++' | awk '{if($9<$10) print$2,$9,$10,$1; else print$2,$10,$9,$1}' OFS="\t" > APOBEC1_exon_2_3_4_5_chicken.bed
#for intron bed (neo@172.30.1.174)
cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/Anas_platyrhynchos/aswin/APOBEC1/
bedtools getfasta -fi ../../genome/GCF_015476345.1_ZJU1.0_genomic.fna -bed <(awk 'NR>1{print$1,$3,p,"intron_"++i,"1","-"}{p=$2}' OFS="\t" test.out.bed) -name -s | tr -d ")(-" > apobec1_duck_intron.fa
cd /media/sagar/disk4/Gallus_gallus_PacBio_bam/Gal6_pacbio/aswin
blastn -task blastn -db Gallus_gallus_gca000002315v5.GRCg6a.dna.toplevel.fa -query apobec1_duck_intron.fa -evalue 0.001 -outfmt 6 | awk '!a[$1]++' | awk '{if($9<$10) print$2,$9,$10,$1; else print$2,$10,$9,$1}' OFS="\t" > APOBEC1_intron_chicken.bed

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Create bed for long reads

#Subset bam to APOBEC1 & it's syntenic region
chr=`awk '$3=="gene" && $0~"protein_coding"' Gallus_gallus_gca000002315v5.GRCg6a.dna.toplevel.gtf | grep "nanog" -iA1 | awk '{print$1}' | sort -u`
cord=`awk '$3=="gene" && $0~"protein_coding"' Gallus_gallus_gca000002315v5.GRCg6a.dna.toplevel.gtf | grep "nanog" -iA1 | awk '{print$4,$5}' OFS="\n" | sort -n | sed -n '1p;$p' | paste -s -d "-"`
samtools view -b ../Gal6_packbio_merged.sorted.bam $(echo -e "$chr:$cord") > Pacbio_apobec1_subseted.bam
#Convert bam to bed
bedtools bamtobed -i Pacbio_apobec1_subseted.bam | awk '{print$1,$2,$3,$4,$3-$2}' OFS="\t" > Pacbio_apobec1_subseted.bed

#Merge bed files
cat APOBEC1_exon_2_3_4_5_chicken.bed | sed '1i track name="APOBEC1 exons" color=255,0,255 visibility=1\n#chrom\tstart\tend\tname' > Pacbio_all_reads_apobec1.bed
cat APOBEC1_intron_chicken.bed | sed '1i track name="APOBEC1 introns" color=255,0,0 visibility=1\n#chrom\tstart\tend\tname' >> Pacbio_all_reads_apobec1.bed
cat Pacbio_apobec1_subseted.bed | sed '1i track name="PacBio chicken" color=0,0,255 visibility=1\n#chrom\tstart\tend\tname\tlength' >> Pacbio_all_reads_apobec1.bed

cat APOBEC1_exon_2_3_4_5_chicken.bed | sed '1i track name="APOBEC1 exons" color=255,0,255 visibility=1\n#chrom\tstart\tend\tname' > Pacbio_reads_longer_than_10kb_apobec1.bed
cat APOBEC1_intron_chicken.bed | sed '1i track name="APOBEC1 introns" color=255,0,0 visibility=1\n#chrom\tstart\tend\tname' >> Pacbio_reads_longer_than_10kb_apobec1.bed
awk '$5>9999' Pacbio_apobec1_subseted.bed | sed '1i track name="PacBio chicken" color=0,0,255 visibility=1\n#chrom\tstart\tend\tname\tlength' >> Pacbio_reads_longer_than_10kb_apobec1.bed

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#UCSC visualization & labelling

#Take screenshot in :
#	1) APOBEC1 with flanking genes : chr1:75,592,831-75,638,416
#	2) APOBEC1 exonic region       : chr1:75,613,821-75,623,539

#plot read distribution online
https://statscharts.com/bar/histogram

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Transfer all output to neo

cd /media/sagar/disk4/Gallus_gallus_PacBio_bam/Gal6_pacbio/aswin
scp -r APOBEC1 neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/Long_read/Gallus_gallus/PacBio

###########################################################################################################################################################################################################################################################################################################
#Nanopore Long read mapping of Gallus gallus 

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download long read data

cat filereport_read_run_SRR13494714_tsv.txt | awk 'FS="\t", OFS="\t" { gsub("ftp.sra.ebi.ac.uk", "era-fasp@fasp.sra.ebi.ac.uk:"); print }' | awk 'NR>1{print$NF}' | awk 'OFS="\n" {print "time /home/ceglab8/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab8/.aspera/connect/etc/asperaweb_id_dsa.openssh" " " $1 " ."}'
time /home/ceglab8/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab8/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR154/042/SRR15421342/SRR15421342_1.fastq.gz .
time /home/ceglab8/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab8/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR154/043/SRR15421343/SRR15421343_1.fastq.gz .
time /home/ceglab8/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab8/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR154/045/SRR15421345/SRR15421345_1.fastq.gz .
time /home/ceglab8/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab8/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR134/014/SRR13494714/SRR13494714_1.fastq.gz .
time /home/ceglab8/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab8/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR134/013/SRR13494713/SRR13494713_1.fastq.gz .
time /home/ceglab8/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab8/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR154/044/SRR15421344/SRR15421344_1.fastq.gz .

esearch -db sra -query "Gallus_gallus[ORGN] AND Oxford Nanopore AND genome" | esummary > r
cat r | grep "SRR[0-9]\+" -o | xargs -n1 sh -c 'sra-meta -id $0'

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Create fasta file of Nanopore
time cat SRR13494714_1.fastq | sed -n '1~4s/^@/>/p;2~4p' >> SRR15421342_3_5.fa

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Map to Gal Galgal6 assembly

#mapping (100m29.334s)
cd /media/sagar/disk4/Gallus_gallus_PacBio_bam/Gal6_pacbio
time minimap2 -t 32 -ax map-ont Gallus_gallus.GRCg6a.dna_sm.toplevel.fa SRR15421342_3_4_5_SRR13494713_4.fa > SRR15421342_3_4_5_SRR13494713_4.sam
#Convert sam to bam (1172m38.775s)
samtools view -bS SRR15421342_3_4_5_SRR13494713_4.sam > SRR15421342_3_4_5_SRR13494713_4.bam
#sorting (242m47.230s)
time samtools sort SRR15421342_3_4_5_SRR13494713_4.bam -o SRR15421342_3_4_5_SRR13494713_4_sorted.bam
rm SRR15421342_3_4_5_SRR13494713_4.bam
#index the sorted bam file (16m11.183s)
time samtools index SRR15421342_3_4_5_SRR13494713_4_sorted.bam SRR15421342_3_4_5_SRR13494713_4_sorted.bam.bai
#Convert bam to tdf file for quick visualization in IGV (217m54.359s)
time /media/aswin/programs/IGV_2.12.3/igvtools count SRR15421342_3_4_5_SRR13494713_4_sorted.bam SRR15421342_3_4_5_SRR13494713_4_sorted.tdf Gallus_gallus.GRCg6a.dna_sm.toplevel.fa

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Check APOBEC1 region

mkdir -p /media/aswin/bird_database/Gallus_gallus_Nanopore_db/aswin/APOBEC1
cd /media/aswin/bird_database/Gallus_gallus_Nanopore_db/aswin/APOBEC1
#Add genome & annotation files inside genome folder

#Send exon (final gallus apobec1 consensus sequence) & intron files from neo
ssh -Y neo@172.30.1.174
cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Gallus_gallus/aswin/APOBEC1
scp manual_consensus.fa ceglab25@172.30.1.131:/media/aswin/bird_database/Gallus_gallus_Nanopore_db/aswin/APOBEC1/
cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/Anas_platyrhynchos/aswin/APOBEC1
scp apobec1_duck_intron.fa ceglab25@172.30.1.131:/media/aswin/bird_database/Gallus_gallus_Nanopore_db/aswin/APOBEC1/

#Create bed files for focal gene
#NOTE: confirm the genome & annotation used for blast is exactly equal to the genome in UCSC at base pair level

#for exon bed
blastn -task blastn -db ../genome/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa -query manual_consensus.fa -evalue 0.01 -outfmt 6 | awk '!a[$1]++' | awk '{if($9<$10) print$2,$9,$10,$1; else print$2,$10,$9,$1}' OFS="\t" > APOBEC1_exon_2_3_4_5_chicken.bed
#for intron bed
blastn -task blastn -db ../genome/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa -query apobec1_duck_intron.fa -evalue 0.01 -outfmt 6 | awk '!a[$1]++' | awk '{if($9<$10) print$2,$9,$10,$1; else print$2,$10,$9,$1}' OFS="\t" > APOBEC1_intron_chicken.bed

#Create bed for long reads
#Subset bam to APOBEC1 & it's syntenic region
chr=`awk '$3=="gene" && $0~"protein_coding"' ../genome/Gallus_gallus.GRCg6a.97.gtf | grep "nanog" -iA1 | awk '{print$1}' | sort -u`
cord=`awk '$3=="gene" && $0~"protein_coding"' ../genome/Gallus_gallus.GRCg6a.97.gtf | grep "nanog" -iA1 | awk '{print$4,$5}' OFS="\n" | sort -n | sed -n '1p;$p' | paste -s -d "-"`
samtools view -b ../../SRR15421342_3_4_5_SRR13494713_4_sorted.bam $(echo -e "$chr:$cord") > Nanopore_apobec1_subseted.bam
samtools index Nanopore_apobec1_subseted.bam Nanopore_apobec1_subseted.bam.bai
#Convert bam to bed
bedtools bamtobed -i Nanopore_apobec1_subseted.bam | awk '{print$1,$2,$3,$4,$3-$2}' OFS="\t" > Nanopore_apobec1_subseted.bed

#Merge bed files
cat APOBEC1_exon_2_3_4_5_chicken.bed | sed '1i track name="APOBEC1 exons" color=255,0,255 visibility=1\n#chrom\tstart\tend\tname' > Nanopore_all_reads_apobec1.bed
cat APOBEC1_intron_chicken.bed | sed '1i track name="APOBEC1 introns" color=255,0,0 visibility=1\n#chrom\tstart\tend\tname' >> Nanopore_all_reads_apobec1.bed
cat Nanopore_apobec1_subseted.bed | sed '1i track name="Nanopore chicken reads" color=0,0,255 visibility=1\n#chrom\tstart\tend\tname\tlength' >> Nanopore_all_reads_apobec1.bed

cat APOBEC1_exon_2_3_4_5_chicken.bed | sed '1i track name="APOBEC1 exons" color=255,0,255 visibility=1\n#chrom\tstart\tend\tname' > Nanopore_reads_longer_than_10kb_apobec1.bed
cat APOBEC1_intron_chicken.bed | sed '1i track name="APOBEC1 introns" color=255,0,0 visibility=1\n#chrom\tstart\tend\tname' >> Nanopore_reads_longer_than_10kb_apobec1.bed
awk '$5>9999' Nanopore_apobec1_subseted.bed | sed '1i track name="Nanopore chicken reads longer than 10kb" color=0,0,255 visibility=1\n#chrom\tstart\tend\tname\tlength' >> Nanopore_reads_longer_than_10kb_apobec1.bed

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#UCSC visualization & labelling

#Take screenshot in :
#	1) APOBEC1 with flanking genes : chr1:75,592,831-75,638,416
#	2) APOBEC1 exonic region       : chr1:75,613,821-75,623,539

#plot read distribution online
https://statscharts.com/bar/histogram

#Transfer all output to neo
cd /media/aswin/bird_database/Gallus_gallus_Nanopore_db/aswin
scp -r APOBEC1 neo@172.30.1.174:/home/neo/bird_db1/aswin/APOBEC1/Long_read/Gallus_gallus/Nanopore

###########################################################################################################################################################################################################################################################################################################
#Map to Galgal 7 assembly

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Mapping using minimap2

cd /media/aswin/bird_database/Gallus_gallus_Nanopore_db
time minimap2 -t 32 -ax map-ont GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna SRR15421342_3_4_5_SRR13494713_4.fa > SRR15421342_3_4_5_SRR13494713_4.sam

#Convert sam to bam (130m38.332s)
samtools view -bS SRR15421342_3_4_5_SRR13494713_4.sam > SRR15421342_3_4_5_SRR13494713_4.bam
#sorting (244m2.818s)
samtools sort SRR15421342_3_4_5_SRR13494713_4.bam -o SRR15421342_3_4_5_SRR13494713_4_sorted.bam
#index the sorted bam file (16m14.023s)
samtools index SRR15421342_3_4_5_SRR13494713_4_sorted.bam SRR15421342_3_4_5_SRR13494713_4_sorted.bam.bai
#Convert bam to tdf file for quick visualization in IGV (217m54.359s)
/media/aswin/programs/IGV_2.12.3/igvtools count SRR15421342_3_4_5_SRR13494713_4_sorted.bam SRR15421342_3_4_5_SRR13494713_4_sorted.tdf GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna

###########################################################################################################################################################################################################################################################################################################
#Leptosomus discolor

#In neo
scp ~/soft_links/Leptosomus_discolor/genome/GCF_000691785.1_ASM69178v1_genomic.fna /home/neo/soft_links/Leptosomus_discolor/genome/GCF_000691785.1_ASM69178v1_genomic.gff ceglab25@172.30.1.131:/media/aswin/Long_read/Leptosomus_discolor/

#In ceglab25
mkdir /media/aswin/Long_read/Leptosomus_discolor
cd /media/aswin/Long_read/Leptosomus_discolor

#Download hifi reads using aws
aws s3 --no-sign-request sync s3://genomeark/species/Leptosomus_discolor/bLepDis1/genomic_data/pacbio_hifi/ .
#Extract and Combine fasta
time zcat m64330e_230511_180019.bc1010--bc1010.hifi_reads.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' | sed '/>/ s!/!_!g' > m64330e_230511_180019.bc1010--bc1010.hifi_reads.fa
time zcat m64055e_230502_084539.bc1010--bc1010.hifi_reads.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' | sed '/>/ s!/!_!g' > m64055e_230502_084539.bc1010--bc1010.hifi_reads.fa
time cat m64330e_230511_180019.bc1010--bc1010.hifi_reads.fa m64055e_230502_084539.bc1010--bc1010.hifi_reads.fa > combined_hifi_reads.fa

#Map using minimap 2 (37m58.498s)
time minimap2 -t 32 -ax map-pb GCF_000691785.1_ASM69178v1_genomic.fna combined_hifi_reads.fa > combined_hifi_reads.sam
#Convert sam to bam (42m1.571s)
time samtools view -bS combined_hifi_reads.sam > combined_hifi_reads.bam
#sorting (66m21.541s)
time samtools sort combined_hifi_reads.bam -o combined_hifi_reads_sorted.bam
rm SRR15421342_3_4_5_SRR13494713_4.bam
#index the sorted bam file (5m54.995s)
time samtools index combined_hifi_reads_sorted.bam combined_hifi_reads_sorted.bam.bai

###########################################################################################################################################################################################################################################################################################################
#Struthio_camelus_australis

#In neo
scp ~/soft_links/Leptosomus_discolor/genome/GCF_000691785.1_ASM69178v1_genomic.fna /home/neo/soft_links/Leptosomus_discolor/genome/GCF_000691785.1_ASM69178v1_genomic.gff ceglab25@172.30.1.131:/media/aswin/Long_read/Leptosomus_discolor/
scp /home/neo/soft_links/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/test.out.bed ceglab25@172.30.1.131:/media/aswin/Long_read/Struthio_camelus/A1_from_cgblast.bed
scp /home/neo/soft_links/Struthio_camelus_australis/aswin/APOBEC1/2nd_gblast/APOBEC1_duplicate/cgblast/APOBEC1_like_Struthio_camelus_1_as_query/Struthio_camelus/test.out.bed ceglab25@172.30.1.131:/media/aswin/Long_read/Struthio_camelus/A1_like_from_cgblast.bed
scp /home/neo/bird_db1/aswin/APOBEC1/TOGA/Struthio_camelus/toga_output/for_vizualization/* ceglab25@172.30.1.131:/media/aswin/Long_read/Struthio_camelus/

#In ceglab25 (37m58.498s)
mkdir /media/aswin/Long_read/Struthio_camelus
cd /media/aswin/Long_read/Struthio_camelus

#Download unmapped PacBio Hifi reads
time aws s3 --no-sign-request cp s3://genomeark/species/Struthio_camelus/bStrCam1/genomic_data/pacbio_hifi/m84091_240223_205537_s4.hifi_reads.bc1017.fastq.gz .	#10m0.324s
time aws s3 --no-sign-request cp s3://genomeark/species/Struthio_camelus/bStrCam1/genomic_data/pacbio_hifi/m84091_240227_220334_s1.hifi_reads.bc1017.fastq.gz .	#12m5.704s

#Extract and Combine fasta
zcat m84091_240223_205537_s4.hifi_reads.bc1017.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' | sed '/>/ s!/!_!g' > m84091_240223_205537_s4.hifi_reads.bc1017.fa	#11m21.822s
zcat m84091_240227_220334_s1.hifi_reads.bc1017.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' | sed '/>/ s!/!_!g' > m84091_240227_220334_s1.hifi_reads.bc1017.fa	#13m1.145s
cat m84091_240223_205537_s4.hifi_reads.bc1017.fa m84091_240227_220334_s1.hifi_reads.bc1017.fa > combined_hifi_reads.fa	#21m58.713s

#Map using minimap 2 (54m49.567s)
time minimap2 -t 32 -ax map-pb GCF_000698965.1_ASM69896v1_genomic.fna combined_hifi_reads.fa > combined_hifi_reads.sam
#Convert sam to bam (72m2.814s)
time samtools view -bS combined_hifi_reads.sam > combined_hifi_reads.bam
#sorting (107m20.051s)
time samtools sort combined_hifi_reads.bam -o combined_hifi_reads_sorted.bam
#index the sorted bam file (9m42.097s)
time samtools index combined_hifi_reads_sorted.bam combined_hifi_reads_sorted.bam.bai

#Subset
samtools view -b combined_hifi_reads_sorted.bam $(echo -e "NW_009272133.1:1-35907") > combined_hifi_reads_sorted_subsetted.bam
#Convert bam to bed
bedtools bamtobed -i combined_hifi_reads_sorted_subsetted.bam | awk '{print$1,$2,$3,$4,$3-$2}' OFS="\t" > combined_hifi_reads_sorted_subsetted.bed
awk '{if(($3-$2)>=2000) print$0}' OFS="\t" combined_hifi_reads_sorted_subsetted.bed > combined_hifi_reads_sorted_subsetted_atleast_2000bp.bed

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Map to updated genome

#NOTE!!!: the pacbio data download was deleted & the data online was also deleted!!!!

cd /media/aswin/Long_read/Struthio_camelus/updated_genome

#Map using minimap 2 (54m49.567s)
time minimap2 -t 32 -ax map-pb GCF_040807025.1_bStrCam1.hap1_genomic.fna ../combined_hifi_reads.fa > combined_hifi_reads.sam
#Convert sam to bam (72m2.814s)
time samtools view -bS combined_hifi_reads.sam > combined_hifi_reads.bam
#sorting (107m20.051s)
time samtools sort combined_hifi_reads.bam -o combined_hifi_reads_sorted.bam
#index the sorted bam file (9m42.097s)
time samtools index combined_hifi_reads_sorted.bam combined_hifi_reads_sorted.bam.bai

#Subset
samtools view -b combined_hifi_reads_sorted.bam $(echo -e "NW_009272133.1:1-35907") > combined_hifi_reads_sorted_subsetted.bam
#Convert bam to bed
bedtools bamtobed -i combined_hifi_reads_sorted_subsetted.bam | awk '{print$1,$2,$3,$4,$3-$2}' OFS="\t" > combined_hifi_reads_sorted_subsetted.bed
awk '{if(($3-$2)>=2000) print$0}' OFS="\t" combined_hifi_reads_sorted_subsetted.bed > combined_hifi_reads_sorted_subsetted_atleast_2000bp.bed



###########################################################################################################################################################################################################################################################################################################
#Malaclemys_terrapin

mkdir /media/aswin/Long_read/Malaclemys_terrapin
cd /media/aswin/Long_read/Malaclemys_terrapin

#Download PacBio
time aws s3 --no-sign-request cp s3://genomeark/species/Malaclemys_terrapin/rMalTer1/genomic_data/pacbio_hifi/m64330e_220320_051323.demultiplex.bc1022--bc1022.hifi_reads.fastq.gz .
time aws s3 --no-sign-request cp s3://genomeark/species/Malaclemys_terrapin/rMalTer1/genomic_data/pacbio_hifi/m64055e_220429_090631.demultiplex.bc1011--bc1011.hifi_reads.fastq.gz .
time aws s3 --no-sign-request cp s3://genomeark/species/Malaclemys_terrapin/rMalTer1/genomic_data/pacbio_hifi/m64334e_220408_155659.demultiplex.bc1011--bc1011.hifi_reads.fastq.gz .

time zcat m64330e_220320_051323.demultiplex.bc1022--bc1022.hifi_reads.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' | sed '/>/ s!/!_!g' > m64330e_220320_051323.demultiplex.bc1022--bc1022.hifi_reads.fa
time zcat m64055e_220429_090631.demultiplex.bc1011--bc1011.hifi_reads.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' | sed '/>/ s!/!_!g' > m64055e_220429_090631.demultiplex.bc1011--bc1011.hifi_reads.fa
time zcat m64334e_220408_155659.demultiplex.bc1011--bc1011.hifi_reads.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' | sed '/>/ s!/!_!g' > m64334e_220408_155659.demultiplex.bc1011--bc1011.hifi_reads.fa
time cat m64330e_220320_051323.demultiplex.bc1022--bc1022.hifi_reads.fa m64055e_220429_090631.demultiplex.bc1011--bc1011.hifi_reads.fa m64334e_220408_155659.demultiplex.bc1011--bc1011.hifi_reads.fa > combined_hifi_reads.fa

time datasets download genome accession GCF_027887155.1 --include gff3,genome

#Map using minimap 2 (54m49.567s)
time minimap2 -t 32 -ax map-pb GCF_027887155.1_rMalTer1.hap1_genomic.fna combined_hifi_reads.fa > combined_hifi_reads.sam
#Convert sam to bam (72m2.814s)
time samtools view -bS combined_hifi_reads.sam > combined_hifi_reads.bam
#sorting (107m20.051s)
time samtools sort combined_hifi_reads.bam -o combined_hifi_reads_sorted.bam
#index the sorted bam file (9m42.097s)
time samtools index combined_hifi_reads_sorted.bam combined_hifi_reads_sorted.bam.bai

#View in IGV
cd /media/aswin/Long_read/Malaclemys_terrapin
samtools faidx GCF_027887155.1_rMalTer1.hap1_genomic.fna
#Visualize the apobec1 region with syntenic genes & note the co-ordinates:
#NC_071505.1:135,198,538-135,412,777
samtools view -b combined_hifi_reads_sorted.bam $(echo -e "NC_071505.1:135198538-135412777") > combined_hifi_reads_sorted_subsetted.bam
#Convert bam to bed
bedtools bamtobed -i combined_hifi_reads_sorted_subsetted.bam | awk '{print$1,$2,$3,$4,$3-$2}' OFS="\t" > combined_hifi_reads_sorted_subsetted.bed
awk '{if(($3-$2)>=2000) print$0}' OFS="\t" combined_hifi_reads_sorted_subsetted.bed > combined_hifi_reads_sorted_subsetted_atleast_2000bp.bed

###########################################################################################################################################################################################################################################################################################################
#Alligator_mississippiensis

mkdir /media/aswin/Long_read/Alligator_mississippiensis
cd /media/aswin/Long_read/Alligator_mississippiensis

#Download PacBio
time aws s3 --no-sign-request cp s3://genomeark/species/Alligator_mississippiensis/rAllMis1/genomic_data/pacbio_hifi/m54306Ue_221121_033625.bc1008--bc1008.hifi_reads.fastq.gz .
time aws s3 --no-sign-request cp s3://genomeark/species/Alligator_mississippiensis/rAllMis1/genomic_data/pacbio_hifi/m54306Ue_221128_200454.bc1008--bc1008.hifi_reads.fastq.gz .
time aws s3 --no-sign-request cp s3://genomeark/species/Alligator_mississippiensis/rAllMis1/genomic_data/pacbio_hifi/m64055e_221208_043940.bc1008--bc1008.hifi_reads.fastq.gz .
time aws s3 --no-sign-request cp s3://genomeark/species/Alligator_mississippiensis/rAllMis1/genomic_data/pacbio_hifi/m64055e_221016_050500.bc1008--bc1008.hifi_reads.fastq.gz .

time aws s3 --no-sign-request cp s3://genomeark/species/Alligator_mississippiensis/rAllMis2/genomic_data/pacbio_hifi/m64055e_220630_180815.hifi_reads.bc1018--bc1018.fastq.gz .
time aws s3 --no-sign-request cp s3://genomeark/species/Alligator_mississippiensis/rAllMis2/genomic_data/pacbio_hifi/m64055e_220624_055853.hifi_reads.bc1018--bc1018.fastq.gz .
time aws s3 --no-sign-request cp s3://genomeark/species/Alligator_mississippiensis/rAllMis2/genomic_data/pacbio_hifi/m64055e_220625_165432.hifi_reads.bc1018--bc1018.fastq.gz .
time aws s3 --no-sign-request cp s3://genomeark/species/Alligator_mississippiensis/rAllMis2/genomic_data/pacbio_hifi/m64055e_220622_190437.hifi_reads.bc1018--bc1018.fastq.gz .
time aws s3 --no-sign-request cp s3://genomeark/species/Alligator_mississippiensis/rAllMis2/genomic_data/pacbio_hifi/m54306Ue_220713_051408.demultiplex.bc1018--bc1018.hifi_reads.fastq.gz .
time aws s3 --no-sign-request cp s3://genomeark/species/Alligator_mississippiensis/rAllMis2/genomic_data/pacbio_hifi/m64055e_220613_163322.demultiplex.bc1018--bc1018.hifi_reads.fastq.gz .

for i in $(ls *.gz)
do
i1=$(echo $i | sed 's/fastq.gz/fa/g')
echo $i1
time zcat $i | sed -n '1~4s/^@/>/p;2~4p' | sed '/>/ s!/!_!g' > $i1
unset i1
done

time cat *.fa > combined_hifi_reads.fa

#Download genome
#datasets download genome accession GCF_030867095.1 --include gff3,genome,seq-report
#cp ncbi_dataset/data/GCF_030867095.1/GCF_030867095.1_rAllMis1_genomic.fna .
#mv ncbi_dataset/data/GCF_030867095.1/genomic.gff GCF_030867095.1_rAllMis1_genomic.gff

#Map using minimap 2 (54m49.567s)
time minimap2 -t 32 -ax map-pb GCF_030867095.1_rAllMis1_genomic.fna combined_hifi_reads.fa > combined_hifi_reads.sam
#Convert sam to bam (72m2.814s)
time samtools view -bS combined_hifi_reads.sam > combined_hifi_reads.bam
#sorting (107m20.051s)
time samtools sort combined_hifi_reads.bam -o combined_hifi_reads_sorted.bam
#index the sorted bam file (9m42.097s)
time samtools index combined_hifi_reads_sorted.bam combined_hifi_reads_sorted.bam.bai


###########################################################################################################################################################################################################################################################################################################
#Alligator_sinensis

mkdir /media/aswin/Long_read/Alligator_sinensis
cd /media/aswin/Long_read/Alligator_sinensis

#Download long read data (In ceglab8)
time prefetch --max-size 100000000 ERR12708335
fasterq-dump --fasta ERR12708335.sra
mv ERR12708335.fasta ERR12708335.fa

#Get genome
scp neo@172.30.1.174:/home/neo/non_mammals_except_birds_db1/amniota/Crocodylia/Alligator_sinensis/genome/GCF_000455745.1_ASM45574v1_genomic.fna .
scp neo@172.30.1.174:/home/neo/non_mammals_except_birds_db1/amniota/Crocodylia/Alligator_sinensis/genome/GCF_000455745.1_ASM45574v1_genomic.gff .
samtools faidx GCF_000455745.1_ASM45574v1_genomic.fna

#Map using minimap 2 (229m46.205s)
time minimap2 -t 32 -ax map-ont GCF_000455745.1_ASM45574v1_genomic.fna ERR12708335.fa > ERR12708335.sam
#Convert sam to bam (41m15.215s)
time samtools view -bS ERR12708335.sam > ERR12708335.bam
#sorting (51m8.463s)
time samtools sort ERR12708335.bam -o ERR12708335_sorted.bam
#index the sorted bam file (4m7.799s)
time samtools index ERR12708335_sorted.bam ERR12708335_sorted.bam.bai



##

