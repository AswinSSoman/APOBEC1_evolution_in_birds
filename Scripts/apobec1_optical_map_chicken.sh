##################################################################################################################################################################################################################################################################################################
#APOBEC1 gene loss : Optical mapping in chicken
##################################################################################################################################################################################################################################################################################################

#Subset genomic region to a single chromosome & keep the chromosome header as simple as possible
cd /media/aswin/optical_maps/APOBEC1/galgal6
awk '/^1\>/ {print">"$0}' RS=">" ../../genome/Gallus_gallus_gca000002315v5.GRCg6a.dna.toplevel.fa | sed '/>/ s/ .*//g' > chr1_chicken.fa
awk '$1~/^#/ {print$0;next} {if($1=="1") print}' ../../genome/Gallus_gallus_gca000002315v5.GRCg6a.dna.toplevel.gtf > chr1_chicken.gtf
grep -P "\tgene\t" chr1_chicken.gtf | cut -f1,4,5,7,9 | sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$6,".",$4,$10,$12,$14 }' | sort -k1,1 -k2,2n | awk '{if ($7=="ensembl") print$1,$2,$3,$4,$5,$6,$8; else print$1,$2,$3,$7,$5,$6,$9}' OFS="\t" > chr1_chicken.bed_tmp

#Convert reference to optical map
#check the nicing enzyme & recognition site
cd /media/aswin/optical_maps
grep "Nickase" *bnx
#check instrument model
grep "#" *.bnx --color=always | grep "InstrumentSerial" -A2 --color=always | sed 's/ /_/g'

#RECOMMENDED : For Saphry/DLE data : 
#Increase the size of RAM - however don't increase ram a lot otherwise the job gets killed (MAX RAM can be mentioned using "-Xmx"; however it is better to not use this since in several tries mentioning this incresed running time)
#Decrease the number of threads - lowering the threads actually increased the running time 
cd ~/aswin/optical_maps/APOBEC1/galgal6/
#running time 5 sec
java -jar /home/sagar/aswin/programs/OMTools-1.4.1a/OMTools.jar FastaToOM --fastain chr1_chicken.fa --enzymestring CTTAAG --refmapout chr1_chicken.ref

#optical map blast
cd /media/aswin/optical_maps/APOBEC1/galgal6
time java -jar /media/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1_chicken.ref --optmapin ../../bGalGal1_Saphyr_DLE1_3273218.bnx --optresout chr1_chicken_3273218.omd --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 32
time java -jar /media/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1_chicken.ref --optmapin ../../bGalGal4_Saphyr_DLE1_3680518.bnx --optresout chr1_chicken_3680518.omd --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 32
time java -jar /media/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1_chicken.ref --optmapin ../../bGalGal5_Saphyr_DLE1_3680519.bnx --optresout chr1_chicken_3680519.omd --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 32 

#Find region of interest
cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Gallus_gallus/aswin/APOBEC1
scp manual_consensus.fa ceglab25@172.30.1.131:/media/aswin/optical_maps/APOBEC1/galgal6/APOBEC1_Gallus_gallus_manual_consensus.fa
makeblastdb -in chr1_chicken.fa -out chr1_chicken.fa -dbtype nucl
blastn -task blastn -evalue 0.01 -db chr1_chicken.fa -query APOBEC1_Gallus_gallus_manual_consensus.fa -num_threads 32 -outfmt "6 qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand"| sed '1i Query\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\n' | column -t > APOBEC1_blast.outfmt6

paste -d ":" <(cat APOBEC1_blast.outfmt6 | awk '!a[$1]++' | awk 'NR>1{print$2}' OFS="\n" | sort -u) \
<(cat APOBEC1_blast.outfmt6 | awk '!a[$1]++' | awk 'NR>1{print$7,$8}' OFS="\n" | sort -n | sed -n '1p;$p' | paste -s -d "-")

#Add annotation of lost gene (else econdition is not written here)
cat APOBEC1_blast.outfmt6 | awk '!a[$1]++' | awk '{if($NF=="plus") a[NR]=$2; b[NR]=$7; c[NR]=$8} END{print a[3],b[2],c[5],"APOBEC1_remnants"}' OFS="\t" > APOBEC1_gene.bed
cat APOBEC1_gene.bed chr1_chicken.bed_tmp | awk '{print$1,$2,$3,$4}' OFS="\t" | sort -k1,1 -k2,2n > chr1_chicken.bed
#rm chr1_chicken.bed_tmp

region=`paste -d ":" <(grep "apobec1" chr1_chicken.bed -iC1 | awk '{print$1}' | sort -u) \
 <(grep "apobec1" chr1_chicken.bed -iC1 | awk '{print$2,$3}' OFS="\n" | sort -n | sed -n '1p;$p' | paste -s -d " " | awk '{print$1,$2+3000}' OFS="-")`

for i in `ls *.omd`
do
j=`echo $i | cut -f3 -d "_"`
java -jar /media/aswin/programs/OMTools-1.4.1a/OMTools.jar OMView --viewrefin chr1_chicken.ref --viewresin $i --viewannoin chr1_chicken.bed --viewregion $region --dnaratio 21 --zoom 0.8 --viewsave . --viewsaveformat OMView_APOBEC1_mapped_"$j".jpg
java -jar /media/aswin/programs/OMTools-1.4.1a/OMTools.jar OMView --viewrefin chr1_chicken.ref --viewresin $i --viewannoin chr1_chicken.bed --viewregion $region --dnaratio 21 --zoom 0.8 --viewunmap true --viewsave . --viewsaveformat OMView_APOBEC1_mapped_unamapped_"$j".jpg
unset j
done

unset region

##################################################################################################################################################################################################################################################################################################
#Map to chicken galgal7

#Subset genomic region to a single chromosome
ssh -Y neo@172.30.1.174
cd /media/neo/5e4dad81-4707-4b68-ac02-d35a20881d06/home/ceglab26/sagar/bird_database/galliforms_database/Gallus_gallus/genome
awk '$1 ~ /^#/ {print $0;next} {if ($1 == "NC_052532.1") print}' GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff | sed 's/NC_052532.1/1/g' > chr1.gff
sortBed -i chr1.gff | gff2bed > chr1.bed
awk '/^NC_052532.1\y/ {print">"$0}' RS=">" GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna | awk '{if($1~/>/) print">1"; else printf$0}' > chr1_chicken.fna
scp chr1.gff chr1_chicken.fna sagar@172.30.1.172:~/aswin/optical_maps/APOBEC1

#Convert 
cd ~/aswin/optical_maps/
time java -jar ~/OMTools-1.4.1a/OMTools.jar FastaToOM --fastain chr1_chicken.fna --enzymestring CTTAAG --refmapout chr1.ref

time java -jar /home/sagar/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1.ref --optmapin bGalGal1_Saphyr_DLE1_3273218.bnx --optresout chr1_3273218.omd --filtermode 1 --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 32
time java -jar /home/sagar/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1.ref --optmapin bGalGal4_Saphyr_DLE1_3680518.bnx --optresout chr1_3680518.omd --filtermode 1 --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 32
time java -jar /home/sagar/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1.ref --optmapin bGalGal5_Saphyr_DLE1_3680519.bnx --optresout chr1_3680519.omd --filtermode 1 --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 32
time java -jar /home/sagar/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1.ref --optmapin bGalGal1_Saphyr_DLE1_3310977.cmap --optresout chr1_3310977.omd --filtermode 1 --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 32
time java -jar /home/sagar/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1.ref --optmapin bGalGal5_Saphyr_DLE1_3680543.cmap --optresout chr1_3680543.omd --filtermode 1 --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 32
time java -jar /home/sagar/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1.ref --optmapin bGalGal4_Saphyr_DLE1_3680542.cmap --optresout chr1_3680542.omd --filtermode 1 --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 32

#View 
java -jar ~/aswin/programs/OMTools-1.4.1a/OMTools.jar OMView --viewrefin chr1.ref --viewresin chr1_3273218.omd chr1_3680518.omd chr1_3680519.omd --viewannoin apobec1_syntenic_genes.bed --viewregion 1:88441465-88637498 --dnaratio 60 --zoom 1.5
java -jar ~/aswin/programs/OMTools-1.4.1a/OMTools.jar OMView --viewrefin chr1.ref --viewresin chr1_3273218.omd --viewannoin apobec1_syntenic_genes.bed --viewregion 1:88441465-88637498


#From other scripts
java -jar OMTools/OMTools.jar OMBlastMapper -refmapin combined.cmap -optmapin combined.cmap -optresout query6.oma -filtermode 1 -alignmentjoinmode 0 -thread 1 -exactmatch false -writeinfo false -writeunmap false -maxalignitem -1 -minconf 0 -fpp 2 -fnp 2 -ear 0.05 -meas 500 -minjoinscore 50 --allowequalrefquery false
java -jar OMTools/OMTools.jar OMBlastMapper -refmapin $input_ref --optmapin $folder/$input_cmap --optresout $folder/$input --writeunmap false --filtermode 1 --alignmentjoinmode 2 --thread 72


1:75878249-75920587 (nanog - aicda)

1:88441465-88637498 (nectin - plcxd2)

##################################################################################################################################################################################################################################################################################################
#Draft scripts

#ceglab25 : ~98 Gb RAM is used with 8 threads - 12.05 hours
time java -Xmx40400000000 -jar /media/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1_chicken.ref --optmapin ../../bGalGal1_Saphyr_DLE1_3273218.bnx --optresout chr1_chicken_3273218.omd --filtermode 1 --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 28 
time java -Xmx40400000000 -jar /media/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1_chicken.ref --optmapin ../../bGalGal4_Saphyr_DLE1_3680518.bnx --optresout chr1_chicken_3680518.omd --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 28 

#sagar : job might get killer if "404000000000" RAM is used
time java -Xmx404000000000 -jar /home/sagar/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1_chicken.ref --optmapin ../../bGalGal1_Saphyr_DLE1_3273218.bnx --optresout chr1_chicken_3273218.omd --filtermode 1 --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 36
time java -Xmx404000000000 -jar /home/sagar/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1_chicken.ref --optmapin ../../bGalGal4_Saphyr_DLE1_3680518.bnx --optresout chr1_chicken_3680518.omd --filtermode 1 --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 36
time java -Xmx404000000000 -jar /home/sagar/aswin/programs/OMTools-1.4.1a/OMTools.jar OMBlastMapper --refmapin chr1_chicken.ref --optmapin ../../bGalGal5_Saphyr_DLE1_3680519.bnx --optresout chr1_chicken_3680519.omd --filtermode 1 --alignmentjoinmode 1 --overlapmergemode 2 --writeunmap false --thread 36








