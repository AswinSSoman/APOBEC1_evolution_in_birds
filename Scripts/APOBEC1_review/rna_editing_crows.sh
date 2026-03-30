################################################################################################################################################################################################################################################
#Comment 10: Line 273, the title is too strong. It is not shown here unequivocally that Apobec1 does not edit RNA in general.
################################################################################################################################################################################################################################################

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download data for Corvus

#Based on : https://doi.org/10.1038/ncomms13195


mkdir /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus
cd /media/aswin/gene_loss/APOBEC1/RNA_editing/reditools/Corvus

#Get metadata from NCBI & ENA browsers
csvformat -T ncbi_SRP022901.csv | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' > ncbi_SRP022901.tsv
csvformat -T ncbi_SRP342045.csv | tr " " "_" | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' > ncbi_SRP342045.tsv

#Choose SRA Runs to download
awk '{print$1,$2,$4/1000/1000/1000,$6,$7/1024/1024/1024,$15,$18,$19,$22,$23,$27,$33,$32,$42,$45}' ncbi_SRP022901.tsv | sort -k8,8 -k10,10 | egrep -v "Corvus_brachyrhynchos|Corvus_frugilegus|Coloeus_monedula" \
 | sort -k8,8 -k10,10 -k9,9 | grep S_Up_H32 | awk '$2=="WGS"' | awk '{print$1}' > dna_ids
awk '{print$1,$2,$4/1000/1000/1000,$6,$7/1024/1024/1024,$15,$18,$19,$22,$23,$27,$33,$32,$42,$45}' ncbi_SRP022901.tsv | sort -k8,8 -k10,10 | egrep -v "Corvus_brachyrhynchos|Corvus_frugilegus|Coloeus_monedula" \
	| sort -k8,8 -k10,10 -k9,9 | grep S_Up_H32 | awk '$2=="RNA-Seq"' | awk '{print$1}' > rna_ids

#Create aspera links
grep -f dna_ids ena_SRP022901.tsv | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | tr " " "_" | awk '{print$20}' | tr ";" "\n"  | awk '{print"time /home/ceglab25/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab25/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@"$0,"."}' > dna_ascp_links.sh
grep -f rna_ids ena_SRP022901.tsv | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' | tr " " "_" | awk '{print$20}' | tr ";" "\n"  | awk '{print"time /home/ceglab25/.aspera/connect/bin/ascp -k2 -QT -l 300m -P33001 -i /home/ceglab25/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@"$0,"."}' > rna_ascp_links.sh

#Download data using aspera
chmod +x dna_ascp_links.sh rna_ascp_links.sh
#RNA ()
time ./rna_ascp_links.sh
#DNA ()
time ./dna_ascp_links.sh

#Uncompress ()
time for i in $(ls SRR*.fastq.gz); do echo ">"$i; time gzip -d $i; done

#Download genome


