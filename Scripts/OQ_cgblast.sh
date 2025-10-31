#1. Genome blast, ORF checking - 11m43.169s

#Run tBlastx if necessary because it is time consuming

qcount=0

while read i
do
j=`echo $i | awk '{print$1}'`

qcount=$((qcount + 1))

echo "  ("$qcount") " $j "(G)"
cd $(echo $i | awk '{print$2}')

#Create folders for results
mkdir -p aswin/$1
cd aswin/$1

#Get the query
find $3/validated_sequences -name "*$2*.fa" -exec cp {} . \;
find $3/validated_sequences_with_flanking_regions/ -name "*$2*_flanking.fa" -exec cp {} . \;
#Query used
q=`find $3/validated_sequences -name "*$2*.fa" | awk -F "/" '{print$NF}'`
#Get the subject & annotation file
g=`find ../../genome/ -maxdepth 1 -name "*genomic.fna"`
a=`find ../../genome/ -maxdepth 1 -name "*genomic.gff"`

#Running Genome BLAST
if [[ $a == "" ]]
then
#blastn
gblast_short $g $q -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank
else
#blastn
gblast_short $g $q -evalue=0.05 -word_size=11 -extend_query -tblastx=yes -eflank $a
fi

#Unset temporary variables
unset j q g a
#Go back to home directory of Gene of interest - customize
cd $3
printf '_%.s' {1..150}; printf '\n\n'

done < <(grep $2 gblast_query_subject_set | awk '{print$2}' | xargs -I {} grep {} /home/neo/bird_db1/aswin/database_details/all_genome_paths)

unset qcount
