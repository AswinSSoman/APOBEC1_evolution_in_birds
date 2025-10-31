#Script to do quality check of validated exon-wise sequences from OrthoQuest pipeline

#Positional arguments :
#1. path to all exon-wise sequences
#2. path to all exon-wise sequences with flanking regions
#3. Exon-coding set file containing  info about how much exons each species have
#4. gene name
#5. path to main gene results folder

echo "QC of validated exon-wise sequences" | GREP_COLORS="mt=04;33" grep "." --color=always | sed 's/^/\n/g' | sed 's/$/\n/g' > $6

#Visualize the length distribution of each exons of each validated genes

for z in `ls $1/*.fa`; do  awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' $z | awk NF | sed "s|^|$z|g" | sed 's/^.*APOBEC1_//g' | sed 's/.fa/_/g'; done > exon_lengt_dist1
for z in `ls $2/*.fa`; do  awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' $z | awk NF | sed "s|^|$z|g" | sed 's/^.*APOBEC1_//g' | sed 's/.fa/_/g'; done > exon_lengt_dist2
for i in `ls $3 | cut -f1 -d "_"  | xargs -I {} seq {}`
do
j=`echo "exon_"$i`
echo "- - - "$j" -" | GREP_COLORS="mt=07;33" grep "exon_[0-9]\+" --color=always
#There is a limitaion in termgraph plotting : i.e, if all the data have exactly same value then "--width" parameter can't be applied; hence to introduce a difference add a fake data "aswin" & then remove it
termgraph <(grep "$j\b" exon_lengt_dist1 | sed 's/_exon_[0-9]\+//g' | sed '1i aswin\t4') --format '{:.1f}' --width 100 | grep -v "aswin" | sed "s/_$j\b//g" | awk NF | cat -n | sed 's/:/ :/g' | column -t | sed 's/^/   /g'
unset j
done > plot1

for i in `ls *coding_exons_set | cut -f1 -d "_"  | xargs -I {} seq {}`
do
j=`echo "exon_"$i`
echo "- - - "$j" -" | GREP_COLORS="mt=07;33" grep "exon_[0-9]\+" --color=always
#There is a limitaion in termgraph plotting : i.e, if all the data have exactly same value then "--width" parameter can't be applied; hence to introduce a difference add a fake data "aswin" & then remove it
termgraph <(grep "$j\b" exon_lengt_dist2 | sed 's/_exon_[0-9]\+//g' | sed '1i aswin\t4') --format '{:.1f}' --width 100 | grep -v "aswin" | sed "s/_$j\b//g" | awk NF | cat -n | sed 's/:/ :/g' | column -t | sed 's/^/   /g'
unset j
done > plot2

paste plot1 plot2 | column -t | sed "/exon/ s/^/\n$(printf '—%.s' {1..310})\n/g" | awk NF | sed 1d | sed '1i \\n1) Exon length distribution: CDS (left) & CDS with flanking region (right) \n' | GREP_COLORS="mt=36" grep "^1. Exon.*\|$" --color=always >> $6
rm exon_lengt_dist* plot*

#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#

#Align & QC all validated sequences

#DNA alignment
ls $1/*.fa | xargs -n1 sh -c 'echo $0 | cut -f2 -d "/" | cut -f2- -d "_" | sed "s/.fa//g" | sed "s/^/>/g"; grep -v ">" $0 | paste -s -d ""' | sed '/>/! s/[a-z]/\U&/g' > msa.fa
clustalo -i msa.fa -t dna --wrap=270 --outfmt=clu --resno | grep -v "CLUSTAL" | awk NF | sed '1i \\n2) DNA: Multiple Sequence Alignment\n' | GREP_COLORS="mt=04;36" grep "^2) DNA.*\|$" --color=always >> $6

#Outlier detection in DNA alignment
#aligning again because outlier detection script "OD-seq" doesn't support clustal clustal format generated in the above step
muscle -in msa.fa -out msa.aln
~/programmes/OD-Seq/OD-seq -i msa.aln -o msa_outler.out -s 3
cat <(echo -e "Outliers in DNA alignment : ") <(grep ">" msa_outler.out | tr -d ">" | paste -s -d "," | awk '{if($0=="") print "Nothing"; else print$0}') | paste -s -d " " | sed 's/^/ ➤ /g' | sed 's/$/\n/g' | grep "." -z --color=always >> $6

#Proetin alignment
transeq <(ls $1/*.fa | xargs -n1 sh -c 'echo $0 | cut -f2 -d "/" | cut -f2- -d "_" | sed "s/.fa//g" | sed "s/^/>/g"; grep -v ">" $0 | paste -s -d ""') -auto -stdout | sed '/>/! s/[a-z]/\U&/g' | sed '/>/ s/_1$//g' > msa.aa
clustalo -i msa.aa -t protein --wrap=270 --outfmt=clu --resno | grep -v "CLUSTAL" | awk NF | sed '1i 3) Protein: Multiple Sequence Alignment\n' | GREP_COLORS="mt=36" grep "^3) Protein.*\|$" --color=always >> $6
rm msa.fa msa.aln msa_outler.out msa.aa

#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#

#Visualize the length distribution of each introns of each validated genes

#Print intron length of validated queries
while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
cd $(echo $i | awk '{print$2}')/aswin/$4
if [[ -d 2nd_gblast ]]
then
k=`awk '{print$2-a}{a=$3}' 2nd_gblast/test.out.bed | sed 1d | tr -d "-" | paste -s -d " "`
else
k=`awk '{print$2-a}{a=$3}' test.out.bed | sed 1d | tr -d "-" | paste -s -d " "`
fi
echo $j $k
unset j k
cd $5
done < <(grep -if <(cut -f2- -d "_" total_validated_queries | sed 's/.fa//g') /home/neo/bird_db1/aswin/database_details/all_genome_paths) | column -t > intron_length

#Plot intron distributions
intron=1
for i in `awk '{print NF}' intron_length | uniq | xargs -I {} seq 2 {}`
do
echo "- - - intron_"$intron " - " > "intron_"$intron
termgraph <(awk -v n="$i" '{print$1,$n}' intron_length) --format '{:.1f}' --width 100 | awk NF | cat -n | sed 's/:/ :/g' | column -t | sed 's/^/   /g' >> "intron_"$intron
intron=$((intron + 1))
done
ls intron_? | paste -d " " - - | xargs -n2 sh -c 'paste $0 $1' | column -t | sed "/intron_/ s/^/\n$(printf '—%.s' {1..310})\n/g" | awk NF | sed 1d | GREP_COLORS="mt=07;33" grep "intron_[0-9]\+" -z --color=always | awk NF | sed '1i 4) Intron length distribution\n' | GREP_COLORS="mt=36" grep -a "^4) Intron.*\|$" --color=always >> $6
unset intron
rm intron_*


#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------#

#Visualize exon distance from 1st exon

while read i
do
j=`echo $i | awk '{print$1}' | cut -f1,2 -d "_"`
cd $(echo $i | awk '{print$2}')/aswin/$4
if [[ -d 2nd_gblast ]]
then
f=`awk 'NR==1{print$2}' 2nd_gblast/test.out.bed`
k=`awk -v d="$f" '{if($NF=="minus") print d-$2; else print $2-d}' 2nd_gblast/test.out.bed | paste -s -d " "`
else
f=`awk 'NR==1{print$2}' test.out.bed`
k=`awk -v d="$f" '{if($NF=="minus") print d-$2; else print $2-d}' test.out.bed | paste -s -d " "`
fi
echo $j $k
unset j f k
cd $5
done < <(grep -if <(cut -f2- -d "_" total_validated_queries | sed 's/.fa//g') /home/neo/bird_db1/aswin/database_details/all_genome_paths) | column -t > exon_dist

exon=1
for i in `awk '{print NF}' exon_dist | uniq | xargs -I {} seq 2 {}`
do
echo "- - - exon_"$exon " - " > "exon_"$exon
termgraph <(awk -v n="$i" '{print$1,$n}' exon_dist) --format '{:.1f}' --width 100 | awk NF | cat -n | sed 's/:/ :/g' | column -t | sed 's/^/   /g' | awk '{if($NF=="0.0") print$1,$2,$3,"|",$4; else print $0}' | column -t >> "exon_"$exon
exon=$((exon + 1))
done
#ls intron_? | paste -d " " - - | xargs -n2 sh -c 'paste $0 $1' | column -t
ls exon_? | paste -d " " - - | xargs -n2 sh -c 'paste $0 $1' | column -t | sed "/exon_/ s/^/\n$(printf '—%.s' {1..310})\n/g" | awk NF | sed 1d | GREP_COLORS="mt=07;33" grep "exon_[0-9]\+" -z --color=always | awk NF | sed '1i \\n5) Exon positioning: Reltive position of exons w.r.t 1st exon\n' | GREP_COLORS="mt=36" grep -a "^5) Exon.*\|$" --color=always >> $6
unset exon
rm exon_dist* exon_[0-9]

