#APOBEC1 identify events : Genomic alignments

#From RTF format extract just alignment data : limitation : alignment width is difficult to change
libreoffice --headless --convert-to txt Duck_lastz.rtf
cat Duck_lastz.txt | tr ":" " " | awk '{if($6!="") print$1,$3,$4,$6; else if ($2!="") print$1,"NA",$2,"NA" ; else print $0}' | sed 's/Anas[^ ]\+/Duck/g' | sed 's/Gallus[^ ]\+/Chicken/g' | column -t | sed '/^Chicken/ s/$/\n/g' > Duck_chicken1.aln

#From clustal format : limitation can't have genomic coordinates & can't reverse complement
cat <(awk '/^anas/ {print$2}' Duck_lastz.aln | paste -s -d "") <(awk '/^gallus/ {printf$2}' Duck_lastz.aln | paste -s -d "") > tmp1
cat <(awk '/^anas/ {print$2}' Duck_lastz.aln | paste -s -d "") <(awk '/^gallus/ {printf$2}' Duck_lastz.aln | paste -s -d "") | sed 's/.\{250\}/& /g' > tmp
for i in `seq 1 125`; do awk -v c="$i" '{print $c}' tmp; done | awk NF | sed -e '1~2s/^/Duck /g' | sed '2~2s/^/Chicken /g' | column -t | sed '/^Chicken/ s/$/\n/g' > Duck_chicken2.aln
rm tmp

#Try
wget -q --header='Content-type:application/json' 'https://rest.ensembl.org/alignment/region/Anas_platyrhynchos_platyrhynchos/1:78929136-78929314?species_set_group=sauropsid'  -O -
wget -q --header='Content-type:application/json' 'https://rest.ensembl.org/alignment/region/Anas_platyrhynchos_platyrhynchos/1:78929136-78929314?method=LASTZ_NET;species_set=Anas_platyrhynchos_platyrhynchos;species_set=gallus_gallus'  -O -
wget -q --header='Content-type:text/x-phyloxml' 'https://rest.ensembl.org/alignment/region/Anas_platyrhynchos_platyrhynchos/1:78929136-78929314?species_set=Anas_platyrhynchos_platyrhynchos;species_set=gallus_gallus' -O -

wget -q --header='Content-type:text/xml' 'https://rest.ensembl.org/alignment/region/Anas_platyrhynchos_platyrhynchos/1:78929136-78929314?method=LASTZ_NET;species_set=Anas_platyrhynchos_platyrhynchos;species_set=gallus_gallus'  -O -
--header='Content-type:text/x-fasta'


#Final
wget -q --header='Content-type:text/xml' 'https://rest.ensembl.org/alignment/region/Anas_platyrhynchos_platyrhynchos/1:78923207-78931771:-1?method=LASTZ_NET;species_set=Anas_platyrhynchos_platyrhynchos;species_set=gallus_gallus'  -O - > f1
cat f1 | egrep "location|mol_seq" | sed 's/<location>//g' | sed 's!</location>!!g' | sed 's/<mol_seq is_aligned="1">//g' | sed 's!</mol_seq>!!g' | column -t | paste -d " " - - | tr ":" " " | sed 's/-/ /1' | awk '{print$1":"$2,$4,$1":"$3}'


#including scientific name
cat f1 | egrep "scientific_name|location|mol_seq" | sed 's/^[ ]\+//g' | sed '/scientific_name/ s/ /_/g' | sed 's/<scientific_name>//g' | sed 's!</scientific_name>!!g' | sed 's/<location>//g' | sed 's!</location>!!g' | sed 's/<mol_seq is_aligned="1">//g' | sed 's!</mol_seq>!!g' | paste -d " " - - - | tr ":" " " | sed 's/-/ /1' | awk '{print$1,$2,$3,$5,$4}'

cat <(sort -k1,1 -k3,3 f2 | awk '/^Anas/ {print$4}' | paste -s -d "") <(sort -k1,1 -k3,3 f2 | awk '/^Gallus/ {print$4}' | paste -s -d "")

cat <(sort -k1,1 -k3,3 f2 | grep "^Anas" | awk 'NR==1{print$1} {printf$4}') <(sort -k1,1 -k3,3 f2 | grep "^Gallus" | awk 'NR==1{print"\n"$1} {printf$4}') 

cat <(sort -k1,1 -k3,3 f2 | grep "^Anas" | awk 'NR==1{print$1} {printf$4} END{print"\n"}') <(sort -k1,1 -k3,3 f2 | grep "^Gallus" | awk 'NR==1{print$1} {printf$4} END{print"\n"}')









AATCCTACCCAACATGCTGAAGTCAACTTC
AATTCTACCCAACACGCTGAAGTCAACTTC

AGGCACTGCTGATAC-TTTATTTCAGGTGGAAGATGCAGCCAGATGACTTCAAATCAAACTATTTGCCTTCCCAACACCCAAAGGTCGTGTACTTGATGTATGCAATCAAGTGGCGCAGAGGTACCATCTGGAGGGGTTGGTGCTCAAAC
ACATA-TGCTGATACTTTTGTTTCAGATGGAAGATTCAGCCAAATGACTTCAAAATGAACTATTTGCCTGACCAACAACTGAAGGTGTTGTACCTTTTCTAAGAACTGAGGTGCAGCACAGGTACTATCTGGTGGAACTGGTGCTAAAAC

AGGCACTGCTGATAC‑TTTATTTCAGGTGGAAGATGCAGCCAGATGACTTCAAATCAAACTATTTGCCTTCCCAACACCCAAAGGTCGTGTACTTGATGTATGCAATCAAGTGGCGCAGAGGTACCATCTGGAGGGGTTGGTGCTCAAACAATCCTACCCAACATGCTGAAGTCAACTTC
ACATA‑TGCTGATACTTTTGTTTCAGATGGAAGATTCAGCCAAATGACTTCAAAATGAACTATTTGCCTGACCAACAACTGAAGGTGTTGTACCTTTTCTAAGAACTGAGGTGCAGCACAGGTACTATCTGGTGGAACTGGTGCTAAAACAATTCTACCCAACACGCTGAAGTCAACTTC



