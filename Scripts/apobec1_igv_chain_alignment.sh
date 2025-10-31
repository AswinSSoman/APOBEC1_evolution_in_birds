########################################################################################################################################################################################################################################################################################
#IGV snapshot
########################################################################################################################################################################################################################################################################################

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Before running for all make sure all species have the subsetted chain file in psl format

mkdir ~/aswin/APOBEC1/TOGA/Gallus_gallus/GRCg7b/chain/test_chain
cd ~/aswin/APOBEC1/TOGA/Gallus_gallus/GRCg7b/chain/test_chain
cp ../chr1_duck_masked.chr1_chicken_masked.allfilled.chain.gz .
gzip -d chr1_duck_masked.chr1_chicken_masked.allfilled.chain.gz
sed 's/NC_051772_1/NC_051772.1/g' chr1_duck_masked.chr1_chicken_masked.allfilled.chain -i
p1=`awk '$3=="gene" && $0~"protein_coding"' ~/aswin/APOBEC1/TOGA/Gallus_gallus/Duck/GCF_015476345.1_ZJU1.0_genomic.gff | grep "=APOBEC1;" | awk '{print$4,$5}' OFS="\n" | sort -n | sed -n '1p'`
p2=`awk '$3=="gene" && $0~"protein_coding"' ~/aswin/APOBEC1/TOGA/Gallus_gallus/Duck/GCF_015476345.1_ZJU1.0_genomic.gff | grep "=APOBEC1;" | awk '{print$4,$5}' OFS="\n" | sort -n | sed -n '$p'`
~/aswin/programmes/chainFilter -tOverlapStart=$p1 -tOverlapEnd=$p2 chr1_duck_masked.chr1_chicken_masked.allfilled.chain > subsetted.chain
~/aswin/programmes/chainToPslBasic subsetted.chain subsetted.psl

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Summary of subsetted chains

while read line
do
a=`echo "$line" | awk '{print$16}'`
b=`echo "$line" | awk '{print$17}'`
#convert to bed12
~/aswin/programmes/pslToBed <(echo "$line") tmp.bed12
#numer of bed entries intersecting 
c=`bedtools intersect -a tmp.bed12 -b /home/morpheus/aswin/APOBEC1/TOGA/Gallus_gallus/Duck/duck.bed -split | wc -l`
if [[ $c > 0 ]]
then
c1="✅"
a1=`echo $a | awk '{print$1-1}'`
a2=`echo $a | awk '{print$1+1}'`
b1=`echo $b | awk '{print$1-1}'`
b2=`echo $b | awk '{print$1+1}'`
#number of chains
d=`~/aswin/programmes/chainFilter -tStartMin=$a1 -tStartMax=$a2 -tEndMin=$b1 -tEndMax=$b2 subsetted.chain | grep chain | wc -l`
#chain IDs
e=`~/aswin/programmes/chainFilter -tStartMin=$a1 -tStartMax=$a2 -tEndMin=$b1 -tEndMax=$b2 subsetted.chain | grep chain | awk '{print$NF}' | paste -s -d ","`
#chain size
f=`~/aswin/programmes/chainFilter -tStartMin=$a1 -tStartMax=$a2 -tEndMin=$b1 -tEndMax=$b2 subsetted.chain | grep chain | awk '{print$7-$6}' | paste -s -d ","`
else
c1="❌"
d="-"
e="-"
f="-"
fi
echo $a $b $c1 $d $e $f
unset a b c c1 a1 a2 b1 b2 d e f
rm tmp.bed12
done < subsetted.psl | sed '1i Psl_target_Start Psl_target_End intersect #_of_chains Chain_ID Chain_size' | column -t > psl_chain_focal_gene_intersect_summary

cd ~/aswin/APOBEC1/TOGA/Gallus_gallus
cp -r GRCg7b/chain .

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cd ~/aswin/APOBEC1/TOGA/Chlamydotis_macqueenii/chain/test_chain
cp whole_chain.psl Chlamydotis_macqueenii.psl

########################################################################################################################################################################################################################################################################################
#Transfer the subsetted psl files for IGV visualization

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Subset chains if there are very short chains (<1000 bases long)
mkdir ~/aswin/APOBEC1/TOGA/IGV_snapshots

for i in `grep -if all_toga_subjects_in_phylogenetic_order <(ls -d */) | sed 's!/!!g' | egrep -v "Chlamydotis_macqueenii"`
do
echo $i
cd $i/chain/test_chain/
echo -n > chain_subset_for_igv.chain

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Filter chains lesser than 1000bp 
for x in `cat psl_chain_focal_gene_intersect_summary | awk '$NF>1000' | grep -v "Chain_ID" | awk '{print$5}'`
do
~/aswin/programmes/chainFilter -id=$x subsetted.chain >> chain_subset_for_igv.chain
unset x
done
~/aswin/programmes/chainToPslBasic chain_subset_for_igv.chain $i".psl"
cp $i".psl" /home/morpheus/aswin/APOBEC1/TOGA/IGV_snapshots/
cd ~/aswin/APOBEC1/TOGA
done

########################################################################################################################################################################################################################################################################################
#Create batch file for IGV snapshots

counter=0
echo "new" > igv_APOBEC1_chain_align.batch
echo -e "genome /home/morpheus/aswin/APOBEC1/TOGA/Gallus_gallus/Duck/GCF_015476345.1_ZJU1.0_genomic.fna\n" >> igv_APOBEC1_chain_align.batch

while read i
do
echo -e "new\nmaxPanelHeight 2000\ngoto NC_051772.1:79,840,796-79,847,607\nload /home/morpheus/aswin/APOBEC1/TOGA/Gallus_gallus/Duck/APOBEC1_exons.bed\nexpand APOBEC1_exons.bed\nsnapshotDirectory /home/morpheus/aswin/APOBEC1/TOGA/IGV_snapshots/\n"
i1=`echo $i | awk '{print$1}'`
i2=`echo $i | awk '{print$2}'`
for x in `cat all_toga_subjects_in_phylogenetic_order | xargs -n1 bash -c 'ls IGV_snapshots/*.psl | grep $0' | grep -v "Chlamydotis_macqueenii" | awk -F "/" '{print$NF}' | sed -n "$i1, $i2 p"`
do
xp=`find /home/morpheus/aswin/APOBEC1/TOGA/IGV_snapshots/ -name "$x"`
echo "load "$xp
counter=$((counter + 1))
col=`sed -n "$counter p" rgb_colors`
echo 'setColor '$col $x
echo "setTrackHeight 17" $x
echo "squish " $x
unset xp col
done
echo -e "snapshot igv_APOBEC1_chain_alignment.png\n"
done < <(echo 1 42) >> igv_APOBEC1_chain_align.batch

#Due to large number of tracks take screenshot using portrait mode
~/aswin/programmes/IGV_Linux_snapshot/igv.sh -b igv_APOBEC1_chain_align.batch

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cd ~/aswin/APOBEC1/TOGA/Gallus_gallus/Duck
awk '/##gff/,/##species/' GCF_015476345.1_ZJU1.0_genomic.sorted.gff > APOBEC1_region.gff
grep -if <(awk '$3=="gene" && $0~"protein_coding"' GCF_015476345.1_ZJU1.0_genomic.sorted.gff | grep "APOBEC1" -C1 | awk '{print$9}' | awk -F";" '{print$1}' | sed 's/ID=gene-//g' | sed 's/$/;/g') GCF_015476345.1_ZJU1.0_genomic.sorted.gff >> APOBEC1_region.gff

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Manually choose isoforms to exclude & include
cat APOBEC1_region.gff | awk '$3=="exon"' | awk '{print$9}' | awk -F ";" '{print$1}' | cut -f2 -d "-" | sort -u | egrep -v "XM_005012419.5|XM_038173084.1|XM_027449862.2" | xargs -n1 sh -c 'echo $0; sed -i "/$0/d" APOBEC1_region.gff'
sed 's/LOC101795050/NANOG/g' APOBEC1_region.gff -i

echo "new" > igv_APOBEC1_synteny_chain_align.batch
echo -e "genome /home/morpheus/aswin/APOBEC1/TOGA/Gallus_gallus/Duck/GCF_015476345.1_ZJU1.0_genomic.fna\n" >> igv_APOBEC1_synteny_chain_align.batch

counter=0
while read i
do
echo -e "new\nmaxPanelHeight 2000\ngoto NC_051772.1:79,825,134-79,883,876\nload /home/morpheus/aswin/APOBEC1/TOGA/Gallus_gallus/Duck/APOBEC1_region.gff\nexpand APOBEC1_region.gff\nsnapshotDirectory /home/morpheus/aswin/APOBEC1/TOGA/IGV_snapshots/\n"
i1=`echo $i | awk '{print$1}'`
i2=`echo $i | awk '{print$2}'`
for x in `cat all_toga_subjects_in_phylogenetic_order | xargs -n1 bash -c 'ls IGV_snapshots/*.psl | grep $0' | grep -v "Chlamydotis_macqueenii" | awk -F "/" '{print$NF}' | sed -n "$i1, $i2 p"`
do
xp=`find /home/morpheus/aswin/APOBEC1/TOGA/IGV_snapshots/ -name "$x"`
echo "load "$xp
counter=$((counter + 1))
col=`sed -n "$counter p" rgb_colors`
echo 'setColor '$col $x
echo "setTrackHeight 17" $x
echo "squish " $x
unset xp col
done
echo -e "snapshot igv_APOBEC1_synteny_chain_alignment.png\n"
done < <(echo 1 42) >> igv_APOBEC1_synteny_chain_align.batch

~/aswin/programmes/IGV_Linux_snapshot/igv.sh -b igv_APOBEC1_synteny_chain_align.batch

unset counter





