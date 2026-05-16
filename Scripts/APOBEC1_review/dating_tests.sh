#Run paml

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Check scripts used in COA1 gene

git clone --filter=blob:none --no-checkout https://github.com/ceglabsagarshinde/COA1_GENE.git
cd COA1_GENE
git sparse-checkout init --cone
git sparse-checkout set Geneloss_timing
git checkout 01b43da88117ed8249d073f9e1b78ddb4732ba83

cd ~/bird_db1/aswin/APOBEC1/Dating/paml/COA1_GENE
find Geneloss_timing/Galliformes/ -name "*.ctl" -type f | sort -V | paste -d " " - - | xargs -n2 sh -c 'echo ">"$0 $1; paste $0 $1 | nl' | less
find Geneloss_timing/Galliformes/ -name "*.ctl" -type f | sort -V | paste -d " " - - | xargs -n2 sh -c 'echo ">"$0 $1; paste $0 $1 -d "\t" | tr " " "_"' | sed 's!Geneloss_timing/Galliformes/!!g' | colnum.sh
find Geneloss_timing/Galliformes/ -name "*.ctl" -type f | sort -V | paste -d " " - - | xargs -n2 sh -c 'echo ">"$0 $1; paste $0 $1' | sed 's!Geneloss_timing/Galliformes/!!g' > all_galliformes_control_files
cat Geneloss_timing/Galliformes/AllPseudogene_AllMix_AllFunctional/Galliformes_F1x4.ctl| cut -f1 -d "=" | tr -d " " > ctl_params

for p in $(cat ctl_params | egrep -v "seqfile|treefile|outfile" | sort -V)
do
echo ">"$p
grep -w $p all_galliformes_control_files | sort | uniq -c
done | less

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run codeml

mkdir -p ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/F1X4
cd ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix/F1X4
cat /home/neo/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_labeled.nwk | sed 's/:[0-9.]\+//g' | sed 's/#/ #/g' > apobec1_final_align_NT_unroot_labeled_branch_labels_removed.nwk

#Labelling: terminal only of all independent single branch loss as #1, clade loss as #2 & rest unlabelled
BASE=~/bird_db1/aswin/APOBEC1/Dating
TREE="$BASE/tree/readd/apobec1_final_align_NT_unroot_labeled.nwk"
ALN="$BASE/alignment/readd/readd_macse/apobec1_final_align_NT.aln.phy"
CTL="$BASE/paml/COA1_GENE/Geneloss_timing/Galliformes/AllPseudogene_AllMix_AllFunctional"
CODEML=~/programmes/paml-4.10.10-linux-x86_64/bin/codeml
OUT="$BASE/paml/indL1_grpL2_intU"

# make branch-label removed tree
mkdir -p "$OUT"
sed 's/:[0-9.]\+//g; s/#/ #/g' "$TREE" > "$OUT/apobec1_final_align_NT_unroot_labeled_nobranch.nwk"

for m in F1x4 F3x4; do
  for t in "$TREE" "$OUT/apobec1_final_align_NT_unroot_labeled_nobranch.nwk"; do
    [[ "$t" == *nobranch* ]] && suf="branch_labels_removed_$m" || suf="$m"
    d="$OUT/$suf"; mkdir -p "$d"; cd "$d"
    cp "$t" "$ALN" .
    cp "$CTL/Galliformes_${m^}.ctl" "$m.ctl"
    sed -i -e 's|Galliformes.aln|apobec1_final_align_NT.aln.phy|g' -e "s|Galliformes.nwk|$(basename "$t")|g" -e "s|omega_mix_functional_${m^}|paml_out_$m|g" "$m.ctl"
    time "$CODEML" "$m.ctl" > run.stdout &
  done
done
wait


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run codeml based on tree manually lablled using hyphy tool later converted to paml

/home/neo/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk

cd ~/bird_db1/aswin/APOBEC1/Dating/paml

BASE=~/bird_db1/aswin/APOBEC1/Dating
TREE="$BASE/tree/readd/apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk"
ALN="$BASE/alignment/readd/readd_macse/apobec1_final_align_NT.aln.phy"
CTL="$BASE/paml/COA1_GENE/Geneloss_timing/Galliformes/AllPseudogene_AllMix_AllFunctional"
CODEML=~/programmes/paml-4.10.10-linux-x86_64/bin/codeml
OUT="$BASE/paml/mixL1_psL2_intU"

# make branch-label removed tree
mkdir -p "$OUT"
sed 's/:[0-9.]\+//g' "$TREE" > "$OUT/apobec1_final_align_NT_unroot_labeled_nobranch.nwk"

#Labelling: terminal & internal: mixed as #2, pseudo as #1, rest unlablled (52.15 mins)
start_time=$(date +%s)
for m in F1x4 F3x4; do
  for t in "$TREE" "$OUT/apobec1_final_align_NT_unroot_labeled_nobranch.nwk"; do
    [[ "$t" == *nobranch* ]] && suf="branch_labels_removed_$m" || suf="$m"
    d="$OUT/$suf"; mkdir -p "$d"; cd "$d"
    cp "$t" "$ALN" .
    cp "$CTL/Galliformes_${m^}.ctl" "$m.ctl"
    sed -i -e 's|Galliformes.aln|apobec1_final_align_NT.aln.phy|g' -e "s|Galliformes.nwk|$(basename "$t")|g" -e "s|omega_mix_functional_${m^}|paml_out_$m|g" "$m.ctl"
    time "$CODEML" "$m.ctl" > run.stdout &
  done
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calculate gene inactivation times
/home/neo/bird_db1/aswin/APOBEC1/Dating/paml/Mammal_ADH_IV/Dating/codeml/F3X4_model/calculate_gene_inactivation.sh


cd ~/bird_db1/aswin/APOBEC1/Dating/paml/all_mix
cp ~/bird_db1/aswin/APOBEC1/Dating/tree/readd/{fg1.txt,fg2.txt,bg.txt} .

cd ~/bird_db1/aswin/APOBEC1/Dating/paml/mixL1_psL2_intU
for i in $(cat ../funtional)
do
~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh Gallus_gallus $i F1x4/paml_out_F1x4 F3x4/paml_out_F3x4 -wp=1 ~/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk -s | grep -v "Mixed_branch_length"
done | sed '1i Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' | column -t > gallus_gallus_gene_loss_date_wrt_diff_functional_branches.out

time for sp in $(cat ../all_lost)
do
gr=$(grep $sp ~/bird_db1/aswin/taxonomy/orders_all_birds | awk '{print$2}')
~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh $sp -f F1x4/paml_out_F1x4 F3x4/paml_out_F3x4 -wp=1 ~/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk -s | grep -v Mixed_branch_length | sed "s/^/$gr\t/g"
unset gr
done | sed '1i Group Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' > all_gene_loss_dates.tsv

cd ~/bird_db1/aswin/APOBEC1/Dating/paml/mixL1_psL2_intU
time for sp in $(cat ../all_lost)
do
gr=$(grep $sp ~/bird_db1/aswin/taxonomy/orders_all_birds | awk '{print$2}')
~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh $sp -f F1x4/paml_out_F1x4 F3x4/paml_out_F3x4 -wp=1 ~/bird_db1/aswin/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_hyphy_labelled_converted_to_paml.nwk -s | grep -v Mixed_branch_length | sed "s/^/$gr\t/g"
unset gr
done | sed '1i Group Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' > all_gene_loss_dates.tsv

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#In neo
cd ~/bird_db1/aswin/APOBEC1/Dating/tree
time scp -r ../../Dating/ ceglab25@172.28.65.125:/media/aswin/gene_loss/APOBEC1/

#In ceglab25
cd /media/aswin/gene_loss/APOBEC1/Dating

#codeml file
/media/aswin/programs/paml-4.10.10-linux-x86_64/bin/codeml

cd /media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted
cp ../../tree/readd/apobec1_final_align_NT.nwk .
cp ../Mixed_as_label1_Pseudo_as_label1_Inatct_as_unlabelled/F3x4/apobec1_final_align_NT.aln.phy .
cp ../Mixed_as_label1_Pseudo_as_label1_Inatct_as_unlabelled/F3x4/F3x4.ctl
#Label
nano apobec1_final_align_NT_labelled.nwk
sed 's/{mi}/ #1/g' apobec1_final_align_NT_hyphy_labelled.nwk | sed 's/{ps}/ #2/g' > apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk
sed 's/:[0-9.]\+//g' apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk > apobec1_final_align_NT_hyphy_labelled_converted_to_paml_branch_labels_removed.nwk

base=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted
tree1=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk
tree2=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/apobec1_final_align_NT_hyphy_labelled_converted_to_paml_branch_labels_removed.nwk
codeml=/media/aswin/programs/paml-4.10.10-linux-x86_64/bin/codeml
ctl=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/F3x4.ctl
aln_path=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/apobec1_final_align_NT.aln.phy
aln=$(echo $aln_path | awk -F "/" '{print$NF}')

cd $base
for i in F1X4 F3X4
do
for j in $tree1 $tree2
do
n=$(echo $j | awk -F "/" '{print$NF}')
mkdir $i"_"$n && cd "$i"_"$n"
if [[ "$i" == "F1X4" ]]
then
sed -e "s/seqfile =.*/seqfile = $aln/g" -e "s/treefile =.*/treefile = $n/g" -e "s/outfile =.*/outfile = paml_out/g" -e "s/CodonFreq =.*/CodonFreq = 1/g" $ctl > control.ctl
elif [[ "$i" == "F3X4" ]]
then
sed -e "s/seqfile =.*/seqfile = $aln/g" -e "s/treefile =.*/treefile = $n/g" -e "s/outfile =.*/outfile = paml_out/g" -e "s/CodonFreq =.*/CodonFreq = 2/g" $ctl > control.ctl
else :
fi
cp "$aln_path" "$j" .
#$codeml control.ctl > run.stdout
cd $base
done
done

#Run codeml (51.6667)
start_time=$(date +%s)
for i in F1X4 F3X4
do
  for j in $tree1 $tree2
  do
    n=$(basename "$j")
    (
      cd "${i}_${n}" || exit
      "$codeml" control.ctl &> run.stdout
    ) &
  done
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Estimate gene loss timing based on Meredith formula

#Using tree without branch labels
time for sp in $(cat ../all_lost)
do
gr=$(grep $sp ~/bird_db1/aswin/taxonomy/orders_all_birds | awk '{print$2}')
	~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh $sp -f F1X4_apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk/paml_out F3X4_apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk/paml_out -wp=1 \
F1X4_apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk/apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk -s | grep -v Mixed_branch_length | sed "s/^/$gr\t/g"
unset gr
done | sed '1i Group Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' > all_gene_loss_dates.tsv

#Using tree with branch labels
time for sp in $(cat ../all_lost)
do
gr=$(grep $sp ~/bird_db1/aswin/taxonomy/orders_all_birds | awk '{print$2}')
	~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh $sp -f F1X4_apobec1_final_align_NT_hyphy_labelled_converted_to_paml_branch_labels_removed.nwk/paml_out F3X4_apobec1_final_align_NT_hyphy_labelled_converted_to_paml_branch_labels_removed.nwk/paml_out -wp=1
	F1X4_apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk/apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk -s | grep -v Mixed_branch_length | sed "s/^/$gr\t/g"
unset gr
done | sed '1i Group Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' > branch_labels_removed_all_gene_loss_dates.tsv

################################################################################################################################################################################################################################################################################################

cd /media/aswin/gene_loss/APOBEC1/Dating/tree/readd
/media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree  apobec1_final_align_NT_unroot.nwk --label fu --list lost_palaeognathae --output apobec1_final_align_NT_unroot_palaeognathe_as_label1.nwk
/media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree  apobec1_final_align_NT_unroot_palaeognathe_as_label1.nwk --label ps --list intact_apobec1_filtered --output apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3.nwk

mkdir -p /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2.nwk
sed 's/{mi}/ #2/g' apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2.nwk | sed 's/{ps}/ #3/g' | sed 's/{fu}/ #1/g' > apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk
sed 's/:[0-9.]\+//g' apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk > apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk

base=/media/aswin/gene_loss/APOBEC1/Dating/paml/paleoL1_mixL2_intL2_unrooted
tree1=/media/aswin/gene_loss/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk
tree2=/media/aswin/gene_loss/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk
codeml=/media/aswin/programs/paml-4.10.10-linux-x86_64/bin/codeml
ctl=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/F3x4.ctl
aln_path=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/apobec1_final_align_NT.aln.phy
aln=$(echo $aln_path | awk -F "/" '{print$NF}')

cd $base
for i in F1X4 F3X4
do
for j in $tree1
do
n=$(echo $j | awk -F "/" '{print$NF}')
mkdir $i"_"$n && cd "$i"_"$n"
if [[ "$i" == "F1X4" ]]
then
sed -e "s/seqfile =.*/seqfile = $aln/g" -e "s/treefile =.*/treefile = $n/g" -e "s/outfile =.*/outfile = paml_out/g" -e "s/CodonFreq =.*/CodonFreq = 1/g" $ctl > control.ctl
elif [[ "$i" == "F3X4" ]]
then
sed -e "s/seqfile =.*/seqfile = $aln/g" -e "s/treefile =.*/treefile = $n/g" -e "s/outfile =.*/outfile = paml_out/g" -e "s/CodonFreq =.*/CodonFreq = 2/g" $ctl > control.ctl
else :
fi
cp "$aln_path" "$j" .
#$codeml control.ctl > run.stdout
cd $base
done
done

#Run codeml (47.6333)
start_time=$(date +%s)
for i in F1X4 F3X4
do
  for j in $tree1 $tree2
  do
    n=$(basename "$j")
    (
      cd "${i}_${n}" || exit
      "$codeml" control.ctl &> run.stdout
    ) &
  done
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#Summary
cd ~/bird_db1/aswin/APOBEC1/Dating/paml/paleoL1_mixL2_intL2_unrooted
time while read species
do
~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh $species Alectura_lathami \
 F1X4_apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk/paml_out \
 F3X4_apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk/paml_out \
 -wp=1,1 -wm=2 -wf=3 \
  F1X4_apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk/apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk \
  F1X4_apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk/apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk -s | grep -v "Functional_branch"
done < ~/bird_db1/aswin/APOBEC1/Dating/paml/lost_galliformes \
 | sed '1i Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' | sed 's/[ \t]\+/\t/g' > all_galliformes_gene_loss_dates.tsv

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Estimate gene loss timing based on Meredith formula

#Using tree without branch labels
cd ~/bird_db1/aswin/APOBEC1/Dating/paml/paleoL1_mixL2_intL2_unrooted
time for sp in $(cat ../all_lost)
do
gr=$(grep $sp ~/bird_db1/aswin/taxonomy/orders_all_birds | awk '{print$2}')
	~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh $sp -f F1X4_apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk/paml_out F3X4_apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk/paml_out -wp=1
	F3X4_apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk/apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk -s | grep -v Mixed_branch_length | sed "s/^/$gr\t/g"
unset gr
done | sed '1i Group Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' > all_gene_loss_dates.tsv

#Using tree with branch labels
time for sp in $(cat ../all_lost)
do
gr=$(grep $sp ~/bird_db1/aswin/taxonomy/orders_all_birds | awk '{print$2}')
	~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh $sp -f F1X4_apobec1_final_align_NT_hyphy_labelled_converted_to_paml_branch_labels_removed.nwk/paml_out F3X4_apobec1_final_align_NT_hyphy_labelled_converted_to_paml_branch_labels_removed.nwk/paml_out -wp=1
	F1X4_apobec1_final_align_NT_hyphy_labelled_converted_to_paml_branch_labels_removed.nwk/apobec1_final_align_NT_hyphy_labelled_converted_to_paml_branch_labels_removed.nwk -s | grep -v Mixed_branch_length | sed "s/^/$gr\t/g"
unset gr
done | sed '1i Group Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' > branch_labels_removed_all_gene_loss_dates.tsv

################################################################################################################################################################################################################################################################################################
#One event: Galliformes Vs intact

mkdir -p /media/aswin/gene_loss/APOBEC1/Dating/paml/galliformes_single_event
cd /media/aswin/gene_loss/APOBEC1/Dating/paml/galliformes_single_event

#Get alignment
cp /media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/apobec1_final_align_NT.aln.phy .

#Label trees
/media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot.nwk \
 --label ps \
 --list /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/lost_galliformes_Meleagris_gallopavo_removed \
 --output apobec1_final_align_NT_unroot_galliformes_as_label1.nwk
/media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree  apobec1_final_align_NT_unroot_galliformes_as_label1.nwk \
 --label fu \
 --list /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/intact_apobec1_filtered \
 --output apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3.nwk

#Manually edit tree to make 1 event for galliformes loss: change branch leading to split of Penelope & other galliformes to mixed
sed 's/{mi}/ #2/g' apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2.nwk | sed 's/{ps}/ #1/g' | sed 's/{fu}/ #3/g' > apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk
awk -iinplace '{while(match($0, /[0-9]+(\.[0-9]+)?[eE][-+]?[0-9]+/)) {val = substr($0, RSTART, RLENGTH)
  $0 = substr($0,1,RSTART-1) sprintf("%.10f", val) substr($0,RSTART+RLENGTH)
} print}' apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk
sed -e 's/Node[0-9]\+//g' -e 's/:[0-9.]\+//g' apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk > apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk

#Set variables
base=/media/aswin/gene_loss/APOBEC1/Dating/paml/galliformes_single_event
tree1=/media/aswin/gene_loss/APOBEC1/Dating/paml/apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2.nwk/apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk
codeml=/media/aswin/programs/paml-4.10.10-linux-x86_64/bin/codeml
ctl=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/F3x4.ctl
aln_path=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/apobec1_final_align_NT.aln.phy
aln=$(echo $aln_path | awk -F "/" '{print$NF}')

cd $base
for i in F1X4 F3X4
do
for j in $tree1
do
n=$(echo $j | awk -F "/" '{print$NF}')
mkdir $i"_"$n && cd "$i"_"$n"
if [[ "$i" == "F1X4" ]]
then
sed -e "s/seqfile =.*/seqfile = $aln/g" -e "s/treefile =.*/treefile = $n/g" -e "s/outfile =.*/outfile = paml_out/g" -e "s/CodonFreq =.*/CodonFreq = 1/g" $ctl > control.ctl
elif [[ "$i" == "F3X4" ]]
then
sed -e "s/seqfile =.*/seqfile = $aln/g" -e "s/treefile =.*/treefile = $n/g" -e "s/outfile =.*/outfile = paml_out/g" -e "s/CodonFreq =.*/CodonFreq = 2/g" $ctl > control.ctl
else :
fi
cp "$aln_path" "$j" .
#$codeml control.ctl > run.stdout
cd $base
done
done

#Run codeml (39.4833)
start_time=$(date +%s)
for i in F1X4 F3X4
do
  for j in $tree1
  do
    n=$(basename "$j")
    (
      cd "${i}_${n}" || exit
      "$codeml" control.ctl &> run.stdout
    ) &
  done
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

################################################################################################################################################################################################################################################################################################
#2 events : Galliformes Vs intact

mkdir -p /media/aswin/gene_loss/APOBEC1/Dating/paml/galliformes_two_events
cd /media/aswin/gene_loss/APOBEC1/Dating/paml/galliformes_two_events

#Get alignment
cp /media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/apobec1_final_align_NT.aln.phy .

#Label trees
/media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot.nwk \
 --label ps \
 --list /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/lost_galliformes_Meleagris_gallopavo_removed \
 --output apobec1_final_align_NT_unroot_galliformes_as_label1.nwk
/media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree  apobec1_final_align_NT_unroot_galliformes_as_label1.nwk \
 --label fu \
 --list /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/intact_apobec1_filtered \
 --output apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3.nwk

#Manually edit tree to make 2 separate events for galliformes loss: change Penelope to mixed, branch leading to other galliformes as well as mixed
nano apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2.nwk
sed 's/{mi}/ #2/g' apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2.nwk | sed 's/{ps}/ #1/g' | sed 's/{fu}/ #3/g' > apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk
awk -iinplace '{while(match($0, /[0-9]+(\.[0-9]+)?[eE][-+]?[0-9]+/)) {val = substr($0, RSTART, RLENGTH)
  $0 = substr($0,1,RSTART-1) sprintf("%.10f", val) substr($0,RSTART+RLENGTH)
} print}' apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk
sed -e 's/Node[0-9]\+//g' -e 's/:[0-9.]\+//g' apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk > apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk

#Set variables
base=/media/aswin/gene_loss/APOBEC1/Dating/paml/galliformes_two_events
tree1=/media/aswin/gene_loss/APOBEC1/Dating/paml/galliformes_two_events/apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk
codeml=/media/aswin/programs/paml-4.10.10-linux-x86_64/bin/codeml
ctl=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/F3x4.ctl
aln_path=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/apobec1_final_align_NT.aln.phy
aln=$(echo $aln_path | awk -F "/" '{print$NF}')

cd $base
for i in F1X4 F3X4
do
for j in $tree1
do
n=$(echo $j | awk -F "/" '{print$NF}')
mkdir $i"_"$n && cd "$i"_"$n"
if [[ "$i" == "F1X4" ]]
then
sed -e "s/seqfile =.*/seqfile = $aln/g" -e "s/treefile =.*/treefile = $n/g" -e "s/outfile =.*/outfile = paml_out/g" -e "s/CodonFreq =.*/CodonFreq = 1/g" $ctl > control.ctl
elif [[ "$i" == "F3X4" ]]
then
sed -e "s/seqfile =.*/seqfile = $aln/g" -e "s/treefile =.*/treefile = $n/g" -e "s/outfile =.*/outfile = paml_out/g" -e "s/CodonFreq =.*/CodonFreq = 2/g" $ctl > control.ctl
else :
fi
cp "$aln_path" "$j" .
#$codeml control.ctl > run.stdout
cd $base
done
done

#Run codeml (37.65)
start_time=$(date +%s)
for i in F1X4 F3X4
do
  for j in $tree1
  do
    n=$(basename "$j")
    (
      cd "${i}_${n}" || exit
      "$codeml" control.ctl &> run.stdout
    ) &
  done
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#Summary
cd ~/bird_db1/aswin/APOBEC1/Dating/paml/two_events
time while read species
do
~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh $species Alectura_lathami \
 F1X4_apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk/paml_out \
 F3X4_apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk/paml_out \
 -wp=1,1 -wm=2 -wf=3 \
 F1X4_apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk/apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk \
 apobec1_final_align_NT_unroot_galliformes_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk -s | grep -v "Functional_branch"
done < ~/bird_db1/aswin/APOBEC1/Dating/paml/lost_galliformes \
 | sed '1i Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' | sed 's/[ \t]\+/\t/g' > all_galliformes_gene_loss_dates.tsv



 ################################################################################################################################################################################################################################################################################################
 #Label whole tree: Loss as #1, Mix as #2 & intact as 3

 mkdir -p /media/aswin/gene_loss/APOBEC1/Dating/paml/lossL1_mixL2_intL3
 cd /media/aswin/gene_loss/APOBEC1/Dating/paml/lossL1_mixL2_intL3

 #Get alignment
 cp /media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/apobec1_final_align_NT.aln.phy .

cp /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/lost_* .

 #Label trees
 /media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot.nwk \
  --label ps \
  --list /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/lost_galliformes_Meleagris_gallopavo_removed \
  --output apobec1_final_align_NT_unroot_galliformes_as_label1.nwk

/media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree apobec1_final_align_NT_unroot_galliformes_as_label1.nwk \
 --label ps \
 --list /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/lost_palaeognathae \
 --output apobec1_final_align_NT_unroot_galliformes_as_label2.nwk

 /media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree apobec1_final_align_NT_unroot_galliformes_as_label2.nwk \
  --label fu \
  --list /media/aswin/gene_loss/APOBEC1/Dating/tree/readd/intact_apobec1_filtered \
  --output lossL1_mixL2_intL3.nwk

 #Manually edit tree to make 2 separate events for galliformes loss: change Penelope to mixed, branch leading to other galliformes as well as mixed
 nano lossL1_mixL2_intL3_hyphy.nwk
 sed 's/{mi}/ #2/g' lossL1_mixL2_intL3_hyphy.nwk | sed 's/{ps}/ #1/g' | sed 's/{fu}/ #3/g' > lossL1_mixL2_intL3_hyphy_converted_to_paml.nwk
 awk -iinplace '{while(match($0, /[0-9]+(\.[0-9]+)?[eE][-+]?[0-9]+/)) {val = substr($0, RSTART, RLENGTH)
   $0 = substr($0,1,RSTART-1) sprintf("%.10f", val) substr($0,RSTART+RLENGTH)
 } print}' lossL1_mixL2_intL3_hyphy_converted_to_paml.nwk
 sed -e 's/Node[0-9]\+//g' -e 's/:[0-9.]\+//g' lossL1_mixL2_intL3_hyphy_converted_to_paml.nwk > lossL1_mixL2_intL3_hyphy_converted_to_paml_branch_labels_removed.nwk

 #Set variables
 base=/media/aswin/gene_loss/APOBEC1/Dating/paml/lossL1_mixL2_intL3
 tree1=/media/aswin/gene_loss/APOBEC1/Dating/paml/lossL1_mixL2_intL3/lossL1_mixL2_intL3_hyphy_converted_to_paml_branch_labels_removed.nwk
 codeml=/media/aswin/programs/paml-4.10.10-linux-x86_64/bin/codeml
 ctl=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/F3x4.ctl
 aln_path=/media/aswin/gene_loss/APOBEC1/Dating/paml/mixL1_psL2_intU_rooted/apobec1_final_align_NT.aln.phy
 aln=$(echo $aln_path | awk -F "/" '{print$NF}')

cd $base
for i in F1X4 F3X4
do
for j in $tree1
do
n=$(echo $j | awk -F "/" '{print$NF}')
mkdir $i"_"$n && cd "$i"_"$n"
if [[ "$i" == "F1X4" ]]
then
sed -e "s/seqfile =.*/seqfile = $aln/g" -e "s/treefile =.*/treefile = $n/g" -e "s/outfile =.*/outfile = paml_out/g" -e "s/CodonFreq =.*/CodonFreq = 1/g" $ctl > control.ctl
elif [[ "$i" == "F3X4" ]]
then
sed -e "s/seqfile =.*/seqfile = $aln/g" -e "s/treefile =.*/treefile = $n/g" -e "s/outfile =.*/outfile = paml_out/g" -e "s/CodonFreq =.*/CodonFreq = 2/g" $ctl > control.ctl
else :
fi
cp "$aln_path" "$j" .
#$codeml control.ctl > run.stdout
cd $base
done
done

#Run codeml (37.65)
start_time=$(date +%s)
for i in F1X4 F3X4
do
  for j in $tree1
  do
    n=$(basename "$j")
    (
      cd "${i}_${n}" || exit
      "$codeml" control.ctl &> run.stdout
    ) &
  done
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#Estimate time
cd ~/bird_db1/aswin/APOBEC1/Dating/paml/lossL1_mixL2_intL3

while read s
do
  s1=$(echo $s | awk '{print$1}')
  s2=$(echo $s | awk '{print$2}')
  s3=$(echo $s | awk '{print$3}')
d=$(~/bird_db1/aswin/APOBEC1/Dating/scripts/calculate_gene_inactivation.sh $s1 $s3 F1X4_lossL1_mixL2_intL3_hyphy_converted_to_paml_branch_labels_removed.nwk/paml_out F3X4_lossL1_mixL2_intL3_hyphy_converted_to_paml_branch_labels_removed.nwk/paml_out -wp=1,1 -wm=2 -wf=3 F1X4_lossL1_mixL2_intL3_hyphy_converted_to_paml_branch_labels_removed.nwk/lossL1_mixL2_intL3_hyphy_converted_to_paml_branch_labels_removed.nwk lossL1_mixL2_intL3_hyphy_converted_to_paml.nwk -s | grep -v "Mixed_branch_length")
echo $s2 "$d"
unset s1 s2 s3
done < ~/bird_db1/aswin/APOBEC1/Dating/paml/all_lost_groups | sed '1i Group Species Functional_branch Mixed_branch_length 1dS_F1X4_Wm 1dS_F1X4_Wf 1dS_F1X4_Wp 1dS_F1X4_Tp 1dS_F3X4_Wm 1dS_F3X4_Wf 1dS_F3X4_Wp 1dS_F3X4_Tp 1dS_Mean_Tp 2dS_F1X4_Tp 2dS_F3X4_Tp 2dS_Mean_Tp' | sed 's/[ \t]/\t/g' > all_gene_loss_dates.out


################################################################################################################################################################################################################################################################################################
#DRAFT

for i in F1X4 F3X4
do
for j in $tree1 $tree2
do
n=$(echo $j | awk -F "/" '{print$NF}')
cd "$i"_"$n"
echo ">"$i $j
$codeml control.ctl &> run.stdout &
cd $base
done
done



cd ~/bird_db1/aswin/APOBEC1/Dating/paml
cat <(cat lost_palaeognathae | awk '{print$1,"palaeognathae"}') <(cat lost_galliformes | awk '{print$1,"galliformes"}') <(cat lost_independent | awk '{print$1,"independent"}') > all_lost_groups
#Manually create the expected timin of loss in each groups & sort them based on time
