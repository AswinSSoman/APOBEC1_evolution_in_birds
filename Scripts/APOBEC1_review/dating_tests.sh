
#In neo
cd ~/bird_db1/aswin/APOBEC1/Dating/tree
time scp -r ../../Dating/ ceglab25@172.28.65.125:/media/aswin/gene_loss/APOBEC1/

#In ceglab25
cd /media/aswin/gene_loss/APOBEC1/Dating

#codeml file
/media/aswin/programs/paml-4.10.10-linux-x86_64/bin/codeml

cd /media/aswin/gene_loss/APOBEC1/Dating/paml/rooted_tree_with_Mixed_as_label1_Pseudo_as_label2_Intact_as_unlabelled
cp ../../tree/readd/apobec1_final_align_NT.nwk .
cp ../Mixed_as_label1_Pseudo_as_label1_Inatct_as_unlabelled/F3x4/apobec1_final_align_NT.aln.phy .
cp ../Mixed_as_label1_Pseudo_as_label1_Inatct_as_unlabelled/F3x4/F3x4.ctl
#Label
nano apobec1_final_align_NT_labelled.nwk
sed 's/{mi}/ #1/g' apobec1_final_align_NT_hyphy_labelled.nwk | sed 's/{ps}/ #2/g' > apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk
sed 's/:[0-9.]\+//g' apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk > apobec1_final_align_NT_hyphy_labelled_converted_to_paml_branch_labels_removed.nwk

base=/media/aswin/gene_loss/APOBEC1/Dating/paml/rooted_tree_with_Mixed_as_label1_Pseudo_as_label2_Intact_as_unlabelled
tree1=/media/aswin/gene_loss/APOBEC1/Dating/paml/rooted_tree_with_Mixed_as_label1_Pseudo_as_label2_Intact_as_unlabelled/apobec1_final_align_NT_hyphy_labelled_converted_to_paml.nwk
tree2=/media/aswin/gene_loss/APOBEC1/Dating/paml/rooted_tree_with_Mixed_as_label1_Pseudo_as_label2_Intact_as_unlabelled/apobec1_final_align_NT_hyphy_labelled_converted_to_paml_branch_labels_removed.nwk
codeml=/media/aswin/programs/paml-4.10.10-linux-x86_64/bin/codeml
ctl=/media/aswin/gene_loss/APOBEC1/Dating/paml/rooted_tree_with_Mixed_as_label1_Pseudo_as_label2_Intact_as_unlabelled/F3x4.ctl
aln_path=/media/aswin/gene_loss/APOBEC1/Dating/paml/rooted_tree_with_Mixed_as_label1_Pseudo_as_label2_Intact_as_unlabelled/apobec1_final_align_NT.aln.phy
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

#Run codeml
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

################################################################################################################################################################################################################################################################################################

cd /media/aswin/gene_loss/APOBEC1/Dating/tree/readd
/media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree  apobec1_final_align_NT_unroot.nwk --label fu --list lost_palaeognathae --output apobec1_final_align_NT_unroot_palaeognathe_as_label1.nwk
/media/aswin/programs/hyphy-2.5.70/hyphy /media/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree  apobec1_final_align_NT_unroot_palaeognathe_as_label1.nwk --label ps --list intact_apobec1_filtered --output apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3.nwk

/media/aswin/gene_loss/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2.nwk
sed 's/{mi}/ #2/g' apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2.nwk | sed 's/{ps}/ #3/g' | sed 's/{fu}/ #1/g' > apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk
sed 's/:[0-9.]\+//g' apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk > apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk

base=/media/aswin/gene_loss/APOBEC1/Dating/paml/unrooted_tree_with_palaeognathe_as_label1_intact_as_label3_mi_as_label2
tree1=/media/aswin/gene_loss/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml.nwk
tree2=/media/aswin/gene_loss/APOBEC1/Dating/tree/readd/apobec1_final_align_NT_unroot_palaeognathe_as_label1_intact_as_label3_mi_as_label2_converted_to_paml_branch_labels_removed.nwk
codeml=/media/aswin/programs/paml-4.10.10-linux-x86_64/bin/codeml
ctl=/media/aswin/gene_loss/APOBEC1/Dating/paml/rooted_tree_with_Mixed_as_label1_Pseudo_as_label2_Intact_as_unlabelled/F3x4.ctl
aln_path=/media/aswin/gene_loss/APOBEC1/Dating/paml/rooted_tree_with_Mixed_as_label1_Pseudo_as_label2_Intact_as_unlabelled/apobec1_final_align_NT.aln.phy
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









