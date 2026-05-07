#Run codeml (37.65)
codeml=~/programmes/paml-4.10.10-linux-x86_64/bin/codeml
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

