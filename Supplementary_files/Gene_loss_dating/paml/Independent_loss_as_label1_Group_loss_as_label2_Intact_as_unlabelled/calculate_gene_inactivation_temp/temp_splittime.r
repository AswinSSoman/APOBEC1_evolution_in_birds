library(ape)
library(reshape2)
a<-read.tree("/home/neo/bird_db1/aswin/APOBEC1/Dating/paml/Mammal_ADH_IV/PGLS/Mammal_ADH_tree_revisions.nwk")
cophenetic(a)->v
subset(melt(v), value!=0)->v1
v1$value= round(v1$value/2,digit=2)
write.table(v1,file="calculate_gene_inactivation_temp/temp_pairwise.txt",quote=F,sep="\t",row.names=F,col.names=F)
