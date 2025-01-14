library(ape)

### STEP 0 - preparation of all input data ###
setwd('/path/to/working/folder')
tr<-read.tree('reference_tree.tre') #reference tree which will serve as starting point for testing alternative topologies - e.g. ASTRAL species tree
names<-read.csv('names_for_reference_tree.csv') #table for tips renaming if necessary (it usually contains two columns - ref_tree_name: containing the names used in reference tree; gen_tree_name: containing the corresponding names used in gene trees)
rownames(names)<-names$ref_tree_name 
names<-names[tr$tip.label,] #reorder the table according to tips in reference tree
tr$tip.label<-as.character(names$gen_tree_name) #rename the reference's tree tips to names used in gene trees
tr<-makeNodeLabel(tr,prefix="N_") #naming the nodes
all.sub<-subtrees(tr) #list of all subtrees in reference tree

# simple table with all nodes in the reference tree and appropriate taxa corresponding to each node
taxa.in.subtr<-data.frame(matrix(ncol=3,nrow=length(all.sub)))
colnames(taxa.in.subtr)<-c("Node","taxa_in_node","N_of_taxa")
for (i in 1:length(all.sub)){
  klk<-all.sub[[i]]
  taxa.in.subtr[i,"Node"]<-klk$node.label[1]
  taxa.in.subtr[i,"taxa_in_node"]<-paste0(klk$tip.label,collapse=" ")
  taxa.in.subtr[i,"N_of_taxa"]<-length(klk$tip.label)
}

# selection of nodes that represents all taxa of desired taxonomic level (i.e. genera and/or separate lineages in our case)
incl.nodes<-c("N_2","N_20","N_29","N_42","N_43","N_50","N_58","N_66","N_80","N_82")
names(incl.nodes)<-c("Andinia","Andreettaea","Dryadella","Specklinia1","Specklinia2","Specklinia3","Sarcinula","Scaphosepalum","Teagueia","Platystele")

# if some separate lineages are represented by single tip in reference tree, they have to be listed here
fixed<-data.frame(matrix(ncol=3,nrow=3))
colnames(fixed)<-c("Node","Taxon","Genus")
fixed[1,]<-c("N_78a","Platystele_aurea_P835","Rubellia")
fixed[2,]<-c("N_65a","Specklinia_acanthodes_P297","P297")
fixed[3,]<-c("N_49a","Specklinia_cucumeris_P1141","P1141")

# list of all gene trees
a<-list.files(path="/path/to/folder/with/gene/trees",pattern="*.tre",full.names=T)

### STEP 1 - simplification of gene trees to genera/lineage representatives ###
# simple loop for simplification of all gene trees to nsim selections of above defined node representatives
# this step randomly choose one specimen from each node (plus fixed taxa), pruned all gene trees to this selection, and write simplified gene trees and taxa selection to files
nsim<-50
for(z in 1:nsim){
  sel.tax<-data.frame(matrix(nrow=length(incl.nodes),ncol=3))
  colnames(sel.tax)<-c("Node","Taxon","Genus")
  for(j in 1:length(incl.nodes)){
    sel.tax[j,"Node"]<-incl.nodes[j]
    sel.tax[j,"Genus"]<-names(incl.nodes[j])
    sel.tax[j,"Taxon"]<-sample(scan(text=taxa.in.subtr[taxa.in.subtr$Node==incl.nodes[j],"taxa_in_node"],what=" ",quiet=T),1)
  }
  sel.tax<-rbind(sel.tax,fixed)
  rm(tr.gen)
  for(x in a){
    tr.q<-read.tree(x)
    tr.run<-drop.tip(tr.q,setdiff(tr.q$tip.label,as.character(sel.tax$Taxon)))
    sel.tax<-sel.tax[match(tr.run$tip.label,sel.tax$Taxon),]
    tr.run$tip.label<-as.character(sel.tax$Genus)
    if(exists("tr.gen")==FALSE){
      tr.gen<-tr.run
    }else{
      tr.gen<-c(tr.gen,tr.run)
    }
  }
  write.tree(tr.gen,paste0("gene_trees_",z,".nwk"))
  write.csv(sel.tax,paste0("taxa_selection_",z,".csv"))
}

### STEP 2 - recursive view on gene trees via selected taxa ###
# this step looks for congruent topologies irrespective of selection of representatives across all gene trees
# it helps to define alternative scenarios of evolution that are actually present in gene trees

q<-list.files(pattern="taxa_selection",full.names=T)
rm(cons_tre,cons_tre_null,noncons_tre)
for(i in a){
  tr.run<-read.tree(i)
  rm(tr.new.fin)
  for(ts in q){
    sel.i<-read.csv(ts,row.names = 1)
    tr.run2<-drop.tip(tr.run,setdiff(tr.run$tip.label,as.character(sel.i$Taxon)))
    sel.i<-sel.i[match(tr.run2$tip.label,sel.i$Taxon),]
    tr.run2$tip.label<-as.character(sel.i$Genus)
    tr.new<-root(tr.run2,"Andinia")
    if(exists("tr.new.fin")==FALSE){
      tr.new.fin<-tr.new
    }else{
      tr.new.fin<-c(tr.new.fin,tr.new)
    }
  }
  tr.cons<-ape::consensus(tr.new.fin,p=0.5)
  tr.cons2<-root(tr.cons,"Andinia")
  if(tr.cons2$Nnode==nrow(sel.i)-1){
    if(exists("cons_tre")==FALSE){
      cons_tre<-tr.cons2
      tr.cons2$node.label<-NULL
      cons_tre_null<-tr.cons2
    }else{
      cons_tre<-c(cons_tre,tr.cons2)
      tr.cons2$node.label<-NULL
      cons_tre_null<-c(cons_tre_null,tr.cons2)
    }
    write.tree(tr.cons2,paste0("consensus_",unlist(strsplit(i,"/"))[11],".tre")) # use it if you want to write consensual tree for each gene separately
  }else{
    if(exists("noncons_tre")==FALSE){
      noncons_tre<-tr.cons2
    }else{
      noncons_tre<-c(noncons_tre,tr.cons2)
    }
    write.tree(tr.cons2,paste0("nonconsensus_",unlist(strsplit(i,"/"))[11],".tre")) # use it if you want to write non-consensual tree for each gene separately
  }
}
write.tree(cons_tre,"consensus.tre")
write.tree(noncons_tre,"non_consensus.tre")

poss.sc2<-unique.multiPhylo(cons_tre_null)
write.tree(poss.sc2,"possible_scenarios.tre")

# checking all possible scenarios for missing important ones (i.e. ASTRAL and other topologies with a meaningful pattern)
# scenarios are ordered based on decreasing sum of node supports
poss.sc1<-unique.multiPhylo(cons_tre)
hlp<-sapply(poss.sc1,function(x) sum(x$node.label))
hlp2<-order(hlp,decreasing = T)
pdf("scenario_inspection.pdf",height=6,width=6)
par(mfrow=c(4,3),mar=c(1,1.5,1,1.5))
for(y in hlp2){
  plot(poss.sc1[[y]])
  nodelabels((poss.sc1[[y]]$node.label))
}
dev.off()
# if the list of scenarios doesn't contain the important ones, you have to add them manually to "possible_scenarios.tre" file in newick format

### STEP 3 - run bash script for PhyParts calculations for all taxa selections and all scenarios - i.e. 50 selections x nrow of scenarios in "possible_scenarios.tre" ###
# you need to get java PhyParts tool from https://bitbucket.org/blackrim/phyparts.git and user bash script to run it for all defined scenarios "script_gen_tree_analyses.sh"
# it provides a crucial file *.concon.tre for each selection of taxa and each evolutionary scenario

### STEP 4 - calculation of mean support for each node ###
# it calculates node's support across all selections of representatives and each scenario
# it also calculates variation in supports for each node
# in the last part it sorts all alternative scenarios by overall scores (mean support of selected nodes) and print all and twelve the most probable scenarios

b<-list.files(path="./phyparts",pattern=".concon.tre",full.names=T)
vysl.all<-data.frame(matrix(ncol=read.tree(b[1])[[1]]$Nnode,nrow=length(readLines("possible_scenarios.tre"))))
colnames(vysl.all)<-c(paste0("N.",1:ncol(vysl.all)))
rownames(vysl.all)<-c(paste0("scenario.",1:nrow(vysl.all)))
vari.vysl<-vysl.all
for(sc in 1:nrow(vysl.all)){
  c<-grep(paste0("\\.",sc,".concon"),b,value=T)
  vysl<-data.frame(matrix(nrow=length(c),ncol=read.tree(c[1])[[1]]$Nnode))
  row.names(vysl)<-c
  colnames(vysl)<-c(paste0("N.",1:read.tree(c[1])[[1]]$Nnode))
  for(w in 1:length(c)){
    q1<-read.tree(c[w])[[1]]
    q2<-read.tree(c[w])[[2]]
    for(t in 1:q1$Nnode){
      vysl[w,t]<-paste0(q1$node.label[t],"/",q2$node.label[t])
    }
  }
  for(u in 1:ncol(vysl)){
    nnn<-as.numeric(unlist(strsplit(vysl[,u],"/")))
    vysl.all[sc,u]<-paste0(mean(nnn[c(TRUE,FALSE)]),"/",mean(nnn[c(FALSE,TRUE)]))
    vari.vysl[sc,u]<-paste0(min(nnn[c(TRUE,FALSE)]),"-",max(nnn[c(TRUE,FALSE)]))
  }
}

scor.tab<-data.frame(matrix(ncol=ncol(vysl.all),nrow=nrow(vysl.all)))
colnames(scor.tab)<-c(paste0("N.",1:ncol(vysl.all)))
rownames(scor.tab)<-rownames(vysl.all)[1:nrow(vysl.all)]
for(gz in rownames(scor.tab)){
  for(ug in 1:ncol(scor.tab)){
    nnn<-as.numeric(unlist(strsplit(vysl.all[gz,ug],"/")))
    scor.tab[gz,ug]<-nnn[1]
  }
}
vysl.all$score<-vari.vysl$score<-rowMeans(scor.tab[,2:ncol(scor.tab)],na.rm=T)
vysl.all<-vysl.all[order(vysl.all$score,decreasing = T),]
vari.vysl<-vari.vysl[order(vari.vysl$score,decreasing = T),]

### the list of possible graphical output ###

# all checked scenarios with decreasing node support
pdf("possible_topology_scenarios.pdf",height=6,width=6.4)
par(mfrow=c(4,3),mar=c(1,0,1.5,1.5),oma=c(0,2,0,0))
for(sc in rownames(vysl.all)){
  tr.sc<-read.tree(paste0("tree_mod_",unlist(strsplit(sc,"[.]"))[2],".tre"))
  plot(tr.sc,cex=0.5,main=paste0(sc," / score:",round(vysl.all[sc,"score"],2)),cex.main=0.5)
  q<-sapply(vysl.all[sc,1:ncol(vysl.all)-1],function(x) unlist(strsplit(x,"/",fixed=T)))
  q[1,1]<-q[2,1]<-""
  nodelabels(round(as.numeric(q[1,]),0),cex=0.4,frame="none",col="red",adj=c(1.2,-0.5))
  nodelabels(round(as.numeric(q[2,]),0),cex=0.4,frame="none",col="red",adj=c(1.2,1.5))
}
dev.off()

# the same output but with variation in node support
pdf("possible_topology_scenarios-variability.pdf",height=6,width=6.4)
par(mfrow=c(4,3),mar=c(1,0,1.5,1.5),oma=c(0,2,0,0))
for(sc in rownames(vari.vysl)){
  tr.sc<-read.tree(paste0("tree_mod_",unlist(strsplit(sc,"[.]"))[2],".tre"))
  plot(tr.sc,cex=0.5,main=paste0(sc," / score:",round(vari.vysl[sc,"score"],2)),cex.main=0.5)
  q<-sapply(vari.vysl[sc,1:ncol(vari.vysl)-1],function(x) unlist(strsplit(x,"-",fixed=T)))
  q[1,1]<-q[2,1]<-""
  nodelabels(q[2,],cex=0.4,frame="none",col="red",adj=c(1.2,-0.5))
  nodelabels(q[1,],cex=0.4,frame="none",col="red",adj=c(1.2,1.5))
}
dev.off()

# the best twelfe scenarios with highest node support
pdf("selected_topology_scenarios.pdf",height=6,width=6.4)
par(mfrow=c(4,3),mar=c(1,0,1.5,1.5),oma=c(0,2,0,0))
for(sc in rownames(vysl.all)[1:12]){
  tr.sc<-read.tree(paste0("tree_mod_",unlist(strsplit(sc,"[.]"))[2],".tre"))
  plot(tr.sc,cex=0.5,main=paste0(sc," / score:",round(vysl.all[sc,"score"],2)),cex.main=0.5)
  q<-sapply(vysl.all[sc,1:ncol(vysl.all)-1],function(x) unlist(strsplit(x,"/",fixed=T)))
  q[1,1]<-q[2,1]<-""
  nodelabels(round(as.numeric(q[1,]),0),cex=0.4,frame="none",col="red",adj=c(1.2,-0.5))
  nodelabels(round(as.numeric(q[2,]),0),cex=0.4,frame="none",col="red",adj=c(1.2,1.5))
}
dev.off()
