fixtree = read.tree("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/RAxML_unconstrained/RAxML_bestTree.unconstrained")
fixtree2 = chronos(fixtree, lambda=0.01)
fixtree3 = fixtree2
fixtree3$edge.length<-
  fixtree3$edge.length/max(nodeHeights(fixtree3)[,2])*35
write.tree(fixtree3, "~/Desktop/starting_tree.tre")