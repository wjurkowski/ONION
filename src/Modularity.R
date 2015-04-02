options(digits = 3, quote = F, scientific = F, stringsAsFactors = F)
library(igraph)
args <- commandArgs(trailingOnly = TRUE)


#read data
t = read.table(args[1], header = F, sep = "\t")
g = graph.edgelist(as.matrix(t[,1:2]))

#get clusters
cl<-clusters(g)
clsID<- which(cl$csize > as.numeric(args[3]))
#create subgraph out of all large clusters
#lclu<-which(cl$membership %in% clsID)
#create subgraph for each cluster separately
results<-data.frame()
for(c in clsID){
  #get list of vertices
  lclu<-which(cl$membership == c)

  # create a sub-network composed by ONLY the nodes in subv and the edges between them
  sg <- induced.subgraph(graph=g,vids=lclu)
  n=V(sg)$name
  el<-get.edgelist(sg)
  nam<-paste(c,"subnetwork.txt",sep="-")
  write.table(el,nam,col.names=F,quote=F,row.names=F,sep="\t")
  # create a sub-network composed by the nodes in subv and, if some of them is connected to other nodes (even if not in subv), 
  # take also them (and of course include all the edges among this bunch of nodes)
  #dg <- decompose.graph(g,mode="weak")
  #neighverts <- unique(unlist(sapply(sg1,FUN=function(s){if(any(V(s)$name %in% lclu)) V(s)$name else NULL})))
  #sgp <- induced.subgraph(graph=g,vids=neighverts)

  #for given modularization type
  if(args[2] == 'fgc'){
  #fast gc (Newman-Girvan algorithm)
    ug = as.undirected(sg,"collapse",edge.attr.comb = getIgraphOpt("edge.attr.comb"))
    sug<-simplify(ug)
    fgc = fastgreedy.community(sug, merges=TRUE, modularity=TRUE)
    maxid = which(fgc$modularity == max(fgc$modularity))
    cat("Cluster:",c,"; Modularization:",maxid,"; modularity = ", max(fgc$modularity),"\n")
    ctm = community.to.membership(sug, fgc$merges, maxid)
    d=data.frame(c,lclu,n,fgc$membership)
    names(d)=c("Cluster","VertexID","VertexName","ModuleMembership")
    results<-rbind(results,d)

    ##VISUALIZE
    #lfg = layout.fruchterman.reingold(g)
    #V(g)$color = colors()[ctm$membership+2]
    #plot.igraph(g, layout = lfg)
    #print(ctm)
  } else if(args[2] == 'ebc'){
  #edge betweenness community
    ebc=edge.betweenness.community (sug, directed = TRUE, edge.betweenness = TRUE, merges = TRUE, bridges = TRUE, labels = TRUE)
    maxid = which(ebc$modularity == max(ebc$modularity))
    cat("Cluster:",c,"; Modularization:",maxid,"; modularity = ", max(ebc$modularity),"\n")
    com=community.to.membership(sug, ebc$merges, maxid, membership=TRUE, csize=TRUE)
    d=data.frame(c,lclu,n,com$membership)
    names(d)=c("Cluster","VertexID","VertexName","ModuleMembership")
    results<-rbind(results,d)
    #sapply(neigh,function(x) sum(V(g)[x[-1]]$color==V(g)[x[1]])/lenght(x-1))
    #V(g)$color
  }
}

#get vertices neighborhood 
ne=neighborhood(g,1)
te<-paste(V(g)$name,unlist(lapply(ne, paste, collapse=" ")),sep="\t")
writeLines(te,"Neighbors.txt")
#writeLines(unlist(lapply(ne, paste, collapse=" ")),"Neighbors.txt")
#save membership table
write.table(results, "ClusterMemberships.txt", row.names=F, col.names=T, quote=F, sep="\t")

#wygeneruj siec
#g = barabasi.game(50)

