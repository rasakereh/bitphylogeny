
# Run BitPhylogeny analysis-----------------------------------------------------
bitphyloR <- function(fin, fout, contains_true_label=FALSE, n=100, b=10, t=1,
                      mode = "methylation", rand_seed = 1234) {
    python.exec('from bitphylogeny.run import run_analysis')
    python.call('run_analysis', fin, fout, contains_true_label, n, b, t, 
                mode, rand_seed)
}
#-------------------------------------------------------------------------------

# Post process BitPhylogeny results --------------------------------------------
vmeasureR <- function(x, y){
    python.exec('from sklearn import metrics')
    return( python.call('metrics.homogeneity_completeness_v_measure',
                        x, y) )
}

# Compute v-measure
compute_vmeasures <- function(fout, fin){
    write_mpear_label(fout)
    python.exec('from bitphylogeny.run import write_vmeasures')
    python.call('write_vmeasures', fout,fin)
}

load_true_label <- function(ft){
  x <- read.csv(ft)
  true_label <- x[,dim(x)[2]]
}

write_vmeasures <- function(filepath, ft, fn){
  labels <- load_label_traces(filepath)
  true_label <- load_true_label(ft)
  vmeasures <- compute_vmeasures(labels, true_label)
  colnames(vmeasures) <- c('homogeneity','completeness', 'v_measure')
  write.csv(vmeasures, file = paste(filepath, fn, sep='/'), row.names=F)
}

load_vmeasures<- function( fp, key ){
  fn <- get_path(fp, key)
  x <- read.csv(fn, header=T, as.is = T)
  return(x)
}

# mpear label
compute_mpear_label <- function(label_traces){
  ltmat <- as.matrix(label_traces)
  class(ltmat) <- 'integer'
  ltmat <- ltmat + 1
  psm <- comp.psm(ltmat)
  mpear <- maxpear(psm)
  mpear_label <- mpear$cl
  return(mpear_label)
}

# wrappers
load_label_traces <- function(filepath){
  filename <- get_path(filepath, 'label_traces')
  label_traces <- read.csv(filename, header = F)
  return(label_traces)
}

write_mpear_label <- function(filepath){
  label_traces <- load_label_traces(filepath)
  mpear_label <- compute_mpear_label(label_traces)
  outfile <- 'mpear_label.csv'
  write.csv(mpear_label,
            file = paste(get_path(filepath,'mcmc-traces'),outfile,sep='/'),
            row.names = F)
}

#-------------------------------------------------------------------------------

# Baseline clustering analysis -------------------------------------------------

get_label_hc <- function(x, K){

  dis <- dist(x, 'binary')

  hc_cand <- lapply(K, function(ii) cutree( hclust( dis ), ii ) )

  hc_silhouette_res <- sapply(1:length(K),
                              function(ii)
                                summary( silhouette(hc_cand[[ii]] ,
                                                    dis) )$avg.width )
  idx <- which.max( hc_silhouette_res )

  hc_label <- hc_cand[[idx]]

  clone <- sapply(unique(hc_label), function(i) which(hc_label==i) )

  n <- length(clone)

  genotype <- matrix(0, n, dim(x)[2])

  for (i in 1:n){
    idx <- clone[[i]]
    if ( length(idx)==1 ){
      genotype[i,] = as.matrix(x[idx,])
    }else{
      genotype[i,] = as.numeric( colMeans(as.matrix(x[idx,])) > 0.5 )
    }
  }

  return(list(hc_label = hc_label, hc_genotype = genotype))
}

get_label_kc <- function( x, K ){

  dis <- dist(x, 'binary')

  kc_cand <- lapply(K, function(ii) pam( dis, ii) )

  kc_silhouette_res <- sapply(1:length(K),
                              function(ii)
                                summary( silhouette(kc_cand[[ii]]$clustering,
                                                    dis) )$avg.width )
  idx <- which.max( kc_silhouette_res )

  kc_label <- kc_cand[[idx]]$clustering
  kc_genotype <- x[kc_cand[[idx]]$medoids,]

  return(list(kc_label = kc_label, kc_genotype = kc_genotype))
}

#-------------------------------------------------------------------------------

# tree building ----------------------------------------------------------------
get_mst <- function(genotype){
  n <- dim(genotype)[1]
  adj <- hamming.distance(as.matrix(genotype))
  gr <- graph.adjacency( adjmatrix= adj, weighted= T, mode='undirected')
  mst <- minimum.spanning.tree(gr, algorithm= 'prim')
}

plot_mst <- function(genotype, label, flag= FALSE){
  reads <- sapply(unique(unlist( label ) ),
                    function(ii) length( which( label==ii )))
  mst <- get_mst(genotype)
  l = layout.reingold.tilford(graph=mst,
                              root=which.min(rowSums(genotype)))
  nodes = matrix(0, dim(genotype)[1], 3 )
  for (ii in 1:dim(nodes)[1]){
    geno = c()
    for (jj in 1:dim(genotype)[2]){
      geno <- paste( geno, genotype[ii,jj] )
    }
    nodes[ii,1] <- ii
    nodes[ii,2] <- reads[ii]
    nodes[ii,3] <- geno
  }
  colnames(nodes) = c('Nodes','Read Counts','Genotype')

  if (flag){
      pdf( paste(p2, '/', fns[i], '.pdf',sep=''),  width=15, height=9)
  }
  par(mfrow = c(1, 2), oma = c(3,3,0,0) + 0.1,
      mar = c(0,0,1,0.5) )
  plot( mst, layout=l )
  textplot(nodes, show.rownames = F, cex=1.25)

  if (flag){
      dev.off()
  }

}

plot_mst_from_dir <- function(p1,p2, flag=T){
  fns <- dir(p1, pattern='genotype')
  fns1 <- dir(p1, pattern='label')
  n <- length(fns1)
  for ( i in 1:n ){
    genotype <- read.csv(paste(p1,fns[i],sep='/'))
    label <- read.csv(paste(p1,fns1[i],sep='/'))
    plot_mst( genotype, label, flag )
  }
}

plot_sankey<-function(df){
  nodemat <- df$nodemat
  edgemat <- df$edgemat
  genomat <- df$genomat
  nn <- nodemat[,1]
  g <- graph.empty() + vertices(nn)
  ee <- unlist(t(edgemat[,1:2]))
  g <- g + edges(sapply(ee, function(i) which(nn==i)))
  ll <- layout.reingold.tilford(g)
  nodemat[,'y'] = ll[,1]
  river <- makeRiver(nodemat, edgemat)
  style <- list(edgecol= "col")
  par(mfrow = c(1, 2))
  riverplot(river, srt=90, lty=1, default_style=style)
  textplot(genomat, show.rownames = F)
}

plot_sankey_mft <- function(fh){
  treefreq <- read.csv(fh)
  mft <- treefreq[which.max(treefreq[,'freq']), 1]
  fn <- paste(dirname(fh), '/nodes-', mft, '.gdl', sep='')
  df <- gdl2df(fn)
  plot_sankey(df)
}

#-------------------------------------------------------------------------------

# Utilities --------------------------------------------------------------------

gdl2df <- function( file_gdl ){
  filetext <- read.csv(file_gdl, as.is  = T, quote="")
  filetext <- filetext[-nrow(filetext),]
  # merge node
  mergenode <- paste(filetext[1], filetext[2], sep='')
  for (ii in 3:length(filetext) ){
    tmp = unlist(strsplit(filetext[ii], ' '))
    if (tmp[1] == 'node:'){
      tmp1 <- paste(filetext[ii], filetext[ii+1], sep='')
      mergenode <- c(mergenode, tmp1)
    }
    if (tmp[1] == 'edge:'){
        mergenode <- c(mergenode, filetext[ii])
    }
  }
  filetext <- mergenode

  # head node
  ftsplit = unlist(strsplit(filetext[1], ' '))
  ID = unlist(strsplit(ftsplit[23],':'))[2]
  ID = unlist(strsplit(ID,'"'))[2]
  layer = length(ID)
  x = as.numeric(ftsplit[18])
  y = 0
  col = 0
  nodes = c(ID, x, y,col,layer)
  edges = c()

  for (ft in 2:length(filetext)) {

    node_sp = unlist(strsplit(filetext[ft], ' '))
    if (node_sp[1] == "node:") {
      # nodes
      ID = unlist(strsplit(node_sp[24],':'))[2]
      ID = unlist(strsplit(ID,'"'))[2]
      layer = length(unlist(strsplit(ID,'-')))
      x = as.numeric(node_sp[19])
      y = 0
      col = 0
      nodes = rbind( nodes,c(ID, x, y, col, layer) )

      # edges
      edge_sp = unlist(strsplit(filetext[ft+1], ' '))
      if (edge_sp[1] == "edge:") {
        ID1 = unlist(strsplit(edge_sp[3],':'))[2]
        ID1 = unlist(strsplit(ID1,'"'))[2]
        ID2 = unlist(strsplit(edge_sp[4],':'))[2]
        ID2 = unlist(strsplit(ID2,'"'))[2]
        node_label = unlist(strsplit(filetext[ft], ' '))
        value = unlist(strsplit(node_label[4], ':'))[2]
        value = unlist(strsplit(value, '"'))[2]
        direction = "A"
        edges = rbind (edges, c(ID1, ID2, value, direction))
      }

    }

  }

  colnames(nodes) = c("ID", "x", "y", "col", "layer")
  rownames(nodes) = NULL
  colnames(edges) = c("N1", "N2", "Value", "direction")

  # change branch length
  for (k in 1:nrow(nodes)) {
    code = unlist(strsplit(nodes[k,"ID"],'-'))

    if (length(code) >=3) {
      index = which(nodes[, "layer" ] == (length(code)-1))
      for (id in index) {
        code1 = unlist(strsplit(nodes[id,"ID"],'-'))
        tt = sapply(2:length(code1), function(i) code1[i] == code[i])
        if ( sum(tt) == (length(code1) - 1) ) {
          nodes[k, "x"] = as.numeric(nodes[k,"x"]) +
            as.numeric(nodes[id, "x"])
          break
        }
      }
    }
  }

  col_pa = c("coral2", "deepskyblue","lightgreen", "pink2", "black", 'yellow')
  for (i in unique(nodes[,"layer"]) ) {
    i = as.numeric(i)
    nodes[nodes[,"layer"] == as.character(i),"col"] = col_pa[i]
    if (i != 1) {
      no_i = sum(nodes[, "layer"] == i)
      lower = as.numeric(nodes[1,"y"]) - i*0.1*as.numeric(nodes[1,"y"])
      upper = as.numeric(nodes[1,"y"]) + i*0.1*as.numeric(nodes[1,"y"])
      nodes[nodes[, "layer"] == i, "y"] =  seq(lower,upper, length = no_i)
    }
  }
  nodes = nodes[, -ncol(nodes)]

  nodes = as.data.frame(nodes)
  nodes[,1] = as.character(nodes[,1])
  nodes[,2] = as.numeric(as.character(nodes[,2]))
  nodes[,3] = as.numeric(as.character(nodes[,3]))
  nodes[,4] = as.character(nodes[,4])
  edges = as.data.frame(edges)
  edges[,1] = as.character(edges[,1])
  edges[,2] = as.character(edges[,2])
  edges[,3] = as.numeric(as.character(edges[,3]))
  edges = cbind(edges, "gray90")
  edges[,5] = as.character(edges[,5])
  colnames(edges)[5] = "col"


  genomat <- c()
  for (i in 1:length(filetext)) {
    ftsplit = unlist(strsplit(filetext[i], ' '))
    index = which(ftsplit == "Genotype")
    if (length(index)!=0) {
      geno_single = c()
      for (k in (index+2):(index+9)) {
        geno_single = paste(geno_single,ftsplit[k])
      }
      genotype = geno_single
      reads= unlist(strsplit(ftsplit[index-2], '"'))[2]
      genomat <- rbind(genomat, c(reads, genotype))
    }
  }

  genomat <- cbind( nodes[,1], genomat)
  colnames(genomat) <- c('Nodes','Reads Count', 'Genotype')
  return(list(nodemat = nodes, edgemat = edges, genomat=genomat))
}

get_path <- function(folder, key){
  trace <- dir(folder, recursive=T, include.dirs =T, pattern = key)
  filepath <- paste(folder, trace, sep='/')
  return(filepath)
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}