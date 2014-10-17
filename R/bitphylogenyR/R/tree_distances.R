df2igraph <- function(df) {
  nodemat <- df$nodemat
  edgemat <- df$edgemat
  nn <- nodemat[,1]
  g <- graph.empty() + vertices(nn)
  ee <- unlist(t(edgemat[,1:2]))
  g <- g + edges(sapply(ee, function(i) which(nn==i)))
  return(g)
}

get_genotype_intercept <- function(g1, g2) {
  
  tt = vector("list", dim(g1)[1])
  for (i in 1:dim(g1)[1]) {
    for (j in 1:dim(g2)[1]) {
      if (sum(g1[i, ] == g2[j, ]) == dim(g1)[2]) {
        tt[[i]] = c(tt[[i]],j)
      }
    }
  }
  
  ind1 <- c()
  for (i in 1:length(tt)) {
    if (!is.null(tt[[i]])) {
      ind1 <- c(ind1, i)
    }
  }
  
  ind2 <- c()
  for (i in 1:length(tt)) {
    if (!is.null(tt[[i]])) {
      ind2 <- c(ind2, tt[[i]][1])
    }
  }
  return(list(ind1 = ind1, ind2 = ind2))
}

get_genotype <- function(tree2) {
  if (is.null(tree2$nodemat)) {
    g2 <- get_mst(tree2)
    t2 <- as.matrix(tree2)
  } else {
    g2 <- df2igraph(tree2)
    genomat <- tree2$genomat 
    t2 <- unlist(strsplit(genomat[1,3], " "))
    for (i in 2:dim(genomat)[1]){
      t2 <- rbind(t2, unlist(strsplit(genomat[i,3], " ")))
    }
    t2 = t2[,-1]
    rownames(t2) <- NULL
    class(t2) <- "numeric"
    t2 <- t2 > 0.5
  }
  class(t2) <- "numeric"
  
  return(list(t = t2, g=g2))
}

get_tree_distance_pairwise <- function(tree1, tree2) {
  
  tmp <- get_genotype(tree1)
  t1 <- tmp$t
  g1 <- tmp$g
  
  tmp <- get_genotype(tree2)
  t2 <- tmp$t
  g2 <- tmp$g
  
  ind <- get_genotype_intercept(t1, t2)
  
  ind1 <- ind$ind1
  ind2 <- ind$ind2
  
  mat1 <- shortest.paths(g1)
  mat2 <- shortest.paths(g2)
  
  submat1 <- mat1[ind1, ind1]
  submat2 <- mat2[ind2, ind2]
  
  m = sum(abs(submat1[lower.tri(submat1)] - submat2[lower.tri(submat2)]))
  
  return(list(m = m, n = length(ind1)))
} 

get_tree_distance_common <- function(tree1, tree2, tree3, tree4) {
  
  tmp <- get_genotype(tree1)
  t1 <- tmp$t
  g1 <- tmp$g
  
  tmp <- get_genotype(tree2)
  t2 <- tmp$t
  g2 <- tmp$g
  
  tmp <- get_genotype(tree3)
  t3 <- tmp$t
  g3 <- tmp$g
  
  tmp <- get_genotype(tree4)
  t4 <- tmp$t
  g4 <- tmp$g
  
  ind1 <- try(get_genotype_intercept(t1, t2))
  ind2 <- try(get_genotype_intercept(t1, t3))
  ind3 <- try(get_genotype_intercept(t1, t4))
  
  ind <- Reduce(intersect, list(ind1$ind1, ind2$ind1, ind3$ind1))
  
  tt <- t1[ind, ]
  
  ind1 <- get_genotype_intercept(tt, t2)
  ind2 <- get_genotype_intercept(tt, t3)
  ind3 <- get_genotype_intercept(tt, t4)
  
  mat1 <- shortest.paths(g1)
  mat2 <- shortest.paths(g2)
  mat3 <- shortest.paths(g3)
  mat4 <- shortest.paths(g4)
  
  submat1 <- mat1[ind, ind]
  submat2 <- mat2[ind1$ind2, ind1$ind2]
  submat3 <- mat3[ind2$ind2, ind2$ind2]
  submat4 <- mat4[ind3$ind2, ind3$ind2]
  
  m1 <- sum(abs(submat1[lower.tri(submat1)] - submat2[lower.tri(submat2)]))
  m2 <- sum(abs(submat1[lower.tri(submat1)] - submat3[lower.tri(submat3)]))
  m3 <- sum(abs(submat1[lower.tri(submat1)] - submat4[lower.tri(submat4)]))
  
  return(list(m = c(m1, m2, m3), n = length(ind)))
} 

get_pairwise_tree_distance_df <- function(f1, f2, f3, treeType, f4){
  
  folders <- dir(f1, pattern = "mutmat", 
                 full.names = TRUE, recursive = TRUE, 
                 include.dirs = TRUE)
  if (treeType == "star-clone") {
    folders <- folders[c(2,4,6,8)]
  }
  
  folders = folders[c(4,1,2,3)] 
  
  pathToSynth <- "extdata/synthetic/"
  
  treeType <- treeType
  
  pathToFiles <- system.file(paste(pathToSynth, treeType, sep = ""), 
                             package = "bitphylogenyR")
  
  files <- dir(pathToFiles, pattern = "mutmat", full.names = T)
  
  true <- gdl2df(file_gdl = f4)
  hc_genotype_files <- dir(path = f2,
                           pattern = "hc_genotype", 
                           recursive = T, full.names = T)
  kc_genotype_files <- dir(path = f3,
                           pattern = "kc_genotype", 
                           recursive = T, full.names = T)
  errors = c("0.01", "0.02", "0.05", "0")
  
  bitphylo_trees = vector("list", length = 4)
  for (i in 1:4) {
    fh <- paste(folders[i], "/treescripts/", sep = "")
    treefreq <- read.csv(dir(path = fh, pattern = "tree-freq", recursive = T,
                             full.names = T))
    mft <- treefreq[which.max(treefreq[, "freq"]), 1]
    fn <- paste(fh, "/nodes-", mft, ".gdl", sep = "")
    bitphylo_trees[[i]] <- gdl2df(fn)
  }
  
  bit_distance = sapply(1:4, function(i)
    get_tree_distance_pairwise(true, bitphylo_trees[[i]]))
  
  kc_distance = sapply(1:4, function(i)
    get_tree_distance_pairwise(true, read.csv(kc_genotype_files[i]))) 
  
  hc_distance = sapply(1:4, function(i)
    get_tree_distance_pairwise(true, read.csv(hc_genotype_files[i])))  
  
  metrics <- data.frame(distance = c(bit_distance$m, kc_distance, hc_distance), 
                        method = c(rep("BitPhylogeny", 4), 
                                   rep("k-centroids", 4),
                                   rep("hierarchical clustering", 4)),
                        tree_type = treeType, 
                        error = errors)
  return(metrics)
  
}

get_common_tree_distance_df <- function(f1, f2, f3, treeType, f4){
  
  folders <- dir(f1, pattern = "mutmat", 
                 full.names = TRUE, recursive = TRUE, 
                 include.dirs = TRUE)
  if (treeType == "star-clone") {
    folders <- folders[c(2,4,6,8)]
  }
  
  folders = folders[c(4,1,2,3)] 
  
  pathToSynth <- "extdata/synthetic/"
  
  treeType <- treeType
  
  pathToFiles <- system.file(paste(pathToSynth, treeType, sep = ""), 
                             package = "bitphylogenyR")
  
  files <- dir(pathToFiles, pattern = "mutmat", full.names = T)
  
  true <- gdl2df(file_gdl = f4)
  hc_genotype_files <- dir(path = f2,
                           pattern = "hc_genotype", 
                           recursive = T, full.names = T)
  kc_genotype_files <- dir(path = f3,
                           pattern = "kc_genotype", 
                           recursive = T, full.names = T)
  errors = c("0", "0.01", "0.02", "0.05")
  
  bitphylo_trees = vector("list", length = 4)
  for (i in 1:4) {
    fh <- paste(folders[i], "/treescripts/", sep = "")
    treefreq <- read.csv(dir(path = fh, pattern = "tree-freq", recursive = T,
                             full.names = T))
    mft <- treefreq[which.max(treefreq[, "freq"]), 1]
    fn <- paste(fh, "/nodes-", mft, ".gdl", sep = "")
    bitphylo_trees[[i]] <- gdl2df(fn)
  }
  
  distance = vector("numeric", length = 12)
  matched_clone = distance
  for (i in 1:4) {
    res <- get_tree_distance_common(true, 
                                    bitphylo_trees[[i]],
                                    read.csv(kc_genotype_files[i]),
                                    read.csv(hc_genotype_files[i]))
    distance[((i-1)*3+1):(i*3)] = res$m
    matched_clone[((i-1)*3+1):(i*3)] = res$n
  }
  
  metrics <- data.frame(distance = distance, 
                        method = c("BitPhylogeny", 
                                   "k-centroids", 
                                   "hierarchical clustering"),
                        tree_type = treeType, 
                        matched_clone = matched_clone,
                        error = c(rep("0", 3),
                                  rep("0.01", 3),
                                  rep("0.02", 3),
                                  rep("0.05", 3))
  )
  return(metrics)
  
}