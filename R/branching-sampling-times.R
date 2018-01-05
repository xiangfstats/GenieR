
branching.sampling.times <- function(phy){
  phy = new2old.phylo(phy)
  if (class(phy) != "phylo")
    stop("object \"phy\" is not of class \"phylo\"")
  tmp <- as.numeric(phy$edge)
  nb.tip <- max(tmp)
  nb.node <- -min(tmp)
  xx <- as.numeric(rep(NA, nb.tip + nb.node))
  names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
  xx["-1"] <- 0
  for (i in 2:length(xx)) {
    nod <- names(xx[i])
    ind <- which(phy$edge[, 2] == nod)
    base <- phy$edge[ind, 1]
    xx[i] <- xx[base] + phy$edge.length[ind]
  }
  depth <- max(xx)
  branching.sampling.times <- depth - xx
  return(branching.sampling.times)
}
