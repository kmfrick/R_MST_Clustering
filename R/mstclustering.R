#' Generate edge list from a distance matrix
#' Duplicates are not deleted, because they will not be counted
#' by Kruskal's algorithm
#' If a check is O(1), this only adds an O(E) overhead, which is negligible
#' @param mat The distance matrix.
#' @return A data frame with three columns: 'from', 'to', 'dist'.
#' @import reshape2
#' @export
gen.edge.list <- function(mat) {
  mat <- as.matrix(mat)
  edge.list <- melt(mat)
  colnames(edge.list) <- c('from', 'to', 'dist')
  edge.list <- edge.list[edge.list$from != edge.list$to,]
  edge.list
}

#' UFDS
#' @param p parent of an item, initially the item itself
#' @param rank upper bound to tree depth excluding root (initially 0)
#' @export
setRefClass("ufds",
            fields=list(
              rank="vector",
              p="vector"
            )
)

#' Initialize UFDS
#' @param m Number of elements.
#' @return An instance of the class ufds containing m disjoint sets.
#' @import methods
#' @export
reset.ufds <- function(m) {
  ufds <- new("ufds", rank=rep(0,m), p=as.numeric(1:m))
  ufds
}

#' Find the set an elements belongs to.
#' @param i The element to check.
#' @param ufds An instance of the class ufds.
#' @export
find.set <- function (i, ufds) {
  if (ufds$p[i] == i) {
    return(i)
  } else {
    # Optimization: path compression
    return(ufds$p[i] <- find.set(ufds$p[i], ufds))
  }
}

#' Check if two elements are in the same set
#' @param i The first element in the tuple.
#' @param j The second element in the tuple.
#' @param ufds An instance of the class ufds.
#' @return TRUE if the elements are in the same set, FALSE otherwise.
#' @export
is.same.set <- function(i, j, ufds) {
  set.i <- find.set(i, ufds)
  set.j <- find.set(j, ufds)
  return(set.i == set.j)
}


#' Join the sets the two elements passed as arguments belong to.
#' @param i The first element in the tuple.
#' @param j The second element in the tuple.
#' @param ufds An instance of the class ufds.
#' @export
union.set <- function(i, j, ufds) {
  if (!is.same.set(i, j, ufds)) {
    set.i <- find.set(i, ufds)
    set.j <- find.set(j, ufds)
    # Optimization: use rank to keep the tree short
    if (ufds$rank[set.i] > ufds$rank[set.j]) {
      ufds$p[set.i] <- set.j
    } else {
      ufds$p[set.j] <- set.i
    }
    if (ufds$rank[set.i] == ufds$rank[set.j]) {
      ufds$rank[set.j] <- ufds$rank[set.j] + 1
    }
  }
}

#' Kruskal's algorithm
#' @param edge.list A data frame with columnns 'from', 'to', 'dist'.
#' @param m Number of nodes.
#' @return A list of edges in the MST, in the same format as the input argument edge.list.
#' @export
kruskal <- function(edge.list, m) {
  mst.cost <- 0
  ufds <- reset.ufds(m)
  edge.list.ordered <- edge.list[order(edge.list$dist), ]
  mst.edge.list <- data.frame(row.names = names(edge.list))
  for (i in 1:nrow(edge.list)) {
    u <- edge.list.ordered[i,]
    if (!is.same.set(u$to, u$from, ufds)) {
      mst.cost <- mst.cost + u$dist
      union.set(u$to, u$from, ufds)
      mst.edge.list <- rbind(mst.edge.list, u)
    }
  }
  cat("MST cost = ", mst.cost)
  mst.edge.list
}

#' Generate an adjacency list
#' @param clust.edge.list The return value of the kruskal() function.
#' @param m Number of nodes.
#' @return An adjacency list in the form of a list of vectors.
#' @export
gen.child.list.mst <- function(clust.edge.list, m) {
  child.list.mst <- vector("list", m)
  for (i in 1:m) {
    child.list.mst[[m]] <- vector()
  }
  for (i in 1:nrow(clust.edge.list)) {
    r <- clust.edge.list[i,]
    child.list.mst[[r$from]] <- c(child.list.mst[[r$from]] , r$to)
    child.list.mst[[r$to]] <- c(child.list.mst[[r$to]] , r$from)
  }
  child.list.mst
}

#' Connected components DFS
#' @export
setRefClass("ClusterList",
            fields=list(
              visited="vector",
              clust.mst="vector"
            )
)

#' Run clustering using MST.
#' Before calling this function, remove some edges from the MST, for example the  k-1 heaviest.
#' For example:
#' ```
#' k <- 3
#' n.edges <- nrow(mst.edge.list)
#' child.list.mst <- child.list.mst[1:(n.edges - (k - 1)),]
#' ```
#' @param child.list.mst The return value of the gen.child.list.mst() function with k-1 edges removed.
#' @param m Number of nodes.
#' @param k The number of clusters.
#' @return A vector whose k-th element is the cluster the k-th point belongs to.
#' @import methods
#' @export
mst.cluster <- function(child.list.mst, m, k) {
  cluster.list.top <- new("ClusterList", visited=rep(FALSE, m), clust.mst=rep(0, m))
  dfs <- function(u, clust, cluster.list=cluster.list.top) {
    if (cluster.list$visited[u]) {
      return()
    } else {
      cluster.list$visited[u] <- TRUE
      cluster.list$clust.mst[u] <- clust
      for (v in child.list.mst[[u]]) {
        dfs(v, clust)
      }
    }
  }
  for (clust in 1:k) {
    dfs.started <- FALSE
    for (u in 1:m) {
      if (!cluster.list.top$visited[u]) {
        dfs.started <- TRUE
        dfs(u, clust)
      }
      if (dfs.started) {
        break
      }
    }
  }
  cluster.list.top$clust.mst
}
