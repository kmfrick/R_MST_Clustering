#' @title gen.edge.list
#' @description Generate edge list from a distance matrix
#' Duplicates are not deleted, because they will not be counted
#' by Kruskal's algorithm
#' If a check is O(1), this only adds an O(E) overhead, which is negligible
#' @param mat The distance matrix.
#' @return A data frame with three columns: 'from', 'to', 'dist'.
#' @importFrom reshape2 melt
#' @export
gen.edge.list <- function(mat) {
  mat <- as.matrix(mat)
  edge.list <- reshape2::melt(mat)
  colnames(edge.list) <- c('from', 'to', 'dist')
  edge.list <- edge.list[edge.list$from != edge.list$to,]
  edge.list
}


#' @title reset.ufds
#' @description Initialize UFDS
#' @param m Number of elements.
#' @return A data table containing a 'rank' column and a 'parent' column.
#' @importFrom data.table data.table
#' @export
reset.ufds <- function(m) {
  ufds <- data.table(rank=rep(0,m), p=as.integer(1:m))
  ufds
}

#' @title find.set
#' @description Find the set an element belongs to.
#' @param i The element to check.
#' @param ufds A data.table representing a UFDS.
#' @return An integer: the root node of the set the element belongs to.
#' @importFrom data.table set
#' @export
find.set <- function (i, ufds) {
  if (ufds$p[i] == i) {
    return(i)
  } else {
    # Optimization: path compression
    set(ufds, i, "p", find.set(ufds$p[i], ufds))
    return(ufds$p[i])
  }
}

#' @title is.same.set
#' @description Check if two elements are in the same set
#' @param i The first element in the tuple.
#' @param j The second element in the tuple.
#' @param ufds A data.table representing a UFDS.
#' @return TRUE if the elements are in the same set, FALSE otherwise.
#' @export
is.same.set <- function(i, j, ufds) {
  set.i <- find.set(i, ufds)
  set.j <- find.set(j, ufds)
  return(set.i == set.j)
}


#' @title union.set
#' @description Join the sets the two elements passed as arguments belong to.
#' @param i The first element in the tuple.
#' @param j The second element in the tuple.
#' @param ufds A data.table representing a UFDS.
#' @return No return value, called for side effects on rank and p.
#' @importFrom data.table set
#' @export
union.set <- function(i, j, ufds) {
  if (!is.same.set(i, j, ufds)) {
    set.i <- find.set(i, ufds)
    set.j <- find.set(j, ufds)
    # Optimization: use rank to keep the tree short
    if (ufds$rank[set.i] > ufds$rank[set.j]) {
      set(ufds, set.i, "p", set.j)
    } else {
      set(ufds, set.j, "p", set.i)
    }
    if (ufds$rank[set.i] == ufds$rank[set.j]) {
      set(ufds, set.j, "rank", ufds$rank[set.j] + 1)
    }
  }
}

#' @title kruskal
#' @description Kruskal's algorithm for MST computation.
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
  cat("MST cost = ", mst.cost, "\n")
  mst.edge.list
}

#' @title gen.child.list.mst
#' @description Generate an adjacency list
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

#' @title mst.cluster
#' @description Run clustering using MST.
#' Before calling this function, remove some edges from the MST, for example the  k-1 heaviest.
#' @param child.list.mst The return value of the gen.child.list.mst() function with k-1 edges removed.
#' @param m Number of nodes.
#' @param k The number of clusters.
#' @return A vector whose k-th element is the cluster the k-th point belongs to.
#' @importFrom data.table set
#' @importFrom data.table data.table
#' @examples
#' iris.clean <- iris[,-5]
#' iris.dist <- as.matrix(dist(iris.clean))
#' iris.edge.list <- gen.edge.list(iris.dist)
#' m <- nrow(iris.dist)
#' iris.mst.edge.list <- kruskal(iris.edge.list, m)
#' k <- 3
#' n.edges <- nrow(iris.mst.edge.list)
#' iris.mst.edge.list <- iris.mst.edge.list[1:(n.edges - (k - 1)),]
#' iris.child.list.mst <- gen.child.list.mst(iris.mst.edge.list, m)
#' iris.clust.mst <- mst.cluster(iris.child.list.mst, m, k)
#' @export
mst.cluster <- function(child.list.mst, m, k) {
  cluster.list.top <- data.table(visited=rep(FALSE, m), clust.mst=rep(0, m))
  dfs <- function(u, clust, cluster.list=cluster.list.top) {
    if (cluster.list$visited[u]) {
      return()
    } else {
      set(cluster.list, u, "visited", TRUE)
      set(cluster.list, u, "clust.mst", clust)
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
