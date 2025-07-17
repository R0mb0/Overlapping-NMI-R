# Implementation of Overlapping NMI library for R

[![Codacy Badge](https://app.codacy.com/project/badge/Grade/995b6a66d00f4901872d566d2ef1d4a6)](https://app.codacy.com/gh/R0mb0/Overlapping-NMI-R/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/R0mb0/Overlapping-NMI-R)
[![Open Source Love svg3](https://badges.frapsoft.com/os/v3/open-source.svg?v=103)](https://github.com/R0mb0/Overlapping-NMI-R)
![License](https://img.shields.io/badge/license-AGPLv3-blue.svg?style=plastic)

<details>
<summary>
   
## 👉 Read the original documentation

</summary>

An implementation of a Normalized Mutual Information (NMI) measure for sets of overlapping clusters.

Fully described in:
   "Normalized Mutual Information to evaluate overlapping community finding algorithms"
   by Aaron F. McDaid, Derek Greene, Neil Hurley
   http://arxiv.org/abs/1110.2515

Our method is based on the method described in Appendix B at the end of:
  "Detecting the overlapping and hierarchical community structure in complex networks"
  by Andrea Lancichinetti, Santo Fortunato and János Kertész
  http://iopscience.iop.org/1367-2630/11/3/033015/


== Usage ==

  make onmi
  onmi FILE1 FILE2

The filesnames record the sets of communities. A typical use case is to have
the "true" communities in one file and and those found by your algorithm
in the other file. One line per community. The nodes are
separated by whitespace, and any non-whitespace characters may be used in the
node names.


== Contact ==

Send any comments or queries or requests to aaronmcdaid@gmail.com


== Citation for Bibtex ==

@article{McDaidNMI,
    abstract = {Given the increasing popularity of algorithms for overlapping clustering, in
particular in social network analysis, quantitative measures are needed to
measure the accuracy of a method. Given a set of true clusters, and the set of
clusters found by an algorithm, these sets of clusters must be compared to see
how similar or different the sets are. A normalized measure is desirable in
many contexts, for example assigning a value of 0 where the two sets are
totally dissimilar, and 1 where they are identical. A measure based on
normalized mutual information, [1], has recently become popular. We demonstrate
unintuitive behaviour of this measure, and show how this can be corrected by
using a more conventional normalization. We compare the results to that of
other measures, such as the Omega index [2].},
    archivePrefix = {arXiv},
    author = {McDaid, Aaron F. and Greene, Derek and Hurley, Neil},
    citeulike-article-id = {9896732},
    citeulike-linkout-0 = {http://arxiv.org/abs/1110.2515},
    citeulike-linkout-1 = {http://arxiv.org/pdf/1110.2515},
    day = {11},
    eprint = {1110.2515},
    month = oct,
    posted-at = {2011-10-13 02:42:56},
    priority = {0},
    title = {Normalized Mutual Information to evaluate overlapping community finding algorithms},
    url = {http://arxiv.org/abs/1110.2515},
    year = {2011}
}

</details>

<details>
<summary>

## 📝 Exercise solved using the R library
   
</summary>

### Where the exercise has been taken 

The exercise is: 38.9.1/2 from this [book](https://www.networkatlas.eu/files/sna_book.pdf)  

### Exercise text 

Use the k-clique algorithm to find overlapping communities in the network at [http://www.networkatlas.eu/exercises/38/1/data.
txt](http://www.networkatlas.eu/exercises/38/1/data.txt). Test how many nodes are part of no community for k equal to 3, 4, and 5.

Compare the k-clique results with the coverage in [http://www.networkatlas.eu/exercises/38/2/comms.txt](http://www.networkatlas.eu/exercises/38/2/comms.txt), by using any variation of overlapping NMI from [https://github.com/aaronmcdaid/Overlapping-NMI](https://github.com/aaronmcdaid/Overlapping-NMI). For which value of k do you get the best performance?

### Library 

```R
# Calculating the Overlapping Normalized Mutual Information (NMI) between two covers
# Robust and compatible with lists of character vectors (each a community)

get_membership_matrix <- function(communities, all_nodes) {
  mat <- matrix(0, nrow=length(all_nodes), ncol=length(communities))
  rownames(mat) <- all_nodes
  for (j in seq_along(communities)) {
    idx <- match(communities[[j]], all_nodes)
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0) {
      mat[idx, j] <- 1
    }
  }
  mat
}

NMI <- function(cover1, cover2) {
  all_nodes <- sort(unique(c(unlist(cover1), unlist(cover2))))
  X <- get_membership_matrix(cover1, all_nodes)
  Y <- get_membership_matrix(cover2, all_nodes)
  n <- length(all_nodes)
  safe_log2 <- function(x) ifelse(x > 0, log2(x), 0)
  cond_entropy <- function(A, B) {
    kA <- ncol(A)
    kB <- ncol(B)
    H <- 0
    for (i in 1:kA) {
      minH <- Inf
      for (j in 1:kB) {
        Nij <- sum(A[,i] & B[,j])
        if (Nij == 0) next
        Ni <- sum(A[,i])
        Nj <- sum(B[,j])
        pij <- Nij / n
        pi <- Ni / n
        pj <- Nj / n
        Hij <- 0
        if (pij > 0 && (pi * pj) > 0)
          Hij <- Hij - (pij) * safe_log2(pij / (pi * pj))
        if ((Ni-Nij) > 0 && (pi * (1-pj)) > 0)
          Hij <- Hij - ((Ni-Nij)/n) * safe_log2(((Ni-Nij)/n) / (pi*(1-pj)))
        if ((Nj-Nij) > 0 && ((1-pi)*pj) > 0)
          Hij <- Hij - ((Nj-Nij)/n) * safe_log2(((Nj-Nij)/n) / ((1-pi)*pj))
        if ((n-Ni-Nj+Nij) > 0 && ((1-pi)*(1-pj)) > 0)
          Hij <- Hij - ((n-Ni-Nj+Nij)/n) * safe_log2(((n-Ni-Nj+Nij)/n) / ((1-pi)*(1-pj)))
        if (!is.nan(Hij) && Hij < minH)
          minH <- Hij
      }
      if (is.finite(minH)) H <- H + minH
    }
    H / kA
  }
  H_XY <- cond_entropy(X, Y)
  H_YX <- cond_entropy(Y, X)
  NMI_value <- 1 - 0.5 * (H_XY + H_YX)
  NMI_value
}
```

### Exercise solution

```R
# Compare the k-clique results with the coverage in http://www.
# networkatlas.eu/exercises/38/2/comms.txt, by using any variation
# of overlapping NMI from https://github.com/aaronmcdaid/
# Overlapping-NMI. For which value of k do you get the best performance?

library(here)
library(igraph)

# Defining the k_clique_communities function for finding k-clique communities
k_clique_communities <- function(graph, k) {
  all_cliques <- cliques(graph, min = k, max = k)
  if (length(all_cliques) == 0) return(list())
  clique_graph <- make_empty_graph(n = length(all_cliques))
  for (i in seq_along(all_cliques)) {
    for (j in seq_len(i-1)) {
      if (length(intersect(all_cliques[[i]], all_cliques[[j]])) == (k-1)) {
        clique_graph <- add_edges(clique_graph, c(i, j))
      }
    }
  }
  comps <- components(clique_graph)
  communities <- lapply(seq_len(comps$no), function(comp_id) {
    idx <- which(comps$membership == comp_id)
    unique(unlist(all_cliques[idx]))
  })
  lapply(communities, function(x) V(graph)$name[x])
}

source(here("OverlappingNMI.R"))

# Loading the edge list and building the graph
edges <- read.table(here("data.txt"))
colnames(edges) <- c("from", "to", "weight")
g <- graph_from_data_frame(edges, directed=FALSE)

# Solution 

# Removing weights for clique percolation
g_unweighted <- delete_edge_attr(g, "weight")

# Loading the ground truth communities from comms.txt
comms_lines <- readLines(here("comms.txt"))
gt_communities <- lapply(comms_lines, function(line) strsplit(line, " ")[[1]])
gt_communities <- lapply(gt_communities, as.character)

# Defining a function for extracting k-clique communities as lists of character vectors
get_kc_comms <- function(graph, k) {
  kc <- k_clique_communities(graph, k)
  lapply(kc, as.character)
}

# Initializing a results data frame
results <- data.frame(k=integer(), nmi=numeric())

# Calculating the NMI for k = 3, 4, 5
for (k in 3:5) {
  cat(sprintf("Calculating for k = %d\n", k))
  detected <- get_kc_comms(g_unweighted, k)
  nmi_value <- NMI(gt_communities, detected)
  results <- rbind(results, data.frame(k=k, nmi=nmi_value))
  cat(sprintf("k = %d | NMI = %.4f\n", k, nmi_value))
}

cat("\nSummary of k-clique NMI results:\n")
print(results, row.names = FALSE)

best_k <- results$k[which.max(results$nmi)]
cat(sprintf("\nBest performance for k = %d (NMI = %.4f)\n", best_k, max(results$nmi)))
```
### Referenced repository 

[R0mb0/The_Atlas_for_the_Aspiring_Network_Scientist_exercises_in_R](https://github.com/R0mb0/The_Atlas_for_the_Aspiring_Network_Scientist_exercises_in_R)

</details>

## Library interface 

```R
# Calculating the Overlapping Normalized Mutual Information (NMI) between two covers
# Robust and compatible with lists of character vectors (each a community)

NMI <- function(cover1, cover2)
```
