library(tidyverse)
library("CliquePercolation")
library("bootnet")


#load data
data(Obama)

#estimate network
net <- bootnet::estimateNetwork(Obama, default = "EBICglasso", missing = "pairwise")

#plot network
graph <- plot(net, layout = "spring")

## Run CP algorithm

#run clique percolation algorithm with optimal k and I
cpk3I.14 <- cpAlgorithm(graph, k = 3, I = 0.14, method = "weighted")

#print results overview
cpk3I.14

# Save list of assigned communities (numbers)
community_assignments_list <- cpk3I.14$list.of.communities.numbers

# Use lapply to associate each sublist with its parent index
CP_community_assignment_df <- do.call(rbind, lapply(seq_along(community_assignments_list), function(i) {
  data.frame(
    Module = i,
    Node_INdex = unlist(community_assignments_list[[i]])
  )
}))
