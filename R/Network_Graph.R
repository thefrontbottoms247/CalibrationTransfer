#' @title Short title
#' @description What the function does
#' @param arg Description
#' @return What it outputs
#' @export
Network_Graph <- function(){
  
  g <- cm.slave_cor[,,1]
reduced_graph <- function(graph){
  ng <- c()
  for (class in 1:sqrt(length(graph))){
    if (sum(graph[class,],graph[,class])/2 > graph[class,class]){
      ng <- c(ng, class)}}
  g <<- graph[ng,ng]}
graph <- graph_from_adjacency_matrix(g, mode = 'directed', weighted = T, diag = F)
oweights <- t(g*as.matrix(as_adjacency_matrix((graph))))
oweights <- as.vector(t(oweights))[as.vector(t(oweights))!=0]
edge_w <- rescale(oweights, to = c(0.5, 5))
vsize <- rescale(degree(graph), to = c(8, 25))

lay <- layout_with_fr(graph)
lay <- layout_with_fr(graph, niter = 1000, grid = "nogrid")
lay <- layout_with_kk(graph)
lay <- layout_with_drl(graph, options = list(simmer.attraction = 0))
lay <- layout_in_circle(graph)

plot(graph,                                   # plot the igraph object
     layout = lay,                            # use precomputed layout for node positions
     vertex.size = vsize,                     # vertex size scaled by node degree
     vertex.shape = "circle",                 # simple circular vertices for clarity
     vertex.color = "navy",                   # fill color of vertices
     vertex.frame.color = "white",            # border color around each vertex
     vertex.frame.width = 1.5,                # border thickness
     vertex.label = colnames(g),              # label vertices with column names
     vertex.label.color = "white",            # label text color
     vertex.label.family = "mono",            # monospace font for alignment
     vertex.label.cex = 0.8,                  # label size
     vertex.label.dist = 0.5,                 # small offset from vertex center
     vertex.label.degree = 0,                 # label angle (0 = outward)
     edge.width = edge_w,                     # line thickness scaled by weight
     edge.color = alpha("maroon", 0.5),       # semi-transparent edge color
     edge.lty = 1,                            # solid line type
     edge.curved = 0.2,                       # slight curvature for readability
     edge.arrow.size = 1.5,                   # arrowhead size
     edge.arrow.width = 0.3,                  # arrowhead width
     edge.label = ifelse(edge_w > 3, oweights, NA),  # show labels only for strong edges
     edge.label.color = "black",              # edge label color
     edge.label.family = "mono",              # monospace font for edge labels
     edge.label.cex = 0.7,                    # edge label size
     edge.loop.angle = 0)                     # disable self-loop angle distortion
g[g <= 50] <- 0
graph <- graph_from_adjacency_matrix(g, mode = "directed", weighted = TRUE, diag = FALSE)

}