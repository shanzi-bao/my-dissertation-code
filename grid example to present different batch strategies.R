# Create a 3x3 grid visualization with distances marked and calculate distance-based batch pairs

# Create the grid points
grid_points <- expand.grid(x = 0:2, y = 0:2)
n_points <- nrow(grid_points)


grid_points$id <- 1:n_points

# Calculate all pairwise distances
distances <- matrix(NA, n_points, n_points)
for (i in 1:n_points) {
  for (j in 1:n_points) {
    distances[i, j] <- sqrt((grid_points$x[i] - grid_points$x[j])^2 + 
                              (grid_points$y[i] - grid_points$y[j])^2)
  }
}

# Find unique distances (excluding 0)
unique_distances <- sort(unique(round(distances[upper.tri(distances)], 3)))
unique_distances

# Create a list to store pairs for each distance batch
distance_batches <- list()
for (d in unique_distances) {
  pairs <- which(round(distances, 3) == d, arr.ind = TRUE)
  # Ensure we don't count (i,j) and (j,i) as separate pairs
  pairs <- pairs[pairs[,1] < pairs[,2], ]
  distance_batches[[as.character(d)]] <- pairs
}


batch_counts <- sapply(distance_batches, nrow)
batch_counts


pdf("grid_visualization.pdf", width = 10, height = 7)
par(mar = c(2, 2, 2, 2))

# Plot 1: Grid with point IDs
plot(grid_points$x, grid_points$y, type = "n", xlim = c(-0.5, 2.5), ylim = c(-0.5, 2.5),
     xlab = "", ylab = "", main = "3x3 Grid with Point IDs", xaxt = "n", yaxt = "n")
points(grid_points$x, grid_points$y, pch = 16, cex = 2)
text(grid_points$x, grid_points$y, labels = grid_points$id, pos = 3, offset = 0.8)
axis(1, at = 0:2)
axis(2, at = 0:2)
grid()


draw_connection <- function(p1, p2, dist_label, color = "black", lwd = 1, lty = 1) {
  lines(c(grid_points$x[p1], grid_points$x[p2]), 
        c(grid_points$y[p1], grid_points$y[p2]), 
        col = color, lwd = lwd, lty = lty)
  
  mid_x <- (grid_points$x[p1] + grid_points$x[p2]) / 2
  mid_y <- (grid_points$y[p1] + grid_points$y[p2]) / 2
  
  text(mid_x, mid_y, labels = dist_label, col = color, cex = 0.8)
}

# Plot 2: Examples of each distance batch
plot(grid_points$x, grid_points$y, type = "n", xlim = c(-0.5, 2.5), ylim = c(-0.5, 2.5),
     xlab = "", ylab = "", main = "Grid with distance presentation", xaxt = "n", yaxt = "n")
points(grid_points$x, grid_points$y, pch = 16, cex = 2)
text(grid_points$x, grid_points$y, labels = grid_points$id, pos = 3, offset = 0.8)
axis(1, at = 0:2)
axis(2, at = 0:2)
grid()


colors <- c("red", "blue", "green", "purple", "orange")
line_types <- c(1, 2, 3, 4, 5)

for (i in 1:length(unique_distances)) {
  d <- unique_distances[i]
  pairs <- distance_batches[[as.character(d)]]
  
  if (nrow(pairs) > 0) {

    p1 <- pairs[1, 1]
    p2 <- pairs[1, 2]
    
    draw_connection(p1, p2, paste0("d", i, "=", d), 
                    color = colors[i], lwd = 2, lty = line_types[i])
  }
}


legend("bottomright", 
       legend = paste("d", 1:length(unique_distances), "=", unique_distances, 
                      " (", batch_counts, " pairs)", sep = ""),
       col = colors, lty = line_types, lwd = 2, cex = 0.8)

# Plot 3: All pairs with the same distance color
plot(grid_points$x, grid_points$y, type = "n", xlim = c(-0.5, 2.5), ylim = c(-0.5, 2.5),
     xlab = "", ylab = "", main = "All Pairs in Each Distance Batch", xaxt = "n", yaxt = "n")
points(grid_points$x, grid_points$y, pch = 16, cex = 2)
text(grid_points$x, grid_points$y, labels = grid_points$id, pos = 3, offset = 0.8)
axis(1, at = 0:2)
axis(2, at = 0:2)
grid()

for (i in 1:length(unique_distances)) {
  d <- unique_distances[i]
  pairs <- distance_batches[[as.character(d)]]
  
  for (j in 1:nrow(pairs)) {
    p1 <- pairs[j, 1]
    p2 <- pairs[j, 2]
    
    lines(c(grid_points$x[p1], grid_points$x[p2]), 
          c(grid_points$y[p1], grid_points$y[p2]), 
          col = colors[i], lwd = 1, lty = line_types[i])
  }
}

legend("bottomright", 
       legend = paste("d", 1:length(unique_distances), "=", unique_distances, 
                      " (", batch_counts, " pairs)", sep = ""),
       col = colors, lty = line_types, lwd = 2, cex = 0.8)

dev.off()

# Create a table with batch summary
batch_table <- data.frame(
  Batch = paste0("d", 1:length(unique_distances)),
  Distance = unique_distances,
  Number_of_Pairs = batch_counts,
  Description = c(
    "Adjacent points (horizontally or vertically)",
    "Diagonal points",
    "Points 2 units apart (horizontally or vertically)",
    "Knight's move pattern",
    "Opposite corners"
  )
)

print(batch_table)

# Alternative approach using igraph for better visualization
if(require(igraph)) {
  # Create a graph where edges represent pairs in each batch
  g <- graph.empty(n = n_points, directed = FALSE)
  V(g)$name <- 1:n_points
  V(g)$x <- grid_points$x
  V(g)$y <- grid_points$y
  
  # Add edges for each distance batch
  for (i in 1:length(unique_distances)) {
    d <- unique_distances[i]
    pairs <- distance_batches[[as.character(d)]]
    
    for (j in 1:nrow(pairs)) {
      g <- add_edges(g, c(pairs[j, 1], pairs[j, 2]), 
                     weight = d, 
                     color = colors[i],
                     width = 6 - i, # Thicker lines for shorter distances
                     lty = line_types[i])
    }
  }
  
  # Plot the graph
  pdf("grid_visualization_graph.pdf", width = 10, height = 10)
  
  # Plot with all edges
  plot(g, layout = as.matrix(grid_points[, c("x", "y")]),
       vertex.size = 20, 
       vertex.label = V(g)$name,
       vertex.label.cex = 1.2,
       edge.curved = 0.2,
       main = "3x3 Grid with All Distance Batches")
  
  # Add a legend
  legend("bottomright", 
         legend = paste("d", 1:length(unique_distances), "=", unique_distances, 
                        " (", batch_counts, " pairs)", sep = ""),
         col = colors, lty = line_types, lwd = 2, cex = 0.8)
  
  # Separate plots for each distance batch
  par(mfrow = c(2, 3))
  for (i in 1:length(unique_distances)) {
    # Create a subgraph with only edges of this distance
    sub_g <- delete_edges(g, which(E(g)$weight != unique_distances[i]))
    
    plot(sub_g, layout = as.matrix(grid_points[, c("x", "y")]),
         vertex.size = 20, 
         vertex.label = V(sub_g)$name,
         vertex.label.cex = 1.2,
         edge.curved = 0.2,
         main = paste("Distance Batch d", i, "=", unique_distances[i], 
                      " (", batch_counts[i], " pairs)", sep = ""))
  }
  
  dev.off()
}

