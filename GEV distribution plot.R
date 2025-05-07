# Visualization of Univariate Extreme Value Theory Distributions


library(evd)       # For extreme value distributions
library(ggplot2)   
library(reshape2)  

# Set random seed for reproducibility
set.seed(123)


# Generalized Extreme Value (GEV) Distribution Family


# Create data points
x <- seq(-3, 6, by = 0.01)

# Generate GEV densities with different shape parameters
gumbel <- dgev(x, loc = 0, scale = 1, shape = 0)     # Gumbel (shape = 0)
frechet1 <- dgev(x, loc = 0, scale = 1, shape = 0.2) # Fréchet (shape > 0)
frechet2 <- dgev(x, loc = 0, scale = 1, shape = 0.5) # More extreme Fréchet
weibull1 <- dgev(x, loc = 0, scale = 1, shape = -0.2) # Weibull (shape < 0)
weibull2 <- dgev(x, loc = 0, scale = 1, shape = -0.5) # More extreme Weibull

# Create dataframe for ggplot
df <- data.frame(
  x = rep(x, 5),
  density = c(gumbel, frechet1, frechet2, weibull1, weibull2),
  Distribution = factor(rep(c("Gumbel (ξ = 0)", 
                              "Fréchet (ξ = 0.2)", 
                              "Fréchet (ξ = 0.5)",
                              "Weibull (ξ = -0.2)",
                              "Weibull (ξ = -0.5)"), each = length(x)))
)

# Plot GEV family distributions
p1 <- ggplot(df, aes(x = x, y = density, color = Distribution)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Generalized Extreme Value Distribution Family",
       subtitle = "Comparison of different shape parameters (ξ)",
       x = "x", 
       y = "Density") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("red", "blue", "darkblue", "green", "darkgreen"))

# Display plot
print(p1)





# Tail Behavior Visualization

# Create data points focusing on the tail
x_tail <- seq(1, 10, by = 0.01)

# Generate right tails for different GEV shape parameters
gumbel_tail <- dgev(x_tail, loc = 0, scale = 1, shape = 0)
frechet_tail <- dgev(x_tail, loc = 0, scale = 1, shape = 0.5)
weibull_tail <- dgev(x_tail, loc = 0, scale = 1, shape = -0.5)


df_tail <- data.frame(
  x = rep(x_tail, 3),
  density = c(gumbel_tail, frechet_tail, weibull_tail),
  Distribution = factor(rep(c("Gumbel (ξ = 0) - Light tail", 
                              "Fréchet (ξ > 0) - Heavy tail", 
                              "Weibull (ξ < 0) - Bounded tail"), each = length(x_tail)))
)

# Plot tail behavior
p2 <- ggplot(df_tail, aes(x = x, y = density, color = Distribution)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Tail Behavior of Extreme Value Distributions",
       subtitle = "Density functions on log scale",
       x = "x", 
       y = "Density (log scale)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("red", "blue", "green")) +
  scale_y_log10() # Log scale to highlight tail behavior differences


print(p2)
