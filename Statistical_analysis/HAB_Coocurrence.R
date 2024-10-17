###############
#HAB Co-occurrence
##############

#---library
library(ggplot2)
library(igraph)
library(ggraph)
library(dplyr)
library(tidyr)

#---data
data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

#Relative abundance
V9data[, 3:(ncol(V9data)-7)] <- apply(V9data[, 3:(ncol(V9data)-7)], 2, function(x) x / sum(x))

#---Tag HAB species
#First, the IOC-UNESCO Taxonomic Reference List of Harmful Micro Algae
#https://www.marinespecies.org/hab/
HABs_taxlist_20240712 <- read_excel("HABs_taxlist_20240712.xlsx")

#Get the unique accepted names
HAB_list = unique(HABs_taxlist_20240712$ScientificName_accepted)

#Load red tide sighting records from previous study or from HKEPD, worms corrected
Red_tide_and_WCZ_02162024 <- read.csv("~/Prelim_phyto_analysis/Red_tide_and_WCZ_02162024.csv")
HAB_list = append(HAB_list, unique(Red_tide_and_WCZ_02162024$Species))

HAB_list = unique(HAB_list)

#Tag species that are HAB species
V9data$HAB_status <- ifelse(V9data$Species %in% HAB_list, "yes", "no")

V9data = subset(V9data, V9data$HAB_status == "yes") #106 unique ASVs

#names back
colnames(V9data) =  c(colnames(V9data)[1],sub("^[^.]+\\.", "", colnames(V9data)[-1]))

# Group by the ASV number and check the Species values
new_data <- V9data %>%
  group_by(ASV.number, Species) %>%
  summarize(across(where(is.numeric), sum), .groups = "drop") %>%
  group_by(ASV.number) %>%
  mutate(unique_species = n_distinct(Species)) %>%
  # If there are multiple Phylums, set the Species to "Ambiguous"
  mutate(Species = ifelse(unique_species > 1, ASV.number, Species)) %>%
  # Keep only the first row for each Species
  slice(1)
new_data$unique_species = NULL

#Combine by Species 
new_data = new_data[-1]

merged_data <- new_data %>%
  group_by(Species) %>%
  summarize(across(where(is.numeric), sum)) %>%
  mutate(total = rowSums(select(., -Species))) %>%
  arrange(desc(total))

merged_data$total = NULL
merged_data = as.data.frame(t(merged_data))
colnames(merged_data) = merged_data[1, ]
merged_data = merged_data[-1,]

merged_data_numeric <- as.data.frame(lapply(merged_data, as.numeric))

wide2 = merged_data_numeric
colnames(wide2) = colnames(merged_data)


#Observation matrix
m_obs = as.matrix(wide2)
co_occurrence <- t(m_obs) %*% m_obs

# Create a graph from the co-occurrence matrix
graph <- graph_from_adjacency_matrix(co_occurrence, mode = "undirected", weighted = TRUE)

# Set the vertex names as species labels
V(graph)$name <- colnames(merged_data)

# Determine the abundance of each species
abundance <- colSums(wide2)
names(abundance) = colnames(wide2)

# Calculate Pearson correlation coefficients
correlation <- cor(m_obs)

# Calculate p-values for correlations
p_values <- matrix(0, ncol = ncol(m_obs), nrow = ncol(m_obs))
for (i in 1:(ncol(m_obs) - 1)) {
  for (j in (i + 1):ncol(m_obs)) {
    p_values[i, j] <- p_values[j, i] <- cor.test(m_obs[, i], m_obs[, j])$p.value
  }
}

# Define the threshold for significance
cor_threshold <- 0.6
p_value_threshold <- 0.05

# Filter edges based on significance criteria
correlated_edges <- E(graph)[abs(correlation[as.integer(E(graph))]) > cor_threshold]
significant_edges <- E(graph)[p_values[as.integer(E(graph))] < p_value_threshold]

# Create a subgraph with only significant edges
subgraph <- delete_edges(graph, E(graph)[!E(graph) %in% correlated_edges])
subgraph <- delete_edges(subgraph, E(subgraph)[!E(subgraph) %in% significant_edges])

subgraph = delete.vertices(subgraph , which(degree(subgraph)==0))#delete points with no connections
abundance = abundance[names(abundance) %in% V(subgraph)$name]

# Positive and negative correlations
neg_edges <- E(subgraph)[correlation[as.integer(E(subgraph))] < 0]
pos_edges <- E(subgraph)[correlation[as.integer(E(subgraph))] > 0]


# Assign different colors to each species
species_colors <- rainbow(length(V(subgraph)))

# Create a color vector for edges based on their sign (positive/negative)
edge_colors <- ifelse(E(subgraph) %in% neg_edges, "negative", "positive") #The colours can't be set as they are already in use by scale_colour manual actually, so there is o pink or green.
Correlations = edge_colors

# Plot the co-occurrence network with significant correlations 
#width = abs(E(subgraph)$weight)
#alpha = abs(E(subgraph)$weight)

ggraph(subgraph, layout = 'dh') +
  geom_edge_arc(strength = 0.2, aes(color = Correlations, alpha = abs(E(subgraph)$weight), width = abs(E(subgraph)$weight))) +
  geom_node_point(aes(size = abundance, color = name), alpha = 0.5) +
  scale_size_continuous(range = c(2, 10)) +
  scale_color_manual(values = species_colors) +
  scale_edge_width(range = c(1, 5), guide = "none") +
  scale_edge_alpha(range = c(0.5, 1), guide = "none") + 
  geom_node_text(aes(label = name), col = "darkblue", size = 3, repel = T, force = 0.004) +
  labs(size = "Relative\nAbundance") + 
  theme_void() +
  guides(color = "none", width = "none")

num_pink_lines <- length(E(subgraph)[correlation[as.integer(E(subgraph))] < 0])
num_green_lines <- length(E(subgraph)[correlation[as.integer(E(subgraph))] > 0])
num_pink_lines
num_green_lines
