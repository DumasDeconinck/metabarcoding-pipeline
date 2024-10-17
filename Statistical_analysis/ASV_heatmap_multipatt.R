#############
#Distribution
#############

#---library
library(ggplot2)
library(vegan)
library(ape)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(ggdendro)
library(grid)
library(tidyr)
library(dplyr)
library(stringr)
library(ggtext)
library(indicspecies)

#---load the data
data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

#Keep onll Phylum and ASV number
V9data = V9data[, !grepl("Sequence|Kindom|Domain|Phylum|Class|Order|Family|Genus|Species", names(V9data))]

#Remove duplicate ASVs
V9data <- V9data %>%
  distinct(ASV.number, .keep_all = TRUE)

colnames(V9data) =  c(colnames(V9data)[1],sub("^[^.]+\\.", "", colnames(V9data)[-1]))

#---Convert to relative abundance
V9data[, 2:ncol(V9data)] <- apply(V9data[, 2:ncol(V9data)], 2, function(x) x / sum(x))

#---Only for universal ASVs
#V9data[V9data == 0] <- NA
#V9data = V9data[complete.cases(V9data),]

#---Long Format, as per usual Plus prepare for ggplot
long_data <- V9data %>%
  pivot_longer(cols = -ASV.number, names_to = "Sample", values_to = "Value")

# Reorder the columns
long_data = data.frame(long_data$Sample, long_data$ASV.number, long_data$Value)
colnames(long_data) = c("Sample", "ASV.number", "Value")

#---dendrogram
tV9data <- as.data.frame(t(V9data))
colnames(tV9data) = V9data$ASV.number
try = tV9data[-1,]
try = V9data[,-1]


# Convert columns to numeric if possible
for (col in names(try)) {
  if (is.character(try[[col]])) {
    # Try to convert the column to numeric
    try[[col]] <- as.numeric(try[[col]])
  }
}

#Scale to remove zeros
df_scaled = try
df_scaled[, c(1:ncol(try))] <- scale(try[, 1:ncol(try)])

#run clustering
df_matrix <- as.matrix(df_scaled)
rownames(df_matrix) = V9data$ASV.number
hclust_obj <- hclust(d = dist(x = df_matrix))
df_dendro <- as.dendrogram(hclust_obj)

# Create dendro
dendro_plot <- ggdendrogram(data = df_dendro, rotate = 90)

#Analyze the clustering
cutree_hclust <- cutree(hclust_obj, k = 4)
df_list <- lapply(1 : 4, function(x) df_matrix[which(cutree_hclust == x), ])

#names(df_list[[1]][,1]) #Doesn't work for singlets
#names(df_list[[2]][,1])
#names(df_list[[3]][,1])
#names(df_list[[4]][,1])

g1 = V9data$ASV.number[which(cutree_hclust == 1)]
g2 = V9data$ASV.number[which(cutree_hclust == 2)]
g3 = V9data$ASV.number[which(cutree_hclust == 3)]
g4 = V9data$ASV.number[which(cutree_hclust == 4)]

g1 = names(cutree_hclust[which(cutree_hclust == 1)])
g2 = names(cutree_hclust[which(cutree_hclust == 2)])
g3 = names(cutree_hclust[which(cutree_hclust == 3)])
g4 = names(cutree_hclust[which(cutree_hclust == 4)])

plot(hclust_obj)
rect.hclust(hclust_obj, k = 4, border = 2:6)
rplot = recordPlot()

colours = brewer.pal(4, "Dark2")

dend_data = dendro_data(df_dendro, type = "rectangle")

names(dend_data)

head(dend_data$segments)

head(dend_data$labels)

dend_data$group <- ifelse(dend_data$labels$label %in% g1, "Group 1",
                          ifelse(dend_data$labels$label %in% g2, "Group 2",
                                 ifelse(dend_data$labels$label %in% g3, "Group 3",
                                        ifelse(dend_data$labels$label %in% g4, "Group 4", NA))))

p <- ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label, colour = dend_data$group),
            hjust = 1, angle = 90, size = 3)
print(p)

ggdendrogram(dend_data) + geom_text(data = dend_data$labels, aes(x, y, label = label, colour = dend_data$group),
                                    hjust = 1, angle = 90, size = 3)

clusters = df_scaled

rownames(clusters)
#rownames(clusters) = V9data$ASV.number

clusters$Association <- ifelse(rownames(clusters) %in% g1, "Group 1",
                               ifelse(rownames(clusters) %in% g2, "Group 2",
                                      ifelse(rownames(clusters) %in% g3, "Group 3",
                                             ifelse(rownames(clusters) %in% g4, "Group 4", NA))))

abundance = clusters[1:ncol(clusters)-1]

association = clusters$Association

abundance[is.na(abundance)] <- 0


# Call the 'multipatt' function with the converted 'abundance' variable
#indicator_r.g <- multipatt(abundance, association, func = "r.g", control = how(nperm = 9999))
#summary(indicator_r.g)

#order
df_order <- order.dendrogram(df_dendro)

#---Original clusters multipatt analysis
tV9data <- as.data.frame(t(V9data))
try = tV9data[-1,]
colnames(try) = tV9data[1,]

# Convert columns to numeric if possible
for (col in names(try)) {
  if (is.character(try[[col]])) {
    # Try to convert the column to numeric
    try[[col]] <- as.numeric(try[[col]])
  }
}

#Scale to remove zeros
df_scaled = try
df_scaled[, c(1:ncol(try))] <- scale(try[, 1:ncol(try)])

clusters = df_scaled
abundance = clusters[1:ncol(clusters)-1]
association = c("East Bottom", "East Surface", "East Bottom", "East Surface","East Bottom", "East Surface","East Bottom", "East Surface",
                "West Bottom", "West Surface", "West Bottom", "West Surface", "West Bottom", "West Surface", "West Bottom", "West Surface", "West Bottom", "West Surface", "West Bottom", "West Surface")
indicator_r.g <- multipatt(abundance, association, func = "r.g", control = how(nperm = 9999))
summary(indicator_r.g)


# Order the levels according to their position in the cluster
long_data$ASV.number = factor(x = long_data$ASV.number, levels = hclust_obj$labels, ordered = TRUE)
long_data$Sample = factor(x = long_data$Sample, levels = c("S1S", "S2S", "S3S", "S4S", "S5S", "S6S",
                                                            "S1B", "S2B", "S3B", "S4B", "S5B", "S6B",
                                                            "P1S", "P2S", "P3S", "P4S",
                                                            "P1B", "P2B", "P3B", "P4B"), ordered = TRUE)

long_data

# Plot the data as a heatmap
plot <- ggplot() +
  geom_tile(data = long_data, aes(x = Sample, y = ASV.number, fill = Value)) +
  scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(11, "RdBu")[2])) +
  labs(x = "Sample", y = "ASV", fill = "Relative\nAbundance") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.2, "lines"),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA))
plot

#try cleaner approach (can't stop it from mirroring)

library(tidyverse)

# Create the exponential sequence
exp_seq <- exp(seq(from = 0, to = 10, length.out = 582))

# Repeat each number 20 times
result <- rev(rep(exp_seq, each = 20))

long_data$height = result

plot <- ggplot() +
  geom_tile(data = long_data, aes(x = Sample, y = ASV.number, fill = Value, height = height)) +
  scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(11, "RdBu")[2])) +
  labs(x = "Sample", y = "ASV", fill = "Relative\nAbundance") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.2, "lines"),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA))
plot

#For top ASVs
top = 20
plot <- ggplot() +
  geom_tile(data = long_data[1:(top*20),], aes(x = Sample, y = ASV.number, fill = Value)) +
  scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(11, "RdBu")[2])) +
  labs(x = "Sample", y = "ASV", fill = "Relative\nAbundance") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.2, "lines"),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA))
plot

#--- univresal ASVs
#---load the data
data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

#Keep onll Phylum and ASV number
V9data = V9data[, !grepl("Sequence|Kindom|Domain|Phylum|Class|Order|Family|Genus|Species", names(V9data))]

#Remove duplicate ASVs
V9data <- V9data %>%
  distinct(ASV.number, .keep_all = TRUE)

colnames(V9data) =  c(colnames(V9data)[1],sub("^[^.]+\\.", "", colnames(V9data)[-1]))

#---Convert to relative abundance
V9data[, 2:ncol(V9data)] <- apply(V9data[, 2:ncol(V9data)], 2, function(x) x / sum(x))

#---Only for universal ASVs
V9data[V9data == 0] <- NA
V9data = V9data[complete.cases(V9data),]

#---Long Format, as per usual Plus prepare for ggplot
long_data <- V9data %>%
  pivot_longer(cols = -ASV.number, names_to = "Sample", values_to = "Value")

# Reorder the columns
long_data <- long_data %>%
  select(Sample, ASV.number, Value)

#---dendrogram
tV9data <- as.data.frame(t(V9data))
colnames(tV9data) = V9data$ASV.number
try = tV9data[-1,]
try = V9data[,-1]


# Convert columns to numeric if possible
for (col in names(try)) {
  if (is.character(try[[col]])) {
    # Try to convert the column to numeric
    try[[col]] <- as.numeric(try[[col]])
  }
}

#Scale to remove zeros
df_scaled = try
df_scaled[, c(1:ncol(try))] <- scale(try[, 1:ncol(try)])

#run clustering
df_matrix <- as.matrix(df_scaled)
rownames(df_matrix) = V9data$ASV.number
hclust_obj <- hclust(d = dist(x = df_matrix))

long_data$ASV.number = factor(x = long_data$ASV.number, levels = hclust_obj$labels, ordered = TRUE)
long_data$Sample = factor(x = long_data$Sample, levels = c("S1S", "S2S", "S3S", "S4S", "S5S", "S6S",
                                                           "S1B", "S2B", "S3B", "S4B", "S5B", "S6B",
                                                           "P1S", "P2S", "P3S", "P4S",
                                                           "P1B", "P2B", "P3B", "P4B"), ordered = TRUE)

long_data$Station <- sub("^[^.]+\\.", "", long_data$Sample)
long_data$Region <- ifelse(substr(long_data$Station, 1, 1) == "P", "East",
                               ifelse(substr(long_data$Station, 1, 1) == "S", "West", "Unknown"))
long_data$Depth <- ifelse(substr(long_data$Station, nchar(long_data$Station), nchar(long_data$Station)) == "S", "Surface",
                              ifelse(substr(long_data$Station, nchar(long_data$Station), nchar(long_data$Station)) == "B", "Bottom", "Unknown"))
long_data$Combined = paste(long_data$Region, long_data$Depth)


ggplot() +
  geom_point(data = long_data, aes(x = Sample, y = ASV.number, size = Value, color = Combined)) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(11, "Dark2")[4],RColorBrewer::brewer.pal(11, "Dark2")[3],RColorBrewer::brewer.pal(11, "Dark2")[2],RColorBrewer::brewer.pal(11, "Dark2")[1])) +
  theme_minimal() +
  labs(size = "Relative\nAbundance", color = NULL) +
  theme(axis.title = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(color = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
        panel.border = element_rect(color = "black", fill = NA),
        legend.key.size = unit(1.2, "lines"),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(0, "lines"))

