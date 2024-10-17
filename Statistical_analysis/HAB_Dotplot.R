#############
#HAB dotplot
############

#---Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(ggdendro)

#library(ggdendro)
library(indicspecies)

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

#---Long Format, as per usual Plus prepare for ggplot
long_data <- merged_data %>%
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Value")

# Reorder the columns
long_data <- long_data %>%
  select(Sample, Species, Value)

#
long_data$Species = factor(x = long_data$Species, levels = unique(merged_data$Species))
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
  geom_point(data = long_data, aes(x = Sample, y = Species, size = Value, , color = Combined)) +
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

dot = ggplot() +
  geom_point(data = long_data %>% filter(Value != 0), 
             aes(x = Sample, y = Species, size = Value, color = Combined),
             shape = 16) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(11, "Dark2")[4],
                                RColorBrewer::brewer.pal(11, "Dark2")[3],
                                RColorBrewer::brewer.pal(11, "Dark2")[2],
                                RColorBrewer::brewer.pal(11, "Dark2")[1])) +
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
dot

#---multipat?
#Scale to remove zeros
df_scaled = as.data.frame(merged_data)[-1]
df_scaled[, c(1:ncol(df_scaled))] <- scale(df_scaled[, 1:ncol(df_scaled)])

clusters = t(df_scaled)
abundance = clusters
abundance = as.data.frame(abundance)
colnames(abundance) = merged_data$Species
association = c("East Bottom", "East Surface", "East Bottom", "East Surface","East Bottom", "East Surface","East Bottom", "East Surface",
                "West Bottom", "West Surface", "West Bottom", "West Surface", "West Bottom", "West Surface", "West Bottom", "West Surface", "West Bottom", "West Surface", "West Bottom", "West Surface")
indicator_r.g <- multipatt(abundance, association, func = "r.g", control = how(nperm = 9999))
summary(indicator_r.g)

#run clustering
df_matrix <- as.matrix(df_scaled)
rownames(df_matrix) = merged_data$Species
hclust_obj <- hclust(d = dist(x = df_matrix))
df_dendro <- as.dendrogram(hclust_obj)

# Create dendro
dendro_plot <- ggdendrogram(data = df_dendro, rotate = 90)

#---ANOVA of HAB relative abundance
long_data_summarized <- long_data %>%
  group_by(Sample, Station, Region, Depth, Combined) %>%
  summarize(Value = sum(Value)) %>%
  ungroup()

box = ggplot() +
  geom_boxplot(data = long_data_summarized, 
               aes(x = Combined, y = Value, fill = Combined),
               outlier.size = 2) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(11, "Dark2")[4],
                               RColorBrewer::brewer.pal(11, "Dark2")[3],
                               RColorBrewer::brewer.pal(11, "Dark2")[2],
                               RColorBrewer::brewer.pal(11, "Dark2")[1])) +
  ylab("Relative HAB species abundance") +
  theme_minimal() +
  labs(fill = NULL) +
  theme_classic() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10))
box

#Anova
anova_results <- aov(Value ~ Depth * Region, data = long_data_summarized)
summary(anova_results)

#Tukey
TukeyHSD(anova_results)


#---plot

