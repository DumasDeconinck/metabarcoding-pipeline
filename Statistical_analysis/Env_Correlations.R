###########
#ASV_Environmental correlations
###########

#---libraries
library(vegan)
library(dplyr)
library(tidyr)
library(readxl)

library(RColorBrewer)
library(ggplot2)
library(ggrepel)

library(linkET)

#---data
data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

#Keep onll Phylum and ASV number
V9data = V9data[, !grepl("Sequence|Kindom|Domain|Phylum|Class|Order|Family|Genus|Species", names(V9data))]
V9data <- V9data %>%
  distinct(ASV.number, .keep_all = TRUE)

colnames(V9data) =  c(colnames(V9data)[1],sub("^[^.]+\\.", "", colnames(V9data)[-1]))

#Relative abundance
V9data[, 2:ncol(V9data)] <- apply(V9data[, 2:ncol(V9data)], 2, function(x) x / sum(x))

#Remove Bottom data
V9data<- V9data %>%
  select(!ends_with("B"))

df_numeric = as.data.frame(t(V9data))
colnames(df_numeric) = df_numeric[1, ]
df_numeric = df_numeric[-1,]

df_numeric <- as.data.frame(lapply(df_numeric, as.numeric))

#Environmental data
HK_env  = read_excel("Environmental data.xlsx")
HK_env$`EC (µS/cm)` = NULL

parameter_names = colnames(HK_env)[(length(colnames(HK_env))-5): length(colnames(HK_env))]

colnames(HK_env) = c("Station", "Depth", "Depth2", "Temperature", "DO", "Salinity", "pH", "Turbidity", "Chl")

HK_env = subset(HK_env, Depth2 == "Surface")
HK_env = as.data.frame(HK_env)
HK_env$Depth = NULL
HK_env$Depth2 = NULL
rownames(HK_env) = HK_env$Station
HK_env$Station = NULL

#---CCA
cc3 = cca(df_numeric ~ ., data = HK_env)

#---Plot
##ggplot
#extracting the data as data frame; env data
veg_1 = as.data.frame(cc3$CCA$biplot)
veg_1["env"] = row.names(veg_1)

#extracting the data; species
veg_2 = as.data.frame(cc3$CCA$u)
veg_2["sites"] = rownames(HK_env)

#Add Grouping variables
rownames(veg_2) = veg_2$sites
veg_2$sites
veg_2$Region <- ifelse(substr(rownames(veg_2), 1, 1) == "P", "East",
                             ifelse(substr(rownames(veg_2), 1, 1) == "S", "West", "Unknown"))

## add polygon for Region
#Prepare polygon
grp.a <- subset(veg_2, Region == unique(veg_2$Region)[1])
grp.b <- subset(veg_2, Region == unique(veg_2$Region)[2])

grp.a <- veg_2[veg_2$Region == unique(veg_2$Region)[1], ][chull(veg_2[veg_2$Region == 
                                                                  unique(veg_2$Region)[1], c("CCA1", "CCA2")]), ]  # hull values
grp.b <- veg_2[veg_2$Region == unique(veg_2$Region)[2], ][chull(veg_2[veg_2$Region == 
                                                                  unique(veg_2$Region)[2], c("CCA1", "CCA2")]), ]  

hull.data <- rbind(grp.a, grp.b)  # combine grp.a and grp.b
hull.data

veg_1$env = c("Temperature (°C)", "DO saturation (%)", "Salinity (PSU)", "pH", "Turbidity (FNU)", "Chlorphyl a (µg/L)")

#Plot
color_palette <- brewer.pal(8, "Dark2")  # Adjust the number of colors as needed
graph <- ggplot() +
  geom_polygon(data = hull.data, aes(x = CCA1, y = CCA2, fill = Region, group = Region), alpha = 0.3) +
  geom_point(data = veg_2, aes(x = CCA1, y = CCA2), size = 1) +
  geom_segment(
    data = veg_1,
    aes(
      x = 0,
      y = 0,
      xend = CCA1*3,#note the times 3, so inaccurate, but makes it legible.
      yend = CCA2*3
    ), colour = "blue",
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_text_repel(
    data = veg_1,
    aes(x = CCA1*3, y = CCA2*3, label = veg_1$env),
    nudge_y = -0.05,
    color = "red",
    size = 3
  ) +
  geom_text_repel(data = veg_2, aes(x = CCA1, y = CCA2, label = veg_2$sites)) +
  scale_colour_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  labs(x = "CCA1", y = "CCA2") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Display the graph
print(graph)

#--cca analysis
cc3
# Anova
AovRes = anova(cc3,by="mar",permutations=999) #analyse significances of marginal effects (“TypeIII effects”)
AovRes

#Analysis of cca
#Monte carlo
perm_test <- anova.cca(cc3, parallel = TRUE, permutations = 999)
perm_test


#---Mantle Plots
#--Phylum all ASV
#data again
data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

# Group by the ASV number and check the Phylum values
new_data <- V9data %>%
  group_by(ASV.number, Phylum) %>%
  summarize(across(where(is.numeric), sum), .groups = "drop") %>%
  group_by(ASV.number) %>%
  mutate(unique_phylums = n_distinct(Phylum)) %>%
  # If there are multiple Phylums, set the Phylum to "Ambiguous"
  mutate(Phylum = ifelse(unique_phylums > 1, "Ambiguous", Phylum)) %>%
  # Keep only the first row for each ASV.number
  slice(1)
new_data$unique_phylums = NULL

new_data$ASV.number = NULL
V9data = new_data

colnames(V9data) =  c(colnames(V9data)[1],sub("^[^.]+\\.", "", colnames(V9data)[-1]))

#Relative abundance
V9data[, 2:ncol(V9data)] <- apply(V9data[, 2:ncol(V9data)], 2, function(x) x / sum(x))

#Removes Phyla that weren't identified
V9data = subset(V9data, Phylum != "Ambiguous") 

#Remove Bottom data
V9data<- V9data %>%
  select(!ends_with("B"))

#Aggregate the Phyla together
V9data_aggregated <- V9data %>%
  group_by(Phylum) %>%
  summarise(across(where(is.numeric), sum))

df_numeric = as.data.frame(t(V9data_aggregated))
colnames(df_numeric) = df_numeric[1, ]
df_numeric = df_numeric[-1,]

df_numeric <- as.data.frame(lapply(df_numeric, as.numeric))

#Mantel
#Mantel test becomes illegible with too much species. Additionally, it takes forever, and it is nonsensical to look at species that don't appear much
#
column_names <- colnames(df_numeric)

column_list <- as.list(column_names)

# Get the number of columns in df_env excluding the first column
num_env_columns <- ncol(HK_env)

# Repeat the species column names based on the number of df_env columns
repeated_column_names <- rep(column_names, each = num_env_columns)

mantel <- mantel_test(df_numeric, HK_env,
                      spec_select = column_list) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

mantel$spec = repeated_column_names #this is to ensure the species names are displayed on the graph

mantela = qcorrplot(correlate(HK_env), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature(),
              nudge = .1,
              label.size = 3) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position='left',
        legend.key.size = unit(1.2, "lines"),
        legend.background = element_rect(fill = "white", color = "black"))
mantela

cora = correlate(df_numeric, HK_env) %>% 
  qcorrplot() +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        legend.position='left',
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(1.2, "lines"),
        legend.background = element_rect(fill = "white", color = "black"))
cora

#--Top HAB species
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

#Remove Bottom data
merged_data<- merged_data %>%
  select(!ends_with("B"))

merged_data = as.data.frame(t(merged_data))
colnames(merged_data) = merged_data[1, ]
merged_data = merged_data[-1,]

merged_data_numeric <- as.data.frame(lapply(merged_data, as.numeric))

df_numeric = merged_data_numeric[1:15]
colnames(df_numeric) = gsub("\\."," " ,colnames(df_numeric))

##Mantel
column_names <- colnames(df_numeric)

column_list <- as.list(column_names)

# Get the number of columns in df_env excluding the first column
num_env_columns <- ncol(HK_env)

# Repeat the species column names based on the number of df_env columns
repeated_column_names <- rep(column_names, each = num_env_columns)

mantel <- mantel_test(df_numeric, HK_env,
                      spec_select = column_list) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

mantel$spec = repeated_column_names #this is to ensure the species names are displayed on the graph

mantelb = qcorrplot(correlate(HK_env), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature(),
              nudge = .1,
              label.size = 3) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position='left',
        legend.key.size = unit(1.2, "lines"),
        legend.background = element_rect(fill = "white", color = "black"))
mantelb

colnames(merged_data_numeric) = gsub("\\."," " ,colnames(merged_data_numeric))
corb = correlate(merged_data_numeric, HK_env) %>% 
  qcorrplot() +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        legend.position='left',
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(1.2, "lines"),
        legend.background = element_rect(fill = "white", color = "black"))
corb
