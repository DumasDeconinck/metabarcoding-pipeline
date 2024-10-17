######
#HAB
#####

#---libraries
library(tidyr)
library(dplyr)
library(readxl)

#---datatidytree#---data
data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

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

#---Create dataframe with HAB species
#-one asv multi spec
# Find the duplicate ASV.number values
duplicate_ASVs <- V9data$ASV.number[duplicated(V9data$ASV.number)]

# Extract the rows with duplicate ASV.number values
data_with_duplicates <- V9data[V9data$ASV.number %in% duplicate_ASVs, ]

length(unique(data_with_duplicates$ASV.number))
length(unique(data_with_duplicates$Species))

#-multi asv one spec
#Remove the Multiple_ASV-One_species
data2 = V9data[ !(V9data$ASV.number %in% data_with_duplicates$ASV.number), ]

#All species that are duplicates in the remaining dataset
duplicate_species <- data2$Species[duplicated(data2$Species)]
data_with_duplicate_species <- data2[data2$Species %in% duplicate_species, ]

length(unique(data_with_duplicate_species$ASV.number))
length(unique(data_with_duplicate_species$Species))

#-One_ASV-One_Species
data3 = data2[ !(data2$ASV.number %in% data_with_duplicate_species$ASV.number), ]
nrow(data3)

#-dataframe
#One-One
data3$ASVs = 1
data3b <- data.frame("Group" = 1, data3$Species, data3$ASVs, data3$ASV.number, data3$Phylum, data3$Class)
colnames(data3b) = c("Group", "Species", "ASV numbers", "ASV_ID", "Phylum", "Class")

#Multi-One
data2_unique <- data_with_duplicate_species %>%
  group_by(Species) %>%
  summarize(ASV_count = n_distinct(ASV.number),
            Phylum = first(Phylum),
            Class = first(Class))

?summarize

data2_unique$Species <- factor(data2_unique$Species, levels = unique(data_with_duplicate_species$Species))

# Order the dada2_unique dataframe based on the ordered factor variable
data2_unique <- data2_unique[order(data2_unique$Species), ]

#Create
data2b <- data.frame("Group" = 2, data2_unique$Species, data2_unique$ASV_count, "ASV_ID" = "-", data2_unique$Phylum, data2_unique$Class)
colnames(data2b) = c("Group", "Species", "ASV numbers", "ASV_ID", "Phylum", "Class")

##one-multi
data1_unique = data_with_duplicates %>%
  group_by(Species) %>%
  summarize(ASV_count = n_distinct(ASV.number),
            Phylum = first(Phylum),
            Class = first(Class))

data1_unique$Species <- factor(data1_unique$Species, levels = unique(data_with_duplicates$Species))

# Order the dada2_unique dataframe based on the ordered factor variable
data1_unique <- data1_unique[order(data1_unique$Species), ]

data1b <- data.frame("Group" = 3, data1_unique$Species, data1_unique$ASV_count, "ASV_ID" = "-", data1_unique$Phylum, data1_unique$Class)
colnames(data1b) = c("Group", "Species", "ASV numbers", "ASV_ID", "Phylum", "Class")

#---Match against red tide sightings
Redtides = read.csv("Red_tide_and_WCZ_02162024.csv")
redtides = unique(Redtides$Species)

#G1
G1sighted = data3b$Species %in% redtides
data3b$HK_Red_tide = G1sighted
table(G1sighted)
subset(data3b, data3b$HK_Red_tide == TRUE)

#G2
G2sighted = data2b$Species %in% redtides
data2b$HK_Red_tide = G2sighted
table(G2sighted)
subset(data2b, data2b$HK_Red_tide == TRUE)

#G3
G3sighted = data1b$Species %in% redtides
data1b$HK_Red_tide = G3sighted
table(G3sighted)
subset(data1b, data1b$HK_Red_tide == TRUE)

#Unique G3
certain = unique(c(as.vector(data2b$Species), data3b$Species))
which(data1b$Species %in% certain)

certain2 = unique(c(as.vector(subset(data2b, data2b$HK_Red_tide == TRUE)$Species), subset(data3b, data3b$HK_Red_tide == TRUE)$Species))
which(subset(data1b, data1b$HK_Red_tide == TRUE)$Species %in% certain2)

#-Bundle
HAB_Species = rbind.data.frame(data3b, data2b, data1b)
unique(HAB_Species$Species)
colnames(HAB_Species) = c("Group", "Species", "ASV numbers", "ASV ID", "Phylum", "Class", "Red tide")

write.csv(HAB_Species, "ASV_HABs.csv", row.names = FALSE)
