##############
#ASV_inspection
###########


#---Libraries
library(readxl)
#---Read file
data = read.csv("Classified_species_phyto.csv")

#-Only inspect V4 and V9
# Create a new dataframe without the columns to be removed
data = data[, !grepl("Ex|FC|LC", names(data))]
V4data = data[, !grepl("Ex|FC|LC|V9", names(data))]
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]

#Only keep rows with counts
data <- data[rowSums(data[, 3:42] != 0) > 0, ]
V4data <- V4data[rowSums(V4data[, 3:22] != 0) > 0, ]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

#To inspect only one primer pair
#data = V4data
#data = V9data

#---Information
length(unique(data$Species))
length(unique(data$ASV.number))

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
data$HAB_status <- ifelse(data$Species %in% HAB_list, "yes", "no")

#---Phytoplankton xASV-Xspecies
#-one asv multi spec
# Find the duplicate ASV.number values
duplicate_ASVs <- data$ASV.number[duplicated(data$ASV.number)]

# Extract the rows with duplicate ASV.number values
data_with_duplicates <- data[data$ASV.number %in% duplicate_ASVs, ]

length(unique(data_with_duplicates$ASV.number))
length(unique(data_with_duplicates$Species))

#Multi asv-one spec
#Remove the Multiple_ASV-One_species
data2 = data[ !(data$ASV.number %in% data_with_duplicates$ASV.number), ]

#All species that are duplicates in the remaining dataset
duplicate_species <- data2$Species[duplicated(data2$Species)]
data_with_duplicate_species <- data2[data2$Species %in% duplicate_species, ]

length(unique(data_with_duplicate_species$ASV.number))
length(unique(data_with_duplicate_species$Species))

#-One_ASV-One_Species
data3 = data2[ !(data2$ASV.number %in% data_with_duplicate_species$ASV.number), ]
nrow(data3)

#---HAB xASV-Xspecies
#First, create a subset for only habspecies
HAB_data = subset(data, HAB_status == "yes")
length(unique(HAB_data$ASV.number))
length(unique(HAB_data$Species))


#-one asv multi spec
# Find the duplicate ASV.number values
duplicate_ASVs <- HAB_data$ASV.number[duplicated(HAB_data$ASV.number)]

# Extract the rows with duplicate ASV.number values
data_with_duplicates <- HAB_data[HAB_data$ASV.number %in% duplicate_ASVs, ]

length(unique(data_with_duplicates$ASV.number))
length(unique(data_with_duplicates$Species))

#-multi asv one spec
#Remove the Multiple_ASV-One_species
data2 = HAB_data[ !(HAB_data$ASV.number %in% data_with_duplicates$ASV.number), ]

#All species that are duplicates in the remaining dataset
duplicate_species <- data2$Species[duplicated(data2$Species)]
data_with_duplicate_species <- data2[data2$Species %in% duplicate_species, ]

length(unique(data_with_duplicate_species$ASV.number))
length(unique(data_with_duplicate_species$Species))

#-One_ASV-One_Species
data3 = data2[ !(data2$ASV.number %in% data_with_duplicate_species$ASV.number), ]
nrow(data3)

