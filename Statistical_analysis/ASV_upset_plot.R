#############
#Upsets
#############

#---library
library(UpSetR)
library(ggupset)
library(ComplexUpset)
library(ggplot2)

#---data
data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

#Keep onll Phylum and ASV number
V9data = V9data[, !grepl("Sequence|Kindom|Domain|Phylum|Class|Order|Family|Genus|Species", names(V9data))]
V9data <- V9data %>%
  distinct(ASV.number, .keep_all = TRUE)

colnames(V9data) =  c(colnames(V9data)[1],sub("^[^.]+\\.", "", colnames(V9data)[-1]))

#Transform into a list
V9data_list <- vector("list", ncol(V9data) - 1)
names(V9data_list) <- names(V9data)[-1] # Use the sample names as the list element names

for (i in 2:ncol(V9data)) { # Iterate over the sample columns (skip the first 'ASV.number' column)
  sample_name <- names(V9data)[i]
  sample_asv <- V9data$ASV.number[V9data[, i] > 0] # Get the ASV numbers with non-zero counts for that sample
  V9data_list[[sample_name]] <- as.vector(sample_asv) # Add the list of ASV numbers to the list, using the sample name as the element name
}

sample_order <- c("S1S", "S2S", "S3S", "S4S", "S5S", "S6S",
                  "S1B", "S2B", "S3B", "S4B", "S5B", "S6B",
                  "P1S", "P2S", "P3S", "P4S",
                  "P1B", "P2B", "P3B", "P4B")

V9data_list <- vector("list", length(sample_order))
names(V9data_list) <- sample_order

for (sample_name in sample_order) {
  sample_asv <- V9data$ASV.number[V9data[, sample_name] > 0] # Get the ASV numbers with non-zero counts for that sample
  V9data_list[[sample_name]] <- as.vector(sample_asv) # Add the list of ASV numbers to the list, using the sample name as the element name
}

#---upset plot
UpSetR::upset(fromList(V9data_list), order.by = "freq", keep.order = T,
              nsets= 20, nintersects = 44,
              mainbar.y.label = "Number of shared ASVs", 
              sets.x.label = "ASVs per Sample", show.numbers = "yes",
              sets = rev(sample_order), point.size = 1, text.scale = 0.8, matrix.dot.alpha = 0,
              line.size = 1)

#---intersects Regions
Reduce(intersect, V9data_list)

# Create a new list for S and P (West and East)
s_list <- V9data_list[names(V9data_list)[startsWith(names(V9data_list), "S")]]
p_list  = V9data_list[names(V9data_list)[startsWith(names(V9data_list), "P")]]

# Find the common ASV numbers
s_common_asvs <- Reduce(intersect, s_list)
p_common_asvs <- Reduce(intersect, p_list)

setdiff(s_common_asvs, unique(unlist(p_list))) #list of ASVS in every western sample, but in none of the eastern samples
setdiff(p_common_asvs, unique(unlist(s_list))) #ASV in every eastern sample, but none of the western samples.

#---intersects Depth
Reduce(intersect, V9data_list)

# Create a new list for S and P (West and East)
b_list <- V9data_list[names(V9data_list)[endsWith(names(V9data_list), "B")]]
s_list  = V9data_list[names(V9data_list)[endsWith(names(V9data_list), "S")]]


# Find the common ASV numbers
b_common_asvs <- Reduce(intersect, b_list)
s_common_asvs <- Reduce(intersect, s_list)

setdiff(b_common_asvs, unique(unlist(s_list))) #list of ASVS in every western sample, but in none of the eastern samples
setdiff(s_common_asvs, unique(unlist(b_list))) #ASV in every eastern sample, but none of the western samples.



#---simplified
# Define a function to rename the columns
rename_columns <- function(column_name) {
  if (startsWith(column_name, "P") && endsWith(column_name, "B")) {
    return("East_Bottom")
  } else if (startsWith(column_name, "P") && endsWith(column_name, "S")) {
    return("East_Surface")
  } else if (startsWith(column_name, "S") && endsWith(column_name, "B")) {
    return("West_Bottom")
  } else if (startsWith(column_name, "S") && endsWith(column_name, "S")) {
    return("West_Surface")
  } else {
    return(column_name) # Keep the original name for columns that don't fit the criteria
  }
}

# Apply the renaming function to the column names
colnames(V9data) <- sapply(colnames(V9data), rename_columns)

library(dplyr)
library(tidyr)

long_data <- V9data %>%
  pivot_longer(cols = -ASV.number, names_to = "Sample", values_to = "Value")
long_sum = long_data %>%
  group_by(ASV.number,Sample) %>% summarize(Total = sum(Value))
V9data = long_sum %>% pivot_wider(names_from = "Sample",values_from="Total")

#Transform into a list
V9data_list <- vector("list", ncol(V9data) - 1)
names(V9data_list) <- names(V9data)[-1] # Use the sample names as the list element names

for (i in 2:ncol(V9data)) { # Iterate over the sample columns (skip the first 'ASV.number' column)
  sample_name <- names(V9data)[i]
  sample_asv <- V9data$ASV.number[V9data[, i] > 0] # Get the ASV numbers with non-zero counts for that sample
  V9data_list[[sample_name]] <- as.vector(sample_asv) # Add the list of ASV numbers to the list, using the sample name as the element name
}

sample_order <- c("West_Surface", "West_Bottom", "East_Surface", "East_Bottom")

V9data_list <- vector("list", length(sample_order))
names(V9data_list) <- sample_order

for (sample_name in sample_order) {
  sample_asv <- V9data$ASV.number[V9data[, sample_name] > 0] # Get the ASV numbers with non-zero counts for that sample
  V9data_list[[sample_name]] <- as.vector(sample_asv) # Add the list of ASV numbers to the list, using the sample name as the element name
}

#---upset plot
UpSetR::upset(fromList(V9data_list), order.by = "freq", keep.order = T,
              nsets= 20, nintersects = NA,
              mainbar.y.label = "Number of shared ASVs", 
              sets.x.label = "ASVs per region", show.numbers = "yes",
              sets = rev(sample_order))

#---intersects
Reduce(intersect, V9data_list)

# Create a new list with only the 'S' names
s_list <- V9data_list[names(V9data_list)[startsWith(names(V9data_list), "West")]]
p_list  = V9data_list[names(V9data_list)[startsWith(names(V9data_list), "East")]]

# Find the common ASV numbers
s_common_asvs <- Reduce(intersect, s_list)
p_common_asvs <- Reduce(intersect, p_list)

setdiff(s_common_asvs, unique(unlist(p_list))) #list of ASVS in every western sample, but in none of the eastern samples
setdiff(p_common_asvs, unique(unlist(s_list))) #ASV in every eastern sample, but none of the western samples.
