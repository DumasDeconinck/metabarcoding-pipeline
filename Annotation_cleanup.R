#########################################################
#Annotation data pipeline                               #
#Dumas Deconinck                                        #
#27/06/2024                                             #
#Version:  R version 4.3.2 (2023-10-31 ucrt) Eye Holes  #
#########################################################

#---libraries
library(dplyr)
library(seqRFLP)
library(worms)

#---Remove low read ASVs
abundancedata = read.csv("read_abundance_table.csv")
abundancedata <- abundancedata[,-1]

#---V4 and V9 filter
# Create a new dataframe without the columns to be removed
V4data = abundancedata[, !grepl("Ex|FC|LC|V9", names(abundancedata))]
V9data = abundancedata[, !grepl("Ex|FC|LC|V4", names(abundancedata))]
dadanonegs = abundancedata[, !grepl("Ex|FC|LC", names(abundancedata))]

#Only keep rows with counts
V4data <- V4data[rowSums(V4data[, 3:22] != 0) > 0, ]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]
dadanonegs = dadanonegs[rowSums(dadanonegs[, 3:42] != 0) > 0, ]

length(unique(abundancedata$ASV))
length(unique(V4data$ASV))
length(unique(V9data$ASV))
length(unique(dadanonegs$ASV))

#abundancedata = V4data
#abundancedata = V9data
abundancedata = dadanonegs

# Assuming the dataframe has n columns
# Keep columns 2 to (n-9) as the count columns
count_cols <- 3:(ncol(abundancedata)-9)

# Create a logical vector indicating rows with at least 2 counts in any of the count columns
keep_rows <- rowSums(abundancedata[, count_cols] >= 2) >= 1

# Subset the dataframe to keep only the rows that meet the criteria
abundancedata_nlr <- abundancedata[keep_rows, ]

nrow(abundancedata)
nrow(abundancedata_nlr)
nrow(abundancedata_nlr)/nrow(abundancedata)

length(unique(abundancedata_nlr$ASV))

#---Worms data 
#Remove NA species, we will blast those later
abundancedata_nlr_classified = na.omit(abundancedata_nlr)

length(unique(abundancedata_nlr_classified$ASV))

#Compare to the worms database.

w <- wormsbynames(gsub("_", " ", abundancedata_nlr_classified$Species), match = TRUE)
#w <- wormsbynames(abundancedata_nlr_classified$Species, match = TRUE)
#This does two extra things. 1) it checks if the species name from our blast results is in the database.
##if it is, and it's accepted, the name will stay the same
##If it is, and it's not accepted, an accepted name is provided
##If it isn't in the worms database, it isn't a recognized marine species, and thus can't be marine phytoplankton. We can ignore these.
##additionally, it provides the fully taxonomy for the accepted species names.

PR2_dada_worms = data.frame(abundancedata_nlr_classified[,1:(length(abundancedata_nlr_classified)-9)],  w$kingdom, w$phylum, w$class, w$order, w$family, w$genus, w$valid_name)
colnames(PR2_dada_worms) = c(names(abundancedata_nlr_classified)[1:(length(abundancedata_nlr_classified)-9)], "Kindom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(PR2_dada_worms)[1] = "ASV number"
PR2_dada_worms = na.omit(PR2_dada_worms)
write.csv(PR2_dada_worms, file = "PR2_dada_worms.csv")


#---Blast species NA
#extract unclassified species
unclassified = subset(abundancedata_nlr, is.na(abundancedata_nlr$Species))

#Remove ASVs that aren't even Eukaryota or Bacteria (Bacteria because of cyanobacteria)
unclassified = subset(unclassified, Domain %in% c("Eukaryota", "Bacteria"))

#Turn the AVS and sequences into a fasta file.
df.fasta = dataframe2fas(unclassified[1:2], file="PR2_dada2_unclassified.fasta")

#Run the line below on blastn in linux/command centre.
#Apps/ncbi-blast-2.15.0+/bin/blastn -query PR2_dada2_unclassified.fasta -db SSU_eukaryote_rRNA -out PR2_dada2_blasted.out -evalue 1e-6 -outfmt "6 qseqid sscinames bitscore score qcovs evalue pident length" -max_target_seqs 5 -perc_identity 99 
PR2_dada2_blasted <- read.delim("~/Preliminary_Phytoplankton_Metabarcoding/PR2_dada2_blasted.out", header=FALSE)
colnames(PR2_dada2_blasted) = c("ASV number", "Species", "bitscore", "score", "query coverage", "evalue", "per. ident", "length")
PR2_dada2_blasted = PR2_dada2_blasted[,1:2]
PR2_dada2_blasted = unique(PR2_dada2_blasted)

length(unique(PR2_dada2_blasted$`ASV number`))

#sum(unique(PR2_dada2_blasted$`ASV number`) %in% unique(unclassified$ASV), na.rm = TRUE)

#This will link the species names to the worms database
w <- wormsbynames(PR2_dada2_blasted$Species, match = TRUE)
#This does two extra things. 1) it checks if the species name from our blast results is in the database.
##if it is, and it's accepted, the name will stay the same
##If it is, and it's not accepted, an accepted name is provided
##If it isn't in the worms database, it isn't a recognized marine species, and thus can't be marine phytoplankton. We can ignore these.
##additionally, it provides the fully taxonomy for the accepted species names.

PR2_dada2_blasted_worms = data.frame(PR2_dada2_blasted$`ASV number`, w$kingdom, w$phylum, w$class, w$order, w$family, w$genus, w$valid_name)
colnames(PR2_dada2_blasted_worms) = c("ASV number", "Kindom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#match ASV numbers back with the read counts and sequences
matchdf = abundancedata_nlr[1:(length(abundancedata_nlr)-9)]

PR2_dada2_blasted_worms_with_matches <- left_join(PR2_dada2_blasted_worms, matchdf, by = c("ASV number" = "ASV"))

PR2_dada2_blasted_worms_with_matches_reordered <- PR2_dada2_blasted_worms_with_matches %>%
  select(1, 9:ncol(PR2_dada2_blasted_worms_with_matches), 2:8)

PR2_dada2_blasted_worms = PR2_dada2_blasted_worms_with_matches_reordered
PR2_dada2_blasted_worms = na.omit(PR2_dada2_blasted_worms)

write.csv(PR2_dada2_blasted_worms, file = "PR2_dada2_blasted_worms.csv")

#---Stick the PR2 dada classified dataframe and the blasted data frame back together

Classified_species = rbind(PR2_dada_worms, PR2_dada2_blasted_worms)
Classified_species$number = substr(Classified_species$`ASV number`, start = 5, stop = nchar(Classified_species$`ASV number`))
Classified_species = Classified_species[order(as.numeric(Classified_species$number)),]
Classified_species$number = NULL

write.csv(Classified_species, file = "Classified_species.csv", row.names=FALSE)

#---Keep only phytoplankton 
Classified_species_phyto = Classified_species  %>%
  filter(Phylum %in% c("Ochrophyta", "Cryptophyta", "Chlorophyta", "Rhodophyta", "Haptophyta", "Cyanobacteria") |
           Class %in% c("Dinophyceae", "Dinoflagellata incertae sedis", "Bacillariophyceae") |
           Family == "Katablepharidaceae"
         )

#write.csv(Classified_species_phyto, file = "Classified_species_phyto.csv", row.names=FALSE)
colnames(Classified_species_phyto)

length(unique(Classified_species_phyto$`ASV number`))

Classified_species_phyto2 <- Classified_species_phyto %>%
  filter(str_count(Species, " ") >= 1)

length(unique(Classified_species_phyto2$`ASV number`))

#Remove ASVs that snuck in that were only identified at the genus level.
write.csv(Classified_species_phyto2, file = "Classified_species_phyto.csv", row.names=FALSE)

