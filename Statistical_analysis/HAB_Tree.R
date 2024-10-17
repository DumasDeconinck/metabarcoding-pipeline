######
#HAB tree
#####

#---libraries
library(tidyr)
library(dplyr)
library(readxl)

#---Libraries
library(devtools)
library(ape)
library(phangorn)
library(ggtree)
library(RColorBrewer)
library(ggplot2)
library(DECIPHER)
library(ips)
library(knitr)

#---data
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

#Keep only Species and ASV number
V9data = V9data[, !grepl("Kindom|Domain|Class|Order|Family|Genus", names(V9data))]

# Group by the ASV number and check the Species values
new_data <- V9data %>%
  group_by(ASV.number, Species, Sequence, Phylum) %>%
  summarize(across(where(is.numeric), sum), .groups = "drop") %>%
  group_by(ASV.number) %>%
  mutate(
    unique_species = n_distinct(Species)
  ) %>%
  # If there are multiple Species, set the Species to "Ambiguous"
  mutate(
    Species = ifelse(unique_species > 1, "Ambiguous", Species)
  ) %>%
  # Detect ambiguous Phyla
  group_by(ASV.number, Species) %>%
  mutate(
    ambiguous_phylum = if(n_distinct(Phylum) > 1) "Ambiguous" else Phylum
  ) %>%
  # Keep only the first row for each ASV.number
  slice(1)

new_data <- V9data %>%
  group_by(ASV.number, Species, Sequence, Phylum) %>%
  summarize(across(where(is.numeric), sum), .groups = "drop") %>%
  group_by(ASV.number) %>%
  mutate(
    unique_species = n_distinct(Species)
  ) %>%
  # If there are multiple Species, set the Species to "Ambiguous"
  mutate(
    Species = ifelse(unique_species > 1, "Ambiguous", Species)
  ) %>%
  # Detect ambiguous Phyla
  group_by(ASV.number, Species) %>%
  mutate(
    ambiguous_phylum = if(n_distinct(Phylum) > 1) "Ambiguous" else Phylum
  ) %>%
  # Keep only the unique rows
  ungroup() %>%
  distinct()

#Creating a tag
new_data$Species = paste(new_data$ASV.number, new_data$Species, sep = "-")
new_data$ASV.number = NULL

V9data = new_data %>% select(Species, Sequence, Phylum)

ASV_df = as.data.frame(V9data)

#---Make the fasta
ASV_fasta = seqRFLP::dataframe2fas(ASV_df[,1:2], "df.fasta")
ASV_fasta = readDNAStringSet("df.fasta")

#---align
Testdata = msa::msaClustalW(ASV_fasta, verbose = TRUE)
ASV_phydat = msa::msaConvert(Testdata, "phangorn::phyDat")

#---Create a distnce matrix
mt = modelTest(ASV_phydat)
dna_dist = dist.ml(ASV_phydat)

#---create and unroot tree
UPGMA_test = upgma(dna_dist)
NJ_test = NJ(dna_dist)
plot(UPGMA_test, main="UPGMA")
plot(NJ_test, main = "Neighbor Joining")

#---Find the best fit
parsimonUPGMA = parsimony(data = ASV_phydat, tree = UPGMA_test)
parsimonNJ = parsimony(data = ASV_phydat, tree = NJ_test)
min(parsimonUPGMA, parsimonNJ)

if(parsimonUPGMA < parsimonNJ){BestFit = UPGMA_test} else {BestFit = NJ_test}

#test_optim <- optim.parsimony(BestFit, ASV_phydat)
#test_pratchet <- pratchet(ASV_phydat)
#plot(test_optim)
#plot(test_pratchet)

#---Maximum Likelihood
fit <- pml(BestFit, ASV_phydat)
print(fit)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0)) #For a final graph, I recommend setting bs to 1000 (bootstrap), but it takes very long
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")

#---Transform into ggtree
# Convert the maximum likelihood tree to a ggtree-compatible format
ggtree_tree <- fitJC$tree
ggtree_tree = plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")

# Plot the maximum likelihood tree using ggtree
#PhylumInfo = split(ggtree_tree$tip.label, ASV_df$ASV_df.Phylum) #THIS DOES NOT WORK, UPDATED LINE BELOW
Splitinfo <- ASV_df$Phylum[match(ggtree_tree$tip.label, ASV_df$Species)]
PhylumInfo = split(ggtree_tree$tip.label, Splitinfo)
PhylumTree = groupOTU(ggtree_tree, PhylumInfo) #this separates the ASVS by phylum
#saveRDS(PhylumTree, file = "PhylumTree_HAB.RData")
#PhylumTree = readRDS("PhylumTree_HAB.RData")


# Define a color palette using RColorBrewer
color_palette <- brewer.pal(n = length(unique(ASV_df$Phylum))+1, name = "Dark2")
color_palette = c(#"white", 
                  brewer.pal(8, "Dark2")[3], 
                  brewer.pal(8, "Dark2")[4], 
                  brewer.pal(8, "Dark2")[5], 
                  brewer.pal(8, "Dark2")[6], 
                  brewer.pal(8, "Dark2")[7])

g = ggtree(PhylumTree, aes(color=group), key_glyph = "rect", branch.length = 0.5) + 
  geom_tiplab(size=2.5, align=TRUE) +
  scale_color_manual(values = color_palette) +
  labs(title = NULL,
       x = NULL,
       y = NULL) + 
  xlim(0, 7) +
  scale_color_manual(name="Phyla", values = color_palette, guide = guide_legend(override.aes = list(shape = 21, size = 1)))
g + theme_tree2()

g = ggtree(PhylumTree, aes(color=group),layout='circular', key_glyph = "rect") + 
  geom_tiplab(size=2, align=TRUE) +
  scale_color_manual(values = color_palette) +
  labs(title = NULL,
       x = NULL,
       y = NULL) + 
  xlim(0, 1.2) +
  scale_color_manual(name="Phyla", values = color_palette, guide = guide_legend(override.aes = list(shape = 21, size = 1)))
  
g 

get_taxa_name(g)


