#############
#Phylogenetic tree
#############

#---Libraries
library(ape) 
library(phangorn) #
library(ggtree) #
library(RColorBrewer) #
library(ggplot2) #
library(DECIPHER)
library(ips)
library(knitr)
#library(msa)
library(Biostrings)
library(dplyr)

data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

V4data = data[, !grepl("Ex|FC|LC|V9", names(data))]
V4data <- V4data[rowSums(V4data[, 3:22] != 0) > 0, ]

#---Load, prepare and transform the data
ASV_df = V9data
#ASV_df = V4data

# ASV_df$Phylum = ASV_df$Class #to investigate at a class level
ASV_df_phylum = data.frame(ASV_df$ASV.number, ASV_df$Sequence, ASV_df$Phylum)
ASV_df_phylum = unique(ASV_df_phylum[c("ASV_df.ASV.number", "ASV_df.Sequence", "ASV_df.Phylum")])

length(unique(ASV_df_phylum$ASV_df.ASV.number))
length(unique(ASV_df_phylum$ASV_df.Phylum))
nrow(ASV_df_phylum)

# Find the duplicate ASV.number values and creat a dataframe without them.
#duplicates <- ASV_df_phylum$ASV_df.ASV.number[duplicated(ASV_df_phylum$ASV_df.ASV.number)]

#data_with_duplicates <- ASV_df_phylum[ASV_df_phylum$ASV_df.ASV.number %in% duplicates, ]

#length(unique(data_with_duplicates$ASV_df.ASV.number))
#length(unique(data_with_duplicates$ASV_df.Phylum)) #In this example, 11 ASVS still are a tossup between multiple phyla

#data_without_duplicates = ASV_df_phylum[!(ASV_df_phylum$ASV_df.ASV.number %in% duplicates), ]

#length(unique(data_without_duplicates$ASV_df.ASV.number))
#length(unique(data_without_duplicates$ASV_df.Phylum))
#nrow(data_without_duplicates) #nrow and 

# Assuming your data frame is named ASV_df_phylum
ASV_df_phylum <- ASV_df_phylum %>%
  group_by(ASV_df.ASV.number, ASV_df.Sequence) %>%
  mutate(ASV_df.Phylum = ifelse(n() > 1, "Ambiguous", ASV_df.Phylum)) %>%
  distinct() %>%
  ungroup()

which(ASV_df_phylum == "Ambiguous")
data_without_duplicates = ASV_df_phylum

#---Summary---Not used (used later in ASV composition including ambiguous reads.)
ASV_df = as.data.frame(data_without_duplicates[, 1:2])
#kable(sort(table(ASV_df$ASV_df.Phylum),decreasing = TRUE))
#kable(sort(table(ASV_df$ASV_df.Phylum)/length(ASV_df$ASV_df.Phylum)*100, decreasing = TRUE))


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
#Splitinfo <- ASV_df$ASV_df.Phylum[match(ggtree_tree$tip.label, ASV_df$ASV_df.ASV.number)]
Splitinfo <- data_without_duplicates$ASV_df.Phylum[match(ggtree_tree$tip.label, ASV_df$ASV_df.ASV.number)]
PhylumInfo = split(ggtree_tree$tip.label, Splitinfo)
PhylumTree = groupOTU(ggtree_tree, PhylumInfo) #this separates the ASVS by phylum
#saveRDS(PhylumTree, file = "PhylumTree2.RData")
#PhylumTree = readRDS("PhylumTree2.RData")


# Define a color palette using RColorBrewer
#color_palette <- brewer.pal(n = length(unique(ASV_df$ASV_df.Phylum)), name = "Dark2")
color_palette <- brewer.pal(n = length(unique(data_without_duplicates$ASV_df.Phylum)), name = "Dark2")

g = ggtree(PhylumTree, aes(color=group), layout='circular', branch.length = "none", key_glyph = "rect") + 
  #geom_tiplab(size=2, aes(angle=angle)) +
  scale_color_manual(values = color_palette) +
  labs(title = NULL,
       x = NULL,
       y = NULL) + 
  scale_color_manual(name="Phyla", values = color_palette, guide = guide_legend(override.aes = list(shape = 21, size = 1)))

g
