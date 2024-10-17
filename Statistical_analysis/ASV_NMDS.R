#############
#NMDS
#############
library(vegan)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pairwiseAdonis)


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

#---nmds
#modify the data
df_numeric = as.data.frame(t(V9data))
colnames(df_numeric) = df_numeric[1, ]
df_numeric = df_numeric[-1,]

df_numeric <- as.data.frame(lapply(df_numeric, as.numeric))

NMDS3 <- metaMDS(df_numeric, k = 2, trymax = 100, trace = F, autotransform = FALSE, distance="bray")

##GGPLOT
data.scores <- as.data.frame(scores(NMDS3, "sites")) #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- colnames(V9data[-1])  # create a column of site names, from the rownames of data.scores
head(data.scores)  #look at the data
rownames(data.scores) = data.scores$site

species.scores <- as.data.frame(scores(NMDS3, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)

#ef.scores = as.data.frame(scores(ef, "vectors"))
#ef.scores$factors = rownames(ef.scores)
#head(ef.scores)

#geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2),size=3) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

data.scores$Region <- ifelse(substr(rownames(data.scores), 1, 1) == "P", "East",
                           ifelse(substr(rownames(data.scores), 1, 1) == "S", "West", "Unknown"))
data.scores$Depth <- ifelse(substr(rownames(data.scores), nchar(rownames(data.scores)), nchar(rownames(data.scores))) == "S", "Surface",
                          ifelse(substr(rownames(data.scores), nchar(rownames(data.scores)), nchar(rownames(data.scores))) == "B", "Bottom", "Unknown"))

#Prepare polygon
grp.a <- subset(data.scores, Region == "East")
grp.b <- subset(data.scores, Region == "West")

grp.a <- data.scores[data.scores$Region == "East", ][chull(data.scores[data.scores$Region == 
                                                                               "East", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$Region == "West", ][chull(data.scores[data.scores$Region == 
                                                                               "West", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

hull.data <- rbind(grp.a, grp.b)  # combine grp.a and grp.b
hull.data

color_palette <- brewer.pal(2, "Dark2")  # Adjust the number of colors as needed

plot = ggplot() +
  geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, fill = Region, group = Region), alpha = 0.3) +
  ggrepel::geom_text_repel(data = data.scores, aes(x = NMDS1, y = NMDS2, label = site), alpha = 0.5, nudge_y = -0.05) +
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, colour = Depth, shape = Depth),size=2) + # add the point markers
  scale_fill_manual(values = color_palette) +
  scale_colour_manual(values = c(brewer.pal(4, "Dark2")[3], brewer.pal(4, "Dark2")[4])) +
  scale_shape_manual(values = c(1:12)) +
  labs(x = "NMDS1", y = "NMDS2") +
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
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA))
plot

NMDS3$stress #< 0.10	Good ordination with no real risk of drawing false inferences and 0.05 almost perfect
gof = goodness(object = NMDS3) #mlarge value would indicate a sample ruines the goodness of fit.
gof
plot(NMDS3, display = "sites", type = "none")
points(NMDS3, display = "sites", cex = 2*gof/mean(gof)) #same plot but larger circles account for more of the variation
plot(NMDS3$diss, NMDS3$dist)#Evaluate fit, Sheppard plot 
stressplot(object = NMDS3, lwd = 5) #non metric fit of 0.997 and linear fit of 0.992

#https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/nmds/
sample_coords <- as.data.frame(NMDS3$points)
species_scores <- as.data.frame(NMDS3$species)
NMDS3$dist
prop_var_explained <- 1 - (NMDS3$stress^2)
prop_var_explained

#PERMANOVA
permanova <- adonis2(df_numeric ~ Region * Depth, data = data.scores, method = "bray", permutations = 999)
permanova

#PERMANOVA only EAST
permanova <- adonis2(df_numeric[1:8,] ~ Depth, data = data.scores[1:8,], method = "bray", permutations = 999)
permanova

#Pairwise permanova (don't see a difference)
pairpermanova = pairwise.adonis2(df_numeric ~ Region + Depth, data = data.scores, method = "bray", permutations = 999)
pairpermanova

# Calculate the distance matrix
d <- vegdist(df_numeric, method = "bray")

# Perform the BETADISPER analysis
betadisper_model <- betadisper(d, data.scores$Region, type = "centroid")

# Test for differences in dispersion between groups
anova(betadisper_model)

# Perform pairwise comparisons of dispersion
TukeyHSD(betadisper_model)


svg("NMDS.svg", width = 760, height = 547, onefile = FALSE)
print(plot)
dev.off()
