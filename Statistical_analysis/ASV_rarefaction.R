#############
#Rarefaction
#############

#---Libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)
library(purrr)

#---load the data
data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

#keep only ASV
V9data = V9data[, !grepl("Sequence|Kindom|Domain|Phylum|Class|Order|Family|Genus|Species", names(V9data))]

#Remove duplicate ASVs
V9data <- V9data %>%
  distinct(ASV.number, .keep_all = TRUE)

#---Quick rarefaction
raredf = t(V9data[,-1])
p = rarecurve(raredf, step = 10, cex = 0.6)

#---Turn into ggplot
# rarefraction curves
out <- rarecurve(t(V9data[,-1]), step = 10, sample = 400)
names(out) <- colnames(V9data[-1])
as_tibble_rc <- function(x){
 # convert rarecurve() output to data.frame
 # bind_rows doesn't work because items are different lengths
 # also need to extract sample sizes from attribute
 # Allocate result dataframe
 nsamples <- map_int(x, length)
 total_samples <- sum(nsamples)
 if(!is.null(names(x))){
   sites <- names(x)
 } else {
   sites <- as.character(1:length(nsamples))
 }
 result <- data_frame(Site = rep("", total_samples),
                      Sample_size = rep(0, total_samples),
                      Species = rep(0, total_samples))
 start <- 1
 for (i in 1:length(nsamples)){
   result[start:(start + nsamples[i]-1), "Site"] <- sites[i]
   result[start:(start + nsamples[i]-1), "Sample_size"] <- attr(x[[i]],
"Subsample")
   result[start:(start + nsamples[i]-1), "Species"] <- x[[i]]
   start <- start + nsamples[i]
 }
 result
}
out <- as_tibble_rc(out)

# add grouping variable
sitedata <- data_frame(Site = as.character(1:50),
                      Type = sample(LETTERS[1:2], 50, replace = TRUE))
alldata <- left_join(out, sitedata, by = "Site")

#ggplot
ggplot(data = alldata,
      mapping = aes(x = Sample_size, y = Species, color = Site)) +
 geom_line()

last_points <- alldata %>%
  group_by(Site) %>%
  slice(n())

ggplot(data = alldata, mapping = aes(x = Sample_size, y = Species, group = Site)) +
  geom_line() +
  ggrepel::geom_label_repel(
    data = last_points,
    aes(x = Sample_size, y = Species, label = Site),
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "grey50",
    size = 3,
    direction = "x",
    hjust = 1.1,
    max.overlaps = 20
  ) +
  theme_classic() +
  ylab("ASVs") +
  xlab("Reads") +
  xlim(min(alldata$Sample_size) * 0.9, max(alldata$Sample_size) * 1.1)

#---Transform the dataframe to long format
long_data <- V9data %>%
  pivot_longer(cols = -ASV.number, names_to = "Sample", values_to = "Value")

# Reorder the columns
long_data <- long_data %>%
  select(Sample, ASV.number, Value)

# View the resulting dataframe
long_data
V9data = long_data

ggplot(V9data, aes(fill=ASV.number, y=Value, x=Sample)) + #quick visualisation
  geom_bar(position="stack", stat="identity") +
  ggtitle("V9 Relative phytoplankton ASV richness") +
  theme(legend.position = "none")


#---rarefaction manual Looks good, takes way too long
###example of list of identified ASVs
  
    x = subset(V9data, V9data$Sample == unique(V9data$Sample)[1])
    ASV_list = rep(x$ASV.number, x$Value)
    
    #This is one rarefied example. Essentially, from the species list, you subsample 20. The larger the subsamble, the more accurate it should be.
    length(unique(sample(ASV_list, 20, replace = F)))
    
    #Then, you want to do this for larger subsizes, here we define how much the subsamble increases per step.
    step_size = 10 #define the step size (the size of the rarefied dataset)
    
    num_iterations <- 2# Define the number of iterations (This is to average out the number of species detect per subsample to avoid stochastic effects, essentially smoothes out the graph)
    
    result_list <- list() # Create an empty list to store the results
    
    for (iteration in 1:num_iterations) {
      Identified_Species <- c()
      Sample_Size <- c()
      i <- 0
      step <- step_size
      
      while (i < (length(ASV_list)-step_size)){
        i <- i + step
        Sample_Size <- append(Sample_Size, i)
        Identified_Species <- append(Identified_Species, length(unique(sample(ASV_list, i, replace = FALSE))))
      }
      
      raredf <- data.frame("Species" = Identified_Species, "Sample size" = Sample_Size)
      result_list[[iteration]] <- raredf
    }
    
    # Calculate the averages for raredf
    avg_raredf <- Reduce("+", result_list) / num_iterations
    
    ggplot(data = avg_raredf, aes(x = Sample.size, y = Species)) +
      geom_line() +
      geom_text(aes(x = max(Sample.size), y = max(Species), label = unique(V9data$Sample)[1]), 
                hjust = 0, vjust = 0.5, size = 2.5) +
      scale_x_continuous(limits = c(0, max(avg_raredf$Sample.size) * 1.1)) +
      theme_classic() +
      xlab("Reads") +
      ylab("ASVs")
    

# Initialize an empty dataframe to store the results
all_raredf <- data.frame()

for (k in 1:length(unique(V9data$Sample))) {
  x <- subset(V9data, V9data$Sample == unique(V9data$Sample)[k])
  ASV_list <- rep(x$ASV.number, x$Value)
  
  # This is one rarefied example. Essentially, from the species list, you subsample 20. The larger the subsample, the more accurate it should be.
  length(unique(sample(ASV_list, 20, replace = F)))
  
  # Then, you want to do this for larger subsizes, here we define how much the subsample increases per step.
  step_size <- 10  # Define the step size (the size of the rarefied dataset)
  num_iterations <- 1  # Define the number of iterations (This is to average out the number of species detect per subsample to avoid stochastic effects, essentially smoothes out the graph)
  
  for (iteration in 1:num_iterations) {
    Identified_Species <- c()
    Sample_Size <- c()
    i <- 0
    step <- step_size
    
    while (i < (length(ASV_list) - step_size)) {
      i <- i + step
      Sample_Size <- append(Sample_Size, i)
      Identified_Species <- append(Identified_Species, length(unique(sample(ASV_list, i, replace = FALSE))))
    }
    
    raredf <- data.frame("Sample" = unique(V9data$Sample)[k], "Sample.size" = Sample_Size, "Species" = Identified_Species)
    all_raredf <- rbind(all_raredf, raredf)
    
    print(paste(unique(V9data$Sample)[k], "Iteration", iteration))
  }
  print(paste("finished", unique(V9data$Sample)[k]))
}



all_raredf %>%
  mutate(label = if_else(Sample.size == max(Sample.size), as.character(Sample), NA_character_)) %>%
  ggplot(aes(x = Sample.size, y = Species, color = Sample)) +
    geom_line() +
    scale_x_continuous(limits = c(0, max(all_raredf$Sample.size) * 1.1)) +
    theme_classic() +
    xlab("Reads") +
    ylab("ASVs") +
    labs(color = "Sample") +
    geom_label_repel(aes(label = label),
                   nudge_x = 1,
                   na.rm = TRUE)

library(directlabels)

all_raredf = read.csv("all_raredf.csv", sep = ";")

ggplot(all_raredf, aes(x = Sample.size, y = Species, colour = Sample)) +
  geom_line() +
  scale_colour_discrete(guide = "none") +
  scale_x_continuous(limits = c(0, max(all_raredf$Sample.size) * 1.1)) +
  #scale_x_log10() +
  geom_dl(aes(label = Sample), method = list(cex = 0.6, dl.combine("last.points"))) +
  theme_classic() +
  xlab("Reads") +
  ylab("ASVs")



