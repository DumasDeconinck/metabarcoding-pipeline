#############
#Diversity Indices
#############

###Libraries
library(dplyr)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(tidyr)
library(fossil)
library(ggpubr)
###

#---load the data
data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

#keep only ASV
V9data = V9data[, !grepl("Sequence|Kindom|Domain|Phylum|Class|Order|Family|Genus|Species", names(V9data))]

#Remove duplicate ASVs
V9data <- V9data %>%
  distinct(ASV.number, .keep_all = TRUE)

#---Transform the dataframe to long format
long_data <- V9data %>%
  pivot_longer(cols = -ASV.number, names_to = "Sample", values_to = "Value")

# Reorder the columns
long_data <- long_data %>%
  select(Sample, ASV.number, Value)

Samples = colnames(V9data)[-1]

#---Chao distances
#Function
calculate_chao1 <- function(x) {# Define a function to calculate Chao1 diversity
  observed <- sum(x > 0)
  singletons <- sum(x == 1)
  doubletons <- sum(x == 2)
  chao1_estimate <- observed + (singletons*(singletons-1))/(2*(doubletons+1))
  return(chao1_estimate)
}

Dist_Chao1 <- long_data %>% # Calculate Chao1 diversity for each sample type
  group_by(Sample) %>%
  summarize(Chao1_Diversity = calculate_chao1(Value))

Dist_Chao1# View the resulting diversity values

#Add to dataframe
Richness = data.frame(Samples, "Chao1" = Dist_Chao1[2])

#---Shannon Diversity
diversity_values <- long_data %>% # Calculate Shannon diversity for each sample type
  group_by(Sample) %>%
  summarize(Shannon_Diversity = vegan::diversity(Value, index = "shannon"))

diversity_values #View the resulting diversity values

Richness = data.frame(Richness, diversity_values[2])

#---ACE
#Formula
calculate_ace <- function(x) {
  S_abund <- sum(x >= 10)
  S_rare <- sum(x < 10)
  N_rare <- sum(x[x < 10])
  
  a <- table(x[x >= 10])
  a1 <- a[1]
  
  C_ace <- 1 - a1/N_rare
  gamma_squared <- max(S_rare/C_ace * sum((1:length(a)) * (1:length(a) - 1) * a) / (N_rare/(N_rare - 1)) - 1, 0)
  
  ACE_estimate <- S_abund + S_rare / C_ace + a1 / C_ace * gamma_squared
  return(ACE_estimate)
}

diversity_values <- long_data %>%
  group_by(Sample) %>%
  summarize(ACE_Diversity = ACE(Value))

diversity_values

Richness = data.frame(Richness, diversity_values[2])

#---Simpson Diversity
diversity_values <- long_data %>%
  group_by(Sample) %>%
  summarize(Simpson_Diversity = vegan::diversity(Value, index = "simpson"))
diversity_values

Richness = data.frame(Richness, diversity_values[2])


#---Pielou's Eveness
#Function
calculate_pielou_evenness <- function(x) {
  n <- sum(x)  # Total number of individuals
  p <- x / n  # Proportions of each species
  p[p == 0] <- NA  # Replace 0s with NA
  ln_p <- log(p)  # Natural logarithm of proportions
  H <- -sum(p * ln_p, na.rm = TRUE)  # Shannon entropy (excluding NAs)
  J <- H / log(length(x))  # Pielou's evenness
  return(J)
}

evenness_values <- long_data %>%
  group_by(Sample) %>%
  summarize(Pielou_Evenness = calculate_pielou_evenness(Value))

evenness_values

Richness = data.frame(Richness, evenness_values[2])

#---Good's Coverage
#Function
calculate_goods_coverage <- function(x) {
  F1 <- sum(x == 1)  # Number of singleton OTUs
  N <- sum(x)  # Total number of individuals
  coverage <- 1 - (F1 / N)  # Good's coverage
  return(coverage)
}

coverage_values <- long_data %>%
  group_by(Sample) %>%
  summarize(Goods_Coverage = calculate_goods_coverage(Value))

coverage_values

Richness = data.frame(Richness, coverage_values[2])

## Create a Graph with the different indices!
Richness
colnames(Richness)= c("Samples", "Chao1", "Shannon", "ACE", "Simpson", "Pielou Eveness", "Good's Coverage")
Richness_long = melt(Richness, id.vars="Samples")

# Determine the number of variables
num_variables <- length(unique(Richness_long$variable))

# Define a color palette that can update with more variables
color_palette <- brewer.pal(num_variables, "Dark2")

# Convert "Samples" to a factor with the desired order
Richness_long$Samples <- factor(Richness_long$Samples, levels = unique(Richness$Samples))

# Draw the graph with improved aesthetics
ggplot(data = Richness_long) +
  geom_line(aes(x = Samples, y = value, colour = variable, group = variable),
            size = 1.2) +
  geom_point(aes(x = Samples, y = value, colour = variable, shape = variable),
             size = 3) +
  scale_y_log10()

Richness_long$Primer <- sub("\\..*", "", Richness_long$Samples)
Richness_long$Station <- sub("^[^.]+\\.", "", Richness_long$Samples)
Richness_long$Region <- ifelse(substr(Richness_long$Station, 1, 1) == "P", "East",
                               ifelse(substr(Richness_long$Station, 1, 1) == "S", "West", "Unknown"))
Richness_long$Depth <- ifelse(substr(Richness_long$Station, nchar(Richness_long$Station), nchar(Richness_long$Station)) == "S", "Surface",
                              ifelse(substr(Richness_long$Station, nchar(Richness_long$Station), nchar(Richness_long$Station)) == "B", "Bottom", "Unknown"))

Richness_East = subset(Richness_long, Region == "East")
Richness_East_Bottom = subset(Richness_East, Depth == "Bottom")
Richness_East_Surface = subset(Richness_East, Depth == "Surface")

Richness_West = subset(Richness_long, Region == "West")
Richness_West_Bottom = subset(Richness_West, Depth == "Bottom")
Richness_West_Surface = subset(Richness_West, Depth == "Surface")

Richness_long2 = rbind(Richness_West_Surface, Richness_West_Bottom, Richness_East_Surface, Richness_East_Bottom)

#Keep order
Richness_long2$Station <- factor(Richness_long2$Station, levels = unique(Richness_long2$Station))

line = ggplot(data = Richness_long2) +
  geom_line(aes(x = Station, y = value, colour = variable, group = variable),
            size = 1.2) +
  geom_point(aes(x = Station, y = value, colour = variable, shape = variable),
             size = 3) +
  scale_y_log10() +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, angle = 45, hjust = 1)
  ) +
  ylab(NULL) +
  xlab("Sample")
line

#Add a combined column for boxplots
Richness_long2$Combined = paste(Richness_long2$Region, Richness_long2$Depth)

#order
Richness_long2$Combined <- factor(Richness_long2$Combined, levels = unique(Richness_long2$Combined))

box = ggplot(data = Richness_long2) +
  geom_boxplot(aes(x = Combined, y = value, colour = variable)) +
  scale_y_log10() +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
  ) +
  ylab(NULL) +
  xlab("Sample")
box

ggarrange(
  line, box,
  align = "h", labels = c("A", "B"),
  common.legend = TRUE,
  legend = "bottom"
)

#---Data
# Check standard deviation and range for each variable
variable_summary <- data.frame(
  variable = unique(Richness_long$variable),
  std_dev = numeric(length(unique(Richness_long$variable))),
  value_range = numeric(length(unique(Richness_long$variable)))
)

for (i in seq_along(variable_summary$variable)) {
  var <- variable_summary$variable[i]
  var_data <- subset(Richness_long, variable == var)
  variable_summary$std_dev[i] <- sd(var_data$value)
  variable_summary$value_range[i] <- max(var_data$value) - min(var_data$value)
}

# Filter out variables with no variance
variable_summary <- variable_summary[variable_summary$std_dev > 0, ]

# Perform t-tests for each variable
t_test_results <- data.frame(
  variable = variable_summary$variable,
  region_pvalue = numeric(nrow(variable_summary)),
  depth_pvalue = numeric(nrow(variable_summary))
)

for (i in seq_along(t_test_results$variable)) {
  var <- t_test_results$variable[i]
  var_data <- subset(Richness_long, variable == var)
  
  # Region effect
  region_t_test <- t.test(value ~ Region, data = var_data)
  t_test_results$region_pvalue[i] <- region_t_test$p.value
  
  # Depth effect
  depth_t_test <- t.test(value ~ Depth, data = var_data)
  t_test_results$depth_pvalue[i] <- depth_t_test$p.value
}

colnames(t_test_results) = c("Diversity Index", "Region", "Depth")
write.csv(t_test_results, "./Output/Diversity_Index_t_test.csv", row.names = FALSE)


#Twoway-Anova
a <- aov(value ~ Region * Depth, data = subset(Richness_long2, variable == "Chao1"))
summary(a)

# Create an empty data frame to store the ANOVA results
anova_results <- data.frame()

# Loop through each variable and run the ANOVA
for (var in unique(Richness_long2$variable)) {
  a <- aov(value ~ Region * Depth, data = subset(Richness_long2, variable == var))
  anova_table <- summary(a)[[1]]
  Diversity_index <- var
  anova_results <- rbind(anova_results, cbind(anova_table, Diversity_index))
}

# Display the results
anova_results


# Identify significant relationships (p-value < 0.05)
significant_relationships <- anova_results[anova_results$`Pr(>F)` < 0.05, ]
significant_relationships = na.omit(significant_relationships)

library(multcomp)

for (var in unique(significant_relationships$Diversity_index)) {
  a <- aov(value ~ Region * Depth, data = subset(Richness_long2, variable == var ))
  tukey <- TukeyHSD(a)
  print(tukey)
}

parameter = as.vector(sapply(rownames(significant_relationships), function(x) {strsplit(x, "\\s+")[[1]][1]}))

Tukey_output = data.frame()
for (i in 1:nrow(significant_relationships)){
  a = aov(value ~ Region * Depth, data = subset(Richness_long2, variable == significant_relationships$Diversity_index[i]))
  b = TukeyHSD(a)
  if (parameter[i] == "Region"){Tukey_output = rbind(Tukey_output, cbind(b$Region, significant_relationships$Diversity_index[i]))}
  if (parameter[i] == "Depth"){Tukey_output = rbind(Tukey_output, cbind(b$Depth, significant_relationships$Diversity_index[i]))}
  if (parameter[i] == "Region:Depth"){Tukey_output = rbind(Tukey_output, cbind(b$`Region:Depth`, significant_relationships$Diversity_index[i]))} #not sure if works
}
Tukey_output

colnames(Tukey_output) = c("diff", "lwr", "upr", "p adj", "Diversity index")

write.csv(Tukey_output, "./Output/Diversity_Index_Anova_TukeyHSD.csv", row.names = TRUE)

