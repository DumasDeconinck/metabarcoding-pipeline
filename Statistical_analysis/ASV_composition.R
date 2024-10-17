#############
#Composition
#############

#---libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(knitr)

#---load the data
data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

#Keep only Phylum and ASV number
V9data = V9data[, !grepl("Sequence|Kindom|Domain|Class|Order|Family|Genus|Species", names(V9data))]

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

kable(sort(table(V9data$Phylum),decreasing = TRUE))
kable(sort(table(V9data$Phylum)/length(V9data$Phylum)*100, decreasing = TRUE))

x = data.frame(as.data.frame(table(V9data$Phylum)), as.data.frame(table(V9data$Phylum)/length(V9data$Phylum)*100)[2])
x = x[rev(order(x$Freq)),]
colnames(x) = c("Phylum", "Number of ASVs", "Percentage of ASVs")

write.csv(x, "V9_frequency_ASVs.csv", row.names = FALSE)

#---Long Format, as per usual Plus prepare for ggplot
long_data <- V9data %>%
  pivot_longer(cols = -Phylum, names_to = "Sample", values_to = "Value")

# Reorder the columns
long_data <- long_data %>%
  select(Sample, Phylum, Value)

#Aggregate by sample
aggregated_data <- long_data %>%
  group_by(Sample, Phylum) %>%
  summarize(Value = sum(Value), .groups = "drop")

#Extract metadata from name
aggregated_data$Primer <- sub("\\..*", "", aggregated_data$Sample)
aggregated_data$Station <- sub("^[^.]+\\.", "", aggregated_data$Sample)
aggregated_data$Region <- ifelse(substr(aggregated_data$Station, 1, 1) == "P", "East",
                               ifelse(substr(aggregated_data$Station, 1, 1) == "S", "West", "Unknown"))
aggregated_data$Depth <- ifelse(substr(aggregated_data$Station, nchar(aggregated_data$Station), nchar(aggregated_data$Station)) == "S", "Surface",
                              ifelse(substr(aggregated_data$Station, nchar(aggregated_data$Station), nchar(aggregated_data$Station)) == "B", "Bottom", "Unknown"))

Aggregated_East = subset(aggregated_data, Region == "East")
Aggregated_East_Bottom = subset(Aggregated_East, Depth == "Bottom")
Aggregated_East_Surface = subset(Aggregated_East, Depth == "Surface")

Aggregated_West = subset(aggregated_data, Region == "West")
Aggregated_West_Bottom = subset(Aggregated_West, Depth == "Bottom")
Aggregated_West_Surface = subset(Aggregated_West, Depth == "Surface")

Agregated_data2 = rbind(Aggregated_West_Surface, Aggregated_West_Bottom, Aggregated_East_Surface, Aggregated_East_Bottom)

#Maintain order
Agregated_data2$Station <- factor(Agregated_data2$Station, levels = unique(Agregated_data2$Station))

#---Ggplot bars
color_palette <- brewer.pal(n = length(unique(Agregated_data2$Phylum)), name = "Dark2")

Phylum1 = ggplot(Agregated_data2, aes(fill = Phylum, y= Value, x = Station)) +
  geom_bar(position = "stack", stat ="identity") + 
  theme_classic() +
  ylab("Richness") +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )
Phylum1

Phylum2 = ggplot(Agregated_data2, aes(fill = Phylum, y= Value, x = Station)) +
  geom_bar(position = "fill", stat ="identity") + 
  theme_classic() +
  ylab("Relative abundance") +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )
Phylum2

#---ggplot pie
Pie = Agregated_data2[2:3]

Pie <- aggregate(Value ~ Phylum, Pie, sum) 

Pie <- Pie %>% 
  arrange(desc(Phylum)) %>%
  mutate(prop = Value / sum(Pie$Value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

#Pie <- Pie %>%
  #filter(prop > 0.5) 

Phylum_Pie = ggplot(Pie, aes(fill=Phylum, y=prop, x = "")) + 
  geom_bar(position="stack", stat="identity") +
  coord_polar("y", start=0) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  geom_text(y = Pie$ypos, label = paste(round(Pie$prop),"%", sep = ""), color = "white", size=3) +
  labs(x = NULL,
       y = NULL)
Phylum_Pie

#---ttest of read counts
#Ttest
Agregated_data2
Agregated_data3 <- Agregated_data2 %>%
  group_by(Sample,Region, Depth) %>%
  summarize(Value = sum(Value), .groups = "drop")

t.test(Value ~ Region, data = Agregated_data3)
t.test(Value ~ Depth, data = Agregated_data3)

#Anova
model <- aov(Value ~ Region * Depth, data = Agregated_data3)
summary(model)

###########################################################
#---Repeat for class
#---load the data
data = read.csv("Classified_species_phyto.csv")
V9data = data[, !grepl("Ex|FC|LC|V4", names(data))]
V9data <- V9data[rowSums(V9data[, 3:22] != 0) > 0, ]

#Keep onll class and ASV number
V9data = V9data[, !grepl("Sequence|Kindom|Domain|Phylum|Order|Family|Genus|Species", names(V9data))]

#The easy way to do it now is to just change the column name of class to Phylum and run the previous script.
V9data$Phylum <- V9data$Class
V9data <- V9data[, -which(names(V9data) == "Class")]

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

kable(sort(table(V9data$Phylum),decreasing = TRUE))
kable(sort(table(V9data$Phylum)/length(V9data$Phylum)*100, decreasing = TRUE))

x = data.frame(as.data.frame(table(V9data$Phylum)), as.data.frame(table(V9data$Phylum)/length(V9data$Phylum)*100)[2])
x = x[rev(order(x$Freq)),]
colnames(x) = c("Phylum", "Number of ASVs", "Percentage of ASVs")

write.csv(x, "V9_frequency_ASVs_class.csv", row.names = FALSE)

#---Long Format, as per usual Plus prepare for ggplot
long_data <- V9data %>%
  pivot_longer(cols = -Phylum, names_to = "Sample", values_to = "Value")

# Reorder the columns
long_data <- long_data %>%
  select(Sample, Phylum, Value)

#Aggregate by sample
aggregated_data <- long_data %>%
  group_by(Sample, Phylum) %>%
  summarize(Value = sum(Value), .groups = "drop")

#Extract metadata from name
aggregated_data$Primer <- sub("\\..*", "", aggregated_data$Sample)
aggregated_data$Station <- sub("^[^.]+\\.", "", aggregated_data$Sample)
aggregated_data$Region <- ifelse(substr(aggregated_data$Station, 1, 1) == "P", "East",
                                 ifelse(substr(aggregated_data$Station, 1, 1) == "S", "West", "Unknown"))
aggregated_data$Depth <- ifelse(substr(aggregated_data$Station, nchar(aggregated_data$Station), nchar(aggregated_data$Station)) == "S", "Surface",
                                ifelse(substr(aggregated_data$Station, nchar(aggregated_data$Station), nchar(aggregated_data$Station)) == "B", "Bottom", "Unknown"))

Aggregated_East = subset(aggregated_data, Region == "East")
Aggregated_East_Bottom = subset(Aggregated_East, Depth == "Bottom")
Aggregated_East_Surface = subset(Aggregated_East, Depth == "Surface")

Aggregated_West = subset(aggregated_data, Region == "West")
Aggregated_West_Bottom = subset(Aggregated_West, Depth == "Bottom")
Aggregated_West_Surface = subset(Aggregated_West, Depth == "Surface")

Agregated_data2 = rbind(Aggregated_West_Surface, Aggregated_West_Bottom, Aggregated_East_Surface, Aggregated_East_Bottom)

#Maintain order
Agregated_data2$Station <- factor(Agregated_data2$Station, levels = unique(Agregated_data2$Station))

#---Ggplot bars
# Get the "Dark2" palette as a starting point
dark2_palette <- brewer.pal(n = 8, name = "Dark2")

# Create a custom color palette with 22 colors
color_palette <- colorRampPalette(dark2_palette)(length(unique(Agregated_data2$Phylum))) #Creates a colour pallete based on the one provided

Class1 = ggplot(Agregated_data2, aes(fill = Phylum, y= Value, x = Station)) +
  geom_bar(position = "stack", stat ="identity") + 
  theme_classic() +
  ylab("Richness") +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )
Class1

Class2 = ggplot(Agregated_data2, aes(fill = Phylum, y= Value, x = Station)) +
  geom_bar(position = "fill", stat ="identity") + 
  theme_classic() +
  ylab("Relative abundance") +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )
Class2

#---ggplot pie
Pie = Agregated_data2[2:3]

Pie <- aggregate(Value ~ Phylum, Pie, sum) 

Pie <- Pie %>% 
  arrange(desc(Phylum)) %>%
  mutate(prop = Value / sum(Pie$Value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

#Pie <- Pie %>%
  #filter(prop > 0.5) 

Class_Pie = ggplot(Pie, aes(fill=Phylum, y=prop, x = "")) + 
  geom_bar(position="stack", stat="identity") +
  coord_polar("y", start=0) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  geom_text(y = Pie$ypos, label = paste(round(Pie$prop),"%", sep = ""), color = "white", size=3) +
  labs(x = NULL,
       y = NULL)
Class_Pie

#----Throw it all in one plot
Phylum1
Phylum2
Phylum_Pie

x = ggarrange(
  Phylum1, Phylum2, Phylum_Pie,
  align = "h",
  labels = c("A", "B", "C"),
  common.legend = TRUE,
  legend = "bottom",
  ncol = 3
)

Class1
Class2
Class_Pie

y = ggarrange(
  Class1, Class2, Class_Pie,
  align = "h",
  labels = c("D", "E", "F"),
  common.legend = TRUE,
  legend = "bottom",
  ncol = 3
)

z = ggarrange(x,y,align="v", ncol = 1)
z


