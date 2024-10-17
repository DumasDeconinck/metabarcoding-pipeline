#####
#Environmental data
####

library(readxl)
library(tidyverse)
library(car)
library(data.table)
library(ggplot2)
library(corrplot)
library(emmeans)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)

HK_env  = read_excel("Environmental data.xlsx")
HK_env$`EC (µS/cm)` = NULL

parameter_names = colnames(HK_env)[(length(colnames(HK_env))-5): length(colnames(HK_env))]

colnames(HK_env) = c("Station", "Depth", "Depth2", "Temperature", "DO", "Salinity", "pH", "Turbidity", "Chl")

# Create a new column 'Area' based on the first letter of the Station
HK_env$Area = ifelse(substring(HK_env$Station, 1, 1) == "P", "East", "West")

##Graphical representation
HK_env3 = HK_env
HK_env3$Station = NULL
HK_env3$Depth = NULL

HK_env_long = melt(setDT(HK_env3), id.vars = c("Depth2", "Area"), variable.name = "Parameter")

HK_env_long

ggplot(HK_env_long, aes(x = Parameter, y = value, fill = interaction(Area, Depth2))) +
  geom_boxplot(position = position_dodge(0.9)) +
  scale_y_sqrt() +
  scale_x_discrete(labels = parameter_names) +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_manual(values = c(brewer.pal(4, "Dark2")[2], brewer.pal(4, "Dark2")[3], brewer.pal(4, "Dark2")[1]), 
    labels = c("West Bottom", "East Surface", "West Surface")) + #Manually changed
  theme_classic()+
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.box.background = element_rect(color="black", size=0),
        legend.box.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Type III anova for depth (For depth, only look at western samples)
# Perform the MANOVA test
manova_result <- manova(cbind(Temperature, DO, Salinity, pH, Turbidity, Chl) ~ Depth2, data = subset(HK_env, Area == "West"))

# Summarize the MANOVA results
summary(manova_result)

# Perform the MANOVA test for region (only look at surface samples)
manova_result <- manova(cbind(Temperature, DO, Salinity, pH, Turbidity, Chl) ~ Area, data = subset(HK_env, Depth2 == "Surface"))
summary(manova_result)

# Summarize the MANOVA results
summary(manova_result)

# Perform individual ANOVAs for each response factor -  Depth
anova_results <- data.frame(Factor = character(), p_value = numeric())
for (factor in c("Temperature", "DO", "Salinity", "pH", "Turbidity", "Chl")) {
  aov_area <- Anova(aov(formula = as.formula(paste(factor, "~ Depth2")), data = subset(HK_env, Area == "West")), type = "III")
  p_value <- aov_area[["Pr(>F)"]][2]
  anova_results <- rbind(anova_results, data.frame(Factor = factor, p_value = p_value))
}

# Display the results
knitr::kable(anova_results, digits = 3)
significant_factors = anova_results[anova_results['p_value'] < 0.05]
significant_factors = significant_factors[1:(length(significant_factors)/2)]

#FOLLOWUP FOR Depth2
library(emmeans)
for (factor in significant_factors) {
  emm <- emmeans(aov(formula = as.formula(paste(factor, "~ Depth2")), data = subset(HK_env, Area == "West")), "Depth2")
  contrast <- contrast(emm, method = "pairwise")
  print(factor)
  print(contrast)
}

# Perform individual ANOVAs for each response factor -  Area
anova_results <- data.frame(Factor = character(), p_value = numeric())
for (factor in c("Temperature", "DO", "Salinity", "pH", "Turbidity", "Chl")) {
  aov_area <- Anova(aov(formula = as.formula(paste(factor, "~ Area")), data = subset(HK_env, Depth2 == "Surface")), type = "III")
  p_value <- aov_area[["Pr(>F)"]][2]
  anova_results <- rbind(anova_results, data.frame(Factor = factor, p_value = p_value))
}

# Display the results
knitr::kable(anova_results, digits = 3)
significant_factors = anova_results[anova_results['p_value'] < 0.05]
significant_factors = significant_factors[1:(length(significant_factors)/2)]

library(emmeans)
for (factor in significant_factors) {
  emm <- emmeans(aov(formula = as.formula(paste(factor, "~ Area")), data = subset(HK_env, Depth2 == "Surface")), "Area")
  contrast <- contrast(emm, method = "pairwise")
  print(factor)
  print(contrast)
}

######Correlation
cor_matrix <- cor(HK_env[, c("Depth", "Temperature", "Salinity", "Chl")])
corrplot(cor_matrix, method = "circle")

###############################
#MAPS
########################

#---Libraries
library(ggmap)
library(readxl)
library(ggspatial)
library(ggrepel)
library(RColorBrewer)
library(osmdata)
library(sf)
library(dplyr)
library(ggnewscale)
library(rnaturalearth)
library(rnaturalearthhires)
library(MASS)
library(gridExtra)



#StadiaMAP API key
ggmap::register_stadiamaps("9b87c6a6-dbdf-4f7e-8da6-92f2b28c72b0") #You have to create a stadia account and create 

HK_locations <- read_excel("Location_data.xlsx")
HK_env = read_excel("Environmental data.xlsx")
colnames(HK_env) = c("Station", "Depth", "Depth2", "Temperature", "DO", "EC", "Salinity", "pH", "Turbidity", "Chl")
HK_env_surf = subset(HK_env, Depth2 == "Surface") #set to Surface or Bottom
HK_env_bot = subset(HK_env, Depth2 == "Bottom")
HK_locations2 = HK_locations[5:10,]


x_min <- 113.78
x_max <- 114.5101
y_min <- 22.13
y_max <- 22.5667

HK_map <- get_stadiamap(bbox = c(x_min, y_min, x_max, y_max), zoom = 12, maptype = "stamen_terrain_background")

map_plot <- ggmap(HK_map) + 
  geom_point(data = HK_locations, aes(x = Lon, y = Lat, size = HK_env_surf$Temperature), color = brewer.pal(8, "Dark2")[3]) +
  geom_point(data = HK_locations2, aes(x = Lon, y = Lat, size = HK_env_bot$Temperature), shape = 1) +
  coord_map() +
  theme_minimal() + 
  theme(axis.title = element_blank(),
        legend.position = "bottom") +
  labs(size = "Temperature (°C)")

scale_bar <- ggsn::scalebar(
  location = "bottomright",
  dist = 10,
  st.size = 2,
  height = 0.015,
  x.min = x_min,
  x.max = x_max,
  y.min = y_min,
  y.max = y_max,
  transform = TRUE,
  model = "WGS84",
  dist_unit = "km"
) # Add scale bar

north_arrow <- annotation_north_arrow(
  location = "topleft",
  which_north = "true",
  pad_x = unit(0.5, "cm"),
  pad_y = unit(0.5, "cm"),
  width = unit(0.9, "cm") # Adjust the width of the north arrow here
) # Add north arrow

#Add text annotation on the right side of the plot
scale_text <- data.frame(
  x = x_max - 0.035,
  y = y_min + 0.01,
  label = "10 KM"
)

final_plot <- map_plot + scale_bar + north_arrow +
  annotate(
    "text",
    x = scale_text$x,
    y = scale_text$y,
    label = scale_text$label,
    size = 3,
    hjust = 1, # Set hjust to 1 for right alignment
    vjust = 0
  )

print(final_plot)

map_plot2 <- ggmap(HK_map) + 
  geom_point(data = HK_locations, aes(x = Lon, y = Lat, size = HK_env_surf$DO), color = brewer.pal(8, "Dark2")[1]) +
  geom_point(data = HK_locations2, aes(x = Lon, y = Lat, size = HK_env_bot$DO), shape = 1) +
  coord_map() +
  theme_minimal() + 
  labs(size = "DO saturation (%)") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")

final_plot2 <- map_plot2 + scale_bar + north_arrow +
  annotate(
    "text",
    x = scale_text$x,
    y = scale_text$y,
    label = scale_text$label,
    size = 3,
    hjust = 1, # Set hjust to 1 for right alignment
    vjust = 0
  )

print(final_plot2)

map_plot3 <- ggmap(HK_map) + 
  geom_point(data = HK_locations, aes(x = Lon, y = Lat, size = HK_env_surf$Salinity), color = brewer.pal(8, "Dark2")[2]) +
  geom_point(data = HK_locations2, aes(x = Lon, y = Lat, size = HK_env_bot$Salinity), shape = 1) +
  coord_map() +
  theme_minimal() + 
  guides(size = guide_legend(nrow = 1)) +
  labs(size = "Salinity (PSU)") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")

final_plot3 <- map_plot3 + scale_bar + north_arrow +
  annotate(
    "text",
    x = scale_text$x,
    y = scale_text$y,
    label = scale_text$label,
    size = 3,
    hjust = 1, # Set hjust to 1 for right alignment
    vjust = 0
  )

print(final_plot3)

map_plot4 <- ggmap(HK_map) + 
  geom_point(data = HK_locations, aes(x = Lon, y = Lat, size = HK_env_surf$pH), color = brewer.pal(8, "Dark2")[4]) +
  geom_point(data = HK_locations2, aes(x = Lon, y = Lat, size = HK_env_bot$pH), shape = 1) +
  coord_map() +
  theme_minimal() + 
  labs(size = "pH") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")

final_plot4 <- map_plot4 + scale_bar + north_arrow +
  annotate(
    "text",
    x = scale_text$x,
    y = scale_text$y,
    label = scale_text$label,
    size = 3,
    hjust = 1, # Set hjust to 1 for right alignment
    vjust = 0
  )

print(final_plot4)

map_plot5 <- ggmap(HK_map) + 
  geom_point(data = HK_locations, aes(x = Lon, y = Lat, size = HK_env_surf$Turbidity), color = brewer.pal(8, "Dark2")[5]) +
  geom_point(data = HK_locations2, aes(x = Lon, y = Lat, size = HK_env_bot$Turbidity), shape = 1) +
  coord_map() +
  theme_minimal() + 
  labs(size = "Turbidity (FNU)") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")

final_plot5 <- map_plot5 + scale_bar + north_arrow +
  annotate(
    "text",
    x = scale_text$x,
    y = scale_text$y,
    label = scale_text$label,
    size = 3,
    hjust = 1, # Set hjust to 1 for right alignment
    vjust = 0
  )

print(final_plot5)

map_plot6 <- ggmap(HK_map) + 
  geom_point(data = HK_locations, aes(x = Lon, y = Lat, size = HK_env_surf$Chl), color = brewer.pal(8, "Dark2")[6]) +
  geom_point(data = HK_locations2, aes(x = Lon, y = Lat, size = HK_env_bot$Chl), shape = 1) +
  coord_map() +
  theme_minimal() + 
  labs(size = "Chlorphyl a (µg/L)") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")

final_plot6 <- map_plot6 + scale_bar + north_arrow +
  annotate(
    "text",
    x = scale_text$x,
    y = scale_text$y,
    label = scale_text$label,
    size = 3,
    hjust = 1, # Set hjust to 1 for right alignment
    vjust = 0
  )

print(final_plot6)

grid.arrange(
  final_plot + ggtitle("A"),  final_plot2 + ggtitle("B"), final_plot3 + ggtitle("C"), final_plot4 + ggtitle("D"), final_plot5 + ggtitle("E"), final_plot6 + ggtitle("F"),
  map_plot7 + ggtitle("G"),   map_plot8 + ggtitle("H"),  map_plot9 + ggtitle("I"),  map_plot10 + ggtitle("J"), 
  ncol = 2, nrow = 5
)
