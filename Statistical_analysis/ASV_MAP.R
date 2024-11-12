########
#map
######

#---Library
library(dplyr)
library(tidyr)
library(readxl)

library(ggmap)
library(RColorBrewer)
library(ggspatial)
library(ggrepel)
library(osmdata)
library(sf)
library(dplyr)
library(ggnewscale)
library(rnaturalearth)
library(rnaturalearthhires)
library(MASS)
library(gridExtra)


#---read data
HK_locations <- read_excel("Location_data.xlsx")

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

# Long format
long_data <- V9data %>%
  pivot_longer(cols = -ASV.number, names_to = "Sample", values_to = "Value")

# Reorder the columns
long_data <- dplyr::select(long_data, Sample, ASV.number, Value)

#---extract information from name
long_data$Region <- ifelse(substr(long_data$Sample, 1, 1) == "P", "East",
                           ifelse(substr(long_data$Sample, 1, 1) == "S", "West", "Unknown"))
long_data$Depth <- ifelse(substr(long_data$Sample, nchar(long_data$Sample), nchar(long_data$Sample)) == "S", "Surface",
                          ifelse(substr(long_data$Sample, nchar(long_data$Sample), nchar(long_data$Sample)) == "B", "Bottom", "Unknown"))
long_data$Location <- substr(long_data$Sample, 1, nchar(long_data$Sample) - 1)

#Match them with location
long_data <- long_data %>%
  left_join(HK_locations, by = c("Location" = "Sample"))

#---Map
ggmap::register_stadiamaps("YOUR STADIA CODE") #You have to create a stadia account and create 

x_min <- 113.78
x_max <- 114.5101
y_min <- 22.13
y_max <- 22.5667

HK_map <- get_stadiamap(bbox = c(x_min, y_min, x_max, y_max), zoom = 12, maptype = "stamen_terrain_background")

#Setup the extras
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

#+ scale_bar + north_arrow + scale_text

#ASV7
x1 = subset(long_data, ASV.number == "ASV_7")
x2 = subset(x1, Depth == "Surface")
x3 = subset(x1, Depth == "Bottom")

map_plot7 <- ggmap(HK_map) + 
  geom_point(data = x2, aes(x = Lon, y = Lat, size = Value), color = brewer.pal(8, "Dark2")[5]) +
  geom_point(data = x3, aes(x = Lon, y = Lat, size = Value), shape = 1) +
  coord_map() +  theme_minimal() +
  labs(size = "Relative\nAbundance") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")  + scale_bar + north_arrow + 
  annotate(
    "text",
    x = scale_text$x,
    y = scale_text$y,
    label = scale_text$label,
    size = 3,
    hjust = 1, # Set hjust to 1 for right alignment
    vjust = 0
  )

#ASV23
x1 = subset(long_data, ASV.number == "ASV_23")
x2 = subset(x1, Depth == "Surface")
x3 = subset(x1, Depth == "Bottom")

map_plot8 <- ggmap(HK_map) + 
  geom_point(data = x2, aes(x = Lon, y = Lat, size = Value), color = brewer.pal(8, "Dark2")[6]) +
  geom_point(data = x3, aes(x = Lon, y = Lat, size = Value), shape = 1) +
  coord_map() +  theme_minimal() +
  labs(size = "Relative\nAbundance") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")  + scale_bar + north_arrow + 
  annotate(
    "text",
    x = scale_text$x,
    y = scale_text$y,
    label = scale_text$label,
    size = 3,
    hjust = 1, # Set hjust to 1 for right alignment
    vjust = 0
  )

#ASV288 An eastern ASV
x1 = subset(long_data, ASV.number == "ASV_288")
x2 = subset(x1, Depth == "Surface")
x3 = subset(x1, Depth == "Bottom")

map_plot9 <- ggmap(HK_map) + 
  geom_point(data = x2, aes(x = Lon, y = Lat, size = Value), color = brewer.pal(8, "Dark2")[7]) +
  geom_point(data = x3, aes(x = Lon, y = Lat, size = Value), shape = 1) +
  coord_map() +  theme_minimal() +
  labs(size = "Relative\nAbundance") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")  + scale_bar + north_arrow + 
  annotate(
    "text",
    x = scale_text$x,
    y = scale_text$y,
    label = scale_text$label,
    size = 3,
    hjust = 1, # Set hjust to 1 for right alignment
    vjust = 0
  )

#ASV195 the most common western only ASV
x1 = subset(long_data, ASV.number == "ASV_195")
x2 = subset(x1, Depth == "Surface")
x3 = subset(x1, Depth == "Bottom")

map_plot10 <- ggmap(HK_map) + 
  geom_point(data = x2, aes(x = Lon, y = Lat, size = Value), color = brewer.pal(8, "Dark2")[8]) +
  geom_point(data = x3, aes(x = Lon, y = Lat, size = Value), shape = 1) +
  coord_map() +  theme_minimal() +
  labs(size = "Relative\nAbundance") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")  + scale_bar + north_arrow + 
  annotate(
    "text",
    x = scale_text$x,
    y = scale_text$y,
    label = scale_text$label,
    size = 3,
    hjust = 1, # Set hjust to 1 for right alignment
    vjust = 0
  )

Maps = grid.arrange(map_plot7 + ggtitle("A"), 
             map_plot8 + ggtitle("B"),
             map_plot9 + ggtitle("C"), 
             map_plot10 + ggtitle("D"))


#---Sample map
line_segments_TH <- data.frame(
  x = subset(HK_locations, Location == "Tolo Harbour")$Lon[-nrow(subset(HK_locations, Location == "Tolo Harbour"))],
  y = subset(HK_locations, Location == "Tolo Harbour")$Lat[-nrow(subset(HK_locations, Location == "Tolo Harbour"))],
  xend = subset(HK_locations, Location == "Tolo Harbour")$Lon[-1],
  yend = subset(HK_locations, Location == "Tolo Harbour")$Lat[-1]
)

line_segments_W <- data.frame(
  x = subset(HK_locations, Location == "Western")$Lon[-nrow(subset(HK_locations, Location == "Western"))],
  y = subset(HK_locations, Location == "Western")$Lat[-nrow(subset(HK_locations, Location == "Western"))],
  xend = subset(HK_locations, Location == "Western")$Lon[-1],
  yend = subset(HK_locations, Location == "Western")$Lat[-1]
)


map_plot <- ggmap(HK_map) +
  geom_segment(data = line_segments_TH, aes(x = x, y = y, xend = xend, yend = yend), color = brewer.pal(8, "Dark2")[1], size = 1, arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(data = line_segments_W, aes(x = x, y = y, xend = xend, yend = yend), color = brewer.pal(8, "Dark2")[2], size = 1, arrow = arrow(length = unit(0.3, "cm"))) +
  geom_point(data = HK_locations, aes(x = Lon, y = Lat), color = brewer.pal(8, "Dark2")[3], size = 2) +
  geom_text_repel(
    aes(x = Lon, y = Lat, label = Sample), 
    data = HK_locations, 
    hjust = -0.4, 
    vjust = -0.4, 
    size = 6) +
  coord_map() +
  theme_minimal() + 
  theme(axis.title = element_blank(),
        legend.title = element_blank())  + scale_bar + north_arrow + 
  annotate(
    "text",
    x = scale_text$x,
    y = scale_text$y,
    label = scale_text$label,
    size = 3,
    hjust = 1, # Set hjust to 1 for right alignment
    vjust = 0
  )

