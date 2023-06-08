library("ggmap")
library(sf)
library("rnaturalearthdata")
library("rnaturalearth")

library(legendMap)

library("legendMap")

PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")


world <- ne_countries(scale="medium", returnclass="sf")
class(world)

metadf <- PS.TSS@sam_data

coordf <- data.frame(Lon=metadf$Longitude, Lat=metadf$Latitude)

#sites <- st_as_sf(coordf, coords = c("Lon", "Lat"),
#                             crs = 4326, agr = "constant")

library("ggspatial")

sampling <-
    ggplot(data=world)+
#    ggmap(area)+
    geom_sf()+
    geom_point(data=metadf,aes(x=Longitude, y=Latitude, fill=HI), size=4, alpha=0.5, shape=21)+
    coord_sf(xlim=c(7, 16), ylim=c(47, 55), expand=FALSE)+
    scale_fill_gradient2("Hybrid\nindex", high="red", low="navy", mid="white", midpoint=0.5)+
#    xlab("Longitude", ylab="Latitude")+
#    annotation_scale(location = "bl", width_hint = 0.5) +
#    annotation_north_arrow(location = "bl", which_north = "true",,
#                                         style = north_arrow_fancy_orienteering) +
    scale_bar(lon = 8, lat = 54.2, arrow_length = 10, arrow_distance = 50,
                 distance_lon = 50, distance_lat = 7, distance_legend = 20,
                 dist_unit = "km", orientation = FALSE, legend_size = 3)+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

sampling

saveRDS(sampling, "tmp/map.RDS")

#saveRDS(sampling, "tmp/sampling_map.rds")
ggplot2::ggsave(file="fig/Figure1.pdf", sampling, width = 120, height = 120, dpi = 300, units="mm")
