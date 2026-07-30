library(here)
suppressMessages({library(tidyverse); library(ggmap); library(ggrepel)})
ggmap::register_stadiamaps(Sys.getenv("STADIA_MAPS_KEY"))
pc <- readRDS(here("output/comparison/cache/processing_cache.rds"))
cols <- c("station","lat","lon","date","temp","salinity","depth_m")
p <- pc$meta_asgard_p2[, cols]
p <- p[p$station != "incubation", ]
p <- p[!is.na(p$lat)&!is.na(p$lon)&!is.na(p$date), ]
.old <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME","C")
p$date_p <- as.POSIXct(p$date, format="%b %d %Y %H:%M", tz="UTC")
p <- p[!is.na(p$date_p), ]; p <- p[order(p$station, p$depth_m), ]
stn <- p[!duplicated(p$station), ]
stn$date_day <- as.Date(stn$date_p); stn$date_num <- as.numeric(stn$date_day)
track <- stn[order(stn$date_p), ]
dbrk <- pretty(stn$date_day, n=5); dlab <- format(dbrk, "%B %d")
Sys.setlocale("LC_TIME", .old)
cat("N proc stations:", nrow(stn), "\n")
bbox <- make_bbox(lon=stn$lon, lat=stn$lat, f=0.15)
mapz <- get_stadiamap(bbox, maptype="stamen_terrain", zoom=6)
fc <- stn[complete.cases(stn[,c("lon","lat","temp","salinity")]), ]
xr<-range(fc$lon)+c(-0.5,0.5); yr<-range(fc$lat)+c(-0.3,0.3)
kx<-111*cos(mean(fc$lat)*pi/180); ky<-111
g<-expand.grid(lon=seq(xr[1],xr[2],length.out=240), lat=seq(yr[1],yr[2],length.out=240))
D<-sqrt((outer(g$lon,fc$lon,"-")*kx)^2+(outer(g$lat,fc$lat,"-")*ky)^2)
W<-1/(D^2); W[!is.finite(W)]<-max(W[is.finite(W)],na.rm=TRUE)*1e6; nn<-apply(D,1,min)
jet<-colorRampPalette(c("#00007F","#0000FF","#007FFF","#00FFFF","#7FFF7F","#FFFF00","#FF7F00","#FF0000","#7F0000"))
th<-theme(axis.title=element_text(size=16,face="bold"), axis.text=element_text(size=12),
          legend.title=element_text(size=13,face="bold"), legend.text=element_text(size=11))
dscale<-scale_color_viridis_c(name="Date", option="C", breaks=as.numeric(dbrk), labels=dlab, guide=guide_colorbar(order=1))
fieldmap<-function(vals, brks, name, pal, lab_every=1){
  gg<-g; gg$z<-as.vector((W%*%vals)/rowSums(W)); gg$z[nn>30]<-NA
  ggmap(mapz)+
    geom_contour_filled(data=gg[!is.na(gg$z),], aes(lon,lat,z=z,fill=after_stat(level_mid)), breaks=brks, alpha=0.55)+
    scale_fill_stepsn(colours=pal(length(brks)-1), breaks=brks, limits=range(brks), name=name,
                      labels=function(b) ifelse(b %% lab_every == 0, as.character(b), ""),
                      guide=guide_colorsteps(barheight=grid::unit(6,"cm"), show.limits=TRUE, order=2))+
    geom_path(data=track, aes(lon,lat), color="grey30", linewidth=0.4, alpha=0.6)+
    geom_point(data=stn, aes(lon,lat,color=date_num), size=2.8, alpha=0.95)+ dscale+
    geom_text_repel(data=stn, aes(lon,lat,label=station), size=2.4, max.overlaps=Inf, alpha=0.85)+
    labs(x="Longitude", y="Latitude")+th
}
pdf(here("output_p/maps/ASGARD_processing_field_date_zoom6.pdf"), width=12, height=9)
print(fieldmap(fc$salinity, seq(25,33,by=0.5), "Salinity (psu)", viridisLite::viridis, 1))
print(fieldmap(fc$temp, seq(-2,11,by=1), "Temperature (°C)", jet, 2))
dev.off(); cat("DONE\n")
