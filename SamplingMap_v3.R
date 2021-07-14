library(raster)
library(sp)
library(xlsx)
library(spatialEco)
require(ggsci)
require(pavo)
require(gplots)
require(rgdal)
require(rgeos)
require(png)

setwd("~/HornedLarkCamoThermo")

### readin plate ###
larkPNG<-readPNG("Figures/alpestris.png")

col_plot_df<-read.table("Utility/col_plot_df.txt")
col_plot_df

### Read in GADM 0 of USA ###
usa1<-readOGR("~/GIS/GADM/USA_adm/USA_adm1.shp")
usa1_trim<-usa1[!usa1@data$NAME_1 %in% c("Alaska","Hawaii"),]
usa1_trim_simp<-gSimplify(usa1_trim,0.001)

plot(usa1_trim_simp)

vouchers<-read.csv("Data/HOLA_VoucherTable_v4.csv")
rownames(vouchers)<-paste("MVZ",vouchers$catalognumber,sep="")
vouchers$infraspecific <-factor(vouchers$infraspecific,levels=c("strigata","merrilli","lamprochroma","utahensis","rubea","sierrae","actia","insularis","ammophila","leucansiptila","occidentalis","adusta","alpestris","arcticola","enertera","hoyti","aphrasta"))

### Remove articola and alpestris ssp as they area island / outside of US ###
vouchers_nomig<-vouchers[!vouchers$infraspecific %in% c("arcticola","alpestris","hoyti","insularis","enertera"),]
vouchers_rgb<-vouchers[rownames(vouchers) %in% rownames(sc_df_pca$scores),]
vouchers_nomig$state

soil_col_df<-read.csv("Data/HOLA_SoilColor_v1.csv",stringsAsFactors=F,row.names=1)
sum_sc_scores<-apply(soil_col_df,1,sum)
sc_df_pca<-princomp(na.omit(soil_col_df))
soil_color_rgb<-apply(na.omit(soil_col_df),1,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))

soil_stack<-stack("Data/CONUS_005cm.tif")
#soil_proj<-projectRaster(soil_stack,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

### Crop soil stack ###
plotRGB(soil_stack,r=1,g=2,b=3)
soil_stack_crop<-crop(soil_stack,extent(x=c(-2398016.2,-548631),y=c(1003325,2899387)))

### Create sampling map ####
### Read in shape files for different ssp ###
hola_shp<-list.files("Data/ShapeFiles/",pattern=".shp",full.names=T)
hola_shp_short<-list.files("Data/ShapeFiles/",pattern=".shp",full.names=F)

shp_list<-list()
for(i in 1:length(hola_shp)){
	shp_list[[i]]<-readOGR(hola_shp[i])
}
names(shp_list)<-gsub(".shp","",sapply(strsplit(hola_shp_short,"_"),function(x) x[3]))
shp_list_trim<-shp_list[names(shp_list) %in% vouchers_nomig$infraspecific]

### Transform projection to match usgs soil data ###
shp_list_trans<-list()
#plotRGB(soil_stack_crop,r=1,g=2,b=3)
for(i in 1:length(shp_list_trim)){
	print(i)
	shp_list_trans[[i]]<-gIntersection(gBuffer(shp_list_trim[[i]],byid=T,width=0),usa1_trim_simp)
	foo<-try(shp_list_trans[[i]])
	if(class(foo) == "try-error"){next}else{
		names(shp_list_trans)[i]<-names(shp_list_trim)[i]
		shp_list_trans[[i]]<-spTransform(shp_list_trans[[i]],proj4string(soil_stack_crop))
		#plot(shp_list_trans[[i]],col=paste0(col2hex(col_plot_df[names(shp_list_trim)[i],1]),"60"),add=T)
	}
}

### Create sampling map ###
png(file="Figures/Fig1_samplingmap_v5.png",height=3.25,width=3.25,units="in",res=300)

quartz(width=3.25,height=3.25)
plotRGB(soil_stack_crop,r=1,g=2,b=3)

for(i in 1:length(shp_list_trans)){
	plot(shp_list_trans[[i]],col= paste0(col2hex(col_plot_df[names(shp_list_trim)[i],1]),"60"),add=T)
}
plot(spTransform(usa1_trim_simp, proj4string(soil_stack_crop)),add=T)


for(i in 1:length(levels(vouchers$infraspecific))){
	ssp_points<-SpatialPoints(matrix(as.numeric(cbind(vouchers[vouchers$infraspecific== levels(vouchers$infraspecific)[i],]$decimallongitude,vouchers[vouchers$infraspecific==levels(vouchers$infraspecific)[i],]$decimallatitude)),ncol=2))
	proj4string(ssp_points)<-"+proj=longlat +datum=WGS84"
	ssp_points<-spTransform(ssp_points,proj4string(soil_stack_crop))
	plot(ssp_points,add=T,bg= col2hex(col_plot_df[levels(vouchers$infraspecific)[i],1]),pch=21,cex=1,lwd=0.75)
}

rect(-1251114.8,1778182,-651867.3, 2830770,col="#FFFFFF80")

### Legend ###
vouch_tab<-table(vouchers$infraspecific)
vouch_tab[levels(vouchers$infraspecific)]
points(x=rep(-1201232,length(vouch_tab)),y=seq(2783762, 1830182,length.out=(length(vouch_tab)+1))[-1],pch=21,bg=col_plot_df[names(vouch_tab),1],cex=1,lwd=0.75)

text(label= c("Subspecies",paste0(names(vouch_tab)," (",vouch_tab,")")),x=rep(-1152043,length(vouch_tab)+1),y=seq(2783762, 1830182,length.out=(length(vouch_tab)+1)),adj=c(0,0.5),cex=0.50,font=c(1,rep(3,length(levels(vouchers$infraspecific)))))


### Add lark plate ###
rasterImage(larkPNG,xleft=-2372212,ybottom= 1066232,xright=-1854174,ytop= 1428670)

### Add scale bar ###
rect(xleft= -1754174,xright= -1754174-100000,ybottom= 1066232,ytop= 1086232,col="white")
rect(xleft= -1754174-100000,xright= -1754174-200000,ybottom= 1066232,ytop= 1086232,col="black")
text(c("0","100","200","km"),x=rev(c(-1754174+100000,-1754174,-1754174-100000,-1754174-200000)),y=rep(1036232,4),cex=0.5)
dev.off()
