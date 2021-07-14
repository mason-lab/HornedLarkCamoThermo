library(gplots)
library(plotrix)
library(MuMIn)
citation("MuMIn")

library(agricolae)#Provides HSD.test()
citation("agricolae")

setwd("~/HornedLarkCamoThermo")

####################
### READ IN DATA ###
####################
col_plot_df<-read.table("Utility/col_plot_df.txt")

### Read in voucher data ###
vouchers<-read.csv("Data/HOLA_VoucherTable_v4.csv")
rownames(vouchers)<-paste("MVZ",vouchers$catalognumber,sep="")
vouchers$infraspecific <-factor(vouchers$infraspecific,levels=c("strigata","merrilli","lamprochroma","utahensis","rubea","sierrae","actia","insularis","ammophila","leucansiptila","occidentalis","adusta","alpestris","arcticola","enertera","hoyti","aphrasta"))

length(table(vouchers$infraspecific)) #17 subspecies

### Read in pavo plumage data ###
back_light_df<-read.csv("Data/HOLA_PavoColorData_v1.csv",stringsAsFactors=F,row.names=1)

### Read in plumage pattern data ###
pattern<-read.csv("Data/HOLA_BackPattern_v1.csv",stringsAsFactors=F,row.names=1)
for(i in 1:ncol(pattern)){
	pattern[,i]<-as.numeric(pattern[,i])
}

### Read in FullSpectrumData Set ###
hola_ir_dorsal <-read.csv("Data/HOLA_Dorsal_FullSpectrum.csv", stringsAsFactors=F,row.names=1)

# ### Perform PCA on solar reflectance ###
# solar_dorsal_pca<-princomp(hola_ir_dorsal[,1:3])
# solar_dorsal_pca$loadings
# plot(solar_dorsal_pca$scores[,1],solar_dorsal_pca$scores[,2], bg=paste0(col2hex(col_plot_df[as.character(vouchers[rownames(solar_dorsal_pca$scores),]$infraspecific),1]),"95"),pch=21)

### Read in abiotic data ###
bio_data<-read.csv("Data/HOLA_BioClim_v1.csv", stringsAsFactors=F,row.names=1)
for(i in 1:ncol(bio_data)){
	bio_data[,i]<-as.numeric(bio_data[,i])
}

### Perform PCA on bioclim ###
bio_pca<-princomp(bio_data)
bio_pca$scores[,1]<-bio_pca$scores[,1]*-1
bio_pca$loadings

### write out supplementary table with bioclim PCA loadings ###
#write.csv(bio_pca$loadings[,1:3],file="HOLA_BioClimPCALoadings.csv")

### Read in soil composition data ###
soil_df<-read.csv("Data/HOLA_SoilComposition_v1.csv", stringsAsFactors=F,row.names=1)
soil_pca<-princomp(soil_df)
soil_pca$loadings

### Run granularity soil pca ###
grain_pca<-princomp(soil_df[,c(1,2,3,4)])

### write out supplementary table with soil composition PCA loadings ###
#write.csv(soil_pca$loadings[,1:3],file="HOLA_SoilPCALoadings.csv")

### Read in soil color data ###
soil_col_df<-read.csv("Data/HOLA_SoilColor_v1.csv",stringsAsFactors=F,row.names=1)
soil_col_pca<-princomp(na.omit(soil_col_df[,1:3]))
nrow(soil_col_pca$scores)

### write out supplementary table with soil color PCA loadings ###
#write.csv(soil_col_pca$loadings[,1:3],file="HOLA_SoilColorPCALoadings.csv")

### Create RGB Values for plotting ###
soil_color_rgb<-apply(na.omit(soil_col_df),1,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))
length(soil_color_rgb)

nrow(vouchers[!rownames(vouchers) %in% names(soil_color_rgb),])

#####################
### DATA ANALYSIS ###
#####################
### Set up plotting factor and colors from darkest to lightest ###
factor_df<-data.frame(ssp=vouchers$infraspecific,hola_ir_dorsal)
avg_ref<-tapply(factor_df$dorsal_solar_reflectance, factor_df$ssp,mean)

vouchers$infraspecific <-factor(vouchers$infraspecific,levels=names(sort(avg_ref)))
ssp_fact<-vouchers$infraspecific
names(ssp_fact)<-rownames(vouchers)

### Figure 2 ###
### Phenotypic varition results / plot ###
col_data<-data.frame(row.names=rownames(vouchers),ssp=vouchers$infraspecificepithet,sex=factor(vouchers$sex),lumMean=pattern[rownames(vouchers),]$lumMean,r.achieved=back_light_df[rownames(vouchers),]$r.achieved,h.theta= back_light_df[rownames(vouchers),]$h.theta, h.phi=back_light_df[rownames(vouchers),]$h.phi,power= pattern[rownames(vouchers),]$sumPower)
col_data $ssp <-factor(col_data $ssp,levels=c("strigata","merrilli","lamprochroma","utahensis","rubea","sierrae","actia","insularis","ammophila","leucansiptila","occidentalis","adusta","alpestris","arcticola","enertera","hoyti","aphrasta"))

### Ensure numerical characters are stored as numerical (not factor) ###
for(i in 3:7){
	col_data[,i]<-as.numeric(col_data[,i])
}

color_pca<-princomp(col_data[,c(3,4,7)]) #Principal component analysis of luminosity, chroma, and power
color_pca$loadings

### write out supplementary table with color PCA loadings ###
#write.csv(color_pca$loadings[,1:3],file="HOLA_SoilColor_v1.csv")

col_data<-cbind(col_data,color_pca$scores[,1:2])
col_data_sex<-split(col_data,col_data$sex)

### Tukey's post-hoc test split by sex to group ssp variation ###
## Female mean luminosity ##
glm_female_lumMean<-glm(lumMean~ssp,data=col_data_sex$female)
aov_female_lumMean<-aov(glm_female_lumMean)
hsdtest_female_lumMean<-HSD.test(aov_female_lumMean,trt="ssp",group=T,console=T)
hsdtest_female_lumMean$groups[levels(col_data_sex$female$ssp)[levels(col_data_sex$female$ssp) %in% rownames(hsdtest_female_lumMean$groups)],]

## Male mean luminosity ##
glm_male_lumMean<-glm(lumMean~ssp,data=col_data_sex$male)
aov_male_lumMean<-aov(glm_male_lumMean)
hsdtest_male_lumMean<-HSD.test(aov_male_lumMean,trt="ssp",group=T,console=T)
hsdtest_male_lumMean$groups[levels(col_data_sex$male$ssp)[levels(col_data_sex$male$ssp) %in% rownames(hsdtest_male_lumMean$groups)],]

## Female r.achieved ##
glm_female_r.achieved<-glm(r.achieved~ssp,data=col_data_sex$female)
aov_female_r.achieved<-aov(glm_female_r.achieved)
hsdtest_female_r.achieved<-HSD.test(aov_female_r.achieved,trt="ssp",group=T,console=T)
hsdtest_female_r.achieved$groups[levels(col_data_sex$female$ssp)[levels(col_data_sex$female$ssp) %in% rownames(hsdtest_female_r.achieved$groups)],]

## Male r.achieved ##
glm_male_r.achieved<-glm(r.achieved~ssp,data=col_data_sex$male)
aov_male_r.achieved<-aov(glm_male_r.achieved)
hsdtest_male_r.achieved<-HSD.test(aov_male_r.achieved,trt="ssp",group=T,console=T)
hsdtest_male_r.achieved$groups[levels(col_data_sex$male$ssp)[levels(col_data_sex$male$ssp) %in% rownames(hsdtest_male_r.achieved$groups)],]

## Female power ##
glm_female_power<-glm(power~ssp,data=col_data_sex$female)
aov_female_power<-aov(glm_female_power)
hsdtest_female_power<-HSD.test(aov_female_power,trt="ssp",group=T,console=T)
hsdtest_female_power$groups[levels(col_data_sex$female$ssp)[levels(col_data_sex$female$ssp) %in% rownames(hsdtest_female_power$groups)],]

## Male power ##
glm_male_power<-glm(power~ssp,data=col_data_sex$male)
aov_male_power<-aov(glm_male_power)
hsdtest_male_power<-HSD.test(aov_male_power,trt="ssp",group=T,console=T)
hsdtest_male_power$groups[levels(col_data_sex$male$ssp)[levels(col_data_sex$male$ssp) %in% rownames(hsdtest_male_power$groups)],]

### Phenotypic variation summary plot ###
pdf(file="Figures/PhenotypicSummary_v1.pdf",width=6.5,height=4.5)
# quartz(width=6.5,height=4.5)

### Set up plot layout ###
layout(matrix(c(1,1,2,2,2,1,1,3,3,3,1,1,4,4,4,1,1,5,5,5,1,1,6,6,6,1,1,7,7,7,8,8,8,8,8,8,8,8,8,8),byrow=T,nrow=8))
#layout.show(8)

#par(mfrow=c(1,2))

par(mar=c(9,4,1,1))

### Scatterplot of PCA space for females ###
plot(col_data_sex[[1]]$Comp.1, col_data_sex[[1]]$Comp.2,pch=21,bg= paste0(col2hex(col_plot_df[as.character(col_data_sex[[1]]$ssp),1]),"95"),xlab="",ylab="",axes=F,cex=1,lwd=0.5,xlim=range(col_data$Comp.1),ylim=range(col_data$Comp.2))

### Add points of PCA space for males ###
points(col_data_sex[[2]]$Comp.1, col_data_sex[[2]]$Comp.2,pch=22,bg= paste0(col2hex(col_plot_df[as.character(col_data_sex[[2]]$ssp),1]),"95"),cex=1,lwd=0.5)

### Axis labels ###
box()

axis(1,mgp=c(0,0.35,0),tck=-0.025,cex.axis=0.7)
axis(2,mgp=c(0,0.35,0),tck=-0.025,cex.axis=0.7)

mtext(paste("PC1 (",round(color_pca $sdev[1]/sum(color_pca $sdev),4) * 100,"%)",sep=""),side=1,line=1.5,cex=0.75)
mtext(paste("PC2 (",round(color_pca $sdev[2]/sum(color_pca $sdev),4) * 100,"%)",sep=""),side=2,line=1.5,cex=0.75)

text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.05,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.05,label=LETTERS[1],font=2,cex=1.25)

### Annotate PCA axes ###
par(xpd=NA)
#locator(1)
text(x=11500,y=-5000,label="Brighter\nMore Patterning",cex=0.7)
text(x=-6500,y=-5000,label="Darker\nLess Patterning",cex=0.7)

#locator(1)
text(x=-14000,y=5250,label="Brighter\nLess Patterning",cex=0.7,srt=90)
text(x=-14000,y=-2000,label="Darker\nMore Patterning",cex=0.7,srt=90)

### Box plots of male and female PC1, PC2 among subspecies ###
### Three characters, males and females, six different plots ###
par(mar=c(0.5,3,0.5,0.5))
par(xpd=T)

## Female lumMean ##
### Make Bar Plot ###
bp_foo<-boxplot(lumMean~ssp,data=col_data_sex$female,axes=F,lwd=0.5,col=col_plot_df[levels(col_data_sex$female$ssp),1],ylim=c(range(col_data_sex$female$lumMean)[1],diff(range(col_data_sex$female$lumMean))*0.4 + range(col_data_sex$female$lumMean)[2]))

### Add Tukey's labels ###
labs<-toupper(hsdtest_female_lumMean $groups[levels(col_data_sex$female$ssp),2])
labs[is.na(labs)]<-""
text(label=labs,y= diff(range(col_data_sex$female$lumMean))*0.3 + range(col_data_sex$female$lumMean)[2]
,x=1:length(levels(col_data_sex$female$ssp)),cex=0.5)

### Add corner label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.02,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.2,label=LETTERS[2],font=2,cex=1.25)

box()

## Male lumMean ##
### Make Bar Plot ###
boxplot(lumMean~ssp,data=col_data_sex$male,axes=F,lwd=0.5,col=col_plot_df[levels(col_data_sex$male$ssp),1],ylim=c(range(col_data_sex$male$lumMean)[1],diff(range(col_data_sex$male$lumMean))*0.4 + range(col_data_sex$male$lumMean)[2]))

### Add Tukey's labels ###
labs<-toupper(hsdtest_male_lumMean $groups[levels(col_data_sex$male$ssp),2])
labs[is.na(labs)]<-""
text(label=labs,y= diff(range(col_data_sex$male$lumMean))*0.3 + range(col_data_sex$male$lumMean)[2]
,x=1:length(levels(col_data_sex$male$ssp)),cex=0.5)

### Add corner label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.02,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.2,label=LETTERS[3],font=2,cex=1.25)

box()

## Female r.achieved ##
### Make Bar Plot ###
bp_foo<-boxplot(r.achieved~ssp,data=col_data_sex$female,axes=F,lwd=0.5,col=col_plot_df[levels(col_data_sex$female$ssp),1],ylim=c(range(col_data_sex$female$r.achieved)[1],diff(range(col_data_sex$female$r.achieved))*0.4 + range(col_data_sex$female$r.achieved)[2]))

### Add Tukey's labels ###
labs<-toupper(hsdtest_female_r.achieved$groups[levels(col_data_sex$female$ssp),2])
labs[is.na(labs)]<-""
text(label=labs,y= diff(range(col_data_sex$female$r.achieved))*0.3 + range(col_data_sex$female$r.achieved)[2]
,x=1:length(levels(col_data_sex$female$ssp)),cex=0.5)

### Add corner label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.02,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.2,label=LETTERS[4],font=2,cex=1.25)

box()

## Male r.achieved ##
### Make Bar Plot ###
boxplot(r.achieved~ssp,data=col_data_sex$male,axes=F,lwd=0.5,col=col_plot_df[levels(col_data_sex$male$ssp),1],ylim=c(range(col_data_sex$male$r.achieved)[1],diff(range(col_data_sex$male$r.achieved))*0.4 + range(col_data_sex$male$r.achieved)[2]))

### Add Tukey's labels ###
labs<-toupper(hsdtest_male_r.achieved $groups[levels(col_data_sex$male$ssp),2])
labs[is.na(labs)]<-""
text(label=labs,y= diff(range(col_data_sex$male$r.achieved))*0.3 + range(col_data_sex$male$r.achieved)[2]
,x=1:length(levels(col_data_sex$male$ssp)),cex=0.5)

### Add corner label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.02,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.2,label=LETTERS[5],font=2,cex=1.25)

box()

## Female power ##
### Make Bar Plot ###
boxplot(power~ssp,data=col_data_sex$female,axes=F,lwd=0.5,col=col_plot_df[levels(col_data_sex$female$ssp),1],ylim=c(range(col_data_sex$female$power)[1],diff(range(col_data_sex$female$power))*0.4 + range(col_data_sex$female$power)[2]))

### Add Tukey's labels ###
labs<-toupper(hsdtest_female_power$groups[levels(col_data_sex$female$ssp),2])
labs[is.na(labs)]<-""

text(label=labs,y= diff(range(col_data_sex$female$power))*0.3 + range(col_data_sex$female$power)[2]
,x=1:length(levels(col_data_sex$female$ssp)),cex=0.5)

### Add corner label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.02,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.2,label=LETTERS[6],font=2,cex=1.25)

box()

## Male power ##
### Make Bar Plot ###
boxplot(power~ssp,data=col_data_sex$male,axes=F,lwd=0.5,col=col_plot_df[levels(col_data_sex$male$ssp),1],ylim=c(range(col_data_sex$male$power)[1],diff(range(col_data_sex$male$power))*0.4 + range(col_data_sex$male$power)[2]))

### Add Tukey's labels ###
labs<-toupper(hsdtest_male_power $groups[levels(col_data_sex$male$ssp),2])
labs[is.na(labs)]<-""

text(label=labs,y= diff(range(col_data_sex$male$power))*0.3 + range(col_data_sex$male$power)[2]
,x=1:length(levels(col_data_sex$male$ssp)),cex=0.5)

### Add corner label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.02,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.2,label=LETTERS[7],font=2,cex=1.25)

box()

### Annotate plot ###
par(xpd=NA)

for(i in 1:length(levels(col_data_sex$male$ssp))){
	text(label= levels(col_data_sex$male$ssp)[i],x=i,y=2500,srt=90,adj=c(1,0.5),col= col_plot_df[levels(col_data_sex$male$ssp)[i],1],font=3,cex=1.35)
}

### Add additional annotations ###
#mfpoints<-locator(2)
#mfpoints<-mfpoints$y
mfpoints<-c(17391.65, 157321.68)
text(y=seq(mfpoints[1],mfpoints[2],length.out=6),x=-0.5,labels=rep(c("male","female"),4),srt=90)

#charpoints<-locator(2)
#charpoints<-charpoints$y
charpoints<-c(30910.94,142979.33)

text(y=seq(charpoints[1], charpoints[2],length.out=3),x=-1.25,labels=c("Back patterning power","Acheived chroma","Brightness"),srt=90)

dev.off()

### GENERAL LINEAR MODELS ###
### Set up data for GLMs / association tests / Appendix ###
color_dataframe<-data.frame(ssp=vouchers[rownames(soil_col_pca$scores),]$infraspecific ,sex=factor(vouchers[rownames(soil_col_pca$scores),]$sex),soilpc1=soil_col_pca$scores[,1], soilpc2=soil_col_pca$scores[,2],lum= pattern[rownames(soil_col_pca$scores),]$lumMean,r.acheived= as.numeric(back_light_df[rownames(soil_col_pca$scores),]$r.achieved))

soil_dataframe<-data.frame(ssp=vouchers[rownames(grain_pca$scores),]$infraspecific ,sex=vouchers[rownames(grain_pca$scores),]$sex, soil_grain=grain_pca$scores[,1],power= pattern[rownames(grain_pca$scores),]$sumPower)

### GLM tests for lum, PC1, ssp, and sex ###
lum_soilpc1_sex1<-glm(lum~soilpc1+sex,data= color_dataframe)

### Summary statistics for GLM reported in manuscript ###
summary(lum_soilpc1_sex1)

### GLM tests for chroma, PC2, ssp, and sex ###
r_soilpc2_sex1<-glm(r.acheived~soilpc2+sex,data= color_dataframe)

### Summary statistics for GLM reported in manuscript ###
summary(r_soilpc2_sex1)

### GLM tests for granluarity, power, ssp, and sex ###
pow_soilgran_sex1<-glm(power~soil_grain+sex,data= soil_dataframe)

### Summary statistics for GLM reported in manuscript ###
summary(pow_soilgran_sex1)

######################################################
### Three-panel plot figure of background matching ###
######################################################

png(width=6.5,height=2.25,file="Figures/BackgroundMatchingPlot_v7.png",units="in",res=500)
#quartz(width=6.5,height=2.35)

par(mar=c(2.75,2.75,1.5,1))
par(mfrow=c(1,3))

## Panel A Soil and plumage brightness ##
plot(soil_col_pca$scores[,1], pattern[rownames(soil_col_pca$scores),]$lumMean,pch=21,bg= soil_color_rgb[rownames(soil_col_pca$scores)],col=col_plot_df[as.character(ssp_fact[rownames(soil_col_pca$scores)]),1],lwd=0.45,cex.axis=0.85,axes=F,xlab="",ylab="",ylim=c(range(pattern[rownames(soil_col_pca$scores),]$lumMean)[1],diff(range(pattern[rownames(soil_col_pca$scores),]$lumMean))*0.1 + range(pattern[rownames(soil_col_pca$scores),]$lumMean)[2]))

box()
axis(1,tck=-0.025,mgp=c(0,0.25,0),cex.axis= 0.625)
mtext("Soil Brightness",side=1,line=1.5,cex=0.75)

axis(2,tck=-0.025,mgp=c(0,0.5,0),cex.axis=0.625)
mtext("Plumage Brightness",side=2,line=1.5,cex=0.75)

test_a<-cor.test(soil_col_pca$scores[,1], pattern[rownames(soil_col_pca$scores),]$lumMean)

### Stats annotation ###
text(x=par("usr")[2]-diff(c(par("usr")[1],par("usr")[2]))*0.025,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.05,label= paste0("r = ",round(test_a$estimate,2),"; P < 0.001"),font=1,cex=1,adj=c(1,0.5))

### Letter label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.05,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.05,label=LETTERS[1],font=2,cex=1.25)

abline(lm(pattern[rownames(soil_col_pca$scores),]$lumMean~ soil_col_pca$scores[,1]))

## Panel B Soil redness and plumage chroma ##
plot(soil_col_pca $scores[,2], back_light_df[rownames(soil_col_pca$scores),]$r.achieved,pch=21,bg= soil_color_rgb[rownames(soil_col_pca$scores)],col=col_plot_df[as.character(ssp_fact[rownames(soil_col_pca$scores)]),1],lwd=0.45,axes=F,xlab="",ylab="",ylim=c(range(back_light_df[rownames(soil_col_pca$scores),]$r.achieved)[1],diff(range(back_light_df[rownames(soil_col_pca$scores),]$r.achieved))*0.1 + range(back_light_df[rownames(soil_col_pca$scores),]$r.achieved)[2]))

box()
axis(1,tck=-0.025,mgp=c(0,0.25,0),cex.axis= 0.625)
mtext("Soil Redness",side=1,line=1.5,cex=0.75)

axis(2,tck=-0.025,mgp=c(0,0.5,0),cex.axis=0.625)
mtext("Acheived Chroma",side=2,line=1.5,cex=0.75)

test_b<-cor.test(soil_col_pca $scores[,2], as.numeric(back_light_df[rownames(soil_col_pca$scores),]$r.achieved))

### Stats annotation ###
text(x=par("usr")[2]-diff(c(par("usr")[1],par("usr")[2]))*0.025,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.05,label= paste0("r = ",round(test_b$estimate,2),"; P < 0.001"),font=1,cex=1,adj=c(1,0.5))

### Letter label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.05,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.05,label=LETTERS[2],font=2,cex=1.25)

abline(lm(back_light_df[rownames(soil_col_pca$scores),]$r.achieved~ soil_col_pca$scores[,2]))

## Panel C Soil granularity and back patterning ##
plot(grain_pca$scores[,1],pattern$sumPower,pch=21,lwd=0.45,col=col_plot_df[as.character(ssp_fact[rownames(grain_pca$scores)]),1],axes=F,ylim=c(range(pattern$sumPower)[1],diff(range(pattern$sumPower))*0.1 + range(pattern$sumPower)[2]))#bg=soil_color_rgb[rownames(grain_pca$scores)],

test_c<-cor.test(grain_pca$scores[,1],pattern$sumPower)

### Stats annotation ###
text(x=par("usr")[2]-diff(c(par("usr")[1],par("usr")[2]))*0.025,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.05,label= paste0("r = ",round(test_c$estimate,2),"; P < 0.001"),font=1,cex=1,adj=c(1,0.5))

### Letter label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.05,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.05,label=LETTERS[3],font=2,cex=1.25)

abline(lm(pattern$sumPower~grain_pca$scores[,1]))

box()
axis(1,tck=-0.025,mgp=c(0,0.25,0),cex.axis= 0.625)
mtext("Soil Granularity",side=1,line=1.5,cex=0.75)

axis(2,tck=-0.025,mgp=c(0,0.5,0),cex.axis=0.625)
mtext("Back Patterning",side=2,line=1.5,cex=0.75)

dev.off()

### GLMS / ANOVAS ###
### Dorsal reflectance ~ precipitation + soil color ###
dorsal_glm_df<-na.omit(data.frame(lumMean=pattern[rownames(soil_col_pca$scores),]$lumMean,dorsal_solar_reflectance=hola_ir_dorsal[rownames(soil_col_pca$scores),]$dorsal_solar_reflectance,soilbright=soil_col_pca$scores[,1],season= bio_pca$scores[rownames(soil_col_pca$scores),1],precip=bio_pca$scores[rownames(soil_col_pca$scores),2],temp=bio_pca$scores[rownames(soil_col_pca$scores),3],sex=vouchers[rownames(soil_col_pca$scores),]$sex))

lumMean_fm<-glm(lumMean~soilbright+season+precip+temp+sex,data=dorsal_glm_df,na.action = "na.fail")
sink("Output/dorsal_lumMean_glm.txt")
summary(lumMean_fm)
sink()
summary(lumMean_fm)$coefficients

dorsal_fm<-glm(dorsal_solar_reflectance~soilbright+season+precip+temp+sex,data=dorsal_glm_df,na.action = "na.fail")
sink("Output/dorsal_solar_absorptance_glm.txt")
summary(dorsal_fm)
sink()
summary(dorsal_fm)$coefficients

