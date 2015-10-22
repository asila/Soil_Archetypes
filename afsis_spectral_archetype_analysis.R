library(archetypes)
library(lattice)
library(doParallel)
library(ggplot2)
library(dplyr)
registerDoParallel()
getDoParWorkers()
library(readr)
library(reshape2)
	a<-"~/Archetypes/AfSIS"
	#Read mir data
	mir <- read_csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/Calibration_Htsxt_MIR.csv")
	#Read chem data
	chem<-read_csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/AfSIS_reference_data.csv")
	setwd(a)
	set.seed(8970)
	#get number of mir columns
	n<-ncol(mir)
	#fit archetypes
	a <- stepArchetypes (mir[,-c(1,n)],k=1:10,verbose = TRUE, nrep=1)
	#fit robust archetype
	#ra <- robustArchetypes (mir[,-1],k=1:17,verbose = TRUE)
	png(file="Scree plots.png",height=500,width=800)
	screeplot(a)
	dev.off()
	#According to elbow criterion k = 3 or maybe k =6 or 8 are the best number of archetypes
	#Corresponding to Occam'srazor we use 3 archetypes; 
	a3 <- bestModel(a[[3]])
	#Transpose the four archetypes for better readability
	param<-t(parameters(a3))
	#Store the parameters
	write.table(param, file="Parameters.csv",row.names=FALSE)
	#atypes <- apply(coef(a3, "alphas"), 2, which.max)
	#Show simplex plot
	par(mfrow=c(1,1))
	png(file="Simplexplot3.png",height=500,width=1200)
	simplexplot(a3, show_direction = FALSE, show_points =TRUE,radius=400,points_col="grey")
	dev.off()
	#Determine the archetypes
	SSN<-as.vector(mir[,1])
	arch.grps <- as.data.frame(cbind(SSN,paste0("Archetype.",apply(predict(a3,mir[,-c(1,n)]),1,which.max))))
	colnames(arch.grps) <- c("SSN","archetypes" )
	#Use barplot in relation to the original data:
	png(file="Archetype_barplot3.png",height=500,width=1200)
	barplot(a3, mir[,-c(1,n)], percentiles = FALSE)
	dev.off()
	#Or use the original raw spectra to show peaks
	mir.arch<-merge(arch.grps,mir)
	wave<-as.numeric(substr(colnames(mir.arch[,-c(1:2)]),2,19))
	colnames(mir.arch) <- c("SSN","archetypes",wave)
	spec.melt<-melt(mir.arch,id=c("SSN","archetypes"))
	#By spectra
	p<-ggplot(data=spec.melt, aes(x=as.numeric(as.vector(variable)),y=value,group=SSN))+
    geom_line(size=0.34,aes(col=as.numeric(variable)))+scale_colour_gradient(high="red",low="blue")+
    ggtitle("Archetypes for AfSIS reference set raw MIR spectra")+
    xlim(c(4000,600))+
    ylim(c(0,3))+ 
    xlab(expression("Wavenumbers cm"^-1))+
    ylab("Absorbance")+
    #theme with white background
    theme_bw() +
    #eliminates background, gridlines, and chart border
    theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
    )
	p<-p+theme(legend.position = "none")
	p<-p+facet_wrap(~archetypes,ncol=1)
	png(file="Archetypes raw spectra.png")
	p
	dev.off()
	#Aggregate
	mir0<-mir.arch[,-1]
	ag<-aggregate(.~archetypes,data=mir0,mean)
	ag.melt<-melt(ag,id="archetypes")
	p<-ggplot(data=ag.melt,aes(x=as.numeric(as.vector(variable)),y=value,color=archetypes)) + geom_line()+
	ggtitle("Averaged raw MIR spectra by archetype")+
    xlim(c(4000,600))+
    ylim(c(0,2))+ 
    xlab(expression("Wavenumbers cm"^-1))+
    ylab("Absorbance")+
    #theme with white background
    theme_bw() +
    #eliminates background, gridlines, and chart border
    theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
    )
	p<-p+theme(legend.justification=c(0,1),legend.position = c(0,1))
	png(file="Archetypes Averaged raw spectra.png")
	p
	dev.off()
	#Merge arch.grps with chem data
	arch<-merge(chem, arch.grps)
	arch.s<-with(arch,table(paste(Country,Site,sep="."),archetypes))	
	#View somw exploratory plots showing distribution of selected soil properties across the obtained archetypes
	with(arch,bwplot(m3.Al~archetypes))
	with(arch,bwplot(m3.Ca~archetypes))
	with(arch,bwplot(ExAc~archetypes))
	with(arch,bwplot(Na~archetypes))
	with(arch,bwplot(psa.c4sand~archetypes))
	with(arch,bwplot(pH~archetypes))
	with(arch,bwplot(Total.Carbon~archetypes))
	with(arch,bwplot(psa.c4clay~archetypes))
	#which sites dorminate archetype 3
	q<-which(arch.s[,3]>29)
	q
	#Save the spectral archetypes
	write.table(arch.grps, file="AfSIS spec archtypes classes3.csv", sep=",",row.names=FALSE)