##------SENSITIVITY ANALYSIS-----------------
#----------STEP 1: GETTING 10DATASETS FROM OUTDEOPTIM FOR 2015-2021 AND 2010-2021 FOR LOC_FERISSBURG-----------
#LOC_Ferrisburg_2010, is the  name of folder for LOC_ferrisburg 2010-2021 SWAT initialization
#LOC_Ferrisburg_2015, is the  name of folder for LOC_ferrisburg 2015-2021 SWAT initialization
library(SWATmodel)
flowgage_id="04282650" #Little Otter Creek at Ferrisburg, VT.
flowgage=get_usgs_gage(flowgage_id,begin_date = "2010-01-01",end_date= "2022-01-01")
flowgage$flowdata$Qm3ps= flowgage$flowdata$flow/24/3600 #m3/s
######remove extreme datapoint
max(flowgage$flowdata$Qm3ps)
flowgage$flowdata=subset(flowgage$flowdata, flowgage$flowdata$Qm3ps<30)

##moving initialization folder to dev/shm directory for faster run
dir.create("/dev/shm/rojakaveh")
setwd("/dev/shm/rojakaveh")
dir.create("CrackFlowFer2010")
dir.create("CrackFlowFer2015")
file.copy(list.files("~/multibasinpaper/LOC_Ferrisburg_2010/",full.names = TRUE,recursive = TRUE),"/dev/shm/rojakaveh/CrackFlowFer2010/",recursive = TRUE)
file.copy(list.files("~/multibasinpaper/LOC_Ferrisburg_2015/",full.names = TRUE,recursive = TRUE),"/dev/shm/rojakaveh/CrackFlowFer2015/",recursive = TRUE)
####load updated functions
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/SWATmodel/R/readSWAT.R?root=ecohydrology")
save(readSWAT,file="readSWAT.R")
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/EcoHydRology/R/setup_swatcal.R?root=ecohydrology")
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/EcoHydRology/R/swat_objective_function_rch.R?root=ecohydrology")

########load calib params
change_params=""
rm(change_params)
load(paste(path.package("EcoHydRology"), "data/change_params.rda", sep = "/"))
calib_range=c("1999-12-31","2021-12-31")
params_select=c(1,2,3,4,5,6,7,8,9,10,11,14,19,21,23,24,32,33)
calib_params=change_params[params_select,]
#calib_params$current[2]=0.5
calib_params$min[9]=0
calib_params$min[10]=0
calib_params$current[9]=2.5
calib_params$current[10]=2.5
calib_params$min[11]=0.01
calib_params$max[11]=1
calib_params$min[13]=40
calib_params$max[13]=95
calib_params$min[14]=0.3
calib_params$max[14]=3
#calib_params$min[15]=0.3
#calib_params$max[15]=3
calib_params$min[15]=0.3
calib_params$max[15]=3
calib_params$min[16]=0.3
calib_params$max[16]=3
calib_params$max[2]=2
calib_params$max[3]=600
calib_params$max[4]=0.3
calib_params$min[2]=0
calib_params$max[2]=1
calib_params$current[2]=0.5
calib_params[c(9,10),4]=0
calib_params[c(9,10),6]=2.5
calib_params[11,5]=1
calib_params[11,6]=.5
calib_params$min[1]=0


calib_params[1:7]
setup_swatcal(calib_params)
rch=3

# Test calibration
x=calib_params$current
swat_objective_function_rch(x,calib_range,calib_params,flowgage,rch,save_results=F)

#################################
#Loop to have 10datasets of outdeoptim
for (calitts in 1:5){
  print(calitts)
  set.seed(calitts)
  outDEoptim<-DEoptim(swat_objective_function_rch,calib_params$min,calib_params$max,
                      DEoptim.control(strategy = 6,NP = 16,itermax=300,parallelType = 1,
                                      packages = c("SWATmodel")),calib_range,calib_params,flowgage,rch)
  x=outDEoptim$optim$bestmem
  filename=format(Sys.time(), "%Y%m%d%H%M_outDEoptim.RData")
  assign(filename,outDEoptim)
  save(outDEoptim,file=format(Sys.time(), "~/multibasinpaper/%Y%m%d%H%M_outDEoptim.RData"))
}

################################
#####lOAD DATASETS
######################################
rm(list = c("NSE"))
for(i in list.files(pattern = "20211130*")){
  load(i)
  if (!exists("NSE")){ 
    NSE=data.frame(1-outDEoptim$member$bestvalit)
  } else { 
    tmpdf1=data.frame(1-outDEoptim$member$bestvalit)
    NSE=cbind(NSE,tmpdf1)
  }
}

rm(list = c("NSE"))
for(i in list.files(pattern = "20211130*")){
  load(i)
  if (!exists("NSE")){ 
    NSE=data.frame(1-outDEoptim$member$bestvalit)
  } else { 
    tmpdf1=data.frame(1-outDEoptim$member$bestvalit)
    NSE=cbind(NSE,tmpdf1)
  }
}
#######################################
#######chnage the colnames
########################################
varnames=c("Dataset")
for (i in 1:length(NSE)){
  names(NSE)[1:length(NSE)] = paste0(rep(varnames),"_",1:length(NSE))
}
####need to load calib params
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/SWATmodel/R/readSWAT.R?root=ecohydrology")
save(readSWAT,file="readSWAT.R")
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/EcoHydRology/R/setup_swatcal.R?root=ecohydrology")
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/EcoHydRology/R/swat_objective_function_rch.R?root=ecohydrology")


change_params=""
rm(change_params)
load(paste(path.package("EcoHydRology"), "data/change_params.rda", sep = "/"))
calib_range=c("1999-12-31","2021-12-31")
params_select=c(1,2,3,4,5,6,7,8,9,10,11,14,19,21,23,24,32,33)
calib_params=change_params[params_select,]
#calib_params$current[2]=0.5
calib_params$min[9]=0
calib_params$min[10]=0
calib_params$current[9]=2.5
calib_params$current[10]=2.5
calib_params$min[11]=0.01
calib_params$max[11]=1
calib_params$min[13]=40
calib_params$max[13]=95
calib_params$min[14]=0.3
calib_params$max[14]=3
#calib_params$min[15]=0.3
#calib_params$max[15]=3
calib_params$min[15]=0.3
calib_params$max[15]=3
calib_params$min[16]=0.3
calib_params$max[16]=3
calib_params$max[2]=2
calib_params$max[3]=600
calib_params$max[4]=0.3
calib_params$min[2]=0
calib_params$max[2]=1
calib_params$current[2]=0.5
calib_params[c(9,10),4]=0
calib_params[c(9,10),6]=2.5
calib_params[11,5]=1
calib_params[11,6]=.5
calib_params$min[1]=0

calib_params[1:7]
#################################################
######plot  params
###############################################
#setwd("~")
rm(list = c("BestMemit"))
for(i in list.files(pattern = "20211130*")){
  load(i)
  if (!exists("BestMemit")){ 
    BestMemit=data.frame(outDEoptim$member$bestmemit)
  } else { 
    tmpdf2=data.frame(outDEoptim$member$bestmemit)
    BestMemit=cbind(BestMemit,tmpdf2)
  }
}
par(mar= c(4,4,1,1), oma = rep(1,4))
layout(rbind(matrix(1:18, 3), rep(19, 6)), heights = c(rep(3, 3),1))
for (i in 1:18){
  plot(BestMemit[,i], ylab= "",xlab="",ylim=c(calib_params$min[i],calib_params$max[i]),cex=0.3)
  axis(1,font = 2)
  axis(2,font = 2)
  mtext(side=1, line=2, "iteration", col="black",cex=1)
  mtext(side=2, line=3, calib_params$parameter[i], col="black", font=2, cex=1)
  points(BestMemit[,i+18], col=2,cex=0.3)
  points(BestMemit[,i+36], col=3,cex=0.3)
  points(BestMemit[,i+54], col=4,cex=0.3)
  points(BestMemit[,i+72], col=5,cex=0.3)
  points(BestMemit[,i+90], col=6,cex=0.3)
  points(BestMemit[,i+108], col=7,cex=0.3)
  points(BestMemit[,i+126], col=8,cex=0.3)
  points(BestMemit[,i+144], col=9,cex=0.3)
  points(BestMemit[,i+162], col=10,cex=0.3)
  
}
plot.new()
par(mar=c(6,0.25,0.25,0.25))
colornames=c(seq(1:length(NSE)))
legendname= paste0(rep(varnames),"_",1:length(NSE))
legend("bottom", legend =c(as.expression(bquote(bold("Dataset_1"))),as.expression(bquote(bold("Dataset_2"))),as.expression(bquote(bold("Dataset_3"))),as.expression(bquote(bold("Dataset_4"))),as.expression(bquote(bold("Dataset_5"))),as.expression(bquote(bold("Dataset_6"))),as.expression(bquote(bold("Dataset_7"))),as.expression(bquote(bold("Dataset_8"))),as.expression(bquote(bold("Dataset_9"))),as.expression(bquote(bold("Dataset_10")))), col = colornames,lty = 1:2, cex = 1,ncol = 5,bty="n",  inset=c(0,1), xpd=TRUE)
dev.off()
#############################################
#####percent changes
#######################################

rm(list = c("myx","tmpdf3","percentchange"))
for(i in list.files(pattern = "20211130*")){
  load(i)
  if (!exists("myx")){ 
    myx=data.frame(outDEoptim$member$bestmemit)
    myx$NSE=1-outDEoptim$member$bestvalit
    percentchange=((tail(myx,-1)-myx[1:(length(myx[,1])-1),]))/(myx[1:(length(myx[,1])-1),])*100
    #myx=tail(myx,-1)
  } else { 
    tmpdf3=data.frame(outDEoptim$member$bestmemit)
    tmpdf3$NSE=1-outDEoptim$member$bestvalit      
    deltatmp1=((tail(tmpdf3,-1)-tmpdf3[1:(length(tmpdf3[,1])-1),]))/(tmpdf3[1:(length(tmpdf3[,1])-1),])*100
    #tmpdf3=tail(tmpdf3,-1)
    myx=rbind(myx,tmpdf3)
    percentchange=rbind(percentchange,deltatmp1)
  }
}

deltanames=c("%Changes in")
par(mar= c(4,4,1,1), oma = rep(1,4))
layout(matrix(1:18, 3), heights = c(rep(3, 3)))
#par(mar=c(4,4,1,1))
#par(mfrow=c(7,3))

for (i in 1:18){
  plot(percentchange[,i], percentchange$NSE,ylim=c(0,100),ylab= "",xlab="",cex=0.5)
  #plot(percentchange[,i], percentchange$NSE,ylab= "%changes NSE",xlab=paste0(rep(deltanames)," ",calib_params$parameter[i]),cex=0.5)
  axis(1,font = 2)
  axis(2,font = 2)
  mtext(side=1, line=2, paste0(rep(deltanames)," ",calib_params$parameter[i]),font=2, col="black",cex=0.8)
  mtext(side=2, line=3, "%changes NSE", col="black", font=2, cex=1)
}
dev.off()


#setwd("~")
rm(list = c("myx","tmpdf","deltaparams"))
for(i in list.files(pattern = "20211130*")){
  load(i)
  if (!exists("myx")){ 
    myx=data.frame(outDEoptim$member$bestmemit)
    myx$NSE=1-outDEoptim$member$bestvalit
    deltaparams=(tail(myx,-1)-myx[1:(length(myx[,1])-1),])
    myx=tail(myx,-1)
  } else { 
    tmpdf=data.frame(outDEoptim$member$bestmemit)
    tmpdf$NSE=1-outDEoptim$member$bestvalit      
    deltatmp=(tail(tmpdf,-1)-tmpdf[1:(length(tmpdf[,1])-1),])
    tmpdf=tail(tmpdf,-1)
    myx=rbind(myx,tmpdf)
    deltaparams=rbind(deltaparams,deltatmp)
  }
}

initialvalues=myx[1:(length(myx[,1])-1),]

for (i in colnames(deltaparams)){
  #i="CNavg"
  nam = paste0("junkpar_",i)
  assign(nam,deltaparams[!deltaparams[,fmatch(i,names(deltaparams))]==0,] )
  nam2=paste0("initialval_",i)
  assign(nam2,initialvalues[!deltaparams[,fmatch(i,names(deltaparams))]==0,])
  nam3=paste0("Sr_",i)
  assign(nam3,abs(((get(paste0("junkpar_",i))$NSE)/
                     (get(paste0("junkpar_",i))[fmatch(i,names(deltaparams))]))*
                    (get(paste0("initialval_",i))[fmatch(i,names(deltaparams))]/
                       get(paste0("initialval_",i))$NSE)))
}


#####percent changes
#
rm(list = c("myx","tmpdf3","percentchange"))
for(i in list.files(pattern = "20211130*")){
  load(i)
  if (!exists("myx")){ 
    myx=data.frame(outDEoptim$member$bestmemit)
    myx$NSE=1-outDEoptim$member$bestvalit
    percentchange=((tail(myx,-1)-myx[1:(length(myx[,1])-1),]))/(myx[1:(length(myx[,1])-1),])*100
    #myx=tail(myx,-1)
  } else { 
    tmpdf3=data.frame(outDEoptim$member$bestmemit)
    tmpdf3$NSE=1-outDEoptim$member$bestvalit      
    deltatmp1=((tail(tmpdf3,-1)-tmpdf3[1:(length(tmpdf3[,1])-1),]))/(tmpdf3[1:(length(tmpdf3[,1])-1),])*100
    #tmpdf3=tail(tmpdf3,-1)
    myx=rbind(myx,tmpdf3)
    percentchange=rbind(percentchange,deltatmp1)
  }
}
#percentchange=percentchange[,!names(percentchange)%in%c("par12","par13")]

deltanames=c("%Changes in")

####################################
#########################################################
######## Sensitive Params
##################################
rm(list = c("tmpdf","outDEoptim","deltatmp","i"))
# Finding Columns/Params that can  be removed
junk= deltaparams[(rowSums(deltaparams[,1:18]>0)>0)&deltaparams$NSE==0,]
rmparams=colnames(junk[colSums(junk)!=0])
calib_params$parameter[colSums(junk)==0]
X=myx[deltaparams$NSE!=0,!(names(myx) %in% rmparams)]
X=X[1:(length(X)-1)]
sensparamcols=which(colSums(junk)==0)[1:(length(which(colSums(junk)==0))-1)]
sensitivity_params=calib_params[sensparamcols,]

rm("myx")
for(i in list.files(pattern = "20211130*")){
  load(i)
  if (!exists("myx")){
    myx=data.frame((outDEoptim$member$bestmemit))
    #    myx$NSE=1-outDEoptim$optim$bestval
  } else {
    tmpdf=data.frame((outDEoptim$member$bestmemit))
    #    tmpdf$NSE=1-outDEoptim$optim$bestval
    myx=rbind(myx,tmpdf)
  }
}
myx=myx[c(1,2,3,6,7,8,9,10,11,12,13,14,15,16,17,18)]
calib_params$parameter[c(1,2,3,6,7,8,9,10,11,12,13,14,15,16,17,18)]

colnames(myx)=calib_params$parameter[c(1,2,3,6,7,8,9,10,11,12,13,14,15,16,17,18)]

summary(myx)
RelSensi=base::as.data.frame(stack(list(ALPHA_BF=Sr_par2,SMFMX=Sr_par9,CN2=Sr_par13,Ksat=Sr_par16,TIMP=Sr_par11,
                                        Depth=Sr_par14,GW_DELAY=Sr_par1,AW=Sr_par15,ESCO=Sr_par17,SURLAG=Sr_par12,
                                        RCHRG_DP=Sr_par6,GWQMN=Sr_par3,EPCO=Sr_par18,SMFMN=Sr_par10,SMTMP=Sr_par8,
                                        SFTMP=Sr_par7
)))


# load the library
library(forcats)
library(ggplot2)
library(dplyr)
#

###########################################################################
#### violin + boxplot + mean 
#######################
names(RelSensi)[2]="Parameter"
library(ggplot2)
# Basic violin plot
p <- ggplot(RelSensi, aes(x=Parameter, y=values, color=Parameter)) + 
  coord_cartesian(ylim = c(0,7))+
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme(legend.position = "None")

p

p + stat_summary(fun.y=mean, geom="point", shape=8, size=4,col="black") +
  xlab("Parameter")+
  ylab("Relative Sensitivity")+
  theme(axis.text = element_text(size =15))+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.title = element_text(size =17))+
  theme(axis.title = element_text(face="bold"))
geom_boxplot(width=0.1)

dev.off()
