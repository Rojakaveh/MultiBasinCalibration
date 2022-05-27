
Sys.unsetenv("http_proxy"); Sys.unsetenv("https_proxy")

############################################## #LITTLE OTTER CK AT MONKTON RD.
flowgage_id="04282629" #LITTLE OTTER CK AT MONKTON RD.
flowgage=get_usgs_gage(flowgage_id,begin_date = "2010-01-01",end_date
                       = "2022-01-01")
flowgage$flowdata$Qmm=(flowgage$flowdata$flow)/(flowgage$area*1000)
max(flowgage$flowdata$Qmm)
flowgage$flowdata=subset(flowgage$flowdata, flowgage$flowdata$Qmm<16)




#############################################USGS 04282586 WEST BR DEAD CREEK NEAR ADDISON, VT
flowgage_id="04282586"
flowgage=get_usgs_gage(flowgage_id,begin_date = "2010-01-01",end_date
                       = "2022-01-01")
flowgage$flowdata$Qmm=(flowgage$flowdata$flow)/(flowgage$area*1000)
max(flowgage$flowdata$Qmm)
flowgage$flowdata=subset(flowgage$flowdata, flowgage$flowdata$Qmm<48)
######################################## USGS 04282581 EAST BR DEAD CREEK NEAR BRIDPORT, VT
load("/groups/bse5304g/Roja/EBDCflowdata.RData")#View(flowgage$flowdata$Qm3psRC)
max(flowgage$flowdata$Qmm)
flowgage$flowdata=subset(flowgage$flowdata, flowgage$flowdata$Qmm<20)


###############################
###SBR 2010-2021
############################
x=c(6.25,0.91,9.31,0.33,1.13,0.77,4.10,2.05,0.66,0.16,1.05,2.22,2.69,0.78,0.90,65,84,94)
######################3
##SBR 2015-2021
#####################################
x=c(3.35,0.52,7.58,0.18,1.20,0.74,2.08,1.21,0.58,0.22,0.90,2.16,2.47,0.70,0.90,66,84,94)
######################################
##SMC
##SMC
#####################################
x=c(25.84,0.94,188,0.31,3.49,0.52,4.66,1.50,0.91,0.57,2.93,2.47,0.57,0.68,0.78,65,84,95)
#############################################

pacman::p_load(httr,EcoHydRology,curl,GSODR,data.table,DEoptim,tidyverse)

change_params=""
rm(change_params)
load(paste(path.package("EcoHydRology"), "data/change_params.rda", sep = "/"))
##############################################################
############now for calibration we will add the details of each sensitive parameters
#####for example we know that CN in *.mgt is a sensitieve parameter so we will alter the file type for each hru in 
###a way that hru#1 is in topslope with the range of 65-75
#####hru#2 is in midslope with the range of 75-85
#########hrue#3 nearest to the stream with the range of 85-95
###########################################################
filetype <-c("*1.mgt","*2.mgt","*3.mgt")
parameter<-c("CN2","CN2","CN2")
alter_type <-c("new","new","new")
min <- c(65,75,85)
max<-c(75,85,95)
current <-c(70,80,90)
multi <- c("FALSE","FALSE","FALSE")
frformat <- c("F20,A80","F20,A80","F20,A80")
fwformat <- c("%20.5f,%s","%20.5f,%s","%20.5f,%s")
newdataframe <- data.frame(filetype,parameter,alter_type,min,max,current,multi,frformat,fwformat)
change_params <- rbind(change_params,newdataframe)
##############
library(SWATmodel)
#########################################################
#  Load updated functions
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/SWATmodel/R/readSWAT.R?root=ecohydrology")
save(readSWAT,file="readSWAT.R")
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/EcoHydRology/R/setup_swatcal.R?root=ecohydrology")
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/EcoHydRology/R/swat_objective_function_rch.R?root=ecohydrology")
calib_range=c("1999-12-31","2021-12-31")

######################################
#####here we should enter the parameters that are sensitive (not all of the params that are below) !!!be carful
########################################
params_select=c(1,2,3,6,7,8,9,10,11,14,21,23,24,32,33,42,43,44)
calib_params=change_params[params_select,]
#View(calib_params)
calib_params$min[1]=0
calib_params$min[2]=0
calib_params$max[2]=1
calib_params$current[2]=0.5
calib_params$max[3]=600

calib_params$min[7]=0
calib_params$min[8]=0
calib_params$current[7]=2.5
calib_params$current[8]=2.5
calib_params$min[9]=0.01
calib_params$max[9]=1
calib_params$max[11]=3
calib_params$min[11]=0.5
calib_params$current[12]=2.25
calib_params$max[12]=2.5
calib_params$min[12]=2
calib_params$min[13]=0.5
calib_params$max[13]=3
calib_params$min[15]=0.1



calib_params[1:7]
setup_swatcal(calib_params)
rch=3
swat_objective_function_rch=function (x, calib_range, calib_params, flowgage, rch,save_results=F)
{
  calib_params$current <- x
  tmpdir=as.character(as.integer((runif(1)+1)*10000))
  tmpdir=paste(c(format(Sys.time(), "%s"),tmpdir,Sys.getpid()),sep="",collapse="")
  print(tmpdir)
  dir.create(tmpdir)
  file.copy(list.files(),tmpdir)
  setwd(tmpdir)
  file.remove(list.files(pattern="output."))
  alter_files(calib_params)
  libarch = if (nzchar(base::version$arch)) paste("libs", base::version$arch, sep = "/") else "libs"
  swatbin <- "rswat2012.exe"
  junkout=system(shQuote(paste(path.package("SWATmodel"), libarch, swatbin, sep = "/")),intern = TRUE)
  start_year = read.fortran(textConnection(readLines("file.cio")[9]), "f20")
  load("readSWAT.R")
  outdata = readSWAT("rch",".")
  outdata$FLOW_OUTmm=(outdata$FLOW_OUTcms*24*3600)/(outdata$AREAkm2*1000)
  test2 = subset(outdata, outdata$RCH == rch)
  test3 = merge(flowgage$flowdata, test2, all = F)
  NS = topmodel::NSeff(test3$Qmm, test3$FLOW_OUTmm)
  print(NS)
  if(save_results){file.copy(list.files(),"../")}
  file.remove(list.files())
  setwd("../")
  file.remove(tmpdir)
  return(abs(NS - 1))
}


swat_objective_function_rch(x,calib_range,calib_params,flowgage,rch,save_results=F)
calib_params$current <- x
tmpdir=as.character(as.integer((runif(1)+1)*10000))
tmpdir=paste(c(format(Sys.time(), "%s"),tmpdir,Sys.getpid()),sep="",collapse="")
print(tmpdir)
dir.create(tmpdir)
file.copy(list.files(),tmpdir)
setwd(tmpdir)
file.remove(list.files(pattern="output."))
alter_files(calib_params)
libarch = if (nzchar(base::version$arch)) paste("libs", base::version$arch, sep = "/") else "libs"
swatbin <- "rswat2012.exe"
junkout=system(shQuote(paste(path.package("SWATmodel"), libarch, swatbin, sep = "/")),intern = TRUE)
start_year = read.fortran(textConnection(readLines("file.cio")[9]), "f20")
load("readSWAT.R")
outdata = readSWAT("rch",".")
outdata$FLOW_OUTmm=(outdata$FLOW_OUTcms*24*3600)/(outdata$AREAkm2*1000)
test2 = subset(outdata, outdata$RCH == rch)
test3 = merge(flowgage$flowdata, test2, all = F) 
NSeff(test3$Qmm, test3$FLOW_OUTmm)


Out <- c()
test3$ID=c(1:length(test3$mdate))
id <- unique(test3$ID)
for (i in id){
  #i=1
 
  fit1 <- try(lm(formula=FLOW_OUTmm~Qmm, data = test3[test3$ID != i,]),silent=TRUE)
  Out[i] <- if (inherits(fit1, "lm")) sim = predict(fit1, newdata=test3[test3$ID==i,])
}
NSeff(Out, test3$FLOW_OUTmm)
#########################################
#specify the cross-validation method
ctrl <- trainControl(method = "LOOCV")
model <- train(FLOW_OUTmm~Qmm, data = test3, method = "lm", trControl = ctrl)
print(model)
