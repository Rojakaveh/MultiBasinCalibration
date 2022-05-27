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
View(calib_params)
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

# Test calibration
x=calib_params$current
##-------------getting the observed data to make observed dataframe--------------
######################################## USGS 04282581 EAST BR DEAD CREEK NEAR BRIDPORT, VT
load("/groups/bse5304g/Roja/EBDCflowdata.RData")#View(flowgage$flowdata$Qm3psRC)
max(flowgage$flowdata$Qmm)
flowgage$flowdata=subset(flowgage$flowdata, flowgage$flowdata$Qmm<48)

flowgage_id3="04282581"
flowgage3=get_usgs_gage(flowgage_id3,begin_date = "2010-01-01",end_date
                        = "2022-01-01")
############################################## #LITTLE OTTER CK AT MONKTON RD.
flowgage_id1="04282629" #LITTLE OTTER CK AT MONKTON RD.
flowgage1=get_usgs_gage(flowgage_id1,begin_date = "2010-01-01",end_date
                        = "2022-01-01")
flowgage1$flowdata$Qmm=(flowgage1$flowdata$flow)/(flowgage1$area*1000)
max(flowgage1$flowdata$Qmm)
flowgage1$flowdata=subset(flowgage1$flowdata, flowgage1$flowdata$Qmm<16)

#############################################USGS 04282586 WEST BR DEAD CREEK NEAR ADDISON, VT
flowgage_id2="04282586"
flowgage2=get_usgs_gage(flowgage_id2,begin_date = "2010-01-01",end_date
                        = "2022-01-01")
flowgage2$flowdata$Qmm=(flowgage2$flowdata$flow)/(flowgage2$area*1000)
max(flowgage2$flowdata$Qmm)
flowgage2$flowdata=subset(flowgage2$flowdata, flowgage2$flowdata$Qmm<48)

###########################################Observed
observed=rbind(flowgage1$flowdata,flowgage2$flowdata,flowgage$flowdata)


unlink("/dev/shm/rojakaveh/", recursive = TRUE)
basedirectory="~/multibasinpaper/" 
calibdirectory="/dev/shm/rojakaveh/"
dir.create(calibdirectory)
setwd(calibdirectory)
swat_project=c("LOC_Monkton","WBDC","EBDC" ) ##SWAT model initializations
site_no=c("4282629","4282586","04282581")
names(site_no)=swat_project
for (i in (swat_project)){
  dir.create(i)
  file.copy(list.files(paste0(basedirectory,i),full.names = TRUE,recursive = TRUE),paste0(calibdirectory,i),recursive = TRUE)
}

#-------------Multibasin Objective function-----------------------------
swat_objective=function(x, calib_range, calib_params, observed, rch,swat_project,calibdirectory,site_no,save_results=F)
{
  for (i in swat_project){
    #i="LOC_Monkton"
    setwd(paste0(calibdirectory,i))
    calib_params$current <- x
    tmpdir=as.character(as.integer((runif(1)+1)*10000))
    tmpdir=paste(c(format(Sys.time(), "%s"),tmpdir,Sys.getpid()),sep="",collapse="")
    print(tmpdir)
    dir.create(tmpdir)
    file.copy(list.files(),tmpdir)
    setwd(tmpdir)
    file.remove(list.files(pattern="output."))
    setup_swatcal(calib_params)
    alter_files(calib_params)
    libarch = if (nzchar(base::version$arch)) paste("libs", base::version$arch, sep = "/") else "libs"
    swatbin <- "rswat2012.exe"
    junkout=system(shQuote(paste(path.package("SWATmodel"), libarch, swatbin, sep = "/")),intern = TRUE)
    start_year = read.fortran(textConnection(readLines("file.cio")[9]), "f20")
    load("readSWAT.R")
    output_rch=readSWAT("rch", ".")
    output_rch=output_rch[RCH==rch]
    output_rch=output_rch[,c(5,6,7,52)]
    output_rch=tibble::add_column(output_rch, site_no=site_no[i])
    output_rch$FLOW_OUTmm=(output_rch$FLOW_OUTcms*24*3600)/(output_rch$AREAkm2*1000)
    assign(paste0("output_rch","_",i),output_rch)
    modeled=rbindlist(mget(ls(pattern = "output_rch_")))
    if(save_results){file.copy(list.files(),"../")}
    file.remove(list.files())
    setwd("../")
    file.remove(tmpdir)
  }
  test=merge(observed,modeled,by=c("mdate","site_no"))
  NS=NSeff(test$Qmm, test$FLOW_OUTmm)
  print(NS)
  return(abs(NS - 1))
}

swat_objective(x, calib_range, calib_params, observed, rch,swat_project,calibdirectory,site_no,save_results=F)
##DEoptim
NewoutDEoptimmulti<-DEoptim(swat_objective,calib_params$min,calib_params$max,
                            DEoptim.control(strategy = 6,NP = 16,itermax=300,parallelType = 1,
                                            packages = c("SWATmodel")),calib_range,calib_params,observed,rch,swat_project,calibdirectory,site_no)
