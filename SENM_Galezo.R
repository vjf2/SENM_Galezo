#SENM for Ali Galezo
#Created February 18, 2018
#Vivienne Foroughirad
#Modified May 16, 2018

library(adehabitatHR)
library(rgdal)
library(rgeos)
library(maptools)
library(spatstat)
library(raster)
library(parallel)
library(foreach)
library(doParallel)
library(igraph)
# library(devtools)
# install_github("vjf2/SocGen")

library(SocGen)

options(stringsAsFactors = FALSE)

source("SENM_helper_functions.R")

all_surveys<-read.csv("clean_all_surveys.csv")

#count up the number of sightings for each individual

all_surveys$sightings<-sapply(all_surveys$dolphin_id, function(i) length(all_surveys$dolphin_id[which(all_surveys$dolphin_id==i)]))

#Take the subset of surveys for individuals that have at least 15 sightings 

mod_surveys<-all_surveys[all_surveys$sightings>=15,]

#select focal juveniles 
#take all surveys from postweaning to age 10

juvs<-mod_surveys$dolphin_id[which(mod_surveys$life_stage=="postweaning" & mod_surveys$age<=10)]

focal_juvs<-data.frame(juv_id=unique(juvs), sightings=numeric(length(unique(juvs))))
focal_juvs$sightings<-sapply(focal_juvs$juv_id, function(i) length(juvs[which(juvs==i)]))
focal_juvs$sex<-life_history_lookup$sex[match(focal_juvs$juv_id, life_history_lookup$dolphin_id)]

#get juveniles which ares sexed and have at least 15 sightings as juveniles

focal_juvs<-focal_juvs[which(focal_juvs$sightings>=15
                             & focal_juvs$sex!=""),]  

#remove RAB since he died before age 4
focal_juvs<-focal_juvs[!focal_juvs$juv_id=="RAB",]

nj<-nrow(focal_juvs)

#Save list of focals
#write.csv(focal_juvs, "focal_juvs.csv", row.names=FALSE)

##run a spatially explicit null model for all individuals with more than 15 sightings
##individuals are available in the model from birth to death/last sighting

#project long lat to UTM with origin in shark bay
xydata<-cbind(mod_surveys$gps_east,mod_surveys$gps_south)
xydata2<-as.data.frame(project(xydata, "+proj=tmerc +lat_0=-25 +lon_0=113 +k=0.99999 +x_0=50000 +y_0=100000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
xydata3<-cbind(mod_surveys[,c("observation_date.x", "observation_id", "dolphin_id")],trunc(xydata2,0))
colnames(xydata3)<-c("Date", "observation_id", "dolphin_id","X","Y")

#create a grid on which to model animal home ranges
grid_buffer=5000

x <- seq(min(xydata3[,"X"])-grid_buffer,max(xydata3[,"X"])+grid_buffer,by=100) # where resolution is the pixel size you desire. 100 is the smallest i would go, if you make it larger you'll get coarser resolution, but faster runtimes
y <- seq(min(xydata3[,"Y"])-grid_buffer,max(xydata3[,"Y"])+grid_buffer,by=100)

xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE

#create UDs for each animal and extract the href smoothing parameter (need to manually select h for boundary method)

hrxydata<-SpatialPointsDataFrame(xydata3[,c("X","Y")],xydata3["dolphin_id"])

uds_href<-kernelUD(hrxydata[,1],grid=xy)

#create simplified coastline for the boundary, length of segments must be greater than 3*h
bound <- structure(list(x = c(122000,122000,116500,110000,108000), y = c(1000,10500,14500,20800,29800)), .Names = c("x", "y"))
bound <- do.call("cbind",bound)
Slo1 <- Line(bound)
Sli1 <- Lines(list(Slo1), ID="frontier1")
barrier <- SpatialLines(list(Sli1))

#if any smoothing parameters are too large for the boundary, set them to max allowed
maxh<-trunc(min(dist(bound))/3,0)

#pull out a list of individual smoothing parameters

hvalues<-list()
for (i in 1:length(uds_href)) {
  h<-uds_href[[i]]@h$h
  h<-ifelse(h>maxh,maxh,h)
  id<-names(uds_href)[[i]]
  hvalues[[i]]<-c(h, id)
}

h<-as.data.frame(do.call("rbind", hvalues), stringsAsFactors=FALSE)

names(h)<-c("h_opt", "dolphin_id")

h$h_opt<-as.numeric(h$h_opt)

#remove uds_href and do a memory clean-up

rm("uds_href");gc()

#recalculate the UDs with the boundary this time
optud<-list()

for (i in 1:dim(h)[1]){
  
  cdol<-hrxydata[hrxydata$dolphin_id==h$dolphin_id[i],]
  hopt<-h$h_opt[i]
  uds_man<-kernelUD(cdol,h=hopt,grid=xy, boundary=barrier)
  optud[[i]]<-uds_man
  cat(i)
}

##this will return some warnings, but they are safe to ignore for now

uddf<-unlist(optud)
class(uddf)<-"estUDm"

#read in a polygon on shark bay to do a fine scale trimming of the UDs
coast_polygon<-readOGR("coastpolygon", "coastpolygon")

coast_polygon<-spTransform(coast_polygon, CRS("+proj=tmerc +lat_0=-25 +lon_0=113 +k=0.99999 +x_0=50000 +y_0=100000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

udsgdf <- as(estUDm2spixdf(uddf),"SpatialGridDataFrame")

#use coast polygon as mask for spatial grid

rgrid <- raster(udsgdf)
#plot(rgrid)
rgrid_msk <- mask(rgrid,coast_polygon, inverse=TRUE)

grid_ae <- as(rgrid_msk, 'SpatialGridDataFrame')
gridded(grid_ae) <- TRUE
grid_ae[[1]] <- as.numeric(!is.na(grid_ae[[1]])) 

#multiply each column of udsgdf by the mask, and restandardize the percentages so that each set of probabilities sums to 1

resu <- lapply(1:ncol(udsgdf), function(i) {udsgdf[[i]] * grid_ae[[1]] / sum(udsgdf[[i]] * grid_ae[[1]]) })

resu <- as.data.frame(resu) 
names(resu) <- names(udsgdf@data) 

# and define it as data slot for udsgdf 
udsgdf@data <- resu 

#now we have our home ranges, create polygons for daily survey effort from survey coordinates
#mcps have 1km buffers, and always include the launch point as one of the vertices.

groupings<-unique(xydata3[,c("Date", "observation_id", "X", "Y")])

days<-split(groupings, groupings$Date)

#add launch to each days surveys (4 so that each day has a least 5 points)

launch<-data.frame(Date=as.Date(rep("2001-01-01",4)),
                   observation_id=rep("launch",4),
                   X=rep(122241,4),
                   Y=seq(11963,11966,1))

days1<-lapply(days, function(x) add_launch(x)) ##add launch in digidolph helper functions

survey_days<-as.data.frame(do.call("rbind", days1))

survey_days[,c("X","Y")]<-apply(survey_days[,c("X","Y")],2, as.numeric)

daily_xydata<-SpatialPointsDataFrame(survey_days[,c("X","Y")],survey_days["Date"])

mcps<-mcp(daily_xydata[,1], percent=100, unin=c("m"), unout=c("m2"))

#Add buffer, make sure whole area is covered

buff_days<-gBuffer(mcps, byid=TRUE,width=1000)

#Number of animals in study
n<-length(unique(xydata3$dolphin_id))

#Number of survey days
d<-length(unique(xydata3$Date))
dates<-sort(unique(xydata3$Date))

#Get availability matrix which has the birthdate and deathdate / last sighting for each dolphin

dolphins<-sort(unique(xydata3$dolphin_id))

matnames<-list(as.character(dates),dolphins)

##create availability list of all individuals 

fast_avail<-data.frame(dolphin_id=dolphins, entry=rep(NA, length(dolphins)), depart=rep(NA, length(dolphins)))

fast_avail$entry<-as.Date(merge(fast_avail, life_history_lookup, by="dolphin_id")$birth_date)

fast_avail$death<-as.Date(merge(fast_avail, life_history_lookup, by="dolphin_id")$death_date)

#If no death date add date of last sighting

lsd_dol<-split(xydata3, as.factor(xydata3$dolphin_id))

silent<-lapply(lsd_dol, function(x) {fast_avail$depart[match(x$dolphin_id[1], fast_avail$dolphin_id)]<<-as.character(max(x$Date))})

fast_avail$death<-as.character(fast_avail$death)

fast_avail$depart<-ifelse(is.na(fast_avail$death), fast_avail$depart, fast_avail$death)

fast_avail$depart<-as.Date(fast_avail$depart, format="%Y-%m-%d")

#if no death then add 6 months to depart

fast_avail$depart<-ifelse(is.na(fast_avail$death), 
                          fast_avail$depart+(6*30.5),
                          fast_avail$depart)

fast_avail$depart<-as.Date(fast_avail$depart, origin="1970-01-01")

fast_avail<-fast_avail[,c(1:3)]

#Give the model a schedule with TRUE/FALSE for each dolphin's avaiability on each survey day

fast_avail$entry<-as.numeric(fast_avail$entry)
fast_avail$depart<-as.numeric(fast_avail$depart)

alive<-Vectorize(FUN=function(r,c) 
  isTRUE(r>=fast_avail$entry[which(fast_avail$dolphin_id==c)] 
         & r<=fast_avail$depart[which(fast_avail$dolphin_id==c)]))

schedule<-outer(as.numeric(dates), dolphins, FUN=alive) #takes 2 min to run

dimnames(schedule)<-matnames

dolphin_density_per_km<-dim(xydata3)[1]/(sum(area(buff_days))/1000000)

#Now run the model, first we'll make a list of simulated allsurveys
#I would run ~3 simulations, fit the gprox, then go back and run a full set
#Note it takes about 4 minutes to run one full 30-year simulation

fullgrid(udsgdf)<-FALSE #for efficient subsetting

#cleanup environment

keep<-c("d", "buff_days", "udsgdf", "dolphin_density_per_km", "schedule",
        "dolphins", "xydata3", "coast_polygon", "cmp_fast_random_points",
        "focal_juvs", "life_history_lookup", "dolphins", "nj")

rm(list=setdiff(ls(),keep))

gc()

# save.image(file="simready.RData")
# load("simready.RData")

#Calculate number of dolphins seen each day

areakm<-area(buff_days)/1000000
numdol<-round(areakm*dolphin_density_per_km)
numdol<-ifelse(numdol<=1, 2, numdol) 


num_sim=1000 #number of simulations to run

#Set up cluster for parallelization
#Timing depends on number of cores available
#7 threads - 1000 sims in 8 hours
source("SENM_helper_functions.R")
library(pbapply)

cl<-makeCluster(detectCores()-1)
clusterEvalQ(cl, library(maptools))
clusterExport(cl, c("d", "buff_days", "udsgdf", "schedule", "num_sim","cmp_fast_random_points", "numdol"))

starttime<-Sys.time()

  nest_days<-pblapply(seq_len(d), FUN=function(i){
    
    bound<-buff_days[i,]
    nd<-numdol[i]
    dailygrid<-udsgdf[bound,,drop=TRUE] 
    probweights<-colSums(dailygrid@data, na.rm=TRUE)
    probweights<-probweights[names(probweights) %in% colnames(schedule)[schedule[i,]==TRUE]]
    dc<-coordinates(dailygrid)
    dgdf<-dailygrid@data
    holder<-replicate(num_sim, cmp_fast_random_points(probweights = probweights, 
                                                  nd=nd, 
                                                  dc=dc,
                                                  dgdf=dgdf), 
                      simplify=FALSE)
    return(holder)
  }, cl=cl)

endtime<-Sys.time()

stopCluster(cl)

endtime-starttime #check run time

#rearrange list so that days are listed within iteration 
sim_surveys<-sapply(1:num_sim, function(i) lapply(nest_days, "[[", i), simplify = FALSE)

#save(sim_surveys, file="100juvs.RData")

rm(nest_days)

#Calculate mean group size in real data (xydata3)

mean_group_size<-mean(table(xydata3$observation_id))

#Group dolphins together using k-means clustering to get same average group size in real data

#First use the real data to estimate how many groups we see each day

dpd<-sapply(unique(xydata3$Date), function(i) length(xydata3$Date[which(xydata3$Date==i)]))

gpd<-sapply(unique(xydata3$Date), function(i) length(unique(xydata3$observation_id[which(xydata3$Date==i)])))

prob_per_day<-data.frame(gpd,dpd)

prob_per_day$gs<-prob_per_day$dpd/prob_per_day$gpd

mod<-lm(gpd~dpd, data=prob_per_day)

#Set up cluster to calculate groups, (should take about 35 min for 1000 sims of full data)

bid<-buff_days$id

cl<-makeCluster(detectCores()-1)
clusterExport(cl, c("sim_surveys", "mod", "bid"))
starttime<-Sys.time()

kfinal<-pblapply(seq_len(num_sim), FUN=function(w) {
  single_sim<-sim_surveys[[w]]
  counter=1
  each_days_assoc<-lapply(single_sim, function(x1) {

    num_clust<-(mod$coefficients[1]+(nrow(x1)*mod$coefficients[2]))+sample(mod$residuals, 1)
    num_clust<-ifelse(num_clust>nrow(x1), nrow(x1)-1, num_clust)
    num_clust<-ifelse(num_clust<1, 1, num_clust)
    num_clust<-ifelse(nrow(x1[!duplicated(x1[,1:2]),])<num_clust,nrow(x1[!duplicated(x1[,1:2]),]),num_clust)
    
    xk<-kmeans(dist(x1[,1:2]),num_clust)
    # gprox<-mean(xk$centers[xk$cluster])
    dayAssocK<-data.frame(IDs=x1[,3],Group=as.numeric(xk$cluster))
    dayAssocK$Permutation<-c(rep(as.character(bid[counter]),nrow(dayAssocK)))
    dayAssocK$Permutation<-as.Date(dayAssocK$Permutation, format=c("%Y-%m-%d"))
    counter <<- counter + 1
    return(dayAssocK)
  })
  eda<-do.call("rbind", each_days_assoc)
  
  eda$id<-paste0(eda$Permutation,"_",eda$Group)
  
  return(eda)
}, cl=cl)

endtime<-Sys.time()

endtime-starttime #check run time

stopCluster(cl)

random_group_sizes<-lapply(kfinal, function(x) mean(table(x$id)))

#Create a table of availability dates for focals

availability_ego<-focal_juvs
availability_ego$entry<-life_history_lookup$weaning_date[match(availability_ego$juv_id,life_history_lookup$dolphin_id)]
availability_ego$depart<-life_history_lookup$birth_date[match(availability_ego$juv_id,life_history_lookup$dolphin_id)]+(10*365.25)
names(availability_ego)[1]<-"dolphin_id"

#Create a table of availability dates for nonfocals (alters)

availability_alter<-data.frame(dolphin_id=dolphins)
availability_alter$entry<-life_history_lookup$birth_date[match(availability_alter$dolphin_id,life_history_lookup$dolphin_id)]+(4*365.25)
availability_alter$depart<-life_history_lookup$birth_date[match(availability_alter$dolphin_id,life_history_lookup$dolphin_id)]+(10*365.25)

availability_ego$entry<-as.numeric(availability_ego$entry)
availability_ego$depart<-as.numeric(availability_ego$depart)
availability_alter$entry<-as.numeric(availability_alter$entry)
availability_alter$depart<-as.numeric(availability_alter$depart)

xydata3$Date<-as.numeric(xydata3$Date)

#Calculate association indices for indiviudals in real data, need
##to calculate the matrix once per focal to allow individuals
##to have different availability ranges as egos vs alters

ai_egos<-list()

lookup<-as.data.frame(t(combn(dolphins, 2)))

for (n in 1:nrow(availability_ego)) {
  
  ego<-as.character(availability_ego$dolphin_id[n])
  start<-availability_ego$entry[n]
  end<-availability_ego$depart[n]
  
  availability_subset<-availability_alter
  availability_subset$entry[availability_subset$dolphin_id==ego]<-start
  availability_subset$depart[availability_subset$dolphin_id==ego]<-end
  
  lookup$start<-availability_subset$entry[match(lookup[,1],availability_subset$dolphin_id)]
  lookup$end<-availability_subset$depart[match(lookup[,1],availability_subset$dolphin_id)]
  lookup$start2<-availability_subset$entry[match(lookup[,2],availability_subset$dolphin_id)]
  lookup$end2<-availability_subset$depart[match(lookup[,2],availability_subset$dolphin_id)]
  lookup$hstart<-with(lookup, ifelse(start>start2, start, start2))
  lookup$hend<-with(lookup, ifelse(end<end2, end, end2))
  lookup$tp<-with(lookup, hend-hstart)
  
  lookup_run<-lookup[which(lookup$tp>1),]
  

  ego_network<-xydata3[xydata3$Date>=start & xydata3$Date<=end,]
  ego_survey_ids<-ego_network$observation_id[ego_network$dolphin_id==ego]
  
  alters<-as.character(ego_network$dolphin_id[ego_network$observation_id %in% ego_survey_ids])
  
  lookup_run_ego<-lookup_run[which(lookup_run$V1 %in% alters & lookup_run$V2 %in% alters),]
  
  lookup_run_ego<-lookup_run_ego[which(lookup_run_ego$hstart<end & lookup_run_ego$hend>start),]
  
  network_ego<-cmp_ai_filtered(sightings=ego_network,
                               group_variable="observation_id", 
                               dates="Date", 
                               IDs="dolphin_id", 
                               lookup_run_ego=lookup_run_ego)
  
  ai_egos[[n]]<-network_ego
  names(ai_egos)[n]<-ego
  cat(n)
}

real_ai_egos<-ai_egos

save(real_ai_egos, file="real_ai_egos.RData")

#Read in genomic and matrilineal relatedness data files

relatedness<-read.csv("max_likelihood_relatedness.csv")
kindat_pos<-read.csv("kindat_pos.csv")

#Calculate network metrics for the real data

nj<-length(real_ai_egos)

network_metrics<-data.frame(ego=character(nj), 
                            mixstrength=numeric(nj),
                            mixdegree=numeric(nj),
                            mixcc=numeric(nj),
                            ss_strength=numeric(nj),
                            ss_cc=numeric(nj),
                            ss_degree=numeric(nj),
                            os_strength=numeric(nj),
                            os_cc=numeric(nj),
                            os_degree=numeric(nj),
                            close_kin=numeric(nj),
                            known_non_kin=numeric(nj),
                            available_kin=numeric(nj),
                            sex=character(nj))

for (i in 1:length(ai_egos)) {
  
  ego<-names(ai_egos[i])
  network_metrics[i,"ego"]<-ego
  focal_sex<-life_history_lookup$sex[match(ego, life_history_lookup$dolphin_id)]
  network_metrics[i,"sex"]<-focal_sex
  
  m<-as.matrix(ai_egos[[i]])
  
  m[lower.tri(m)]=t(m)[lower.tri(m)]
  
  m[is.nan(m)]<-0 #if neither was sighted during overlap set to 0
  
  #remove inds that are NA for ego (and add ego back in)
  inds<-sort(c(colnames(m)[which(!is.na(m[ego,]))], ego))
  
  m<-m[inds, inds]
  
  g<-graph.adjacency(m, mode="undirected", weighted=TRUE, diag=FALSE)
  
  eg<-make_ego_graph(g, order=1, nodes=ego)[[1]] #is a list, return element 1
  
  if(length(E(eg))==0){next} else{
    
    #add sexes of individuals
    
    V(eg)$sex<-life_history_lookup$sex[match(V(eg)$name, life_history_lookup$dolphin_id)]
    
    network_metrics[i, "mixdegree"]<-degree(eg, ego)
    network_metrics[i, "mixstrength"]<-strength(eg, ego)
    network_metrics[i, "mixcc"]<-transitivity(eg, type="local", vids=ego)
    
    #pull out just same sex
    
    egss<-induced_subgraph(eg, vids=V(eg)$name[which(V(eg)$sex==focal_sex)])
    
    network_metrics[i, "ss_degree"]<-degree(egss, ego)
    network_metrics[i, "ss_strength"]<-strength(egss, ego)
    network_metrics[i, "ss_cc"]<-transitivity(egss, type="local", vids=ego)
    
    #pull out just opposite sex (and ego)
      
    opposite_sex_names <- c(ego, V(eg)$name[which(V(eg)$sex!=focal_sex)])
    egops<-induced_subgraph(eg, vids=V(eg)$name %in% opposite_sex_names)
      
    network_metrics[i, "os_degree"]<-degree(egops, ego)
    network_metrics[i, "os_strength"]<-strength(egops, ego)
    network_metrics[i, "os_cc"]<-transitivity(egops, type="local", vids=ego)
    
    #assign relatedness to edges, and an unknown number
    
    seq_edges<-as.numeric(E(eg) [ from(ego) ])
    
    edges<-sapply(seq_edges, function(x) V(eg)[inc(x)]$name)
    el<-as.data.frame(t(edges))
    
    el<-merge_pairs(el, relatedness[,c("xID1", "ID2", "relatedness")], 
                    "V1", "V2", "xID1", "ID2", all.x=TRUE, all.y=FALSE)
    
    el<-merge_pairs(el, kindat_pos, "V1", "V2", "ID1", "ID2", all.x=TRUE, all.y=FALSE)
    
    el$relatedness<-ifelse(is.na(el$relatedness), el$matpedR, el$relatedness)
    
    close_kin<-length(na.omit(el$relatedness[el$relatedness>=0.125]))
    network_metrics[i,"close_kin"]<-close_kin
    known_non_kin<-length(el$relatedness[!is.na(el$relatedness)])
    network_metrics[i,"known_non_kin"]<-known_non_kin
    
  }
}

write.csv(network_metrics, "network_metrics.csv", row.names = FALSE)

real_network_metrics<-read.csv("network_metrics.csv")

#Repeat for the results of the random model

#load("kfinal1000.RData")

#Format the columns for faster and more consistent subsetting

kfinal<-lapply(kfinal, function(x) {names(x)<-c("dolphin_id", "DayGroup", "Date", "observation_id");x})
kfinal<-lapply(kfinal, function(x) {x[,"observation_id"]<-as.numeric(as.factor(x[,"observation_id"]));x})
kfinal<-lapply(kfinal, function(x) {x[,"Date"]<-as.numeric(x[,"Date"]);x})
kfinal<-lapply(kfinal, function(x) {x[,"dolphin_id"]<-as.character(x[,"dolphin_id"]);x})

lookup<-as.data.frame(t(combn(dolphins, 2)))

#Set up a cluster to run in parallel
#This is the longest section, it may take up to 8 hours on a laptop

cl<-makeCluster(detectCores()-1)
clusterExport(cl, c("kfinal", "availability_ego", "availability_alter", "cmp_ai_filtered"))
registerDoParallel(cl)

starttime<-Sys.time()

ai_egos_rand<-foreach (n=1:nrow(availability_ego), .errorhandling = "pass") %dopar% {
  
  ego<-as.character(availability_ego$dolphin_id[n])
  start<-availability_ego$entry[n]
  end<-availability_ego$depart[n]
  
  availability_subset<-availability_alter
  availability_subset$entry[availability_subset$dolphin_id==ego]<-start
  availability_subset$depart[availability_subset$dolphin_id==ego]<-end
  
  lookup$start<-availability_subset$entry[match(lookup[,1],availability_subset$dolphin_id)]
  lookup$end<-availability_subset$depart[match(lookup[,1],availability_subset$dolphin_id)]
  lookup$start2<-availability_subset$entry[match(lookup[,2],availability_subset$dolphin_id)]
  lookup$end2<-availability_subset$depart[match(lookup[,2],availability_subset$dolphin_id)]
  lookup$hstart<-with(lookup, ifelse(start>start2, start, start2))
  lookup$hend<-with(lookup, ifelse(end<end2, end, end2))
  lookup$tp<-with(lookup, hend-hstart)
  
  lookup_run<-lookup[which(lookup$tp>1),]
  ai_egos<-list()
  
  for (j in 1:num_sim) {
    
    
    random1<-kfinal[[j]]
    ego_network<-random1[random1$Date>=start & random1$Date<=end,]
    ego_survey_ids<-ego_network$observation_id[ego_network$dolphin_id==ego]
    
    alters<-as.character(ego_network$dolphin_id[ego_network$observation_id %in% ego_survey_ids])
    
    lookup_run_ego<-lookup_run[which(lookup_run$V1 %in% alters & lookup_run$V2 %in% alters),]
    
    lookup_run_ego<-lookup_run_ego[which(lookup_run_ego$hstart<end & lookup_run_ego$hend>start),]
    
    if(nrow(lookup_run_ego)==0) {ai_egos[[j]]<-NA
      next}
    
    network_ego<-cmp_ai_filtered(sightings=ego_network,
                                 group_variable="observation_id", 
                                 dates="Date", 
                                 IDs="dolphin_id", 
                                 lookup_run_ego=lookup_run_ego)
    
    ai_egos[[j]]<-network_ego
    
  }
   
  ai_egos
  
}


stopCluster(cl)
endtime<-Sys.time()

endtime-starttime #check run time

names(ai_egos_rand)<-availability_ego$dolphin_id

save(ai_egos_rand, file="ai_egos_rand.RData")

#Calculate network metrics for the random data

random_network_metrics<-list()

for (k in 1:length(ai_egos_rand)){
  
  network_metrics<-data.frame(ego=character(num_sim), 
                              mixstrength=numeric(num_sim),
                              mixdegree=numeric(num_sim),
                              mixcc=numeric(num_sim),
                              ss_strength=numeric(num_sim),
                              ss_cc=numeric(num_sim),
                              ss_degree=numeric(num_sim),
                              os_strength=numeric(num_sim),
                              os_cc=numeric(num_sim),
                              os_degree=numeric(num_sim),
                              close_kin=numeric(num_sim),
                              known_non_kin=numeric(num_sim),
                              available_kin=numeric(num_sim),
                              sex=character(num_sim),
                              iteration=numeric(num_sim))
  
  ai_egos1<-ai_egos_rand[[k]]

    for (i in 1:length(ai_egos1)) {
    
    ego<-names(ai_egos_rand[k])
    network_metrics[i,"ego"]<-ego
    focal_sex<-life_history_lookup$sex[match(ego, life_history_lookup$dolphin_id)]
    network_metrics[i,"sex"]<-focal_sex
    network_metrics[i,"iteration"]<-i
    
    m<-as.matrix(ai_egos1[[i]])
    m[lower.tri(m)]=t(m)[lower.tri(m)]
    m[is.nan(m)]<-0 #if neither was sighted during overlap set to 0
    
    if(all(is.na(m))){next}
    
    #remove inds that are NA for ego (and add ego back in)
    inds<-sort(c(colnames(m)[which(!is.na(m[ego,]))], ego))
    
    m<-m[inds, inds]
    
    g<-graph.adjacency(m, mode="undirected", weighted=TRUE, diag=FALSE)
    
    eg<-make_ego_graph(g, order=1, nodes=ego)[[1]] #is a list, return element 1
    
    if(length(E(eg))==0){next} else{
      
    #add sexes of individuals
    
    V(eg)$sex<-life_history_lookup$sex[match(V(eg)$name, life_history_lookup$dolphin_id)]
    
    network_metrics[i, "mixdegree"]<-degree(eg, ego)
    network_metrics[i, "mixstrength"]<-strength(eg, ego)
    network_metrics[i, "mixcc"]<-transitivity(eg, type="local", vids=ego)
    
    #pull out just same sex
    
    egss<-induced_subgraph(eg, vids=V(eg)$name[which(V(eg)$sex==focal_sex)])
    
    network_metrics[i, "ss_degree"]<-degree(egss, ego)
    network_metrics[i, "ss_strength"]<-strength(egss, ego)
    network_metrics[i, "ss_cc"]<-transitivity(egss, type="local", vids=ego)
      
    #pull out just opposite sex (and ego)
    
    opposite_sex_names <- c(ego, V(eg)$name[which(V(eg)$sex!=focal_sex)])
    egops<-induced_subgraph(eg, vids=V(eg)$name %in% opposite_sex_names)
    
    network_metrics[i, "os_degree"]<-degree(egops, ego)
    network_metrics[i, "os_strength"]<-strength(egops, ego)
    network_metrics[i, "os_cc"]<-transitivity(egops, type="local", vids=ego)
    
    #assign relatedness to edges, and an unknown number
    
    seq_edges<-as.numeric(E(eg) [ from(ego) ])
    
    edges<-sapply(seq_edges, function(x) V(eg)[inc(x)]$name)
    el<-as.data.frame(t(edges))
    
    el<-merge_pairs(el, relatedness[,c("xID1", "ID2", "relatedness")], 
                    "V1", "V2", "xID1", "ID2", all.x=TRUE, all.y=FALSE)
    
    el<-merge_pairs(el, kindat_pos, "V1", "V2", "ID1", "ID2", all.x=TRUE, all.y=FALSE)
    
    el$relatedness<-ifelse(is.na(el$relatedness), el$matpedR, el$relatedness)
    
    close_kin<-length(na.omit(el$relatedness[el$relatedness>=0.125]))
    network_metrics[i,"close_kin"]<-close_kin
    known_non_kin<-length(el$relatedness[!is.na(el$relatedness)])
    network_metrics[i,"known_non_kin"]<-known_non_kin
    
    }
}

  random_network_metrics[[k]]<-network_metrics
  cat(k)
  
}

all_random_metrics<-do.call("rbind", random_network_metrics)

write.csv(all_random_metrics, "all_random_metrics_kmeanswJUI.csv", row.names = FALSE)

#Calculate available kin in population, available non-kin, and available individuals of unknown relatedness

focals<-availability_ego[,"dolphin_id"]
alters<-availability_alter[,"dolphin_id"]

###change to rlookup

rlookup<-data.frame(Var1=rep(focals, each=length(alters)), Var2=rep(alters, length(focals)))
#remove self pairs

rlookup[,2]<-ifelse(rlookup[,1]==rlookup[,2], NA, rlookup[,2])
rlookup<-rlookup[complete.cases(rlookup),]

rlookup$start<-availability_ego$entry[match(rlookup[,1],availability_ego$dolphin_id)]
rlookup$end<-availability_ego$depart[match(rlookup[,1],availability_ego$dolphin_id)]
rlookup$start2<-availability_alter$entry[match(rlookup[,2],availability_alter$dolphin_id)]
rlookup$end2<-availability_alter$depart[match(rlookup[,2],availability_alter$dolphin_id)]
rlookup$hstart<-as.Date(with(rlookup, ifelse(start>start2, start, start2)),origin="1970-01-01")
rlookup$hend<-as.Date(with(rlookup, ifelse(end<end2, end, end2)),origin="1970-01-01")
rlookup$tp<-with(rlookup, hend-hstart)

rlookup_run<-rlookup[which(rlookup$tp>1),]

#Add relatedness data to list of all possible pairs

akin<-merge_pairs(rlookup_run[,c("Var1", "Var2", "tp")], relatedness[,c("xID1", "ID2", "relatedness")], "Var1", "Var2", "xID1", "ID2", all.x=TRUE, all.y=FALSE)

akin<-merge_pairs(akin, kindat_pos, "Var1", "Var2", "ID1", "ID2", all.x=TRUE, all.y=FALSE)

akin$relatedness<-ifelse(is.na(akin$relatedness), akin$matpedR, akin$relatedness)

dolphs<-split(akin, akin$Var1)

available_close_kin<-unlist(lapply(dolphs, function (x) length(na.omit(x$relatedness[x$relatedness>=0.125]))))

available_non_kin<-unlist(lapply(dolphs, function (x) length(na.omit(x$relatedness[x$relatedness<0.125]))))

available_unknown<-unlist(lapply(dolphs, function (x) length(x$relatedness[is.na(x$relatedness)])))

#Add available kin to real data

real_network_metrics$available_kin<-available_close_kin[match(real_network_metrics$ego, names(available_close_kin))]
real_network_metrics$available_non_kin<-available_non_kin[match(real_network_metrics$ego, names(available_non_kin))]
real_network_metrics$available_unknown<-available_unknown[match(real_network_metrics$ego, names(available_unknown))]

real_network_metrics$percent_close_kin<-real_network_metrics$close_kin/real_network_metrics$available_kin

#Add available kin to random data as well 

all_random_metrics$available_kin<-real_network_metrics$available_kin[match(all_random_metrics$ego, real_network_metrics$ego)]
all_random_metrics$available_non_kin<-real_network_metrics$available_non_kin[match(all_random_metrics$ego, real_network_metrics$ego)]
all_random_metrics$available_unknown<-real_network_metrics$available_unknown[match(all_random_metrics$ego, real_network_metrics$ego)]

all_random_metrics$percent_close_kin<-all_random_metrics$close_kin/all_random_metrics$available_kin

#write out updated network files

write.csv(real_network_metrics, "real_network_metrics.csv")
write.csv(all_random_metrics, "all_random_metrics.csv")

####See Figure Plotting for figures and aggregating results

save(kfinal, file="kfinal1000.RData")

