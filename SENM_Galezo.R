#SENM for Ali Galezo
#Created February 18, 2018
#Vivienne Foroughirad
#Modified January 25, 2019

library(adehabitatHR)
library(rgdal)
library(rgeos)
library(raster)
library(parallel)
library(foreach)
library(doParallel)
library(igraph)
# library(devtools)
# install_github("vjf2/SocGen")

library(SocGen)

options(stringsAsFactors = FALSE)

all_surveys<-read.csv("clean_all_surveys.csv")

#read and format life history data

life_history<-read.delim("Raw_input_files/LifeHistory_20180507.txt", header=TRUE, sep="\t")
life_history$birth_date<-as.Date(life_history$birth_date, format="%Y-%m-%d")
life_history$weaning_date<-as.Date(life_history$weaning_date, format="%Y-%m-%d")
life_history$death_date<-as.Date(life_history$death_date, format="%Y-%m-%d")
life_history_lookup<-life_history[!duplicated(life_history$dolphin_id),]

#count up the number of sightings for each individual

sightings<-table(all_surveys$dolphin_id)

#Take the subset of surveys for individuals that have at least 15 sightings 

mod_surveys<-all_surveys[all_surveys$dolphin_id %in% names(sightings)[which(sightings>=15)],]

#select focal juveniles 
#take all surveys from postweaning to age 10

juvs<-mod_surveys$dolphin_id[which(mod_surveys$life_stage=="postweaning" & mod_surveys$age<=10)]

focal_juvs<-as.data.frame(table(juvs))

focal_juvs$sex<-life_history_lookup$sex[match(focal_juvs$juvs, life_history_lookup$dolphin_id)]

#get juveniles which ares sexed and have at least 15 sightings as juveniles

focal_juvs<-focal_juvs[which(focal_juvs$Freq>=15
                             & focal_juvs$sex!=""),]  

#remove RAB since he died before age 4
focal_juvs<-focal_juvs[!focal_juvs$juvs=="RAB",]

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

#add a set of launch site locations to each days surveys (4 so that each day has a least 5 points)

launch<-data.frame(Date=as.Date(rep("2001-01-01",4)),
                   observation_id=rep("launch",4),
                   X=rep(122241,4),
                   Y=seq(11963,11966,1))

days1<-lapply(days, function(x) {
                launch[1]<-x[1,1]
                x<-rbind(x, launch)
                return(x)}) 

survey_days<-as.data.frame(do.call("rbind", days1))

survey_days[,c("X","Y")]<-apply(survey_days[,c("X","Y")],2, as.numeric)

daily_xydata<-SpatialPointsDataFrame(survey_days[,c("X","Y")],survey_days["Date"])

mcps<-mcp(daily_xydata[,1], percent=100, unin=c("m"), unout=c("m2"))

#Add buffer, make sure whole area is covered

buff_days<-gBuffer(mcps, byid=TRUE,width=1000)
buff_days<-gSimplify(buff_days, tol=50)

#Number of animals in study
n<-length(unique(xydata3$dolphin_id))

#Number of survey days
d<-length(unique(xydata3$Date))
dates<-sort(unique(xydata3$Date))

#Get availability matrix which has the birthdate and deathdate / last sighting for each dolphin

dolphins<-sort(unique(xydata3$dolphin_id))

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

schedule<-schedulize(fast_avail, id="dolphin_id", start="entry", end="depart", 
                     dates=dates, format="sim")

dolphin_density_per_km<-dim(xydata3)[1]/(sum(area(buff_days))/1000000)

fullgrid(udsgdf)<-FALSE #convert to SpatialPixels for more efficient subsetting

#cleanup environment

keep<-c("d", "buff_days", "udsgdf", "schedule",
        "dolphins", "xydata3", "coast_polygon","focal_juvs", 
        "life_history_lookup", "dolphins", "nj", "dates")

rm(list=setdiff(ls(),keep))

gc()

# save.image(file="simready.RData")
# load("simready.RData")

# areakm<-area(buff_days)/1000000
# numdol<-round(areakm*dolphin_density_per_km)
# numdol<-ifelse(numdol<=1, 2, numdol) 

numdol<-c(table(xydata3$Date))

num_sim=1000 #number of simulations to run
gridrad<-udsgdf@grid@cellsize[1]/2

#Set up cluster for parallelization
#Timing depends on number of cores available
#7 threads - 1000 sims in 8 hours

cl<-makeCluster(detectCores()-1)
clusterEvalQ(cl, {library(sp);library(SocGen)})
clusterExport(cl, c("d", "buff_days", "udsgdf", "schedule", "num_sim", "numdol", "gridrad"))

starttime<-Sys.time()

  nest_days<-parLapplyLB(cl=cl, seq_len(d), fun=function(i){
    
    bound<-buff_days[i,]
    nd<-numdol[i]
    dailygrid<-udsgdf[bound,,drop=TRUE] 
    probweights<-colSums(dailygrid@data, na.rm=TRUE)
    probweights<-probweights[names(probweights) %in% colnames(schedule)[schedule[i,]==TRUE]]
    dc<-coordinates(dailygrid)
    dgdf<-dailygrid@data
    holder<-replicate(num_sim, fast_random_points(probweights = probweights,
                                                  nd=nd,
                                                  dc=dc,
                                                  dgdf=dgdf,
                                                  gridrad=gridrad),
                      simplify=FALSE)
    return(holder)
  })

endtime<-Sys.time()

stopCluster(cl)

endtime-starttime #check run time

#rearrange list so that days are listed within iteration 
sim_surveys<-sapply(1:num_sim, function(i) lapply(nest_days, "[[", i), simplify = FALSE)

#save(sim_surveys, file="100juvs.RData")

rm(nest_days)

#Calculate mean group size in real data (xydata3)

mean_group_size<-mean(table(xydata3$observation_id))

#Group dolphins together using hclust clustering to get same average group size in real data

groupperday<-table(xydata3$Date[!duplicated(xydata3$observation_id)])

sim_surveys<-lapply(sim_surveys, function(i) lapply(1:length(i), function(q) {
  names(i[[q]])<-c("y", "x", "id")
  return(i[[q]])}))

kfinal<-group_assign(data=sim_surveys[1:10], id="id", xcoord ="x", ycoord="y",
                     time = names(groupperday),group_vector=groupperday, method="hclust")

save(kfinal, file="kfinal1000.RData")

rm(sim_surveys)

random_group_sizes<-lapply(kfinal, function(x) mean(table(x$observation_id)))

#Create a table of availability dates for focals

availability_ego<-focal_juvs
availability_ego$entry<-life_history_lookup$weaning_date[match(availability_ego$juvs,life_history_lookup$dolphin_id)]
availability_ego$depart<-life_history_lookup$birth_date[match(availability_ego$juvs,life_history_lookup$dolphin_id)]+(10*365.25)
names(availability_ego)[1]<-"dolphin_id"

#Create a table of availability dates for nonfocals (alters)

availability_alter<-data.frame(dolphin_id=dolphins)
availability_alter$entry<-life_history_lookup$birth_date[match(availability_alter$dolphin_id,life_history_lookup$dolphin_id)]+(4*365.25)
availability_alter$depart<-life_history_lookup$birth_date[match(availability_alter$dolphin_id,life_history_lookup$dolphin_id)]+(12*365.25)

#Calculate association indices for individuals in real data, need
##to calculate the matrix once per focal to allow individuals
##to have different availability ranges as egos vs alters

xydata3$Date<-as.Date(xydata3$Date, origin="1970-01-01")
xydata3<-xydata3[order(xydata3$Date),]

ai_mask<-schedulize(availability_alter, dates=dates, format="mask")

ai_egos<-vector(mode = "list", nrow(availability_ego))

for (n in 1:nrow(availability_ego)) {
  
  ego<-as.character(availability_ego$dolphin_id[n])
  start<-availability_ego$entry[n]
  end<-availability_ego$depart[n]
  
  ego_network<-xydata3[xydata3$Date>=start & xydata3$Date<=end,]
  
  eid<-unique(ego_network$observation_id[which(ego_network$dolphin_id==ego)])
  dolls<-unique(ego_network$dolphin_id[which(ego_network$observation_id %in% eid)])
  ego_network<-ego_network[which(ego_network$dolphin_id %in% dolls),]
  
  ego_mask<-ai_mask
  
  ego_mask[ego,]<-1
  
  network_ego<-simple_ratio(sightings=ego_network,
                            group_variable="observation_id", 
                            dates="Date", 
                            IDs="dolphin_id", 
                            symmetric=FALSE, 
                            mask=ego_mask)
  
  ai_egos[[n]]<-network_ego
  names(ai_egos)[n]<-ego
  cat(n)
}

real_ai_egos<-ai_egos

save(real_ai_egos, file="real_ai_egos.RData")

library(igraph)

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
                            ss_strength_kin=numeric(nj),
                            ss_strength_kin_prop=numeric(nj),
                            n_ss_kin=numeric(nj),
                            ss_cc=numeric(nj),
                            ss_degree=numeric(nj),
                            os_strength=numeric(nj),
                            os_strength_kin=numeric(nj),
                            os_strength_kin_prop=numeric(nj),
                            n_os_kin=numeric(nj),
                            os_cc=numeric(nj),
                            os_degree=numeric(nj),
                            sex=character(nj))

for (i in 1:length(ai_egos)) {
  
  ego<-names(ai_egos[i])
  network_metrics[i,"ego"]<-ego
  focal_sex<-life_history_lookup$sex[match(ego, life_history_lookup$dolphin_id)]
  network_metrics[i,"sex"]<-focal_sex
  
  m<-as.matrix(ai_egos[[i]])
  
  m[lower.tri(m)]=t(m)[lower.tri(m)]
  
  m[is.nan(m)]<-0 #if neither was sighted during overlap set to 0
  
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
    
    #pull out animals with relatedness >= 0.125
    relatives <- unique(c(el[which(el$relatedness >= 0.125),]$V1, el[which(el$relatedness >= 0.125),]$V2))
    
    #same-sex strength, kin only:
    egsskin<-induced_subgraph(eg, vids=V(eg)$name[which(V(eg)$sex==focal_sex & V(eg)$name %in% relatives)])
    if(length(V(egsskin)$name) == 0){
      network_metrics[i, "ss_strength_kin"]<- NA
      network_metrics[i, "n_ss_kin"]<- 0
      network_metrics[i, "ss_strength_kin_prop"] <- NA
    } else {
      if(degree(egsskin, ego) == 0){
        network_metrics[i, "ss_strength_kin"]<- NA
        network_metrics[i, "n_ss_kin"]<- 0
        network_metrics[i, "ss_strength_kin_prop"] <- NA
      } else {
        network_metrics[i, "ss_strength_kin"]<-strength(egsskin, ego)
        network_metrics[i, "n_ss_kin"]<-degree(egsskin, ego)
        network_metrics[i, "ss_strength_kin_prop"] <- network_metrics[i, "ss_strength_kin"] / network_metrics[i, "ss_strength"]
      }
    }
    
    #opposite-sex strength, kin only:
    egopskin<-induced_subgraph(eg, vids=V(eg)$name %in% opposite_sex_names[opposite_sex_names %in% relatives])
    if(length(V(egopskin)$name) == 0){
      network_metrics[i, "os_strength_kin"]<- NA
      network_metrics[i, "n_os_kin"]<- 0
      network_metrics[i, "os_strength_kin_prop"] <- NA
    } else {
      if(degree(egopskin, ego) == 0){
        network_metrics[i, "os_strength_kin"]<- NA
        network_metrics[i, "n_os_kin"]<- 0
        network_metrics[i, "os_strength_kin_prop"] <- NA
      } else {
        network_metrics[i, "os_strength_kin"]<-strength(egopskin, ego)
        network_metrics[i, "n_os_kin"]<-degree(egopskin, ego)
        network_metrics[i, "os_strength_kin_prop"] <- network_metrics[i, "os_strength_kin"] / network_metrics[i, "os_strength"]
      }
    }

  }
}

write.csv(network_metrics, "real_network_metrics.csv", row.names = FALSE)

#Repeat for the results of the random model

library(foreach)
library(doParallel)

date_lookup<-sort(xydata3$Date)

starttime<-Sys.time()

cl<-makeCluster(detectCores()-1, outfile="../ids_completed.txt")
clusterEvalQ(cl, {library(SocGen); library(igraph)})
clusterExport(cl, c("kfinal", "availability_ego", "ai_mask", "date_lookup", "life_history_lookup"))
registerDoParallel(cl)

random_network_metrics<-foreach (n=1:nrow(availability_ego), .errorhandling='pass') %dopar% {
  
  ego<-as.character(availability_ego$dolphin_id[n])
  start<-availability_ego$entry[n]
  end<-availability_ego$depart[n]
  ego_sex<-life_history_lookup$sex[match(ego, life_history_lookup$dolphin_id)]
  ego_mask<-ai_mask
  
  ego_mask[ego,]<-1
  
  num_sim=length(kfinal)
  
  network_metrics<-data.frame(ego=rep(ego, num_sim), 
                              mixstrength=numeric(num_sim),
                              mixdegree=numeric(num_sim),
                              mixcc=numeric(num_sim),
                              ss_strength=numeric(num_sim),
                              ss_cc=numeric(num_sim),
                              ss_degree=numeric(num_sim),
                              os_strength=numeric(num_sim),
                              os_cc=numeric(num_sim),
                              os_degree=numeric(num_sim),
                              sex=rep(ego_sex, num_sim),
                              iteration=1:num_sim)
  
  for (j in 1:num_sim) {
    
    random1<-kfinal[[j]]
    ego_network<-random1[which(date_lookup>=start & date_lookup<=end),]
    ego_network$dates<-date_lookup[which(date_lookup>=start & date_lookup<=end)]
    eid<-unique(ego_network$group[which(ego_network$id==ego)])
    dolls<-unique(ego_network$id[which(ego_network$group %in% eid)])
    ego_network<-ego_network[which(ego_network$id %in% dolls),]
    
    network_ego<-simple_ratio(sightings=ego_network,
                              group_variable="group", 
                              dates="dates", 
                              IDs="id", 
                              symmetric=TRUE, 
                              mask=ego_mask)
    
    if(!is.matrix(network_ego)){next}
    
    network_ego[is.nan(network_ego)]<-0
    
    g<-graph.adjacency(network_ego, mode="undirected", weighted=TRUE, diag=FALSE)
    
    eg<-make_ego_graph(g, order=1, nodes=ego)[[1]] #is a list, return element 1
    
    if(length(E(eg))==0){next} else{
      
      #add sexes of individuals
      
      V(eg)$sex<-life_history_lookup$sex[match(V(eg)$name, life_history_lookup$dolphin_id)]
      
      network_metrics[j, "mixdegree"]<-degree(eg, ego)
      network_metrics[j, "mixstrength"]<-strength(eg, ego)
      network_metrics[j, "mixcc"]<-transitivity(eg, type="local", vids=ego)
      
      #pull out just same sex
      
      egss<-induced_subgraph(eg, vids=V(eg)$name[which(V(eg)$sex==ego_sex)])
      
      network_metrics[j, "ss_degree"]<-degree(egss, ego)
      network_metrics[j, "ss_strength"]<-strength(egss, ego)
      network_metrics[j, "ss_cc"]<-transitivity(egss, type="local", vids=ego)
      
      #pull out just opposite sex (and ego)
      
      opposite_sex_names <- c(ego, V(eg)$name[which(V(eg)$sex!=ego_sex)])
      egops<-induced_subgraph(eg, vids=V(eg)$name %in% opposite_sex_names)
      
      network_metrics[j, "os_degree"]<-degree(egops, ego)
      network_metrics[j, "os_strength"]<-strength(egops, ego)
      network_metrics[j, "os_cc"]<-transitivity(egops, type="local", vids=ego)
      # ego_sims[[j]]<-network_ego
    }
    
  }
  cat(paste0(n, " networks complete for ", ego, "\n"))
  network_metrics
}

stopCluster(cl)
endtime<-Sys.time()

endtime-starttime #check run time (4.4 hours on last run)

all_random_metrics<-do.call("rbind", random_network_metrics)

write.csv(all_random_metrics, "all_random_metrics.csv", row.names = FALSE)

####See Figure Plotting for figures and aggregating results



