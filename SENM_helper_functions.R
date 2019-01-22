#SENM helper functions
library(compiler)

#This feature adds the launch point and extra points so a sampling area can always be calculated

add_launch<-function(x) {
  launch[1]<-x[1,1]
  x<-rbind(x, launch)
  return(x)
}

#fast random points
# fast_random_points<-function(probweights, nd, dc, dgdf) {
#   daily_dolphins<-sample(names(probweights),size=nd,replace=FALSE,prob=probweights)
#   df<-dgdf[,which(colnames(dgdf) %in% daily_dolphins)]
#   df[is.na(df)]<-0
#   cells<-apply(df, 2, function(x) sample(seq_len(nrow(df)),1,prob=x))
#   res<-data.frame(dc[cells,], names(cells))
# }

# cmp_fast_random_points<-cmpfun(fast_random_points)


#This is a speed optimized version to calculate the simple ratio index with a one-day sampling period
#Takes a dataframe of start and end times for each dyad and a dataframe of sightings


ai_filtered<-function(sightings=sightings, 
                      group_variable=group_variable, 
                      dates=dates, 
                      IDs=IDs,
                      lookup_run_ego=lookup_run_ego){
  z<-sightings
  focals<-as.character(unique(c(lookup_run_ego[,1], lookup_run_ego[,2])))
  nf<-length(focals)
  matnames<-list(focals,focals)
  currentmat<-matrix(c(rep(NA,nf^2)),nrow=nf,dimnames=matnames)
  
  for (i in 1:nrow(lookup_run_ego)){
    
    ego<-as.character(lookup_run_ego[i,1])
    alter<-as.character(lookup_run_ego[i,2])
    
    all_ego<-z[which(z[,IDs]==ego),]
    all_alter<-z[which(z[,IDs]==alter),]
    
    all_ego<-all_ego[(all_ego[,dates]>=lookup_run_ego[i,7]),] 
    all_ego<-all_ego[(all_ego[,dates]<=lookup_run_ego[i,8]),]
    
    all_alter<-all_alter[(all_alter[,dates]>=lookup_run_ego[i,7]),] 
    all_alter<-all_alter[(all_alter[,dates]<=lookup_run_ego[i,8]),]
    
    sample<-all_ego[(all_ego[,group_variable] %in% all_alter[,group_variable]),]
    
    #Numerator
    X<-length(unique(sample[,dates]))
    #Ego without alter
    Ya<-length(setdiff(all_ego[,dates], all_alter[,dates]))
    #Alter without ego
    Yb<-length(setdiff(all_alter[,dates], all_ego[,dates]))
    #Both seen but not together
    Yab<-length(intersect(all_ego[,dates], all_alter[,dates]))-X
    
    currentmat[ego,alter]<-(X/(X+Ya+Yb+Yab))
  }
  return(currentmat)
}

cmp_ai_filtered<-cmpfun(ai_filtered)
