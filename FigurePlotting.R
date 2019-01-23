#SENM Galezo Figure Plotting

options(stringsAsFactors = FALSE)

#Real data
real_network_metrics<-read.csv("real_network_metrics.csv")
real_network_metrics$ss_cc[is.na(real_network_metrics$ss_cc)]<-0

rmixdegree<-aggregate(mixdegree~sex,data=real_network_metrics, mean)
rmixstrength<-aggregate(mixstrength~sex,data=real_network_metrics, mean)
rmixcc<-aggregate(mixcc~sex,data=real_network_metrics, mean)
rss_strength<-aggregate(ss_strength~sex,data=real_network_metrics, mean)
rss_degree<-aggregate(ss_degree~sex,data=real_network_metrics, mean)
rss_cc<-aggregate(ss_cc~sex,data=real_network_metrics, mean)

#Simulated data
all_random_metrics<-read.csv("all_random_metrics.csv")
all_random_metrics$ss_cc[is.nan(all_random_metrics$ss_cc)]<-0

randmixdegree<-aggregate(mixdegree~sex+iteration,data=all_random_metrics, mean)
randmixstrength<-aggregate(mixstrength~sex+iteration, data=all_random_metrics, mean)
randmixcc<-aggregate(mixcc~sex+iteration,data=all_random_metrics, mean)
randss_degree<-aggregate(ss_degree~sex+iteration,data=all_random_metrics, mean)
randss_strength<-aggregate(ss_strength~sex+iteration,data=all_random_metrics, mean)
randss_cc<-aggregate(ss_cc~sex+iteration,data=all_random_metrics, mean)

#Plot real data vs simulated data

####Degreee
{
windows()
  #pdf(file="rand_degree.pdf")
par(mfrow=c(2,2), mar=c(4.1, 4.1, 3.1, 1.1))
hist(randmixdegree$mixdegree[which(randmixdegree$sex=="FEMALE")], xlim=c(15,30), main="Female Degree", xlab="", col="grey",
     border="darkgrey")  

segments(x0=rmixdegree$mixdegree[which(rmixdegree$sex=="FEMALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
arrows(x0=rmixdegree$mixdegree[which(rmixdegree$sex=="FEMALE")], y0=1, y1=0, length=0.12,
       lwd=2, col="darkblue")

hist(randmixdegree$mixdegree[which(randmixdegree$sex=="MALE")], xlim=c(15,30), main="Male Degree", xlab="", col="grey",
     border="darkgrey", ylab="")  

segments(x0=rmixdegree$mixdegree[which(rmixdegree$sex=="MALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
arrows(x0=rmixdegree$mixdegree[which(rmixdegree$sex=="MALE")], y0=1, y1=0, length=0.12,
       lwd=2, col="darkblue")


hist(randss_degree$ss_degree[which(randss_degree$sex=="FEMALE")], xlim=c(5,16), main="Female Same Sex Degree", xlab="Number of Associates",
     col="grey",
     border="darkgrey")  

segments(x0=rss_degree$ss_degree[which(rss_degree$sex=="FEMALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
arrows(x0=rss_degree$ss_degree[which(rss_degree$sex=="FEMALE")], y0=1, y1=0, length=0.12,
       lwd=2, col="darkblue")

hist(randss_degree$ss_degree[which(randss_degree$sex=="MALE")], xlim=c(5,16), main="Male Same Sex Degree", xlab="Number of Associates",
     col="grey",
     border="darkgrey", ylab="")  
segments(x0=rss_degree$ss_degree[which(rss_degree$sex=="MALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
arrows(x0=rss_degree$ss_degree[which(rss_degree$sex=="MALE")], y0=1, y1=0, length=0.12,
       lwd=2, col="darkblue")


dev.off()
}

####Strength
{
windows()
  #pdf(file="rand_strength.pdf")
par(mfrow=c(2,2), mar=c(4.1, 4.1, 3.1, 1.1))
hist(randmixstrength$mixstrength[which(randmixstrength$sex=="FEMALE")], xlim=c(0.5,1.6), main="Female Strength", xlab="", col="grey",
     border="darkgrey")  
segments(x0=rmixstrength$mixstrength[which(rmixstrength$sex=="FEMALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
arrows(x0=rmixstrength$mixstrength[which(rmixstrength$sex=="FEMALE")], y0=1, y1=0, length=0.12,
       lwd=2, col="darkblue")


hist(randmixstrength$mixstrength[which(randmixstrength$sex=="MALE")], xlim=c(0.5,1.6), main="Male Strength", xlab="", col="grey",
     border="darkgrey", ylab="")  
segments(x0=rmixstrength$mixstrength[which(rmixstrength$sex=="MALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
arrows(x0=rmixstrength$mixstrength[which(rmixstrength$sex=="MALE")], y0=1, y1=0, length=0.12,
       lwd=2, col="darkblue")


hist(randss_strength$ss_strength[which(randss_strength$sex=="FEMALE")], xlim=c(0.2,1), main="Female Same Sex Strength", xlab="Cumulative Weighted Degree",
     col="grey",
     border="darkgrey")  
segments(x0=rss_strength$ss_strength[which(rss_strength$sex=="FEMALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
arrows(x0=rss_strength$ss_strength[which(rss_strength$sex=="FEMALE")], y0=1, y1=0, length=0.12,
       lwd=2, col="darkblue")


hist(randss_strength$ss_strength[which(randss_strength$sex=="MALE")], xlim=c(0.2,1), main="Male Same Sex Strength", xlab="Cumulative Weighted Degree",
     col="grey",
     border="darkgrey", ylab="")  
segments(x0=rss_strength$ss_strength[which(rss_strength$sex=="MALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
arrows(x0=rss_strength$ss_strength[which(rss_strength$sex=="MALE")], y0=1, y1=0, length=0.12,
       lwd=2, col="darkblue")



dev.off()
}

####Clustering
{
  windows()
  #pdf(file="rand_clustering.pdf")
  par(mfrow=c(2,2), mar=c(4.1, 4.1, 3.1, 1.1))
  hist(randmixcc$mixcc[which(randmixcc$sex=="FEMALE")], xlim=c(0.2,0.6), main="Female Clustering", xlab="", col="grey",
       border="darkgrey")  
  segments(x0=rmixcc$mixcc[which(rmixcc$sex=="FEMALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
  arrows(x0=rmixcc$mixcc[which(rmixcc$sex=="FEMALE")], y0=1, y1=0, length=0.12,
         lwd=2, col="darkblue")
  
  hist(randmixcc$mixcc[which(randmixcc$sex=="MALE")], xlim=c(0.2,0.6), main="Male Clustering", xlab="", col="grey",
       border="darkgrey", ylab="")  
  segments(x0=rmixcc$mixcc[which(rmixcc$sex=="MALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
  arrows(x0=rmixcc$mixcc[which(rmixcc$sex=="MALE")], y0=1, y1=0, length=0.12,
         lwd=2, col="darkblue")
  
  hist(randss_cc$ss_cc[which(randss_cc$sex=="FEMALE")], xlim=c(0.2,0.6), main="Female Same Sex Clustering", xlab="Clustering Coefficient",
       col="grey",
       border="darkgrey")  
  segments(x0=rss_cc$ss_cc[which(rss_cc$sex=="FEMALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
  arrows(x0=rss_cc$ss_cc[which(rss_cc$sex=="FEMALE")], y0=1, y1=0, length=0.12,
         lwd=2, col="darkblue")
  
  hist(randss_cc$ss_cc[which(randss_cc$sex=="MALE")], xlim=c(0.2,0.6), main="Male Same Sex Clustering", xlab="Clustering Coefficient",
       col="grey",
       border="darkgrey", ylab="")  
  segments(x0=rss_cc$ss_cc[which(rss_cc$sex=="MALE")], y0=0, y1=par("usr")[4], lty=2, lwd=2, col="darkblue")
  arrows(x0=rss_cc$ss_cc[which(rss_cc$sex=="MALE")], y0=1, y1=0, length=0.12,
         lwd=2, col="darkblue")
  
  dev.off()
  }

##Male and Female colors

mcol <- rgb(54, 100, 139, max = 255, alpha = 175, names = "blue50")
fcol <- rgb(0, 200, 55, max = 255, alpha = 100, names = "green50")

#Plot real data for individual males and females

####Male vs Female Degree
{
  #windows()
  pdf(file="Degree.pdf", width=7, height=4)
  par(mfrow=c(1,2), mai=c(1.02, 0.82, 0.82, 0.22))
  hist(real_network_metrics$mixdegree[which(real_network_metrics$sex=="MALE")], 
       col=mcol, 
       xlim=c(0,45),ylim=c(0,22), 
       main="Mixed Sex Degree", 
       xlab="Number of Associates", 
       lty="blank")
  hist(real_network_metrics$mixdegree[which(real_network_metrics$sex=="FEMALE")], 
       col=fcol,
       lty="blank",
       add=TRUE)
  abline(v=rmixdegree$mixdegree[which(rmixdegree$sex=="FEMALE")], 
         lty=2,
         lwd=2, 
         col="darkgreen")
  abline(v=rmixdegree$mixdegree[which(rmixdegree$sex=="MALE")], 
         lty=3, 
         lwd=2, 
         col="darkblue")
  legend(x=28,
         y=22,
         pch=15,
         pt.cex=1.25,
         col=c(mcol, fcol), 
         bty="n",
         legend=c(" Male", " Female"))
  legend(x=25.5,
         y=22,
         lty=c(3,2),
         col=c("darkblue", "darkgreen"),
         legend=c(NA, NA),
         bty="n",
         lwd=2)
  
#Male vs Female Same Sex Degree
  hist(real_network_metrics$ss_degree[which(real_network_metrics$sex=="MALE")], 
       col=mcol, 
       xlim=c(0,22),
       ylim=c(0,15), 
       main="Same Sex Degree", 
       xlab="Number of Same Sex Associates", 
       nclass=8,
       lty="blank",
       ylab=NULL)
  hist(real_network_metrics$ss_degree[which(real_network_metrics$sex=="FEMALE")], 
       lty="blank", 
       col=fcol, 
       add=TRUE)
  abline(v=rss_degree$ss_degree[which(rss_degree$sex=="FEMALE")], 
         lty=2, 
         lwd=2, 
         col="darkgreen")
  abline(v=rss_degree$ss_degree[which(rss_degree$sex=="MALE")], 
         lty=3, 
         lwd=2, 
         col="darkblue")
  dev.off()
}

####Male vs Female Strength
{
  pdf(file="Strength.pdf", width=7, height=4)
  par(mfrow=c(1,2), mai=c(1.02, 0.42, 0.82, 0.42))
  hist(real_network_metrics$mixstrength[which(real_network_metrics$sex=="MALE")], 
       col=mcol, 
       xlim=c(0,3.5),
       ylim=c(0,30), 
       main="Strength", 
       xlab="Cumulative Weighted Degree (SRI)", 
       lty="blank") 
  hist(real_network_metrics$mixstrength[which(real_network_metrics$sex=="FEMALE")], 
       lty="blank", 
       col=fcol, 
       add=TRUE)
  abline(v=rmixstrength$mixstrength[which(rmixstrength$sex=="FEMALE")], 
         lty=2,
         lwd=2, 
         col="darkgreen")
  abline(v=rmixstrength$mixstrength[which(rmixstrength$sex=="MALE")], 
         lty=3, 
         lwd=2, 
         col="darkblue")
  legend(x=2,
         y=30,
         pch=15,
         pt.cex=1.25,
         col=c(mcol, fcol), 
         bty="n",
         legend=c("  Male", "  Female"))
  legend(x=1.85,
         y=30,
         lty=c(3,2),
         col=c("darkblue", "darkgreen"),
         legend=c(NA, NA),
         bty="n",
         lwd=2)


##Male vs Female Same Sex Strength

  hist(real_network_metrics$ss_strength[which(real_network_metrics$sex=="MALE")], 
       col=mcol, 
       xlim=c(0,2.5),
       ylim=c(0,40), 
       main="Same Sex Strength", 
       xlab="Same Sex Cumulative Weighted Degree (SRI)", 
       lty="blank")
  hist(real_network_metrics$ss_strength[which(real_network_metrics$sex=="FEMALE")], 
       lty="blank", 
       col=fcol,
       nclass=4,
       add=TRUE)
  abline(v=rss_strength$ss_strength[which(rss_strength$sex=="FEMALE")], 
         lty=2,
         lwd=2, 
         col="darkgreen")
  abline(v=rss_strength$ss_strength[which(rss_strength$sex=="MALE")], 
         lty=3, 
         lwd=2, 
         col="darkblue")
  # legend("topright", 
  #        fill=c(mcol, fcol), 
  #        lty=c(3,2),
  #        legend=c("Male", "Female"))
  dev.off()
  }

####Male vs Female Clustering
{
  pdf(file="Clustering.pdf", width=7, height=4)
  par(mfrow=c(1,2), mai=c(1.02, 0.42, 0.82, 0.42))
  hist(real_network_metrics$mixcc[which(real_network_metrics$sex=="MALE")], 
       col=mcol, 
       xlim=c(0,1.1),
       ylim=c(0,26), 
       main="Clustering Coefficient", 
       xlab="Local Transitivity", 
       lty="blank") 
  hist(real_network_metrics$mixcc[which(real_network_metrics$sex=="FEMALE")], 
       lty="blank", 
       col=fcol, 
       add=TRUE)
  abline(v=rmixcc$mixcc[which(rmixcc$sex=="FEMALE")], 
         lty=2,
         lwd=2, 
         col="darkgreen")
  abline(v=rmixcc$mixcc[which(rmixcc$sex=="MALE")], 
         lty=3, 
         lwd=2, 
         col="darkblue")
  legend(x=0.65,
         y=25,
         pch=15,
         pt.cex=1.25,
         col=c(mcol, fcol), 
         bty="n",
         legend=c("  Male", "  Female"))
  legend(x=0.60,
         y=25,
         lty=c(3,2),
         col=c("darkblue", "darkgreen"),
         legend=c(NA, NA),
         bty="n",
         lwd=2)

  hist(real_network_metrics$ss_cc[which(real_network_metrics$sex=="MALE")], 
       col=mcol, 
       xlim=c(0,1.1),
       ylim=c(0,26), 
       main="Same Sex Clustering Coefficient", 
       xlab="Same Sex Local Transitivity", 
       lty="blank") 
  hist(real_network_metrics$ss_cc[which(real_network_metrics$sex=="FEMALE")], 
       lty="blank", 
       col=fcol, 
       add=TRUE)
  abline(v=rss_cc$ss_cc[which(rss_cc$sex=="FEMALE")], 
         lty=2,
         lwd=2, 
         col="darkgreen")
  abline(v=rss_cc$ss_cc[which(rss_cc$sex=="MALE")], 
         lty=3, 
         lwd=2, 
         col="darkblue")
  dev.off()
}


#Add expected values to observed data

real_network_metrics$expected_mixstrength<-aggregate(mixstrength~ego, data=all_random_metrics, mean)[,2]

real_network_metrics$expected_mixdegree<-aggregate(mixdegree~ego, data=all_random_metrics, mean)[,2]

real_network_metrics$expected_mixcc<-aggregate(mixcc~ego, data=all_random_metrics, mean)[,2]

real_network_metrics$expected_ss_strength<-aggregate(ss_strength~ego, data=all_random_metrics, mean)[,2]

real_network_metrics$expected_ss_cc<-aggregate(ss_cc~ego, data=all_random_metrics, mean)[,2]

real_network_metrics$expected_ss_degree<-aggregate(ss_degree~ego, data=all_random_metrics, mean)[,2]

expected_percent_close_kin<-aggregate(percent_close_kin~ego, data=all_random_metrics, mean)

real_network_metrics$expected_percent_close_kin<-expected_percent_close_kin$percent_close_kin[match(real_network_metrics$ego, expected_percent_close_kin$ego)]

write.csv(real_network_metrics, "observed_expected_network_metrics.csv")










