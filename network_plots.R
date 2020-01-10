#Plot examples of female and male networks

library(igraph)

options(stringsAsFactors = FALSE)

load("real_ai_egos.RData")

lh<-read.delim("Raw_input_files/LifeHistory_20180507.txt")
lh<-lh[!duplicated(lh$dolphin_id),]

male<-"AGA" #or try STA
female<-"SHO" #or try CAB KIY SHO SQL

male_ego<-real_ai_egos[[male]] 
female_ego<-real_ai_egos[[female]] 

male_ego[is.na(male_ego)]<-0
female_ego[is.na(female_ego)]<-0

g<-graph.adjacency(male_ego, mode="undirected", weighted=TRUE,diag=FALSE)
gf<-graph.adjacency(female_ego, mode="undirected", weighted=TRUE,diag=FALSE)

g<-make_ego_graph(g, order=1, nodes=male)[[1]] #is a list, return element 1
gf<-make_ego_graph(gf, order=1, nodes=female)[[1]] 

V(g)$sex<-lh$sex[match(V(g)$name,lh$dolphin_id)]
V(g)$color<-ifelse(V(g)$sex=="MALE", "blue3", "red3")
V(g)$color[which(V(g)$name==male)]<-"black"
V(g)$age<-substr(lh$birth_date,1,4)[match(V(g)$name, lh$dolphin_id)]
V(g)$size<-max(as.numeric(V(g)$age))-as.numeric(V(g)$age)
V(g)$shape<-ifelse(V(g)$sex=="MALE", "square", "circle")

set.seed(7)  #to get same network layout

minC <- rep(-Inf, vcount(g))
maxC <- rep(Inf, vcount(g))
minC[which(V(g)$name==male)]<-maxC[which(V(g)$name==male)]<-0
male_layout <- layout_with_fr(g, minx=minC, maxx=maxC,
                     miny=minC, maxy=maxC)

V(gf)$sex<-lh$sex[match(V(gf)$name,lh$dolphin_id)]
V(gf)$color<-ifelse(V(gf)$sex=="MALE", "blue3", "red3")
V(gf)$color[which(V(gf)$name==female)]<-"black"
V(gf)$age<-substr(lh$birth_date,1,4)[match(V(gf)$name, lh$dolphin_id)]
V(gf)$size<-max(as.numeric(V(gf)$age))-as.numeric(V(gf)$age)
V(gf)$shape<-ifelse(V(gf)$sex=="MALE", "square", "circle")

minC <- rep(-Inf, vcount(gf))
maxC <- rep(Inf, vcount(gf))
minC[which(V(gf)$name==female)]<-maxC[which(V(gf)$name==female)]<-0
female_layout <- layout_with_fr(gf, minx=minC, maxx=maxC,
                              miny=minC, maxy=maxC)

#clustering coefficients
gall<-graph.adjacency(male_ego, mode="undirected", weighted=TRUE,diag=FALSE)
gfall<-graph.adjacency(female_ego, mode="undirected", weighted=TRUE,diag=FALSE)


gsize<-transitivity(gall, type="local", vids=V(gall)[V(gall)$name %in% V(g)$name])*100
gsize[is.na(gsize)]<-1

gfsize<-transitivity(gfall, type="local", vids=V(gfall)[V(gfall)$name %in% V(gf)$name])*100


#windows()
pdf(file="networks20200109_2.pdf", width=7, height=4.33)

par(mar=c(0,0,0,0), mfrow=c(1,2))

plot(g,
     # vertex.size=((degree(g)/max(degree(g)))+0.25)*100,
     vertex.size=gsize,
     vertex.color=adjustcolor(V(g)$color,alpha.f=1),
     vertex.shapes=V(g)$shape,
     edge.width = edge_attr(g)$weight*12,
     edge.curved = rep(-.2,length(edge_attr(g)$weight)),
     edge.arrow.size = 0, 
     vertex.label = NA,
     edge.color = adjustcolor("black", alpha.f = 0.4), 
     layout=male_layout, 
     rescale=FALSE,
     xlim=c(-9,9), 
     ylim=c(-10,10)
     )
text(x=-6, y=8, "a", cex=1.5)

plot(gf,
     # vertex.size=((degree(gf)/max(degree(gf)))+0.25)*100, 
     vertex.size=gfsize,
     vertex.color=adjustcolor(V(gf)$color,alpha.f=1),
     vertex.shapes=V(gf)$shape,
     edge.width = edge_attr(gf)$weight*12,
     edge.curved = rep(-.2,length(edge_attr(gf)$weight)),
     edge.arrow.size = 0, 
     vertex.label = NA,
     edge.color = adjustcolor("black", alpha.f = 0.4), 
     layout=female_layout, 
     rescale=FALSE,
     xlim=c(-9,9), 
     ylim=c(-10,10)
     )
text(x=-6, y=8, "b", cex=1.5)

dev.off()
