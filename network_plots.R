#Plot examples of female and male networks

library(igraph)

options(stringsAsFactors = FALSE)

load("../ai_egos_for_plotting.RData")

lh<-read.delim("../Raw_input_files/LifeHistory_20180507.txt")
lh<-lh[!duplicated(lh$dolphin_id),]

male<-"AGA" #or try STA
female<-"SQL" #or try CAB KIY SHO

male_ego<-ai_egos[[male]] 
female_ego<-ai_egos[[female]] 

male_ego[is.na(male_ego)]<-0
female_ego[is.na(female_ego)]<-0

g<-graph.adjacency(male_ego, mode="undirected", weighted=TRUE,diag=FALSE)
gf<-graph.adjacency(female_ego, mode="undirected", weighted=TRUE,diag=FALSE)

g<-make_ego_graph(g, order=1, nodes=male)[[1]] #is a list, return element 1
gf<-make_ego_graph(gf, order=1, nodes=female)[[1]] 

V(g)$sex<-lh$sex[match(V(g)$name,lh$dolphin_id)]
V(g)$color<-ifelse(V(g)$sex=="MALE", "blue", "red")
V(g)$color[which(V(g)$name==male)]<-"black"
V(g)$age<-substr(lh$birth_date,1,4)[match(V(g)$name, lh$dolphin_id)]
V(g)$size<-max(as.numeric(V(g)$age))-as.numeric(V(g)$age)

V(gf)$sex<-lh$sex[match(V(gf)$name,lh$dolphin_id)]
V(gf)$color<-ifelse(V(gf)$sex=="MALE", "blue", "red")
V(gf)$color[which(V(gf)$name==female)]<-"black"
V(gf)$age<-substr(lh$birth_date,1,4)[match(V(gf)$name, lh$dolphin_id)]
V(gf)$size<-max(as.numeric(V(gf)$age))-as.numeric(V(gf)$age)

set.seed(3) #not sure if this controls fr algorithm but give it a whirl

# windows()
pdf(file="../networks.pdf")
layout(matrix(c(rep(1, 6), rep(2, 4)), ncol=2, byrow = TRUE))
par(mar=c(1,1,1,1))
plot(g,
     vertex.size=degree(g)+3, 
     vertex.color=adjustcolor(V(g)$color,alpha.f=0.7),
     edge.width = edge_attr(g)$weight*10,
     edge.curved = rep(-.2,length(edge_attr(g)$weight)),
     edge.arrow.size = 0, 
     vertex.label = NA,
     edge.color = adjustcolor("black", alpha.f = 0.4), 
     layout=layout_with_fr(g)
     )


plot(gf,
     vertex.size=degree(gf)+3, 
     vertex.color=adjustcolor(V(gf)$color,alpha.f=0.7),
     edge.width = edge_attr(gf)$weight*10,
     edge.curved = rep(-.2,length(edge_attr(gf)$weight)),
     edge.arrow.size = 0, 
     vertex.label = NA,
     edge.color = adjustcolor("black", alpha.f = 0.4), 
     layout=layout_with_fr(gf)
     )
dev.off()


