rm(list = ls())
library(igraph)
library(igraphdata)
library(dplyr)
library(ggplot2)
library(visNetwork)
library(hrbrthemes)
library(viridis)
set.seed(42)

data <- read.csv("international_currencies_data.txt", header=TRUE, sep="\t",
                 skip=9)
edges <- data[,c(1:2,4,6:8)] %>% filter(quote1900==1) %>%
  select(-quote1900)
vertices <- read.csv("vertices.csv", header=TRUE, sep=";")

#Визуализация
g <- graph_from_data_frame(d=edges, vertices=vertices, directed=T)
plot(g)

cols <- c("lightcoral", "lightsalmon3", "green", "yellow", "dodgerblue", "light blue")
z <- rep(NA, length(V(g)$region))
for (reg in unique(V(g)$region)){
  z[which(V(g)$region==reg)] <- cols[which(unique(V(g)$region)==reg)]
}
V(g)$color <- z

plot(g, layout = layout_with_kk, vertex.size=V(g)$rgdp*0.75,
     edge.arrow.size=0.5, edge.arrow.width=0.75, vertex.color=V(g)$color,
     edge.color=ifelse(E(g)$colony==1,"gold","grey"))

# Степени вершин
degree_in <- degree(g, mode="in")
degree_out <- degree(g, mode="out")
p <- edge_density(g)

tmp <- data.frame(degree_in, degree_out, row.names=names(degree_in))
tmp <- tmp[order(tmp$degree_in, decreasing=T),]
barplot(t(as.matrix(tmp)),beside=T,las=2, ylim=c(0,45), col=c("dodger blue","lightcoral"),
        main="Центральность по степени")
legend("right",rownames(t(as.matrix(tmp))),cex =0.8,fill=c("dodger blue","lightcoral"),title="Degree mode")

par(mfrow=c(1,2))
g.random <- erdos.renyi.game(n=gorder(g), p.or.m=p, type="gnp", directed=TRUE)
plot(g.random, edge.arrow.size=0.5, edge.arrow.width=0.75, vertex.size=5, 
     vertex.label=NA)
g.BA <- sample_pa(n=gorder(g), power=1, directed=TRUE)
plot(g.BA, edge.arrow.size=0.5, edge.arrow.width=0.75, vertex.size=5, 
     vertex.label=NA)

#################################################
# Сравнение метрик графов
#################################################

c1 <- rgb(30,144,255, max=255, alpha=80)
c2 <- rgb(240,128,128, max=255, alpha=80)
par(mfrow=c(1,2))
for (gr in list(g.random, g.BA)){
  h1 <- hist(degree(gr, mode="in"), breaks = 0:max(degree(gr, mode="in")), prob = TRUE, xlab="Degree", col=c1, alpha = 80,
             main="")
  h2 <- hist(degree(gr, mode="out"), breaks = 0:max(degree(gr, mode="out")), prob = TRUE, xlab="Degree", col=c2, alpha = 80, add=TRUE,
             main="")
  legend("right",c("in","out"),cex =0.8,fill=c(c1,c2),title="Degree mode")
}

graphs <- list(g, g.random, g.BA)
ed <- unlist(lapply(graphs, edge_density))
df <- data.frame(edge_density=ed, row.names=c("Исходный граф","Случайный граф",
                         "Барабаши-Альберт"))
df$diameter <- unlist(lapply(graphs, diameter))
df$mean_distance <- unlist(lapply(graphs, mean_distance))
df$transitivity <- unlist(lapply(graphs,
                                 function(x){transitivity(x,type="global")}))
df$reciprocity <- unlist(lapply(graphs, reciprocity))
df$cohesion <- unlist(lapply(graphs, cohesion))
denom1 <- gorder(g)*(gorder(g)-1)/2
df$dyad_mut <- unlist(lapply(graphs, 
                             function(x){dyad_census(x)$mut}))/denom1
df$dyad_asym <- unlist(lapply(graphs, 
                             function(x){dyad_census(x)$asym}))/denom1
df$dyad_null <- unlist(lapply(graphs, 
                             function(x){dyad_census(x)$null}))/denom1
df$triangles <- unlist(lapply(graphs, 
                              function(x){triad_census(x)[16]}))
df <- as.data.frame(round(t(df),3))
# write.csv(df,"models_metrics.csv",row.names=T)

#################################################
# Центральности
#################################################

par(mfrow=c(1,1))
clos <- sort(closeness(g, mode="in"), decreasing=T)
barplot(clos,las=2, col="dodger blue", main="Центральность по близости")

betw <- sort(betweenness(g, normalized=T), decreasing=T)
barplot(betw,las=2, col="dodger blue", main="Центральность по кратчайшему пути")

eig_dir <- eigen_centrality(g, directed=T)$vector
eig_undir <- eigen_centrality(g, directed=F)$vector
tmp <- data.frame(eig_undir, eig_dir, row.names=names(eig_undir))
tmp <- tmp[order(tmp$eig_undir, tmp$eig_dir, decreasing=T),]
barplot(t(as.matrix(tmp)),beside=T,las=2, col=c("dodger blue","lightcoral"),
        main="Центральность по собственному значению")
legend("right",c("undirected","directed"),cex =0.8,
       fill=c("dodger blue","lightcoral"),title="Degree mode")

#################################################
# Предпочтительное присоединение
#################################################

assortativity_degree(g)
assortativity_nominal(g, types=as.numeric(as.factor(V(g)$region)))

edgeList <- as_edgelist(g)
region <- matrix(c(V(g)$name, V(g)$region), nrow = gorder(g), ncol = 2)
FromLabel <- region[match(edgeList[,1],region[,1]),2]
ToLabel <- region[match(edgeList[,2],region[,1]),2]

tmp <- c("Dyadicity","Heterophilicity","number of vertices")
for (reg in unique(V(g)$region)){
  same_label <- sum((FromLabel==reg) & (ToLabel==reg))
  diff_label1 <- sum((FromLabel==reg) & (ToLabel!=reg))
  diff_label2 <- sum((FromLabel!=reg) & (ToLabel==reg))
  diff_label <- diff_label1 + diff_label2
  
  nc <- sum(V(g)$region==reg)
  nnc <- sum(V(g)$region!=reg)
  
  exp_dyad <- nc*(nc-1)*p
  exp_hetero <- nnc*nc*p*2
  
  tmp <- cbind(tmp,c(same_label/exp_dyad,diff_label/exp_hetero,nc))
}
tmp <- as.data.frame(tmp[,-1])
tmp <- round(sapply(tmp, as.numeric),3)
colnames(tmp) <- unique(V(g)$region)
rownames(tmp) <- c("Dyadicity","Heterophilicity","number of vertices")
tmp <- t(tmp)[order(t(tmp)[,1], decreasing=T),]
tmp
# write.csv(tmp, "dyad_hetero.csv", row.names=T)

#################################################
# Кластеризация
#################################################

V(g)$n_reg <- as.numeric(as.factor(V(g)$region))
g_undir <- as.undirected(g)

par(mfrow=c(1,2))
colrs <- c("gray50", "tomato", "gold", "yellowgreen", "dodger blue",
           "lightsalmon", "white")

ceb <- cluster_edge_betweenness(g_undir)
plot(ceb, g_undir, main="Edge betweenness")
sizes(ceb)
modularity(ceb)

cfg <- cluster_fast_greedy(g_undir)
plot(cfg, g_undir, main="Fast greedy")
membership(cfg) 
sizes(cfg)
modularity(cfg)

cl <- cluster_louvain(g_undir)
plot(cl, g_undir, main="Louvain")
membership(cl) 
sizes(cl)
modularity(cl)

A <- as_adjacency_matrix(g_undir, names = TRUE)
A <- as.matrix(A)
S <- cor(A)
D <- 1-S
d <- as.dist(D)
cc <- hclust(d)
# plot(cc)
cls <- cutree(cc, k = 7)
cor(cls, V(g_undir)$n_reg)
plot(g_undir, vertex.color = colrs[cls], main="Hierarchical")

compare(ceb,cl)
compare(ceb,cfg)
compare(cfg,cl)

#################################################
# Модерирование
#################################################

library(statnet)
library(ergm)

dd <- data %>% select_if(~!any(is.na(.)))
dd$country_A
logit <- glm(quote1900 ~ log(rgdp) + log(rgdp_B) + log(dist) + log(1+bitrade) +
               log(1+coverage) + as.factor(colony) + as.factor(gold) + as.factor(gold_B),
             family = 'binomial', data = dd)
summary(logit)

ergm_g <- as.network(edgeList, vertices=vertices)
v_id <- network.vertex.names(ergm_g)
ergm_g%v%"rgdp" <- vertices[match(v_id, vertices[,1]),]$rgdp
ergm_g%v%"region" <- vertices[match(v_id, vertices[,1]),]$region
list.vertex.attributes(ergm_g)

fit <- ergm(ergm_g ~ edges + mutual + nodematch("region") +
              absdiff("rgdp") + nodefactor("region")) 
summary(fit)




