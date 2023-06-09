library("raster")
library("dplyr")
library("stargazer") 
library("ggplot2")
library("openintro")
library("gdata")
library("ggmap")
library("here")
library("readr")
library("geosphere")
library("stargazer")
library("sandwich")
library("osmdata")

set.seed(42)
register_google(key = "AIzaSyBf14UoBj_BiZITyKTGpEIne_d0swiSfBo")

data <- read.csv(here("MosRentUrbanMasteR.csv"))
MoscowCenterMap <- get_map(location = c("Moscow"), zoom = 10, maptype="roadmap")
center <- geocode(c("Kremlin"))

# Длина МКАД - 108.9 км -> радиус примерно 17.33
radius <- 108.9 / 2 / pi
data$DistToCenter <- as.numeric(distm(data[21:22], center)/1000)
inner <- data[data$DistToCenter < radius,]

# Пункт 1.
ggmap(MoscowCenterMap) +
  geom_point(aes(x=longitude, y=latitude, color=factor(ethnic_discrim)), size=1.5, data=inner, alpha=0.7) +
  scale_color_manual("Ethnic discrimination", values=c("green3", "red"))
ggsave(here("UrbanMap1_EfimovAlexandr.pdf"))

# Пункт 3.
data$inside <- as.numeric(data$DistToCenter < radius)
ggmap(MoscowCenterMap) +
  geom_point(aes(x=longitude, y=latitude, color=factor(inside)), size=1.5, data, alpha=0.7) +
  scale_color_manual("House is inside MKAD", values=c("red", "blue"))
ggsave(here("UrbanMap2_EfimovAlexandr.pdf"))

# km.res <- kmeans(inner[,c("longitude","latitude")], 100, nstart=50)
# inner$cluster <- km.res$cluster
# ggmap(MoscowCenterMap) +
#   geom_point(aes(x=longitude, y=latitude, color=factor(cluster)), size=1.5, inner, alpha=0.8)

# Пункт 4.
q <- opq(bbox = "Moscow") %>%
  add_osm_feature(key = "amenity", value = "university")
educ <- osmdata_sf(q)$osm_polygons
educ <- educ[!(is.na(educ$name)),]
educ <- educ[educ$amenity=="university",]

educ$centroid_lon <- NA
educ$centroid_lat <- NA
for(i in 1:nrow(educ)){
  input <- unname(unlist(educ[i,"geometry"]))
  if (!(is.null(input))){
    centr <- centroid(matrix(input, ncol=2, byrow=FALSE))
    educ[i,"centroid_lon"] <- centr[1]
    educ[i,"centroid_lat"] <- centr[2]
  }
}
educ <- educ[!is.na(educ$centroid_lon),]
educ <- educ[,c("osm_id","name","centroid_lon","centroid_lat")]
educ <- educ[educ$name!="Общежитие МФТИ",]

DistToUni <- distm(inner[,c("longitude","latitude")], educ[,c("centroid_lon","centroid_lat")])/1000
inner$DistToUni <- NA
inner$ClosestUni <- NA
for(i in 1:nrow(inner)){
  inner[i,"DistToUni"] <- min(DistToUni[i,])
  inner[i,"ClosestUni"] <- which.min(DistToUni[i,])
}

ggmap(MoscowCenterMap) +
  geom_point(aes(x=longitude, y=latitude, color=factor(ClosestUni)), size=1.5, inner, alpha=0.7) +
  theme(legend.position = "none") +
  geom_point(aes(x=centroid_lon, y=centroid_lat), size=1.5, educ, alpha=1, fill="black", shape=23)
  # labs(color="Кластеризация по ближайшему ВУЗу")
ggsave(here("UrbanMap3_EfimovAlexandr.pdf"))

# Пункт 5.
DistInner <- distm(inner[,c("longitude","latitude")], inner[,c("longitude","latitude")])/1000
for(i in 1:nrow(DistInner)){
  DistInner[i,i] <- NA
}
inner$mean_rent <- NA
inner$mean_discr <- NA
for(i in 1:nrow(DistInner)){
  inner[i,"mean_rent"] <- mean(inner[which(DistInner[,i] < 1),]$rent_amount)
  inner[i,"mean_discr"] <- mean(inner[which(DistInner[,i] < 1),]$ethnic_discrim)
}

m5.1 <- lm(log(rent_amount)~log(mean_rent), inner)
m5.2 <- glm(ethnic_discrim~mean_discr, inner, family="binomial")

stargazer(list(m5.1,m5.2), type="html", df=FALSE,
          out=here("UrbanTable1_EfimovAlexandr.html"))

# Пункт 6.
fm1 <- log(rent_amount)~ethnic_discrim
fm2 <- log(rent_amount)~ethnic_discrim+log(total_area)+one_room+two_rooms+three_rooms+
  firstfloor+log(floor_number)+repair+parking+blochnii1+kirpichnii3+
  monolitno_kirpichnii4+monolitnii5+panelnii6+stalinskii7+starii_fond8
fm3 <- log(rent_amount)~ethnic_discrim+log(total_area)+log(DistToCenter)+log(DistToUni)
fm4 <- log(rent_amount)~ethnic_discrim+log(total_area)+log(DistToCenter)+
  log(DistToUni)+as.factor(ClosestUni)
fm5 <- log(rent_amount)~ethnic_discrim+log(total_area)+log(DistToCenter)+
  log(DistToUni)+as.factor(ClosestUni)+one_room+two_rooms+three_rooms+
  firstfloor+log(floor_number)+repair+parking+blochnii1+kirpichnii3+
  monolitno_kirpichnii4+monolitnii5+panelnii6+stalinskii7+starii_fond8
m6.1 <- lm(fm1, inner)
m6.2 <- lm(fm2, inner)
m6.3 <- lm(fm3, inner)
m6.4 <- lm(fm4, inner)
m6.5 <- lm(fm5, inner)

clse = function(reg) {
  rob = sqrt(diag(vcovCL(reg, cluster=~ClosestUni, type = "HC1")))
  return(rob)
}

stargazer(list(m6.1, m6.2, m6.3, m6.4, m6.5), se=list(clse(m6.1),clse(m6.2),clse(m6.3),clse(m6.4),clse(m6.5)),
          type="html", align=TRUE, df=FALSE, out=here("UrbanTable2_EfimovAlexandr.html"),
          keep=c("ethnic_discrim","log","room","floor","repair","parking","blochnii1",
                 "kirpichnii","monolitnii5","panelnii6","stalinskii7","starii_fond8","Const"),
          order=c("ethnic_discrim","total_area","Dist","number"),
          add.lines=list(c("Индивидуальные эффекты ВУЗов","Нет","Нет","Нет","Да","Да"),
                         c("","","","","","")))

# Пункт 7.
psm <- glm(ethnic_discrim~log(total_area)+one_room+two_rooms+three_rooms+
             firstfloor+log(floor_number)+repair+parking+blochnii1+kirpichnii3+
             monolitno_kirpichnii4+monolitnii5+panelnii6+stalinskii7+starii_fond8+
             log(DistToCenter)+
             log(DistToUni)+as.factor(ClosestUni), inner, family="binomial")
inner$ps <- psm$fitted.values

# ggplot(data=inner, aes(x=ps, group=as.factor(ethnic_discrim), fill=as.factor(ethnic_discrim))) +
#   geom_density(adjust=1.5, alpha=.4) + 
#   guides(fill=guide_legend(title="Этническая дискриминация")) +
#   scale_fill_manual(values=c("green3", "red"))

inner$X <- 1:nrow(inner)
treated <- inner[inner$ethnic_discrim==1,]
control <- inner[inner$ethnic_discrim==0,]
treat.vs.control.dist <- DistInner[treated$X, control$X]

match.data <- function(n){
  score.diffs <- c()
  matched.control.idx <- c()
  treat.no.match <- c()
  for(i in 1:nrow(treated)){
    neighbors <- which(treat.vs.control.dist[i,] < 1)
    if (length(neighbors) >= n){
      ps.diff <- abs(control[neighbors,"ps"] - treated[i,"ps"])
      if (min(ps.diff) <= 0.05){
        score.diffs <- c(score.diffs,min(ps.diff))
        matched.control.idx <- c(matched.control.idx,neighbors[order(ps.diff)[1:n]])
      } else{
        treat.no.match <- c(treat.no.match,i)
      }
      
    } else{
      treat.no.match <- c(treat.no.match,i)
    }
  }
  return(rbind(treated[-treat.no.match,],control[unique(matched.control.idx),]))
}

matched.data1 <- match.data(n=1)
matched.data3 <- match.data(n=3)
m7.1 <- lm(log(rent_amount)~ethnic_discrim, matched.data3)
m7.2 <- lm(log(rent_amount)~ethnic_discrim+log(DistToCenter)+log(DistToUni)+
             as.factor(ClosestUni), matched.data3)
m7.3 <- lm(log(rent_amount)~ethnic_discrim, matched.data1)
m7.4 <- lm(log(rent_amount)~ethnic_discrim+log(DistToCenter)+log(DistToUni)+
             as.factor(ClosestUni), matched.data1)

stargazer(list(m7.1, m7.2, m7.3, m7.4), se=list(clse(m7.1),clse(m7.2),clse(m7.3),clse(m7.4)),
          type="html", out=here("UrbanTable3_EfimovAlexandr.html"),
          column.labels=c("N = 3","N = 1"), column.separate=c(2,2),
          keep=c("ethnic_discrim","log","log","Const"), align=TRUE, df=FALSE,
          add.lines=list(c("Индивидуальные эффекты ВУЗов","Нет","Да","Нет","Да"),
                         c("","","","","")))
