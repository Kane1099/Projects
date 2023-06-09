library(raster)
library(readxl)
library(dplyr)
library(spgwr)
library(stargazer)
library(rgdal)
library(maptools)
library(spdep)
library(sf)
library(spatialreg)
library(ggplot2)
library(tibble)
library(ggthemes)
library(plm)
library(lmtest)
library(splm)
library(panelr)
library(writexl)

geo_data <- raster::getData('GADM', country = 'Russia', level = 1)
geo_data@data$order<-1:dim(geo_data@data)[1]
      
geo_data@data$NL_NAME_1[is.na(geo_data@data$NL_NAME_1)] <- 'г. Москва'
df <- read_excel('project_data_final.xlsx')

#summary statistics
summary(df[c('crime', 'alcohol', 'divorce', 'grp', 'stud', 'unemp')])


#panel summary statistics
panel_data(df, id=region_geo, wave = year)


df<-merge(df, geo_data@data[,c("NL_NAME_1","order")], by.x="region_geo",by.y="NL_NAME_1")
df<-df %>% arrange(order,year)

df$crime_on1000<-df$crime/df$population*1000
df$alco_on1000<-df$alcohol/df$population*1000
df$grp_on<-df$grp/1000000
df$stud[df$stud==0]<-1
df$stud_on1000<-df$stud/1000



data_2016<-df %>% filter(df$year==2016)
data_2020<-df %>% filter(df$year==2020)

shp_2016 <- merge(geo_data,data_2016,by.x='NL_NAME_1',by.y='region_geo',all.x=FALSE)
shp_2020<-merge(geo_data,data_2020,by.x='NL_NAME_1',by.y='region_geo')



for (i in 1:83) {
  for (j in 1:length(slot(slot(shp_2016, "polygons")[[i]], "Polygons"))) {
    slot(slot(slot(shp_2016, "polygons")[[i]], "Polygons")[[j]], "coords")[,1] <- ifelse(
      slot(slot(slot(shp_2016, "polygons")[[i]], "Polygons")[[j]], "coords")[,1] < 0,
      slot(slot(slot(shp_2016, "polygons")[[i]], "Polygons")[[j]], "coords")[,1] + 360,
      slot(slot(slot(shp_2016, "polygons")[[i]], "Polygons")[[j]], "coords")[,1]
    )
  }
}

for (i in 1:83) {
  for (j in 1:length(slot(slot(shp_2020, "polygons")[[i]], "Polygons"))) {
    slot(slot(slot(shp_2020, "polygons")[[i]], "Polygons")[[j]], "coords")[,1] <- ifelse(
      slot(slot(slot(shp_2020, "polygons")[[i]], "Polygons")[[j]], "coords")[,1] < 0,
      slot(slot(slot(shp_2020, "polygons")[[i]], "Polygons")[[j]], "coords")[,1] + 360,
      slot(slot(slot(shp_2020, "polygons")[[i]], "Polygons")[[j]], "coords")[,1]
    )
  }
}


###карта для 2016 года

brks <- round(quantile(shp_2016$crime_on1000, probs=seq(0,1,0.2)), digits=2)
colours <- c("salmon1", "salmon2", "red3", "brown", "black")  
a_2016<-findInterval(shp_2016$crime_on1000,brks)
a_2016[a_2016==6]<-5
a_2016[a_2016==0]<-1
plot(shp_2016,xlim=c(30,180),ylim=c(30,100),col = colours[a_2016])  
legend( 'topleft', legend=leglabs(brks), fill = colours, bty="n")
title(main=paste("Количество преступлений на 1000 чел. в 2016 году"))

###карта для 2020 года
brks1 <- round(quantile(shp_2020$crime_on1000, probs=seq(0,1,0.2)), digits=2)
colours1 <- c("salmon1", "salmon2", "red3", "brown", "black")
a_2020<-findInterval(shp_2020$crime_on1000,brks)
a_2020[a_2020==6]<-5
a_2020[a_2020==0]<-1
plot(shp_2020,xlim=c(30,180),ylim=c(30,100),col = colours1[a_2020]) 
legend( 'topleft', legend=leglabs(brks), fill = colours1, bty="n")
title(main=paste("Количество преступлений на 1000 чел. в 2020 году"))

###индексы пространственной зависимости
coords_2016<- coordinates(shp_2016)
nb_knn_2016 <- knearneigh(coords_2016, k = 3) %>% knn2nb()
w_2016 <-nb2listw(nb_knn_2016, style="W", zero.policy=TRUE)
moran.test(shp_2016$crime_on1000,w_2016) #значим, <1 значит есть положительная пространственная зависимость
str(geary(shp_2016$crime_on1000, w_2016, length(nb_knn_2016), length(nb_knn_2016)-1,Szero(w_2016)))#меньше 1, положительная пространственная корреляция
moran.plot(shp_2016$crime_on1000,w_2016) #много внизу слева и справа сверху 

lmi<-localmoran_perm(shp_2016$crime_on1000,w_2016,na.action = na.exclude) #для каждого региона pv мало, есть пространственная зависимость
xd <- as.data.frame(lmi)
brks <- round(quantile(xd$Ii, probs=seq(0,1,0.2)), digits=2)
a_local<-findInterval(xd$Ii, brks)
a_local[a_local==6]<-5
a_local[a_local==0]<-1
plot(shp_2016,xlim=c(30,180),ylim=c(30,100),col = colours[a_local]) 
legend('topleft', legend=leglabs(brks), fill = colours, bty="n")
title(main=paste("Локальный индекс Морана в 2016 году"))


##panels
#pooling
model<-log(crime_on1000)~log(stud_on1000)+log(alco_on1000)+unemp+log(divorce)+log(grp_on)
gp <- plm(model,data = df, model = "pooling")
summary(gp)


#FE 
one_way <- plm(model,data = df, model = "within", effect = "individual")


clse = function(reg) { 
  G = length(unique(index(reg,"id")))
  N = length(index(reg,"id"))
  dfa = (G/(G - 1)) 
  rob = sqrt(diag(dfa*vcovHC(reg, method="arellano", type = "HC1", 
                             cluster = "group")))
  return(rob)
}

wi_time <- plm(model,data = df, model = "within", effect = "time")
wi <- plm(model,data = df, model = "within", effect = "twoways")

re <- plm(model,data = df, model = "random", effect = "individual")

stargazer(gp,one_way,wi_time,wi,re,se=list(sqrt(diag(gp$vcov)),clse(one_way),clse(wi_time),clse(wi),clse(re)),type='html',out='models.html')

pooltest(gp,one_way)#h1 принимается, лучше использовать FE 
pFtest(one_way, gp) #fe лучше
pFtest(wi_time, gp) #лучше не добавлять временной тренд! 
pFtest(wi, gp) #twoways лучше 
phtest(one_way, re) #fe лучше
phtest(wi, re) #fe лучше

#-->панели FE individual

coords<- coordinates(geo_data)
nb_knn <- knearneigh(coords, k = 3) %>% knn2nb()
W <- nb2mat(nb_knn, zero.policy=T)


pcdtest(one_way, test = c("cd")) #CDtest есть пространственная зависимость между регионами
slmtest(one_way, listw=W, test="lme") #есть пространственная корреляция в остатках 

within_ind_lag <- spml(model, data = df, listw = spdep::mat2listw(W),
                   model="within", spatial.error="none", Hess = FALSE, effect = 'individual',lag=TRUE)
summary(within_ind_lag) #значимы разводы + пространс лаг crime #SAR
impac1 <- spatialreg::impacts(within_ind_lag, listw=spdep::mat2listw(W, style="W"), time=5)
summary(impac1, short=T)


within_ind_lag_er <- spml(model, data = df, listw = spdep::mat2listw(W),
                       model="within", spatial.error="b", Hess = FALSE, effect = 'individual',lag=TRUE)
summary(within_ind_lag_er) #значимая пространст кор в ошибках,crime разводы  #SARAR
impac2 <- spatialreg::impacts(within_ind_lag_er, listw=spdep::mat2listw(W, style="W"), time=5)
summary(impac2, short=T, zstats=T)

within_ind_er <- spml(model, data = df, listw = spdep::mat2listw(W),
                          model="within", spatial.error="b", Hess = FALSE, effect = 'individual',lag=FALSE)
summary(within_ind_er) #значимая пространст кор в ошибках, разводы #SEM

one_way_lag <- plm(update(model,.~.+slag(unemp, listw=spdep::mat2listw(W)))
                          ,data = df, model = "within", effect = "individual") #SLX с безработицей
summary(one_way_lag) #пространственный лаг безработицы тоже значим

one_way_lag_div <- plm(update(model,.~.+slag(log(divorce), listw=spdep::mat2listw(W)))
                   ,data = df, model = "within", effect = "individual") #SLX с разводами
summary(one_way_lag_div) #пространственный лаг разводов тоже значим

# within_ind_lag_y <- spml(update(model,.~.+slag(log(divorce), listw=spdep::mat2listw(W))),listw=spdep::mat2listw(W),data=df,
#                        model="random", spatial.error="none", Hess = FALSE, effect = 'individual',lag=TRUE)
# summary(within_ind_lag_y)

within_time_lag <- spml(model, data = df, listw = spdep::mat2listw(W),
                       model="within", spatial.error="none", Hess = FALSE, effect = 'time',lag=TRUE)
summary(within_time_lag) #пространс лаг crime не значим

