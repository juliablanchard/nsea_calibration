#read saved data

f_history<-as.matrix(readRDS("data/FmatWeightedInterpolated.rds"))[1:73,]
row.names(f_history)<-f_history[,1]
f_history<-f_history[as.character(1947:2019),-1]


# check first value 
firstRecordedEffort <- NULL
for(iSpecies in 1:ncol(f_history))
  firstRecordedEffort <- c(firstRecordedEffort,f_history[which(!is.na(f_history[,iSpecies]))[1],iSpecies]) # first value of F for each species (different years)
# Gurnard has no data so inputing same value as whiting
firstRecordedEffort[8] <- firstRecordedEffort[6]


indexFirstEffort<-NULL
for(iSpecies in 1:ncol(f_history))indexFirstEffort<-c(indexFirstEffort,which(!is.na(f_history[,iSpecies]))[1])
indexFirstEffort <-indexFirstEffort+100

pre_effort <- matrix(NA, byrow = TRUE, nrow = 100,
                     ncol = ncol(f_history), dimnames = list(1847:1946))

effort <- as.array(rbind(pre_effort, f_history))
colnames(effort)<-gear_params$species
rownames(effort)<-1847:2019

# Historical effort

# Use logistic curve (fit to Fmort data through time), with a descending limb
# starting from 1850

# f(y) ~ F/(1+ exp(-y0(y-y1))) * F/(1 + exp(-y2(y-y3))
# where:
# y1, the year value of the sigmoid's midpoint;
# F, the curve's maximum value;
# k, the logistic growth rate or steepness of the curve


# Move This to another file...
# visually fit double logistic curves to fishing moratlit rates ( see Svuwalski et al PNAS)

feffort<- effort

yr<-1847:2019

feffort[,"Sprat"]<-(1.5/(1+ exp(-0.5*(yr-1970))))
#*(2/(1 + exp(0.1*(yr-2005))))
plot(1847:2019,effort[,"Sprat"])
points(1847:2019,feffort[,"Sprat"],col="red")

feffort[,"Sandeel"]<-(1/(1+ exp(-1*(yr-1980))))*(0.3/(1 + exp(0.1*(yr-2005))))
plot(1847:2019,effort[,"Sandeel"])
points(1847:2019,feffort[,"Sandeel"],col="red")


feffort[,"N.pout"]<-(1.3/(1+ exp(-0.9*(yr-1980))))*(1.3/(1 + exp(0.1*(yr-1995))))
feffort["2019","N.pout"]<-1.2
plot(1847:2019,effort[,"N.pout"])
points(1847:2019,feffort[,"N.pout"],col="red")

feffort[,"Herring"]<-0.1 +(2/(1+ exp(-0.3*(yr-1963))))*(1/(1 + exp(0.1*(yr-1970))))
plot(1847:2019,effort[,"Herring"])
points(1847:2019,feffort[,"Herring"],col="red")

feffort[,"Dab"]<-0.1 + (1/(1+ exp(-1*(yr-2002))))
plot(1847:2019,effort[,"Dab"])
points(1847:2019,feffort[,"Dab"],col="red")

feffort[,"Whiting"]<-0.1 +(0.9/(1+ exp(-0.2*(yr-1976))))*(1.05/(1 + exp(0.2*(yr-2000))))
plot(1847:2019,effort[,"Whiting"])
points(1847:2019,feffort[,"Whiting"],col="red")

feffort[,"Sole"]<- 0.1 + (0.7/(1+ exp(-0.2*(yr-1963))))*(0.9/(1 + exp(0.1*(yr-2015))))
plot(1847:2019,effort[,"Sole"])
points(1847:2019,feffort[,"Sole"],col="red")

feffort[,"Gurnard"]<-(1.5/(1+ exp(-0.5*(yr-1978))))*(1.3/(1 + exp(0.1*(yr-2000))))
plot(1847:2019,effort[,"Gurnard"])
points(1847:2019,feffort[,"Gurnard"],col="red")

feffort[,"Plaice"]<-0.1+(0.7/(1+ exp(-0.1*(yr-1970))))*(1/(1 + exp(0.2*(yr-2002))))
plot(1847:2019,effort[,"Plaice"])
points(1847:2019,feffort[,"Plaice"],col="red")


feffort[,"Haddock"]<-0.1+(1/(1+ exp(-0.2*(yr-1968))))*(1/(1 + exp(0.1*(yr-2005))))
plot(1847:2019,effort[,"Haddock"])
points(1847:2019,feffort[,"Haddock"],col="red")


feffort[,"Cod"]<-0.1+(1/(1+ exp(-0.5*(yr-1975))))
#*(1.2/(1 + exp(0.09*(yr-2003))))
plot(1847:2019,effort[,"Cod"])
points(1847:2019,feffort[,"Cod"],col="red")

feffort[,"Saithe"]<-0.1 + (0.8/(1+ exp(-0.2*(yr-1973))))*(0.9/(1 + exp(0.09*(yr-2003))))
plot(1847:2019,effort[,"Saithe"])
points(1847:2019,feffort[,"Saithe"],col="red")

# this takes the extrapolated feffort to fill in historical period, otherwise use feffort as a smoother function
for(iSpecies in 1:12)
  effort[c(1:indexFirstEffort[iSpecies]),iSpecies]<-feffort[c(1:indexFirstEffort[iSpecies]),iSpecies]


## no fishing during war years
# effort[as.character(1939:1945),]<-0
# effort[as.character(1914:1918),]<-0

saveRDS(effort,"effortTime.RDS")
saveRDS(feffort,"logistic_effortTime.RDS")




# ######### Below are the links to most recent ICES North Sea stock assessment data that I used for this model:
# 
# Overall database:
#   https://standardgraphs.ices.dk/stockList.aspx
# 
# 
# Sprat: https://standardgraphs.ices.dk/ViewSourceData.aspx?key=13322
# 
# Sandeel : https://standardgraphs.ices.dk/ViewCharts.aspx?key=13303
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13301
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13298
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13304
# 
# N.pout - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13166
# 
# Herring - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13422
# 
# Dab - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13184
# 
# Whiting - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13525
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13506
# 
# Sole - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13743
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13495
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13828
# 
# Gurnard - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13493
# 
# Plaice - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13484
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13744
# 
# Haddock - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13204
# 
# Cod - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13740
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13838
# 
# Saithe - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13511


# add benthic resource
# use thermizer


