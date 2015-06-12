library(mgcv)
library(BiodiversityR)
library(MASS)
library(raster)


##################
#                #
#   Mapeo PTM    #
#                #
##################


setwd("/home/dupas/Bureau/restauration/media/dupas/1To/IRD/ARTICLES/ECOFOR/pheromones/")
setwd("/media/dupas/1To/IRD/ARTICLES/ECOFOR/pheromones")
setwd("/media/1To/IRD/ARTICLES/ECOFOR/pheromones")
setwd("/home/legs/InsectVirus/")


xmin<--81;xmax<--69;ymin<--5;ymax<-13
resolution <- 0.008333333333333333333
Offset=c(round((90-ymax)/resolution,0),round((180+xmin)/resolution,0))
Region.dim=c(round((ymax-ymin)/resolution,0),round((xmax-xmin)/resolution,0))

insect <- read.table("datosFeromonas.csv",header=TRUE)
predictors <- stack(c("alt.bil",paste("bio_", 1:19, ".bil",sep="")))
predictors <- crop(predictors,extent(c(xmin,xmax,ymin,ymax)))
bioclim_data <- extract(predictors,insect[,c("X","Y")])
insect_bioclim <- cbind(insect,bioclim_data)
insect_bioclim <- insect_bioclim[which(abs(insect_bioclim[,c("Alt")]-insect_bioclim[,c("alt")])<200),]


dep="Phth"
pht <- na.omit(insect_bioclim[,-which(colnames(insect_bioclim)%in%c("Tec","Sym"))]);dim(pht)
formula.glm2 <- as.formula(paste(dep,"~Dias+Northern_Ecuador+Central_Colombia+Venezuela+",paste(c("alt",grep("bio_",colnames(insect_bioclim),value=TRUE)),sep="",collapse="+"),"+","I(",paste(c("alt",grep("bio_",colnames(insect_bioclim),value=TRUE)),sep="",collapse="^2)+I("),"^2)"))
fullmodelPht <- glm.nb(formula.glm2,data=insect_bioclim)
nullmodelPht <- glm.nb(paste(dep,"~1",sep=""),data=insect_bioclim)
resultPht2 <- step(fullmodelPht,k=log(length(fullmodelPht$residuals)))
#resultPht_f <- step(nullmodelPht,direction="forward",scope=list(lower=~1,upper=as.formula(paste("~",paste(c("alt",grep("bio",colnames(insect_bioclim),value=TRUE)),sep="",collapse="+"),"+","I(",paste(c("alt",grep("bio",colnames(insect_bioclim),value=TRUE)),sep="",collapse="^2)+I("),"^2)"))))
#resultPht$aic
resultPht2$aic
summary(resultPht2)
#summary(resultPht_f)
pht$phtPred <- predict(resultPht2,type="response")
glm.nb.PredObsPhth <- glm.nb(Phth~log(phtPred),data=pht)
jpeg("pht_obs_pred.jpg"); plot(pht$Phth,pht$phtPred, main = paste("Phthorimmaea operculella"),sub = paste ("R = ",cor(pht$Phth,pht$phtPred), "p = ",summary(glm.nb.PredObsPhth)$coefficients["log(phtPred)","Pr(>|z|)"]));dev.off()


dep="Tec"
tec <- na.omit(insect_bioclim[,-which(colnames(insect_bioclim)%in%c("Phth","Sym"))]);dim(tec)
formula.glm2 <- as.formula(paste(dep,"~Dias+Northern_Ecuador+Central_Colombia+Venezuela+",paste(c("alt",grep("bio",colnames(tec),value=TRUE)),sep="",collapse="+"),"+","I(",paste(c("alt",grep("bio",colnames(tec),value=TRUE)),sep="",collapse="^2)+I("),"^2)"))
fullmodelTec <- glm.nb(formula.glm2,data=tec)
resultTec2 <- step(fullmodelTec,k=log(length(fullmodelTec$residuals)))
#resultTec <- step(fullmodelTec)
#nullmodelTec <- glm.nb(paste(dep,"~1",sep=""),data=insect_bioclim)
#resultTec_f <- step(nullmodelTec,direction="forward",scope=list(lower=~1,upper=as.formula(paste("~",paste(c("alt",grep("bio",colnames(insect_bioclim),value=TRUE)),sep="",collapse="+"),"+","I(",paste(c("alt",grep("bio",colnames(insect_bioclim),value=TRUE)),sep="",collapse="^2)+I("),"^2)"))))
#resultTec$aic
resultTec2$aic
#summary(resultTec)
summary(resultTec2)
#summary(resultTec_f)
tec$tecPred <- predict(resultTec2,type="response")
glm.nb.PredObsTec <- glm.nb(Tec~log(tecPred),data=tec)
jpeg("tecia_obs_pred.jpg"); plot(tec$Tec,tec$tecPred, main = paste("Tecia solanivora"),sub = paste ("R = ",cor(tec$Tec,tec$tecPred), "p = ",summary(glm.nb.PredObsTec)$coefficients["log(tecPred)","Pr(>|z|)"]));dev.off()

dep="Sym"
sym <- na.omit(insect_bioclim[,-which(colnames(insect_bioclim)%in%c("Phth","Tec","Sitio","Provincia","Pais","X","Y"))]);dim(sym)
formula.glm2 <- as.formula(paste(dep,"~Dias+Northern_Ecuador+",paste(c("alt",grep("bio",colnames(sym),value=TRUE)),sep="",collapse="+"),"+","I(",paste(c("alt",grep("bio",colnames(sym),value=TRUE)),sep="",collapse="^2)+I("),"^2)"))
formula.glm2 <- as.formula(paste(dep,"~Dias+Northern_Ecuador+alt+bio_1+bio_2+bio_3+bio_4+bio_6+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18+bio_19+I(alt^2)+I(bio_1^2)+I(bio_2^2)+I(bio_3^2)+I(bio_4^2)+I(bio_6^2)+I(bio_7^2)+I(bio_8^2)+I(bio_9^2)+I(bio_10^2)+I(bio_11^2)+I(bio_12^2)+I(bio_13^2)+I(bio_14^2)+I(bio_15^2)+I(bio_16^2)+I(bio_17^2)+I(bio_18^2)+I(bio_19^2)"))
fullmodelSym <- glm.nb(formula.glm2,data=sym)
resultSym2 <- step(fullmodelSym,k=log(fullmodelSym$df.null+1))
#resultSym <- step(fullmodelSym)
nullmodelSym <- glm.nb(paste(dep,"~1",sep=""),data=insect_bioclim)
#resultSym_f <- step(nullmodelSym,direction="forward",scope=list(lower=~1,upper=as.formula(paste("~",paste(c("alt",grep("bio",colnames(insect_bioclim),value=TRUE)),sep="",collapse="+"),"+","I(",paste(c("alt",grep("bio",colnames(insect_bioclim),value=TRUE)),sep="",collapse="^2)+I("),"^2)"))))
#resultSym$aic
resultSym2$aic
#summary(resultSym)
summary(resultSym2)
#summary(resultSym_f)
sym$symPred <- predict(resultSym2,type="response")
glm.nb.PredObsSym <- glm.nb(Sym~log(symPred),data=sym)
jpeg("sym_obs_pred.jpg"); plot(sym$Sym,sym$symPred, main = paste("Symetrischemma tangolias"),sub = paste ("R = ",cor(sym$Sym,sym$symPred), "p = ",summary(glm.nb.PredObsSym)$coefficients["log(symPred)","Pr(>|z|)"]));dev.off()

##
##
## Prediction GLM
##
##

PTMActGLM <- crop(stack(c("bio_1.bil","bio_1.bil","bio_1.bil")),extent(c(xmin,xmax,ymin,ymax)) )
datosPTM <- as.data.frame(cbind(Central_Ecuador=1/4, Northern_Ecuador=1/4, Central_Colombia=1/4, Venezuela=1/4,Dias=21,values(predictors)))
values(PTMActGLM) <-  as.matrix(data.frame(Phth=predict(resultPht2,datosPTM,type="response"),
                                           Tec=predict(resultTec2,datosPTM,type="response"),
                                           Sym=predict(resultSym2,datosPTM,type="response")))
writeRaster(PTMActGLM,"PTMActGLM.tif",driver="GTiff",overwrite=TRUE)
PTMActGLMSymNotInvadedColombia <- PTMActGLM
values(PTMActGLMSymNotInvadedColombia)[which(coordinates(PTMActGLM)[,"y"]>=1.25),3] <- 0  # Symetrischemma is absent in central colombia
writeRaster(PTMActGLMSymNotInvadedColombia,"PTMActGLMSymNotInvadedColombia.tif",overwrite=TRUE)

phtMap <- subset(PTMActGLMSymNotInvadedColombia,"Phth")
tecMap <- subset(PTMActGLMSymNotInvadedColombia,"Tec")
symMap <- subset(PTMActGLMSymNotInvadedColombia,"Sym")

writeRaster(phtMap,"phtActGLM.tif",driver="GTiff",overwrite=TRUE)
writeRaster(tecMap,"tecActGLM.tif",driver="GTiff",overwrite=TRUE)
writeRaster(symMap,"symActGLMSymNotInvadedColombia.tif",driver="GTiff",overwrite=TRUE)

phtMap <- raster("phtActGLM.tif")
tecMap <- raster("tecActGLM.tif")
symMap <- raster("symActGLM.tif")
#symMap <- tecMap; values(symMap) <- values(PTMMap)[,"PTMActGLM.3"] ; names(symMap)  <- "symActGLM"
raster("symActGLMSymNotInvadedColombia.tif")
#PTMMap <- stack("PTMActGLM.tif")

#plot(raster("phtActGLM.tif"))

###############
# Mapeo Virus #
###############
################
#              #

#   Mapeo p89  #
#              #
################


p89 <- read.table("virusecofor.txt")
dep="p89"
bioclim_p89 <- extract(predictors,p89[,c("X","Y")])
bioclim_p89 <- extract(subset(predictors,c("alt","bio_1","bio_6","bio_9","bio_12","bio_18")),p89[,c("X","Y")])
bioclim_p89 <- cbind(p89,bioclim_p89)
formula.glm <- as.formula(paste(dep,"~",paste(c("alt",grep("bio",colnames(bioclim_p89),value=TRUE)),sep="",collapse="+")))
formula.glm2 <- as.formula(paste(dep,"~",paste(c("alt",grep("bio",colnames(bioclim_p89),value=TRUE)),sep="",collapse="+"),"+","I(",paste(c("alt",grep("bio",colnames(bioclim_p89),value=TRUE)),sep="",collapse="^2)+I("),"^2)"))
fullmodelp89 <- glm(formula.glm,data=bioclim_p89,family=binomial)
fullmodelp89 <- glm(formula.glm2,data=bioclim_p89,family=binomial)
nullmodelp89 <- glm("p89~1",data=bioclim_p89,family=binomial)
resultp89_f <- step(nullmodelp89,direction="forward",scope=list(lower=~1,upper=as.formula(paste("~",paste(c("alt",grep("bio",colnames(bioclim_p89),value=TRUE)),sep="",collapse="+"),"+","I(",paste(c("alt",grep("bio",colnames(bioclim_p89),value=TRUE)),sep="",collapse="^2)+I("),"^2)"))))
resultp89_b <- step(resultp89_f,k=log(resultp89_f$df.null+1))
resultp89_2 <- step(fullmodelp89,k=log(fullmodelp89$df.null+1))
p89plot <- predict(predictors,resultp89_2,type="response")
plot(p89plot)
writeRaster(p89plot,"p89plot.tif",driver="GTiff",overwrite=TRUE)
p89plot_f <- predict(predictors,resultp89_f,type="response")
writeRaster(p89plot_f,"p89plot_f.tif",driver="GTiff",overwrite=TRUE)


#
# Correlations
#

# ranges

library(matrixStats)

insect_bioclimPTM <- insect_bioclim[-c(grep("Neg",insect_bioclim[,1]),grep("abs",insect_bioclim[,1])),]

bio_1p89range <- range(bioclim_p89[!is.na(bioclim_p89$p89),"bio_1"])
bio_6p89range <- range(bioclim_p89[!is.na(bioclim_p89$p89),"bio_6"])
bio_12p89range <- range(bioclim_p89[!is.na(bioclim_p89$p89),"bio_12"])
p89Var  <- c(grep("alt",names(resultp89_2$coefficients),value=TRUE)[!grepl("I\\(",grep("alt",names(resultp89_2$coefficients),value=TRUE))],
             grep("bio_",names(resultp89_2$coefficients),value=TRUE)[!grepl("I\\(",grep("bio_",names(resultp89_2$coefficients),value=TRUE))])
p89ranges <- colRanges(as.matrix(bioclim_p89[,p89Var])); rownames(p89ranges) <- p89Var


SymVar <- c(grep("alt",names(resultSym2$coefficients),value=TRUE)[!grepl("I\\(",grep("alt",names(resultSym2$coefficients),value=TRUE))],
            grep("bio_",names(resultSym2$coefficients),value=TRUE)[!grepl("I\\(",grep("bio_",names(resultSym2$coefficients),value=TRUE))])
Symranges <- colRanges(as.matrix(insect_bioclimPTM[,SymVar])); rownames(Symranges) <- SymVar

TecVar <- c(grep("alt",names(resultTec2$coefficients),value=TRUE)[!grepl("I\\(",grep("alt",names(resultTec2$coefficients),value=TRUE))],
            grep("bio_",names(resultTec2$coefficients),value=TRUE)[!grepl("I\\(",grep("bio_",names(resultTec2$coefficients),value=TRUE))])
Tecranges <- colRanges(as.matrix(insect_bioclimPTM[,TecVar])); rownames(Tecranges) <- TecVar

PhtVar <- c(grep("alt",names(resultPht2$coefficients),value=TRUE)[!grepl("I\\(",grep("alt",names(resultPht2$coefficients),value=TRUE))],
            grep("bio_",names(resultPht2$coefficients),value=TRUE)[!grepl("I\\(",grep("bio_",names(resultPht2$coefficients),value=TRUE))])
Phtranges <- colRanges(as.matrix(insect_bioclimPTM[,PhtVar])); rownames(Phtranges) <- PhtVar

allRanges <- colRanges(na.omit(as.matrix(datosPTM[,c("alt",grep("bio_",colnames(datosPTM),value=TRUE))])))
colnames(allRanges) <- c("allMin","allMax")
rownames(allRanges) <- c("alt",grep("bio_",colnames(datosPTM),value=TRUE))
allRanges <- cbind(allRanges,allRanges,allRanges,allRanges,allRanges,allRanges)
colnames(allRanges) <- c("allMin","allMax","p89Min","p89Max","phtMin","phtMax","tecMin","tecMax","symMin","symMax","Min","Max")
allRanges[rownames(p89ranges),c("p89Min","p89Max")] <- p89ranges 
allRanges[rownames(Phtranges),c("phtMin","phtMax")] <- Phtranges 
allRanges[rownames(Tecranges),c("tecMin","tecMax")] <- Tecranges
allRanges[rownames(Symranges),c("symMin","symMax")] <- Symranges
allRanges[,"Min"] <- rowMaxs(allRanges[,2*1:4+1])
allRanges[,"Max"] <- rowMins(allRanges[,2*1:4+2])

p89TecRanges <- data.frame(Min=rowMaxs(allRanges[,c("tecMin","p89Min")]),Max=rowMins(allRanges[,c("tecMax","p89Max")]));rownames(p89TecRanges)=rownames(allRanges)

validTecP89 <- (predictors>p89TecRanges[,"Min"])&(predictors<p89TecRanges[,"Max"])
validTecP89 <- prod(validTecP89)
plot(validTecP89)
writeRaster(validTecP89,"validRangeTecP89.tif",driver="GTiff",overwrite=TRUE)

validAll <- (predictors>allRanges[,"Min"])&(predictors<allRanges[,"Max"])
validAll <- prod(validAll)
plot(validAll)
writeRaster(validAll,"validRange.tif",driver="GTiff",overwrite=TRUE)

#
# Correlations
#

Tecp89 <- stack(tecMap,p89plot)
validDatosTecP89 <- Tecp89[validTecP89]
writeRaster(validTecP89,"validTecP89map.tif")
corTecP89 <- cor(na.omit((log(values(Tecp89)[,1])),na.omit(log(values(Tecp89)[,2])/(1-log(values(Tecp89)[,2])))))
pvalTecP89 <- lm(logit(values(Tecp89)[,2])~log(values(TecP89)[,1]))$pvalue
svg("validTecP89plot.svg")
plot(log(values(TecP89)[,1]),logit(values(Tecp89)[,2]), main="Relation between niche model predicted Tecia solanivora abundance and niche model predicted PhoGV prevalence on Tecia solanivora across sampled bioclimatic and altitudinal range", xlab="log(Tecia solanivora)",ylab="logit(phoGV prevalence)",sub=paste("R = ",corTecP89,", p-value = ",pvalTecP89))
dev.off()

symp89 <- stack(symMap,p89plot)
validDatossymP89 <- symp89[validsymP89]
writeRaster(validDatossymP89,"validsymP89map.tif")
svg("validsymP89plot.svg")
corsymP89 <- cor(log(values(symP89)[,1]),logit(values(symp89)[,2]))
pvalsymP89 <- lm(logit(values(symp89)[,2])~log(values(symP89)[,1]))$pvalue
plot(log(values(symP89)[,1]),logit(values(symp89)[,2]), main="Relation between niche model predicted Symatrischemma tangolias abundance and niche model predicted PhoGV prevalence on Tecia solanivora across sampled bioclimatic and altitudinal range", xlab="log(Symetrischemma tangolias)",ylab="logit(phoGV prevalence)",sub=paste("R = ",corsymP89,", p-value = ",pvalsymP89))
dev.off()

phtp89 <- stack(phtMap,p89plot)
validDatosphtP89 <- phtp89[validphtP89]
writeRaster(validDatosphtP89,"validphtP89map.tif")
svg("validphtP89plot.svg")
corphtP89 <- cor(log(values(phtP89)[,1]),logit(values(phtp89)[,2]))
pvalphtP89 <- lm(logit(values(phtp89)[,2])~log(values(phtP89)[,1]))$pvalue
plot(log(values(phtP89)[,1]),logit(values(phtp89)[,2]), main="Relation between niche model predicted Phthorimmaea operculella abundance and niche model predicted PhoGV prevalence on Tecia solanivora across sampled bioclimatic and altitudinal range", xlab="log(Phthorimmaea operculella)",ylab="logit(phoGV prevalence)",sub=paste("R = ",corphtP89,", p-value = ",pvalphtP89))
dev.off()
