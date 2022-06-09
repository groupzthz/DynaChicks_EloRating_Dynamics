getwd()
setwd("//nas-vetsuisse.campus.unibe.ch/vetsuisse/Benutzer/kg19h146/PROJECT2/focalSelection/R-Elo")
#setwd("H:/PROJECT2/focalSelection/R-Elo")

#clear workspace
rm(list = ls())
#install.packages("EloRating")
#install_github("gobbios/EloRating", build_vignettes = TRUE)
library(aniDom) # randomised Elo-rating
library(EloRating) # original Elo rating
library(domstruc) # focus and position (who interacts with whom)
library(data.table) # data.frame better
#library(compete) # needed for aniDom (poisson) calculation
library(ggplot2) # plotting
library(lme4) # mixed models
library(multcomp) # posthoc analysis
#library(ggfortify)
#library(useful)
library(DHARMa) # model diagnistics
library(moments) # normal asumptions
library(MuMIn) # model comparison
library(irr) # ICC calculations
source('helper_functions_Elo.R')
set.seed(42)

##### Loading and preparing Data ###########
liveData <- fread("InteractionsLive.csv",header = TRUE, sep = ";")
videoData <- fread("InteractionsVideo.csv",header = TRUE, sep = ";")
Individ <- fread("Individuals.csv",header = TRUE, sep = ";")
dataMature <- fread("InteractionsUpdate.csv",header = TRUE, sep = ";")

## relevant Data extraction
lD = liveData[, 1:5]
vD = videoData[, 1:5]

dataYoung = rbind(liveData[, 1:5], videoData[, 1:5])

#check Data
str(dataYoung)
str(dataMature)
length(unique(dataYoung$Date)) == 15
length(unique(dataMature$Date)) == 4
unique(dataYoung$Pen)
unique(dataMature$Pen)

#convert Date
dataYoung$Date <- as.Date(as.character(dataYoung$Date), format = "%d.%m.%Y")
dataMature$Date <- as.Date(as.character(dataMature$Date), format = "%d.%m.%Y")

#check correct order
dataYoung = dataYoung[order(Date),]
dataMature = dataMature[order(Date),]
unique(dataYoung$Date)
unique(dataMature$Date)

#split into pens
PenAYt = dataYoung[Pen == "A", c(1, 3, 4,5)]
PenBYt = dataYoung[Pen == "B", c(1, 3, 4,5)]
PenCYt = dataYoung[Pen == "C", c(1, 3, 4,5)]
PenDYt = dataYoung[Pen == "D", c(1, 3, 4,5)]
PenEYt = dataYoung[Pen == "E", c(1, 3, 4,5)]
PenFYt = dataYoung[Pen == "F", c(1, 3, 4,5)]

PenAMt = dataMature[Pen == "A", c(1, 3, 4,5)]
PenBMt = dataMature[Pen == "B", c(1, 3, 4,5)]
PenCMt = dataMature[Pen == "C", c(1, 3, 4,5)]
PenDMt = dataMature[Pen == "D", c(1, 3, 4,5)]
PenEMt = dataMature[Pen == "E", c(1, 3, 4,5)]
PenFMt = dataMature[Pen == "F", c(1, 3, 4,5)]

#check entered Individuals
PenAY = PenAYt[(Winner %in% unique(Individ$ID[Individ$Pen == "A"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "A"])),]
PenBY = PenBYt[(Winner %in% unique(Individ$ID[Individ$Pen == "B"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "B"])),]
PenCY = PenCYt[(Winner %in% unique(Individ$ID[Individ$Pen == "C"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "C"])),]
PenDY = PenDYt[(Winner %in% unique(Individ$ID[Individ$Pen == "D"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "D"])),]
PenEY = PenEYt[(Winner %in% unique(Individ$ID[Individ$Pen == "E"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "E"])),]
PenFY = PenFYt[(Winner %in% unique(Individ$ID[Individ$Pen == "F"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "F"])),]

PenAM = PenAMt[(Winner %in% unique(Individ$ID[Individ$Pen == "A"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "A"])),]
PenBM = PenBMt[(Winner %in% unique(Individ$ID[Individ$Pen == "B"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "B"])),]
PenCM = PenCMt[(Winner %in% unique(Individ$ID[Individ$Pen == "C"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "C"])),]
PenDM = PenDMt[(Winner %in% unique(Individ$ID[Individ$Pen == "D"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "D"])),]
PenEM = PenEMt[(Winner %in% unique(Individ$ID[Individ$Pen == "E"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "E"])),]
PenFM = PenFMt[(Winner %in% unique(Individ$ID[Individ$Pen == "F"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "F"])),]

#create full data sets
PenAall = rbind(PenAY,PenAM)
PenBall = rbind(PenBY,PenBM)
PenCall = rbind(PenCY,PenCM)
PenDall = rbind(PenDY,PenDM)
PenEall = rbind(PenEY,PenEM)
PenFall = rbind(PenFY,PenFM)


#exclusion of animals unidentified after barn switch for full data:
PenAall = PenAall[Winner != "ZA" & Loser != "ZA"]
PenBall = PenBall[Winner != "TU" & Winner != "VA" & Loser != "TU" & Loser != "VA", ]
PenCall = PenCall[Winner != "TX" & Winner != "ZZ" & Winner != "OL" & Loser != "TX" &
                    Loser != "ZZ" & Loser != "OL", ]
#exclusion of double tagged animal during Mature
PenAM = PenAM[Winner != "ZA" & Loser != "ZA"]


rm(liveData, videoData, dataMature, lD, vD, dataYoung, 
   PenAYt, PenBYt, PenCYt, PenDYt, PenEYt, PenFYt,
   PenAMt, PenBMt, PenCMt, PenDMt, PenEMt, PenFMt)

################ Diagnostics ###############################################################

# Dataset analytics:

diagnAY = dataset_diagnostics(PenAY, Individ$ID[Individ$Pen == 'A'])
diagnBY = dataset_diagnostics(PenBY, Individ$ID[Individ$Pen == 'B'])
diagnCY = dataset_diagnostics(PenCY, Individ$ID[Individ$Pen == 'C'])
diagnDY = dataset_diagnostics(PenDY, Individ$ID[Individ$Pen == 'D'])
diagnEY = dataset_diagnostics(PenEY, Individ$ID[Individ$Pen == 'E'])
diagnFY = dataset_diagnostics(PenFY, Individ$ID[Individ$Pen == 'F'])

diagnAM = dataset_diagnostics(PenAM, Individ$ID[Individ$Pen == 'A'])
diagnBM = dataset_diagnostics(PenBM, Individ$ID[Individ$Pen == 'B'])
diagnCM = dataset_diagnostics(PenCM, Individ$ID[Individ$Pen == 'C'])
diagnDM = dataset_diagnostics(PenDM, Individ$ID[Individ$Pen == 'D'])
diagnEM = dataset_diagnostics(PenEM, Individ$ID[Individ$Pen == 'E'])
diagnFM = dataset_diagnostics(PenFM, Individ$ID[Individ$Pen == 'F'])

diagnAall = dataset_diagnostics(PenAall, Individ$ID[Individ$Pen == 'A'])
diagnBall = dataset_diagnostics(PenBall, Individ$ID[Individ$Pen == 'B'])
diagnCall = dataset_diagnostics(PenCall, Individ$ID[Individ$Pen == 'C'])
diagnDall = dataset_diagnostics(PenDall, Individ$ID[Individ$Pen == 'D'])
diagnEall = dataset_diagnostics(PenEall, Individ$ID[Individ$Pen == 'E'])
diagnFall = dataset_diagnostics(PenFall, Individ$ID[Individ$Pen == 'F'])
# maybe add vector containing which data is which age period


##################### Elo rating ################################################################
# RATING
ratingAY = elo_analysis(PenAY, Individ$ID[Individ$Pen == 'A']) #carefull!
ratingBY = elo_analysis(PenBY, Individ$ID[Individ$Pen == 'B']) #carefull!
ratingCY = elo_analysis(PenCY, Individ$ID[Individ$Pen == 'C']) #carefull!
ratingDY = elo_analysis(PenDY, Individ$ID[Individ$Pen == 'D']) #carefull!
ratingEY = elo_analysis(PenEY, Individ$ID[Individ$Pen == 'E'])
ratingFY = elo_analysis(PenFY, Individ$ID[Individ$Pen == 'F'])

ratingAM = elo_analysis(PenAM, Individ$ID[Individ$Pen == 'A']) #carefull!
ratingBM = elo_analysis(PenBM, Individ$ID[Individ$Pen == 'B'])
ratingCM = elo_analysis(PenCM, Individ$ID[Individ$Pen == 'C'])
ratingDM = elo_analysis(PenDM, Individ$ID[Individ$Pen == 'D'])
ratingEM = elo_analysis(PenEM, Individ$ID[Individ$Pen == 'E'])
ratingFM = elo_analysis(PenFM, Individ$ID[Individ$Pen == 'F'])

ratingAall = elo_analysis(PenAall, Individ$ID[Individ$Pen == 'A'], all = TRUE)
ratingBall = elo_analysis(PenBall, Individ$ID[Individ$Pen == 'B'], all = TRUE)
ratingCall = elo_analysis(PenCall, Individ$ID[Individ$Pen == 'C'], all = TRUE)
ratingDall = elo_analysis(PenDall, Individ$ID[Individ$Pen == 'D'], all = TRUE)
ratingEall = elo_analysis(PenEall, Individ$ID[Individ$Pen == 'E'], all = TRUE)
ratingFall = elo_analysis(PenFall, Individ$ID[Individ$Pen == 'F'], all = TRUE)


#######################################################################################

#### EXAMPLE PEN E #################
plotdata = data.table( Names= ratingEM$Individuals$ID,
                      Mature = ratingEM$Individuals$elos,
                      rankM = ratingEM$Individuals$rank,
                      Young = ratingEY$Individuals$elos,
                      rankY = ratingEY$Individuals$rank,
                      All = ratingEall$rating_rand[sort(names(ratingEall$rating_rand))])

plotdata[, Diff := rankM - rankY ]

mean = mean(abs(plotdata$Diff))

#plot the Young ranks vs Mature ranks of Pen E
 g1 = ggplot(data = plotdata, mapping = aes(x = rankM, y = rankY)) + 
  geom_abline(intercept = 0, slope = 1, size = 1, linetype = "dashed", colour = 'red')+
  geom_abline(intercept = mean, slope = 1, size = 1, linetype = "dashed", colour = 'grey')+
  geom_abline(intercept = -mean, slope = 1, size = 1, linetype = "dashed", colour = 'grey')+
  theme_classic(base_size = 18)+
  labs(x = 'Dominance rank Mature', y= 'Dominance rank Young')+
  geom_point() +
  geom_text(
    label=plotdata$Names,
    nudge_x=0.5, nudge_y=0.5, size = 5)


# Elo rating Young
mean.eloY = rowMeans(ratingEY$randElo)
identitiesY = rownames(ratingEY$randElo)
identitiesY <- identitiesY[order(mean.eloY)]
CIsY <- apply(ratingEY$randElo,1,quantile,c(0.025,0.975),na.rm=TRUE)
CIsY <- CIsY[,order(mean.eloY)]
mean.eloY <- mean.eloY[order(mean.eloY)]

plotTable = data.table(number = 1:20,  
                       ranks = mean.eloY,
                       IDs = identitiesY
)

g2 = ggplot(plotTable, aes(x = number, y= ranks))+
  geom_errorbar(aes(x= number, ymin=CIsY[1,], ymax=CIsY[2,]), width=.1)+
  geom_point(shape = 21, size = 5, colour = "white", fill = "white", stroke = 1)+
  
  geom_text(
    label=plotTable$IDs,
    size = 4)+
  labs(y="Elo-rating Young")+
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank())  

#Elo rating Mature
mean.eloM = rowMeans(ratingEM$randElo)
identitiesM = rownames(ratingEM$randElo)
identitiesM <- identitiesM[order(mean.eloM)]
CIsM <- apply(ratingEM$randElo,1,quantile,c(0.025,0.975),na.rm=TRUE)
CIsM <- CIsM[,order(mean.eloM)]
mean.eloM <- mean.eloM[order(mean.eloM)]

plotTable2 = data.table(number = 1:20,  
                       ranks = mean.eloM,
                       IDs = identitiesM
)

g3 = ggplot(plotTable2, aes(x = number, y= ranks))+
  geom_errorbar(aes(x= number, ymin=CIsM[1,], ymax=CIsM[2,]), width=.1)+
  geom_point(shape = 21, size = 5, colour = "white", fill = "white", stroke = 1)+
  scale_y_continuous(limits = c(-450, 650))+
  geom_text(
    label=plotTable2$IDs,
    size = 4)+
  labs(x = "Individuals", y="Elo-rating Mature")+
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  



library(patchwork)
design <- "
  13
  23
    "

fig2 = g2+ g3 + g1 + plot_layout(design = design, width = c(1, 1.4))
fig2[[2]] <- fig2[[2]] + plot_layout(tag_level = 'keep')
  fig2 + plot_annotation(tag_levels = c("A", "1"))&  theme(plot.tag = element_text(size = 18, face = "bold"))

cor.test(plotdata$rankM, plotdata$rankY, method = "spearman")


#############################################################################
##### Behaviour prediction of Mature Elo ############

#dataset
behaviourData = rbind(ratingAM$Individuals, ratingBM$Individuals, ratingCM$Individuals, ratingDM$Individuals,
                                 ratingEM$Individuals, ratingFM$Individuals)
behaviourData[, Pen := factor(c(rep('A', length(ratingAM$Individuals$rank)), 
                           rep('B',length(ratingBM$Individuals$rank)), 
                           rep('C',length(ratingCM$Individuals$rank)),
                           rep('D', 20), rep('E',20), rep('F',20)))]

behaviourData[, Condition := factor(c(rep("large", length(Pen)-60), 
                                 rep("small", 60)))]
behaviourData[, WinsY := c(ratingAY$Individuals$Wins, ratingBY$Individuals$Wins, 
                               ratingCY$Individuals$Wins, ratingDY$Individuals$Wins,
                               ratingEY$Individuals$Wins, ratingFY$Individuals$Wins)]
behaviourData[, LossesY := c(ratingAY$Individuals$Losses, ratingBY$Individuals$Losses, 
                               ratingCY$Individuals$Losses, ratingDY$Individuals$Losses,
                               ratingEY$Individuals$Losses, ratingFY$Individuals$Losses)]
behaviourData = behaviourData[!is.na(rank),]
behaviourData[, scaleElos := as.numeric(scale(elos)), by = Pen]
behaviourData[, scaleWins := as.numeric(scale(WinsY)), by = Pen]
behaviourData[, scaleLosses := as.numeric(scale(LossesY)), by = Pen]

# exclude unidentified animals
behaviourData = behaviourData[!(Pen == "A" & ID == "ZA")&
                                !(Pen == "B" & ID == "TU")&
                                !(Pen == "B" & ID == "VA")&
                                !(Pen == "C" & ID == "OL")&
                                !(Pen == "C" & ID == "ZZ")&
                                !(Pen == "C" & ID == "TX"),]

behaviourData[, LossesCluster2 := ifelse(scaleLosses < 0,"less losses than average Young", "more losses than average Young"), by = Pen]
behaviourData$LossesCluster2 = factor(behaviourData$LossesCluster2, levels = c("less losses than average Young", "more losses than average Young"))
behaviourData[, WinCluster2 := ifelse(scaleWins < 0,"less wins than average Young", "more wins than average Young"), by = Pen]
behaviourData$WinCluster2 = factor(behaviourData$WinCluster2, levels = c("less wins than average Young", "more wins than average Young"))

behaviourData[, LossesCluster := cut(scaleLosses,
                                     breaks=quantile(scaleLosses),
                                     include.lowest=TRUE,
                                     labels=c("<25%","25-50%","50-75%", ">75%"))]
behaviourData$LossesCluster = factor(behaviourData$LossesCluster, levels = c("<25%","25-50%","50-75%", ">75%"))
behaviourData[, WinCluster := cut(scaleWins,
                                  breaks=quantile(scaleWins),
                                  include.lowest=TRUE,
                                  labels=c("<25%","25-50%","50-75%", ">75%"))]
behaviourData$WinCluster = factor(behaviourData$WinCluster, levels = c("<25%","25-50%","50-75%", ">75%"))

behaviourData[, TopLoss := ifelse(scaleLosses>quantile(scaleLosses, 0.75), 1, 0)]
behaviourData[, TopAgress := ifelse(scaleWins>quantile(scaleWins, 0.75), 1, 0)]
behaviourData[, scaleWinsM := as.numeric(scale(Wins)), by = Pen]
behaviourData[, scaleLossesM := as.numeric(scale(Losses)), by = Pen]

behaviourData[Condition == "small", TopRankM := ifelse(rank<3, "TopRank", "Normal")]
behaviourData[Condition == "large", TopRankM := ifelse(rank<13, "TopRank", "Normal")]

behaviourData[Condition == "small", BottomRankM := ifelse(rank>17, "BottomRank", "Normal")]
behaviourData[Condition == "large", BottomRankM := ifelse(rank>107, "BottomRank", "Normal")]

behaviourData[, TopAgressM := ifelse(scaleWinsM>quantile(scaleWinsM, 0.75), 1, 0)]

fwrite(behaviourData, "IndividualData.csv", sep = ";")

# look at outliers
#Agress = rbind(outlier_animals(behaviourData[Pen == 'A']),
#               outlier_animals(behaviourData[Pen == 'B']),
#               outlier_animals(behaviourData[Pen == 'C']),
#               outlier_animals(behaviourData[Pen == 'D']),
#               outlier_animals(behaviourData[Pen == 'E']),
#               outlier_animals(behaviourData[Pen == 'F']))

#relationship Wins and Elos  
ggplot(behaviourData, aes(x = scaleLosses, y = scaleElos))+
  geom_abline(intercept = 0, slope = 0, linetype = "dashed")+
  #geom_point(
  #  data = behaviourData[TopAgress == 1], colour = "blue",
  #  size = 3
  #)+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x, col = "black")+
  facet_grid(.~WinCluster2)+
  labs(x = "scaled number of Losses Young", y = "scaled Elo rating Mature")+
  theme_classic(base_size = 18)


plotWins = behaviourData
plotWins_long = melt(plotWins, id.vars = c("Pen", "ID", "LossesCluster", "TopAgress", "TopRankM"), measure.vars = c("scaleWins", "scaleWinsM"))
plotWins_long[, fullID := paste0(Pen, ID)]

#plotLosses = behaviourData
#plotLosses_long = melt(plotLosses, id.vars = c("Pen", "ID", "LossesCluster", "TopAgress", "TopLoss"), measure.vars = c("scaleLosses", "scaleLossesM"))
#plotLosses_long[, fullID := paste0(Pen, ID)]


ggplot(plotWins_long, aes(x = variable, y = value)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group = fullID), colour = "darkgrey") + 
  geom_point(size = 2, colour = "darkgrey")+
  geom_point(
    data = plotWins_long[TopRankM == "TopRank",], colour = "#0072B2",
    size = 4
  )+
  geom_line(
    data = plotWins_long[TopRankM == "TopRank",],
    aes(group = fullID),colour = "#0072B2",
  )+
  geom_point(
    data = plotWins_long[TopAgress == 1,], 
    aes(colour = ">75%"), colour = "#D55E00",
    size = 3
  )+
  geom_line(
    data = plotWins_long[TopAgress == 1,],
    aes(colour = ">75%",group = fullID),
    colour = "#D55E00"
  )+
  scale_colour_discrete(labels=">75%")+
  facet_grid(.~Pen)+
  labs(y = "scaled number of Wins")+
  scale_x_discrete(labels = c("Young", "Mature"))+
  
  theme_bw(base_size = 18)+
  theme(axis.title.x=element_blank())


behaviourData[ ,diffWin := scaleWins -scaleWinsM]
behaviourData[ ,diffLoss := scaleLosses -scaleLossesM]

ggplot(behaviourData, aes(y = diffWin, x= LossesCluster)) + 
  geom_boxplot(outlier.shape = NA)+
  #geom_line(aes(group = fullID)) + 
  geom_point(size = 2)+
  theme_bw(base_size = 18)

ggplot(behaviourData, aes( x= scaleWins, y = scaleElos)) + 
  #geom_boxplot(outlier.shape = NA)+
  geom_point(size = 2)+
  geom_point(
    data = behaviourData[TopAgress == 1,], 
    aes(colour = ">75%"), colour = "#D55E00",
    size = 3
  )+
  geom_smooth()+
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
  #geom_boxplot(data = Top_agress[TopAgress == 1,], 
  #             colour = "orange",)+
  #geom_point(
  #  data = Top_agress[TopAgress == 1,], 
  #  colour = "orange",
  #  size = 4
  #)+
  facet_grid(.~as.factor(LossesCluster))+
  labs(y = "standardised Elo-ratings", x = "standardised Young Wins")+
  #scale_x_discrete(labels = c("Young", "Mature"))+
  theme_bw(base_size = 18)+
  theme(axis.title.x=element_blank())



# heatmap of wins, losses Young and Elo
ggplot(behaviourData, aes(WinCluster, LossesCluster)) +
  geom_tile(aes(fill = scaleElos), colour = "white") +
  scale_fill_gradient(low = "white", high = "steelblue")


ggplot(behaviourData, aes(WinCluster, y = scaleElos, fill = LossesCluster))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())

####### Model behaviour ###################
# plot for Elo distribution against normal in red
ggplot(behaviourData, aes(x = scaleElos)) +            
  geom_density(alpha = 0.1, position = 'identity')+
  stat_function(fun = dnorm, n = 420, args = list(mean = 0, sd = 1), col = 'red')+ 
  theme_classic(base_size = 18)
skewness(behaviourData$scaleElos) # -> positive skewed (right)
kurtosis(behaviourData$scaleElos) # -> nearly perfect kurtosis (normal = 3)
jarque.test(as.numeric(behaviourData$scaleElos)) # -> significantly different from normal but only slightly


model1 = lmer(scaleElos ~ scaleWins*scaleLosses*Condition + (1|Pen), behaviourData) 
#-> singularity issue -> leave out random effect
fullPredElo = lm(scaleElos ~ scaleWins*scaleLosses*Condition, behaviourData)
empty = lm(scaleElos ~ 1, behaviourData)
resid.Elo<- simulateResiduals(fullPredElo, 1000)
plot(resid.Elo)

options(na.action = "na.fail")
res = dredge(fullPredElo)
subset(res, delta <= 2, recalc.weights=FALSE)
summary(model.avg(res, revised.var=FALSE))
sw(res) # importance of predictors

reduce1 = lm(scaleElos ~ scaleWins*scaleLosses, behaviourData)
resid.Elo<- simulateResiduals(reduce1, 1000)
plot(resid.Elo)
summary(reduce1)
confint(reduce1)


library(sjPlot)
plot_model(reduce1, transform = NULL, show.values = TRUE, show.p = TRUE) +  
  ylim(-0.5, 0.5) +  
  theme_bw(base_size = 18)+
  geom_hline(yintercept= 0)

library(emmeans)
#contrasttrends = emtrends(reduce1, ~scaleWins:scaleLosses, var = "scaleWins", at = list(scaleLosses = c(-1,0,1)))
#contrasttrends2 = emtrends(reduce1, ~scaleWins:scaleLosses, var = "scaleLosses", at = list(scaleWins = c(-1,0,1)))
contrastmeans= emmeans(reduce1, ~scaleLosses:scaleWins, 
                       at = list(scaleWins = c(min(behaviourData$scaleWins),mean(behaviourData$scaleWins),max(behaviourData$scaleWins)), 
                                 scaleLosses = c(min(behaviourData$scaleLosses),mean(behaviourData$scaleLosses),max(behaviourData$scaleLosses))))

###
# plotting effect
library(plotly)
press_grid = seq(-3, 6, by = 0.1)
v_grid <- seq(-3, 6, by = 0.1)
newdat <- expand.grid(press_grid, v_grid)  #the grid results in the same values as the newdat in the OP
colnames(newdat) <- c("scaleWins", "scaleLosses")

pred <- predict.lm(reduce1, newdata = newdat, se=TRUE)

ymin <- pred$fit - 1.96 * pred$se.fit
ymax <- pred$fit + 1.96 * pred$se.fit
fit <- pred$fit 
z <- matrix(fit, length(press_grid))
ci.low <- matrix(ymin, length(press_grid))
ci.up <- matrix(ymax, length(press_grid))
ref = matrix(rep(0, length(ymin)), length(press_grid))

#3d plot
plot_ly(x = press_grid, y = v_grid) %>% 
  add_surface(z = z,
              colorscale = list(c(0,1),c("red","blue"))) %>% 
  add_surface(z = ci.low, opacity = 0.5, showscale = FALSE, colorscale = list(c(0,1),c("grey","grey"))) %>% 
  add_surface(z = ci.up, opacity = 0.5, showscale = FALSE, colorscale = list(c(0,1),c("grey","grey"))) %>%
  layout(zaxis = list(gridcolor = 'black'),
  #add_trace( z = ref, 
            type = "surface", mode = "lines",
            opacity = .8)
  #add_surface(z = ref, opacity = 0.5, showscale = FALSE, mode='lines')


detach("package:plotly", unload = TRUE)

#behaviour in numbers: 
behaviourData[, .N, by = .(TopAgress, Condition, TopRankM)]
behaviourData[, .N, by = .(TopAgress, Condition, BottomRankM)]

######################################################################
##### Method evaluation ####################

methodData= data.table(ID = c(ratingAM$Individuals$ID,
                               ratingBM$Individuals$ID,
                               ratingCM$Individuals$ID,
                               ratingDM$Individuals$ID,
                               ratingEM$Individuals$ID,
                               ratingFM$Individuals$ID),
                        rankM = c(ratingAM$Individuals$rank,
                                  ratingBM$Individuals$rank,
                                  ratingCM$Individuals$rank,
                                  ratingDM$Individuals$rank,
                                  ratingEM$Individuals$rank,
                                  ratingFM$Individuals$rank),
                       eloM = c(ratingAM$Individuals$elos,
                                 ratingBM$Individuals$elos,
                                 ratingCM$Individuals$elos,
                                 ratingDM$Individuals$elos,
                                 ratingEM$Individuals$elos,
                                 ratingFM$Individuals$elos),
                        rankAllorig = c(setNames(length(ratingAall$rating):1,
                                               names(ratingAall$rating))[sort(names(ratingAall$rating))],
                                        setNames(length(ratingBall$rating):1,
                                                 names(ratingBall$rating))[sort(names(ratingBall$rating))],
                                        setNames(length(ratingCall$rating):1,
                                                 names(ratingCall$rating))[sort(names(ratingCall$rating))],
                                        setNames(length(ratingDall$rating):1,
                                                 names(ratingDall$rating))[sort(names(ratingDall$rating))],
                                        setNames(length(ratingEall$rating):1,
                                                 names(ratingEall$rating))[sort(names(ratingEall$rating))],
                                        setNames(length(ratingFall$rating):1,
                                                 names(ratingFall$rating))[sort(names(ratingFall$rating))]),
                        rankAllrand = c(setNames(length(ratingAall$rating_rand):1,
                                                 names(ratingAall$rating_rand))[sort(names(ratingAall$rating_rand))],
                                        setNames(length(ratingBall$rating_rand):1,
                                                 names(ratingBall$rating_rand))[sort(names(ratingBall$rating_rand))],
                                        setNames(length(ratingCall$rating_rand):1,
                                                 names(ratingCall$rating_rand))[sort(names(ratingCall$rating_rand))],
                                        setNames(length(ratingDall$rating_rand):1,
                                                 names(ratingDall$rating_rand))[sort(names(ratingDall$rating_rand))],
                                        setNames(length(ratingEall$rating_rand):1,
                                                 names(ratingEall$rating_rand))[sort(names(ratingEall$rating_rand))],
                                        setNames(length(ratingFall$rating_rand):1,
                                                 names(ratingFall$rating_rand))[sort(names(ratingFall$rating_rand))]),
                       eloAllorig = c(ratingAall$rating[sort(names(ratingAall$rating))],
                                      ratingBall$rating[sort(names(ratingBall$rating))],
                                      ratingCall$rating[sort(names(ratingCall$rating))],
                                      ratingDall$rating[sort(names(ratingDall$rating))],
                                      ratingEall$rating[sort(names(ratingEall$rating))],
                                      ratingFall$rating[sort(names(ratingFall$rating))]),
                       eloAllrand = c(ratingAall$rating_rand[sort(names(ratingAall$rating_rand))],
                                      ratingBall$rating_rand[sort(names(ratingBall$rating_rand))],
                                      ratingCall$rating_rand[sort(names(ratingCall$rating_rand))],
                                      ratingDall$rating_rand[sort(names(ratingDall$rating_rand))],
                                      ratingEall$rating_rand[sort(names(ratingEall$rating_rand))],
                                      ratingFall$rating_rand[sort(names(ratingFall$rating_rand))]))

methodData[, Pen := factor(c(rep('A', length(ratingAM$Individuals$rank)),
                                rep('B',length(ratingBM$Individuals$rank)),
                                rep('C',length(ratingCM$Individuals$rank)),
                                rep('D', 20), rep('E',20), rep('F',20)))]

methodData[, Condition := factor(c(rep('large', length(ratingAM$Individuals$rank)+
                                                         length(ratingBM$Individuals$rank)+
                                                         length(ratingCM$Individuals$rank)),
                                    rep('small', 20*3)))]

methodData = methodData[!(Pen == "A" & ID == "ZA")&
                        !(Pen == "B" & ID == "TU")&
                                !(Pen == "B" & ID == "VA")&
                                !(Pen == "C" & ID == "OL")&
                                !(Pen == "C" & ID == "ZZ")&
                                !(Pen == "C" & ID == "TX"),]
methodData = methodData[!is.na(rankM),]

methodData[, diffRank := abs(rankM - rankAllorig),]
# Spearman rank correlation between mature Elo and all Elo rand & orig

compResults = methodData[, .(corOrig = cor.test(rankM,  rankAllorig, method = "spearman")$estimate,
               corOrigp =round(cor.test(rankM, rankAllorig, method = "spearman")$p.value, 5),
               corRand = cor.test(rankM,  rankAllrand, method = "spearman")$estimate,
               corRandp =round(cor.test(rankM, rankAllrand, method = "spearman")$p.value, 5),
               meanRankDiffOrig = mean(abs(rankM-rankAllorig)),
               sdRankDiffOrig = sd(abs(rankM-rankAllorig)),
               lCIRankDiffOrig = mean(abs(rankM-rankAllorig))- 
                 qt(0.975,df=length(rankM)-1)*sd(abs(rankM-rankAllorig))/sqrt(length(rankM)),
               hCIRankDiffOrig = mean(abs(rankM-rankAllorig))+ 
                 qt(0.975,df=length(rankM)-1)*sd(abs(rankM-rankAllorig))/sqrt(length(rankM)),
               maxRankDiffOrig = max(abs(rankM-rankAllorig)),
               minRankDiffOrig = min(abs(rankM-rankAllorig)),
               meanRankDiffRand = mean(abs(rankM-rankAllrand)),
               sdRankDiffRand = sd(abs(rankM-rankAllrand)),
               maxRankDiffRand = max(abs(rankM-rankAllrand)),
               minRankDiffRand = min(abs(rankM-rankAllrand)),
               iccOrig = icc(data.table(eloM, eloAllorig))$value,
               iccRand = icc(data.table(eloM, eloAllrand))$value
               ),by = Pen]

compResults2 = methodData[, .(
                             meanRankDiffOrig = mean(abs(rankM-rankAllorig)),
                             sdRankDiffOrig = sd(abs(rankM-rankAllorig)),
                             lCIRankDiffOrig = mean(abs(rankM-rankAllorig))- 
                               qt(0.975,df=length(rankM)-1)*sd(abs(rankM-rankAllorig))/sqrt(length(rankM)),
                             hCIRankDiffOrig = mean(abs(rankM-rankAllorig))+ 
                               qt(0.975,df=length(rankM)-1)*sd(abs(rankM-rankAllorig))/sqrt(length(rankM)),
                             maxRankDiffOrig = max(abs(rankM-rankAllorig)),
                             minRankDiffOrig = min(abs(rankM-rankAllorig)),
                             meanRankDiffRand = mean(abs(rankM-rankAllrand)),
                             sdRankDiffRand = sd(abs(rankM-rankAllrand)),
                             maxRankDiffRand = max(abs(rankM-rankAllrand)),
                             minRankDiffRand = min(abs(rankM-rankAllrand)),
                             iccOrig = icc(data.table(eloM, eloAllorig))$value,
                             iccRand = icc(data.table(eloM, eloAllrand))$value
),by = Condition]


methodData[Pen == "A"| Pen == "B"| Pen == "C", mean(diffRank)]
methodData[Pen == "A"| Pen == "B"| Pen == "C", quantile(diffRank, c(0.05, 0.75, 0.95))]

methodData[Pen == "D"| Pen == "E"| Pen == "F", mean(diffRank)]
methodData[Pen == "D"| Pen == "E"| Pen == "F", quantile(diffRank, c(0.05, 0.75,0.95))]



######################################################################
###### logistic regression examples of growth rate ###############
library("viridis")

#TODO: change to new formula!!

plotdata = data.table(RankDiff = 0:19)
plotdata[, `0.05` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/20)*0.05*x))})]
plotdata[, `0.1` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/20)*0.1*x))})]
plotdata[, `0.2` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/20)*0.2*x))})]
plotdata[, `0.3` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/20)*0.3*x))})]
plotdata[, `0.5` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/20)*0.5*x))})]
plotdata[, `1` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/20)*1*x))})]
plotdata[, `2` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/20)*2*x))})]
plotdata[, `15` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/20)*15*x))})]

plotdata_long = melt(plotdata, id.vars = "RankDiff", value.name = "Probab")

 f1 = ggplot(plotdata_long, aes(x = RankDiff, y = Probab, color = variable))+
  geom_line( size =1.5)+
  labs(colour = expression(paste("Steepness ",italic("k"))))+
   scale_color_viridis(discrete = TRUE)+
  theme_bw(base_size = 18)+
   theme(axis.title.x = element_blank(),
         axis.title.y = element_blank())

plotdata2 = data.table(RankDiff = 1:119)
plotdata2[, `0.05` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/120)*0.05*x))})]
plotdata2[, `0.1` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/120)*0.1*x))})]
plotdata2[, `0.2` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/120)*0.2*x))})]
plotdata2[, `0.3` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/120)*0.3*x))})]
plotdata2[, `0.5` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/120)*0.5*x))})]
plotdata2[, `1` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/120)*1*x))})]
plotdata2[, `2` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/120)*2*x))})]
plotdata2[, `15` := sapply(RankDiff, function(x){1/(1+exp(-1*(10/120)*15*x))})]

plotdata2_long = melt(plotdata2, id.vars = "RankDiff", value.name = "Probab")

f2 = ggplot(plotdata2_long, aes(x = RankDiff, y = Probab, color = variable))+
  geom_line( size =1.5)+
  labs(colour = expression(paste("Steepness ",italic("k"))))+
  scale_color_viridis(discrete = TRUE)+
  theme_bw(base_size = 18)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

###### hierarchy example pen E ##########
library(ggpubr)
p1 = ratingEY$hierarchy_shape$plot+ annotate(geom = 'label', x=17, y=0.27, label="Young", size = 6)+scale_y_continuous(limits=c(0.25,1),oob = rescale_none)
p2 = ratingEM$hierarchy_shape$plot+ annotate(geom = 'label', x=17, y=0.27, label="Mature", size = 6)+scale_y_continuous(limits=c(0.25,1),oob = rescale_none)
p3 = ratingEall$hierarchy_shape$plot+ annotate(geom = 'label', x=17, y=0.27, label="Full set", size = 6)+scale_y_continuous(limits=c(0.25,1),oob = rescale_none)


fig1 = ggarrange(p1,p2,p3,ncol = 3, labels = c("A", "B", "C"), 
                 common.legend = TRUE, legend = "right")
annotate_figure(fig1,
                bottom = text_grob("Difference in rank",
                                   size = 18),
                left = text_grob("Probability that higher rank wins", rot = 90, size = 18))
                
library("viridis") 
figSupp = ggarrange(f1 +  annotate(geom = 'label', x=15, y=0.37, label="Group Size = 20", size = 5),
                    f2 + annotate(geom = 'label', x=100, y=0.37, label="Group Size = 120", size = 5), 
                    labels = c("A", "B"), common.legend = TRUE, legend = "right")
annotate_figure(figSupp,
                bottom = text_grob("Difference in rank",
                                   size = 18),
                left = text_grob("Probability that higher rank wins", rot = 90, size = 18))




###############################################################
### PAPER 2 
###### DYNAMICS OF INTERACTIONS #################
## done on mature data


Interact = rbind(ratingAM$Individuals, ratingBM$Individuals, ratingCM$Individuals, ratingDM$Individuals,
                 ratingEM$Individuals, ratingFM$Individuals)
Interact[, Pen := factor(c(rep('A', length(ratingAM$Individuals$rank)), 
                           rep('B',length(ratingBM$Individuals$rank)), 
                           rep('C',length(ratingCM$Individuals$rank)),
                           rep('D', 20), rep('E',20), rep('F',20)))]

Interact[, Condition := factor(c(rep("large", length(Pen)-60), 
                                 rep("small", 60)))]
Interact = na.omit(Interact)

Interact[, ratio := sum/(sum(sum)*0.5), by = Pen]
Interact[order(rank), cumRatio := cumsum(ratio), by = Pen]
Interact[order(-rank), cumRatioLow := cumsum(ratio), by = Pen]
Interact[order(-rank), which(cumRatioLow > 0.5)[1], by = Pen]

ggplot(Interact[Condition == 'small',],aes(x = rank, y = ratio, colour = Pen)) + geom_point()+geom_smooth(se = F)
ggplot(Interact[Condition == 'large',],aes(x = rank, y = ratio, colour = Pen)) + geom_point()+geom_smooth(se = F)

#test = glmer(Interactions ~ poly(Elo,) + (1|Pen), plotdata1, family = 'poisson')
#library(DHARMa)
#resid.df2<- simulateResiduals(test, 1000)
#plot(resid.df2, asFactor = T)

# Plot for Dynamics of interactions by rank (small)
ggplot(data = Interact[Condition == 'small',], mapping = aes(x = rank, y =sum, colour = Pen)) + 
  geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(y = 'Number of Interactions')+
  theme_classic(base_size = 18)

ggplot(data = Interact[Condition == 'small',], mapping = aes(x = rank, y =CumuSum, colour = Pen)) + 
  geom_line()+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(y = 'Number of Interactions')+
  theme_classic(base_size = 18)

# Plot for Dynamics of interactions by rank (large)
ggplot(data = Interact[Condition == 'large',], mapping = aes(x = rank, y = sum, colour = Pen)) + 
  geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(y = 'Number of Interactions')+
  theme_classic(base_size = 18)
###

InteractionSmall = data.table(LoserRank = c(ratingDM$rankMatch$LoserRank, 
                                     ratingEM$rankMatch$LoserRank, 
                                     ratingFM$rankMatch$LoserRank),
                       WinnerRank = c(ratingDM$rankMatch$WinnerRank, 
                                      ratingEM$rankMatch$WinnerRank, 
                                      ratingFM$rankMatch$WinnerRank),
                       Pen = factor(c(rep('D', length(ratingDM$rankMatch$LoserRank)), 
                                      rep('E', length(ratingEM$rankMatch$LoserRank)), 
                                      rep('F', length(ratingFM$rankMatch$LoserRank)))))
InteractionSmall[ ,Condition := "small"]
InteractionLarge = data.table(LoserRank = c(ratingAM$rankMatch$LoserRank, 
                                     ratingBM$rankMatch$LoserRank, 
                                     ratingCM$rankMatch$LoserRank),
                       WinnerRank = c(ratingAM$rankMatch$WinnerRank, 
                                      ratingBM$rankMatch$WinnerRank, 
                                      ratingCM$rankMatch$WinnerRank),
                       Pen = factor(c(rep('A', length(ratingAM$rankMatch$WinnerRank)), 
                                      rep('B',length(ratingBM$rankMatch$WinnerRank)), 
                                      rep('C',length(ratingCM$rankMatch$WinnerRank)))))

InteractionLarge[, Condition := "large"]
allInteractions =rbind(InteractionSmall, InteractionLarge)

allInteractions[, RankDiff := abs(WinnerRank-LoserRank)]
allInteractions[, HighRankWins := WinnerRank < LoserRank]



# Plot for Dynamics of interactions by rank (small)
ggplot(data = InteractionSmall, mapping = aes(x = WinnerRank, y =LoserRank)) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  labs(y = 'Loser rank', x= 'Winner rank')+
  theme_classic(base_size = 18)+
  facet_grid(. ~ Pen)+
  #stat_density_2d(aes(fill = ..level..), geom="polygon")
  #geom_density_2d(aes(colour = Pen), size = 2)
  #geom_density_2d_filled(contour_var = "ndensity")
  #geom_density_2d_filled(contour_var = "count") 
  geom_density_2d_filled(alpha = 0.7, contour_var = "count")+
  geom_abline(intercept = 0 , slope = 1, linetype = "dashed", colour = 'grey')+
  geom_jitter(size = 1)+
  xlim(1,20)+
  ylim(1,20)




# Plot for Dynamics of interactions by rank (large)
ggplot(data = InteractionLarge, mapping = aes(x = WinnerRank, y =LoserRank)) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  #geom_point(size = 2.5) + 
  labs(y = 'Loser rank', x= 'Winner rank')+
  theme_classic(base_size = 18)+
  facet_grid(. ~ Pen)+
  #stat_density_2d(aes(fill = ..level..), geom="polygon")
  #geom_density_2d(aes(colour = Pen), size = 2)
  #geom_density_2d_filled(contour_var = "ndensity")
  #geom_density_2d_filled(contour_var = "count") 
  geom_density_2d_filled(alpha = 0.7, contour_var = "count")+
  geom_abline(intercept = 0 , slope = 1, linetype = "dashed", colour = 'grey')+
  geom_jitter(size = 0.6, colour = 'black')+
  xlim(1,120)+
  ylim(1,120)

fwrite(rbind(InteractionSmall, InteractionLarge), file = "InteractionsRank.csv", sep = ";")


# model propability to win according to difference of rank
Win_data <- allInteractions %>%
  #code buzz cola choice as a binary variable
  mutate(Winner_high = case_when(
    choice == "buzz_cola" ~ 1,
    choice == "slurm" ~ 0
  )) %>%
  #group by combinations and find the proportion of buzz cola choices
  group_by(buzz_cola, slurm) %>%
  summarise(fraction_choose_cola = mean(buzz_cola_choice))

WinProb = allInteractions[, .(highRankWinProb = mean(HighRankWins), Condition = Condition), by = c("Pen", "RankDiff")]

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

#plot the logistic regression on the entire choice data
  ggplot(data = WinProb[Condition == "small"], aes(x = RankDiff, y = highRankWinProb, colour = factor(Pen))) +
  geom_point() +
  binomial_smooth(se = FALSE) +
  #add in some aesthetics
  scale_colour_discrete(name = "Pen") +
  labs(x = "RankDiff",
       y = "Porbability that high ranking individual wins") +
  theme_minimal()

WinProb[RankDiff == 2,]  
  
WinChance.model = glmer(HighRankWins ~ RankDiff + (1|Pen), data = allInteractions[Condition == "small",], family = binomial)
resid.Win<- simulateResiduals(WinChance.model, 100)
plot(resid.Win)
summary(WinChance.model)

##### clustering ############



Interact_clean = Interact[!is.na(rank),]

#PCA on wins losses and elo
interact.pca = prcomp(Interact_clean[,c(2,3, 4)], center = T, scale = TRUE)

summary(interact.pca)
str(interact.pca)

autoplot(interact.pca, data =Interact_clean, colour = 'Pen', size = 3, loadings = T,  loadings.label=T,scale =0) #to see individuals use label= T


set.seed(1)
interact.kmeans = kmeans(Interact_clean[,c(2,3, 4)], centers = 3)
plot(interact.kmeans, data = Interact_clean, size= 3)

Interact_clean[, cluster := as.factor(interact.kmeans$cluster)]
Interact_clean = Interact_clean[order(Pen,rank),]


ggplot(data = Interact_clean, mapping = aes(x = cluster, y = rank, fill = Condition)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), alpha= 0.4, size = 2)+
  theme_classic(base_size = 18)



descdist(Interact_clean$sum, discrete = FALSE)
poisson = fitdist(Interact_clean$sum, 'poisson')
gamma = fitdist(Interact_clean$sum, 'gamma')
plot(gamma)


modelInteract = glmer.nb(sum ~ cluster*Condition + (1|Pen), Interact_clean)

#IDEA: do on exact scaled ELo-rating instead of rank or cluster
# vergleich cluster und Elo-rating (verteilung, skalieren)
modelInteract2 = glmer.nb(sum ~ poly(scaleElos,3)*Condition + (1|Pen), Interact_scale)
# how to interpret this????
resid.df1<- simulateResiduals(modelInteract, 500)
resid.df2<- simulateResiduals(modelInteract2, 500)
plot(resid.df1)
plot(resid.df2)
summary(modelInteract2)
Interact_clean$CC <- interaction(Interact_clean$cluster, Interact_clean$Condition)
modelInteract_CC = glmer.nb(sum ~ CC + (1|Pen), Interact_clean)
#estimates r?ckrechnen e^estimate 

summary(glht(modelInteract_CC,linfct=mcp(CC="Tukey")))


Interact_scale = Interact_clean
Interact_scale[,scaleWins := scale(Wins), by = Pen]
Interact_scale[,scaleLosses := scale(Losses), by = Pen]
Interact_scale[,scaleElos := scale(elos), by =  Pen]

sample <- sample(c(TRUE, FALSE), nrow(Interact_scale), replace=TRUE, prob=c(0.7,0.3))
train <- data.table(Interact_scale[sample, ])
test <- data.table(Interact_scale[!sample, ]) 
lda_model <- lda(cluster~Wins + Losses +elos, data=train)
predicted <- predict(lda_model, test)

mean(predict(lda_model, train)$class == train$cluster)
mean(predicted$class==test$cluster)
ggord(lda_model, train$cluster)

#reducing


######
# survival on 

#model = survfit(Surv(sum(Wins), sum(Losses))~RankDiff, data=DATA)
# summary(model)
# plot(model)
#coxph(Surv(Wins,Losses)~RankDiff, data) how likely that with higher rank you have a higher chance to win
