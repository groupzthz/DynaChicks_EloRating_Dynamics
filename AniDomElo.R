getwd()

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
Individ <- fread("Individuals.csv",header = TRUE, sep = ";")
dataMature <- fread("MatureInteractionsFull.csv",header = TRUE, sep = ";")

#check Data
str(dataMature)
length(unique(dataMature$Date)) == 4
unique(dataMature$Pen)

#convert Date
dataMature$Date <- as.Date(as.character(dataMature$Date), format = "%d.%m.%Y")
dataMature$VideoTime <- as.POSIXct(as.character(dataMature$VideoTime), format = "%H:%M:%S")
#check correct order
dataMature = dataMature[order(Date, Pen, VideoTime),]

#check that data is complete
dataMature[, .(minTime = min(VideoTime), maxTime = max(VideoTime)), by = .(Date,Pen)]

#check order of data

#split into pens
PenAt = dataMature[Pen == "A", c(1, 3, 4,5, 7)]
PenBt = dataMature[Pen == "B", c(1, 3, 4,5, 7)]
PenCt = dataMature[Pen == "C", c(1, 3, 4,5, 7)]
PenDt = dataMature[Pen == "D", c(1, 3, 4,5, 7)]
PenEt = dataMature[Pen == "E", c(1, 3, 4,5, 7)]
PenFt = dataMature[Pen == "F", c(1, 3, 4,5, 7)]

#check entered Individuals

PenA = PenAt[(Winner %in% unique(Individ$ID[Individ$Pen == "A"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "A"])),]
PenB = PenBt[(Winner %in% unique(Individ$ID[Individ$Pen == "B"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "B"])),]
PenC = PenCt[(Winner %in% unique(Individ$ID[Individ$Pen == "C"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "C"])),]
PenD = PenDt[(Winner %in% unique(Individ$ID[Individ$Pen == "D"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "D"])),]
PenE = PenEt[(Winner %in% unique(Individ$ID[Individ$Pen == "E"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "E"])),]
PenF = PenFt[(Winner %in% unique(Individ$ID[Individ$Pen == "F"]))& 
                 (Loser %in% unique(Individ$ID[Individ$Pen == "F"])),]

#exclusion of double tagged animal during Mature
PenA = PenA[Winner != "ZA" & Loser != "ZA"]


rm(dataMature, 
   PenAt, PenBt, PenCt, PenDt, PenEt, PenFt)

################ Diagnostics ###############################################################

# Dataset analytics:

diagnA = dataset_diagnostics(PenA, Individ$ID[Individ$Pen == 'A'])
diagnB = dataset_diagnostics(PenB, Individ$ID[Individ$Pen == 'B'])
diagnC = dataset_diagnostics(PenC, Individ$ID[Individ$Pen == 'C'])
diagnD = dataset_diagnostics(PenD, Individ$ID[Individ$Pen == 'D'])
diagnE = dataset_diagnostics(PenE, Individ$ID[Individ$Pen == 'E'])
diagnF = dataset_diagnostics(PenF, Individ$ID[Individ$Pen == 'F'])

##################### Elo rating ################################################################
# RATING
ratingA = elo_analysis(PenA, Individ$ID[Individ$Pen == 'A']) 
ratingB = elo_analysis(PenB, Individ$ID[Individ$Pen == 'B']) 
ratingC = elo_analysis(PenC, Individ$ID[Individ$Pen == 'C']) 
ratingD = elo_analysis(PenD, Individ$ID[Individ$Pen == 'D']) 
ratingE = elo_analysis(PenE, Individ$ID[Individ$Pen == 'E'])
ratingF = elo_analysis(PenF, Individ$ID[Individ$Pen == 'F'])

printCalculations(ratingA)
printCalculations(ratingB)
printCalculations(ratingC)
printCalculations(ratingD)
printCalculations(ratingE)
printCalculations(ratingF)

#######################################################################################

####### DATASET DESCRIPTIONS ################

PenA[, .N, by = Condition]
PenB[, .N, by = Condition]
PenC[, .N, by = Condition]
PenD[, .N, by = Condition]
PenE[, .N, by = Condition]
PenF[, .N, by = Condition]

PenA[, .N, by = Code]
PenB[, .N, by = Code]
PenC[, .N, by = Code]
PenD[, .N, by = Code]
PenE[, .N, by = Code]
PenF[, .N, by = Code]

###### SUM OF INTERACTIONS #################

Interact = rbind(ratingA$Individuals, ratingB$Individuals, ratingC$Individuals, ratingD$Individuals,
                 ratingE$Individuals, ratingF$Individuals)
Interact[, Pen := factor(c(rep('A', length(ratingA$Individuals$rank)), 
                           rep('B',length(ratingB$Individuals$rank)), 
                           rep('C',length(ratingC$Individuals$rank)),
                           rep('D', 20), rep('E',20), rep('F',20)))]

Interact[, Condition := factor(c(rep("large", length(Pen)-60), 
                                 rep("small", 60)))]
Interact = na.omit(Interact)

# ratio of Inetraction per Individual
Interact[, ratio := sum/(sum(sum)*0.5), by = Pen]

#scaled elos
Interact[, scaleElo := scale(elos), by = Pen]

Interact[order(rank), cumRatio := cumsum(ratio), by = Pen]
Interact[order(-rank), cumRatioLow := cumsum(ratio), by = Pen]
Interact[order(-rank), which(cumRatioLow > 0.5)[1], by = Pen]

#plot ratio of interactions by individual
ggplot(Interact[Condition == 'small',],aes(x = rank, y = ratio, colour = Pen)) + 
  geom_point()+
  geom_smooth(se = F)
ggplot(Interact[Condition == 'large',],aes(x = rank, y = ratio, colour = Pen)) + 
  geom_point()+
  geom_smooth(se = F)

# test = glmer(sum ~ poly(scaleElo,2)*Condition + (1|Pen), Interact, family = 'poisson')
# resid.df2<- simulateResiduals(test, 1000)
# plot(resid.df2)

# Plot for Dynamics of interactions by rank (small)
ggplot(data = Interact[Condition == 'small',], mapping = aes(x = rank, y =sum, colour = Pen)) + 
  geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(x = 'Rank', y = 'Number of Interactions')+
  theme_classic(base_size = 18)


# Plot for Dynamics of interactions by rank (large)
ggplot(data = Interact[Condition == 'large',], mapping = aes(x = rank, y = sum, colour = Pen)) + 
  geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(x = 'Rank', y = 'Number of Interactions')+
  theme_classic(base_size = 18)


#TODO: test this by equalising intervals (min-max scaling) and then Kolgomorv to test if distributions are the same?

#### INTERACTIONS BY RANK ############


RankTable = rbind(ratingA$rankMatch, ratingB$rankMatch, ratingC$rankMatch, 
          ratingD$rankMatch, ratingE$rankMatch, ratingF$rankMatch)

RankTable[, Pen := factor(c(rep('A', length(ratingA$rankMatch$WinnerRank)), 
                            rep('B',length(ratingB$rankMatch$WinnerRank)), 
                            rep('C',length(ratingC$rankMatch$WinnerRank)),
                                rep('D', length(ratingD$rankMatch$LoserRank)), 
                                      rep('E', length(ratingE$rankMatch$LoserRank)), 
                                      rep('F', length(ratingF$rankMatch$LoserRank))))]
RankTable[Pen == "A" | Pen == "B"| Pen == "C" ,Condition := "large"]
RankTable[Pen == "D" | Pen == "E"| Pen == "F" ,Condition := "small"]

RankTable[, RankDiff := abs(WinnerRank-LoserRank)]
RankTable[, HighRankWins := WinnerRank < LoserRank]

RankTable[Code == "Avoidance"| Code == "Threat", AggressLvl := "non_physical"]
RankTable[Code == "Peck"| Code == "Fight", AggressLvl := "physical"]


# Plot for Dynamics of interactions by rank (small)
ggplot(data = RankTable[Condition == "small",], mapping = aes(x = WinnerRank, y =LoserRank)) + 
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

ggplot(data = RankTable[Condition == "small",], mapping = aes(x = WinnerRank, y =LoserRank)) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  labs(y = 'Loser rank', x= 'Winner rank')+
  theme_classic(base_size = 18)+
  facet_grid(AggressLvl ~ Pen)+
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
ggplot(data = RankTable[Condition == "large",], mapping = aes(x = WinnerRank, y =LoserRank)) + 
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

# facet by aggression intensity
ggplot(data = RankTable[Condition == "large",], mapping = aes(x = WinnerRank, y =LoserRank)) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  #geom_point(size = 2.5) + 
  labs(y = 'Loser rank', x= 'Winner rank')+
  theme_classic(base_size = 18)+
  facet_grid(AggressLvl ~ Pen)+
  #stat_density_2d(aes(fill = ..level..), geom="polygon")
  #geom_density_2d(aes(colour = Pen), size = 2)
  #geom_density_2d_filled(contour_var = "ndensity")
  #geom_density_2d_filled(contour_var = "count") 
  geom_density_2d_filled(alpha = 0.7, contour_var = "count")+
  geom_abline(intercept = 0 , slope = 1, linetype = "dashed", colour = 'grey')+
  geom_jitter(size = 0.6, colour = 'black')+
  xlim(1,120)+
  ylim(1,120)


ggplot(data = RankTable[Condition == "small",], aes(x = AggressLvl, y = RankDiff))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  facet_grid(.~Pen)
ggplot(data = RankTable[Condition == "large",], aes(x = AggressLvl, y = RankDiff))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  facet_grid(.~Pen)


#TODO: how to check if dstributions are similar?
#Aggression type by rank difference
ggplot(data = RankTable[Condition == "small",], aes(x = RankDiff, color = AggressLvl))+
  geom_density(size = 2)+
  #facet_grid(.~Pen)+
  theme_bw(base_size = 18)
qqplot(RankTable[Condition == "small" & AggressLvl == "non_physical",RankDiff],
       RankTable[Condition == "small" & AggressLvl == "physical",RankDiff])
abline(a = 0, b = 1, lty = 3)
ks.test(RankTable[Condition == "small" & AggressLvl == "non_physical",RankDiff],
        RankTable[Condition == "small" & AggressLvl == "physical",RankDiff])

ggplot(data = RankTable[Condition == "large",], aes(x = RankDiff, color = AggressLvl))+
  geom_density(size = 2)+
  #facet_grid(.~Pen)+
  theme_bw(base_size = 18)
qqplot(RankTable[Condition == "large" & AggressLvl == "non_physical",RankDiff],
       RankTable[Condition == "large" & AggressLvl == "physical",RankDiff])
abline(a = 0, b = 1, lty = 3)
ks.test(RankTable[Condition == "large" & AggressLvl == "non_physical",RankDiff],
        RankTable[Condition == "large" & AggressLvl == "physical",RankDiff])

# Aggression type by individual rank
# by Winner rank
ggplot(data = RankTable[Condition == "large",], aes(x = WinnerRank, color = AggressLvl))+
  geom_density(size = 2)+
  #facet_grid(.~Pen)+
  theme_bw(base_size = 18)
qqplot(RankTable[Condition == "large" & AggressLvl == "non_physical",WinnerRank],
       RankTable[Condition == "large" & AggressLvl == "physical",WinnerRank])
abline(a = 0, b = 1, lty = 3)
ks.test(RankTable[Condition == "large" & AggressLvl == "non_physical",WinnerRank],
        RankTable[Condition == "large" & AggressLvl == "physical",WinnerRank])

testData = 


test = glmer()


ggplot(data = RankTable[Condition == "large",], aes(x = LoserRank, color = AggressLvl))+
  geom_density(size = 2)+
  #facet_grid(.~Pen)+
  theme_bw(base_size = 18)
qqplot(RankTable[Condition == "large" & AggressLvl == "non_physical",LoserRank],
       RankTable[Condition == "large" & AggressLvl == "physical",LoserRank])
abline(a = 0, b = 1, lty = 3)
ks.test(RankTable[Condition == "large" & AggressLvl == "non_physical",LoserRank],
        RankTable[Condition == "large" & AggressLvl == "physical",LoserRank])


ggplot(data = RankTable[Condition == "small",], aes(x = WinnerRank, color = AggressLvl))+
  geom_density(size = 2)+
  #facet_grid(.~Pen)+
  theme_bw(base_size = 18)
qqplot(RankTable[Condition == "small" & AggressLvl == "non_physical",WinnerRank],
       RankTable[Condition == "small" & AggressLvl == "physical",WinnerRank])
abline(a = 0, b = 1, lty = 3)
ks.test(RankTable[Condition == "small" & AggressLvl == "non_physical",WinnerRank],
        RankTable[Condition == "small" & AggressLvl == "physical",WinnerRank])

ggplot(data = RankTable[Condition == "small",], aes(x = LoserRank, color = AggressLvl))+
  geom_density(size = 2)+
  #facet_grid(.~Pen)+
  theme_bw(base_size = 18)
qqplot(RankTable[Condition == "small" & AggressLvl == "non_physical",LoserRank],
       RankTable[Condition == "small" & AggressLvl == "physical",LoserRank])
abline(a = 0, b = 1, lty = 3)
ks.test(RankTable[Condition == "small" & AggressLvl == "non_physical",LoserRank],
        RankTable[Condition == "small" & AggressLvl == "physical",LoserRank])


#fwrite(rbind(InteractionSmall, InteractionLarge), file = "InteractionsRank.csv", sep = ";")

######## AGGRESSION BY CONDITION ##############
ggplot(data = RankTable[Condition == "small",], aes(x = RankDiff, color = AggressLvl))+
  geom_density(size = 2)+
  facet_grid(.~Situation)+
  theme_bw(base_size = 18)

ggplot(data = RankTable[Condition == "large",], aes(x = RankDiff, color = AggressLvl))+
  geom_density(size = 2)+
  facet_grid(.~Situation)+
  theme_bw(base_size = 18)


ggplot(data = RankTable[Condition == "small",], aes(x = LoserRank, color = AggressLvl))+
  geom_density(size = 2)+
  facet_grid(Pen~Situation)+
  theme_bw(base_size = 18)



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
