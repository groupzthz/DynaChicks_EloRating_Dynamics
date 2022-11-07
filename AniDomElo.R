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
library(parameters) # model parameters
library(emmeans) #model means
library(effects) #model effects visualisation
library(RColorBrewer) # color for plotting
library(EloSteepness) # for steepness measure
library(sjPlot)
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

##################### Bayesian Steepness ##########################################################
#Steepness
steepA = elo_steepness_from_sequence(winner = PenA$Winner,
                                      loser = PenA$Loser,
                                      refresh = 0,
                                      cores = 2)
steepB = elo_steepness_from_sequence(winner = PenB$Winner,
                                     loser = PenB$Loser,
                                     refresh = 0,
                                     cores = 2)

steepC = elo_steepness_from_sequence(winner = PenC$Winner,
                                     loser = PenC$Loser,
                                     refresh = 0,
                                     cores = 2)
steepD = elo_steepness_from_sequence(winner = PenD$Winner,
                                     loser = PenD$Loser,
                                     refresh = 0,
                                     cores = 2)
steepE = elo_steepness_from_sequence(winner = PenE$Winner,
                                     loser = PenE$Loser,
                                     refresh = 0,
                                     cores = 2)
steepF = elo_steepness_from_sequence(winner = PenF$Winner,
                                     loser = PenF$Loser,
                                     refresh = 0,
                                     cores = 2)


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

###### DATA Tables to work with #################

# Overview of data of all individuals 
Interact = rbind(ratingA$Individuals, ratingB$Individuals, ratingC$Individuals, ratingD$Individuals,
                 ratingE$Individuals, ratingF$Individuals)
Interact[, Pen := factor(c(rep('A', length(ratingA$Individuals$rank)), 
                           rep('B',length(ratingB$Individuals$rank)), 
                           rep('C',length(ratingC$Individuals$rank)),
                           rep('D', 20), rep('E',20), rep('F',20)))]

Interact[, Group_Size := factor(c(rep("large", length(Pen)-60), 
                                 rep("small", 60)))]


# ratio of Interaction per Individual
Interact[, ratio := sum/(sum(sum)*0.5), by = Pen]
Interact[, ratioIntens := physAggr/(physAggr+nonphysAggr)]


#fwrite(Interact, "IndividualData.csv", sep = ";")

Interact = Interact[!is.na(elos),]

Interact[, rowNum := 1:.N]

####

#All single Interactions with updated ranks and elos
RankTable = rbind(ratingA$rankMatch, ratingB$rankMatch, ratingC$rankMatch, 
                  ratingD$rankMatch, ratingE$rankMatch, ratingF$rankMatch)

RankTable[, Pen := factor(c(rep('A', length(ratingA$rankMatch$WinnerRank)), 
                            rep('B',length(ratingB$rankMatch$WinnerRank)), 
                            rep('C',length(ratingC$rankMatch$WinnerRank)),
                            rep('D', length(ratingD$rankMatch$LoserRank)), 
                            rep('E', length(ratingE$rankMatch$LoserRank)), 
                            rep('F', length(ratingF$rankMatch$LoserRank))))]
RankTable[Pen == "A" | Pen == "B"| Pen == "C" ,Group_Size := "large"]
RankTable[Pen == "D" | Pen == "E"| Pen == "F" ,Group_Size := "small"]

RankTable[, RankDiff := abs(WinnerRank-LoserRank)]
RankTable[, HighRankWins := WinnerRank < LoserRank]
RankTable[, EloDiff := WinnerElo - LoserElo]
RankTable[, rowNum := 1:.N]
RankTable[, Situation := as.factor(Situation)]
RankTable[, Group_Size := as.factor(Group_Size)]

RankTable[Code == "Avoidance"| Code == "Threat", AggressLvl := "non_physical"]
RankTable[Code == "Peck"| Code == "Fight", AggressLvl := "physical"]
RankTable[, AggressBool := ifelse(AggressLvl == "non_physical", 0, 1)]


### Sum of Interactions #####################################

#rough overview plot of ratio of interactions by individual
ggplot(Interact[Group_Size == 'small',],aes(x = rank, y = ratio, colour = Pen)) + 
  geom_point()+
  geom_smooth(se = F)
ggplot(Interact[Group_Size == 'large',],aes(x = rank, y = ratio, colour = Pen)) + 
  geom_point()+
  geom_smooth(se = F)

#old model
#model for the number of interactions by Elo
# sum.model = glmer(sum ~ poly(scaleElos,2)*Group_Size + (1|Pen), Interact, family = 'poisson')
# resid.df2<- simulateResiduals(sum.model, 1000)
# plot(resid.df2) # overdispersion not good fit
# plotResiduals(resid.df2, form = Interact$Condition)
# plotResiduals(resid.df2, form = Interact$scaleElos)
# #try negative binomial fit to compensate overdispersal
# sum.model = glmer.nb(sum ~ poly(scaleElos,2)*Group_Size + (1|Pen), Interact)
# resid.df2<- simulateResiduals(sum.model, 1000)
# plot(resid.df2) #looks great
# plotResiduals(resid.df2, form = Interact$Condition)
# plotResiduals(resid.df2, form = Interact$scaleElos)
# 
# sum.model.null = glmer.nb(sum ~ 1 + (1|Pen), Interact) 
# anova(sum.model, sum.model.null, test = "Chisq")
# 
# #take out 2-way
# drop1(sum.model, test = "Chisq")
# sum.model.red = glmer.nb(sum ~ poly(scaleElos,2) + Group_Size + (1|Pen), Interact)
# resid.df2<- simulateResiduals(sum.model.red, 1000)
# plot(resid.df2)
# plotResiduals(resid.df2, form = Interact$Condition)
# plotResiduals(resid.df2, form = Interact$scaleElos)
# 
# summary(sum.model.red)
# summary(emtrends(sum.model.red, ~ scaleElos, "scaleElos", max.degree = 2, type = "response"), infer = c(T,T))
# summary(allEffects(sum.model.red))
#calculate in percent from hand
# round((negbinom$family$linkinv(estimate)-1)*100)
# round((negbinom$family$linkinv(estimate-sd)-1)*100)
# round((negbinom$family$linkinv(estimate+sd)-1)*100)
#https://stats.stackexchange.com/questions/365623/how-to-report-negative-binomial-regression-results-from-r

#Interact[, PredictSum := (predict(sum.model.red, type = "response"))]
#how to get CI for plot?? se.fit in predicted doesn't work for glmer
#Interact[, LCISum := exp(predict(sum.model.red)- 1.96 * se.fit)]
#Interact[, UCISum := exp(predict(sum.model.red))]
#Yamenah will check if possible, might try without


#Who accounts for how much percentage of Interactions
InteractSum = Interact[order(Pen, rank),]
InteractSum[, CumSum := cumsum(sum)/(sum(sum)*0.5), by = Pen]
InteractSum[, which(CumSum>0.8)[1]/.N, by = Pen]

#including Condition in model?

#ERST BERECHNEN WIE DIE SUMME ALLER HQ ETC IS
InteractSum[, HQ_all := HQ + HQRec]
InteractSum[, Feed_all := Feed + FeedRec]
InteractSum[, Normal_all := Normal + NormalRec]
InteractSum[, HQ_min := HQ_all/40]
InteractSum[, Feed_min := Feed_all/30]
InteractSum[, Normal_min := Normal_all/30]


InteractWide = melt(InteractSum, id.vars = c("ID","Pen", "Group_Size", "scaleElos"), 
                    measure.vars = c("HQ_all", "Feed_all", "Normal_all"),
                    variable.name = "Situation", 
                    value.name = "Sum")

InteractWideMin = melt(InteractSum, id.vars = c("ID","Pen", "Group_Size", "scaleElos"), 
                    measure.vars = c("HQ_min", "Feed_min", "Normal_min"),
                    variable.name = "Situation", 
                    value.name = "Sum")

#include observation time to use as offset in model to account for time observed
InteractWide[Situation == "HQ_all", Minutes := 40]
InteractWide[Situation == "Feed_all", Minutes := 30]
InteractWide[Situation == "Normal_all", Minutes := 30]

ggplot(InteractWide, mapping = aes(x = scaleElos, y =Sum))+#, colour = Pen)) + 
  geom_smooth()+#method = glmer.nb, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(x = 'scaled Elo rating', y = 'Number of Interactions')+
  #facet_grid(Situation~ Group_Size) + 
  theme_bw(base_size = 18)+
  scale_color_manual(values=c(largeCol, smallCol))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

hist(InteractWide$Sum)
InteractWide[, rowNum := 1:.N]
#start with max two-way interactions & poly not interacting
#sum.model = glmer(Sum ~ poly(scaleElos,2) + scaleElos*Group_Size+ Situation*Group_Size + Situation*scaleElos+(1|Pen), InteractWide, family = 'poisson')
#resid.df2<- simulateResiduals(sum.model, 1000)
#plot(resid.df2) # overdispersion not good fit -> negative binomial
sum.model = glmer.nb(Sum ~ poly(scaleElos,2) + scaleElos*Group_Size+ Situation*Group_Size + Situation*scaleElos+ offset(log(Minutes))+
                       (1|Pen), InteractWide)
resid.df2<- simulateResiduals(sum.model, 1000)
plot(resid.df2) 
plotResiduals(resid.df2, form = InteractWide$Group_Size)
plotResiduals(resid.df2, form = InteractWide$Situation) #problematic?
plotResiduals(resid.df2, form = InteractWide$scaleElos)

sum.model.null = glmer.nb(Sum ~ 1 + offset(log(Minutes))+(1|Pen), InteractWide) 
anova(sum.model, sum.model.null, test = "Chisq")

#take out 2-way
drop1(sum.model, test = "Chisq")
sum.model.red1 = glmer.nb(Sum ~ poly(scaleElos,2)+Group_Size +Situation+offset(log(Minutes))+(1|Pen), InteractWide)
anova(sum.model.red1, sum.model.null, test = "Chisq")
resid.df2<- simulateResiduals(sum.model.red1, 1000)
plot(resid.df2)
plotResiduals(resid.df2, form = InteractWide$Group_Size)
plotResiduals(resid.df2, form = InteractWide$Situation) #problematic? no homoscedacity not assumed in negative binomial
plotResiduals(resid.df2, form = InteractWide$scaleElos)


#check if situation makes model better
sum.model.red2 = glmer.nb(Sum ~ poly(scaleElos,2)+Group_Size+offset(log(Minutes))+(1|Pen), InteractWide)
anova(sum.model.red1, sum.model.red2)
sum.model.red2 = glmer.nb(Sum ~ Situation+Group_Size+offset(log(Minutes))+(1|Pen), InteractWide)
anova(sum.model.red1, sum.model.red2)

#check against poisson
#m3 <- glmer(Sum ~ poly(scaleElos,2)+Group_Size +Situation+offset(log(Minutes))+(1|Pen), family = "poisson", data = InteractWide)
#pchisq(2 * (logLik(sum.model.red1) - logLik(m3)), df = 1, lower.tail = FALSE)

#results makes sense (see plot) take results as are
ggplot(InteractWide, aes(x = Situation, y = Sum/Minutes))+
  geom_boxplot()+
  facet_grid(.~Group_Size)

summary(sum.model.red1)
summary(allEffects(sum.model.red1))
test = emmeans(sum.model.red1, ~ pairwise ~Situation, type = "response", offset = log(InteractWide$Minutes))
confint(test, adjust = "bonferroni", level = 0.95)
plot(test, comparison = TRUE) +theme_bw()
tab_model(sum.model.red1) #still very high
plot(predictorEffects(sum.model.red1), lines=list(multiline=TRUE))

#WHAT DOES THIS DO? not helpful I think
emt <- emtrends(sum.model.red1, ~ degree | scaleElos, var = "scaleElos", 
                max.degree = 2, offset = log(InteractWide$Minutes),
                type = "response", at = list(scaleElos = c(min(Interact$scaleElos), 0, max(Interact$scaleElos))))
summary(emt, infer = c(TRUE, TRUE))

InteractWide[, PredictSum := (predict(sum.model.red1, type = "response"))]
plotSums = dcast(InteractWide, ID +Pen + scaleElos + Group_Size ~ Situation, value.var = c("PredictSum", "Sum"))
plotSums[, Sum_HQ_all := (Sum_HQ_all/40)*30]
plotSums[, Sum := Sum_HQ_all + Sum_Feed_all + Sum_Normal_all]
plotSums[, PredictSum := PredictSum_HQ_all + PredictSum_Feed_all + PredictSum_Normal_all]

largeCol = brewer.pal(n = 8, name = "Blues")[c(4,6,8)]
smallCol = brewer.pal(n = 8, name = "OrRd")[c(4,6,8)]

#effect plot of Elo and number of interactions split by group

ggplot(data = plotSums, mapping = aes(x = scaleElos, y =Sum, colour = Pen)) + 
  geom_smooth(aes(x = scaleElos, y = PredictSum), se= FALSE)+#method = glmer.nb, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(x = 'scaled Elo rating', y = 'Number of Interactions')+
  facet_grid(.~ Group_Size) + 
  theme_bw(base_size = 18)+
  scale_color_manual(values=c(largeCol, smallCol))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

#effect plot of Elo and number of interactions split by Situation
ggplot(data = InteractWide, mapping = aes(x = scaleElos, y =Sum, colour = Pen)) + 
  geom_smooth(aes(x = scaleElos, y = PredictSum), se= FALSE)+#method = glmer.nb, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(x = 'scaled Elo rating', y = 'Number of Interactions')+
  facet_grid(Group_Size~ Situation) + 
  theme_bw(base_size = 18)+
  scale_color_manual(values=c(largeCol, smallCol))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())



#plot Effects of Situation
plotData <- as.data.table(emmeans(sum.model.red1, ~ pairwise ~Situation, type = "response", offset = log(InteractWide$Minutes))$emmeans)
plotData2 = InteractWide
plotData2[Situation == "HQ_all", Sum := (Sum/40)*30]

ggplot(plotData, aes(x = Situation, y = response))+
  #geom_pointrange(data = plotData2[, .(median = median(Sum), 
  #                                   q1 = quantile(Sum, 0.25),
  #                                   q2 = quantile(Sum, 0.75)),by = Situation], 
  #           aes(x = Situation, y= median, ymin = q1, ymax = q2), size = 1, colour = "grey")+
  geom_jitter(data = plotData2, aes(x = Situation, y= Sum, colour = Pen), width=0.15)+
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL), size = 1, colour = "black", shape=23, fill="yellow")+
  scale_color_manual(values=c(largeCol, smallCol))+
  labs(x = "Situation", y= "Number of interactions by individuals")+
  scale_x_discrete(labels=c("Grape", "Feed", "Normal"))+
  theme_classic(base_size = 18)

ggplot(plotData, aes(x = Situation, y = response))+
  geom_pointrange(data = plotData2[, .(median = median(Sum), 
                                     q1 = quantile(Sum, 0.25),
                                     q2 = quantile(Sum, 0.75)),by = Situation], 
             aes(x = Situation, y= median, ymin = q1, ymax = q2), size = 1, colour = "grey")+
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL), size = 1)+
  scale_color_manual(values=c(largeCol, smallCol))+
  labs(x = "Situation", y= "Predicted number of interactions by individuals")+
  scale_x_discrete(labels=c("Grape", "Feed", "Normal"))+
  theme_classic(base_size = 18)

#### Patterns of aggression ############

#strategy
diagnA$strategy
diagnD$strategy
diagnC$strategy
diagnD$strategy
diagnE$strategy
diagnF$strategy

#strategy plots
dom_plot_strategy(diagnA$focus_pos, diagnA$blur, show_data_ci = T)
dom_plot_strategy(diagnB$focus_pos, diagnB$blur, show_data_ci = T)
dom_plot_strategy(diagnC$focus_pos, diagnC$blur, show_data_ci = T)
dom_plot_strategy(diagnD$focus_pos, diagnD$blur, show_data_ci = T)
dom_plot_strategy(diagnE$focus_pos, diagnE$blur, show_data_ci = T)
dom_plot_strategy(diagnF$focus_pos, diagnF$blur, show_data_ci = T)


# Plot for Dynamics of interactions by rank (small)
ggplot(data = RankTable[Group_Size == "small",], mapping = aes(x = WinnerRank, y =LoserRank)) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  labs(y = 'Loser rank', x= 'Winner rank')+
  theme_classic(base_size = 18)+
  theme(legend.position = "none")+
  facet_grid(. ~ Pen)+
  geom_density_2d_filled(alpha = 0.7, contour_var = "ndensity")+
  geom_abline(intercept = 0 , slope = 1, linetype = "dashed", colour = 'grey')+
  geom_jitter(size = 1.2, width = 0.4)+
  xlim(0.5,20.5)+
  ylim(0.5,20.5)



# Plot for Dynamics of interactions by rank (large)
ggplot(data = RankTable[Group_Size == "large",], mapping = aes(x = WinnerRank, y =LoserRank)) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  #geom_point(size = 2.5) + 
  labs(y = 'Loser rank', x= 'Winner rank')+
  theme_classic(base_size = 18)+
  theme(legend.position = "none")+
  facet_grid(. ~ Pen)+
  #stat_density_2d(aes(fill = ..level..), geom="polygon")
  #geom_density_2d(aes(colour = Pen), size = 2)
  #geom_density_2d_filled(contour_var = "ndensity")
  #geom_density_2d_filled(contour_var = "count") 
  geom_density_2d_filled(alpha = 0.7, contour_var = "ndensity")+
  geom_abline(intercept = 0 , slope = 1, linetype = "dashed", colour = 'grey')+
  geom_jitter(size = 0.6, width = 0.4)+
  xlim(0.5,120)+
  ylim(0.5,120)

# Example plot for Dynamics of interactions by rank 
ggplot(data = RankTable[Pen == "E",], mapping = aes(x = WinnerRank, y =LoserRank)) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  labs(y = 'Loser rank', x= 'Winner rank')+
  theme_classic(base_size = 18)+
  geom_density_2d_filled(alpha = 0.7, contour_var = "ndensity")+
  geom_abline(intercept = 0 , slope = 1, linetype = "dashed", colour = 'grey')+
  geom_jitter(size = 1.5, width = 0.4)+
  xlim(0.5,20.5)+
  ylim(0.5,20.5)


########### Aggression Intensity #################

#facet by intensity
ggplot(data = RankTable[Group_Size == "small",], mapping = aes(x = WinnerRank, y =LoserRank)) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  labs(y = 'Loser rank', x= 'Winner rank')+
  theme_classic(base_size = 18)+
  facet_grid(AggressLvl ~ Pen)+
  geom_density_2d_filled(alpha = 0.7, contour_var = "count")+
  geom_abline(intercept = 0 , slope = 1, linetype = "dashed", colour = 'grey')+
  geom_jitter(size = 1, width = 0.4)+
  xlim(1,20)+
  ylim(1,20)


# facet by aggression intensity
ggplot(data = RankTable[Group_Size == "large",], mapping = aes(x = WinnerRank, y =LoserRank)) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  #geom_point(size = 2.5) + 
  labs(y = 'Loser rank', x= 'Winner rank')+
  theme_classic(base_size = 18)+
  facet_grid(AggressLvl ~ Pen)+
  geom_density_2d_filled(alpha = 0.7, contour_var = "count")+
  geom_abline(intercept = 0 , slope = 1, linetype = "dashed", colour = 'grey')+
  geom_jitter(size = 0.6, colour = 'black')+
  xlim(1,120)+
  ylim(1,120)

# non-physical set as reference 0
intensity.model = glmer( AggressBool ~ WinnerElo+ LoserElo + Situation+ Group_Size + 
                           WinnerElo:Group_Size + Situation:WinnerElo + Group_Size:Situation +
                           LoserElo:Group_Size + LoserElo:WinnerElo + LoserElo:Situation + (1|Pen), data = RankTable, family = binomial)
#singularity
intensity.model = glm( AggressBool ~ WinnerElo+ LoserElo + Situation+ Group_Size + 
                           WinnerElo:Group_Size + Situation:WinnerElo + Group_Size:Situation +
                           LoserElo:Group_Size + LoserElo:WinnerElo + LoserElo:Situation, data = RankTable, family = binomial)
resid.intensity = simulateResiduals(intensity.model)
plot(resid.intensity)
plotResiduals(resid.intensity, form = RankTable$DiffElo) #nearly perfect
plotResiduals(resid.intensity, form = RankTable$Group_Size) #heterogenity but okay
plotResiduals(resid.intensity, form = RankTable$Situation) #good

intensity.null = glm( AggressBool ~ 1 , data = RankTable, family = binomial)
anova(intensity.model, intensity.null, test = "Chisq")
library(car)
vif(intensity.model)
#drop 2-way
drop1(intensity.model, k = 2) #keep situation*Groupsize & Winner*Loser
intensity.model = glm(AggressBool ~ WinnerElo+ LoserElo + Situation+ Group_Size + 
                                         Group_Size:Situation + LoserElo:WinnerElo, data = RankTable, family = binomial)
resid.intensity = simulateResiduals(intensity.model)
plot(resid.intensity)
plotResiduals(resid.intensity, form = RankTable$DiffElo) #nearly perfect
plotResiduals(resid.intensity, form = RankTable$Group_Size) #heterogenity
plotResiduals(resid.intensity, form = RankTable$Situation) #good
intensity.model2 = glm(AggressBool ~ WinnerElo+ LoserElo + Situation+ Group_Size + 
                        Group_Size:Situation, data = RankTable, family = binomial)


anova(intensity.null,intensity.model, test = "Chisq")
anova(intensity.model, intensity.model2, test = "Chisq")


summary(intensity.model)
parameters(intensity.model, exponentiate = T)
group_sit = emmeans(intensity.model, ~ pairwise ~ Group_Size*Situation, type = "response")
emmeans(intensity.model, ~ pairwise ~ WinnerElo*LoserElo, type = "response")
emmeans(intensity.model, ~ WinnerElo, 
        at = list(WinnerElo = c(min(RankTable$WinnerElo),mean(RankTable$WinnerElo),max(RankTable$WinnerElo))), type= "response")
confint(group_sit, adjust = "bonferroni", level = 0.95)
plot(allEffects(intensity.model))

tab_model(intensity.model)
plot_model(intensity.model, show.value = T)

p2 = plot_model(intensity.model, type = "pred", terms = c("WinnerElo[all]", "LoserElo[-2, 0, 2]"))

p2 + theme_classic(base_size = 18) + labs(title = "", x = "Elo rating of the winner", y = "Probability for high intensity aggression",
                                          colour = "Elo rating of \nthe loser")+ geom_line(size = 2)

##
#Effect plots
#Aggression intensity by WinnerElo -> can't be done because part of interaction

RankTable[ ,PredictIntens := predict(intensity.model, type = "response")]

ggplot(data = RankTable, aes(x = WinnerElo, y = PredictIntens, colour = factor(AggressBool)))+
  geom_point()+
  #facet_grid(.~Pen)+
  theme_bw(base_size = 18)

ggplot(data = Interact, aes(x = scaleElos))+
  geom_point(aes(y = physAggr), size = 2, colour = "red")+
  geom_smooth(aes(y = physAggr),colour = "red")+
  geom_point(aes(y = nonphysAggr), size = 2, colour = "blue")+
  geom_smooth(aes(y = nonphysAggr),colour = "blue")+
    #facet_grid(.~Pen)+
  theme_bw(base_size = 18)

ggplot(data = na.omit(Interact), aes(x = scaleElos, y = ratioIntens))+
  geom_point( size = 2)+
  geom_smooth()+
  #facet_grid(.~Pen)+
  theme_bw(base_size = 18)

plot_model(intensity.model, type = "pred", terms = c("WinnerElo[all]"))+
  geom_rug()+
  labs(x = "Elo of winner", y = "predicted probability for physical aggression")+
  theme_classic(base_size = 18)

#Aggression intensity by Situation and Group Size

plot_model(intensity.model, type = "pred", terms = c("Group_Size", "Situation"),)+
  theme_classic(base_size = 18)

  
plotData <- as.data.table(emmeans(intensity.model, ~ pairwise ~ Situation*Group_Size, type = "response")$emmeans)

observedProb = RankTable[, .(prob = sum(AggressBool)/.N), by= .(Group_Size, Situation)]

contrast = as.data.table(emmeans(intensity.model, ~ pairwise ~ Group_Size*Situation, type = "response")$contrasts)


library(ggsignif)
ggplot(plotData, aes(x = Group_Size, y = prob))+
  geom_pointrange(aes(shape = Situation, ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge(width = 0.75), size = 1.5)+
  #geom_point(data = observedProb, aes(shape = Situation), position = position_jitterdodge(jitter.width = 0.3) , size = 3, colour = "grey")+
  labs(x = "Group size", y= "Predicted probability for high intensity aggression")+
  ylim(0.1,1)+
  scale_shape_discrete(labels=c('Feeder', 'Grape', 'Normal'))+
  theme_classic(base_size = 18)
  
# p +   geom_signif( #large Feeder vs. large HQ
#   annotation = "**",#formatC(contrast$p.value[2], digits = 1),
#   y_position = 0.65, xmin = 0.75, xmax = 0.99,
#   tip_length = c(0.2, 0.1),
#   textsize = 8,
#   size = 1
# ) +
#   
#   geom_signif( #large HQ vs. large Normal
#     annotation = "***",#formatC(contrast$p.value[11], digits = 1),
#     y_position = 0.65, xmin = 1.01, xmax = 1.25,
#     tip_length = c(0.1, 0.3),
#     textsize = 8,
#     size = 1
#   )+
#   
#   geom_signif( #large Feeder vs. large Normal
#     annotation = "0.06",#formatC(contrast$p.value[11], digits = 1),
#     y_position = 0.75, xmin = 0.75, xmax = 1.25,
#     tip_length = c(0.1, 0.1),
#     textsize = 5,
#     size = 1
#   )+
# geom_signif( #small Feeder vs. small HQ
#   annotation = "***",#formatC(contrast$p.value[7], digits = 1),
#   y_position = 0.8, xmin = 1.75, xmax = 1.99,
#   tip_length = c(0.2, 0.04),
#   textsize = 8,
#   size = 1
# ) +
#   
#   geom_signif( #small HQ vs. small Normal
#     annotation = "**",#formatC(contrast$p.value[14], digits = 1),
#     y_position = 0.8, xmin = 2.01, xmax = 2.25,
#     tip_length = c(0.04, 0.2),
#     textsize = 8,
#     size = 1
#   )+
#   
#   geom_signif( #small HQ vs. large HQ
#     annotation = "*",#formatC(contrast$p.value[10], digits = 1),
#     y_position = 0.95, xmin = 1, xmax = 2,
#     tip_length = c(0.25, 0.1),
#     textsize = 8,
#     size = 1
#   )


####################### Individual contacts data ######################################

#How many individuals does one hen interact with? 

#create unique identifier
RankTable[, IDWinner := paste0(Pen, Winner)]
RankTable[, IDLoser := paste0(Pen, Loser)]
Interact[, uniqueID := paste0(Pen, ID)]

#extract unique Interaction companions
UniqueContacts1 = RankTable[,.(InteractUniq = unique(IDLoser) ), by = (IDWinner)]
UniqueContacts2 = RankTable[,.(InteractUniq = unique(IDWinner) ), by = (IDLoser)]
#rename column to same name for rbind
colnames(UniqueContacts1)[1] <- "uniqueID" 
colnames(UniqueContacts2)[1] <- "uniqueID"
UniqueContacts = rbind(UniqueContacts1, UniqueContacts2)

#order and delete duplicate entries
UniqueContacts = UniqueContacts[order(uniqueID, InteractUniq), ]
UniqueContacts[, dupl := rleid(uniqueID, InteractUniq)]
UniqueContacts[, duplBool := duplicated(dupl)]
UniqueContacts = UniqueContacts[duplBool == F,]

#summary per individual
UniqueContactsN = UniqueContacts[, .(Contacts = .N), by = "uniqueID"] 

#add to Individual overview table
Interact = Interact[UniqueContactsN, on = .(uniqueID)]

#summaries
Interact[, .(max = max(Contacts), min = min(Contacts), median = median(Contacts)), by = Pen]

# Plot for Dynamics of interactions by Elo
ggplot(data = Interact, mapping = aes(x = scaleElos, y =Contacts, colour = Pen)) + 
  geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(x = 'Rank', y = 'Number of unique Contacts')+
  theme_classic(base_size = 18)

#Do large group high ranking interact with other high ranking?
RankTable[WinnerRank == 1, .(min = min(LoserRank), max = max(LoserRank)), by = Pen]
RankTable[WinnerRank == 2, .(min = min(LoserRank), max = max(LoserRank)), by = Pen]
RankTable[WinnerRank == 3, .(min = min(LoserRank), max = max(LoserRank)), by = Pen]
RankTable[WinnerRank == 4, .(min = min(LoserRank), max = max(LoserRank)), by = Pen]
RankTable[WinnerRank == 5, .(min = min(LoserRank), max = max(LoserRank)), by = Pen]


###### Examples of steepness plots ##############

#examples for hierarchy shape

ratingA$hierarchy_shape_rand$plot

library(ggpubr)
p1 = ratingA$hierarchy_shape_rand$plot+
  scale_y_continuous(limits=c(0.25,1),oob = rescale_none)+
  scale_x_continuous(breaks = c(1,seq(0,120, 10)[2:12], 119),oob = rescale_none)+
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text(angle = 90)
  )+
  labs(x = "Difference in rank", y = "Probability that high rank wins")
p2 = ratingE$hierarchy_shape_rand$plot+
  scale_y_continuous(limits=c(0.25,1),oob = rescale_none)+
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text(angle = 90)
  )+
  labs(x = "Difference in rank", y = "Probability that high rank wins")


##### Examples of hierarchy plots #######

mean.eloM = rowMeans(ratingF$randElo)
identitiesM = rownames(ratingF$randElo)
identitiesM <- identitiesM[order(mean.eloM)]
CIsM <- apply(ratingF$randElo,1,quantile,c(0.025,0.975),na.rm=TRUE)
CIsM <- CIsM[,order(mean.eloM)]
mean.eloM <- mean.eloM[order(mean.eloM)]

plotTable1 = data.table(number = 1:20,  
                        ranks = mean.eloM,
                        IDs = identitiesM
)
plotTable1[, CIl := CIsM[1,]]
plotTable1[, CIu := CIsM[2,]]

ggplot(plotTable1, aes(x = number, y= ranks))+
  geom_errorbar(aes(x= number, ymin=CIl, ymax=CIu), width=.1)+
  geom_point(size = 2)+
  #geom_point(shape = 21, size = 4.7, colour = "white", fill = "white", stroke = 1)+
  #scale_y_continuous(limits = c(-450, 650))+
  #geom_text(
  #  label=plotTable$IDs,
  #  size = 4)+
  labs(x = "Individuals", y="Dominance rank")+
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())  


mean.eloM = rowMeans(ratingB$randElo)
identitiesM = rownames(ratingB$randElo)
identitiesM <- identitiesM[order(mean.eloM)]
CIsM <- apply(ratingB$randElo,1,quantile,c(0.025,0.975),na.rm=TRUE)
CIsM <- CIsM[,order(mean.eloM)]
mean.eloM <- mean.eloM[order(mean.eloM)]

plotTable2 = data.table(number = 1:length(mean.eloM),  
                       ranks = mean.eloM,
                       IDs = identitiesM
)
plotTable2[, CIl := CIsM[1,]]
plotTable2[, CIu := CIsM[2,]]

ggplot(plotTable2, aes(x = number, y= ranks))+
  geom_errorbar(aes(x= number, ymin=CIl, ymax=CIu), width=.1)+
  geom_point( size = 2)+
  #scale_y_continuous(limits = c(-450, 650))+
  #geom_text(
  #  label=plotTable$IDs,
  #  size = 4)+
  labs(x = "Individuals", y="Elo-rating")+
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 


plotTable1[, Group_Size := "small"]
plotTable2[, Group_Size := "large"]

plotTable = rbind(plotTable1, plotTable2)
plotTable$Group_Size = factor(plotTable$Group_Size, levels=c("small", "large"), labels = c("small group", "large group"))
BCol = brewer.pal(n = 8, name = "Blues")[6]
FCol = brewer.pal(n = 8, name = "OrRd")[8]

pp = ggplot(plotTable, aes(x = number, y= ranks, colour = Group_Size) )+
  geom_errorbar(aes(x= number, ymin=CIl, ymax=CIu), width=.1)+
  geom_point( size = 2)+
  #scale_y_continuous(limits = c(-450, 650))+
  #geom_text(
  #  label=plotTable$IDs,
  #  size = 4)+
  facet_grid( . ~ Group_Size, scales = "free")+
  scale_color_manual(values=c(FCol, BCol))+
  labs(x = "Individuals", y="Dominance rank")+
  scale_y_continuous(breaks=c(-400, 700),
                     labels=c("low rank", "high rank"))+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(angle=90)) 
ggsave("HierarchyEx.tiff", pp, "tiff", width = 26, height = 10, units = "cm", dpi = 300)


##### plots for poster ####
plotHelper = Interact

plotHelper$Group_Size = factor(Interact$Group_Size, levels=c("small", "large"), labels = c("small group", "large group"))

p = ggplot(data =plotHelper, mapping = aes(x = scaleElos, y =sum, colour = Pen)) + 
  geom_smooth(aes(x = scaleElos, y = PredictSum), se= FALSE)+#method = glmer.nb, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(x = 'Dominance rank', y = 'Number of interactions')+
  facet_grid(.~ Group_Size) + 
  theme_bw(base_size = 18)+
  scale_color_manual(values=c(largeCol, smallCol))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(-2, 2.5),
                   labels=c("low rank", "high rank"))+
  theme(axis.ticks.x=element_blank()) 
ggsave("NumberInteractions.tiff", p, "tiff", width = 20, height = 10, units = "cm", dpi = 300)
