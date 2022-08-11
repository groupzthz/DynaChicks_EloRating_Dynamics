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

Interact = na.omit(Interact)

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


############Sum of Interactions #####################################

#rough overview plot of ratio of interactions by individual
ggplot(Interact[Group_Size == 'small',],aes(x = rank, y = ratio, colour = Pen)) + 
  geom_point()+
  geom_smooth(se = F)
ggplot(Interact[Group_Size == 'large',],aes(x = rank, y = ratio, colour = Pen)) + 
  geom_point()+
  geom_smooth(se = F)


#model for the number of interactions by Elo
sum.model = glmer(sum ~ poly(scaleElos,2)*Group_Size + (1|Pen), Interact, family = 'poisson')
resid.df2<- simulateResiduals(sum.model, 1000)
plot(resid.df2) # overdispersion not good fit
plotResiduals(resid.df2, form = Interact$Condition)
plotResiduals(resid.df2, form = Interact$scaleElos)
#try negative binomial fit to compensate overdispersal
sum.model = glmer.nb(sum ~ poly(scaleElos,2)*Group_Size + (1|Pen), Interact)
resid.df2<- simulateResiduals(sum.model, 1000)
plot(resid.df2) #looks great
plotResiduals(resid.df2, form = Interact$Condition)
plotResiduals(resid.df2, form = Interact$scaleElos)

sum.model.null = glmer.nb(sum ~ 1 + (1|Pen), Interact)#throws weird error 
anova(sum.model, sum.model.null, test = "Chisq")

#take out 2-way
drop1(sum.model, test = "Chisq")
sum.model.red = glmer.nb(sum ~ poly(scaleElos,2)+ Group_Size + (1|Pen), Interact)
resid.df2<- simulateResiduals(sum.model.red, 1000)
plot(resid.df2)
plotResiduals(resid.df2, form = Interact$Condition)
plotResiduals(resid.df2, form = Interact$scaleElos)

summary(sum.model.red)
parameters(sum.model.red, exponentiate = T)

Interact[, PredictSum := (predict(sum.model.red, type = "response"))]
#TODO: how to get CI for plot?? se.fit in predicted doesn't work for glmer
#Interact[, LCISum := exp(predict(sum.model.red)- 1.96 * se.fit)]
#Interact[, UCISum := exp(predict(sum.model.red))]


largeCol = brewer.pal(n = 8, name = "Blues")[c(4,6,8)]
smallCol = brewer.pal(n = 8, name = "OrRd")[c(4,6,8)]

#effect plot of Elo and number of interactions split by group
ggplot(data = Interact, mapping = aes(x = scaleElos, y =sum, colour = Pen)) + 
       geom_smooth(aes(x = scaleElos, y = PredictSum), se= FALSE)+#method = glmer.nb, formula = y ~ splines::bs(x, 3), se = FALSE)+
       geom_point(size = 2.5) + 
       labs(x = 'scaled Elo rating', y = 'Number of Interactions')+
       facet_grid(.~ Group_Size) + 
      theme_bw(base_size = 18)+
   scale_color_manual(values=c(largeCol, smallCol))+
     theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

#Who accounts for how much percentage of Interactions
InteractSum = Interact[order(Pen, rank),]
InteractSum[, CumSum := cumsum(sum)/(sum(sum)*0.5), by = Pen]
InteractSum[, which(CumSum>0.75)[1]/.N, by = Pen]


#### Patterns of aggression ############


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


########### AGGRESSION Intensity #################

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
plotResiduals(resid.intensity, form = RankTable$Group_Size) #heterogenity
plotResiduals(resid.intensity, form = RankTable$Situation) #good
#TODO adjust heterogneity correct with 1|rowNum?? shouts signularity issue now but hetero gone...
intensity.model = glmer( AggressBool ~ WinnerElo+ LoserElo + Situation+ Group_Size + 
                           WinnerElo:Group_Size + Situation:WinnerElo + Group_Size:Situation +
                           LoserElo:Group_Size + LoserElo:WinnerElo + LoserElo:Situation + (1|rowNum), data = RankTable, family = binomial)

intensity.null = glmer( AggressBool ~ 1+ (1|rowNum), data = RankTable, family = binomial)
anova(intensity.model, intensity.null, test = "Chisq")

#drop 2-way
drop1(intensity.model, test = "Chisq") #keeo situation*Groupsize & Winner*Loser
intensity.model = glmer(AggressBool ~ WinnerElo+ LoserElo + Situation+ Group_Size + 
                                         Group_Size:Situation + LoserElo:WinnerElo + (1|rowNum), data = RankTable, family = binomial)
resid.intensity = simulateResiduals(intensity.model)
plot(resid.intensity)
plotResiduals(resid.intensity, form = RankTable$DiffElo) #nearly perfect
plotResiduals(resid.intensity, form = RankTable$Group_Size) #heterogenity
plotResiduals(resid.intensity, form = RankTable$Situation) #good

anova(intensity.model, intensity.null, test = "Chisq")

summary(intensity.model)
parameters(intensity.model, exponentiate = T)
emmeans(intensity.model, ~ pairwise ~ Group_Size | Situation, type = "response")
emmeans(intensity.model, ~ pairwise ~ Situation | Group_Size, type = "response")

plot(allEffects(intensity.model))
library(sjPlot)
tab_model(intensity.model)
plot_model(intensity.model, show.value = T)

plot_model(intensity.model, type = "pred", terms = c("WinnerElo[all]", "LoserElo[-2, 0, 2]"))

##
#Effect plots
#Aggression intensity by WinnerElo 
#TODO: how can the effect be so strong if not visible?
RankTable[ ,PredictIntens := predict(intensity.model, type = "response")]

ggplot(data = RankTable, aes(x = WinnerElo, y = PredictIntens, colour = factor(AggressBool)))+
  geom_point()+
  #facet_grid(.~Pen)+
  theme_bw(base_size = 18)

ggplot(data = Interact, aes(x = scaleElos))+
  geom_point(aes(y = physAggr), size = 2, colour = "red")+
  geom_smooth(aes(x = scaleElos, y = physAggr),colour = "red")+
  geom_point(aes(y = nonphysAggr), size = 2, colour = "blue")+
  geom_smooth(aes(x = scaleElos, y = nonphysAggr),colour = "blue")+
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
  theme_classic(base_size = 18)+

  
plotData <- as.data.table(emmeans(intensity.model, ~ pairwise ~ Situation | Group_Size, type = "response")$emmeans)



# DIFF IN ELO


ggplot(data = RankTable, aes(x = DiffElo, y = AggressBool))+
  geom_point()+
  geom_smooth(formula = y~x)

ggplot(data = RankTable[, .N, by = .(Condition, AggressBool)], aes(x = Condition, fill = as.factor(AggressBool), y = N))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)

ggplot(data = RankTable, aes(x = DiffElo, y = AggressBool))+
  geom_point()+
  geom_smooth(formula = y~x)+
  facet_grid(.~Condition)


# WINNER AND LOSER ELO

ggplot(data = RankTable, aes(x = WinnerElo, y = AggressBool))+
  geom_point()+
  geom_smooth(formula = y~x, method="glm",method.args=list(family="binomial"))
ggplot(data = RankTable, aes(x = LoserElo, y = AggressBool))+
  geom_point()+
  geom_smooth(formula = y~x, method="glm",method.args=list(family="binomial"))

ggplot(data = RankTable, aes(x = WinnerElo, y = AggressBool))+
  geom_point()+
  geom_smooth(formula = y~x, method="glm",method.args=list(family="binomial"))+
  facet_grid(.~Condition)
ggplot(data = RankTable, aes(x = LoserElo, y = AggressBool))+
  geom_point()+
  geom_smooth(formula = y~x, method="glm",method.args=list(family="binomial"))+
  facet_grid(.~Condition)

# non-physical set as reference 0
test = glmer( AggressBool ~ WinnerElo*LoserElo*Condition + (1|Pen), data = RankTable, family = binomial)
resid.test = simulateResiduals(test)
plot(resid.test)
plotResiduals(resid.test, form = RankTable$WinnerElo) #okay?
plotResiduals(resid.test, form = RankTable$LoserElo) 
plotResiduals(resid.test, form = RankTable$Condition) 
test_null = glmer( AggressBool ~ 1+ (1|Pen), data = RankTable, family = binomial)
anova(test, test_null, test = "Chisq")

#drop 3-way
drop1(test, test = "Chisq")
test = glmer( AggressBool ~ WinnerElo*Condition + WinnerElo*LoserElo + LoserElo*Condition + (1|Pen), data = RankTable, family = binomial)
resid.test = simulateResiduals(test)
plot(resid.test)
plotResiduals(resid.test, form = RankTable$WinnerElo) #okay?
plotResiduals(resid.test, form = RankTable$LoserElo) 
plotResiduals(resid.test, form = RankTable$Condition) 

#drop 2-way
drop1(test, test = "Chisq")
test = glmer( AggressBool ~ WinnerElo*LoserElo+ Condition+ (1|Pen), data = RankTable, family = binomial)
resid.test = simulateResiduals(test)
plot(resid.test)
plotResiduals(resid.test, form = RankTable$WinnerElo) #okay?
plotResiduals(resid.test, form = RankTable$LoserElo) 
plotResiduals(resid.test, form = RankTable$Condition) 


anova(test, test_null, test = "Chisq")
summary(test)
parameters(test, exponentiate = T)
emmeans(test, ~ Condition, type = "response")
emmeans(test, ~ WinnerElo, type = "response")
plot(allEffects(test))
#significant but meaningless effect?? (see magnitude)




#other plots
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


##### Individual contacts data ###########

RankTable[, IndividWin := paste0(Pen, Winner)]
RankTable[, IndividLos := paste0(Pen, Loser)]
Interact[, Individual := paste0(Pen, ID)]

UniqueContacts1 = RankTable[,.(InteractUniq = unique(IndividLos) ), by = (IndividWin)]
UniqueContacts1 = RankTable[,.(InteractUniq = unique(IndividWin) ), by = (IndividLos)]
colnames(UniqueContacts1)[1] <- "Individual" 
colnames(UniqueContacts2)[1] <- "Individual"

UniqueContacts = rbind(UniqueContacts1, UniqueContacts2)
UniqueContacts = UniqueContacts[order(Individual, InteractUniq), ]
UniqueContacts[, dupl := rleid(Individual, InteractUniq)]
UniqueContacts[, duplBool := duplicated(dupl)]
UniqueContacts = UniqueContacts[duplBool == F,]

UniqueContactsN = UniqueContacts[, .(Contacts = .N), by = "Individual"] 

Interact = Interact[UniqueContactsN, on = .(Individual)]

Interact[, .(max = max(Contacts), min = min(Contacts), median = median(Contacts)), by = Pen]

# Plot for Dynamics of interactions by rank (small)
ggplot(data = Interact[Condition == 'small',], mapping = aes(x = rank, y =Contacts, colour = Pen)) + 
  geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(x = 'Rank', y = 'Number of unique Contacts')+
  theme_classic(base_size = 18)


# Plot for Dynamics of interactions by rank (large)
ggplot(data = Interact[Condition == 'large',], mapping = aes(x = rank, y = Contacts, colour = Pen)) + 
  geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(x = 'Rank', y = 'Number of unique Contacts')+
  theme_classic(base_size = 18)

Contacts.model = glmer(Contacts ~ poly(scaleElos,2)*Condition + (1|Pen) + (1|rowNum), Interact, family = 'poisson')
resid.df2<- simulateResiduals(Contacts.model, 1000)
plot(resid.df2)
plotResiduals(resid.df2, form = Interact$Condition)
plotResiduals(resid.df2, form = Interact$scaleElos)
#heterogenity and uniformity issues

Contacts.model.null = glmer(Contacts ~ 1 + (1|Pen) + (1|rowNum), Interact, family = 'poisson')
anova(Contacts.model, Contacts.model.null, test = "Chisq")

#check negative binomial fit
Contacts.model2 = glmer.nb(Contacts ~ poly(scaleElos,2)*Condition + (1|Pen), Interact)
resid.df2<- simulateResiduals(Contacts.model2, 1000) # heterogenity
plot(resid.df2)
plotResiduals(resid.df2, form = Interact$Condition) # heterogenity
plotResiduals(resid.df2, form = Interact$scaleElos)

library("glmmTMB")
Contacts.nbm0 = glmmTMB(Contacts~poly(scaleElos,2)*Condition + (1|Pen),disp=~Condition, Interact, family=nbinom2,
                        control = glmmTMBControl(optimizer = optim, 
                                                 optArgs = list(method="BFGS",iter.max=1e50,eval.max=1e50)))
resid.df2<- simulateResiduals(Contacts.nbm0, 1000) # heterogenity
plot(resid.df2)
plotResiduals(resid.df2, form = Interact$Condition) # heterogenity
plotResiduals(resid.df2, form = Interact$scaleElos)

drop1(Contacts.nbm0, test = "Chisq") # nope no drop
summary(Contacts.nbm0)
# important though????? -> obvious that large groups can have more contacts, median more interesting

#Do large group high ranking inetract with other high ranking?

###### Steepness plots ##############

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


##### Hierrachy plots #######

mean.eloM = rowMeans(ratingE$randElo)
identitiesM = rownames(ratingE$randElo)
identitiesM <- identitiesM[order(mean.eloM)]
CIsM <- apply(ratingE$randElo,1,quantile,c(0.025,0.975),na.rm=TRUE)
CIsM <- CIsM[,order(mean.eloM)]
mean.eloM <- mean.eloM[order(mean.eloM)]

plotTable = data.table(number = 1:20,  
                        ranks = mean.eloM,
                        IDs = identitiesM
)

ggplot(plotTable, aes(x = number, y= ranks))+
  geom_errorbar(aes(x= number, ymin=CIsM[1,], ymax=CIsM[2,]), width=.1)+
  geom_point(shape = 21, size = 4.7, colour = "white", fill = "white", stroke = 1)+
  #scale_y_continuous(limits = c(-450, 650))+
  geom_text(
    label=plotTable$IDs,
    size = 4)+
  labs(x = "Individuals", y="Elo-rating")+
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  


mean.eloM = rowMeans(ratingB$randElo)
identitiesM = rownames(ratingB$randElo)
identitiesM <- identitiesM[order(mean.eloM)]
CIsM <- apply(ratingB$randElo,1,quantile,c(0.025,0.975),na.rm=TRUE)
CIsM <- CIsM[,order(mean.eloM)]
mean.eloM <- mean.eloM[order(mean.eloM)]

plotTable = data.table(number = 1:length(mean.eloM),  
                       ranks = mean.eloM,
                       IDs = identitiesM
)

ggplot(plotTable, aes(x = number, y= ranks))+
  geom_errorbar(aes(x= number, ymin=CIsM[1,], ymax=CIsM[2,]), width=.1)+
  geom_point( size = 2)+
  #scale_y_continuous(limits = c(-450, 650))+
  #geom_text(
  #  label=plotTable$IDs,
  #  size = 4)+
  labs(x = "Individuals", y="Elo-rating")+
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 




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
modelInteract2 = glmer.nb(sum ~ poly(scaleElos,2)*Condition + (1|Pen), Interact)
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
