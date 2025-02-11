getwd()

#clear workspace
rm(list = ls())


library(aniDom) # randomised Elo-rating
library(EloRating) # original Elo rating
library(domstruc) # pattern of aggression (who interacts with whom)
library(data.table) # data.frame better
#library(compete) # needed for aniDom (poisson) calculation
library(ggplot2) # plotting
library(lme4) # mixed models
library(multcomp) # posthoc analysis
library(DHARMa) # model diagnistics
library(moments) # normal asumptions
library(MuMIn) # model comparison and pseudo R²
library(irr) # ICC calculations
library(parameters) # model parameters
library(emmeans) #model means
library(effects) #model effects visualisation
library(RColorBrewer) # color for plotting
library(EloSteepness) # for steepness measure
library(ggpubr)#plot combination

#library(sjPlot)
#library(nlstools)
#library(elo)
#library(binom) #confidence intervall for binomial sample
source('helper_functions_Elo.R')
source('plot_utils.R')
set.seed(42)

#### Loading and preparing Data ###########
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

#clean workspace
rm(dataMature, 
   PenAt, PenBt, PenCt, PenDt, PenEt, PenFt)

#### Diagnostics, transitivity ###############################################################

# Dataset analytics:
#number of observations by resource condition
PenA[, .N, by = Condition]
PenB[, .N, by = Condition]
PenC[, .N, by = Condition]
PenD[, .N, by = Condition]
PenE[, .N, by = Condition]
PenF[, .N, by = Condition]

#number of observations by observed interaction type
PenA[, .N, by = Code]
PenB[, .N, by = Code]
PenC[, .N, by = Code]
PenD[, .N, by = Code]
PenE[, .N, by = Code]
PenF[, .N, by = Code]


# were all inidviduals observed?
# total number of interactions
# ratio of interactions per bird 
# proportion of unknown pair-relationships
# expected proportion of known relationships + CI under Poisson distribution 
# actual number of unknown relationships
# triangle transitivity
# focus and position + CIs
# simulated data of downward heuristic reference
# estimated pattern of aggression
diagnA = dataset_diagnostics(PenA, Individ$ID[Individ$Pen == 'A'])
diagnB = dataset_diagnostics(PenB, Individ$ID[Individ$Pen == 'B'])
diagnC = dataset_diagnostics(PenC, Individ$ID[Individ$Pen == 'C'])
diagnD = dataset_diagnostics(PenD, Individ$ID[Individ$Pen == 'D'])
diagnE = dataset_diagnostics(PenE, Individ$ID[Individ$Pen == 'E'])
diagnF = dataset_diagnostics(PenF, Individ$ID[Individ$Pen == 'F'])

#### Elo rating calculation ################################################################

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


#### Bayesian Steepness ##########################################################

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

#example of one output
summary(steepD)
plot_steepness(steepD)
plot_scores(steepD)
plot_steepness_regression(steepD, width_fac = 5)
D_scores <- as.data.table(scores(steepD, quantiles = c(0.055, 0.945)))

#### Data tables for analyses #################
#colors for large and small groups
largeCol = brewer.pal(n = 8, name = "Blues")[c(4,6,8)]
smallCol = brewer.pal(n = 8, name = "OrRd")[c(4,6,8)]

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



#### Sum of Interactions #####################################

#Descriptives
#divide number of interactions to not count them twice
descript = Interact[, .(Sum = sum(sum)/2, N = .N), by = Group_Size]
descript[, rate := (Sum/N/100)*60]
descript = Interact[, .(Sum = sum(sum)/2, N = .N), by = Pen]
descript[, rate := (Sum/N/100)*60]
#Interact[, summary(sum), by = Group_Size]

#rough overview plot of ratio of interactions by individual
ggplot(Interact[Group_Size == 'small',],aes(x = rank, y = ratio, colour = Pen)) + 
  geom_point()+
  geom_smooth(se = F)
ggplot(Interact[Group_Size == 'large',],aes(x = rank, y = ratio, colour = Pen)) + 
  geom_point()+
  geom_smooth(se = F)


#Who accounts for how much percentage of Interactions
InteractSum = Interact[order(Pen, rank),]
InteractSum[, CumSum := cumsum(sum)/(sum(sum)*0.5), by = Pen]
InteractSum[, which(CumSum>0.8)[1]/.N, by = Pen]


#First calculate sums within conditions and add offset  
InteractSum[, HQ_all := HQ + HQRec]
InteractSum[, Feed_all := Feed + FeedRec]
InteractSum[, Normal_all := Normal + NormalRec]
InteractSum[, HQ_min := HQ_all/40]
InteractSum[, Feed_min := Feed_all/30]
InteractSum[, Normal_min := Normal_all/30]

#change to wide format
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

#desriptive for condition
InteractWideMin[, .(Mean = mean(Sum*60/2), SD = sd(Sum*60/2)), by = c("Situation")]
ggplot(InteractWide, aes(x = Situation, y = Sum/Minutes))+
  geom_boxplot()+
  facet_grid(.~Group_Size)

#look at data distribution
ggplot(InteractWide, mapping = aes(x = scaleElos, y =Sum))+#, colour = Pen)) + 
  geom_point(size = 2.5) + 
  labs(x = 'scaled Elo rating', y = 'Number of Interactions')+
  #facet_grid(Situation~ Group_Size) + 
  theme_bw(base_size = 18)+
  scale_color_manual(values=c(largeCol, smallCol))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

#look at data distribution
hist(InteractWide$Sum) #-> poisson model

##### statistical model ####
# poly(raw= T) otherwise impossible estimates with raw = F 
sum.model = glmer(Sum ~ poly(scaleElos,2, raw = TRUE) + scaleElos*Group_Size+ Situation*Group_Size + Situation*scaleElos+ offset(log(Minutes))+
                       (1|Pen), InteractWide, family = "poisson")
resid.df2<- simulateResiduals(sum.model, 1000)
plot(resid.df2) #overdispersal!! -> negative binomial model necessary

sum.model = glmer.nb(Sum ~ poly(scaleElos,2, raw = TRUE) + scaleElos*Group_Size+ Situation*Group_Size + Situation*scaleElos+ offset(log(Minutes))+
                       (1|Pen), InteractWide)
resid.df2<- simulateResiduals(sum.model, 1000)
plot(resid.df2) 
plotResiduals(resid.df2, form = InteractWide$Group_Size) #good
plotResiduals(resid.df2, form = InteractWide$Situation) #problematic? no homoscedacity not assumed in negative binomial
plotResiduals(resid.df2, form = InteractWide$scaleElos) #good

sum.model.null = glmer.nb(Sum ~ 1 +(1|Pen), InteractWide) 
anova(sum.model, sum.model.null) #better than intercept only

#simplify model
#take out 2-way
drop1(sum.model) # only removing polynomial effect would make the model much worse

#reduced final model
sum.model.red1 = glmer.nb(Sum ~ poly(scaleElos,2, raw = TRUE)+Group_Size +Situation+offset(log(Minutes))+(1|Pen), InteractWide)
resid.df2<- simulateResiduals(sum.model.red1, 1000)
plot(resid.df2)
plotResiduals(resid.df2, form = InteractWide$Group_Size) #good
plotResiduals(resid.df2, form = InteractWide$Situation) #problematic? no homoscedacity not assumed in negative binomial
plotResiduals(resid.df2, form = InteractWide$scaleElos) # good
AIC(sum.model.red1, sum.model.null)

#model estimate
summary(sum.model.red1)
parameters(sum.model.red1, exp = TRUE)
tab_model(sum.model.red1)
plot(predictorEffects(sum.model.red1), lines=list(multiline=TRUE))
summary(allEffects(sum.model.red2))

#estimate variance explained
r.squaredGLMM(sum.model.red1, sum.model.null)

#Post-hoc comparison of conditions
postHocCond = emmeans(sum.model.red1, ~ pairwise ~Situation, 
                      type = "response", offset = log(InteractWide$Minutes))
confint(postHocCond, adjust = "bonferroni", level = 0.95)

#test model without resource conditions
testInteract = InteractWide[Situation == "HQ_all"]

hist(testInteract$Sum) #-> poisson model

# poly(raw= T) otherwise impossible estimates with raw = F 
sum.model = glmer.nb(Sum ~ poly(scaleElos,2, raw = TRUE) + scaleElos*Group_Size+(1|Pen), testInteract)
resid.df2<- simulateResiduals(sum.model, 1000)
plot(resid.df2) 
plotResiduals(resid.df2, form = testInteract$Group_Size) #good
plotResiduals(resid.df2, form = testInteract$scaleElos) #good

drop1(sum.model)
sum.model = glmer.nb(Sum ~ poly(scaleElos,2, raw = TRUE) + Group_Size+(1|Pen), testInteract)
resid.df2<- simulateResiduals(sum.model, 1000)
plot(resid.df2) 
plotResiduals(resid.df2, form = testInteract$Group_Size) #good
plotResiduals(resid.df2, form = testInteract$scaleElos) #good

sum.model.null = glmer.nb(Sum ~ 1+(1|Pen), testInteract)
AIC(sum.model.null, sum.model)
summary(sum.model)
parameters(sum.model, exp = TRUE)


##### plots #####

#save model estimates
InteractWide[, PredictSum := (predict(sum.model.red1, type = "response"))]

plotSums = dcast(InteractWide, ID +Pen + scaleElos + Group_Size ~ Situation, value.var = c("PredictSum", "Sum"))
plotSums[, Sum_HQ_all := (Sum_HQ_all/40)*30]
plotSums[, Sum := Sum_HQ_all + Sum_Feed_all + Sum_Normal_all]
plotSums[, PredictSum := PredictSum_HQ_all + PredictSum_Feed_all + PredictSum_Normal_all]




#effect plot of Elo and number of interactions split by group
eloPlot = ggplot(data = plotSums, mapping = aes(x = scaleElos, y =Sum, colour = Pen)) + 
  geom_smooth(aes(x = scaleElos, y = PredictSum), se= FALSE)+#method = glmer.nb, formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_point(size = 2.5) + 
  labs(x = 'Scaled Elo rating', y = 'Sum of interactions')+
  facet_grid(.~ Group_Size,labeller = labeller(Group_Size = c("large" = "large groups", "small" = "small groups"))) + 
  theme_bw(base_size = 18)+
  scale_color_manual(values=c(largeCol, smallCol))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


#plot effects of the condition
plotData <- as.data.table(emmeans(sum.model.red1, ~ pairwise ~Situation*Group_Size, type = "response", offset = log(InteractWide$Minutes))$emmeans)
plotData2 = InteractWide
plotData2[Situation == "HQ_all", Sum := (Sum/40)*30]

#not used in paper
ggplot(plotData, aes(x = Situation, y = response))+
  geom_jitter(data = plotData2, aes(x = Situation, y= Sum, colour = Pen), width=0.27)+
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL, shape = "model estimate"), 
                  size = 1, linewidth = 1.5, colour = "black", 
                  fill="yellow")+
  scale_color_manual(values=c(largeCol, smallCol))+
  facet_grid(.~ Group_Size,labeller = labeller(Group_Size = c("large" = "large groups", "small" = "small groups"))) + 
  labs(x = "Resource condition", y= "Sum of interactions")+
  scale_shape_manual(name = NULL, values = c("model estimate" = 23))+
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2))+
  scale_x_discrete(labels=c("HQ", "FC", "CON"))+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position = "top")

plotData <- as.data.table(emmeans(sum.model.red1, ~ pairwise ~Situation, type = "response", offset = log(InteractWide$Minutes))$emmeans)

#plot used in paper
condPlot = 
  ggplot(plotData, aes(x = Situation, y = response))+
  geom_pointrange(data = plotData2[, .(mean = mean(Sum), 
                                     q1 = quantile(Sum, 0.25),
                                     q2 = quantile(Sum, 0.75)),by = Situation], 
             aes(x = Situation, y= mean, ymin = q1, ymax = q2, 
                 colour = "observed data"),
             size = 1,linewidth = 1.5,
             position = position_nudge(x = -0.025))+
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL,
                      colour = "model estimates"), 
                  size = 1, linewidth = 1.5,
                  position = position_nudge(x = 0.025))+
  scale_color_manual(name = NULL, values=c("observed data" = "grey", "model estimates" = "black"))+
  labs(x = "Resource condition", y= "Sum of interactions")+
  scale_x_discrete(labels=c("HQ", "FC", "CON"))+
  theme_classic(base_size = 18)


figInteractions = ggarrange(eloPlot,condPlot, ncol = 2, labels = c("a)", "b)"), 
                             font.label=list(color="black",size=16),widths = c(1,0.7), legend = "right")

ggsave(path = "plots", "SumInteractions.tiff", figInteractions, "tiff", width = 35, height= 18, units = "cm", dpi = 300)


#### Aggression Intensity #################

#descriptors
descript = RankTable[, .(Sum = sum(AggressBool), N = .N, perc = sum(AggressBool)/.N), by = Pen]


#create dyad identifier
RankTable[, Dyad := paste(sort(c(Winner, Loser)), collapse = " "), by = 1:nrow(RankTable)]
RankTable[, Dyad := paste(Pen, Dyad)]


##### statistical model #####

#to change level orders
RankTable[, Situation := factor(Situation, levels = c("HQ", "Feeder", "Normal"))]
RankTable[, Group_Size := factor(Group_Size, levels = c("small", "large"))]

# non-physical set as reference 0
#binomial model for intensity of aggression
intensity.model = glmer( AggressBool ~ WinnerElo+ LoserElo + Situation+ Group_Size + 
                           WinnerElo:Group_Size + Situation:WinnerElo + Group_Size:Situation +
                           LoserElo:Group_Size + LoserElo:WinnerElo + LoserElo:Situation + (1|Dyad/Pen), data = RankTable, family = binomial)
#failed to converge -> Pen barely any variation -> try without
intensity.model = glmer( AggressBool ~ WinnerElo+ LoserElo + Situation+ Group_Size + 
                           WinnerElo:Group_Size + Situation:WinnerElo + Group_Size:Situation +
                           LoserElo:Group_Size + LoserElo:WinnerElo + LoserElo:Situation + (1|Dyad), 
                         data = RankTable, family = binomial, 
                         glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 100000)))

resid.intensity = simulateResiduals(intensity.model)
plot(resid.intensity)
plotResiduals(resid.intensity, form = RankTable$DiffElo) #nearly perfect
plotResiduals(resid.intensity, form = RankTable$Group_Size) #good
plotResiduals(resid.intensity, form = RankTable$Situation) #good

intensity.null = glmer( AggressBool ~ 1 + (1|Dyad), data = RankTable, family = binomial)
anova(intensity.model, intensity.null) #better than null model

#reduce complexity: test 2-way interactions
drop1(intensity.model) #keep situation*Groupsize & Winner*Loser


intensity.model = glmer(AggressBool ~ WinnerElo+ LoserElo + Situation+ Group_Size + 
                                         Group_Size:Situation + LoserElo:WinnerElo + (1|Dyad), 
                        data = RankTable, family = binomial, 
                        glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
resid.intensity = simulateResiduals(intensity.model)
plot(resid.intensity)
plotResiduals(resid.intensity, form = RankTable$DiffElo) #nearly perfect
plotResiduals(resid.intensity, form = RankTable$Group_Size) #good
plotResiduals(resid.intensity, form = RankTable$Situation) #good

anova(intensity.model, intensity.null) #better than null model

#model estimates
summary(intensity.model)
parameters(intensity.model, exponentiate = TRUE)
tab_model(intensity.model)
plot(predictorEffects(intensity.model), lines=list(multiline=TRUE))
plot(allEffects(intensity.model))
plot_model(intensity.model, show.value = T)

#explained variance
r.squaredGLMM(intensity.model, intensity.null)

#estimated means and post hoc testing

#aggression intensity by winnerelo fixed at min, mean, or max loserelo
emmeans(intensity.model, ~ WinnerElo, 
        at = list(WinnerElo = c(min(RankTable$WinnerElo),mean(RankTable$WinnerElo),max(RankTable$WinnerElo))), type= "response")

winBeta = summary(emtrends(intensity.model, ~ LoserElo, var="WinnerElo", 
                   at=list(LoserElo=c( min(RankTable$LoserElo),
                                       mean(RankTable$LoserElo), 
                                       max(RankTable$LoserElo)))), adjust = "bonferroni", level = 0.95, type = "response")
#if odds ratio is needed exponentiate instead of probability
#winBeta$WinnerElo.trend = exp(winBeta$WinnerElo.trend)
#winBeta$asymp.LCL = exp(winBeta$asymp.LCL)
#winBeta$asymp.UCL = exp(winBeta$asymp.UCL)
winBeta


#aggression intensity by group size and situation    
group_sit = emmeans(intensity.model, ~ pairwise ~ Group_Size*Situation, type = "response")
confint(group_sit, adjust = "bonferroni", level = 0.95)


##### plots ####

#Aggression intensity by WinnerElo*LoserElo

#change to show interaction between winner and loser elo at fix loserElo and range of winner
plotData = as.data.table(summary(emmeans(intensity.model, ~ WinnerElo*LoserElo, 
                          at = list(WinnerElo = seq(from = min(Interact$scaleElo), to = max(Interact$scaleElo), 
                                             by = 0.1),
                             LoserElo = c(min(Interact$scaleElo),mean(Interact$scaleElo),max(Interact$scaleElo))), type= "response")))

plotData[, LoserCat := "mean" ]
plotData[LoserElo %in% min(unique(LoserElo)), LoserCat := "min" ]
plotData[LoserElo %in% max(unique(LoserElo)), LoserCat := "max" ]
plotData[, LoserCat := as.factor(LoserCat)]

#change level order
#RankTable[, Situation := factor(Situation, levels = c("Feeder", "HQ","Normal"))]


intensPlot1 = 
  ggplot(plotData, aes(x = WinnerElo, y = prob)) +
  geom_jitter(data = RankTable, aes(y = AggressBool, colour = Pen),
              height = 0.02, alpha = 0.5)+
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL, fill = LoserCat), alpha = 0.6) +
  geom_line(aes(linetype = LoserCat), linewidth = 1)+
  theme_classic(base_size = 18)+
  labs(x = 'Elo rating of Winners', y = "Probability for high intensity aggression",
       )+
  #scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA])) 
  #                   labels = c("25", "35", "45" , "55"))+
  scale_color_manual( values=c(largeCol, smallCol))+
  scale_fill_manual(name = "Elo rating \nof losers", values=c("grey70", "grey40", "grey90"))+
  scale_linetype_manual(name = "Elo rating \nof losers", values = c("min" = "dotted", "mean" = "solid", "max" = "dashed")) +
  guides(color = guide_legend(order = 1), "Elo rating \nof losers" = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))


#Aggression intensity by Situation and Group Size

plotData <- as.data.table(emmeans(intensity.model, ~ pairwise ~Situation*Group_Size, type = "response")$emmeans)

#get 95% CI for bernoulli data 
plotData2 = RankTable[, .(mean = mean(AggressBool), 
              CI1 = binom.confint(sum(AggressBool), 
                                  length(AggressBool), 
                                  conf.level = 0.95, method = 'logit')$lower,
              CI2 = binom.confint(sum(AggressBool), 
                                  length(AggressBool), 
                                  conf.level = 0.95, method = 'logit')$upper),
          by = .(Situation, Group_Size)]


#change level order
plotData[, Group_Size := factor(Group_Size, levels = c("large", "small"))]
plotData[, Situation := factor(Situation, levels = c("HQ","Feeder","Normal"))]
plotData2[, Group_Size := factor(Group_Size, levels = c("large", "small"))]
plotData2[, Situation := factor(Situation, levels = c("HQ","Feeder","Normal"))]

intensPlot2 = 
  ggplot(plotData, aes(x = Situation, y = prob))+
  geom_pointrange(data = plotData2, 
                  aes(x = Situation, y= mean, ymin = CI1, ymax = CI2, 
                      colour = "observed data"),
                  size = 1,linewidth = 1.5,
                  position = position_nudge(x = -0.03))+
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL,
                      colour = "model estimates"), 
                  size = 1, linewidth = 1.5,
                  position = position_nudge(x = 0.03))+
  scale_color_manual(name = NULL, values=c("observed data" = "grey", "model estimates" = "black"))+
  labs(x = "Resource Condition", y= "Probability for high intensity aggression")+
  facet_grid(.~Group_Size, labeller = labeller(Group_Size = c("large" = "large groups", "small" = "small groups")))+
  scale_x_discrete(labels=c("HQ", "FC", "CON"))+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position = "top")

figIntensity = ggarrange(intensPlot1,intensPlot2, ncol = 2, labels = c("a)", "b)"), 
                            font.label=list(color="black",size=16))

ggsave(path = "plots", "Intensity2.tiff", figIntensity, "tiff", width = 38, height= 15, units = "cm", dpi = 300)


#### Patterns of aggression ############

#strategy
diagnA$strategy
diagnB$strategy
diagnC$strategy
diagnD$strategy
diagnE$strategy
diagnF$strategy

#focus & position
diagnA$focus_pos
diagnB$focus_pos
diagnC$focus_pos
diagnD$focus_pos
diagnE$focus_pos
diagnF$focus_pos

####

#strategy plots
# Two rows, two columns
par(mfrow = c(2, 3))

# Plots

adj_dom_plot_strategy(diagnA$focus_pos, diagnA$blur, show_data_ci = T)
title("Pen A")
legend("topleft", legend = c("Downward heuristic","Bullying","Close competitors","Undefined"),
       col = "black",
       pt.bg = c(
         scales::alpha("white", 0.2),
         scales::alpha("blue", 0.2),
         scales::alpha("red", 0.2),
         scales::alpha("grey", 0.4)
       ),
       pch = 22, cex = 1, pt.cex = 2)

adj_dom_plot_strategy(diagnB$focus_pos, diagnB$blur, show_data_ci = T)
title("Pen B")
adj_dom_plot_strategy(diagnC$focus_pos, diagnC$blur, show_data_ci = T)
title("Pen C")
adj_dom_plot_strategy(diagnD$focus_pos, diagnD$blur, show_data_ci = T)
title("Pen D")
adj_dom_plot_strategy(diagnE$focus_pos, diagnE$blur, show_data_ci = T)
title("Pen E")
adj_dom_plot_strategy(diagnF$focus_pos, diagnF$blur, show_data_ci = T)
title("Pen F")

# Back to the original graphics device
par(mfrow = c(1, 1))

# Plot for Dynamics of interactions by rank (large)
largePattern = 
  ggplot(data = RankTable[Group_Size == "large",], aes(x = WinnerRank, y =LoserRank)) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  #geom_point(size = 2.5) + 
  labs(y = 'Loser rank', x= 'Winner rank')+
  theme_classic(base_size = 18)+
  facet_grid(. ~ Pen)+
  geom_density_2d_filled(alpha = 0.7, contour_var = "ndensity", show.legend = FALSE, bins = 4)+
  scale_fill_manual(values = c("white", "#FFEBCD", "gold", "darkorange"))+
  geom_abline(intercept = 0 , slope = 1, linetype = "dashed", colour = 'grey')+
  geom_jitter(size = 0.6, width = 0.4)+
  xlim(0.5,120)+
  ylim(0.5,120)



# Plot for Dynamics of interactions by rank (small)
smallPattern = 
  ggplot(data = RankTable[Group_Size == "small",], mapping = aes(x = WinnerRank, y =LoserRank, shape = "observations")) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  theme_classic(base_size = 18)+
  labs(y = 'Loser rank', x= 'Winner rank', shape = NULL, fill = "density \nrange")+
  theme(legend.position = "top")+
  facet_grid(. ~ Pen)+
  geom_density_2d_filled(alpha = 0.7, contour_var = "ndensity",  bins = 4)+
  scale_fill_manual(values = c("white", "#FFEBCD", "gold", "darkorange"))+
  geom_abline(intercept = 0 , slope = 1, linetype = "dashed", colour = 'grey')+
  geom_jitter(size = 0.8, width = 0.4)+
  xlim(0.5,20.5)+
  ylim(0.5,20.5)




# Example plot for Dynamics of interactions by rank 
ggplot(data = RankTable[Pen == "E",], mapping = aes(x = WinnerRank, y =LoserRank)) + 
  #geom_smooth(se= FALSE)+#method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
  labs(y = 'Loser rank', x= 'Winner rank', fill = "density \nrange")+
  theme_classic(base_size = 18)+
  geom_density_2d_filled(alpha = 0.7, contour_var = "ndensity", bins = 4)+
  scale_fill_manual(values = c("white", "#FFEBCD", "gold", "darkorange"))+
  geom_abline(intercept = 0 , slope = 1, linetype = "dashed", colour = 'grey')+
  geom_jitter(size = 1.5, width = 0.4)+
  xlim(0.5,20.5)+
  ylim(0.5,20.5)

legend <- get_legend(smallPattern)

figPattern = ggarrange(largePattern,smallPattern, nrow = 2, 
                            font.label=list(color="black",size=16),widths = c(1,0.7), legend = "top", 
                       align = "hv", legend.grob = legend)

ggsave(path = "plots", "Pattern2.tiff", figPattern, "tiff", width = 25, height= 23, units = "cm", dpi = 300)


#### Individual contacts data ######################################

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

#### Plots ##############
##### Examples of steepness plots ##############

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

mean.eloM = rowMeans(ratingE$randElo)
identitiesM = rownames(ratingE$randElo)
identitiesM <- identitiesM[order(mean.eloM)]
CIsM <- apply(ratingE$randElo,1,quantile,c(0.025,0.975),na.rm=TRUE)
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
  labs(x = "Rank order", y="Elo rating")+
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
  #scale_y_continuous(breaks = c(-400, 0, 400))+
  #geom_text(
  #  label=plotTable$IDs,
  #  size = 4)+
  labs(x = "Rank order", y="Elo rating")+
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 


plotTable1[, Group_Size := "small"]
plotTable2[, Group_Size := "large"]

plotTable = rbind(plotTable2, plotTable1)
plotTable$Group_Size = factor(plotTable$Group_Size, levels=c("large","small"), labels = c( "large group","small group"))
BCol = "#2F5597"
FCol = "#E61C45"

pp = ggplot(plotTable, aes(x = number, y= ranks, colour = Group_Size) )+
  geom_errorbar(aes(x= number, ymin=CIl, ymax=CIu), width=.1)+
  geom_point( size = 2)+
  scale_y_continuous(breaks = c(-600, -300, 0, 300, 600))+
  #geom_text(
  #  label=plotTable$IDs,
  #  size = 4)+
  facet_grid( . ~ Group_Size, scales = "free")+
  scale_color_manual(values=c(BCol, FCol))+
  labs(x = "Rank order", y="Elo rating")+
  #scale_y_continuous(breaks=c(-400, 700),
  #                   labels=c("low rank", "high rank"))+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) 
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
