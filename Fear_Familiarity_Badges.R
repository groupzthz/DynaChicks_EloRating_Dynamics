#clear workspace
rm(list = ls())


library(data.table)
library(ggplot2)
library(MuMIn)
options(na.action = "na.fail")
library(lme4)
library(DHARMa)
library(parameters)
library(emmeans)
library(effects)
library(RColorBrewer) # color for plotting
library(MASS) #for negative binomial glm
library(car)
set.seed(1)


#### Loading and preparing Data ###########
arenaTest1 = fread("ArenaTest1.csv")
arenaTest2 = fread("ArenaTest2.csv")
Individ = fread("IndividualDataTests.csv")
Individ[, scaleElo:= as.numeric(scale(elos)), by = Pen]
Individ[, Unique:= paste0(Pen, ID)]
Thermal = fread("ThermalResults.csv")

#merging data
tempData = merge(arenaTest1[, !c("Date", "Time", "Pophole", "Defaecation", "Behaviours", "Comb", "Weight", "Claws", "Notes", "exclusionPot"), with = F], 
                 Individ[, !c("Focals"), with = F], 
                 by = c("Pen", "ID"))

fullData = merge(arenaTest2[,!c("Date", "Start", "Defaecation24", "Behaviours", "Claws26","Comb26", "Weight26", "Notes26", "exclusionPot"), with = F ], 
                 tempData, by = c("Pen", "ID"))
fullData = merge(Thermal, fullData, by = c("Pen", "ID"))

fullData[, "Foot26" := LeftFoot26 + RightFoot26]
fullData[, "Foot" := LeftFoot + RightFoot]
fullData[, "GroupBin" := ifelse(Group_Size == "Large", 1, 0)]
fullData[, FamiliarElo := Individ$scaleElo[match(FamiliarID, Individ$Unique)]]
fullData[, StrangerElo := Individ$scaleElo[match(StrangerID, Individ$Unique)]]
fullData[, Group_Size := as.factor(Group_Size)]
fullData[, Offset := 600]

#change label names for plotting with facets
fullData$plotGroup_Size <- factor(fullData$Group_Size, levels = c("large", "small"),
                                  labels = c("large groups", "small groups"))


#Careful: AKO has no Elo because of wrong labels during observations

#colours for plots
largeCol = brewer.pal(n = 8, name = "Blues")[c(4,6,8)]
smallCol = brewer.pal(n = 8, name = "OrRd")[c(4,6,8)]

#### Focal individuals ####

fullData[, list(Loser = sum(scaleElo <0), Winner = sum(scaleElo >0)), by = c("Group_Size")]
# -> more subordinate tendency animals than dominant... especially in large groups

#distribution of Elo ranking by Pen/Group_Size

#Selected individuals in each pen
ggplot(data = fullData, aes(x = Pen, y = scaleElo, fill = Pen, colour = Pen))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = Individ, aes(x = Pen, y = scaleElo, shape = "non focal"), colour = "black", alpha = 0.5, size = 1)+
  geom_violin(alpha = 0.3)+
  geom_point(alpha = 0.8, size = 3)+
  theme_classic(base_size = 18)+
  scale_fill_manual(values=c(largeCol, smallCol))+
  scale_colour_manual(values=c(largeCol, smallCol))+
  #scale_fill_manual(values = c("orange", "lightblue"))+
  labs(y = "scaled Elo rating", fill = "Pen")+
  scale_shape_manual(name = NULL, values = c("non focal" = 4))+ 
  guides(fill = "none",
         colour = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  coord_flip()+
  scale_x_discrete(limits=rev)

#distribution of Elo ratings by group size
ggplot(data = fullData, aes(x = Group_Size, y = scaleElo))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha = 0.7, size = 3)+
  theme_classic(base_size = 18)



#### Badges of status ########


##### model statistics ########

#evaluated on the weight and comb size measured at 26 WoA 

#Effect of badges of status on Elo
badges_model = lmer(scaleElo ~ scale(Weight2)*scale(CombSize2)+ scale(CombSize2)*Group_Size + scale(Weight2)*Group_Size + (1|Pen), data = fullData)
#singularity warning
badges_model = lm(scaleElo ~ scale(Weight2)*scale(CombSize2)+ scale(CombSize2)*Group_Size + scale(Weight2)*Group_Size, data = fullData)
resid.Badges<- simulateResiduals(badges_model, 1000)
plot(resid.Badges)
plotResiduals(resid.Badges, form = fullData$Group_Size) #heteroscedacity? still ok
plotResiduals(resid.Badges, form = scale(fullData$CombSize2))
plotResiduals(resid.Badges, form = scale(fullData$Weight2)) 

vif(badges_model) # not too high co-linear (values <5)

drop1(badges_model) #better without interactions

badges_model.red = lm(scaleElo ~ scale(Weight2)+scale(CombSize2)+Group_Size, data = fullData)
resid.Badges<- simulateResiduals(badges_model.red, 1000)
plot(resid.Badges)
plotResiduals(resid.Badges, form = fullData$Group_Size) #heteroscedacity? still ok
plotResiduals(resid.Badges, form = scale(fullData$CombSize2)) 
plotResiduals(resid.Badges, form = scale(fullData$Weight2)) 

badges_model_null = lm(scaleElo ~ 1, data = fullData)

AIC( badges_model.red, badges_model_null)
anova(  badges_model_null, badges_model.red)
r.squaredGLMM(badges_model.red, badges_model_null)

summary(badges_model.red)


##### plots ########

#Comb size and Elo rating
ggplot(data= fullData, aes(y = scaleElo)) + 
  geom_point(aes( x = CombSize2))+ 
  geom_smooth(aes( x = CombSize2), method = lm, formula = y~x, colour = "black") + 
  geom_point(aes( x = CombSize1))+ 
  geom_smooth(aes( x = CombSize1), method = lm, formula = y~x) + 
  facet_grid(.~Group_Size)+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Comb Size (cm²)", y = "scaled Elo rating")

#weight and Elo rating
ggplot(data= fullData, aes(y = scaleElo)) + 
  geom_point(aes( x = Weight2))+ 
  geom_smooth(aes( x = Weight2), method = lm, formula = y~x, colour = "black") + 
  geom_point(aes( x = Weight1))+ 
  geom_smooth(aes( x = Weight1), method = lm, formula = y~x) +
  facet_grid(.~plotGroup_Size)+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Weight (g)", y = "scaled Elo rating")


#are weight and comb size correlated?
cor.test(fullData[,CombSize2], fullData[,Weight2]) # not very high r = 0.28
ggplot(data= fullData, aes(y = Weight2, x = CombSize2)) + 
  geom_point()+ 
  geom_smooth(method = lm, formula = y~x) + 
  facet_grid(.~Group_Size)


###
#Effect of badges of the number of Interactions
ggplot(data= fullData, aes(y = sum, x = CombSize2)) + 
  geom_point()+ 
  geom_smooth(method = lm, formula = y~x) + 
  facet_grid(.~plotGroup_Size)+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Comb Size (cm²)", y = "Number of interactions")

ggplot(data= fullData, aes(y = sum, x = Weight2)) + 
  geom_point()+ 
  geom_smooth(method = lm, formula = y~x) + 
  facet_grid(.~Group_Size)


#### Fear test ##########

###### LATENCY ########


Lat.model = lmer( Box ~ scaleElo*Group_Size + offset(Offset) + (1|Pen), data = fullData)
resid.Lat = simulateResiduals(Lat.model, 1000)
plot(resid.Lat)#no good fit

Lat.model = glmer( Box ~ scaleElo*Group_Size + offset(log(Offset))+ (1|Pen), family = poisson, data = fullData)
resid.Lat = simulateResiduals(Lat.model, 1000)
plot(resid.Lat) #even worse

Lat.model = glmer.nb( Box ~ scaleElo*Group_Size + offset(log(Offset))+ (1|Pen), data = fullData)
#singularity fit

Lat.model = glm.nb( Box ~ scaleElo*Group_Size + offset(log(Offset)), data = fullData)
resid.Lat = simulateResiduals(Lat.model, 1000)
plot(resid.Lat)
plotResiduals(resid.Lat, form = fullData$scaleElo)
plotResiduals(resid.Lat, form = fullData$Group_Size)

drop1(Lat.model)

Lat.model.red = glm.nb( Box ~ scaleElo+Group_Size + offset(log(Offset)), data = fullData)
resid.Lat = simulateResiduals(Lat.model.red, 1000)
plot(resid.Lat)
plotResiduals(resid.Lat, form = fullData$scaleElo)
plotResiduals(resid.Lat, form = fullData$Group_Size)

Lat.model.null = glm.nb( Box ~ 1+ offset(log(Offset)), data = fullData)#TODO: why AIC different?
AIC(Lat.model.red, Lat.model.null) #worse than null
anova(Lat.model.null, Lat.model.red, test = "Chisq")
r.squaredGLMM(Lat.model.red, Lat.model.null)


#Latency to leave box by group size and rank
ggplot(data = fullData, aes(colour = Group_Size, y = Box, x = scaleElo))+
  geom_point(size = 2)+
  geom_smooth(method = "lm", formula = y~x, size = 1.5)+
  theme_bw(base_size = 16)+
  scale_color_manual(values=c(largeCol[2], smallCol[3]))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Latency to leave dived by full time given
ggplot(fullData, aes(x = scaleElo, y = Box/Offset))+
  geom_point()+
  geom_smooth(method="glm",method.args=list(family="binomial"), formula = y~x)+
  facet_grid(.~Group_Size)+
  theme_bw(base_size = 16)

# #testing via survival analysis ... old
# library(survival)
# library(survminer)
# #library(dplyr)
# library(coxme)
# #treating non-successfull emergence as censored data = 0
# # censoring occurs if subject has not experienced the event of interest by the end of collection
# fullData[, Emergence := ifelse(Box == 600, 0, 1)]
# surv_object1 <- Surv(time = fullData$Box, event = fullData$Emergence)
# 
# 
# #full model
# fit1 = coxme(surv_object1 ~ scaleElo*Group_Size+ (1|Pen), data = fullData)
# fit2 <- coxph(surv_object1 ~ scaleElo*Group_Size, data=fullData)
# anova(fit1, fit2) # no grouping necessary
# fit3 <- coxph(surv_object1 ~ scaleElo+Group_Size, data=fullData)
# anova(fit2, fit3) # no interaction necessary
# fitnull = coxme(surv_object1 ~ 1+ (1|Pen), data = fullData)
# 
# anova(fit3, fitnull, test = "Chisq") # not better than null model
# fittest = cox.zph(fit3)
# ggcoxzph(fittest)
# ggcoxdiagnostics(fit3,
#                  type = "deviance",
#                  ox.scale = "linear.predictions")
# ggforest(fit3)
# ggadjustedcurves(fit3, data=fullData, variable = "Group_Size")

# ggsurvplot(
#   fit = survfit(surv_object1 ~ Elo_cut, data = fullData), 
#   xlab = "Latency to emerge", 
#   ylab = "Overall emergence probability")


###### FIRST CHOICE #######
fullData[FirstChoice.y == "", FirstChoice.y := "No Choice"]
comp.Choice = fullData[, .N, by = c("FirstChoice.y", "Side", "Group_Size")]

ggplot(data = comp.Choice, aes(x = FirstChoice.y, y = N, fill = Side))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)+
  facet_grid(. ~Group_Size)

fullData[, Side := as.factor(Side)]


# Side bias:
choiceData1 = fullData[FirstChoice.y != "No Choice",]
choiceData1[Side == "A", Side_Social := "Right"]
choiceData1[Side == "S", Side_Social := "Left"]
choiceData1[Side == "A", ChoiceDirection := ifelse(FirstChoice.y == "Social", "Right", "Left")]
choiceData1[Side == "S", ChoiceDirection := ifelse(FirstChoice.y == "Social", "Left", "Right")]
choiceData1[,ChoiceDirection := as.factor(ChoiceDirection)]
choiceData1[,Side_Social := as.factor(Side_Social)]
choiceData1[,ChoiceDirectionBool := ifelse(ChoiceDirection == "Right", 0, 1)]

comp.Choice = choiceData1[, .N, by = c("ChoiceDirection", "Side_Social", "Group_Size")]
ggplot(data = comp.Choice, aes(x = ChoiceDirection, y = N, fill = Side_Social))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)+
  facet_grid(. ~Group_Size)

ggplot(comp.Choice, aes(x = Group_Size, y = N,  fill = ChoiceDirection))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)

#right set as 1 left as 0
choice.model1 = glmer(ChoiceDirectionBool ~ scaleElo*Group_Size+ Group_Size*Side_Social + scaleElo*Side_Social + (1|Pen), data = choiceData1, family = "binomial")
#singularity
choice.model1 = glm(ChoiceDirectionBool ~ scaleElo*Group_Size+ Group_Size*Side_Social + scaleElo*Side_Social, data = choiceData1, family = "binomial")
resid.choice<- simulateResiduals(choice.model1, 1000)
plot(resid.choice)
plotResiduals(resid.choice, form = choiceData1$Side_Social)
plotResiduals(resid.choice, form = choiceData1$scaleElo) #problematic?
plotResiduals(resid.choice, form = choiceData1$Group_Size)

#take out 2-way
drop1(choice.model1)

choice.model1.red = glm(ChoiceDirectionBool ~ Group_Size + Side_Social + scaleElo, data = choiceData1, family = "binomial")
resid.choice<- simulateResiduals(choice.model1.red, 1000)
plot(resid.choice)
plotResiduals(resid.choice, form = choiceData1$Side_Social)
plotResiduals(resid.choice, form = choiceData1$scaleElo) 
plotResiduals(resid.choice, form = choiceData1$Group_Size)

choice.model1.null = glm(ChoiceDirectionBool ~ 1, data = choiceData1, family = "binomial")
AIC(choice.model1.null, choice.model1.red)

summary(choice.model1.red)
parameters(choice.model1.red, exponentiate = TRUE)
r.squaredGLMM(choice.model1.red, choice.model1.null)

emmeans(choice.model.red, ~ Group_Size + Side_Social + scaleElo, type = "response")
emmeans(choice.model.red, ~ Side_Social, type = "response")
emmeans(choice.model.red, ~ Group_Size, type = "response")





###### FEEDING ######
ggplot(fullData, aes(x = scaleElo, y = Feeding))+ 
  geom_point(aes(colour = Group_Size),size = 3)+
  geom_smooth(method="glm",method.args=list(family="binomial"), formula = y~x, colour = "black")+
  #geom_point(data = fullData[Emergence == 0,], colour = "red", size = 3)+
  theme_bw(base_size = 16)

#Latency to leave box by group size and rank
ggplot(data = fullData, aes(colour = Group_Size, y = Feeding, x = scaleElo))+
  geom_point(size = 2)+
  geom_smooth(method="glm",method.args=list(family="binomial"), formula = y~x, size = 1.5)+
  theme_bw(base_size = 16)+
  scale_color_manual(values=c(largeCol[2], smallCol[3]))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(fullData, aes(x = scaleElo, y = Feed/Offset))+
  geom_point()+
  geom_smooth(method="glm",method.args=list(family="binomial"), formula = y~x)+
  facet_grid(.~Group_Size)+
  geom_point(data = fullData[Feeding == 1,], colour = "red", size = 3)+
  theme_bw(base_size = 16)

ggplot(fullData, aes(x = Group_Size, y = Feed))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(data = fullData[Feeding == 0,], colour = "black", size = 3)+
  geom_jitter(data = fullData[Feeding == 1,], colour = "red", size = 3)+ # careful, double plotted because of jitter
  theme_bw(base_size = 16)

Feeding.model = glmer(Feeding ~ scaleElo*Group_Size +(1|Pen), data = fullData, family = binomial)
#singularity
Feeding.model = glm(Feeding ~ scaleElo*Group_Size, data = fullData, family = binomial)
resid.Feeding = simulateResiduals(Feeding.model, 1000)
plot(resid.Feeding)
plotResiduals(resid.Feeding, form = fullData$Group_Size)
plotResiduals(resid.Feeding, form = fullData$scaleElo)

#take out 2-way
drop1(Feeding.model)

Feeding.model.red = glm(Feeding ~ scaleElo + Group_Size, data = fullData, family = binomial)
resid.Feeding<- simulateResiduals(Feeding.model.red, 1000)
plot(resid.Feeding)
plotResiduals(resid.Feeding, form = fullData$Group_Size)
plotResiduals(resid.Feeding, form = fullData$scaleElo)

Feeding.null = glm(Feeding ~ 1, data = fullData, family = binomial)
AIC(Feeding.null, Feeding.model)# only trend to be better than null??
anova(Feeding.null, Feeding.model, test = "Chisq")

summary(Feeding.model.red)
parameters(Feeding.model.red, exponentiate = T)
r.squaredGLMM(Feeding.model.red, Feeding.null)

### plot

plotData = summary(emmeans(Feeding.model.red, ~scaleElo, at = list(scaleElo = seq(from = min(fullData$scaleElo),
                                                               to = max(fullData$scaleElo), by = 0.1)), type = "response"))

ggplot(plotData, aes(x = scaleElo, y = prob))+
  geom_jitter(data = fullData, aes(y = Feeding, colour = Pen), height = 0.02,size = 2.5)+
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), fill = "grey70", alpha = 0.7) +
  geom_line(linewidth = 1)+
  theme_classic(base_size = 18)+
  labs(x = 'scaled Elo ratings', y = "Probability for feeding in the arena")+
  #scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
  #                   labels = c("25", "35", "45" , "55"))+
  scale_color_manual(values=c(largeCol, smallCol))+
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))


###### SOCIAL TIME #####
#PROBLEM WITH THE DATA:
# zero not defined for geometric mean = exp(mean(log(x))) or transformation
# arena 1; zero entries are ignored pulling the means closer to each other than they actually are 
# auf social family binomial/neg binom/beta?

#select all individuals
# Modelling
modelData = fullData[, .(Pen, Group_Size, scaleElo, 
                         Social, Feed, Between, Box,
                         Familiar, Stranger, Between2, FirstEntry,
                         Offset)]

hist(modelData$Social) # many zeros, pattern rather negbinom than poisson
modelData$Social # no value 600, so no censoring needed

model.Social.FULL = lmer(Social ~ Group_Size*scaleElo + (1|Pen),
                   data = modelData) 
resid.Social<- simulateResiduals(model.Social.FULL, 1000)
plot(resid.Social)
testZeroInflation(resid.Social) # zero-inflation

Social.glmm <- glmmTMB(formula = Social ~ scaleElo*Group_Size + (1|Pen), 
                     REML = F,        # needs to be stated since default = ML
                     ziformula = ~1, 
                     data = modelData) 
resid.test = simulateResiduals(test.glmm, 1000)
plot(resid.test)
plotResiduals(resid.test, form = fullData$scaleElo)
plotResiduals(resid.test, form = fullData$Group_Size)

drop1(Social.glmm, test = "Chisq")
Social.glmm <- glmmTMB(formula = Social ~ scaleElo+Group_Size + (1|Pen), 
                       REML = F,        # needs to be stated since default = ML
                       ziformula = ~1, 
                       data = modelData) 
resid.Social<- simulateResiduals(Social.glmm, 1000)
plot(resid.Social)
plotResiduals(resid.Social, form = fullData$scaleElo)
plotResiduals(resid.Social, form = fullData$Group_Size)
Social.glmm.null <- glmmTMB(formula = Social ~ 1 + (1|Pen), 
                       REML = F,        # needs to be stated since default = ML
                       ziformula = ~1, 
                       data = modelData) 

anova(Social.glmm, Social.glmm.null)
summary(Social.glmm) #not sigficantly different from zero


#test Feeding time
test = lmer(Feed ~ Group_Size*scaleElo + (1|Pen),
                         data = modelData)
resid.test<- simulateResiduals(test, 1000)
plot(resid.test)
plotResiduals(resid.test, form = modelData$Group_Size) 
plotResiduals(resid.test, form = modelData$scaleElo) 

library(glmmTMB)
#account for zeroinflation
test.glmm <- glmmTMB(formula = Feed ~ scaleElo*Group_Size + (1|Pen), 
                     REML = F,        # needs to be stated since default = ML
                     ziformula = ~1, 
                     data = modelData) 
resid.test = simulateResiduals(test.glmm, 1000)
plot(resid.test)
plotResiduals(resid.test, form = fullData$scaleElo)
plotResiduals(resid.test, form = fullData$Group_Size)

drop1(test.glmm, test = "Chisq")
test.glmm <- glmmTMB(formula = Feed ~ scaleElo+Group_Size + (1|Pen), 
                     REML = F,        # needs to be stated since default = ML
                     ziformula = ~1, 
                     data = modelData) 
resid.test = simulateResiduals(test.glmm, 1000)
plot(resid.test)
plotResiduals(resid.test, form = fullData$scaleElo)
plotResiduals(resid.test, form = fullData$Group_Size)
summary(test.glmm)

plotData = melt(modelData[, c(1:7)], 
                id.vars = c('Pen', 'Group_Size', "scaleElo"))

ggplot(data= plotData, aes(x = Group_Size, y = value, fill = variable))+
  geom_violin(draw_quantiles = c(0.5))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2))

#### Recognition test #########

###### FIRST CHOICE ######
#Side bias
comp.Choice = fullData[, .N, by = c("FirstChoice.x", "Side_Familiar", "Group_Size")]
ggplot(data = comp.Choice, aes(x = FirstChoice.x, y = N, fill = Side_Familiar))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)+
  facet_grid(. ~Group_Size)


# small group animals choose the left side (significant?) which is the hemisphere they use for recognition (Vallortigara)

fullData[, Side_Familiar := as.factor(Side_Familiar)]

choiceData2 = fullData[!is.na(FirstChoice.x),]
choiceData2[, FirstChoice.x := as.factor(FirstChoice.x)]
choiceData2[Side_Familiar == "Left", ChoiceDirection := ifelse(FirstChoice.x == "Stranger", "Right", "Left")]
choiceData2[Side_Familiar == "Right", ChoiceDirection := ifelse(FirstChoice.x == "Stranger", "Left", "Right")]
choiceData2[,ChoiceDirectionBool := ifelse(ChoiceDirection == "Right", 0, 1)]

#right is set as reference 0
choice.model = glmer(ChoiceDirectionBool ~ scaleElo*Group_Size + scaleElo*Side_Familiar + Group_Size*Side_Familiar + (1|Pen), data = choiceData2, family = "binomial")
#Singularity
choice.model = glm(ChoiceDirectionBool ~ scaleElo*Group_Size + scaleElo*Side_Familiar + Group_Size*Side_Familiar, data = choiceData2, family = "binomial")

resid.choice<- simulateResiduals(choice.model, 1000)
plot(resid.choice)
plotResiduals(resid.choice, form = choiceData2$Side_Familiar)
plotResiduals(resid.choice, form = choiceData2$Group_Size)
plotResiduals(resid.choice, form = choiceData2$scaleElo)

#drop 2-way
drop1(choice.model) #no clear indication that any interactions make the model better

choice.model.red = glm(ChoiceDirectionBool ~ scaleElo+ Side_Familiar + Group_Size, data = choiceData2, family = "binomial")
resid.choice<- simulateResiduals(choice.model.red, 1000)
plot(resid.choice) # problematic?
plotResiduals(resid.choice, form = choiceData2$Side_Familiar)
plotResiduals(resid.choice, form = choiceData2$Group_Size)
plotResiduals(resid.choice, form = choiceData2$scaleElo)

choice.model.null = glm(ChoiceDirectionBool ~ 1, data = choiceData2, family = "binomial")
AIC(choice.model.null, choice.model.red) # no clear difference
anova(choice.model.null, choice.model.red, test = "Chisq")

summary(choice.model.red)
parameters(choice.model.red, exponentiate = TRUE)
emmeans(choice.model.red, ~ Group_Size, type = "response")
vif(choice.model.red)
r.squaredGLMM(choice.model.red, choice.model.null)

#test
choice.model.red2 = glm(ChoiceDirectionBool ~ Group_Size, data = choiceData2, family = "binomial")
AIC(choice.model.red2, choice.model.null)

#plot preferential direction
comp.Choice = choiceData2[, .N, by = c("ChoiceDirection", "Side_Familiar", "Group_Size")]

comp.Choice$Group_Size <- factor(comp.Choice$Group_Size, levels = c("large", "small"),
                                 labels = c("large groups", "small groups")
)
ggplot(data = comp.Choice, aes(x = ChoiceDirection, y = N, fill = Side_Familiar))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)+
  labs(x = "First choice of direction", y = "Number of hens", fill = "Position of familiar")+
  facet_grid(. ~Group_Size)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_discrete(labels = c("Familiar to the left", "Familiar to the right"))

#first choice affects where to stay?
fullData[, DiffFam := scaleElo-FamiliarElo]
fullData[, DiffStr := scaleElo-StrangerElo]
cols = c('Pen', 'ID', 'Group_Size', 'Aggression_Str', 'Aggression_Fam', 'Familiar', 'Stranger', "FirstChoice.x", "scaleElo", "DiffFam", "DiffStr")
meltData = melt(fullData[,.SD, .SDcols = cols], id.vars = c('Pen', 'ID', 'Group_Size', 'Aggression_Str', 'Aggression_Fam', 'FirstChoice.x', "scaleElo", "DiffFam", "DiffStr"))
meltData[, Unique := paste0(Pen, ID)]
meltData$Group_Size <- factor(meltData$Group_Size, levels = c("large", "small"),
                              labels = c("large groups", "small groups")
)

meltData[, ratioDur := value/600]
ggplot(data = na.omit(meltData), aes(x = FirstChoice.x, y = value, fill = variable))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge())+
  facet_grid(.~Group_Size)+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "First choice of field", y = "Time spent in field", fill = "Field")

ggplot(data = na.omit(meltData), aes(x = FirstChoice.x, y = ratioDur, group = variable))+
  geom_point(position=position_dodge(width=0.4))+
  geom_line(aes(group = Unique), colour = "darkgrey")+
  facet_grid(.~Group_Size)+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "First choice of field", y = "Time spent in field", color = "Field")


###### FLOCKMATE TIME ######

#PROBLEM WITH THE DATA:
# zero not defined for geometric mean = exp(mean(log(x))) or transformation
# arena 2: too many entries containing zeros -> only 25 valid entries left

# Modelling
dat = fullData[, .(Pen, Group_Size, scaleElo, ID, 
                   Social, Feed, Between, Box,
                   Familiar, Stranger, Between2, FirstEntry,
                   Offset)]


#for Arena2 test take out animal that did not move at all because it has no first entry
# auf familiar cbind(success, (failure oder total)) family binomial/neg binom/beta?
modelData = dat[!is.na(FirstEntry),]

modelData$Familiar #no 600 so no censoring needed
hist(modelData$Familiar) 
# here comes the gaussian based on the distribution seen before:
# with pen as random:

modelData$Pen <- as.factor(modelData$Pen)
modelData$FirstEntry <- as.factor(modelData$FirstEntry)
model.Familiar = lmer(Familiar ~ Group_Size:scaleElo + Group_Size:FirstEntry + 
                                  scaleElo:FirstEntry + Group_Size + scaleElo + FirstEntry + (1|Pen),
                  data = modelData)
# singular fit issue, most likely due to Pen:
summary(model.Familiar)
# yes, indeed Pen is the trouble maker (no variance explained)
# therefore, lets try out an lm model:
model.Familiar = lm(Familiar ~ Group_Size:scaleElo + Group_Size:FirstEntry + 
                  scaleElo:FirstEntry + Group_Size + scaleElo + FirstEntry,
                data = modelData)
resid.test<- simulateResiduals(model.Familiar, 1000)
plot(resid.test)# very even distribution of residuals :-)
plotResiduals(resid.test, form = modelData$Group_Size) #good
plotResiduals(resid.test, form = modelData$scaleElo) # problem? not really
plotResiduals(resid.test, form = modelData$FirstEntry) # good

#take-out two-way
drop1(model.Familiar) # include GroupSize: FirstEntry

model.Familiar.red = lm(Familiar ~ Group_Size:FirstEntry +Group_Size+scaleElo+FirstEntry, data = modelData)
resid.model<- simulateResiduals(model.Familiar.red, 1000)
plot(resid.model) # good
plotResiduals(resid.model, form = modelData$FirstEntry)#good
plotResiduals(resid.model, form = modelData$Group_Size)#good
plotResiduals(resid.model, form = modelData$scaleElo)#bit weird still

model.Familiar.null = lm(Familiar ~ 1, data = modelData)
AIC(model.Familiar.null, model.Familiar.red) 

summary(model.Familiar.red)
parameters(model.Familiar.red)
r.squaredGLMM(model.Familiar.red, model.Familiar.null)


plot(allEffects(model.Familiar.red))
library(lmerTest)
anova(model.test.3)

emmeans(model.Familiar.red, pairwise ~   FirstEntry*Group_Size) #this is needed to get all comparisons
#do not show p-values from contrasts (no correction!! also only trend too weak for post-hoc because of power)
#mach plots in greyscale



boxplot(modelData$Familiar~modelData$FirstEntry)
# to be on safe side:
interGF<-interaction(modelData$Group_Size, modelData$FirstEntry)
kruskal.test(Familiar ~ interGF,data=modelData)
# still seems to be a trend in the interaction and mainly due to FirstEntry

table(modelData$Pen,modelData$FirstEntry)
sd=31.492*sqrt(83)
summary(lm(Familiar ~ Group_Size+scaleElo+FirstEntry+Group_Size:FirstEntry, data = modelData))
power.t.test(delta=57.055,sd=sd,n=83)
# given this very low power makes it almost impossible to find a significant effect even though there is 
# a true significant effect ==> so your tendency is clearly not by chance but true effect
# hence I would go with the final model as follows:
summary(model.test.3)
# and report the tendency with the effect size of the interaction
# report also the CIs
# no information to be reported on the main effect level.

# Model diagnostics of final model
resid.model<- simulateResiduals(model.test.3, 1000)
plot(resid.model) # good
plotResiduals(resid.model, form = modelData$FirstEntry)#good
plotResiduals(resid.model, form = modelData$Group_Size)#good

modelData$Group_Size = factor(modelData$Group_Size, levels=c("large", "small"), labels = c("large group", "small group"))

plotPoints = as.data.frame(emmeans(model.Familiar.red, pairwise ~ Group_Size | FirstEntry))
plotPoints = as.data.table(plotPoints)[1:4, ]

plotPoints$Group_Size = factor(plotPoints$Group_Size, levels=c("large", "small"), labels = c("large group", "small group"))


#plot of data
pp = ggplot(modelData, aes(x = FirstEntry, y = Familiar))+
  #geom_boxplot(aes(fill = Group_Size), alpha = 0.5)+
  geom_jitter(aes(colour = Pen),width = 0.1,
             size = 2)+
  geom_pointrange(data=plotPoints, mapping=aes(y=emmean, ymin=lower.CL, ymax=upper.CL), 
                  position = position_dodge(width = 0.75), size=1, shape=22)+
    labs(x = "Animal first visited", y = "Time spent with flockmate")+
  scale_x_discrete(labels=c("Flockmate", "Non-\nflockmate"))+
  theme_bw(base_size = 14)+
  scale_color_manual(values=c(largeCol, smallCol))+
  facet_grid(.~Group_Size)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggsave(path = "plots", "FlockmateTime2.tiff", pp, "tiff", width = 16, height = 11, units = "cm", dpi = 300)

#for comparison purposes: show time with familiar vs time with unfamiliar
longData = melt(modelData[, .SD, .SDcols = c("Pen", "ID", "scaleElo", "Group_Size", "Familiar", "Stranger", "FirstEntry")], 
                id.vars = c("Pen", "scaleElo", "Group_Size", "ID", "FirstEntry"),
                variable.name = "Area", value.name = "Time")
longData[, uID := paste0(Pen, ID)]

ggplot(longData, aes(x = Area, y = Time, colour = FirstEntry))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5)+
  facet_grid(.~Group_Size)+
  geom_point(position = position_jitterdodge(),
              size = 2)+
  labs(x = "Area", y = "Time spent in area")+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#descriptives
longData[, mean(Time), by = .(Group_Size, Area, FirstEntry)]

