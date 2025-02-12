#clear workspace
rm(list = ls())


library(data.table) # data.frame better
library(ggplot2) # plotting
library(lme4) # mixed models
library(DHARMa) # model diagnostics
library(parameters) # model parameters output
library(emmeans) # posthoc testing
library(effects) # model parameters output
library(RColorBrewer) # color for plotting
library(MASS) #for negative binomial glm
library(car) # testing colinearity in models
#library(MuMIn)

options(na.action = "na.fail")
set.seed(1)

#colours for plots
largeCol = brewer.pal(n = 8, name = "Blues")[c(4,6,8)]
smallCol = brewer.pal(n = 8, name = "OrRd")[c(4,6,8)]


#### Loading and preparing Data ################################################
arenaTest1 = fread("ArenaTest1.csv")
arenaTest2 = fread("ArenaTest2.csv")
Individ = fread("IndividualsOutput.csv")
Individ[, Unique:= paste0(Pen, ID)]
Thermal = fread("ThermalResults.csv")

#merging data
tempData = merge(arenaTest1[, c("Pen", "ID", "Side", "Box", "FirstChoice", "Feeding", "Weight"), with = F], 
                 Individ[, c("Pen", "ID", "Unique", "Group_Size", "scaleElos"), with = F], 
                 by = c("Pen", "ID"))

fullData = merge(arenaTest2[,c("Pen", "ID", "Side_Flockmate", "Flockmate","NonFlockmate", 
                               "FirstChoice", "FirstEntry", "Weight26"), with = F ], 
                 tempData, by = c("Pen", "ID"))
fullData = merge(Thermal[, c("Pen", "ID", "CombSize1", "CombSize2"), with = F], fullData, by = c("Pen", "ID"))

#fullData[, "GroupBin" := ifelse(Group_Size == "Large", 1, 0)]

#factorise
fullData[, Group_Size := as.factor(Group_Size)]
fullData[, Side := as.factor(Side)]
fullData[, Side_Flockmate := as.factor(Side_Flockmate)]
fullData[,Pen :=  as.factor(Pen)]
fullData[, FirstEntry := as.factor(FirstEntry)]


#test duration of 10 min added as offset
fullData[, Offset := 600]

#change label names for plotting with facets
fullData$plotGroup_Size <- factor(fullData$Group_Size, levels = c("large", "small"),
                                  labels = c("large groups", "small groups"))


#Careful: AKO has no Elo because of wrong labels during observations
rm(Thermal, arenaTest1, arenaTest2, Individ, tempData)

#### Focal individuals ############################################################

fullData[, list(Loser = sum(scaleElos <0), Winner = sum(scaleElos >0)), by = c("Group_Size")]
# -> more subordinate tendency animals than dominant... especially in large groups

#distribution of Elo ranking by Pen/Group_Size

#supplementary figure
#Selected individuals in each pen
ggplot(data = fullData, aes(x = Pen, y = scaleElos, fill = Pen, colour = Pen))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = Individ, aes(x = Pen, y = scaleElos, shape = "non focal"), colour = "black", alpha = 0.5, size = 1)+
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
ggplot(data = fullData, aes(x = Group_Size, y = scaleElos))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha = 0.7, size = 3)+
  theme_classic(base_size = 18)



#### Badges of status #############################################################


##### model statistics ##########################################################

#evaluated on the weight and comb size measured at 26 WoA 

hist(fullData[,Weight26]) # normal
hist(fullData[,CombSize2]) # near-normal

#are weight and comb size correlated?
cor.test(fullData[,CombSize2], fullData[,Weight26]) # not very high r = 0.28
ggplot(data= fullData, aes(y = Weight26, x = CombSize2)) + 
  geom_point()+ 
  geom_smooth(method = lm, formula = y~x) + 
  facet_grid(.~Group_Size)


#Effect of badges of status on Elo
badges.model = lmer(scaleElos ~ scale(Weight26)*scale(CombSize2)+ scale(CombSize2)*Group_Size + scale(Weight26)*Group_Size + (1|Pen), data = fullData)
#singularity warning
badges.model = lm(scaleElos ~ scale(Weight26)*scale(CombSize2)+ scale(CombSize2)*Group_Size + scale(Weight26)*Group_Size, data = fullData)
resid.Badges<- simulateResiduals(badges.model, 1000)
plot(resid.Badges)
plotResiduals(resid.Badges, form = fullData$Group_Size) #heteroscedacity? still ok
plotResiduals(resid.Badges, form = scale(fullData$CombSize2)) # okay
plotResiduals(resid.Badges, form = scale(fullData$Weight26)) # good

#test if comb size and weight too co-linear
vif(badges.model) # not too high co-linear (values <5)


#reduce interactions which don't improve model fit
drop1(badges.model) #best without any interactions

badges.model.red = lm(scaleElos ~ scale(Weight26)+scale(CombSize2)+Group_Size, data = fullData)
resid.Badges<- simulateResiduals(badges.model.red, 1000)
plot(resid.Badges)
plotResiduals(resid.Badges, form = fullData$Group_Size) #heteroscedacity? still ok
plotResiduals(resid.Badges, form = scale(fullData$CombSize2)) 
plotResiduals(resid.Badges, form = scale(fullData$Weight26)) 

#null model
badges.model_null = lm(scaleElos ~ 1, data = fullData)

AIC(badges.model.red, badges.model_null) # null model better
anova(badges.model_null, badges.model.red) # no significant difference

#variance explained
r.squaredGLMM(badges.model.red, badges.model_null)

tab_model(badges.model.red)


##### plots ####################################################################

#not used in paper

#Comb size and Elo rating (showing both measurement time points)
ggplot(data= fullData, aes(y = scaleElos)) + 
  geom_point(aes( x = CombSize2))+ 
  geom_smooth(aes( x = CombSize2), method = lm, formula = y~x, colour = "black") + 
  geom_point(aes( x = CombSize1))+ 
  geom_smooth(aes( x = CombSize1), method = lm, formula = y~x) + 
  facet_grid(.~Group_Size)+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Comb Size (cm²)", y = "scaled Elo rating")

#weight and Elo rating (showing both measurement time points)
ggplot(data= fullData, aes(y = scaleElos)) + 
  geom_point(aes( x = Weight26))+ 
  geom_smooth(aes( x = Weight26), method = lm, formula = y~x, colour = "black") + 
  geom_point(aes( x = Weight))+ 
  geom_smooth(aes( x = Weight), method = lm, formula = y~x) +
  facet_grid(.~plotGroup_Size)+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Weight (g)", y = "scaled Elo rating")


#### Fear test #################################################################

###### Latency #################################################################

fullData$Box # several times 600 -> censoring offset needed
hist(fullData$Box) # difficult distribution, try normal with offset

Lat.model = lmer( Box ~ scaleElos*Group_Size + offset(Offset) + (1|Pen), data = fullData)
resid.Lat = simulateResiduals(Lat.model, 1000)
plot(resid.Lat)#no good fit

#try poisson
Lat.model = glmer( Box ~ scaleElos*Group_Size + offset(log(Offset))+ (1|Pen), family = poisson, data = fullData)
resid.Lat = simulateResiduals(Lat.model, 1000)
plot(resid.Lat) #even worse

#try negative binomial
Lat.model = glmer.nb( Box ~ scaleElos*Group_Size + offset(log(Offset))+ (1|Pen), data = fullData)
#singularity fit

#without Pen as it does not seem to account for any variation
Lat.model = glm.nb( Box ~ scaleElos*Group_Size + offset(log(Offset)), data = fullData)
resid.Lat = simulateResiduals(Lat.model, 1000)
plot(resid.Lat)
plotResiduals(resid.Lat, form = fullData$scaleElos) #decent
plotResiduals(resid.Lat, form = fullData$Group_Size) # decent

#reduce model complexity
drop1(Lat.model) # model is better without interaction

Lat.model.red = glm.nb( Box ~ scaleElos+Group_Size + offset(log(Offset)), data = fullData)
resid.Lat = simulateResiduals(Lat.model.red, 1000)
plot(resid.Lat)
plotResiduals(resid.Lat, form = fullData$scaleElos)
plotResiduals(resid.Lat, form = fullData$Group_Size)

#intercept only model
Lat.model.null = glm.nb( Box ~ 1+ offset(log(Offset)), data = fullData)

AIC(Lat.model.red, Lat.model.null) #intercept model better 
anova(Lat.model.red,Lat.model.null, test = "Chisq") # no significant difference

#variance explained
r.squaredGLMM(Lat.model.red, Lat.model.null)

tab_model(Lat.model.red)

##

#plots

#not used in paper
#Latency to leave box by group size and rank
#as ratio -> dived by full time given (offset)
ggplot(fullData, aes(x = scaleElos, y = Box/Offset))+
  geom_point()+
  geom_smooth(method="glm",method.args=list(family="binomial"), formula = y~x)+
  facet_grid(.~Group_Size)+
  theme_bw(base_size = 16)


###### First Choice Feed vs. Social ##############################################

#deal with null entries
fullData[FirstChoice.y == "", FirstChoice.y := "No Choice"]

#overview of how many of each choice we have
comp.Choice = fullData[, .N, by = c("FirstChoice.y", "Side", "Group_Size")]

ggplot(data = comp.Choice, aes(x = FirstChoice.y, y = N, fill = Side))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)+
  facet_grid(. ~Group_Size)


# Side bias:
# select only animals that made a choice for the analysis 
choiceData1 = fullData[FirstChoice.y != "No Choice",]

#determining whether the animal went left or right depending on starting position
choiceData1[Side == "A", Side_Social := "Right"]
choiceData1[Side == "S", Side_Social := "Left"]
choiceData1[Side == "A", ChoiceDirection := ifelse(FirstChoice.y == "Social", "Right", "Left")]
choiceData1[Side == "S", ChoiceDirection := ifelse(FirstChoice.y == "Social", "Left", "Right")]

#factorise
choiceData1[,ChoiceDirection := as.factor(ChoiceDirection)]
choiceData1[,Side_Social := as.factor(Side_Social)]

#create Boolean with right set as 0 and left as 1
choiceData1[,ChoiceDirectionBool := ifelse(ChoiceDirection == "Right", 0, 1)]

#overview of data
comp.Choice = choiceData1[, .N, by = c("ChoiceDirection", "Side_Social", "Group_Size")]

#appear to mostly go towards the conspecific
ggplot(data = comp.Choice, aes(x = ChoiceDirection, y = N, fill = Side_Social))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)+
  facet_grid(. ~Group_Size)

#no apparent side bias
ggplot(comp.Choice, aes(x = Group_Size, y = N,  fill = ChoiceDirection))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)

###

#model statistics

#right set as 0 left as 1
choice.model1 = glmer(ChoiceDirectionBool ~ scaleElos*Group_Size+ Group_Size*Side_Social + scaleElos*Side_Social + (1|Pen), data = choiceData1, family = "binomial")
#singularity
choice.model1 = glm(ChoiceDirectionBool ~ scaleElos*Group_Size+ Group_Size*Side_Social + scaleElos*Side_Social, data = choiceData1, family = "binomial")
resid.choice<- simulateResiduals(choice.model1, 1000)
plot(resid.choice)
plotResiduals(resid.choice, form = choiceData1$Side_Social) # decent
plotResiduals(resid.choice, form = choiceData1$scaleElos) #decent
plotResiduals(resid.choice, form = choiceData1$Group_Size) # good

#reduce model complexity
drop1(choice.model1) #best without interactions

choice.model1.red = glm(ChoiceDirectionBool ~ Group_Size + Side_Social + scaleElos, data = choiceData1, family = "binomial")
resid.choice<- simulateResiduals(choice.model1.red, 1000)
plot(resid.choice)
plotResiduals(resid.choice, form = choiceData1$Side_Social)#good
plotResiduals(resid.choice, form = choiceData1$scaleElos) #decent
plotResiduals(resid.choice, form = choiceData1$Group_Size)#good

#null model
choice.model1.null = glm(ChoiceDirectionBool ~ 1, data = choiceData1, family = "binomial")

AIC(choice.model1.null, choice.model1.red) #much better than null model

summary(choice.model1.red)
parameters(choice.model1.red, exponentiate = TRUE)
tab_model(choice.model1.red)
plot(predictorEffects(choice.model1.red), lines=list(multiline=TRUE))

#variance explained
r.squaredGLMM(choice.model1.red, choice.model1.null)

#probabilities
emmeans(choice.model1.red, ~ Side_Social, type = "response")



###### Probability to feed #####################################################

#Binomial model
Feeding.model = glmer(Feeding ~ scaleElos*Group_Size +(1|Pen), data = fullData, family = binomial)
#singularity
Feeding.model = glm(Feeding ~ scaleElos*Group_Size, data = fullData, family = binomial)
resid.Feeding = simulateResiduals(Feeding.model, 1000)
plot(resid.Feeding)
plotResiduals(resid.Feeding, form = fullData$Group_Size) #good
plotResiduals(resid.Feeding, form = fullData$scaleElos) #good

#reduce model complexity
drop1(Feeding.model) # better without interaction

Feeding.model.red = glm(Feeding ~ scaleElos + Group_Size, data = fullData, family = binomial)
resid.Feeding<- simulateResiduals(Feeding.model.red, 1000)
plot(resid.Feeding)
plotResiduals(resid.Feeding, form = fullData$Group_Size) #good
plotResiduals(resid.Feeding, form = fullData$scaleElos) # good

#null model
Feeding.null = glm(Feeding ~ 1, data = fullData, family = binomial)

AIC(Feeding.null, Feeding.model) # only slightly better than null model
anova(Feeding.null, Feeding.model, test = "Chisq") # trend to be better than null model

summary(Feeding.model.red) 
parameters(Feeding.model.red, exponentiate = T)
tab_model(Feeding.model.red)

#variance explained
r.squaredGLMM(Feeding.model.red, Feeding.null)

### plot

plotData = summary(emmeans(Feeding.model.red, ~scaleElos, at = list(scaleElos = seq(from = min(fullData$scaleElos),
                                                               to = max(fullData$scaleElos), by = 0.1)), type = "response"))

ggplot(plotData, aes(x = scaleElos, y = prob))+
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


#### Recognition test ##########################################################

###### First Choice Flockmate vs. NonFlockmate ##################################

#Overview of side bias
comp.Choice = fullData[, .N, by = c("FirstChoice.x", "Side_Flockmate", "Group_Size")]
ggplot(data = comp.Choice, aes(x = FirstChoice.x, y = N, fill = Side_Flockmate))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)+
  facet_grid(. ~Group_Size)

#exclude hens which made no choice
choiceData2 = fullData[!is.na(FirstChoice.x),]

#factorise
choiceData2[, FirstChoice.x := as.factor(FirstChoice.x)]

#extract which hens went left and which right dependning on starting position 
choiceData2[Side_Flockmate == "Left", ChoiceDirection := ifelse(FirstChoice.x == "NonFlockmate", "Right", "Left")]
choiceData2[Side_Flockmate == "Right", ChoiceDirection := ifelse(FirstChoice.x == "NonFlockmate", "Left", "Right")]

#right is set as reference 0
choiceData2[,ChoiceDirectionBool := ifelse(ChoiceDirection == "Right", 0, 1)]

#overview of data
comp.Choice = choiceData2[, .N, by = c("ChoiceDirection", "Side_Flockmate", "Group_Size")]
comp.Choice$Group_Size <- factor(comp.Choice$Group_Size, levels = c("large", "small"),
                                 labels = c("large groups", "small groups")
)

#strange side bias in small group?
ggplot(data = comp.Choice, aes(x = ChoiceDirection, y = N, fill = Side_Flockmate))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)+
  labs(x = "First choice of direction", y = "Number of hens", fill = "Position of Flockmate")+
  facet_grid(. ~Group_Size)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_discrete(labels = c("Flockmate to the left", "Flockmate to the right"))

#side bias in small group!
ggplot(comp.Choice, aes(x = Group_Size, y = N,  fill = ChoiceDirection))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_bw(base_size = 18)


###

#model statistics

#binomial model
choice.model2 = glmer(ChoiceDirectionBool ~ scaleElos*Group_Size + scaleElos*Side_Flockmate + Group_Size*Side_Flockmate + (1|Pen), data = choiceData2, family = "binomial")
#Singularity
choice.model2 = glm(ChoiceDirectionBool ~ scaleElos*Group_Size + scaleElos*Side_Flockmate + Group_Size*Side_Flockmate, data = choiceData2, family = "binomial")
resid.choice<- simulateResiduals(choice.model2, 1000)
plot(resid.choice)
plotResiduals(resid.choice, form = choiceData2$Side_Flockmate)#good
plotResiduals(resid.choice, form = choiceData2$Group_Size)#good
plotResiduals(resid.choice, form = choiceData2$scaleElos)#good

#reduce model complexity
drop1(choice.model2) #no clear indication that any interactions make the model better

choice.model2.red = glm(ChoiceDirectionBool ~ scaleElos+ Side_Flockmate + Group_Size, data = choiceData2, family = "binomial")
resid.choice<- simulateResiduals(choice.model2.red, 1000)
plot(resid.choice) 
plotResiduals(resid.choice, form = choiceData2$Side_Flockmate)#good
plotResiduals(resid.choice, form = choiceData2$Group_Size)#good
plotResiduals(resid.choice, form = choiceData2$scaleElos)#good

#null model
choice.model2.null = glm(ChoiceDirectionBool ~ 1, data = choiceData2, family = "binomial")

AIC(choice.model2.null, choice.model2.red) # no clear difference
anova(choice.model2.null, choice.model2.red, test = "Chisq") #no difference

summary(choice.model2.red)
parameters(choice.model2.red, exponentiate = TRUE)
tab_model(choice.model2.red)

#variance explained
r.squaredGLMM(choice.model2.red, choice.model2.null)

#probability
emmeans(choice.model2.red, ~ Group_Size, type = "response")


###### Flockmate time ##########################################################


#for Arena2 test take out animal that did not move at all because it has no first entry
flockmateData = fullData[!is.na(FirstEntry),]

flockmateData$Flockmate #no 600 so no censoring needed
hist(flockmateData$Flockmate) #looks near normal

# here comes the gaussian based on the distribution seen before:
flockmate.model = lmer(Flockmate ~ Group_Size:scaleElos + Group_Size:FirstEntry + 
                                  scaleElos:FirstEntry + Group_Size + scaleElos + FirstEntry + (1|Pen),
                  data = flockmateData)
# singular fit issue, most likely due to Pen:
# therefore, lets try out an lm model:
flockmate.model = lm(Flockmate ~ Group_Size:scaleElos + Group_Size:FirstEntry + 
                  scaleElos:FirstEntry + Group_Size + scaleElos + FirstEntry,
                data = flockmateData)
resid.test<- simulateResiduals(flockmate.model, 1000)
plot(resid.test)# very even distribution of residuals
plotResiduals(resid.test, form = flockmateData$Group_Size) #good
plotResiduals(resid.test, form = flockmateData$scaleElos) # still okay
plotResiduals(resid.test, form = flockmateData$FirstEntry) # good

#reduce model complexity
drop1(flockmate.model) # include GroupSize:FirstEntry

flockmate.model.red = lm(Flockmate ~ Group_Size:FirstEntry +Group_Size+scaleElos+FirstEntry, data = flockmateData)
resid.model<- simulateResiduals(flockmate.model.red, 1000)
plot(resid.model) # good
plotResiduals(resid.model, form = flockmateData$FirstEntry)#good
plotResiduals(resid.model, form = flockmateData$Group_Size)#good
plotResiduals(resid.model, form = flockmateData$scaleElos)#still okay

#null model
flockmate.model.null = lm(Flockmate ~ 1, data = flockmateData)

AIC(flockmate.model.null, flockmate.model.red) #slightly better than null

summary(flockmate.model.red)
parameters(flockmate.model.red)
tab_model(flockmate.model.red)
plot(allEffects(flockmate.model.red))


#varaince explained
r.squaredGLMM(flockmate.model.red, flockmate.model.null)

#post-hoc comparisons
emmeans(flockmate.model.red, pairwise ~   FirstEntry*Group_Size) 


###

#plots

#change level order for plotting
flockmateData$Group_Size = factor(flockmateData$Group_Size, levels=c("large", "small"), labels = c("large group", "small group"))

plotPoints = as.data.frame(emmeans(flockmate.model.red, pairwise ~ Group_Size | FirstEntry))
plotPoints = as.data.table(plotPoints)[1:4, ]

plotPoints$Group_Size = factor(plotPoints$Group_Size, levels=c("large", "small"), labels = c("large group", "small group"))


#plot used in paper
#time spent with flockmate by group size and animal first visited
pp = ggplot(flockmateData, aes(x = FirstEntry, y = Flockmate))+
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

#save plot
ggsave(path = "plots", "FlockmateTime2.tiff", pp, "tiff", width = 16, height = 11, units = "cm", dpi = 300)
