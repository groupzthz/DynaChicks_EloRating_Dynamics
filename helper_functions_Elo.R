#Some of the following functions have been adapted from functions written by 
#Damien Farine & Alfredo Sanchez-Tojar for the AniDom package


#function: dataset_diagnostics
# calculates diagnostics for the input dataset
# Input: dataset = social observation sequence (data.frame or data.table)
#        Individuals = all existing individuals in the group (data.frame or data.table)
# Outputs list of: 
#        'seen_notSeen' = how many individuals were observed, how many were never observed  
#        'tot_interact' = total number of interactions
#        'ratio_interact' =  ratio of number of interactions divided by the amount of individuals observed (d/N)
#        'prop_unknown' = the proportion of dyads in data set for which no interactions have been observed
#        'prop_known' = the proportion of dyads in data set for which interactions have been observed
#        'triang_trans' = Triangle Transitivity  as described in Shizuka & McDonald, 2012
#        'focus_pos' = Focus and Position as described in Hobson et al., 2021
#        'blur' = Downward null model as described in Hobson et al., 2021
#        'strategy' = Pattern of aggression as described in Hobson et al., 2021

dataset_diagnostics <- function(dataset, Individuals){
  
  winners = dataset$Winner
  losers = dataset$Loser
  date = dataset$Date
  cat("New Dataset started\n")
  #create matrix from interactions
  data_matrix = as.data.frame(creatematrix(winner = winners, loser = losers))
  
  #How many animals were observed, how many not
  seen_notSeen = table((Individuals %in% names(data_matrix)))
  
  # How many interactions? 
  tot_interact = sum(data_matrix)
  
  #Ratio of interactions to individuals d/N
  ratio_interact = round(sum(data_matrix)/length(names(data_matrix)), 2)
  
  #how large is the proportion of dyads in data set for which no interactions have
  #been observed
  prop_unknown = prunk(as.matrix(data_matrix))
  
  cat("Calculating expected under Poisson distribution\n")
  #does the observed proportion of known dyads lie within the range of what
  # you would expect under a Poisson process?
  # commented out because of compete package problems under this version of R
  prop_known = proportion_of_known(data_matrix)
  
  cat("Calculating Triangle Transitivity\n")
  # Triangle Transitivity (Shizuka & McDonald 2012)
  ttri = round(transitivity(as.matrix(data_matrix), runs = 1000),7)
  
  cat("Calculating Focus, Position and Strategy\n")
  # based on patterns of aggression (Hobson et al., 2021) "domstruc" package
  
  # Compute bootstrap estimates of confidence intervals for focus and position 
  # given an aggression matrix  
  focus_pos = dom_make_data(as.matrix(data_matrix))
  
  #Compute bootstrap estimates of focus and position for the downward null model 
  #with confidence intervals 
  blur = dom_make_blur_data(as.matrix(data_matrix))
  
  #Categorize group-level strategy -> pattern of aggression
  strategy = dom_categorize_strategy(focus_pos, blur, use_data_ci = TRUE)
  
  
  output = list(seen_notSeen, tot_interact, ratio_interact, prop_unknown, prop_known, 
                ttri,focus_pos, blur, strategy)
  
  names(output) = c('seen_notSeen', 'tot_interact', 'ratio_interact', 'prop_unknown', 'prop_known', 
                    'triang_trans','focus_pos', 'blur', 'strategy')
  
  return(output)
}

#function: proportion_of_known
# calculates whether the observed proportion of known dyads lie within the range of what you would expect under a Poisson process
# Input: input_matrix = social observation matrix
# Outputs list of: 
#        'exp_poiss' = expected proportion of known dyads under a Poisson process
#        'exp_poiss_CI' = confidence inetravl of expected proportion
#        'obs_known' =  the proportion of dyads in data set for which interactions have been observed

proportion_of_known <- function(input_matrix) {
  
  data = input_matrix
  avalues = c(30)
  bvalues = c(5)
  N.inds.values <- c(table(rowSums(data)+colSums(data)!=0))
  N.obs.values <- c(sum(data)/(table(rowSums(data)+colSums(data)!=0)))
  poiss = c(TRUE)
  dombias = c(FALSE)
  db.sim <- data.frame(Ninds=integer(),
                       Nits=integer(),
                       poiss=logical(),
                       dombias=logical(),
                       unknowndyads=numeric(),
                       stringsAsFactors=FALSE)
  
  for (simnum in 1:500){
       # generating interactions using aniDom
      # supress warning that one variable has a column name that is deleted
      output <- suppressWarnings(generate_interactions(N.inds.values,
                                      N.inds.values*N.obs.values,
                                      a=avalues,
                                      b=bvalues,
                                      id.biased=poiss,
                                      rank.biased=dombias))
      # generating sociomatrix and estimating number of
      matrix<-creatematrix(winners = output$interactions$Winner, losers = output$interactions$Loser)
      unknowndyads<-prunk(matrix)[1]
      
      # adding values to db
      db.sim<-rbind(db.sim,c(N.inds.values,N.obs.values,
                             poiss,dombias,
                             unknowndyads))
     }
  
  names(db.sim) <- c("Ninds","N.obs.values",
                     "poiss","dombias",
                     "unknowndyads")
  
  db.sim$knowndyads <- 1-db.sim$unknowndyad
  
  exp_poiss = round(mean(db.sim$knowndyads),2)
  
  exp_poiss_CI = round(quantile(db.sim$knowndyads, probs = c(0.025, 0.975)), 2)
  
  obs_known =  1-prunk(as.matrix(data))[1]
  
  output = list(exp_poiss, exp_poiss_CI, obs_known)
  
  names(output) = c('exp_poiss', 'exp_poiss_CI', 'obs_known')
  return(output)
  
}

#function: elo_analysis
# calculates Elo ratings for the social interaction sequence and estimates uncertainty
# Input: dataset = social observation sequence (data.frame or data.table)
#        Individuals = all existing individuals in the group (data.frame or data.table)
# Outputs list of: 
#       'summary' = provides a summary of the input data
#       'rating' = Elo rating in original order using AniDom package
#       'hierarchy_shape' = plot of the hierarchy shape
#       'uncer_repeat' = repeatability of Elo ratings
#       'uncer_split'= median split performance
#       'uncer_split_rand' = randomised split performance
#       'rating_rand' = Elo rating based on randomised Elo ratings
#       'comp_rating' = Spearman rank correlation based on original Elo ratings and randomised Elo ratings
#       'comp_elos' = ICC of Elo ratings from original order data compared to Elo ratings from randomised order
#       'hierarchy_shape_rand' =  plot of the hierarchy shape based on randomsied ELo ratings
#       'Individuals' = optional data.table with details for the Individuals 
#       'rankMatch' = optional data.table with details for the social interactions
#       'randElo' = randomised Elo rating out in full (elo ratings for each randomisation) 


elo_analysis = function(dataset, Individuals){
  
  winners = dataset$Winner
  losers = dataset$Loser
  date = dataset$Date
  #presence = data.table(Date = seq(from = unique(date[1]), to = unique(date[length(date)]), by = "day"))
  #presence = cbind(presence, as.data.table(matrix(1,length(presence$Date), length(Individuals))))
  #presence[!(Date %in% date),2:(length(Individuals)+1)] = 0
  #names(presence)[2:(length(Individuals)+1)] = Individuals
  cat("New Dataset started\n")
  
  # 1. Original Elo extraction
  # not exactly the same values as in Elo_rating package -> sigmoid.param?
  # does it matter for this paper? I think not?
  
  #origElo = elo.seq(winner = winners, loser = losers, 
  #                 Date = date,runcheck = TRUE, startvalue = 0, k = 200, presence = presence)#k = oK$best$k)
  origElo = elo.seq(winner = winners, loser = losers, 
                    Date = date,runcheck = TRUE, startvalue = 0, k = 200)
  #print summary of input data
  summary(origElo)
  
  #using AniDom elo package
  newElo = elo_scores(winners=winners,losers=losers,
                    identities = Individuals, randomise=FALSE)
  
  rating = newElo[order(newElo[,1]),]
  
  #takes a set of winners and losers from observed interactions and plots the probability
  #of the dominant individual in an interaction winning given the difference in rank to the subordinate
  #in the same interaction.
  group.size = length(rating)
  
  #hierarchy_shape = my_plot_hierarchy_shape(identity= names(rating),
  #                                           rank=length(rating):1,
  #                                           winners=winners,
  #                                           losers=losers,
  #                                          group.size = group.size)
  
  # estimates repeatability of Elo ratings -> 
  # "Simulations suggest that a repeatability score above 0.8 suggests
  # a reasonably robust hierarchy, given a large input dataset 
  #(can be unreliable for small datasets, i.e. < 10 observations per individual), 
  # or for extremely flat hierarchies"
  # from rptr Gaussian
  cat('Calculating uncertainity by repeatability\n')
  uncer_repeat = estimate_uncertainty_by_repeatability(winners, losers, 
                                               identities = Individuals,
                                               n.rands = 1000)
  #the Spearman rank correlation coefficient of the first half and second half
  #of the data.
  #non-randomised
  uncer_split = estimate_uncertainty_by_median_splitting(winners, losers, 
                                                  identities = Individuals,
                                                   randomise = FALSE)
  #randomised
  #the mean and 95% range of the Spearman rank correlation
  #coefficients from two halves of the data with the ordering randomised each time.
  #Our simulations suggest that correlations above 0.5 suggests a
  #robust dominance hierarchy (or low uncertainty).
  cat('Calculating uncertainity by splitting\n')
  uncer_split_rand = estimate_uncertainty_by_splitting(winners, losers, 
                                                  identities = Individuals,
                                                  randomise = TRUE,
                                                  n.rands = 1000)
  
  
  # 2. Randomised Elo
  
  cat('Calculating randomised Elo\n')
  randElo = elo_scores(winners=winners,losers=losers,
                       identities = Individuals, 
                        randomise=TRUE, n.rands = 1000)
  
  rating_rand = rowMeans(randElo)
  rating_rand = rating_rand[order(rating_rand)]
  
  #deals with NA, which can interfere with rankings
  if (any(is.na(rating_rand))){
    ranking_rand = c((length(Individuals)-sum(is.na(rating_rand))):1, rep(NA, sum(is.na(rating_rand))))
    names(ranking_rand) = names(rating_rand)
  }else {
    ranking_rand = length(Individuals):1
    names(ranking_rand) = names(rating_rand)
  }
 
  #hierarchy_shape_rand = my_plot_hierarchy_shape(identity= names(rating_rand),
  #                                           rank=length(rating_rand):1,
  #                                           winners=winners,
  #                                           losers=losers,
  #                                           group.size = group.size)
  
  # plotting ranks after randomised Elo with CI 
  #plot_ranks(randElo,plot.CIs=TRUE, ordered.by.rank = TRUE)
  
  # 3 . Comparison randomised and original
  
  comp_rating = cor.test(rating[sort(names(rating))], rating_rand[sort(names(rating_rand))],  method = 'spearman',na.rm = TRUE)
  comp_elos = icc(data.table(rating[sort(names(rating))], rating_rand[sort(names(rating_rand))]), model = "twoway", type = "agreement")
  
  
  diff_rating = rating[sort(names(rating))]-rating_rand[sort(names(rating_rand))]
  
  
  
  # 4. Individual based diagnostics
  
  # Interactions by individuals

    #Overview by individual
    Individuals = data.table(ID = Individuals)
    Individuals[, Wins := sum(winners == ID), by = ID]
    Individuals[, Losses := sum(losers == ID), by = ID]
    Individuals$elos = rating_rand[sort(names(rating_rand))]
    Individuals$rank = ranking_rand[sort(names(ranking_rand))]
    Individuals[, sum := Wins+Losses]
    Individuals[,scaleElos := scale(elos)] 
    Individuals[, physAggr := sum((dataset$Code == "Peck" |dataset$Code == "Fight") &
                                    winners == ID), by = ID]
    Individuals[, nonphysAggr := sum((dataset$Code == "Threat" |dataset$Code == "Avoidance" )&
                                       winners == ID), by = ID]
    Individuals[, HQ := sum(dataset$Condition == "HQ" &
                            winners == ID), by = ID]
    Individuals[, Feed := sum(dataset$Condition == "Feeder" &
                                winners == ID), by = ID]
    Individuals[, Normal := sum(dataset$Condition == "Normal" & 
                                  winners == ID), by = ID]
    Individuals[, physAggrRec := sum((dataset$Code == "Peck" |dataset$Code == "Fight") &
                                    losers == ID), by = ID]
    Individuals[, nonphysAggrRec := sum((dataset$Code == "Threat" |dataset$Code == "Avoidance" )&
                                          losers == ID), by = ID]
    Individuals[, HQRec := sum(dataset$Condition == "HQ" &
                                 losers == ID), by = ID]
    Individuals[, FeedRec := sum(dataset$Condition == "Feeder" &
                                   losers == ID), by = ID]
    Individuals[, NormalRec := sum(dataset$Condition == "Normal" & 
                                     losers == ID), by = ID]
    
    #each interaction plus end ranks and elos
    rankMatch = data.table(Winner = winners, Loser = losers, Code = dataset$Code, Situation = dataset$Condition)
    rankMatch[, WinnerRank := Individuals$rank[match(winners, Individuals$ID)]]
    rankMatch[, LoserRank := Individuals$rank[match(losers, Individuals$ID)]]
    rankMatch[, WinnerElo := Individuals$scaleElos[match(winners, Individuals$ID)]]
    rankMatch[, LoserElo := Individuals$scaleElos[match(losers, Individuals$ID)]]

  

  output = list(origElo,rating, uncer_repeat,uncer_split,uncer_split_rand,
                rating_rand,comp_rating, comp_elos, Individuals, rankMatch, randElo)
  names(output) = c('Elo_object','rating','uncer_repeat','uncer_split','uncer_split_rand',
                'rating_rand','comp_rating','comp_elos', 'Individuals', 'rankMatch', 'randElo')

  return(output)
}

#function to estimate uncertainty based on a median split
estimate_uncertainty_by_median_splitting <-
  function(winners, losers, identities=NULL, sigmoid.param=1/100, K=200, init.score=0, randomise=FALSE) {
    
    if (is.null(identities)) {
      identities <- unique(c(winners,losers))  
    }
    
    n.inds <- length(identities)
    
    n.observations <- length(winners)
    obs1 <- seq(1,n.observations/2,1)
    obs2 <- seq(max(obs1)+1,n.observations,1)
    
    scores1 <- elo_scores(winners[obs1],losers[obs1], identities, sigmoid.param, K, init.score, randomise= FALSE)
    scores2 <- elo_scores(winners[obs2],losers[obs2], identities, sigmoid.param, K, init.score, randomise = FALSE)
    scores.cor <- cor.test(scores1,scores2,method="spearman")
    # in case we want CIs
    # package RVAideMemoire
    #spearman.ci(var1, var2, nrep = 1000, conf.level = 0.95)
   
    return(scores.cor)
    
  }

#OLD version to calculate steepness -> instead used EloSteepness package
#function to create hierarchy shape plots and to calculate the steepness of the hierarchy
#takes a set of winners and losers from observed interactions and plots the probability
#of the dominant individual in an interaction winning given the difference in rank to the subordinate
#in the same interaction.
# Outputs list of:
#     plot =  plot of the shape
#     data = data.table of data used for plot
#     rate = steepness k of the hierarchy 
#     confint_Rate = 95 % confidence interval for steepness

# my_plot_hierarchy_shape <-
#   function(identity, rank, winners, losers, group.size, standard = 10) {
#     
#     winners.rank <- rank[match(winners,identity)]
#     losers.rank <- rank[match(losers,identity)]
#     xx <- winners.rank-losers.rank
#     x <- 1:(max(abs(xx)))
#     y <- rep(NA,length(x))
#     totInteract <- rep(NA,length(x))
#     sumHighWin <- rep(NA,length(x))
#     CI.upper <- y
#     CI.lower <- y
#     for (i in 1:length(x)) {
#       totInteract[i] = sum(abs(xx)==x[i])
#       sumHighWin[i] = sum(xx==-x[i])
#       y[i] <- sumHighWin[i]/totInteract[i]
#       CI.upper[i] <- y[i] + sqrt(y[i]*(1-y[i])/sum(abs(xx)==x[i])) + 0.5/sum(abs(xx)==x[i])
#       CI.upper[i] <- min(CI.upper[i],1)
#       CI.lower[i] <- y[i] - sqrt(y[i]*(1-y[i])/sum(abs(xx)==x[i])) - 0.5/sum(abs(xx)==x[i])
#       CI.lower[i] <- max(CI.lower[i],0)
#     }
#     CI.upper <- CI.upper[!is.na(y)]
#     CI.lower <- CI.lower[!is.na(y)]
#     x <- x[!is.na(y)]
#     y <- y[!is.na(y)]
#     totInteract <- totInteract[!is.na(totInteract)]
#     
#     a = standard
#     
#     #old formula
#     #logisticModel <- nls(y~1/(1+exp(-1*(a/group.size)*r*x)), # -1 to keep growth rate positive, a/N as  
#     #start=list(r=-1), 
#     #control=list(maxiter=1000, minFactor=.00000000001))
#     
#     #to standardise all rank differences to a scale from 1 to 10
#     #include difference of 0 
#     x = c(0,x)
#     standard_X = (x-min(x))/(max(x)-min(x))*10
#     #baseline-probability that higher ranking individual wins -> b = 0, probability = 0.5
#     b = 0
#     base_prob = 1/(1+exp(b)) 
#     y = c(base_prob, y)
#     logisticModel <- nls(y~1/(1+exp(-1*r*standard_X + b)), # -1 to keep growth rate positive, a/N as  
#                          start=list(r=-1), 
#                          control=list(maxiter=1000, minFactor=.00000000001))
#     
#     
#     growth.rate = coef(logisticModel)
#     library(nlstools)
#     
#     if(sum(y<1)/length(y) > 0.2){
#       confint_Rate = confint2(logisticModel, method = "asymptotic")
#     }else { confint_Rate = "Not computable, more than 80% > 1"}
#     
#     
#     plotdata = data.table(RankDiff = x[-1], ProbWin = y[-1], Predict = predict(logisticModel)[-1])
#     
#     if(length(x) %% 2 == 0){
#       breaks.x = seq(2,length(x),2)
#     }else{breaks.x = seq(1,length(x),2)}
#     
#     library(scales)
#     plot = ggplot(plotdata, aes(x = RankDiff, y= ProbWin))+
#       geom_point(size = 3)+
#       #geom_smooth(aes(col = 'loess fit'), method = 'loess', formula = y~x, size = 1.5, se = FALSE)+
#       geom_errorbar(aes(x= RankDiff, ymin=CI.lower, ymax=CI.upper), width=.1)+
#       geom_line(aes(x=RankDiff, y=Predict, color = "logistic fit"), size =1.5)+
#       labs(colour = "Fit function")+
#       scale_y_continuous(limits=c(0.35,1),oob = rescale_none)+
#       scale_x_continuous(breaks = breaks.x)+
#       theme_bw(base_size = 18)+
#       theme(axis.title.x = element_blank(), axis.title.y = element_blank())
#     
#     invisible(list(plot = plot, 
#                    data = data.table(tot.Interact = totInteract,
#                                      Rank.diff=x[-1],Prob.dom.win=y[-1],
#                                      CI.upper=CI.upper,CI.lower=CI.lower),
#                    rate = growth.rate,
#                    confint_Rate = confint_Rate
#     ))
#     
#   }


# function to have a look at the elo_analysis outcome in a better formatting
# takes the result of elo_analysis and outputs as formatted text in the console

printCalculations <- function(rating){
  cat("Summary:\n")
  summary(rating$Elo_object)
  #cat("\n \n Steepness:\n", rating$hierarchy_shape$rate, "Confint:",rating$hierarchy_shape$confint_Rate)
  cat("\n \n Repeatability:\n", rating$uncer_repeat)
  cat("\n \n Median split:\n")
  print(rating$uncer_split)
  cat("\n \n Randomised Split:\n",rating$uncer_split_rand[1], "Confint:", rating$uncer_split_rand[2:3])
  cat("\n \n Comparison Original and Randomised:\n") 
  print(rating$comp_rating)
  print(rating$comp_elos)
  #cat("\n \n Randomised Elo: \n")
  #cat("\n Steepness:\n", rating$hierarchy_shape_rand$rate, "Confint:",rating$hierarchy_shape_rand$confint_Rate)
}


# Function to plot focus and position plots
# adjusted from https://rdrr.io/github/danm0nster/domstruc/src/R/plotting.R

library(alphahull)
library(scales)

#' Plot focus, position and show strategy regions based on downward null
#' heuristic.
#'
#' @param data A data frame containing one row with columns `focus` and `position`
#' @param blur_data A data frame that is created by `dom_make_blur_data()`
#'
#' @return
#' @export
#'
#' @examples
adj_dom_plot_strategy <- function(data, blur_data, show_data_ci = FALSE, show_legend = FALSE) {
  # TODO: Check input data is in the right format (option to not have ci in data if show_data_ci == FALSE)
  # Get convex hull of blur data with error bars
  ahuld.conv <- convex_hull(blur_data)
  
  polygons <- make_polygons(blur_data)
  bully.poly <- polygons[["bully"]]
  clcomps.poly <- polygons[["clcomps"]]
  undef.poly <- polygons[["undef"]]
  
  ########## plot showing strategy polygons
  # plot real data
  with(data, plot(focus, position,
                  ylim = c(-0.1, 1), xlim = c(-0.1, 1), las = 1,
                  col = "blue", pch = 5,
                  xlab = "focus", ylab = "position", cex = 2
                  
  ))
  
  # Add strategy polygons to bottom plotted layer
  # undefined (light grey)
  polygon(undef.poly, col = scales::alpha("grey", 0.4),
          border = scales::alpha("black", 0.5), lty = 2, lwd = 1.5)
  
  # bully (blue)
  polygon(bully.poly, col = scales::alpha("blue", 0.2),
          border = scales::alpha("black", 0.5), lty = 2, lwd = 1.5)
  
  # close competitors (red)
  polygon(clcomps.poly, col = scales::alpha("red", 0.2),
          border = scales::alpha("black", 0.5), lty = 2, lwd = 1.5)
  
  # downward heuristic (black, plotted over white)
  polygon(ahuld.conv, col = scales::alpha("white", 1), lty = 2, lwd = 1.5) # plot white without border
  polygon(ahuld.conv, col = scales::alpha("black", 0), # plot black border only
          border = scales::alpha("black", 0.4), lwd = 1.5)
  
  
  if(show_legend){
    legend("bottomright",
           c(
             "Downward heuristic",
             "Bullying",
             "Close competitors",
             "Undefined"
           ),
           # bty="n",
           col = "black",
           pt.bg = c(
             scales::alpha("black", 0.2),
             scales::alpha("blue", 0.2),
             scales::alpha("red", 0.2),
             scales::alpha("grey", 0.2)
           ),
           pch = 22, cex = 1, pt.cex = 1.5
    )
  }
  
  # add points and segments to upper layer
  
  # plot line for increasingly blurred downward heuristic, text above top of high position errorbars
  lines(x = blur_data$focus, y = blur_data$position)
  points(x = blur_data$focus, y = blur_data$position)
  text(x = blur_data$focus, y = blur_data$position_ci_hi + 0.025,
       label = blur_data$blur, cex = 0.6)
  
  # horizontal focus error, at blurred position
  graphics::segments(
    x0 = blur_data$focus_ci_lo, x1 = blur_data$focus_ci_hi,
    y0 = blur_data$position, y1 = blur_data$position,
    col = scales::alpha("black", 0.7), lwd = 1.5
  )
  
  # vertical position error, at blurred focus
  graphics::segments(
    x0 = blur_data$focus, x1 = blur_data$focus,
    y0 = blur_data$position_ci_lo, y1 = blur_data$position_ci_hi,
    col = scales::alpha("black", 0.7), lwd = 1.5
  )
  
  # add real data
  if (show_data_ci) {
    # horizontal focus error, at real position
    graphics::segments(
      x0 = data$focus_ci_lo, x1 = data$focus_ci_hi,
      y0 = data$position, y1 = data$position,
      col = scales::alpha("blue", 0.7), lwd = 2
    )
    
    # vertical position error, at real focus
    graphics::segments(
      x0 = data$focus, x1 = data$focus,
      y0 = data$position_ci_lo, y1 = data$position_ci_hi,
      col = scales::alpha("blue", 0.7), lwd = 2
    )
  }
  with(data, points(focus, position, col = "white", bg = "blue",
                    pch = 23, cex = 1.75))
  
}
