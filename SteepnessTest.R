#Test of steepness package
library(EloSteepness)


data = PenB

testseq = elo_steepness_from_sequence(winner = data$Winner,
                                      loser = data$Loser,
                                      refresh = 0,
                                      cores = 1)
summary(testseq)

plot_steepness(testseq)
steepness_precis(testseq)
plot_scores(testseq)
plot_steepness_regression(testseq, width_fac = 0.2)


dataE = PenE
testseqE = elo_steepness_from_sequence(winner = dataE$Winner,
                                      loser = dataE$Loser,
                                      refresh = 0,
                                      cores = 1)
summary(testseqE)

plot_steepness(testseqE)
steepness_precis(testseqE)
plot_scores(testseqE)
plot_scores(testseqE, color = "red", subset_ids = "DV")
plot_steepness_regression(testseqE, width_fac = 0.2)
