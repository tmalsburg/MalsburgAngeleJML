
setwd("/home/malsburg/Documents/Uni/Projekte/AnalyzingETMeasures/eye-movements-false-alarms")

library(lme4)

source("R/new_etmeasures.function.R")
source("R/compute_model_parameters.function.R")

# Parameters ranges from the two experiments above but with subject
# and item numbres that are more realistic for a simple design with
# just two conditions:
parameter.ranges <- matrix(ncol=2, byrow=T,
                           data=c(
                             20, 50,                    # n.subjects
                             20, 50,                    # n.item
                             219.72406370, 232.3518431, # mean.ffd
                             1.30784904, 1.4299300,     # sd.ffd
                             1.10117695, 1.1572831,     # sd.subjects
                             1.07238407, 1.0814842,     # sd.items
                             0.14098156, 0.3234607,     # p.refix
                             0.07158487, 0.4297653,     # p.regr
                             0.19256018, 0.4051862,     # p.reread
                             197.09194420, 204.2718762, # mean.gazediff
                             1.40134663, 1.6932467,     # sd.gazediff
                             312.02408176, 558.2689046, # mean.gopastdiff
                             1.68605892, 1.8390409,     # sd.gopastdiff
                             242.25922238, 291.4012766, # mean.tvtdiff
                             1.54550643, 1.8689263))    # sd.tvtdiff

rownames(parameter.ranges) <- c("n.subjects", "n.items", "mean.ffd",
                                "sd.ffd", "sd.subjects", "sd.items",
                                "p.refix", "p.regr", "p.reread",
                                "mean.gazediff", "sd.gazediff",
                                "mean.gopastdiff", "sd.gopastdiff",
                                "mean.tvtdiff", "sd.tvtdiff")

sample.parameters <- function(parameter.ranges) {

  p <- lapply(1:nrow(parameter.ranges), function(i) {
    runif(1, parameter.ranges[i,1], parameter.ranges[i,2])
  })

  names(p) <- row.names(parameter.ranges)

  p[[1]] <- round(p[[1]])
  p[[2]] <- round(p[[2]])

  p
}

recordWarnings <- function(expr)
{
    warn = NULL
    frame_number <- sys.nframe()
    ans <- withCallingHandlers(expr, warning = function(w)
    {
      assign("warn", w, envir = sys.frame(frame_number))
      invokeRestart("muffleWarning")
    })
    list(ans, warn)
}

# Given a data set this function runs the analyses for all dependent
# measures.  Each analysis is run with a fixed effect for condition
# and without.  Likelihood-ratio tests are used to test significance
# of condition factor.  The result is a list containing the p-values
# for the effect of condition and the parameter estimates of the
# fixed-effects for all models.
analyze <- function(nd) {

  wlmer <- function(...) { recordWarnings(lmer(...)) }

  mer0.ffd <- wlmer(ffd ~ 1    + (1|item) + (1|subj), nd, REML=F)
  mer0.gzd <- wlmer(gzd ~ 1    + (1|item) + (1|subj), nd, REML=F)
  mer0.gpd <- wlmer(gpd ~ 1    + (1|item) + (1|subj), nd, REML=F)
  mer0.tvt <- wlmer(tvt ~ 1    + (1|item) + (1|subj), nd, REML=F)
  
  mer1.ffd <- wlmer(ffd ~ cond + (1|item) + (1|subj), nd, REML=F)
  mer1.gzd <- wlmer(gzd ~ cond + (1|item) + (1|subj), nd, REML=F)
  mer1.gpd <- wlmer(gpd ~ cond + (1|item) + (1|subj), nd, REML=F)
  mer1.tvt <- wlmer(tvt ~ cond + (1|item) + (1|subj), nd, REML=F)

  p.ffd <- anova(mer0.ffd[[1]], mer1.ffd[[1]])[2,8]
  p.gzd <- anova(mer0.gzd[[1]], mer1.gzd[[1]])[2,8]
  p.gpd <- anova(mer0.gpd[[1]], mer1.gpd[[1]])[2,8]
  p.tvt <- anova(mer0.tvt[[1]], mer1.tvt[[1]])[2,8]

  list(p.ffd = p.ffd,
       p.gzd = p.gzd,
       p.gpd = p.gpd,
       p.tvt = p.tvt,
       mer0.ffd = fixef(mer0.ffd[[1]]),
       mer0.gzd = fixef(mer0.gzd[[1]]),
       mer0.gpd = fixef(mer0.gpd[[1]]),
       mer0.tvt = fixef(mer0.tvt[[1]]),
       mer1.ffd = fixef(mer1.ffd[[1]]),
       mer1.gzd = fixef(mer1.gzd[[1]]),
       mer1.gpd = fixef(mer1.gpd[[1]]),
       mer1.tvt = fixef(mer1.tvt[[1]]),
       mer0.ffd.warn = mer0.ffd[[2]],
       mer0.gzd.warn = mer0.gzd[[2]],
       mer0.gpd.warn = mer0.gpd[[2]],
       mer0.tvt.warn = mer0.tvt[[2]],
       mer1.ffd.warn = mer1.ffd[[2]],
       mer1.gzd.warn = mer1.gzd[[2]],
       mer1.gpd.warn = mer1.gpd[[2]],
       mer1.tvt.warn = mer1.tvt[[2]])
}

niter <- 100     # number of iterations

# We cannot just add and substract half the effect on the millisecond
# scale because that would change the geometric mean (i.e., the
# intercept).  Given the gmean m and the effect size d we have to find
# a and b such that a-b=d and gmean(c(a,b)) = m.  Then a is the
# measurement in the condition coded with 0.5 and b the measurement
# for the condition coded with -0.5. The formulae below are derived by
# using gemean(c(a,b))=sqrt(a*b), substituting into a-b=c, and solving
# for a and b.
add.effect <- function(m, d, y, cond) {
  a <- (d + sqrt(d**2 + 4*m**2))/2
  dl <- log(a) - log(m)
  ifelse(cond<0, y-dl, y+dl)
}

system.time(
for (i in 1:niter) {

  message("Running iteration ", i, " out of ", niter, ".")
  
  p <- sample.parameters(parameter.ranges)
  nd <- new.etmeasures(p)
  
  # Analyze data set without effect of condition:

  nd$ffd <- log(nd$ffd)
  nd$gzd <- log(nd$gzd)
  nd$gpd <- log(nd$gpd)
  nd$tvt <- log(nd$tvt)
  results.0 <- analyze(nd)

  ffd <- nd$ffd
  gzd <- nd$gzd
  gpd <- nd$gpd
  tvt <- nd$tvt
  cond <- nd$cond

  # Analyze data set with a 2.5 ms effect of condition:

  nd$ffd <- add.effect(p$mean.ffd,                                                          2.5, ffd, cond)
  nd$gzd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff,                               2.5, gzd, cond)
  nd$gpd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$regr  *p$mean.gopastdiff, 2.5, gpd, cond)
  nd$tvt <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$reread*p$mean.tvtdiff,    2.5, tvt, cond)
  results.2.5 <- analyze(nd)

  # Analyze data set with a 5 ms effect of condition:

  nd$ffd <- add.effect(p$mean.ffd,                                                          5, ffd, cond)
  nd$gzd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff,                               5, gzd, cond)
  nd$gpd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$regr  *p$mean.gopastdiff, 5, gpd, cond)
  nd$tvt <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$reread*p$mean.tvtdiff,    5, tvt, cond)
  results.5 <- analyze(nd)
  
  # Analyze data set with a 10 ms effect of condition:

  nd$ffd <- add.effect(p$mean.ffd,                                                          10, ffd, cond)
  nd$gzd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff,                               10, gzd, cond)
  nd$gpd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$regr  *p$mean.gopastdiff, 10, gpd, cond)
  nd$tvt <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$reread*p$mean.tvtdiff,    10, tvt, cond)
  results.10 <- analyze(nd)
  
  # Analyze data set with a 20 ms effect of condition:

  nd$ffd <- add.effect(p$mean.ffd,                                                          20, ffd, cond)
  nd$gzd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff,                               20, gzd, cond)
  nd$gpd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$regr  *p$mean.gopastdiff, 20, gpd, cond)
  nd$tvt <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$reread*p$mean.tvtdiff,    20, tvt, cond)
  results.20 <- analyze(nd)
  
  # Analyze data set with a 40 ms effect of condition:

  nd$ffd <- add.effect(p$mean.ffd,                                                          40, ffd, cond)
  nd$gzd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff,                               40, gzd, cond)
  nd$gpd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$regr  *p$mean.gopastdiff, 40, gpd, cond)
  nd$tvt <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$reread*p$mean.tvtdiff,    40, tvt, cond)
  results.40 <- analyze(nd)
  
  # Analyze data set with a 80 ms effect of condition:

  nd$ffd <- add.effect(p$mean.ffd,                                                          80, ffd, cond)
  nd$gzd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff,                               80, gzd, cond)
  nd$gpd <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$regr  *p$mean.gopastdiff, 80, gpd, cond)
  nd$tvt <- add.effect(p$mean.ffd + nd$refix*p$mean.gazediff + nd$reread*p$mean.tvtdiff,    80, tvt, cond)
  results.80 <- analyze(nd)

  # Store the results of all analyses in a file:

  results.file <- tempfile("iteration_", tmpdir="./simulation_results")
  save(results.0, results.2.5, results.5, results.10, results.20,
       results.40, results.80, p, file=results.file)
  
})



