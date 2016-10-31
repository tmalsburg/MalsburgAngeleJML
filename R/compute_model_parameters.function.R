#rm(list=ls())

rename_df_columns <- function(g){
  g$subj <- g$subject
  g$ffd <- g$ffd_n1
  g$gzd <- g$gzd_n1
  g$gpd <- g$gopast_n1
  g$tvt <- g$tvt_n1
  g
}


#setwd("~/Documents/Documents/Experiments/eye-movements-false-alarms/data/full_datasets_for_verification")

#g <- rename_df_columns(read.csv("Eureka_1_data_full.csv"))

require(lme4)
require(reshape)

compute.model.parameters <- function(subj, item, cond, ffd, gzd, gpd, tvt) {
  
  g <- data.frame(subj, item, cond, ffd, gzd, gpd, tvt, stringsAsFactors=F)
  
  g$cond <- factor(g$cond)

  # to ensure that the model intercept is the grand mean
  contrasts(g$cond) <- contr.helmert(length(levels(g$cond)))

  g.lm.ffd <- lmer(data = g, log(ffd) ~ cond + (1|subj) + (1|item))

  ps <- list()
  
  # Number of subjects and items:
  ps$n.subjects <- length(unique(g$subj))
  ps$n.items <- length(unique(g$item))
  
  ps$mean.ffd <- exp(coef(summary(g.lm.ffd))[1,1])
  ps$sd.ffd   <- exp(sd(resid(g.lm.ffd)))

  # Random intercept sd for subjects and items:
  ps$sd.subjects <- exp(attr(summary(g.lm.ffd)$varcor$subj, "stddev"))
  ps$sd.items    <- exp(attr(summary(g.lm.ffd)$varcor$item, "stddev"))
  
  g$refixation_duration <- g$gzd - g$ffd
  g$regression_path_only_duration <- g$gpd - g$gzd
  g$rereading_duration <- g$tvt - g$gzd

  refix  <- ifelse(!is.na(g$refixation_duration > 0), g$refixation_duration > 0, FALSE)
  regout <- ifelse(!is.na(g$regression_path_only_duration > 0), g$regression_path_only_duration > 0, FALSE)
  reread <- ifelse(!is.na(g$rereading_duration > 0), g$rereading_duration > 0, FALSE)

  ps$p.refix  <- mean(refix)
  ps$p.regr   <- mean(regout)
  ps$p.reread <- mean(reread)

  # only use cases where a refixation/regression/non-first-pass fixation was actually made
  g[g == 0] <- NA

  g.lm.refixation_duration <- lmer(data = g, log(refixation_duration) ~ cond + (1|subj) + (1|item))
  ps$mean.gazediff <- exp(coef(summary(g.lm.refixation_duration))[1,1])
  ps$sd.gazediff   <- exp(sd(resid(g.lm.refixation_duration)))


  g.lm.regression_path_only_duration <- lmer(data = g, log(regression_path_only_duration) ~ cond + (1|subj) + (1|item))
  ps$mean.gopastdiff <- exp(coef(summary(g.lm.regression_path_only_duration))[1,1])
  ps$sd.gopastdiff   <- exp(sd(resid(g.lm.regression_path_only_duration)))


  g.lm.rereading_duration <- lmer(data = g, log(rereading_duration) ~ cond + (1|subj) + (1|item))
  ps$mean.tvtdiff <- exp(coef(summary(g.lm.rereading_duration))[1,1])
  ps$sd.tvtdiff   <- exp(sd(resid(g.lm.rereading_duration)))

  ps

}

