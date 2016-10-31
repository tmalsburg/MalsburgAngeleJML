### R code from vignette source '/home/malsburg/Documents/Uni/Projekte/AnalyzingETMeasures/MalsburgAngele2016JML/manuscript.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: load-libraries
###################################################

  library(knitr)
  library(dplyr)
  library(ggplot2)
  library(PropCIs)

  source("R/new_etmeasures.function.R")
  source("R/compute_model_parameters.function.R")



###################################################
### code chunk number 2: calculate-parameters
###################################################

  d1 <- read.csv("eyetracking_data/AngeleEtAl2013.csv")
  parameters1 <- with(d1, compute.model.parameters(subject, item, condition, ffd_n1, gzd_n1, gopast_n1, tvt_n1))
  
  d2 <- read.csv("eyetracking_data/MetznerEtAl2016.tsv", sep="\t", head=T)
  parameters2 <- with(d2, compute.model.parameters(subj, item, cond, FFD, FPRT, RPD, TFT))

  d <- data.frame(parameter=names(parameters1),
                  angele=unlist(parameters1),
                  metzner=unlist(parameters2))
  
  d[1,2:3] <- c(20, 50)
  d[2,2:3] <- c(20, 50)


###################################################
### code chunk number 3: table-parameters
###################################################

  kable(d[c(1:2, 5:6, 7:9, 3, 10, 12, 14, 4, 11, 13, 15),], digits=2, row.names=F, col.names=c("Parameter", "Lower bound", "Upper bound"), format="latex", booktabs=T)



###################################################
### code chunk number 4: measure-correlations
###################################################

cor.angele  <- cor(subset(d1, select=c(ffd_n1, gzd_n1, gopast_n1, tvt_n1)), use="complete")
cor.metzner <- cor(subset(d2, select=c(FFD, FPRT, RPD, TFT)), use="complete")

set.seed(1)
d1a <- new.etmeasures(parameters1)
set.seed(1)
d2a <- new.etmeasures(parameters2)

cor.angele.a  <- cor(subset(d1a, select=c(ffd, gzd, gpd, tvt)))
cor.metzner.a <- cor(subset(d2a, select=c(ffd, gzd, gpd, tvt)))

colnames(cor.angele)    <- rownames(cor.angele)    <- 
colnames(cor.metzner)   <- rownames(cor.metzner)   <- 
colnames(cor.angele.a)  <- rownames(cor.angele.a)  <- 
colnames(cor.metzner.a) <- rownames(cor.metzner.a) <- c("FFD", "GZD", "GPD", "TVT")

kable(cbind(cor.angele,  cor.angele.a),  digits=2, format="latex")
kable(cbind(cor.metzner, cor.metzner.a), digits=2, format="latex")


###################################################
### code chunk number 5: figure-scatter-plots-angele
###################################################

set.seed(1)
t <- as.logical(rbinom(nrow(d2), 1, 0.5))
par(mfrow=c(2, 3), mar=c(4,4,1,1))
with(d1[t,],  plot(ffd_n1, gzd_n1,    xlab="FFD", ylab="GZD", xlim=c(0, 800), ylim=c(0, 1500), cex=0.5, pch=20))
with(d1[t,],  plot(gzd_n1, gopast_n1, xlab="GZD", ylab="GPD", xlim=c(0, 800), ylim=c(0, 1500), cex=0.5, pch=20))
with(d1[t,],  plot(gzd_n1, tvt_n1,    xlab="GZD", ylab="TVT", xlim=c(0, 800), ylim=c(0, 1500), cex=0.5, pch=20))

with(d1a[t,], plot(ffd,    gzd,       xlab="Artificial FFD", ylab="Artificial GZD", xlim=c(0, 800), ylim=c(0, 1500), cex=0.5, pch=20))
with(d1a[t,], plot(gzd,    gpd,       xlab="Artificial GZD", ylab="Artificial GPD", xlim=c(0, 800), ylim=c(0, 1500), cex=0.5, pch=20))
with(d1a[t,], plot(gzd,    tvt,       xlab="Artificial GZD", ylab="Artificial TVT", xlim=c(0, 800), ylim=c(0, 1500), cex=0.5, pch=20))



###################################################
### code chunk number 6: figure-scatter-plots-metzner
###################################################

set.seed(1)
t <- as.logical(rbinom(nrow(d2), 1, 0.1))
par(mfrow=c(2, 3), mar=c(4,4,1,1))
with(d2[t,], plot(FFD, FPRT,  xlab="FFD", ylab="GZD", xlim=c(0, 2000), ylim=c(0, 4000), cex=0.5, pch=20))
with(d2[t,], plot(FPRT, RPD, xlab="GZD", ylab="GPD", xlim=c(0, 2000), ylim=c(0, 4000), cex=0.5, pch=20))
with(d2[t,], plot(FPRT, TFT,  xlab="GZD", ylab="TVT", xlim=c(0, 2000), ylim=c(0, 4000), cex=0.5, pch=20))

t <- as.logical(rbinom(nrow(d2a), 1, 0.1))
with(d2a[t,], plot(ffd, gzd, xlab="Artificial FFD", ylab="Artificial GZD", xlim=c(0, 2000), ylim=c(0, 4000), cex=0.5, pch=20))
with(d2a[t,], plot(gzd, gpd, xlab="Artificial GZD", ylab="Artificial GPD", xlim=c(0, 2000), ylim=c(0, 4000), cex=0.5, pch=20))
with(d2a[t,], plot(gzd, tvt, xlab="Artificial GZD", ylab="Artificial TVT", xlim=c(0, 2000), ylim=c(0, 4000), cex=0.5, pch=20))



###################################################
### code chunk number 7: load-simulation-results
###################################################

load("simulation_results/results.Rda")

# Multiply the p-values with the sign of the coefficient for the
# effect of condition if there was a true effect.  In the power
# analysis, we only consider effects that are significant and in the
# correct direction.
sims$ps <- with(sims, p * sign(mer1.cond))



###################################################
### code chunk number 8: one-effect
###################################################
# Rate of simulations in which at least one single-variable model
# shows an effect.  The uncorrected alpha level (0.05) is used to
# determine significance.

# Determine significance:

r <- sims %>%
  mutate(significant    = 0<=p  & p <=0.05,
         hit            = 0<=ps & ps<=0.05,
         significant.bf = 0<=p  & p <=0.05/4,
         hit.bf         = 0<=ps & ps<=0.05/4)

r.ffd <- r %>%
  group_by(effect.size, sim) %>%
  summarize(significant    = sum(significant),
            hit            = sum(hit),
            significant.bf = sum(significant.bf),
            hit.bf         = sum(hit.bf))

r.gzd <- r %>%
  filter(variable != "ffd") %>%
  group_by(effect.size, sim) %>%
  summarize(significant    = sum(significant),
            hit            = sum(hit),
            significant.bf = sum(significant.bf),
            hit.bf         = sum(hit.bf))

r.gpd <- r %>%
  filter(variable == "gpd") %>%
  select(sim, effect.size, significant, hit, significant.bf, hit.bf)

r.tvt <- r %>%
  filter(variable == "tvt") %>%
  select(sim, effect.size, significant, hit, significant.bf, hit.bf)

# Calculate rates and confidence intervals:

one.effect.ffd <- r.ffd %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2],
            hit.rate       = mean(hit>0),
            ht.ci.lower    = exactci(sum(hit>0),         n(), 0.95)$conf.int[1],
            ht.ci.upper    = exactci(sum(hit>0),         n(), 0.95)$conf.int[2])

two.effects.ffd <- r.ffd %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant>1),
            dr.ci.lower = exactci(sum(significant>1), n(), 0.95)$conf.int[1],
            dr.ci.upper = exactci(sum(significant>1), n(), 0.95)$conf.int[2],
            hit.rate = mean(hit>1),
            ht.ci.lower = exactci(sum(hit>1), n(), 0.95)$conf.int[1],
            ht.ci.upper = exactci(sum(hit>1), n(), 0.95)$conf.int[2])

bonferroni.effect.ffd <- r.ffd %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant.bf>0),
            dr.ci.lower = exactci(sum(significant.bf>0), n(), 0.95)$conf.int[1],
            dr.ci.upper = exactci(sum(significant.bf>0), n(), 0.95)$conf.int[2],
            hit.rate = mean(hit.bf>0),
            ht.ci.lower = exactci(sum(hit.bf>0), n(), 0.95)$conf.int[1],
            ht.ci.upper = exactci(sum(hit.bf>0), n(), 0.95)$conf.int[2])

one.effect.gzd <- r.gzd %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2],
            hit.rate       = mean(hit>0),
            ht.ci.lower    = exactci(sum(hit>0),         n(), 0.95)$conf.int[1],
            ht.ci.upper    = exactci(sum(hit>0),         n(), 0.95)$conf.int[2])

two.effects.gzd <- r.gzd %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant>1),
            dr.ci.lower = exactci(sum(significant>1), n(), 0.95)$conf.int[1],
            dr.ci.upper = exactci(sum(significant>1), n(), 0.95)$conf.int[2],
            hit.rate = mean(hit>1),
            ht.ci.lower = exactci(sum(hit>1), n(), 0.95)$conf.int[1],
            ht.ci.upper = exactci(sum(hit>1), n(), 0.95)$conf.int[2])

bonferroni.effect.gzd <- r.gzd %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant.bf>0),
            dr.ci.lower = exactci(sum(significant.bf>0), n(), 0.95)$conf.int[1],
            dr.ci.upper = exactci(sum(significant.bf>0), n(), 0.95)$conf.int[2],
            hit.rate = mean(hit.bf>0),
            ht.ci.lower = exactci(sum(hit.bf>0), n(), 0.95)$conf.int[1],
            ht.ci.upper = exactci(sum(hit.bf>0), n(), 0.95)$conf.int[2])

one.effect.gpd <- r.gpd %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2],
            hit.rate       = mean(hit>0),
            ht.ci.lower    = exactci(sum(hit>0),         n(), 0.95)$conf.int[1],
            ht.ci.upper    = exactci(sum(hit>0),         n(), 0.95)$conf.int[2])

two.effects.gpd <- r.gpd %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant>1),
            dr.ci.lower = exactci(sum(significant>1), n(), 0.95)$conf.int[1],
            dr.ci.upper = exactci(sum(significant>1), n(), 0.95)$conf.int[2],
            hit.rate = mean(hit>1),
            ht.ci.lower = exactci(sum(hit>1), n(), 0.95)$conf.int[1],
            ht.ci.upper = exactci(sum(hit>1), n(), 0.95)$conf.int[2])

bonferroni.effect.gpd <- r.gpd %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant.bf>0),
            dr.ci.lower = exactci(sum(significant.bf>0), n(), 0.95)$conf.int[1],
            dr.ci.upper = exactci(sum(significant.bf>0), n(), 0.95)$conf.int[2],
            hit.rate = mean(hit.bf>0),
            ht.ci.lower = exactci(sum(hit.bf>0), n(), 0.95)$conf.int[1],
            ht.ci.upper = exactci(sum(hit.bf>0), n(), 0.95)$conf.int[2])

one.effect.tvt <- r.tvt %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2],
            hit.rate       = mean(hit>0),
            ht.ci.lower    = exactci(sum(hit>0),         n(), 0.95)$conf.int[1],
            ht.ci.upper    = exactci(sum(hit>0),         n(), 0.95)$conf.int[2])

two.effects.tvt <- r.tvt %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant>1),
            dr.ci.lower = exactci(sum(significant>1), n(), 0.95)$conf.int[1],
            dr.ci.upper = exactci(sum(significant>1), n(), 0.95)$conf.int[2],
            hit.rate = mean(hit>1),
            ht.ci.lower = exactci(sum(hit>1), n(), 0.95)$conf.int[1],
            ht.ci.upper = exactci(sum(hit>1), n(), 0.95)$conf.int[2])

bonferroni.effect.tvt <- r.tvt %>%
  group_by(effect.size) %>%
  summarize(detection.rate = mean(significant.bf>0),
            dr.ci.lower = exactci(sum(significant.bf>0), n(), 0.95)$conf.int[1],
            dr.ci.upper = exactci(sum(significant.bf>0), n(), 0.95)$conf.int[2],
            hit.rate = mean(hit.bf>0),
            ht.ci.lower = exactci(sum(hit.bf>0), n(), 0.95)$conf.int[1],
            ht.ci.upper = exactci(sum(hit.bf>0), n(), 0.95)$conf.int[2])



###################################################
### code chunk number 9: decision-criteria
###################################################

all.criteria <- rbind(one.effect.ffd, two.effects.ffd, bonferroni.effect.ffd,
                      one.effect.gzd, two.effects.gzd, bonferroni.effect.gzd,
                      one.effect.gpd, two.effects.gpd, bonferroni.effect.gpd,
                      one.effect.tvt, two.effects.tvt, bonferroni.effect.tvt)

all.criteria$measure <- factor(rep(c("ffd", "gzd", "gpd", "tvt"), each=21),
                               c("ffd", "gzd", "gpd", "tvt"))

all.criteria$criterion <- factor(rep(c("One effect", "Two effects", "Bonferroni"), each=7),
                                 c("One effect", "Two effects", "Bonferroni"))



###################################################
### code chunk number 10: decision-criteria-table
###################################################
kable(subset(all.criteria, effect.size==0&measure=="ffd", select=c(criterion, detection.rate, dr.ci.lower, dr.ci.upper)), digits=3, format="latex")


###################################################
### code chunk number 11: figure-false-positives-by-parameters
###################################################

r0 <- subset(r.ffd, effect.size==0)[c("sim", "significant")]
s0 <- subset(sims, effect.size==0, select=c(sim, p.refix, p.regr, p.reread, mean.ffd, mean.gazediff, mean.gopastdiff, mean.tvtdiff))
s0 <- unique(s0)
r0 <- merge(r0, s0)

for (parameter.name in colnames(r0)[-(1:2)]) {
    
  r0[[parameter.name]] <- cut(
      r0[[parameter.name]],
      quantile(r0[[parameter.name]], seq(0, 1, 1/2)),
      c("low", "high"),
      T)
    
}

x1 <- r0 %>%
  group_by(p.refix) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2])
x2 <- r0 %>%
  group_by(p.regr) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2])
x3 <- r0 %>%
  group_by(p.reread) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2])
x4 <- r0 %>%
  group_by(mean.ffd) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2])
x5 <- r0 %>%
  group_by(mean.gazediff) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2])
x6 <- r0 %>%
  group_by(mean.gopastdiff) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2])
x7 <- r0 %>%
  group_by(mean.tvtdiff) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2])

colnames(x1)[1] <- colnames(x2)[1] <- colnames(x3)[1] <- colnames(x4)[1] <- colnames(x5)[1] <- colnames(x6)[1] <- colnames(x7)[1] <- "Parametervalue"
x <- rbind(x1, x2, x3, x4, x5, x6, x7)
x$parameter <- rep(colnames(r0)[-(1:2)], each=2)

x$detection.rate <- 100 * x$detection.rate
x$dr.ci.lower <- 100 * x$dr.ci.lower
x$dr.ci.upper <- 100 * x$dr.ci.upper

print(ggplot(x, aes(Parametervalue, detection.rate)) +
    geom_point() +
    geom_pointrange(aes(ymax=dr.ci.lower, ymin=dr.ci.upper)) +
    ylab("False positive rate (%)") +
    xlab("Parameter value") + 
    facet_wrap( ~ parameter, nrow=2, as.table=F) +
    theme_bw())



###################################################
### code chunk number 12: table:false-positives-max-effect
###################################################

xl <- subset(r0, p.regr=="low"  &  p.reread=="low"  &  mean.gopastdiff=="low")
xh <- subset(r0, p.regr=="high" &  p.reread=="high" &  mean.gopastdiff=="high")
xl$params <- "low"
xh$params <- "high"
x <- rbind(xl, xh)

xx <- x %>%
  group_by(params) %>%
  summarize(detection.rate = mean(significant>0),
            dr.ci.lower    = exactci(sum(significant>0), n(), 0.95)$conf.int[1],
            dr.ci.upper    = exactci(sum(significant>0), n(), 0.95)$conf.int[2])

kable(xx, format="latex", digits=3)
    


###################################################
### code chunk number 13: figure-power-by-criterion
###################################################

ggplot(filter(all.criteria, measure=="ffd")) +
  aes(factor(effect.size), hit.rate, group=criterion, linetype=criterion) +
  geom_line() +
  guides(linetype = guide_legend(title="Criterion:")) + 
  scale_x_discrete() +
  ylim(0, 1) +
  xlab("True effect (ms)") +
  ylab("Power") +
  theme_bw()



###################################################
### code chunk number 14: figure-power-by-stage
###################################################

# Plot how the power changes when effects are introduced at later stages:

ggplot(filter(all.criteria)) +
  aes(factor(effect.size), hit.rate, group=criterion, linetype=criterion) +
  geom_line() +
  guides(linetype = guide_legend(title="Criterion:")) + 
  scale_x_discrete() +
  ylim(0, 1) +
  xlab("True effect (ms)") +
  ylab("Power") +
  theme_bw() +
  facet_wrap(~measure, nrow=1, as.table=F)


###################################################
### code chunk number 15: figure-power-by-measure
###################################################
detection.rate.individual.models <- sims %>%
  group_by(effect.size, variable) %>%
  summarize(detection.rate=mean (0<=ps & ps<=0.05),
            ci.lower=exactci(sum(0<=ps & ps<=0.05), n(), 0.95)$conf.int[1],
            ci.upper=exactci(sum(0<=ps & ps<=0.05), n(), 0.95)$conf.int[2])

levels(detection.rate.individual.models$variable) <- toupper(levels(detection.rate.individual.models$variable))

ggplot(detection.rate.individual.models) +
  aes(factor(effect.size), detection.rate, group=variable, linetype=variable) +
  geom_line() +
  guides(linetype = guide_legend(title="Measure:")) + 
  scale_x_discrete() +
  ylim(0, 1) +
  xlab("True effect (ms)") +
  ylab("Power") +
  theme_bw()


