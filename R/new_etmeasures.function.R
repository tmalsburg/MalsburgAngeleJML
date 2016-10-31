
# Some functions that might come in handy when working with the
# log-normal distribution:
gmean  <- function(x, ...) exp(mean(log(x), ...))
gsd    <- function(x, ...) exp(sd(log(x), ...))
rlnorm <- function(n, gmean, gsd) exp(rnorm(n, log(gmean), log(gsd)))

#' Generates simulated eyetracking measures for a fictitious reading
#' experiment.  There are two conditions but no effects.  Also
#' missing: random slopes for items and subjects.  Given an
#' appropriate set of parameters the typical variances and
#' correlations between eyetracking measures are reasonably well
#' reproduced.
#'
#' @title Generate a artificial eyetracking measures for a fictitious reading experiment.
#' @param p a list containing the parameters used for simulating
#' eyetracking measures:
#' \describe{
#'   \item{n.subjects}{the number of subjects}
#'   \item{n.items}{the number of items}
#'   \item{p.refix}{the probability to refixate a word after the first fixation}
#'   \item{p.regr}{the probability to regress from a word during the first pass}
#'   \item{p.reread}{the probability to reread a word after the first pass}
#'   \item{mean.ffd}{mean duration of the first fixation}
#'   \item{mean.gazediff}{mean difference between the first fixation and the gaze duration}
#'   \item{mean.gopastdiff}{mean difference between the go past time (a.k.a. regression path duration) and the gaze duration}
#'   \item{mean.tvtdiff}{mean difference between total viewing time and gaze duration}
#'   \item{sd.ffd}{standard deviation of the duration of the first fixation}
#'   \item{sd.gazediff}{standard deviation of the difference between the first fixation and the gaze duration}
#'   \item{sd.gopastdiff}{standard deviation of the difference between the go past time (a.k.a. regression path duration) and the gaze duration}
#'   \item{sd.tvtdiff}{standard deviation of the difference between total viewing time and gaze duration}
#'   \item{sd.subjects}{standard deviation of the random intercepts for subjects}
#'   \item{sd.items}{standard deviation of the random intercepts for items}
#' }
#' @return data frame containing simulated eyetracking measures This
#' data frame has the following columns:
#'  \item{trial}{the trial id}
#'  \item{subj}{the subject id}
#'  \item{item}{the item id}
#'  \item{cond}{the condition (either -0.5 or 0.5)}
#'  \item{refix}{whether or not a refixation occurred}
#'  \item{regr}{whether or not a first pass regression occurred}
#'  \item{reread}{whether or not the word was reread after the first pass}
#'  \item{ffd}{the first fixation duration}
#'  \item{gzd}{the gaze duration}
#'  \item{gpd}{the go past time}
#'  \item{tvt}{the total viewing time}
#' @keywords eye movements, eyetracking measures, simulated data
#' @export
#' @examples
#' p <- list(n.subjects      = 40,
#'           n.items         = 132,
#'           p.refix         = 0.14,
#'           p.regr          = 0.07,
#'           p.reread        = 0.197,
#'           mean.ffd        = 221,
#'           mean.gazediff   = 194,
#'           mean.gopastdiff = 238,
#'           mean.tvtdiff    = 236,
#'           sd.ffd          = 1.3,
#'           sd.gazediff     = 1.4,
#'           sd.gopastdiff   = 1.8,
#'           sd.tvtdiff      = 1.5,
#'           sd.subjects     = 22,
#'           sd.items        = 21)
#' d <- new.etmeasures(p)
#' head(d)
new.etmeasures <- function(p) {

    n     <- p$n.subjects * p$n.items
    
    subj <- rep(1:p$n.subjects, each=p$n.items)
    item <- rep(1:p$n.items, p$n.subjects)

    # Condition labels following a Latin square.  List is assigned
    # based on subject id (odd or even):
    cond <- rep_len(c(-0.5,0.5), p$n.items)
    cond <- rep(cond, p$n.subj)
    cond <- (subj%%2 * 2 - 1) * cond
    d    <- data.frame(subj=subj, item=item, cond=cond)
    
    # Adding random intercepts and slopes for subjects:
    re.int.subj <- rnorm(p$n.subjects, 0, log(p$sd.subjects))
    re.int.subj <- rep(re.int.subj, each=p$n.items)

    # Adding random intercepts and slopes for items:
    re.int.item <- rnorm(p$n.items, 0, log(p$sd.items))
    re.int.item <- rep(re.int.item, p$n.subjects)

    # Sample ffds from log-normal:
    d$ffd <- rlnorm(n,
                    exp(log(p$mean.ffd) + re.int.subj + re.int.item),
                    p$sd.ffd)
    
    # Generate refixations, regressions, and rereads:
    d$refix  <- rbinom(n, 1, p$p.refix)
    d$regr   <- rbinom(n, 1, p$p.regr)
    d$reread <- rbinom(n, 1, p$p.reread)

    # Generate gaze durations, go-past times, and total viewing times:
    gazediff   <- rlnorm(n, p$mean.gazediff,   p$sd.gazediff)
    gopastdiff <- rlnorm(n, p$mean.gopastdiff, p$sd.gopastdiff)
    tvtdiff    <- rlnorm(n, p$mean.tvtdiff,    p$sd.tvtdiff)

    d$gzd <- d$ffd + ifelse(d$refix,  gazediff,   0)
    d$gpd <- d$gzd + ifelse(d$regr,   gopastdiff, 0)
    d$tvt <- d$gzd + ifelse(d$reread, tvtdiff,   0)

    d$trial <- 1:nrow(d)
    
    # Reorder the columns:
    d[c("trial", "subj", "item", "cond", "refix", "regr", "reread",
        "ffd", "gzd", "gpd", "tvt")]
    
}

