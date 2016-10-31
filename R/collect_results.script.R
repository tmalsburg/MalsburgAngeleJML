
# Script used to combine the results of the individual iterations into
# on big file.  Run this in the directory named simulation_results.

library(data.table)

files <- list.files(".", "^iteration_[^r].*")
nfiles <- length(files)

not.is.null <- function(x) !is.null(x)

sims <- lapply(1:nfiles, function(i) {

  load(files[i])

  message("Processing file ", i, " (", files[i], ") ...")

  x <- lapply(list(results.0, results.2.5, results.5, results.10, results.20, results.40, results.80),
              function (r) {
                t(sapply(c("ffd", "gzd", "gpd", "tvt"),
                       function (m) {
                         c(unlist(r[grep(m, names(r))]),
                           unlist(parameters))
                       }))
              })
  x <- do.call(rbind, x)
  x <- data.frame(x, row.names=NULL)
  x$sim <- files[i]
  x$variable <- c("ffd", "gzd", "gpd", "tvt")
  x$effect.size <- rep(c(0, 2.5, 5, 10, 20, 40, 80), each=4)

  x
})

message("Joining data frames ...")
sims <- data.frame(rbindlist(sims))

message("Renaming colums ...")
colnames(sims)[c(1:4)] <- c("p", "mer0.int", "mer1.int", "mer1.cond")

message("Rearranging columns ...")
sims <- sims[c(20, 21, 22, 1:19)]

message("Convert variable to factor ...")
sims$variable <- factor(sims$variable, levels=c("ffd", "gzd", "gpd", "tvt"))

message("Saving data ...")
save(sims, file="results.Rda")
