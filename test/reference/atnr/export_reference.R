# Reference outputs from ATNr (Gauzens et al.; CRAN, GPL) for FoodWebDynamics' ATN models.
# The package is run, not copied: its initialised parameters, right-hand sides and lsoda
# trajectories are written to CSV, and test_atnr.jl builds ATNModels from those parameters
# and must reproduce the derivatives and trajectories.
#
#   Rscript export_reference.R [LIBRARY]   # LIBRARY: an R library containing ATNr (1.1.2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) .libPaths(c(args[1], .libPaths()))
suppressMessages(library(ATNr))
out <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
cat("ATNr", as.character(packageVersion("ATNr")), "\n")

set.seed(20261116)
S <- 10
nb_b <- 3
repeat {                             # every consumer needs at least one resource
  BM <- sort(10^runif(S, 0, 4))
  fw <- suppressWarnings(create_Lmatrix(BM, nb_b, Ropt = 100))
  if (all(colSums(fw[, (nb_b + 1):S] > 0) > 0)) break
}
write.table(fw, file.path(out, "L.csv"), sep = ",", row.names = FALSE, col.names = FALSE)   # raw link strengths
fw[fw > 0] <- 1                      # rows resources, columns consumers
writeLines(as.character(BM), file.path(out, "BM.csv"))
write.table(fw, file.path(out, "fw.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
y_species <- runif(S, 0.5, 2)

dump_model <- function(name, m, y, times) {
  dir <- file.path(out, name)
  dir.create(dir, showWarnings = FALSE)
  for (n in ls(m)) {
    v <- tryCatch(m[[n]], error = function(e) NULL)
    if (is.numeric(v) && !(n %in% c("dB", "F", "log_BM"))) {
      write.table(as.matrix(v), file.path(dir, paste0(n, ".csv")), sep = ",", row.names = FALSE,
                  col.names = FALSE)
    }
  }
  m$initialisations()                # precomputations ODE needs (lsoda_wrapper does this itself)
  write.table(cbind(y, m$ODE(y, 0)), file.path(dir, "ode.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
  sol <- lsoda_wrapper(times, y, m, rtol = 1e-10, atol = 1e-12)
  write.table(sol, file.path(dir, "traj.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
}

m <- initialise_default_Scaled(create_model_Scaled(S, nb_b, BM, fw))
dump_model("scaled", m, y_species, seq(0, 50, length.out = 6))

m <- initialise_default_Unscaled(create_model_Unscaled(S, nb_b, BM, fw), temperature = 20)
dump_model("unscaled", m, y_species, seq(0, 5e6, length.out = 6))

m <- initialise_default_Unscaled_nuts(create_model_Unscaled_nuts(S, nb_b, 2, BM, fw),
                                      create_Lmatrix(BM, nb_b, Ropt = 100), temperature = 20)
dump_model("nutrients", m, c(runif(2, 5, 10), y_species), seq(0, 50, length.out = 6))
