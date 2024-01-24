### R implementation of PENGUIN(S)
options(stringsAsFactors = FALSE)

### Require dependencies, and install automatically if user doesn't already have them installed
repos <- "https://cloud.r-project.org"
if (!require(data.table)) {
  install.packages("data.table", repos = repos)
  library(data.table)
}
if (!require(optparse)) {
  install.packages("optparse", repos = repos)
  library(optparse)
}
if (!require(GenomicSEM)) {
  if (!require(devtools)) {
    install.packages("devtools", repos = repos)
    library(devtools)
  }
  install_github("GenomicSEM/GenomicSEM")
  library(GenomicSEM)
}

### Read the arguments into R
option_list <- list(
  # required
  make_option("--sumstat_files", action = "store", default = NULL, type = "character"),
  make_option("--type", action = "store", default = "individual", type = "character"),
  make_option("--output", action = "store", default = "./", type = "character"),
  make_option("--N", action = "store", default = NULL, type = "character"),
  make_option("--Ns", action = "store", default = NULL, type = "integer"),
  # optional
  make_option("--individual_files", action = "store", default = NULL, type = "character"),
  make_option("--hm3", action = "store", default = "./w_hm3.snplist", type = "character"),
  make_option("--ld", action = "store", default = "./eur_w_ld_chr/", type = "character"),
  make_option("--wld", action = "store", default = "./eur_w_ld_chr/", type = "character")
)

### Read in user inputs
opt <- parse_args(OptionParser(option_list = option_list))
# Required inputs
sumstat_files <- as.character(unlist(strsplit(opt$sumstat_files, split = ",")))
type <- opt$type
output_path <- opt$output
N <- as.numeric(unlist(strsplit(opt$N, split = ",")))
Ns <- opt$Ns
# Optional inputs
if (type == "individual" && is.null(opt$individual_files)) {
  stop("You didn't provide individual-level data, did you want to use the PENGUIN-S sumstats approach? If so, please pass the flag --type sumstats")
}
if (!is.null(opt$individual_files)) {
  individual_files <- as.character(unlist(strsplit(opt$individual_files, split = ",")))
}
hm3 <- opt$hm3
ld <- opt$ld
wld <- opt$wld
# Munged files path
munged_files <- c(paste0(output_path, "/penguin1"), paste0(output_path, "/penguin2"))

### Run LD score regression for inputed GWAS sumstat files using GenomicSEM
cat("\nMunging input GWAS sumstats...\n")
munge(files = sumstat_files, N = N,
      hm3 = hm3, trait.names = munged_files,
      log.name = paste0(output_path, "/penguin"),
      parallel = TRUE, cores = 2)
## sample prev and population prev currently assuming continuous variables
cat("\nRunning LDSC...\n")
LDSCoutput <- ldsc(traits = paste0(munged_files, ".sumstats.gz")
      , sample.prev = c(NA, NA)
      , population.prev = c(NA, NA)
      , ldsc.log = paste0(output_path, "/penguin")
      , ld = ld, wld = wld)
# Extract data from LDSCoutput that will be used in both PENGUIN and PENGUIN-S
x.var <- LDSCoutput$S[4]
xy.cov <- LDSCoutput$S[2]
k <- nrow(LDSCoutput$S)
se <- matrix(0, k, k)
se[lower.tri(se, diag = TRUE)] <- sqrt(diag(LDSCoutput$V))
se.var <- se[4]
se.cov <- se[2]

### Calculate phenotypic covariance between the exposure and outcome
cat("\nCalculating closed form solution for genetic confounding...\n")
if (type == "individual") {
  ### PENGUIN - Genetic confounding with INDIVIDUAL-LEVEL DATA
  cat("Using PENGUIN - calculating closed form solution with individual-level data.\n")
  # read phenotype and covariates data and merge together
  phen <- fread(individual_files[1])
  colnames(phen)[3] <- "Y"
  phen <- subset(phen, select = -FID)
  covariate <- fread(individual_files[2])
  colnames(covariate)[3] <- "X"
  covariate <- subset(covariate, select = -FID)
  phen <- base::merge(x = phen, y = covariate, by.x = c("IID"), by.y = c("IID"))
  phen <- na.omit(phen)
  covariate_columns <- names(phen)[4:length(phen)]

  # extract x and y residuals from lm fitted with covariates
  y.res.formula <- as.formula(paste("Y ~", paste(covariate_columns, collapse = "+")))
  x.res.formula <- as.formula(paste("X ~", paste(covariate_columns, collapse = "+")))
  y.res <- resid(lm(y.res.formula, data = phen))
  x.res <- resid(lm(x.res.formula, data = phen))
  y.res <- scale(y.res)
  x.res <- scale(x.res)
  dat <- data.frame(Y = y.res, X = x.res)

  # calculate beta
  cov.covxy <- cov(dat$Y , dat$X) - xy.cov
  cov.varx <- var(dat$X) - x.var
  cov.beta <- cov.covxy / cov.varx

  # calculate standard error
  mu_x <- cov.covxy
  mu_y <- cov.varx
  sigma2_x <-  (1 + cov(dat$Y, dat$X)^2) / (nrow(dat) - 1) + (se.cov^2)
  sigma2_y <-  2 / (nrow(dat) - 1) + (se.var^2)
  se.est <- sqrt(sigma2_x / (mu_y^2) + (mu_x^2) * sigma2_y / (mu_y^4))

  # calculate marginal regression results as a baseline to output if using individual level data
  lm.marg <- lm(Y ~ X, data = dat)
  se.marg <- summary(lm.marg)$coefficient[2, 2]
  beta.marg <- as.numeric(coefficients(lm.marg)[2])
  p.marg <- 2 * pnorm(abs(beta.marg / se.marg), lower.tail = FALSE)

  # calculate p-value
  p.cov <- 2 * pnorm(abs(cov.beta / se.est), lower.tail = FALSE)
  # output data
  out <- c("BETA" = cov.beta, "SE" = se.est, "P" = p.cov, "BETA_Marginal_Regression" = beta.marg, "SE_Marginal_Regression" = se.marg, "P_Marginal_Regression" = p.marg)
} else {
  ### PENGUIN-S - Genetic confounding with SUMSTATS DATA
  cat("Using PENGUIN-S - calculating closed form solution with sumstats data.\n")
  # calculate xy covariance from LDSC intercept
  ldsc.xycov <- (LDSCoutput$N[2] / Ns) * LDSCoutput$I[2]

  # calculate beta
  cov.covxy <- ldsc.xycov - xy.cov
  cov.varx <- 1 - x.var
  cov.beta <- cov.covxy / cov.varx

  # calculate standard error
  ldsc.int <- strsplit(system(paste0("grep -E '^Cross trait Intercept:' ", output_path, "/penguin_ldsc.log"),intern = T),split="\\s+")[[1]]
  ldsc.se.int <- as.numeric(gsub(x = ldsc.int[5], pattern = "^\\((.*)\\)", replacement = "\\1"))
  mu_x <- cov.covxy
  mu_y <- cov.varx
  sigma2_x <- (LDSCoutput$N[2]^2) / (Ns^2) * (ldsc.se.int^2) + (se.cov^2)
  sigma2_y <-  2 / (Ns - 1) + (se.var^2)
  se.est <- sqrt(sigma2_x / (mu_y^2) + (mu_x^2) * sigma2_y / (mu_y^4))

  # calculate p-value
  p.cov <- 2 * pnorm(abs(cov.beta / se.est), lower.tail = FALSE)
  # output data
  out <- c("BETA" = cov.beta, "SE" = se.est, "P" = p.cov)
}

cat("\nGenetic confounding sumstats:\n")
print(out)
