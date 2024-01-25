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
  make_option("--sumstats", action = "store", default = NULL, type = "character"),
  make_option("--type", action = "store", default = "individual", type = "character"),
  make_option("--out", action = "store", default = "./", type = "character"),
  make_option("--N", action = "store", default = NULL, type = "character"),
  # required for PENGUIN
  make_option("--phen", action = "store", default = NULL, type = "character"),
  make_option("--phen.name", action = "store", default = NULL, type = "character"),
  make_option("--phen.col", action = "store", default = NULL, type = "character"),
  make_option("--covar", action = "store", default = NULL, type = "character"),
  make_option("--covar.name", action = "store", default = NULL, type = "character"),
  make_option("--covar.col", action = "store", default = NULL, type = "character"),
  # required for binary traits
  make_option("--ss.prev", action = "store", default = NULL, type = "character"),
  make_option("--ind.prev", action = "store", default = NULL, type = "character"),
  # required for PENGUIN-S
  make_option("--Ns", action = "store", default = NULL, type = "integer"),
  # optional - with defaults provided
  make_option("--gc", action = "store_true", default = FALSE),
  make_option("--hm3", action = "store", default = "./w_hm3.snplist", type = "character"),
  make_option("--ld", action = "store", default = "./eur_w_ld_chr/", type = "character"),
  make_option("--wld", action = "store", default = "./eur_w_ld_chr/", type = "character")
)

### Read in user inputs
opt <- parse_args(OptionParser(option_list = option_list))
# Required inputs
sumstat_files <- as.character(unlist(strsplit(opt$sumstats, split = ",")))
type <- opt$type
output_path <- opt$out
N <- as.numeric(unlist(strsplit(opt$N, split = ",")))
# input for PENGUIN - individual approach
if (type == "individual") {
  # read phenotype and covariate data and merge together
  if (is.null(opt$phen)) {
     stop("You didn't provide individual-level data, did you want to use the PENGUIN-S sumstats approach? If so, please pass the flag --type sumstats.")
  }
  # read phenotype data
  phen <- fread(opt$phen)
  if (!is.null(opt$phen.name)) {
    phen_columns <- as.character(unlist(strsplit(opt$phen.name, split = ",")))
    phen_columns <- match(phen_columns, names(phen))
  } else if (!is.null(opt$phen.col)) {
    phen_columns <- as.numeric(unlist(strsplit(opt$phen.col, split = ",")))
  } else {
      stop("You didn't provide the argument phen.name or phen.col. PENGUIN cannot decipher the phen file without this information. Please rerun PENGUIN with one of these arguments provided.")
  }
  colnames(phen)[phen_columns[1]] <- "Y"
  colnames(phen)[phen_columns[2]] <- "X"
  cols <- c("IID", "Y", "X")
  phen <- phen[, ..cols]
  # read covariate data
  if (!is.null(opt$covar)) {
    covariate <- fread(opt$covar)
    if (!is.null(opt$covar.name)) {
      covariate_columns <- as.character(unlist(strsplit(opt$covar.name, split = ",")))
    } else if (!is.null(opt$covar.col)) {
      covar.col <- as.numeric(unlist(strsplit(opt$covar.col, split = ",")))
      covariate_columns <- colnames(covariate)[covar.col]
    } else {
      stop("You provided a covariates file, but didn't provide the argument covar.name or covar.col. PENGUIN cannot decipher the covar file without this information. Please rerun PENGUIN with one of these arguments provided.")
    }
    cols <- c("IID", covariate_columns)
    covariate <- covariate[, ..cols]
    # merge phen and covariate together on IID
    phen <- base::merge(x = phen, y = covariate, by.x = c("IID"), by.y = c("IID"))
    phen <- na.omit(phen)
  } else {
    cat("\nWarning: you provided no covariates. PENGUIN will be run only with phenotype data.\n")
  }
} else {
  # input for PENGUIN-S - sumstats approach
  if (is.null(opt$Ns)) {
    stop("You must provide the argument Ns - the overlapping sample size of the two sumstats in order to run PENGUIN-S.")
  }
  Ns <- opt$Ns
}
gc <- opt$gc
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
save(LDSCoutput, file = paste0(output_path, "/LDSCoutput.RData"))
# Extract data from LDSCoutput that will be used in both PENGUIN and PENGUIN-S
x.var <- LDSCoutput$S[4]
gen.cov <- LDSCoutput$S[2]
k <- nrow(LDSCoutput$S)
se <- matrix(0, k, k)
se[lower.tri(se, diag = TRUE)] <- sqrt(diag(LDSCoutput$V))
se.var <- se[4]
se.cov <- se[2]

out_df <- data.frame(METHOD = character(), BETA = numeric(), SE = numeric(), P = numeric())

### Calculate phenotypic covariance between the exposure and outcome
cat("\nCalculating closed form solution for genetic confounding...\n")
if (type == "individual") {
  ### PENGUIN - Genetic confounding with INDIVIDUAL-LEVEL DATA
  cat("Using PENGUIN - calculating closed form solution with individual-level data.\n")

  # extract x and y residuals from lm fitted with covariates
  if (!is.null(opt$covar)) {
    y.res.formula <- as.formula(paste("Y ~", paste(covariate_columns, collapse = "+")))
    x.res.formula <- as.formula(paste("X ~", paste(covariate_columns, collapse = "+")))
    y.res <- resid(lm(y.res.formula, data = phen))
    x.res <- resid(lm(x.res.formula, data = phen))
  } else {
    y.res <- phen$Y
    x.res <- phen$X
  }
  y.res <- scale(y.res)
  x.res <- scale(x.res)
  dat <- data.frame(Y = y.res, X = x.res)

  # calculate beta
  xy.cov <- cov(dat$Y, dat$X)
  # apply genomic control for the cov(X, Y)
  if (gc) {
    xy.cov <- xy.cov / (LDSCoutput$I[1] * LDSCoutput$I[4])
  }
  cov.covxy <- xy.cov - gen.cov
  cov.varx <- var(dat$X) - x.var
  cov.beta <- cov.covxy / cov.varx

  # calculate standard error
  mu_x <- cov.covxy
  mu_y <- cov.varx
  sigma2_x <-  (1 + cov(dat$Y, dat$X)^2) / (nrow(dat) - 1) + (se.cov^2)
  sigma2_y <-  2 / (nrow(dat) - 1) + (se.var^2)
  se.est <- sqrt(sigma2_x / (mu_y^2) + (mu_x^2) * sigma2_y / (mu_y^4))
  # calculate p-value
  p.cov <- 2 * pnorm(abs(cov.beta / se.est), lower.tail = FALSE)
  penguin <- data.frame(METHOD = "PENGUIN", BETA = cov.beta, SE = se.est, P = p.cov)

  # calculate marginal regression results as a baseline to output if using individual level data
  lm.marg <- lm(Y ~ X, data = dat)
  se.marg <- summary(lm.marg)$coefficient[2, 2]
  beta.marg <- as.numeric(coefficients(lm.marg)[2])
  p.marg <- summary(lm.marg)$coefficient[2, 4]
  marginal <- data.frame(METHOD = "Marginal Reg", BETA = beta.marg, SE = se.marg, P = p.marg)

  # output data
  out_df <- rbind(out_df, penguin, marginal)
} else {
  ### PENGUIN-S - Genetic confounding with SUMSTATS DATA
  cat("Using PENGUIN-S - calculating closed form solution with sumstats data.\n")
  # calculate xy covariance from LDSC intercept
  ldsc.xycov <- (LDSCoutput$N[2] / Ns) * LDSCoutput$I[2]
  # apply genomic control for the cov(X, Y)
  if (gc) {
    ldsc.xycov <- ldsc.xycov / (LDSCoutput$I[1] * LDSCoutput$I[4])
  }
  # calculate beta
  cov.covxy <- ldsc.xycov - gen.cov
  cov.varx <- 1 - x.var
  cov.beta <- cov.covxy / cov.varx

  # calculate standard error
  ldsc.int <- strsplit(system(paste0("grep -E '^Cross trait Intercept:' ", output_path, "/penguin_ldsc.log"),intern = T),split="\\s+")[[1]]
  ldsc.se.int <- as.numeric(gsub(x = ldsc.int[5], pattern = "^\\((.*)\\)", replacement = "\\1"))
  ldsc.se.xycov <- ldsc.se.int * LDSCoutput$N[2] / Ns
  mu_x <- cov.covxy
  mu_y <- cov.varx
  sigma2_x <- (ldsc.se.xycov^2) + (se.cov^2)
  sigma2_y <- 2 / (Ns - 1) + (se.var^2)
  se.est <- sqrt(sigma2_x / (mu_y^2) + (mu_x^2) * sigma2_y / (mu_y^4))

  # calculate p-value
  p.cov <- 2 * pnorm(abs(cov.beta / se.est), lower.tail = FALSE)
  # output data
  penguins <- data.frame(METHOD = "PENGUIN-S", BETA = cov.beta, SE = se.est, P = p.cov)
  out_df <- rbind(out_df, penguins)
}

cat("\nGenetic confounding sumstats:\n")
print(out_df)
cat(paste0("\n\nResults written to ", output_path, "/results.txt"))
fwrite(out_df, paste0(output_path, "/results.txt"))
