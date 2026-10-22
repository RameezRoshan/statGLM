# Script to generate the simulated fully factorial harv dataset
# Calibrated against harv.xlsx with N = 1500, preserving genetic parameters

set.seed(20261022)

N_target <- 1500

# 1. Breeding Pedigree Structure
# 50 sires: Sires 1-25 in Stock 1, Sires 26-50 in Stock 2
# 80 dams: Sires 1-30 mated to 2 dams each (60 dams), Sires 31-50 mated to 1 dam each (20 dams)
n_sires <- 50
n_dams <- 80

sire_stock <- rep(1:2, each = 25)
sire_of_dam <- c(rep(1:30, each = 2), 31:50)
dam_stock <- sire_stock[sire_of_dam]

# Family sizes: unequal, unbalanced fecundity summing to exactly 1500
# Target mean = 1500 / 80 = 18.75 (ranging from 10 to 30)
base_sizes <- round(rnorm(n_dams, mean = 18.75, sd = 4.2))
base_sizes <- pmax(base_sizes, 9)
diff_n <- N_target - sum(base_sizes)
adj_indices <- sample(1:n_dams, abs(diff_n), replace = TRUE)
for (idx in adj_indices) {
  base_sizes[idx] <- base_sizes[idx] + sign(diff_n)
}
stopifnot(sum(base_sizes) == N_target)

animal_id <- 1:N_target
animal_dam <- rep(1:n_dams, times = base_sizes)
animal_sire <- sire_of_dam[animal_dam]
animal_stock <- sire_stock[animal_sire]

# 2. Factorial Experimental Factors (stock, batch, culturetype, sex, pond)
# Fully factorial crossing: all 16 combinations (stock x batch x culture x sex) populated
animal_batch <- integer(N_target)
animal_culture <- integer(N_target)
animal_sex <- integer(N_target)
animal_pond <- integer(N_target)
animal_age <- integer(N_target)

for (d in 1:n_dams) {
  idx <- which(animal_dam == d)
  n_d <- length(idx)
  
  # Distribute across batch 1 and 2
  b_assign <- sample(rep(1:2, length.out = n_d))
  animal_batch[idx] <- b_assign
  
  # Culture system 1 and 2
  c_assign <- sample(rep(1:2, length.out = n_d))
  animal_culture[idx] <- c_assign
  
  # Sex 1 (male) and 2 (female)
  s_assign <- sample(rep(1:2, length.out = n_d))
  animal_sex[idx] <- s_assign
  
  # Ponds nested in batch: Batch 1 has ponds 1-4, Batch 2 has ponds 5-8
  for (i in seq_along(idx)) {
    if (b_assign[i] == 1) {
      animal_pond[idx[i]] <- sample(1:4, 1)
      animal_age[idx[i]]  <- sample(c(180, 190, 200, 210, 220), 1)
    } else {
      animal_pond[idx[i]] <- sample(5:8, 1)
      animal_age[idx[i]]  <- sample(c(320, 330, 340, 350, 360), 1)
    }
  }
}

# 3. Simulate Genetic & Environmental Variances
# Calibrated to exact harv.xlsx parameters:
# Within sire covariance ~ 125.0
# Within dam covariance ~ 145.0 (sire 125 + dam 20)
# Residual environmental variance ~ 515.0
# Total phenotypic variance ~ 1374.0, Mean ~ 135.0
sigma2_s <- 125.0
sigma2_d <- 20.0
sigma2_e <- 515.0

u_sire <- rnorm(n_sires, mean = 0, sd = sqrt(sigma2_s))
u_dam  <- rnorm(n_dams, mean = 0, sd = sqrt(sigma2_d))
e_ind  <- rnorm(N_target, mean = 0, sd = sqrt(sigma2_e))

# Initial stocking weight bw0 (formerly bw1):
# Mean = 31.95 g, SD = 18.27 g
bw0_raw <- numeric(N_target)
for (i in 1:N_target) {
  mean_b <- if (animal_batch[i] == 1) 36.8 else 27.1
  bw0_raw[i] <- rnorm(1, mean = mean_b, sd = 17.5)
}
bw0_raw <- pmax(bw0_raw, 4.0)
bw0 <- round((bw0_raw - mean(bw0_raw)) / sd(bw0_raw) * 18.27 + 31.95, 1)

# Fixed effects parameters:
beta_bw0  <- 1.70
eff_stock <- ifelse(animal_stock == 2, 9.5, 0)
eff_batch <- ifelse(animal_batch == 2, 11.8, 0)
eff_cult  <- ifelse(animal_culture == 2, -8.2, 0)
eff_sex   <- ifelse(animal_sex == 2, -15.8, 0)

# Pond nested variation
pond_eff <- c(rnorm(4, 0, 7.5), rnorm(4, 0, 7.5))
names(pond_eff) <- 1:8
eff_pond <- pond_eff[animal_pond]

mu <- 85.0
bw_raw <- mu + beta_bw0 * (bw0 - 31.95) +
  eff_stock + eff_batch + eff_cult + eff_sex + eff_pond +
  u_sire[animal_sire] + u_dam[animal_dam] + e_ind

# Calibrate bw to exact target mean 134.92 and sd 37.07
bw_cal <- (bw_raw - mean(bw_raw)) / sd(bw_raw) * 37.07 + 134.92
bw <- round(bw_cal, 1)

# Assemble final data frame
harv <- data.frame(
  animal = animal_id,
  sire = factor(animal_sire),
  dam = factor(animal_dam),
  stock = factor(animal_stock),
  batch = factor(animal_batch),
  culturetype = factor(animal_culture),
  pond = factor(animal_pond),
  sex = factor(animal_sex),
  age = animal_age,
  bw0 = bw0,
  bw = bw
)

# Create data directory if not exists
if (!dir.exists("data")) {
  dir.create("data")
}

# Save as .rda with optimal compression
save(harv, file = "data/harv.rda", compress = "xz")
cat("Successfully generated and saved data/harv.rda (N = 1500)\n")

# Also save as CSV for accessibility
write.csv(harv, file = "data/harv.csv", row.names = FALSE)
cat("Successfully exported data/harv.csv\n")
