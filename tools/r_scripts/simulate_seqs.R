# Created: 07/08/2017
# Author: Gustavo Barroso
# This script performs simulations of recombination landscape under diverse scenarios


##############################################
#
# General Parameters 
#
##############################################

sequence_length <- 100e+6
N0 <- 30000
Ne <- N0 # used when simulating fluctuating pop. sizes
little_r <- 1e-8 # recombination rate per site per generation
mean_rho <- 4 * N0 * little_r
little_mu <- 1.5e-8 # mutation rate per site per generation
mean_theta <- 4 * N0 * little_mu
num_haploids <- 68 # sample size (in haploids)


##############################################
#
# Demographic History 
#
##############################################

# For demographic history, the directory stays the same

# Introgression (1-pulse event between 2 populations)
admix <- FALSE
if(admix) {
  admix_proportion <- 0.1 # proportion of genetic material from source to target population
  admix_time <- 0.125 # time (in coalescent units) of secondary contact between populations
  split_time <- 2.0 # time (in coalescent units) of separation between populations
}  

# discretized time intervals used by iSMC
number_of_intervals <- 30
time_boundaries <- rep(0, number_of_intervals)
for(i in 2:number_of_intervals) {
  time_boundaries[i] = - log(1 - ((i - 1) / number_of_intervals))
}

# vector of time points (in coal. units) where demographic changes occur
time_demo_changes <- c(0.22, 1.1, 2.3)
# by how much the pop. size changes (past-wards at each step) at the above time points 
fold_changes <- c(1, 25, 3.3, 10)

# Now, since we have fluctuating pop. sizes, we must ajust N0 to have the (time-average) Ne that we want
# computes total (maximum) number of generations by assuming max. scaled coal. time is equal to max. time in discretisation
max_time <- time_boundaries[number_of_intervals] # can be viewed as ~expected maximum tree height in the ARG
num_gen_total <- 4 * Ne * max_time

# number of generations spent in-between each pop. size change
num_gen_ib <- 4 * Ne * time_demo_changes
num_gen_vec <- num_gen_ib[1]
for(i in 2:length(num_gen_ib)) {
  num_gen_vec[i] <- num_gen_ib[i] - num_gen_ib[i - 1]
}
num_gen_vec <- c(num_gen_vec, num_gen_total - num_gen_vec[length(num_gen_vec)])

# we solve for N0 using the harmonic mean formula:
harmonic_pop_sizes <- 0
for(i in 1:length(num_gen_vec)) {
  harmonic_pop_sizes <- harmonic_pop_sizes + num_gen_vec[i] / fold_changes[i]
}
N0 <- (harmonic_pop_sizes * Ne) / num_gen_total
  

# since SCRM sets theta and rho based on N0, we must update these values
# IMPORTANT: such rates are "genome-wide 'means'" (relative to mu and little_r), but NOT "time-wise 'means'" (relative to Ne). 
mean_theta <- 4 * N0 * little_mu 
mean_rho <- 4 * N0 * little_r

# plots demo. history
coal_times <- c(0.0, time_demo_changes, time_boundaries[number_of_intervals])
generation_times <- 4 * N0 * coal_times
pop_sizes <- c(fold_changes, fold_changes[length(fold_changes)])
pop_sizes <- pop_sizes * N0

pdf("demographic_history.pdf")
plot(x = generation_times, 
     y = pop_sizes, type = "s", lwd = 2.5, xaxt = "n", ylab = "Pop. Size",  xlab = "Time (4Ne Generations)")
# plot x-axis to match time boundaries (best if x-axis is plotted in log-scale):
x_axis <- round(4 * N0 * time_boundaries)
# plot x axis in equally spaced boundaries (best if x-axis is plotted in linear scale)
x_axis <- round(4 * N0 * seq(from = 0.0, to = max(time_boundaries), length.out = 40))
axis(1, at = x_axis, las = 2, cex.axis = 0.5)
dev.off()


# Note: for the case of human demography, we paste the commandline from simulation #4 in the supplemental material of Fu et al 2014:
#ms 68 1 -I 11 20 20 20 1 1 1 1 1 1 1 1 -en 0 1 1 -en 0 2 2.41428571428571 -en 0 3 3.23571428571429 -eg 0 2
#97.6909397920288 -eg 0 3 123.512032930849 -en 0 4 7.14285714285714e-11 -en 0 5 7.14285714285714e-11 -
#  en 0 6 7.14285714285714e-11 -en 0 7 7.14285714285714e-11 -en 0 8 7.14285714285714e-11 -en 0 9
#7.14285714285714e-11 -en 0 10 7.14285714285714e-11 -en 0 11 7.14285714285714e-11 -ej
#0.0321428571428571 5 4 -ej 0.0321428571428571 7 6 -en 0.0321430357142857 4 0.714285714285714 -en
#0.0321430357142857 6 0.714285714285714 -ej 0.0321446428571429 6 4 -en 0.0321464285714286 4
#0.714285714285714 -ej 0.0357142857142857 2 3 -ej 0.0357142857142857 4 3 -en 0.0357160714285714 3
#0.132857142857143 -en 0.0357160714285714 3 0.132857142857143 -es 0.0392857142857143 3 0.97 -en
#0.0392875 12 0.178571428571429 -en 0.0392875 3 0.132857142857143 -ej 0.0428571428571429 9 8 -ej
#0.0428571428571429 11 10 -en 0.0428573214285714 8 0.178571428571429 -en 0.0428573214285714 10
#0.178571428571429 -ej 0.0428589285714286 10 8 -en 0.0428607142857143 8 0.178571428571429 -ej
#0.0535714285714286 3 1 -en 0.0535732142857143 1 1 -ej 0.0714285714285714 12 8 -en 0.0714303571428571
#8 0.178571428571429 -en 0.107142857142857 1 0.521428571428571 -ej 0.214285714285714 8 1 -en
#0.2142875 1 0.714285714285714 -r 28000 50000000 -t 42000 -p 12 -seeds 46 47 48
  

##############################################
#
# Gamma model of spatial variation in Rho
#
##############################################

# Gamma distribution of rho
alpha_rho <- 0.5
beta_rho <- alpha_rho

# how often we change rho values along the genome (inverse of 'g' parameter described in Barroso et al.)
rho_transition_prob <- 1e-5

# starts simulation        
number_of_rho_transitions = 0
rho_transition_points = c(NULL)
current_transition_point = 0
rho_span_vector = c(NULL)

# gets the points where rho values change
while(current_transition_point < sequence_length) {
  number_of_rho_transitions <- number_of_rho_transitions + 1
  rho_span <- rgeom(1, rho_transition_prob)
  rho_span_vector <- c(rho_span_vector, rho_span)
  current_transition_point <- current_transition_point + rho_span
  rho_transition_points[number_of_rho_transitions] <- current_transition_point
}

# if the last transition point happens to fall outside the range of our sequence, we delete it
if(rho_transition_points[number_of_rho_transitions] > sequence_length){
  rho_transition_points <- rho_transition_points[- number_of_rho_transitions]
  number_of_rho_transitions <- number_of_rho_transitions - 1
  rho_span_vector <- rho_span_vector[- number_of_rho_transitions]
}

# gets the rho values from the gamma distribution
rho_values <- rgamma(n = number_of_rho_transitions + 1, shape = alpha_rho, rate = beta_rho)
rho_values <- append(rho_values, rho_values[length(rho_values)], after = length(rho_values))
rho_values <- rho_values * mean_rho # scale by genome-wide average rho to get a meaningful rate
first_rho <- rho_values[1]

# prints and writes rho landscape to file
rho_transition_points <- append(rho_transition_points, 0, after = 0)
rho_transition_points <- append(rho_transition_points, sequence_length, after = length(rho_transition_points))
# writes landscapes scaled back to (time-average) Ne, since this is how iSMC infers it
rho_landscape <- cbind(as.data.frame(rho_values * (Ne / N0)), as.data.frame(rho_transition_points))

pdf("sim_rho_landscape.pdf")
plot(y = rho_landscape$rho_values, x = as.integer(rho_landscape$rho_transition_points / 1e+3),
     type = "s", ylab = expression(rho), xlab = "Position (kb)", lwd = 0.5, main = "Gamma-simulated recombination landscape")
abline(h = mean_rho * (Ne / N0), lty = 2, col = "red")
dev.off()
write.table(rho_landscape, file = "rho_landscape.txt", quote = F, row.names = F, col.names = F)

# writes simulated parameters to file 
# we write mean_theta and mean_rho as a function of Ne, not N0, so we re-scale them
values <- c(mean_theta * (Ne / N0), mean_rho * (Ne / N0), alpha_rho, (1 - rho_transition_prob), little_r, N0, Ne)
names <- c("mean_theta", "mean_rho", "rho.alpha", "r_ii", "little_r", "N0", "Ne") 
sim_params <- cbind(names, values)
write.table(sim_params, file = "sim_params.txt", quote = F, row.names = F, col.names = F, sep = "\t")


##############################################
#
# Hotspot model of spatial variation in Rho
#
##############################################

if(getwd() != root_dir) {
  cat("Attempted to start simulation outside root directory! Moving back...")
  setwd(root_dir)
}
# appends to root directory
dir.create("Hotspot")
setwd("Hotspot")

# using different name because it is not expected to be genome-wide avg. rho in the hotspot case
background_rho <- mean_rho
initial_rho <- background_rho * sequence_length # arbitrarily start landscape in background state
# intensity of the hotspot relative to background
intensity <- 500 
# creates dir based on hotspot intensity
intensity_dir <- paste("intensity_", as.character(intensity), sep = "")
if(!dir.exists(intensity_dir)){
  dir.create(intensity_dir, showWarnings = F)
} else {
  print(paste("Dir", intensity_dir, "already exists! Over-writing..."))
  dir.create(intensity_dir, showWarnings = F)
}
setwd(intensity_dir)

# distribution of hotspots along the sequence
back_2_hot <- 1e-5 # probability of switching from background rho to hotspot
back_2_hot_dir <- paste("avg_bckgnd_span_", as.character((1 / back_2_hot) / 1e+3), "kb", sep = "")
hot_2_back <- 4e-4  # probability of switching from hotspot rho to background
hot_2_back_dir <- paste("avg_hot_span_", as.character((1 / hot_2_back) / 1e+3), "kb", sep = "")
# combines back_2_hot and hot_2_back params in 1 folder name
comb_dir <- paste(back_2_hot_dir, "_", hot_2_back_dir, sep = "")

if(!dir.exists(comb_dir)){
  dir.create(comb_dir, showWarnings = F)
} else {
  print(paste("Dir", comb_dir, "already exists! Over-writing..."))
  dir.create(comb_dir, showWarnings = F)
}
setwd(comb_dir)

# creates directory to store simulation files
dir.create("RhoSim")

# starts simulating landscape
number_of_rho_transitions = 0
rho_transition_points = c(NULL)
current_transition_point = 0
rho_span_vector = c(NULL)
current_state <- -1 # -1 for background, 1 for hotspot

# gets the points where rho values change
while(current_transition_point < sequence_length) {
  number_of_rho_transitions <- number_of_rho_transitions + 1
  if(current_state == -1) {
    rho_span <- rgeom(1, back_2_hot)
    rho_span_vector <- c(rho_span_vector, rho_span)
  }
  else if(current_state == 1) {
    rho_span <- rgeom(1, hot_2_back)
    rho_span_vector <- c(rho_span_vector, rho_span)
  }
  current_transition_point <- current_transition_point + rho_span
  rho_transition_points[number_of_rho_transitions] <- current_transition_point
  current_state = current_state * -1 # switches back and forth between background and hotspot
}

# if the last transition point happens to fall outside the range of our sequence, we delete it
if(rho_transition_points[number_of_rho_transitions] > sequence_length){
  rho_transition_points <- rho_transition_points[- number_of_rho_transitions]
  number_of_rho_transitions <- number_of_rho_transitions - 1
  rho_span_vector <- rho_span_vector[- number_of_rho_transitions]
}

# we start this vector with hotspot since it represents TRANSITION values, and the landscape starts with background
cold_hot <- rep(background_rho, length = length(rho_transition_points))
for(i in 1:length(cold_hot)){
  if(i %% 2 == 1) {
    cold_hot[i] <- background_rho * intensity
  }
}

rho_values <- append(cold_hot, background_rho, after = 0)
rho_values <- append(rho_values, rho_values[length(rho_values)], after = length(rho_values))
first_rho <- rho_values[1]
rho_transition_points <- append(rho_transition_points, 0, after = 0)
rho_transition_points <- append(rho_transition_points, sequence_length, after = length(rho_transition_points))
rho_landscape <- cbind(as.data.frame(rho_values), as.data.frame(rho_transition_points))

# computes mean rho across the genome (background is not the mean)
mean_rho <- 0
for(i in 1:(nrow(rho_landscape) - 1)) {
  span <- rho_landscape[i + 1, 2] - rho_landscape[i, 2] 
  mean_rho <- mean_rho + span * rho_landscape[i, 1]
}
mean_rho <- mean_rho / sequence_length
cat("FYI: mean_rho = ", as.character(mean_rho), sep = "")

# writes rho landscape to file
write.table(rho_landscape, file = "RhoSim/rho_landscape.txt", quote = F, row.names = F, col.names = F)

# writes simulated parameters to file
values <- c(mean_theta, background_rho, back_2_hot, hot_2_back, intensity, little_r, N0)
names <- c("mean_theta", "background_rho", "back_2_hot", "hot_2_back", "intensity", "little_r", "N0")
sim_params <- cbind(names, values)
write.table(sim_params, file = "RhoSim/sim_params.txt", quote = F, row.names = F, col.names = F, sep = "\t")

pdf("RhoSim/sim_rho_landscape.pdf")
plot(y = rho_landscape$rho_values, x = as.integer(rho_landscape$rho_transition_points),
     type = "s", ylab = expression(rho), xlab = "Position (bp)", lwd = 2, main = "Hotspots-simulated recombination landscape")
abline(h = mean_rho, col = "magenta", lty = 2)
dev.off()


##############################################
#
# Heterogeneous Theta landscape 
#
##############################################

# The shape of the distribution 
theta_dist <- "gamma" # or "uniform"

# If Gamma distribution of theta values 
alpha_theta <- 1
beta_theta <- alpha_theta
# If Uniform distribution of theta values
min_unif <- 0.1
max_unif <- 10

# how often we change theta along the sequence (inverse of 'f' parameter described in Barroso et al.)
theta_transition_prob <- 1e-3

number_of_theta_transitions = 0
theta_transition_points = c(NULL)
current_transition_point = 0
theta_span_vector = c(NULL)

# gets the points where theta values change
while(current_transition_point < sequence_length) {
  number_of_theta_transitions <- number_of_theta_transitions + 1
  theta_span <- rgeom(1, theta_transition_prob)
  theta_span_vector <- c(theta_span_vector, theta_span)
  current_transition_point <- current_transition_point + theta_span
  theta_transition_points[number_of_theta_transitions] <- current_transition_point
}

# if last transition != sequence length, we delete it
if(theta_transition_points[number_of_theta_transitions] >= sequence_length) {
  theta_transition_points <- theta_transition_points[- number_of_theta_transitions]
  number_of_theta_transitions <- number_of_theta_transitions - 1
  theta_span_vector <- theta_span_vector[- number_of_theta_transitions]
}

theta_values <- numeric(length = number_of_theta_transitions + 1)
# draws from the distribution chosen above
if(theta_dist == "gamma") {
  theta_values <- rgamma(n = number_of_theta_transitions + 1, shape = alpha_theta, rate = beta_theta)
} else if(theta_dist == "uniform") {
  theta_values <- runif(n = number_of_theta_transitions + 1, min = min_unif, max = max_unif)
} else {
  stop("Mis-specified distribution of theta values!")
}
theta_values <- append(theta_values, theta_values[length(theta_values)], after = length(theta_values))
# transforms to mean 1
#theta_values <- theta_values / mean(theta_values) 
theta_values <- theta_values * mean_theta
first_theta <- theta_values[1]

# assembles the theta landscape for reference
theta_transition_points <- append(theta_transition_points, 0, after = 0)
theta_transition_points <- append(theta_transition_points, sequence_length, after = length(theta_transition_points))
theta_landscape <- as.data.frame(cbind(theta_values * (Ne / N0),theta_transition_points))

# go back
pdf("sim_theta_landscape.pdf")
plot(y = as.numeric(theta_landscape[,1]), x = as.numeric(theta_landscape[,2]), lwd = 2,
     type = "s", ylab = expression(theta), xlab = "Position (bp)", main = "Theta landscape")
abline(h = mean_theta, lty = 2, col = "blue")
dev.off()

write.table(theta_landscape, file = "theta_landscape.txt", quote = F, row.names = F, col.names = F)

# writes simulated parameters to file
if(theta_dist == "gamma") {
  values <- c(alpha_theta, theta_transition_prob, little_mu, N0)
  names <- c("theta.alpha", "t_ij", "little_mu", "N0")
  sim_params <- cbind(names, values)
  write.table(sim_params, file = "sim_params_theta.txt", quote = F, row.names = F, col.names = F, sep = "\t")
} else if(theta_dist == "uniform") {
  values <- c(min_unif, max_unif, theta_transition_prob, little_mu, N0)
  names <- c("theta.unif.min", "theta.unif.max", "t_ij", "little_mu", "N0")
  sim_params <- cbind(names, values)
  write.table(sim_params, file = "sim_params_theta.txt", quote = F, row.names = F, col.names = F, sep = "\t")
} else {
  stop("Mis-specified distribution of theta values!")
}


##############################################
#
# Creating the SCRM command line
#
##############################################

# build an SCRM commandline based on the parameters used so far in the script
sink("scrm_commandline.txt")
if(exists("rho_landscape")) {
  cat(paste(as.character(num_haploids), "1 -r", as.character(first_rho * sequence_length),
            as.character(sequence_length), "", sep = " "))
  for (i in 2:(number_of_rho_transitions + 1)){
    cat("-sr ")
    cat(as.character(rho_transition_points[i]))
    cat(" ")
    cat(as.character(rho_values[i] * sequence_length))
    cat(" ")
  }
} else {
  cat(paste(as.character(num_haploids), "1 -r", as.character(mean_rho * sequence_length),
            as.character(sequence_length), "", sep = " "))
}
if(exists("theta_landscape")) {
  cat(paste("-t", as.character(first_theta * sequence_length), "", sep = " "))
  for(i in 2:(number_of_theta_transitions + 1)) {
    cat("-st ")
    cat(as.character(theta_transition_points[i]))
    cat(" ")
    cat(as.character(theta_values[i] * sequence_length))
    cat(" ")
  }
} else {
  cat(paste("-t", as.character(mean_theta * sequence_length)))
  cat(" ")
}
if(exists("time_demo_changes")) {
  for(i in 1:length(time_demo_changes)){
    cat("-eN ")
    cat(as.character(time_demo_changes[i]))
    cat(" ")
    cat(as.character(fold_changes[i]))
    cat(" ")
  }
}
if(exists("admix_proportion")) {
  cat("-es ")
  cat(as.character(admix_time))
  cat(" ")
  cat("1")
  cat(" ")
  cat(as.character(1 - admix_proportion))
  cat(" ")
  cat("-ej ")
  cat(as.character(split_time))
  cat(" ")
  cat("1")
  cat(" ")
  cat("2")
}
cat("-T") 
cat("\n")
sink()


##############################################
#
# Running SCRM 
#
##############################################

# reads generated SCRM commandline to map mutations into haplotypes
require(scrm)
require(stringr)

# This part of the script is compatible with R version > 3.6.0

command_line <- readLines("scrm_commandline.txt")

# number of ARGs to be generated (replicates) based on the simulated landscape
num_reps <- 10
# creates one dir per replicate
for(i in 1:num_reps) {
  rep_dir <- paste("rep_", as.character(i), sep = "")
  if(!dir.exists(rep_dir)){
    dir.create(rep_dir, showWarnings = F)
  } else {
    print(paste("Dir", rep_dir, "already exists! Over-writing..."))
    dir.create(rep_dir, showWarnings = F)
  }
}

pb <- txtProgressBar(min = 0, max = num_reps, style = 3)
for(i in 1:num_reps) {
  setTxtProgressBar(pb, i)
  full_arg <- scrm(command_line) # runs scrm
  haps <- as.data.frame(t(as.data.frame(full_arg$seg_sites))) # organises haplotypes
  pos <- row.names(haps) 
  pos <- str_replace(pos, "X", "")
  # gets infinite-sites SNPs positions
  num_pos <- as.numeric(str_replace(pos, "e.", "e-")) 
  # converts to finite loci coordinates
  snp_loci <- round(num_pos *  sequence_length, 0)
  # moves multiple hits to the right
  snp_loci[which(duplicated(snp_loci))] <- snp_loci[which(duplicated(snp_loci))] + 1
  # if it happ  ens again, we discard them
  if(length(which(duplicated(snp_loci))) > 0) {
    dup_pos <- which(duplicated(snp_loci))
    snp_loci <- snp_loci[-dup_pos]
    haps <- haps[-dup_pos,]
  }
  haps$pos <- snp_loci
  
  full_haplotyptes <- as.data.frame(matrix(ncol = num_haploids, nrow = sequence_length))
  full_haplotyptes[is.na(full_haplotyptes)] <- 0
  # maps mutations into haplotypes
  for(j in 1:num_haploids) {
    names(full_haplotyptes)[j] <- paste("hap_", as.character(j), sep = "")
    hap_seq <- numeric(length = sequence_length)
    # gets snp coordinates for haplotype j
    snp_pos <- haps[which(haps[,j] == 1), ncol(haps)]
    hap_seq[snp_pos] <- 1
    full_haplotyptes[,j] <- hap_seq
  }
  file_path <- paste(getwd(), "/rep_", as.character(i), "/rep_", as.character(i), sep = "")
  # writes with col. names as a header
  write.table(full_haplotyptes, file = file_path, row.names = F, col.names = T)
}
close(pb)

