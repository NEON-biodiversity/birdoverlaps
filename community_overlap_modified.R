community_overlap <- function(traits, sp, discrete = FALSE, circular = FALSE, normal = TRUE, output = "median", weight_type= "hmean", randomize_weights = FALSE, unique_values = NULL, circular_args = list(), density_args = list(), hypervolume_set_args = list()) {

  # Return error if circular or discrete are specified with multivariate data.
  if (circular & 'matrix' %in% class(traits)) {
    stop("circular data types are not supported with multivariate data.")
  }
  if (discrete & 'matrix' %in% class(traits)) {
    stop("discrete data types are not supported with multivariate data.")
  }

  # Clean input, removing missing values and species with <2 values.
  sp <- as.character(sp)
  dat <- cbind(as.data.frame(traits), sp = sp)
  dat <- dat[stats::complete.cases(dat), ]
  abunds <- table(dat$sp)
  abunds <- abunds[abunds>1]
  dat <- dat[dat$sp %in% names(abunds), ]
  traitlist <- split(dat[, -ncol(dat)], dat$sp)
  uniquespp <- unique(dat$sp)
  nspp <- length(uniquespp)

  # Overlap cannot be calculated if there are less than 2 species with at least 2 individuals each.
  if (nspp < 2) return(NA)

  # Define common grid limits so that all density functions and hypervolumes are estimated across the same domain.
  grid_limits <- apply(dat[, -ncol(dat), drop = FALSE], 2, range)

  # Add a multiplicative factor to the upper and lower end of the range so that the tails aren't cut off.
  extend_grid <- c(-0.5, 0.5) %*% t(apply(grid_limits,2,diff))
  grid_limits <- grid_limits + extend_grid

  # In the discrete case for non-circular data, find the common unique values across which to build the histogram.
  if (discrete & !circular) {
    unique_values <- sort(unique(unlist(dat[, -ncol(dat)])))
  }

  # Create a list of univariate density functions (in univariate case) or hypervolumes (in multivariate case)
  density_list <- lapply(traitlist, function(x) Ostats:::trait_density(x, grid_limits, normal, discrete, circular, unique_values, density_args, circular_args))

  overlaps <- NULL
  abund_pairs <- NULL

  # All possible pairwise combinations.
  combs <- combn(1:nspp, 2)

  for (idx in 1:ncol(combs)) {

    overlaps <- c(overlaps, Ostats:::pairwise_overlap(density_list[[combs[1, idx]]], density_list[[combs[2, idx]]], discrete, density_args, hypervolume_set_args))
    if (weight_type == "hmean")
      abund_pairs <- c(abund_pairs, 2/(1/abunds[combs[1, idx]] + 1/abunds[combs[2, idx]]))
    if (weight_type == "mean")
      abund_pairs <- c(abund_pairs, (abunds[combs[1, idx]] + abunds[combs[2, idx]]))

  }

  if (randomize_weights == TRUE){
    abund_pairs <- sample(abund_pairs)}
  if (output == "median" && weight_type == "none"){
    final_output <- stats::median(overlaps)}
  if (output == "median" && weight_type == "hmean"){
    final_output <- matrixStats::weightedMedian(x = as.vector(overlaps), w = abund_pairs)}
  if (output == "median" && weight_type == "mean"){
    final_output <- matrixStats::weightedMedian(x = overlaps, w = abund_pairs)}
  if (output == "mean" && weight_type == "hmean"){
    final_output <- stats::weighted.mean(x = overlaps, w = abund_pairs)}
  if (output == "mean" && weight_type == "mean"){
    final_output <- stats::weighted.mean(x = overlaps, w = abund_pairs)}
  if (output == "mean" && weight_type == "none"){
    final_output <- mean(overlaps)}

  # Create a data frame of raw overlaps to return as well.
  overlaps_df <- as.data.frame(cbind(t(combs), overlaps))
  overlaps_df[,1] <- uniquespp[overlaps_df[,1]]
  overlaps_df[,2] <- uniquespp[overlaps_df[,2]]
  names(overlaps_df) <- c('species1', 'species2', 'overlap')


  return(list(value = final_output, raw = overlaps_df))
}
