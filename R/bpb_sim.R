# simulation functions for transmission histories for use with phybreak
#' @importFrom stats runif
#' @importFrom  graphics plot
sim.coal.tree.m <- function(tt = NULL,  # nolint: object_name_linter
                            names = NULL, donor = NULL, t_inf = NULL, t_sam = NULL,
                            a = rep(0, dim(tt)[1]), b = rep(3, dim(tt)[1]),
                            rhoD = rep(0, dim(tt)[1]),  # nolint: object_name_linter
                            rhoR = rep(0, dim(tt)[1]),  # nolint: object_name_linter
                            nSamples = NULL,  # nolint: object_name_linter
                            sample.times = NULL,  # nolint: object_name_linter
                            gen.rate = 365 / 1.5,  # nolint: object_name_linter
                            tr_window = rep(Inf, dim(tt)[1]), # can transmit after sampling
                            file_path = "",
                            tree_name = 1,
                            tree_name_digits = nchar(as.character(tree_name)),
                            save_tree = FALSE,
                            plot_tree = FALSE,
                            plot_colors = NULL,
                            leaf_raw_ids = NULL) {
  # multiple coalescent simulations
  prefix <- "[Simualte coalsecent tree]"
  # make sure there is a valid tree or way to build one
  if (!is.null(names) && !is.null(donor) && !is.null(t_inf) && !is.null(t_sam)) {
    nInd <- length(names)  # nolint: object_name_linter
    if (nInd < 2) {
      stop("There must be at least two individuals in the transmission tree")
    } else if (nInd != length(donor) || nInd != length(t_inf) || nInd != length(t_sam)) {
      stop("Lenth of 'names', 'donor', 't_inf', and 't_sam' must all be the same")
    } else if (!is.null(tt)) {
      stop("Conflicting inputs: must only input one of transmission tree or names, donor, t_inf, and t_sam")
    } else {
      tt <- data.frame(names, donor, t_inf, t_sam, stringsAsFactors = FALSE)
    }
  }
  if (is.null(tt)) {
    msg <- paste("No way to build transmission tree: Either a transmission tree",
                 "or names, donor, t_inf, and t_sam must be supplied")
    msg <- paste(prefix, msg)
    stop(msg)
  }


  # number of individuals in transmission tree
  nInd <- dim(tt)[1]  # nolint: object_name_linter

  # make sure inputs have the right number of elements
  if (length(a) != nInd) stop("a must be same length as number of individuals")
  if (length(b) != nInd) stop("b must be same length as number of individuals")
  if (length(rhoD) != nInd) stop("rhoD must be same length as number of individuals")
  if (length(rhoR) != nInd) stop("rhoR must be same length as number of individuals")

  if (any(c(a, b, rhoD, rhoR) < 0)) stop("Values of a, b, rhoD, and rhoR must be non-negative")

  # check whether (forwards time) transmission from recipient back to donor is ever possible
  rhoD_any <- any(rhoD > 0)  # nolint: object_name_linter

  # make sure nSamples and sample.times inputs are consistent
  if (is.null(nSamples) && is.null(sample.times)) {
    # default to one sample per individual if neither are provided
    nSamples <- rep(1, nInd)  # nolint: object_name_linter
  }
  if (is.null(nSamples)) { # sample.times provided but not nSamples
    if (!is.list(sample.times)) {
      stop("sample.times must be given as a list of vectors with the sample times of each sequence for each patient")
    } else if (length(sample.times) != nInd) {
      stop("sample.times must be given as a list of vectors with the sample times of each sequence for each patient")
    } else {
      # number of samples for each individual
      nSamples <- sapply(sample.times, length)  # nolint: object_name_linter
    }
  }
  if (is.null(sample.times)) { # all sample times from each individual
    sample.times <- list()  # nolint: object_name_linter
    for (i in 1:nInd) {
      sample.times[[i]] <- rep(tt$t_sam[i], nSamples[i])  # nolint: object_name_linter
    }
  }
  # make sure it is a list and lengths are correct if they are provded
  if (!is.null(nSamples) && !is.null(sample.times)) {
    if (!is.list(sample.times)) {
      stop("sample.times must be given as a list of vectors with the sample times of each sequence for each patient")
    } else if (length(sample.times) != nInd) {
      stop("sample.times must be given as a list of vectors with the sample times of each sequence for each patient")
    } else {
      for (i in 1:nInd) {
        if (length(sample.times[[i]]) != nSamples[i]) {
          stop("Length of sample.times[[i]] must be equal to nSamples[i]")
        }
      }
    }
  }



  # transmission tree (and other parameters) must be ordered from oldest to most recent infection date
  if (is.unsorted(tt$t_inf)) {
    tt_old <- tt
    new_order <- order(tt$t_inf)
    tt <- tt[new_order, ]
    leaf_raw_ids <- leaf_raw_ids[new_order]
    a <- a[new_order]
    b <- b[new_order]
    rhoR <- rhoR[new_order]  # nolint: object_name_linter
    rhoD <- rhoD[new_order]  # nolint: object_name_linter
    nSamples <- nSamples[new_order]  # nolint: object_name_linter
    sample.times <- sample.times[new_order]  # nolint: object_name_linter
    tr_window <- tr_window[new_order]
    donor_shuff <- tt$donor
    for (i in 1:nInd) {
      if (tt$donor[i] != 0) {
        tt$donor[i] <- which(tt$names == tt_old$names[donor_shuff[i]])
      }
    }
  }

  # find which individuals have at least 1 sample
  nSamples_nonzero <- which(nSamples > 0)  # nolint: object_name_linter

  # sample times as vector
  sample.times.unlist <- unlist(sample.times)  # nolint: object_name_linter
  # find all unique sample times for easier checking of values
  sample.times.unique <- unique(unlist(sample.times))  # nolint: object_name_linter

  # find times when transmission windows open
  tr_window_times <- tt$t_inf + tr_window

  # TODO make sure values are possible

  # starting guess for how many iterations may be needed (can be extended in loop)
  iter_guess <- length(sample.times.unique) + nInd # start guessing no migration events

  # vector for times of events
  t <- vector(length = iter_guess)
  # find oldest sampling time
  t[1] <- max(sample.times.unique)

  # initial lineages present
  lin_init <- -(1:sum(nSamples))
  # initial locations of each lineage (which patient each lineage is in)
  loc_init <- list()
  names <- character(0)
  leaf_times <- unlist(sample.times)
  for (i in 1:nInd) {
    # negative values imply samples were taken from that host at an earlier time
    loc_init[[i]] <- rep(-i, nSamples[i])
    # positive values imply samples are in that host at this time
    loc_init[[i]][which(sample.times[[i]] == t[1])] <- i
    names <- c(names, rep(tt$names[i], nSamples[i]))
  }
  # unlist into vector
  loc_init <- unlist(loc_init)

  lins <- data.frame(lin = lin_init, loc = loc_init)

  # current number of lineages
  k_all <- nSamples
  # lineages in each host (including ones that have already been sampled and can't coalesce)
  k <- integer(nInd) # lineages available to coalesce in each host
  for (i in 1:nInd) {
    k[i] <- length(which(loc_init == i))
  }
  k_tot <- sum(nSamples) # total number of lineages

  # names for leaves on newick tree
  names_newick <- paste(names, -lin_init, sep = "_")

  # declare matrix for internal nodes in tree
  nodes <- matrix(nrow = k_tot - 1, ncol = 2)

  int_node <- 0 # counter for number of internal nodes

  # initialize matrices for RVs
  z <- matrix(Inf, nrow = iter_guess, ncol = nInd)
  mD <- matrix(Inf, nrow = iter_guess, ncol = nInd)  # nolint: object_name_linter
  mR <- matrix(Inf, nrow = iter_guess, ncol = nInd)  # nolint: object_name_linter
  t_new_w <- matrix(Inf, nrow = iter_guess, ncol = nInd)
  t_new_s <- matrix(Inf, nrow = iter_guess, ncol = nInd)

  # initialize other things to keep track
  t_event <- numeric(iter_guess) # time for each event (not cumulative)
  event <- integer(iter_guess)

  # initialize list for migration events
  migrate <- vector("list", iter_guess)

  i <- 1 # index for number of events

  # loop as long as there is more than one lineage in longest infected patient
  # or it is during the time second patient is infected
  while (k_all[1] > 1 || (rhoD_any && t[i] - tt$t_inf[2] > -1e-15) || sum(k_all[2:nInd]) > 0) {
    # double lengths if necessary
    if (i == iter_guess) {
      t <- c(t, numeric(iter_guess))

      z <- rbind(z, matrix(Inf, nrow = iter_guess, ncol = nInd))
      mD <- rbind(mD, matrix(Inf, nrow = iter_guess, ncol = nInd))  # nolint: object_name_linter
      mR <- rbind(mR, matrix(Inf, nrow = iter_guess, ncol = nInd))  # nolint: object_name_linter
      t_new_w <- rbind(t_new_w, matrix(Inf, nrow = iter_guess, ncol = nInd))
      t_new_s <- rbind(t_new_s, matrix(Inf, nrow = iter_guess, ncol = nInd))

      t_event <- c(t_event, numeric(iter_guess))
      event <- c(event, integer(iter_guess))

      migrate <- c(migrate, vector("list", iter_guess))

      iter_guess <- 2 * iter_guess
    }

    k_atleast2 <- which(k >= 2)
    nk_atleast2 <- length(k_atleast2)
    z[i, k_atleast2] <- Fz(runif(nk_atleast2), k[k_atleast2], a[k_atleast2],
                           b[k_atleast2], (t[i] - tt$t_inf[k_atleast2]) * gen.rate) / gen.rate

    # (reverse time) migrations events from donor to recipient
    if (rhoD_any) {
      for (j in 1:nInd) {
        if (rhoD[j] > 0 && k[j] >= 1 &&
          # test to see if any recipients with donor j are in a contact window with j
          any(t[i] - tt$t_inf[which(tt$donor == j)] > 0 &&
            t[i] - tt$t_inf[which(tt$donor == j)] <= tr_window[which(tt$donor == j)])) {
          mD[i, j] <- Fm(runif(1), k[j], rhoD[j], a[j], b[j],  # nolint: object_name_linter
                         (t[i] - tt$t_inf[j]) * gen.rate) / gen.rate
        }
      }
    }

    # (reverse time) migrations events from recipient to donor
    t_t_inf_diffs <- t[i] - tt$t_inf
    k_nonzero <- which(k >= 1)
    for (j in k_nonzero) {
      if (j != 1 && t_t_inf_diffs[j] >= 0 && t[i] <= tr_window_times[j]) {
        if (rhoR[j] > 0) {
          mR[i, j] <- min(t_t_inf_diffs[j],    # nolint: object_name_linter
                          Fm(runif(1), k[j], rhoR[j], a[j], b[j],
                            (t[i] - tt$t_inf[j]) * gen.rate) / gen.rate)
        } else {
          mR[i, j] <- t_t_inf_diffs[j]    # nolint: object_name_linter
        }
      }
    }

    tr_window_times_t_diffs <- t[i] - tr_window_times
    t_new_w[i, tr_window_times_t_diffs > 0] <- tr_window_times_t_diffs[tr_window_times_t_diffs > 0]

    # time when a new sample is taken #i think there is a simpler way to do this
    for (j in nSamples_nonzero) {
      if (any(t[i] > sample.times[[j]])) {
        # is this after (forwards time) any samples have been taken from patient j?
        t_sam_diff <- t[i] - sample.times[[j]] # find difference between current time and sample times
        t_new_s[i, j] <- min(t_sam_diff[t_sam_diff > 0]) # minimum of values that are positive
      }
    }

    # which event happenes first within each class of event
    z_i <- which.min(z[i, ]) # index of first coalescence event
    z_t <- z[i, z_i] # time of first coalescence event
    mD_i <- which.min(mD[i, ])  # nolint: object_name_linter
    mD_t <- mD[i, mD_i]  # nolint: object_name_linter
    mR_t <- min(mR[i, ])  # nolint: object_name_linter
    mR_i <- which(mR[i, ] == mR_t)  # nolint: object_name_linter

    t_new_w_i <- which.min(t_new_w[i, ])
    t_new_w_t <- t_new_w[i, t_new_w_i]
    t_new_s_i <- which.min(t_new_s[i, ])
    t_new_s_t <- t_new_s[i, t_new_s_i]


    # which kind of event happens first
    event[i] <- which.min(c(z_t, mD_t, mR_t, t_new_w_t, t_new_s_t))
    # time to first event
    t_event[i] <- c(z_t, mD_t, mR_t, t_new_w_t, t_new_s_t)[event[i]]

    lins_prev <- lins # make copy of old lineages

    # first event is coalescence somewhere
    if (event[i] == 1) {
      int_node <- int_node + 1
      # coalescence takes place in individual z_i
      possible_lins <- lins_prev$lin[which(lins_prev$loc == z_i)] # lineages in z_i that can coalesce
      # choose two lineages to coalesce
      coal <- sample(possible_lins, 2, replace = FALSE)
      nodes[int_node, ] <- coal # add coalescence event to nodes in tree
      # add new internal lineage
      lins <- rbind(lins, c(int_node, z_i))
      # remove coalesced lineages
      lins <- lins[-which(lins_prev$lin == coal[1] | lins_prev$lin == coal[2]), ]
    } else if (event[i] == 3) {
      # was this the initial transmission event, where all lineages go to donor?
      if (any(t_event[i] == t[i] - tt$t_inf[mR_i])) {
        # move migrating lineages (reverse time) from recipient to donor
        for (j in 1:length(mR_i)) {  # nolint:seq_linter
          lins$loc[which(lins_prev$loc == mR_i[j])] <- tt$donor[mR_i[j]]
        }
      } else {
        possible_lins <- lins$lin[which(lins_prev$loc == mR_i)]
        migrate <- possible_lins[sample.int(length(possible_lins), 1)]
        # move migrating lineage (reverse time) from recipient to donor
        lins$loc[which(lins$lin == migrate)] <- tt$donor[mR_i]
      }
    }
    # first event is new transmission window opening/new sample
    if (event[i] == 4 || event[i] == 5 || any(abs(t[i] - sample.times.unique - t_event[i]) <= 1e-15)) {

      # see if event is new sample
      if (any(abs(t[i] - sample.times.unique - t_event[i]) <= 1e-15)) {
        new_sam_index <- which((t[i] - sample.times.unlist - t_event[i]) <= 1e-15)
        new_sam <- lin_init[new_sam_index]
        lins$loc[which(lins$lin %in% new_sam)] <- abs(lins$loc[which(lins$lin %in% new_sam)])
      }
    }

    # reduce current time by elapsed time during event
    t[i + 1] <- t[i] - t_event[i]

    k <- tabulate(lins$loc, nbins = nInd)
    k_all <- k + tabulate(-lins$loc, nbins = nInd)

    i <- i + 1
  }

  # if(dim(lins[[i]])[1] != 1){
  if (dim(lins)[1] != 1) {
    stop("Tree Not Fully Coalesced")
  }


  # times of coalescence events
  coal_times <- t[which(event == 1) + 1]

  # declare character vector for building up parts of newick tree
  branch <- character(int_node)

  # build newick tree
  for (i in 1:int_node) {
    if (all(nodes[i, ] < 0)) { # both edges go to leaves
      branch[i] <- paste0(
        "(", names_newick[-nodes[i, 1]], ":", leaf_times[-nodes[i, 1]] - coal_times[i], ",",
        names_newick[-nodes[i, 2]], ":", leaf_times[-nodes[i, 2]] - coal_times[i], ")"
      )
    } else if (nodes[i, 1] < 0 && nodes[i, 2] > 0) { # first edge goes to leaf, second goes to internal node
      branch[i] <- paste0(
        "(", names_newick[-nodes[i, 1]], ":", leaf_times[-nodes[i, 1]] - coal_times[i], ",",
        branch[nodes[i, 2]], ":", -coal_times[i] + coal_times[nodes[i, 2]], ")"
      )
    } else if (nodes[i, 1] > 0 && nodes[i, 2] < 0) { # first edge goes to internal node, second goes to leaf
      branch[i] <- paste0(
        "(", branch[nodes[i, 1]], ":", -coal_times[i] + coal_times[nodes[i, 1]], ",",
        names_newick[-nodes[i, 2]], ":", leaf_times[-nodes[i, 2]] - coal_times[i], ")"
      )
    } else if (all(nodes[i, ] > 0)) { # both edges go to internal nodes
      branch[i] <- paste0(
        "(", branch[nodes[i, 1]], ":", -coal_times[i] + coal_times[nodes[i, 1]], ",",
        branch[nodes[i, 2]], ":", -coal_times[i] + coal_times[nodes[i, 2]], ")"
      )
    }
  }

  # add semicolon to finish tree
  tree_newick <- paste0(branch[int_node], ";")


  # Node signs denote leaf or internal.
  # If negative, node is leaf. If positive, is internal.
  # Build geneology for internal use

  new_nodes <- nodes
  leaf_times <- leaf_times * 12  # Convert from months to years
  coal_times <- coal_times * 12  # Convert from months to years

  status <- sign(nodes)  # Get +1/-1 entries
  new_nodes[new_nodes > 0] <- new_nodes[new_nodes > 0] + nInd  # Shift forward
  # Flip forward the node indexes so they correspond to the leaf times
  new_nodes[new_nodes < 0] <- new_nodes[new_nodes < 0] * (-1)

  phylogeny <- matrix(data = 0, nrow = 2 * nInd, ncol = 7)  # Merge tree and distances
  # The 6th column is a flag for "leaf node"
  # The 7th column is the sample index (from input)
  # Set the index into the matrix
  phylogeny[1:(2 * nInd - 1)] <- 1:(2 * nInd - 1)

  for (row in 1:(dim(nodes)[1])) { # Loop over the rows
    for (col in 1:(dim(nodes)[2])) { # Loop over the columns. Should be of length 2
      index <- new_nodes[row, col]
      # Place the re-indexed merge in the tree.
      phylogeny[index, 2] <- row + nInd
      # Store the (input) sample index into column 7
      if (!is.na(leaf_raw_ids[index]) && leaf_raw_ids[index] > 0) {
        phylogeny[index, 7] <- leaf_raw_ids[index]
        # If we recorded a raw index, then we are a leaf node.
        phylogeny[index, 6] <- 1
      }
      outer_node_status <- status[row, col]  # Get a +1 or -1
      if (outer_node_status == -1) {
        outer_distance <- leaf_times[index]
      } else {
        outer_distance <- coal_times[index - nInd]
      }  # The inner node is always an interior node
      inner_distance <- coal_times[row]  # Get the coalescent time
      # Store the delta between absolute times
      phylogeny[index, 3] <- outer_distance
      phylogeny[index, 4] <- outer_distance - inner_distance
    }
  }

  phylogeny[, 5] <- phylogeny[, 4] / sum(phylogeny[, 4])  # normalize.
  # tip states
  tip_states <- vector()
  for (i in 1:nInd) {
    tip_states <- c(tip_states, rep(i - 1, nSamples[i]))
  }

  # make text object with newick tree and tip names/states
  result <- list("tree_newick" = tree_newick, "names_newick" = t(names_newick),
                 "tip_states" = t(tip_states), "phylogeny" = phylogeny)


  return(result)
}


# functions for inverse cdf sampling for coalescence and migration
# coalescence
Fz <- function(u, k, a, b, t) {  # nolint: object_name_linter
  Fz <- (1 - (1 - u)^(b / choose(k, 2))) * (a + b * t) / b # nolint: object_name_linter
  return(Fz)
}

# migration
Fm <- function(u, k, rho, a, b, t) {  # nolint: object_name_linter
  Fm <- (1 - (1 - u)^(b / (k * rho))) * (a + b * t) / b # nolint: object_name_linter
  return(Fm)
}
