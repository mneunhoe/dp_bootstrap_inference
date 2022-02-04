gaussian_scale <- function(epsilon, sensitivity, cdp = TRUE) {
  if (cdp) {
    # Adding N(0, sigma^2) noise satisfies GS/(2sigma^2)-zCDP
    # so N(0, (1/eps)^2) and GS=1 satisfies eps^2/2-zDP = rho-zCDP
    return(sensitivity / epsilon)
  } else {
    return(sensitivity * sqrt(2.0 * log(1.26 / delta)) / epsilon)
  }
}

### DP Functions ###
dpHistogram <-
  function(x,
           n,
           lower_bound,
           upper_bound,
           epsilon,
           bins) {
    # true histogram
    breaks <- seq(lower_bound, upper_bound, length.out = bins + 1)
    true_hist <- .Call(graphics:::C_BinCount, x, breaks, T, T)
    #true_hist = np.histogram(x, bins=bins, range=(lower_bound, upper_bound))
    # calculate laplace noise
    # noise = np.random.laplace(0, 2.0/epsilon, size=len(true_hist[0])) #histogram queries have sensitivity 2
    # caculate gaussian noise
    cdf_sensitivity <- 2.0
    noise <-
      rnorm(bins, 0, gaussian_scale(epsilon, cdf_sensitivity, cdp = TRUE))
    # noise = np.random.normal(0, common.gaussian_scale(epsilon, cdf_sensitivity, cdp=True), size=len(true_hist[0]))
    # add noise to bins
    #out = ([i+j for (i,j) in zip(true_hist[0], noise)], true_hist[1])
    out <- list(true_hist + noise, breaks)
    return(out)
  }


dpTree <-
  function(x,
           n,
           lower_bound,
           upper_bound,
           epsilon,
           depth,
           cdp) {
    # initialize list of number of bins at each level of tree
    #bins = [2**i for i in range(1,depth)]
    bins <- 2 ^ (1:depth)
    # divide up epsilon
    if (cdp) {
      eps <- epsilon / sqrt(depth)
    } else {
      eps <- epsilon / depth
    }
    # build noisy histogram at each level of tree
    #tree = [([n],(lower_bound, upper_bound))] + [dpHistogram(x, n, lower_bound, upper_bound, eps, bins[i]) for i in range(depth-1)]
    tree <-
      append(list(list(n, c(
        lower_bound, upper_bound
      ))), lapply(bins, function(b)
        dpHistogram(x, n, lower_bound, upper_bound, eps, b)))
    
    
    return(tree)
  }



inverseVariance <- function(epsilon) {
  # histogram has sensitivity 2
  b <- epsilon / 2
  return(2 * (b ^ 2))
}

# Create an array of elements that are adjacent e.g. adjacentElements([1,2,3,4]) returns [2,1,4,3]
adjacentElements <- function(ls) {
  adj <- vector(mode = "numeric", length = length(ls))
  i <- 0
  while (i < length(ls)) {
    if (i %% 2 == 1) {
      adj[i + 1] <- ls[(i + 1) - 1]
    } else {
      adj[i + 1] = ls[(i + 1) + 1]
    }
    i <- i + 1
  }
  return(adj)
}

wBelow <- function(tree) {
  #  """
  #    Recursive weight estimate from below. This assumes that the variance is the same for every node at in the tree.
  #    :param tree: Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
  #    :return: Single weight for each level of the tree.
  #    """
  weights <-
    vector(mode = "numeric", length = length(tree))  # initialize with one weight per level of the tree
  i <- length(tree)
  while (i >= 1) {
    if (i == length(tree)) {
      weights[i] = 1
    } else{
      prev = weights[i + 1]
      weights[i] = (2.0 * prev) / (2 * prev + 1)
    }#coerce to floats with 2.0
    i <- i - 1
  }
  return(weights)
}


countBelow <- function(tree, wBelows) {
  #  """
  #    Recursively compute counts from below
  #    :param tree: Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
  #    :param wBelows: Array of weights of same length as tree, where the ith weight corresponds to the ith weight calculated from below.
  #    :return: List of counts for each node of the tree.
  #    """
  counts = vector(mode = "list", length = length(tree))
  i = length(tree)
  while (i >= 1) {
    if (i == length(tree)) {
      counts[[i]] = tree[[i]][[1]]
    } else {
      w = wBelows[i]
      child = tree[[i + 1]][[1]]
      #childSum = [sum(x) for x in zip(child[0:len(child)-1],child[1:len(child)])] # sum all pairs of children counts
      #childSum = childSum[0:len(childSum):2] # pick out the sums that are children with the same parents
      childSum <-
        (child[1:(length(child) - 1)] + child[2:(length(child))])[seq(1, (length(child) -
                                                                            1), 2)]
      weightedT <- w * tree[[i]][[1]] #weigh parent counts
      weightedChild = (1 - w) * childSum #weigh child counts
      #counts[i] = [sum(x) for x in zip(weightedT, weightedChild)]
      counts[[i]] <- weightedT + weightedChild
    }#sum parent and child counts
    i  <- i - 1
  }
  return(counts)
}


wAbove <- function(tree, wBelows) {
  #  """
  #    Recursive weight estimation from above
  #    :param tree: Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
  #    :param wBelows: Array of weights of same length as tree, where the ith weight corresponds to the ith weight calculated from below.
  #    :return: Single weight for each level of the tree.
  #    """
  weights = vector(mode = "numeric", length = length(tree))
  i = 1
  while (i <= length(tree)) {
    if (i == 1) {
      weights[i] = 1
    } else {
      prevAbove = weights[i - 1]
      prevBelow = wBelows[i]
      weights[i] = 1.0 / (1.0 + (prevAbove + prevBelow) ^ (-1))
    }
    i <- i + 1
  }
  return(weights)
}

countAbove <- function(tree, countsBelow, wAboves) {
  #  """
  #    Recursively compute counts from above
  #    :param tree: Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
  #    :param countsBelow: Array of counts of same length as tree, as calculated by countBelow function
  #    :param wAboves: Weights computed from above, assuming each node in tree has same variance.
  #    :return: List of counts for each node of the tree.
  #    """
  counts = vector(mode = "list", length = length(tree))
  i = 1
  while (i <= length(tree)) {
    if (i == 1) {
      counts[[i]] = tree[[i]][[1]]
    } else {
      w = wAboves[i]
      #parents = [val for val in counts[i-1] for _ in (0, 1)] #replicating parent counts so dimensions match
      parents <- rep(counts[[i - 1]], each = 2)
      adjacents <- adjacentElements(tree[[i]][[1]])
      parentAdjDiff <-
        parents - adjacents # get difference between parent and adjacent node
      weightedT <- w * tree[[i]][[1]]  #weight current node count
      weightedPA = (1 - w) * parentAdjDiff #weighted parent - adjacent count
      counts[[i]] = weightedT + weightedPA
    }
    i <- i + 1
  }
  return(counts)
}


optimalCount <- function(tree, wA, cA, cB) {
  #  """
  #    Optimal counts for nodes of tree
  #    :param tree: Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
  #    :param wA: Array of weights calculated from above with wAbove
  #    :param cA: Array of counts calculated from below with countBelow
  #    :param cB: Array of counts calculated from above with countAbove
  #    :return: Optimized tree
  #    """
  counts = vector(mode = "list", length = length(tree))
  i = 1
  while (i <= length(tree)) {
    if (i == 1) {
      counts[[i]] = tree[[i]][[1]]
    } else {
      w = wA[i]
      parents = rep(cA[[i - 1]], each = 2) #replicating parent counts so dimensions match
      adjacents = adjacentElements(cB[[i]])
      parentAdjDiff = parents - adjacents # get difference between parent and adjacent node
      weightedT = w * cB[[i]] #weight current node count
      weightedPA = (1 - w) * parentAdjDiff #weighted parent - adjacent count
      counts[[i]] = weightedT + weightedPA
    }
    i <- i + 1
  }
  
  for (i in seq_along(tree))
    tree[[i]][[1]] <- counts[[i]]
  
  return(tree)
}

optimalPostProcess <- function(tree, epsilon) {
  #  """
  #     Optimal Post Processing
  #
  #     Wrapper function that generates optimal tree from noisy tree generated with the Laplace mechanism. The general idea is that
  #     you can leverage the fact that the child counts at each node should sum to the parent count, and a node's count should be
  #     equal to its parent count minus the count of the adjacent child node to generate less noisy counts at every node of the tree.
  #
  #     You can think of the leveraging of child node's information to get the parent node's counts as a recursive process that
  #     begins at the leaves of a tree, which here we refer to with the tag "Below" in the helper function names. Similarly, leveraging
  #     a parent node and the adjacent child nodes can be thought of as a recursive process that begins with the root node, which
  #     is referred to here with the tage "Above" in the helper function names. A new count at every node in the tree can then be calculated
  #     using the counts that are generated in this way, which each contribute an amount according to some weight, which are calculated
  #     by the wBelow and wAbove functions respectively.
  #
  #    The theory behind this is explored in detail in the Honaker paper, whose implementation here is described in extra_docs/tree-post-processing.
  #     The implementation here assumes that the variance of the noise added is equal at every level of the tree, which also means that the weights
  #     wBelow and wAbove are the same for nodes at the same level of the tree. Honaker provides a more general implementation in his work.
  
  #     reference: Honaker, James. "Efficient Use of Differentially Private Binary Trees." (2015).
  
  #    :param tree: Differentially private tree generated from dpTree$release method
  #    :param epsilon: The epsilon value used for the noise addition at each node (note this is not the same as the global epsilon value.)
  #    """
  wB = wBelow(tree) # calculate the weights of counts from below
  wA = wAbove(tree, wB) # calculate the weights of counts from above
  cB = countBelow(tree, wB) # calculate counts from below
  cA = countAbove(tree, cB, wA) # calculate counts from above
  
  return(list(optimalCount(tree, wA, cA, cB), wA, wB))
}


### Post-processed CDF
postProcessedCDF <- function(tree, epsilon, monotonic = F) {
  # Commenting out node effects code because we can pre-process that
  # smallest granularity cdf possible uses leaf buckets
  
  # tree1 = list(tree)
  tree1 = tree
  vals = tree1[[length(tree1)]][[2]]
  counts <- NULL
  
  i = 1
  
  # generate trees of node effects
  depth = length(tree1)
  # node_effect_trees = node_effects_efficient(depth, epsilon)
  
  while (i <= length(vals)) {
    # build tree to keep track of the effects of each node on the calculation
    # effects_tree = empty_tree(len(tree1))
    # start out with min and max values at top of tree
    mM = tree1[[1]][[2]]
    # initialize count for i
    count = 0
    # iterate through layers of tree
    index = 0
    j = 1
    while (j <= length(tree1)) {
      # determine if should traverse tree to left or right
      mid = mM[1] + (mM[2] - mM[1]) / 2
      # if looking at leftmost node in the tree, we know empirical cdf should evaluate
      # to 0 not to bin size.
      if (i == 1) {
        break
      }
      # if you don't need higher granularity, stop traversal
      if (vals[i] == mM[2]) {
        count <- count + tree1[[j]][[1]][index + 1]
        # effects = effects_on_node(depth, [j, index], node_effect_trees)
        # effects_tree = sum_trees(effects_tree, effects)
        break
      }
      # if at leaves of tree, record the count there
      if (j == length(tree1)) {
        count <- count + tree1[[j]][[1]][index + 1]
        # effects = effects_on_node(depth, [j, index], node_effect_trees)
        # effects_tree = sum_trees(effects_tree, effects)
        # if traversing left
      } else if (vals[i] <= mid) {
        # reset max value to the mid, don't add to the count
        mM[2] = mid
        # set next index of node to look at in next layer
        index = index * 2
      } else{
        # if at end of tree, record count at that node
        #if traversing right,
        
        # reset min value to the mid
        mM[1] = mid
        # set to next index of node to look at in next layer
        index = index * 2 + 1
        count <-
          count +  tree1[[j + 1]][[1]][index] # add the node's left child to the count
        # effects = effects_on_node(depth, [j+1, index-1], node_effect_trees)
        # effects_tree = sum_trees(effects_tree, effects)
      }
      j <- j + 1
    }
    counts <- c(counts, count)
    
    
    
    ## Calculate variance of the count
    # variance = variance_from_effects(effects_tree, epsilon)
    i <- i + 1
  }
  n = tree1[[1]][[1]][1] # pull public n from root of tree
  #Prevent negative counts
  # for(i in 2:length(counts)){
  #   counts[i] <- ifelse(counts[i]<counts[i-1], counts[i-1], counts[i])
  # }
  percents = counts / n
  
  percents[percents < 0] <- 0
  percents[percents > 1] <- 1
  #percents = counts/max(counts) # normalize counts by n
  #percents <- percents/max(percents)
  # variances.append(variance/n**2) # scaling by the normalization
  # variances.append(variance)
  
  return(percents)
}





dpCDF <-
  function(x,
           lower_bound,
           upper_bound,
           epsilon,
           granularity,
           cdp = TRUE,
           num_trials = 1) {
    #  """
    #    :param x: List or numpy array of real numbers
    #    :param lower_bound: Lower bound on values in x.
    #    :param upper_bound: Upper bound on values in x.
    #    :param epsilon: Desired value of epsilon.
    #    :param hyperparameters: dictionary that contains:
    #        :param depth: depth of the tree desired.
    #    :param num_trials: Number of fresh DP median estimates to generate using dataset x.
    #    :return: List of 1/2 eps-DP estimates for the median of x.
    #    """
    
    n_bins = (upper_bound - lower_bound) %/% granularity + 1
    depth = log2(n_bins) %/% 1
    
    results = vector("list", length = num_trials)
    for (i in 1:num_trials) {
      t = dpTree(x, length(x), lower_bound, upper_bound, epsilon, depth, cdp)
      tOpt = optimalPostProcess(t, epsilon)[[1]]
      
      cdf = postProcessedCDF(tOpt, epsilon)
      
      results[[i]] <- list(cdf, t[[length(t)]][[2]])
    }
    return(results)
  }


single_value_tree <- function(depth, one_loc) {
  # initialize list of number of bins at each level of tree
  #bins = [2**i for i in range(1,depth)]
  bins <- 2 ^ (1:depth)
  
  
  tree <-
    append(list(list(0)), lapply(bins, function(b)
      list(rep(0, b))))
  
  
  tree[[one_loc[1]]][[1]][one_loc[2]] <- 1
  
  
  
  return(tree)
}

empty_tree <- function(depth) {
  # initialize list of number of bins at each level of tree
  #bins = [2**i for i in range(1,depth)]
  bins <- 2 ^ (1:depth)
  
  
  tree <-
    append(list(list(0)), lapply(bins, function(b)
      list(rep(0, b))))
  
  
  
  return(tree)
}



sum_trees <- function(tree_1, tree_2) {
  lapply(1:length(tree_1), function(x)
    list(tree_1[[x]][[1]] + tree_2[[x]][[1]]))
}


flip_level <- function(level, index, node_vals) {
  bin_rep <- rev(as.integer(intToBits(index)))
  
  if (sum(bin_rep) > 0) {
    st <- paste0(bin_rep[min(which(bin_rep == 1)):32], collapse = "")
  } else {
    st <- paste0(0)
  }
  
  index_bin_string <- stringr::str_pad(st, level, pad = "0")
  
  start <- 0
  end <- length(node_vals)
  
  for (c in strsplit(index_bin_string, "")[[1]]) {
    mid <- as.integer((start + end) / 2)
    if (c == 1) {
      flipped_vals <- node_vals[(mid + 1):end]
      flipped_vals <- c(flipped_vals, node_vals[(start + 1):mid])
      for (i in seq(length(node_vals))) {
        if (i > start & i <= end) {
          node_vals[i] <- flipped_vals[i - start]
        } else {
          node_vals[i] <- node_vals[i]
        }
      }
      start <- mid
    } else {
      end <- mid
    }
  }
  return(node_vals)
}

node_effects <- function(depth, epsilon) {
  # returns a list of trees that describe the effect that a node at each level of the tree effects the measurements at any other node.
  # figure out how many nodes there are at each level of the tree
  nodes_per_level = 2 ^ (0:depth)
  
  # build all the trees of 0's with single 1s
  all_input_trees = list()
  for (i in 1:depth) {
    input_trees_at_level = list()
    for (j in 1:nodes_per_level[i]) {
      input_trees_at_level <-
        append(input_trees_at_level, list(single_value_tree(depth, c(i, j))))
    }
    all_input_trees <-
      append(all_input_trees, list(input_trees_at_level))
  }
  
  
  node_effect_trees = list()
  for (level in all_input_trees) {
    node_effect_trees_at_level = list()
    
    for (input_tree in level) {
      tOpt = optimalPostProcess(input_tree, epsilon)[[1]]
      node_effect_trees_at_level <-
        append(node_effect_trees_at_level, list(tOpt))
    }
    node_effect_trees <-
      append(node_effect_trees, list(node_effect_trees_at_level))
  }
  
  return(node_effect_trees)
  
  
}


effects_on_node <-
  function(depth, node_location, node_effect_trees) {
    # tree that displays how much a particular node is affected by all other node's counts
    # node location is entered as [depth, location inside level from the left]
    
    # figure out how many nodes there are at each level of the tree
    nodes_per_level <- 2 ^ (0:depth)
    # figure out how many nodes there are at each level of the tree
    effects_on_node <- empty_tree(depth)
    
    # go through all the node_effect_trees and read off the effects
    #
    for (i in 1:depth) {
      for (j in 1:nodes_per_level[i]) {
        t <- node_effect_trees[[i]][[j]]
        effect = t[[node_location[1]]][[1]][node_location[2]]
        effects_on_node[[i]][[1]][j] = effect
      }
    }
    return(effects_on_node)
  }



gaussian_scale <- function(epsilon, sensitivity, cdp = TRUE) {
  if (cdp) {
    # Adding N(0, sigma^2) noise satisfies GS/(2sigma^2)-zCDP
    # so N(0, (1/eps)^2) and GS=1 satisfies eps^2/2-zDP = rho-zCDP
    return(sensitivity / epsilon)
  } else {
    return(sensitivity * sqrt(2.0 * log(1.26 / delta)) / epsilon)
  }
}

variance_from_effects <- function(effects_tree, epsilon) {
  # calculate the variance of the noise added originally
  # variance = 2.0*(2.0/epsilon)**2 # why times 2.0?
  cdf_sensitivity <- 2.0
  variance <- (gaussian_scale(epsilon, cdf_sensitivity, cdp = TRUE)) ^ 2
  
  # square and sum all the values in the effects_tree
  total_effect <- 0
  for (i in 1:length(effects_tree)) {
    for (j in 1:length(effects_tree[i])) {
      total_effect <- total_effect + effects_tree[[i]][[1]][[j]] ^ 2
      net_variance = total_effect * variance
    }
  }
  
  return(net_variance)
}

## Pre-process Post-processed variances (not data-dependent)
getVariances <- function(tree, epsilon, cdp = TRUE) {
  tree1 = tree
  vals = tree1[[length(tree)]][[2]]
  variances = NULL
  i = 1
  
  # generate trees of node effects
  depth = length(tree1)
  node_effect_trees = node_effects(depth, epsilon)
  max_eps_per_val <-
    if (cdp)
      epsilon / sqrt(depth)
  else
    epsilon / depth
  
  
  while (i <= length(vals)) {
    mM <- tree1[[1]][[2]]
    # build tree to keep track of the effects of each node on the calculation
    effects_tree = empty_tree(depth)
    # iterate through layers of tree
    index = 0
    j = 1
    while (j <= depth) {
      # determine if should traverse tree to left or right
      mid = mM[1] + (mM[2] - mM[1]) / 2.0
      # if looking at leftmost node in the tree, we know empirical cdf should evaluate
      # to 0 not to bin size.
      if (i == 0) {
        break
      }
      # if you don't need higher granularity, stop traversal
      if (vals[i] == mM[2]) {
        effects = effects_on_node(depth, c(j, index + 1), node_effect_trees)
        effects_tree = sum_trees(effects_tree, effects)
        break
      }
      # if at leaves of tree, record the count there
      if (j == depth) {
        effects = effects_on_node(depth, c(j, index + 1), node_effect_trees)
        effects_tree = sum_trees(effects_tree, effects)
      } else if (vals[i] <= mid) {
        # if traversing left
        
        # reset max value to the mid, don't add to the count
        M = mid
        # set next index of node to look at in next layer
        index = index * 2
      } else {
        # if at end of tree, record count at that node
        #if traversing right,
        
        # reset min value to the mid
        m = mid
        # set to next index of node to look at in next layer
        index = index * 2 + 1
        effects = effects_on_node(depth, c(j + 1, index), node_effect_trees)
        effects_tree = sum_trees(effects_tree, effects)
      }
      j <- j + 1
    }
    ## Calculate variance of the count and std of the quantile
    variance = variance_from_effects(effects_tree, max_eps_per_val)
    variances <- c(variances, variance)
    
    i <- i + 1
  }
  return(variances)
}





a_binsearch <-
  function(n,
           qm_probs,
           noise_scale,
           a_lower,
           a_upper,
           alpha,
           epsilon) {
    a_gran = alpha / 10.0
    num_iters = as.integer(log2((a_upper - a_lower) / a_gran))
    # Very rough heuristic for scaling up iterations based on epsilon
    if (epsilon < 1.0) {
      num_iters = num_iters
    } else if (epsilon < 3.0) {
      num_iters = num_iters + 1
    } else if (epsilon < 5.0) {
      num_iters = num_iters + 2
    } else {
      num_iters = num_iters + 3
    }
    ret = a_upper
    # print("num_iters:", num_iters)
    for (i in 1:num_iters) {
      a_mid = (a_lower + a_upper) / 2.0
      terms = NULL
      for (j in 1:n) {
        qm = (j) / n
        qm_prob = qm_probs[j]
        normal_prob = 1.0 - pnorm(a_mid, mean = qm, sd = noise_scale)
        terms <- c(terms, qm_prob * normal_prob)
      }
      res_prob = sum(terms)
      # print("a_lower, a_mid, a_upper, res_prob:", a_lower, a_mid, a_upper, res_prob)
      if (res_prob < alpha) {
        a_upper = a_mid
        ret = a_mid
      } else{
        a_lower = a_mid
        ret = a_upper
      }
    }
    return(ret)
  }

n <- 100
lower_bound <- -4
upper_bound <- 4
granularity <- 0.4
epsilon <- 1
preProcessCDF <-
  function(n,
           lower_bound,
           upper_bound,
           granularity,
           epsilon,
           alpha,
           cdp = TRUE,
           variances = NULL) {
    if (is.null(variances)) {
      x = seq(lower_bound, upper_bound, length.out = n)
      n_bins = (upper_bound - lower_bound) %/% granularity + 1
      depth = as.integer(log2(n_bins) %/% 1)
      tt = dpTree(x, length(x), lower_bound, upper_bound, epsilon, depth, cdp =
                    cdp)
      tOpt = optimalPostProcess(tt, epsilon)[[1]]
      
      variances = getVariances(tOpt, epsilon, cdp = cdp)
    }
    # Compute all the binomial/normal approx coefficients
    qm_probs = dbinom(0:n, n, 0.5)
    
    # Compute the a_is
    a_lower = 0.0
    a_upper = 1.0
    a_list_lower = NULL
    a_list_upper = NULL
    for (i in 1:length(variances)) {
      noise_scale = sqrt(variances[i]) / (n) # variances are for counts, we want std for percents
      # print("noise_scale:", noise_scale)
      a_i_upper = a_binsearch(n,
                              qm_probs,
                              noise_scale,
                              a_lower,
                              a_upper,
                              alpha / 2.0,
                              epsilon)
      
      a_list_upper <- c(a_list_upper, a_i_upper)
      a_list_lower <- c(a_list_lower, 1.0 - a_i_upper)
    }
    cat("preprocessing finished!\n")
    return(list(a_list_lower, a_list_upper, variances))
  }

cdfFancyPostProcess <-
  function(vals, cdf, a_list_lower, a_list_upper) {
    # for lower interval, look at upper cis
    distances = cdf - a_list_lower
    indices = which(distances >= 0)
    if (length(indices) > 0 & min(indices) > 1) {
      i = min(indices) - 1
    } else {
      i = 1
    }
    lower_val = vals[i]
    
    # vice versa
    distances = cdf - a_list_upper
    indices = which(distances <= 0)
    if (length(indices) > 0 & max(indices) < length(vals) - 1) {
      i = max(indices) + 1
    } else{
      i = length(vals)
    }
    # print("upper measurement", cdf[i], "corresponding a", a_list_upper[i])
    upper_val = vals[i]
    
    return(c(lower_val, upper_val))
  }




# Sample from mixture of normals
rmixnorm <-
  function(n = 1000,
           weights = c(0.5, 0.5),
           component_means = c(-2, 2),
           component_sds = c(0.5, 0.5)) {
    # Sample component ids based on component weights
    component_id <-
      sample(
        1:length(weights),
        size = n,
        prob = weights,
        replace = TRUE
      )
    
    # Vectorized sampling from rnorm based on component id
    samples <-
      rnorm(n = n, mean = component_means[component_id], sd = component_sds[component_id])
    
    return(samples)
    
  }


# Get analytical cdf from mixture of normals
pmixnorm <- function(x,
                     weights = c(0.5, 0.5),
                     component_means = c(-2, 2),
                     component_sds = c(0.5, 0.5)) {
  q <- 0
  for (i in 1:length(weights)) {
    q <- q + weights[i] * pnorm(x, component_means[i], component_sds[i])
  }
  
  return(q)
  
}

# Get analytical density from mixture of normals
dmixnorm <- function(x,
                     weights = c(0.5, 0.5),
                     component_means = c(-2, 2),
                     component_sds = c(0.5, 0.5)) {
  q <- 0
  for (i in 1:length(weights)) {
    q <- q + weights[i] * dnorm(x, component_means[i], component_sds[i])
  }
  
  return(q)
  
}

# Make sure it's the same discretization as for the input
cdf_bootstrap <-
  function(cdf,
           B = 1000,
           n_data,
           lower_bound,
           upper_bound,
           epsilon,
           granularity,
           cdp,
           projection_step = TRUE) {
    probs <- c(cdf[[1]][1], diff(cdf[[1]]))
    
    
    boot_dist <- vector("list", B)
    
    for (i in 1:B) {
      samp_x <-   rmultinom(1, size = n_data, probs)
      x_boot <- NULL
      for (j in 1:length(probs)) {
        x_boot <- c(x_boot, rep(cdf[[2]][j], samp_x[j]))
        
      }
      boot_cdf <-
        dpCDF(x_boot,
              lower_bound,
              upper_bound,
              epsilon,
              granularity,
              cdp,
              num_trials = 1)
      
      if (projection_step) {
        boot_dist[[i]] <-
          isotone::gpava(boot_cdf[[1]][[2]], boot_cdf[[1]][[1]])$x
      } else {
        boot_dist[[i]] <- boot_cdf[[1]][[1]]
      }
      
    }
    
    return(boot_dist)
  }


# minimize F(X) >= q

get_quantile <- function(cdf,
                         q,
                         cdf_x = NULL,
                         cdf_y = NULL) {
  if (is.null(cdf_x) & is.null(cdf_y)) {
    cdf_x <- cdf[[2]]
    cdf_y <- cdf[[1]]
  }
  
  min(cdf_x[cdf_y >= q])
}

summarize_results <-
  function(res_list,
           true_value = 0,
           n = 1000,
           dataset = "Mixture of Normals",
           index = 1) {
    ci_coverage <-
      sapply(lapply(res_list, function(x)
        do.call("rbind", x[[index]])), function(x)
          mean(x[, 1] <= true_value & x[, 2] >= true_value))
    
    plot(
      ci_coverage,
      xaxt = "n",
      bty = "n",
      las = 1,
      ylab = "CI Coverage",
      xlab = "Epsilon",
      main = paste0("Coverage of CIs\n", dataset, " n = ", n),
      pch = 19,
      ylim = c(0.5, 1)
    )
    axis(1,
         at = 1:length(ci_coverage),
         labels = names(ci_coverage))
    abline(h = 0.95, lty = "dashed")
    legend(
      "bottomleft",
      legend = c("Nominal Coverage", "Empirical Coverage"),
      pch = c(NA, 19),
      lty = c("dashed", NA),
      bty = "n"
    )
    
    ci_len <-
      sapply(lapply(res_list, function(x)
        do.call("rbind", x[[index]])), function(x)
          mean(apply(x, 1, diff)))
    
    plot(
      ci_len,
      xaxt = "n",
      bty = "n",
      las = 1,
      ylab = "CI Length",
      xlab = "Epsilon",
      main = paste0("Length of CIs\n", dataset, "; N = ", n),
      pch = 19
    )
    axis(1, at = 1:length(ci_len), labels = names(ci_len))
    
  }

dp_ci <-
  function(data,
           B = 1000,
           quantile_of_interest = 0.5,
           alpha = 0.05,
           lower_bound = -4,
           upper_bound = 4,
           epsilon = 0.1,
           granularity = 0.01,
           cdp = TRUE,
           projection_step = TRUE) {
    n_data <- length(data)
    # Generate dp cdf
    private_cdf <-
      dpCDF(
        data,
        lower_bound,
        upper_bound,
        epsilon = epsilon,
        granularity = granularity,
        cdp = cdp,
        num_trials = 1
      )
    
    if (projection_step) {
      pava <- isotone::gpava(private_cdf[[1]][[2]], private_cdf[[1]][[1]])
      private_cdf[[1]][[1]] <- pava$x
    }
    
    boot_dist <-
      cdf_bootstrap(
        cdf = private_cdf[[1]],
        n_data = n_data,
        B = B,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        epsilon = epsilon,
        granularity = granularity,
        cdp = cdp,
        projection_step = projection_step
      )
    
    
    ci <- quantile(sapply(boot_dist, function(x)
      get_quantile(
        q = quantile_of_interest,
        cdf_x = private_cdf[[1]][[2]],
        cdf_y = x
      )),
      c(alpha / 2, 1 - alpha / 2))
    
    return(ci)
  }


draw_sample <- function(P, n = 1000) {
  if (is.function(P)) {
    P(n)
  }
  else if (is.numeric(P)) {
    sample(P, n)
  } else {
    stop(
      "P needs to be either a function (e.g. to sample from an analytical distribution) or a numeric vector"
    )
  }
}


## Exact CI for median from DescTools package SignTest.default
## See also ttp://www.stat.umn.edu/geyer/old03/5102/notes/rank.pdf,
## http://de.scribd.com/doc/75941305/Confidence-Interval-for-Median-Based-on-Sign-Test
cimed <- function(x, alpha=0.05, na.rm=FALSE) {
  if(na.rm) x <- x[! is.na(x)]
  n <- length(x)
  k <- qbinom(p=alpha / 2, size=n, prob=0.5, lower.tail=TRUE)
  ## Actual CL: 
  # 1 - 2 * pbinom(k - 1, size=n, prob=0.5)
  # 1 - 2 * pbinom(k - 1, size=n, prob=0.5) >= 1 - alpha
  sort(x)[c(k, n - k + 1)]
}
