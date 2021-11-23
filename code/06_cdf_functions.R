gaussian_scale<- function(epsilon, sensitivity, cdp=TRUE){
  if (cdp){
    # Adding N(0, sigma^2) noise satisfies GS/(2sigma^2)-zCDP
    # so N(0, (1/eps)^2) and GS=1 satisfies eps^2/2-zDP = rho-zCDP
    return(sensitivity / epsilon)
  } else {
    return(sensitivity * sqrt(2.0*log(1.26/delta)) / epsilon)
  }
}

### DP Functions ###
dpHistogram <- function(x, n, lower_bound, upper_bound, epsilon, bins){
  # true histogram
  breaks <- seq(lower_bound, upper_bound, length.out = bins+1)
  true_hist <- .Call(graphics:::C_BinCount, x, breaks, T, T)
  #true_hist = np.histogram(x, bins=bins, range=(lower_bound, upper_bound))
  # calculate laplace noise
  # noise = np.random.laplace(0, 2.0/epsilon, size=len(true_hist[0])) #histogram queries have sensitivity 2
  # caculate gaussian noise
  cdf_sensitivity <- 2.0
  noise <- rnorm(bins, 0, gaussian_scale(epsilon, cdf_sensitivity, cdp = TRUE))
  # noise = np.random.normal(0, common.gaussian_scale(epsilon, cdf_sensitivity, cdp=True), size=len(true_hist[0]))
  # add noise to bins
  #out = ([i+j for (i,j) in zip(true_hist[0], noise)], true_hist[1])
  out <- list(true_hist + noise, breaks)
  return(out)
}


dpTree <- function(x, n, lower_bound, upper_bound, epsilon, depth, cdp){
  # initialize list of number of bins at each level of tree
  #bins = [2**i for i in range(1,depth)]
  bins <- 2^(1:depth)
  # divide up epsilon
  if(cdp){
    eps <- epsilon/sqrt(depth)
  } else {
    eps <- epsilon/depth
  }
  # build noisy histogram at each level of tree
  #tree = [([n],(lower_bound, upper_bound))] + [dpHistogram(x, n, lower_bound, upper_bound, eps, bins[i]) for i in range(depth-1)]
  tree <- append(list(list(n, c(lower_bound, upper_bound))), lapply(bins, function(b) dpHistogram(x, n, lower_bound, upper_bound, eps, b)))
  
  
  return(tree)
}



inverseVariance <- function(epsilon){
  # histogram has sensitivity 2
  b <- epsilon/2
  return(2*(b^2))
}

# Create an array of elements that are adjacent e.g. adjacentElements([1,2,3,4]) returns [2,1,4,3]
adjacentElements <- function(ls){
  adj <- vector(mode = "numeric", length = length(ls)) 
  i <- 0
  while(i < length(ls)){
    if (i%%2 == 1){
      adj[i+1] <- ls[(i+1)-1]
    } else {
      adj[i+1] = ls[(i+1)+1]
    }
    i <- i + 1 
  }
  return(adj)
}

wBelow <- function(tree){
  #  """
  #    Recursive weight estimate from below. This assumes that the variance is the same for every node at in the tree.
  #    :param tree: Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.  
  #    :return: Single weight for each level of the tree.
  #    """
  weights <- vector(mode = "numeric", length = length(tree))  # initialize with one weight per level of the tree
  i <- length(tree)
  while(i >= 1){
    if(i == length(tree)){
      weights[i] = 1
    } else{
      prev = weights[i+1]
      weights[i] = (2.0*prev)/(2*prev+1)
    }#coerce to floats with 2.0
    i <- i - 1
  }
  return(weights)
}


countBelow <- function(tree, wBelows){
  #  """
  #    Recursively compute counts from below
  #    :param tree: Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
  #    :param wBelows: Array of weights of same length as tree, where the ith weight corresponds to the ith weight calculated from below.
  #    :return: List of counts for each node of the tree.
  #    """
  counts = vector(mode = "list", length = length(tree)) 
  i = length(tree) 
  while(i >= 1){
    if (i == length(tree) ) {
      counts[[i]] = tree[[i]][[1]]
    } else {
      w = wBelows[i]
      child = tree[[i+1]][[1]]
      #childSum = [sum(x) for x in zip(child[0:len(child)-1],child[1:len(child)])] # sum all pairs of children counts
      #childSum = childSum[0:len(childSum):2] # pick out the sums that are children with the same parents
      childSum <- (child[1:(length(child)-1)] + child[2:(length(child))])[seq(1,(length(child)-1),2)]
      weightedT <- w*tree[[i]][[1]] #weigh parent counts
      weightedChild = (1-w)*childSum #weigh child counts
      #counts[i] = [sum(x) for x in zip(weightedT, weightedChild)]
      counts[[i]] <- weightedT + weightedChild
    }#sum parent and child counts
    i  <- i - 1
  }
  return(counts)
}


wAbove <- function(tree, wBelows){
  #  """
  #    Recursive weight estimation from above
  #    :param tree: Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
  #    :param wBelows: Array of weights of same length as tree, where the ith weight corresponds to the ith weight calculated from below.
  #    :return: Single weight for each level of the tree.
  #    """
  weights = vector(mode = "numeric", length = length(tree)) 
  i = 1
  while(i <= length(tree)) {
    if(i == 1){
      weights[i] = 1
    } else {
      prevAbove = weights[i-1]
      prevBelow = wBelows[i]
      weights[i] = 1.0/(1.0 + (prevAbove + prevBelow)^(-1))
    }
    i <- i + 1 
  }
  return(weights)
}

countAbove <- function(tree, countsBelow, wAboves){
  #  """
  #    Recursively compute counts from above
  #    :param tree: Tree, formatted as a list of arrays where the contents of the ith array in the list is the ith level of the tree.
  #    :param countsBelow: Array of counts of same length as tree, as calculated by countBelow function
  #    :param wAboves: Weights computed from above, assuming each node in tree has same variance.
  #    :return: List of counts for each node of the tree.
  #    """
  counts = vector(mode = "list", length = length(tree)) 
  i = 1
  while(i <= length(tree)){
    if(i == 1){
      counts[[i]] = tree[[i]][[1]]
    } else {
      w = wAboves[i]
      #parents = [val for val in counts[i-1] for _ in (0, 1)] #replicating parent counts so dimensions match
      parents <- rep(counts[[i-1]], each = 2)
      adjacents <- adjacentElements(tree[[i]][[1]])
      parentAdjDiff <- parents - adjacents # get difference between parent and adjacent node
      weightedT <- w*tree[[i]][[1]]  #weight current node count
      weightedPA = (1-w)*parentAdjDiff #weighted parent - adjacent count
      counts[[i]] = weightedT + weightedPA
    }
    i <- i + 1
  }
  return(counts)
}


optimalCount <- function(tree, wA, cA, cB){
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
  while(i <= length(tree)){
    if(i == 1){
      counts[[i]] = tree[[i]][[1]]
    } else {
      w = wA[i]
      parents = rep(cA[[i-1]], each = 2) #replicating parent counts so dimensions match
      adjacents = adjacentElements(cB[[i]])
      parentAdjDiff = parents - adjacents # get difference between parent and adjacent node
      weightedT = w*cB[[i]] #weight current node count
      weightedPA = (1-w)*parentAdjDiff #weighted parent - adjacent count
      counts[[i]] = weightedT + weightedPA
    }
    i <- i + 1
  }
  
  for(i in seq_along(tree)) tree[[i]][[1]] <- counts[[i]]
  
  return(tree)
}

optimalPostProcess <- function(tree, epsilon){
  
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
postProcessedCDF <- function(tree, epsilon, monotonic=F){
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
  
  while(i <= length(vals)){
    # build tree to keep track of the effects of each node on the calculation
    # effects_tree = empty_tree(len(tree1))
    # start out with min and max values at top of tree
    mM = tree1[[1]][[2]]
    # initialize count for i
    count = 0
    # iterate through layers of tree
    index = 0
    j = 1
    while(j<=length(tree1)){
      # determine if should traverse tree to left or right
      mid = mM[1] + (mM[2]-mM[1])/2
      # if looking at leftmost node in the tree, we know empirical cdf should evaluate
      # to 0 not to bin size.
      if(i == 1){
        break
      }
      # if you don't need higher granularity, stop traversal
      if(vals[i] == mM[2]){
        count <- count + tree1[[j]][[1]][index+1]
        # effects = effects_on_node(depth, [j, index], node_effect_trees)
        # effects_tree = sum_trees(effects_tree, effects)
        break
      }
      # if at leaves of tree, record the count there
      if(j == length(tree1)){
        count <- count + tree1[[j]][[1]][index+1]
        # effects = effects_on_node(depth, [j, index], node_effect_trees)
        # effects_tree = sum_trees(effects_tree, effects)
        # if traversing left
      } else if (vals[i] <= mid){
        # reset max value to the mid, don't add to the count
        mM[2] = mid
        # set next index of node to look at in next layer
        index = index*2
      } else{
        # if at end of tree, record count at that node
        #if traversing right, 
        
        # reset min value to the mid
        mM[1] = mid
        # set to next index of node to look at in next layer
        index = index*2 + 1
        count <- count +  tree1[[j+1]][[1]][index] # add the node's left child to the count
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
  percents = counts/n
  
  percents[percents<0] <- 0
  percents[percents>1] <- 1
  #percents = counts/max(counts) # normalize counts by n
  #percents <- percents/max(percents)
  # variances.append(variance/n**2) # scaling by the normalization
  # variances.append(variance)
  
  return(percents)
}





dpCDF <- function(x, lower_bound, upper_bound, epsilon, granularity, cdp = TRUE, num_trials = 1){
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
  
  n_bins = (upper_bound-lower_bound) %/% granularity + 1
  depth = log2(n_bins)%/%1 
  
  results = vector("list", length = num_trials)
  for(i in 1:num_trials){
    t = dpTree(x, length(x), lower_bound, upper_bound, epsilon, depth, cdp)
    tOpt = optimalPostProcess(t, epsilon)[[1]]
    
    cdf = postProcessedCDF(tOpt, epsilon)
    
    results[[i]] <- list(cdf, t[[length(t)]][[2]])
  }
  return(results)
}

