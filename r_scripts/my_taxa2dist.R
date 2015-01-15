# WTD functions

test_fun <- function(x,y) {
#  if (x == "unclassified") {x = NA}
#  if (y == "unclassified") {y = NA}
  if (xor(x == "unclassified", y == "unclassified")) {
    # only one is missing
    return(1)
  } else if (x != y) {
    return(2)
  } else {
    return(0)
  }
}

test_fun_vec <- Vectorize(test_fun) # vectorize test_fun

# slightly modified implementation of the taxa2dist function to handle partial taxonomies and the
# weighted taxonomic distance (WTD)
my_taxa2dist <- function (x, varstep = FALSE, check = FALSE, wtd= FALSE, labels) 
{
  # count the number of unique members in Genus Family, Order, etc.
  rich <- apply(x, 2, function(taxa) length(unique(taxa)))
  S <- nrow(x) # number of unique taxa/species
  
  # removes classes that have only one taxa or classes that are all constant
  if (check) {
    keep <- rich < S & rich > 1
    rich <- rich[keep]
    x <- x[, keep]
  }
  # orders columns by their number of unique entries (taxa)
  # TODO: Figure out when this is nessisary.
  # i <- rev(order(rich))
  # x <- x[, i]
  # rich <- rich[i] # do the same for richness
  
  if (varstep) {
    # vary step lengths relative to their proportinal loss of the number of distinct classes
    add <- -diff(c(nrow(x), rich, 1)) # create a set of iterated differences
    add <- add/c(S, rich) # scale results based on number of species in each class
    add <- add/sum(add) * 100 # scale between 0 and 100
  }
  else if (wtd) {
    add <- 1/(2^((ncol(x)-1):-1))
    add <- add *100 / sum(add)
    # add <- rep(1.0, length(add))
    # add <- add/sum(wtd_add) * 100 # scale between 0 and 100
  }
  else {
    # uniform scaling
    add <- rep(100/(ncol(x) + check), ncol(x) + check)
  }

  if (!is.null(names(add))) 
    names(add) <- c("Base", names(add)[-length(add)])
  if (!check) 
    add <- c(0, add)
  
  # calculate distances
  # out <- matrix(add[1], nrow(x), nrow(x))
  out <- matrix(0, nrow(x), nrow(x))
  for (i in 1:ncol(x)) {
    # adds the distance to the appropreate cells of the matrix
    out <- out + add[i+1] * outer(x[, i], x[, i], test_fun_vec)
  }
  
  out <- as.dist(out) # transforms the matrix to a distance object
  attr(out, "method") <- "taxa2dist" # add some annotation
  attr(out, "steps") <- add # add the step-vector
  
  # attach the labels
  if (missing(labels)) {
    attr(out, "Labels") <- rownames(x)
  } else {
    if (length(labels) != nrow(x)) 
      warning("Labels are wrong: needed ", nrow(x), " got ", 
              length(labels))
    attr(out, "Labels") <- as.character(labels)
  }
  if (!check && any(out <= 0)) 
    warning("you used 'check=FALSE' and some distances are zero -- was this intended?")
  out
}