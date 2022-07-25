#### Function to compute probability of presence at larger spatial or hierarchical scale ####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

aggreg_prob = function(x, na.rm)
{ 
  # Case with all NA
  if (all(is.na(x)))
  { 
    y <- NA 
  } else {
    # Case with some NA but not all
    x <- x[!is.na(x)]
    
    # General case
    y <- 1-prod(1-x) # Probability of presence of species = probability of presence of at least one OMU = opposite of probability of absence of all OMU
  }
  
  # Output
  return(y) 
}