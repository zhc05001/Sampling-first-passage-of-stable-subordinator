# Exact sampling of the first passage of a stable subordinator
This R package is for exact sampling of the first passage event of a stable subordinator across a boundary.  The algorithms implemented in the package are developed in 

- Chi, Z. (2024). *Complexity of exact sampling of the first passage of a stable subordinator*. [arXiv:xxxxx](http://merlot.stat.uconn.edu/~zhc05001/)

Currently the package only samples a certain transformation of the undershoot at the first pasage across constant level 1.  However, using the package, the entire event across any constant level or non-constant regular boundary can be easily sampled using the package; see Algorithm 2.1 in Chi (2024).  A short R function is provided below.

To sample `n` first passage events aross constant level 1 by a "standard" stable subordinator of index `alpha`, use the function `sample.fp` in the file `ar-fp.R` as follows
```R
source("ar-fp.R")
X=sample.fp(n, alpha)
```
`X` is a list of the following variables, each being an array of length `n`:`X$z`, `X$y`, `X$log.y`, and `X$theta`.  It is a little complicated to explain what these variables are.  However, `X$z` is a random parameter critical for the sample, `X$y` is a transformation of the undershoot, `X$log.y` is the logarithm of `X$y`, and `X$theta` is a value sampled otherwise with `X$y` from a bivariate distribution.

Now, suppose you want to sample the first passage event across a non-decreasing regular boundary.  By definition, the boundary is a function $`b(t)`$ of $`t\ge0`$ that is positive at least for $`t`$ small enough.  In addition to define $`b(t)`$, you also need to define $`b'(t)`$ and the inverse of $`B(t) = t^{-1/\alpha} b(t)`$:
```R
diff.b <- function(t) {
## Define the deriviative of b(t) here
...
}

inv.B <- function(t) {
## Define the inverse of B(t) here
...
}
```

