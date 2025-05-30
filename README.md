# Exact sampling of the first passage of a stable subordinator
This R package is for exact sampling of the first passage event of a stable subordinator across a boundary.  The algorithms implemented in the package are developed in 

- Chi, Z. (2024). *Complexity of exact sampling of the first passage of a stable subordinator*. [arXiv:xxxxx](http://merlot.stat.uconn.edu/~zhc05001/)

Let $`b(t)>0`$ be a non-increasing differentable function of $`t>0`$.  To sample $`n`$ first passage events aross $`b(t)`$ a "standard" stable subordinator of index $`\alpha`$, use the function `sample.fp` in the file `ar-fp.R` as follows.

First, supply the definitions of $`b(t)`$, its derivative $`b'(t)`$, and $`\log [B^{-1}(s)]`$, where $`B(t)=t^{-1/\alpha} b(t)`$.  For example, suppose $`b(t) = (M - t^{1/\alpha})_+`$, where $`M>0`$ is a parameter that you want to be able to adjust, then $`b'(t) = -(1/\alpha) t^{1/\alpha-1} I\{t<M^\alpha\}`$ and $`B^{-1}(s) = [M/(s+1)]^\alpha`$.  
```R
    b <- function(t,a,M) {
        x=rep(0,length(t));
        I = which(t<M^a); x[I]=M-t[I]^(1/a);
        x
    }
    diff.b <- function(t,a,M) {
        x=rep(0,length(t));
        I = which(t<M^a); x[I]=-(1/a)*t[I]^(1/a-1);
        x
    }
    linv.B <- function(s,a,M) { a*(log(M) - log1p(s)) }  ## log1p(s) = log(1+s).
```
aaa
```


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

