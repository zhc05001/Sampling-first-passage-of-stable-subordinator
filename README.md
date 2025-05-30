# Exact sampling of the first passage of a stable subordinator
This R package is for exact sampling of the first passage event of a stable subordinator across a boundary.  The algorithms implemented in the package are developed in 

- Chi, Z. (2024). *Complexity of exact sampling of the first passage of a stable subordinator*. [arXiv:xxxxx](http://merlot.stat.uconn.edu/~zhc05001/)

Let $b(t)>0$ be a non-increasing differentable function of $t>0$.  To sample $n$ first passage events aross $`b(t)`$ a "standard" stable subordinator of index $`\alpha`$, use the function `sample.fp` in the file `ar-fp.R` as follows.

First, supply the definitions of $b(t)$, its derivative $b'(t)$, and $\log [B^{-1}(s)]$, where $B(t)=t^{-1/\alpha} b(t)$.  For example, suppose $b(t) = (M - t^{1/\alpha})_+$, where $M>0$ is a parameter that you want to be able to adjust, then $b'(t) = -(1/\alpha) t^{1/\alpha-1} I{t<M^\alpha}$ and $B^{-1}(s) = [M/(s+1)]^\alpha$.  To allow vectorized computation, the following R code can be used.  Note that the definition of $`\log B^{-1}(s)`$ instead of $`B^{-1}(s)`$ has to be supplied.
```R
    b <- function(t,a,M) {  ## t can be an array of positive numbers
        x=rep(0,length(t))
        I = which(t<M^a)
        x[I]=M-t[I]^(1/a)
        x
    }
    diff.b <- function(t,a,M) {
        x=rep(0,length(t))
        I = which(t<M^a)
        x[I]=-(1/a)*t[I]^(1/a-1)
        x
    }
    linv.B <- function(s,a,M) { ## s can be an array of positive numbers
        a*(log(M) - log1p(s))   ## log1p(s) = log(1+s)
    } 
```
After the functions are defined, the first passage events can be done as follows
```
source("ar-fp.R")
X=sample.fp(n, alpha, b, diff.b, linv.B, M)
```
`X` is a list of the following variables, each being an array of length `n`.  The ones directly related to the first passage event are
- **`X$t`**, `X$log.y`, and `X$y`.
