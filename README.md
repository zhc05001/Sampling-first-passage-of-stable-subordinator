# Exact sampling of the first passage of a stable subordinator
This R package is for exact sampling of the first passage event of a stable subordinator across a boundary.  The algorithms implemented in the package are developed in 

- Chi, Z. (2025). *Complexity of exact sampling of the first passage of a stable subordinator*. [arXiv:xxxxx](http://merlot.stat.uconn.edu/~zhc05001/)

Let $b(t)>0$ be a non-increasing differentable function of $t>0$.  To sample $n$ first passage events aross $`b(t)`$ by a "standard" stable subordinator $`S_t`$ of index $`0<\alpha<1`$, whose moment generating function is
$`E(e^{-\lambda S_t}) = \exp(-t\lambda^\alpha)`$, $`t>0`$, $`\lambda>0`$, proceed as follows.

First, supply the definitions of $`b(t)`$, its derivative $`b'(t)`$, and $`\ln [B^{-1}(s)]`$, where $`B(t)=t^{-1/\alpha} b(t)`$.  The general format is as follows
```R
    b <- function(t, alpha, par1, par2, ...) {
        ## t and alpha must be the first two input arguments, par1, par2, ... are additional parameters if needed
        ## there can be zero number of additional parameters, so that the only input arguments are t and alpha
        ## the definition should allow t to be an array of positive numbers
        ...
    }
    diff.b <- functon(t, alpha, par1, par2, ...) {
        ## derivatives of b
        ## same requirements on the input arguments as function b
        ...
    }
    linv.B <- function(s, alpha, par1, par2, ...) {
        ## logarithm of the inverse of B
        ## the definitin should allow s to be an array of numbers
        ...
    }
```
After the functions are defined, the first passage events can be sampled as follows, 
```
source("ar-fp.R")
X=sample.fp(n, alpha, b, diff.b, linv.B, par1, par2, ...)
```
`X` is a list of named items, each being an array of length `n`.  The ones directly related to the first passage event are
- **`t`:** time $`\tau`$ of the first passage event
- **`y`:** the value such that give $`\tau=t`$, $`b(t)(1+y)^{1-1/\alpha}`$ is the undershoot of the first passage, i.e., the value $`S_{t-}`$ of the subordinator just before the passage.  Given the undershoot, the jump of the subordinator at the first passage can be sampled by $`[b(t)-S_{t-}] V^{-1/\alpha}`$, where $`V`$ is uniformly distributed on $(0,1)$.  Note that for the above barrier, the subordinator may cross it by creeping, i.e., moving continously instead of jumping.  In the case of creeping, the value of $`y`$ is zero and $`S_{t-} = S_t = b(t)`$.
- **`log.y`:** the logarithm of `X$y`.  The function `sample.fp` actually samples $`\ln y`$, and then exponentiates the result to get $`y`$.  The reason for choosing this approach is that when $`\alpha`$ is close to 1, the sample value of $`y`$ is often extremely small, causing numerical issues.  In contrast, $`\ln y`$ is more numerically stable to sample.

`X` also has the following named items, which store sample valuse of intermediate random variables.
- **`z`** the value of a random variable dentoed $`z`$ in Chi (2025).  It is the first random variable to be sampled by `sample.fp`.  The time $`\tau`$ of the first passage is a deterministic function of $`z`$: $`\tau = B^{-1}(s)`$, where $`s = \alpha[(1-\alpha)/z]^{1/\alpha-1}`$.

### Example
Suppose $`b(t) = M - t^{1/\alpha}`$ if $`0\leq t\leq M^\alpha`$ and 0 otherwise, where $M>0$ is a parameter.  Then $`b'(t) = -(1/\alpha) t^{1/\alpha-1} I\{t<M^\alpha\}`$ and $B^{-1}(s) = [M/(s+1)]^\alpha$.  To allow vectorized computation, the following R code can be used.  Note that the definition of $`\log B^{-1}(s)`$ instead of $`B^{-1}(s)`$ has to be supplied.
```R
    b1 <- function(t,alpha,M) {  ## 1st barrier function to be tested.
        x=rep(0,length(t))
        I = which(t<M^alpha)
        x[I]=M-t[I]^(1/alpha)
        x
    }
    diff.b1 <- function(t,alpha,M) {
        x=rep(0,length(t))
        I = which(t<M^alpha)
        x[I]=-(1/alpha)*t[I]^(1/alpha-1)
        x
    }
    linv.B1 <- function(s,alpha,M) {
        alpha*(log(M) - log1p(s))   ## log1p(s) = log(1+s)
    } 
```
After the functions are defined, the first passage events can be done as follows, where $`n=1000`$ samples are drawn, the index of the subordinator is $`\alpha=0.9`$, and $`M=100`$ in the barrier function,
```
source("ar-fp.R")
X=sample.fp(1000, 0.9, b1, diff.b1, linv.B1, 100)
```
### Citation
If you find this code useful, please cite it using the following BibTeX entry:
```bibtex
@software{chi:fp:2025,
    author = {Chi, Zhiyi},
     month = June,
     title = {Exact sampling of the first passage of a stable subordinator},
       url = {https://github.com/zhc05001/Exact-sampling-of-the-first-passage-of-a-stable-subordinator},
   version = {1.0.0},
      year = 2025,
}
```

