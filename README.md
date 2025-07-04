The R code in this project is for exact sampling of the first passage event of a stable subordinator across a boundary.  The algorithms implemented are developed in 

- Chi, Z. (2025). [*Complexity of exact sampling of the first passage of a stable subordinator*](https://arxiv.org/abs/2506.03047/).

Let $b(t)>0$ be a non-increasing differentable function of $t>0$.  To sample $n$ first passage events aross $`b(t)`$ by a "standard" stable subordinator $`S_t`$ of index $`0<\alpha<1`$, whose moment generating function is
$`E(e^{-\lambda S_t}) = \exp(-t\lambda^\alpha)`$, $`t>0`$, $`\lambda>0`$, proceed as follows.

First, supply the definitions of $`b(t)`$, its derivative $`b'(t)`$, and $`\ln [B^{-1}(s)]`$, where $`B(t)=t^{-1/\alpha} b(t)`$.  The general format is as follows
```R
    b <- function(t, alpha, par1, par2, ...) {
        ## t and alpha must be the first two input arguments, par1, par2, ... are additional
        ## parameters
        ## The function has to allow t to be an array of positive numbers
        ...
    }
    diff.b <- functon(t, alpha, par1, par2, ...) {
        ## derivatives of b
        ## same requirements on the input arguments as function b
        ...
    }
    linv.B <- function(s, alpha, par1, par2, ...) {
        ## logarithm of the inverse of B
        ## The function has to allow s to be an array of numbers
        ...
    }
```
After the functions are defined, the first passage events can be sampled as follows, 
```
source("ar-fp.R")
X=sample.fp(n, alpha, b, diff.b, linv.B, par1, par2, ...)
```
`X` is a list of named items, each being an array of length `n`.  The ones directly related to the first passage event are
- **`t`:** the time of the first passage event
- **`y`:** the value such that give the first passage occurs at time $`t`$, $`b(t)(1+y)^{1-1/\alpha}`$ is the undershoot of the first passage, i.e., the value $`S_{t-}`$ of the subordinator just before the passage.  Given the undershoot, the jump of the subordinator at the first passage can be sampled by $`[b(t)-S_{t-}] V^{-1/\alpha}`$, where $`V`$ is uniformly distributed on $(0,1)$.  Note that if $`b'(t)<0`$, the subordinator may cross the boundary by creeping, i.e., moving continously instead of jumping, the probability of which is $`-b'(t)/[-b'(t) + b(t)/(\alpha t)].`$  In the case of creeping, the value of $`y`$ is zero and $`S_{t-} = S_t = b(t)`$.
- **`log.y`:** the logarithm of `y`.  The function `sample.fp` actually samples $`\ln y`$, and then exponentiates the result to get $`y`$.  The reason for choosing this approach is that when $`\alpha`$ is close to 1, the sample value of $`y`$ is often extremely small, causing numerical issues.  In contrast, $`\ln y`$ is more numerically stable to sample.

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
After the functions are defined, to sample $`n=1000`$ first passage events of a subordinator with index $`\alpha=0.9`$ across barrier $`b(t) = \max(100 - t^{1/0.9},0)`$,
```
source("ar-fp.R")
X=sample.fp(1000, 0.9, b1, diff.b1, linv.B1, 100)
```

`X` also has the following named items, which store sample valuse of intermediate random variables.
- **`z`:** the value of a random variable denoted $`z`$.  Given the value of $`z`$, the time of the first passage is equal to $`B^{-1}(s)`$, where $`s = \alpha[(1-\alpha)/z]^{1/\alpha-1}`$.
- **`log.z`:** the logarithm of $`z`$
- **`theta`:** Given the value of $`z`$, instead of being sampled along, $`y`$ is jointly sampled with another random variable $`\theta`$ from a bivariate distribution parameterized by $`\alpha`$ and $`z`$.  In the case of creeping, $`\theta`$ is not defined.  The hard part is how to sample $`(y,\theta)`$ given that the subordinator jumps across $`b(t)`$.  The bulk of the package is for this task.  The package contains three functions to sample from the bivariate p.d.f. for different regions of $`(\alpha,z)`$.  They are `sample.chi.alpha.z.3.1`, `sample.chi.alpha.z.3.2`, and `sample.chi.alpha.z.5.1`, which are named after the algorithms in Chi (2025) to sample from a bivariate p.d.f. which is in proportin to a certain function $`\chi_{\alpha,z}(y,\theta)`$.

### Citation
If you find this code useful, please cite it using the following BibTeX entry:
```bibtex
@software{chi:25:github,
    author = {Chi, Zhiyi},
     month = June,
     title = {Sample first passage of stable subordinator},
       url = {https://github.com/zhc05001/Sampling-first-passage-of-stable-subordinator},
   version = {1.0.0},
      year = 2025,
}
```

