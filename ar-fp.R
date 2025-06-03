## $Id: ar-fp.R,v 1.1 2025/06/03 15:30:19 zchi Exp zchi $

library("stabledist");
library("pracma"); ### Must be loaded before "rwfec" so that its sinc() is
                   ### masked by the same named function in rwfec
library('rwfec')   ### For sinc()
library('crayon');
log1p.recip <- function(t)
### Calcultate log(1+1/t).
### When t is extremely small, log1p(1/t) does not work well.
### For example, if t = 10^-320, then numerically t is nonzero.  However
### in R, log1p(1/t)=Inf, while using log1p(1/t) = log1p(t) - log(t),
### the value is 736.83.
###
### Input t can be an array.
{
    eps = 10^-6;
    val = rep(NaN, length(t));
    I = which(abs(t)>1 | (t<=1 & t>=eps));
    val[I] = log1p(1/t[I]);
    I = which(t<eps&t>=0);
    val[I] = log1p(t[I]) -log(t[I]);
    return(val);
} ## log1p.recip

sample.mat <- function(w.mat)
### The input argument w.mat is an array or matrix.  If an array, it is
### treated as a one-column matrix.  All the entries of w.mat must be
### nonnegative and each of its column sums must be positive.  The number
### of outputs is the same as the number of columns, such that, if there are
### N rowns in w.mat, the i-th output is sampled from 1 to N according to
### the weights in the i-th column in w.mat
{
    if (length(w.mat)==0) return(numeric(0));
    if (!is.matrix(w.mat)) {
        w.mat = matrix(w.mat, nc=1);
    }

    N = nrow(w.mat);
    probs = apply(w.mat, 2, cumsum);
    s = rowSums(runif(ncol(w.mat))*probs[N,]>t(probs))+1;
    return(s);
} ## sample.mat

log.sinc <- function(x, difforder=0)
### Logarithm of sinc or the derivative of log of sinc.  The argument x
### can be an array.
###
### Formulas used: 25.8.8, DLMF NIST for the expansion of log(sinc(x));
### section 25.6(i) for the expression of zeta(2*n);
### section 24.1 and 24.6 for B_{2n} (Bernoulli numbers)
{
    eps = 10^-6;
    if (length(x)==0) return(numeric(0));
    if (sum(difforder == 0:2)==0) {
        stop("Derivative can be only calculated up to order 2");
    }

    val = x;
    I = which(abs(x)>=eps);
    J = which(abs(x)<eps);
    x2=x[J]^2; x4=x2*x2; x6=x2*x4;
    
    if (difforder==0) {
        val[I] = log(sinc(x[I]))
        val[J] = -x2/6*(1 + x2/30 + 2*x4/945 + x6/8100)
    } else if (difforder==1) {
        val[I] = cos(x[I])/sin(x[I]) - 1/x[I];
        val[J] = -x[J]/3*(1 + x2/15 + 2*x4/315 + x6/2025);
    } else {
        val[I] = 1/x[I]^2 - 1/(sin(x[I]))^2;
        val[J]= -(1+ x2/5 + 2*x4/63 + 7*x6/2025)/3
    } 
    return(val);
} ## log.sinc

x.over.expm1  <- function(x, series=FALSE)
### The value of x/(exp(x)-1).  The input argument x can be an array
{
    eps = 10^-4;
    if (length(x)==0) return(numeric(0));
    val = x;
    I = which(abs(x)>eps)
    J = which(abs(x)<=eps);

    val[I] = x[I]/expm1(x[I]);
    
    z2=x[J]^2; z4 = z2*z2; z6 = z2*z4; z8 = z4*z4;
    val[J] = 1 - x[J]/2 + z2/12 - z4/720 + z6/30240 - z8/1209600
    return(val);
} ## x.over.expm1

log.expm1 <- function(x)
### The value of log(exp(x)-1).  The input argument x can be an array
{
    eps = 10^-4; M=10^2;
    if (length(x)==0) return(numeric(0));
    val = rep(NaN, length(x));

    I = which(x>eps & x <=M);
    val[I] = log(expm1(x[I]));
    
    I = which(x<=eps & x>=0);
    val[I] = log(x[I])-log(x.over.expm1(x[I]));

    I = which(x>M);
    val[I] = x[I] + log(-expm1(-x[I]));

    return(val);
} ## log.expm1

log.rgamma.small.shape  <- function(n, shape, scale=1)
### A modified version of rgamss.R.  It only returns log of the sampled
### values.
###
### Goal: Sample small-shape gamma random variables via accept-reject           
### Paper: arXiv:1302.1884                                                      
### Authors: C. Liu, R. Martin (www.math.uic.edu/~rgmartin), and N. Syring
{
    if (shape > 0.2) return(log(rgamma(n, shape, scale=scale)));

    lambda = 1/shape - 1;
    w = shape/exp(1)/(1-shape);
    
    eta <- function(x) {
        y = rep(0, length(x));
        I = which(x>=0);
        y[I] = exp(-x[I]);
        
        I = which(x<0);
        y[I] = w*lambda*exp(lambda*x[I]);
        
        return(y);
    }
    
    h <- function(x) {
        exp(-x - exp(-x/shape));
    }
    
    z = numeric(n);
    I = 1:n;
    while(length(I)>0) {
        z[I] = rexp(length(I));
        
        k = sample(c(0,1), length(I), replace=TRUE, prob=c(1,w));
        J = which(k!=0);
        z[I[J]] = -z[I[J]]/lambda;
        
        I = I[which(eta(z[I])*runif(length(I))>= h(z[I]))];
    }
    
    oo = log(scale) - z/shape;
    return(oo);
} ## log.rgamma.small.shape

Xi.alpha <- function(theta, alpha)
### Calculate $\Xi_\alpha(\theta)$ in Prop. 2.2.
### Input argument theta is an array, and alpha must be a scalar in (0,1).

{
    if (length(theta)==0) return(numeric(0));
    if (length(alpha)!=1 | alpha<=0 | alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in (0,1)");
    }
    delta = 1-alpha;
    lsx0 = log.sinc(delta*theta);
    lsx1 = log.sinc(pi - alpha*theta);
    lsx=log.sinc(pi - theta);
    return(exp(lsx0 - lsx + (lsx1-lsx)*alpha/delta));
} ## Xi.alpha

K.alpha  <- function(theta, alpha)
### Calculate $K_\alpha(\theta)$ in Eq. (2.11).  Same requirement on
### the dimensions of the input arguments as Xi.alpha(theta, alpha)
{
    return(Xi.alpha(theta,alpha)*theta/(pi/alpha-theta));
} ## K.alpha

log.H.alpha <- function(theta, alpha)
### Logarithm of $H_\alpha(\theta)$ in Eq. (2.10).  The input argument
### theta is an array of numbers in [0, pi] and alpha must be a scalar in
### [0,1]
{
    if (length(theta)==0) return(numeric(0));
    if (any(theta<0 | theta>pi)) {
        stop(paste("All entries of", bold(red("theta")), "must be in [0,pi]"));
    }
    if (length(alpha)!=1 | alpha<0 | alpha>1 ) {
        stop(bold(red("alpha")), " must be a scalar in [0,1]");
    }

    delta = 1-alpha;
    thresh=10^-4;            

    val = rep(0, length(theta));
    if (alpha>0) {
        val[which(theta==pi)] = Inf;
    }
    if (alpha==1) {
        I = which(theta>0 & theta < pi);
        sx = sinc(theta[I]);
        lsx = log.sinc(theta[I]);
        val[I] = 1 - cos(theta[I])/sx - lsx;
    } else if (alpha>0) {
        I = which(theta > 0 & theta < (1-thresh)*pi);
        x = sinc(delta*theta[I]) / sinc(theta[I]);
        y = sinc(alpha*theta[I]) / sinc(theta[I]);
        ## val[I] = log1p(x-1) + alpha/delta*log1p(y-1);
        val[I] = log(x) + alpha/delta*log(y);
        
        I = which(theta >= (1-thresh)*pi & theta < pi);
        x = K.alpha(theta[I], alpha);
        y = delta*pi / alpha/(pi - theta[I]);
        val[I] = log(x) + 1/delta*log1p(y);
    }
    return(val);
} ## log.H.alpha

dlog.Xi.alpha <- function(theta, alpha)
### Derivative of $\log [\Xi_\alpha(\theta)]$.  Same requirements on the
### input arguments as Xi.alpha(theta, alpha), except that the input
### argument alpha must be a scalar in (0,1).
{
    if (length(theta)==0) return(numeric(0));
    if (length(alpha)!=1 | alpha<=0 | alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in (0,1)");
    }

    delta = 1 - alpha;
    d = delta * log.sinc(delta*theta, difforder=1) -
        (alpha^2 * log.sinc(pi - alpha*theta, difforder=1) -
         log.sinc(pi-theta, difforder=1))/delta;
    return(d);
} ## dlog.Xi.alpha

dlog.H.alpha <- function(theta, alpha)
### Derivative of $\log[H_\alpha(\theta)]$.  Same requirements on the input
### arguments as Xi.alpha(theta, alpha).  When theta is close to 0 or
### pi, the derivative is evaluated using the last formula in Prop. 2.4.
{
    if (length(theta)==0) return(numeric(0));
    if (any(theta<0 | theta>pi)) {
        stop("All entries of ", bold(red("theta")), " must be in [0,pi]");
    }
    if (length(alpha)!=1 | alpha<0 | alpha>1) {
        stop(bold(red("alpha")), " must be a scalor in [0,1]");
    }
    
    ## From Eq. (2.8) and the fact that ln H_alpha/alpha is symmetric about 1/2
    if (alpha==1) {
        val = theta + (cos(theta)/sinc(theta)-1) *
            log.sinc(theta,difforder=1);
        return(val);
    } else if (alpha==0) {
        val = rep(0, length(theta));
        return(val);
    }

    delta = 1-alpha;
    thresh = 10^-4;

    val = rep(0, length(theta));
    val[which(theta==pi)] = Inf;
    
    I = which(theta > 0 & theta < (1-thresh)*pi);
    val[I] = delta * log.sinc(delta * theta[I], difforder=1) +
        (alpha^2 * log.sinc(alpha * theta[I], difforder = 1) -
         log.sinc(theta[I], difforder = 1)) / delta;
    
    I = which(theta >= (1-thresh)*pi & theta < pi);
    d = pi^2 / ((pi-theta[I]) * theta[I] * (pi - alpha*theta[I]));
    val[I] = dlog.Xi.alpha(theta[I], alpha) + d;

    return(val);
} ## dlog.H.alpha

sample.z <- function(n, alpha, log=FALSE)
### Sample from the distribution in Eq. (2.7).  The input argument n is the
### size of the sample, and alpha must be a scalar in [0,1]
{
    if (n<=0) return(numeric(0));
    if (length(alpha)!=1 | alpha<0 | alpha>1) {
        stop(bold(red("alpha")), " must be a scalar in [0,1]");
    }
    if (!log) {
        return(rexp(n)/exp(log.H.alpha(pi*runif(n), alpha)));
    } else {
        return(log(rexp(n)) - log.H.alpha(pi*runif(n), alpha));
    }
} ## sample.z

c.alpha  <- function(alpha)
### $c_\alpha$ in Eq. (4.1)
{
    delta=1-alpha;
    (alpha/delta)^alpha
} ## c.alpha

sample.chi.alpha.z.3.1  <- function(z, alpha)
### Algorithm 3.1 to sample from the normalized $\chi_{\alpha,z}$.  The
### input argument z must be an array of positive numbers, and alpha must
### be a scaler in (0,1).  Note: the logarithms of the sample values from
### the distribution are also returned 
{
    cat("Entering sample.chi.alpha.z.3.1...\n");

    nz=length(z);
    if (nz==0) return(numeric(0));
    if (length(alpha)!=1 | alpha<=0 | alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in (0,1)");
    }
    if (any(z<=0)) {
        stop("All the entries of ", bold(red("z")), " must be positive");
    }
    
    sample.m <- function(s) {
        ## Sample from the normalized $m_s(x) = e^{-sx^2/2}I\{0<x<\pi\}$
        ## defined in the display below Eq. (3.4); s is an array
        theta = rep(0, length(s)); r = sqrt(s);

        I = which(r<=1/2);
        while (length(I)>0) {
            theta[I] = runif(length(I))*pi;
            I = I[which(runif(length(I)) > exp(-s[I]*theta[I]^2/2))];
            ## I is updated to the *rejected* entries
        }

        I = which(r>1/2);
        while (length(I)>0) {
            theta[I] = abs(rnorm(length(I)))/r[I];
            I = I[which(theta[I]>pi)];
        }
        return(theta);
    }

    delta = 1-alpha;  daratio = delta/alpha;
    az = alpha*z;
    c.z = (z+delta)*z^(-delta)*pmax(1 + alpha * pi^2/2, 1/z)*exp(-z);
    ## c.z is $r_{\alpha,z} * e^{-z}$, where $r_{\alpha,z}$ is defined in
    ## the display below Eq. (3.4).
    
    log.z = log(z);
    y = rep(Inf, nz); log.y = y; theta = y; tau = y;
    
    I = 1:nz;  ## Indices of entries not yet sampled
    while (length(I)>0) {
        J = I;
        while (length(J)>0) { ## Inner loop of Algorithm 3.1
            theta[J] = sample.m(az[J]);
            log.tau = log.z[J] + log.H.alpha(theta[J], alpha);
            tau[J] = exp(log.tau);

            ## Keep entries that are *rejected*
            J = J[which(tau[J]==Inf | runif(length(J)) * c.z[J] >
                        exp(az[J] * theta[J]^2/2 + alpha*log.tau - tau[J] +
                            log1p(delta/tau[J])))];
        } ## End innter loop

        log.y[I] = log.rgamma.small.shape(length(I), delta, 1/tau[I]);
        y[I] = exp(log.y[I]);
        
        sub=I[which(runif(length(I))<=delta/(tau[I] + delta))];
        if (length(sub)>0) {
            y[sub] = y[sub] + rexp(length(sub))/tau[sub];
            log.y[sub] = log(y[sub]);
        }

        R = x.over.expm1(-daratio*log1p(y[I]))/x.over.expm1(log1p(y[I]));
        
        ## Keep entries that are *rejected*
        I= I[which((y[I] + 1) * runif(length(I)) >= R^alpha)];
    }
    one.minus.x = -expm1(-daratio*log1p(y));  ## one minus undershoot

    return(list(log.y=log.y, y=y, theta=theta, undershoot.gap=one.minus.x,
                undershoot=1-one.minus.x, alpha=alpha));
} ## sample.chi.alpha.z.3.1

sample.chi.alpha.z.3.2 <- function(z, alpha)
### Algorithm 3.2 to sample from the normalized $\chi_{\alpha,z}$.  The
### input argument z must be an array of positive numbers, and alpha
### must be a scaler in (0,1).
{
    cat("Entering sample.chi.alpha.z.3.2...\n");
    nz = length(z);
    if (nz==0) return(numeric(0));
    if (length(alpha)!=1 | alpha<=0 | alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in (0,1)");
    }
    if (any(z<=0)) {
        stop("All the entries of ", bold(red("z")), " must be positive");
    }

    delta = 1-alpha
    g.delta = gamma(delta);
    daratio = delta/alpha;

    ## cc is c_{2,alpha}/c_alpha on line 6 of Algorithm 4.2
    cc = 1/c.alpha(min(delta, alpha));

    log.z = log(z);
    y = rep(Inf, nz); log.y = y; theta = y; tau = y;
    one.minus.x = y;  ## one minus undershoot
    
    I = 1:nz; ## Indices of entries not yet sampled

    while (length(I)>0) {
        J = I;

        while (length(J)>0) {    ## Inner loop of Algorithm 3.2
            theta[J] = runif(length(J))*pi;
            tau[J] = exp(log.z[J] + log.H.alpha(theta[J], alpha));

            
            acp = rep(0, length(J));
            .aj = which(tau[J]<Inf);
            if (length(.aj)>0) {
                acp[.aj] = exp(alpha*log(tau[J[.aj]]) + log(g.delta + tau[J[.aj]]^(-alpha)) -tau[J[.aj]])/(g.delta+1);
            }
            
            ## Keep entries that are *rejected*
            J=J[runif(length(J))>=acp];
        }
        i = (runif(length(I)) > 1/(g.delta*tau[I]^alpha + 1));

        sub = I[which(i)];
        if (length(sub)>0) {
            log.y[sub] = log.rgamma.small.shape(length(sub), delta, 1/tau[sub]);
            y[sub] = exp(log.y[sub]);
        }
        sub = I[which(!i)];
        if (length(sub)>0) {
            y[sub] = rexp(length(sub))/tau[sub];
            log.y[sub] = log(y[sub]);
        }
        one.minus.x[I] = -expm1(-daratio*log1p(y[I]));

        R = x.over.expm1(-daratio*log1p(y[I]))/x.over.expm1(log1p(y[I]));
        I = I[which(runif(length(I))*cc*(1+y[I]^alpha)>R^alpha)];
    }

    return(list(y=y, log.y=log.y, theta=theta, undershoot.gap=one.minus.x,
                undershoot=1-one.minus.x, alpha=alpha));
} ## sample.chi.alpha.z.3.2

truncated.exp <- function(x, t)
### The value of the pdf $(t/2)*\min(1, e^{1-tx})I\{x>0\}$ in Eq. (4.1).
### The input arguments x and t can be arrays.  If the two arrays have
### different length, the shorter one internally are padded with its
### repeats to have the same length as the longer one.  Each entry of t
### must be nonnegative.
{
    if (length(x)==0 | length(t)==0) return(numeric(0));
    if (any(t<0)) {
        stop("All entries of ", bold(red("t")), " must be nonnegative");
    }
    max.l = max(length(x), length(t));
    x = x + rep(0, max.l);
    t = t + rep(0, max.l);
    
    val = rep(0, max.l);
    sub = which(x>=0);
    val[sub] = (t[sub]/2)*pmin(1, exp(1-t[sub]*x[sub]));

    return(val);
} ## truncated.exp

sample.truncated.exp  <- function(t)
### Sample the pdf $(t/2)*min(1, e^{1-t*x})I\{x>0\}$ in Eq. (4.1).  The
### input t is an array of positive numbers.  The output has the same
### length as t, such that its i-th entry is sampled is from the truncated
### exponential with parameter value t[i].
{
    if (length(t)==0) return(numeric(0));
    if (any(t<=0)) {
        stop("All entries of ", bold(red("t"))," must be positive");
    }
    y = runif(length(t))*2;
    sub = which(y>1);
    y[sub] = 1-log1p(1-y[sub]);
    return(y/t);
} ## sample.truncated.exp

M.abc <- function(a,b,c)
### Evaluate $M(a,b,c)$ in Eq. (4.6).  The input arguments a and b are
### arrays whose entries must be nonnegative and c must be a scalar in
### (0,1).  If a and b have different lengths, the shorter one internally
### are padded with its repeats to have the same length as the longer one.
### Because for large $b$, $M(a,b,c)$ numerically can be Inf, the log
### of $M(a,b,c)$ is also output.
{
    if (length(a)==0 | length(b)==0) return(numeric(0));
    if (length(c)!=1 | c<0 | c>=1) {
        stop(bold(red("c"))," must be a scalar in (0,1)");
    }
    if (any(c(a,b)<0)) {
        stop("All entries of ", bold(red("a"))," and ", bold(red("b")), " must be nonnegative");
    }
    max.l = max(length(a), length(b));
    a = a + rep(0, max.l);
    b = b + rep(0, max.l);
    
    g <- function(x,y) { ## y^(x-1)*(e^y - 1 - y), x > -1, y >=0
        eps=10^-10; y0=300;
        val = rep(0, length(y)); log.val=log(val);
        I = which(y>eps & y<=y0);
        val[I] = y[I]^x * (expm1(y[I])/y[I]-1);
        log.val[I] = log(val[I]);
        I = which(y<=eps);
        val[I] = y[I]^(x+1)*(0.5  + y[I]/6 + y[I]^2*(1 + y[I]/5)/24);
        log.val[I] = log(val[I]);
        I = which(y>y0);
        log.val[I]=y[I]+(x-1)*log(y[I]) + log(-expm1(-y[I] + log(1+y[I])));
        val[I] = exp(log.val[I]);
        return(list(val=val, log.val=log.val));
    }

    val = rep(0, length(a)); log.val=log(val);
    I = which(a<b);
    if (c>10^-10) {
        .b=g(c,b[I]); .a = g(c, a[I]); .v = (b[I]^c - a[I]^c)/c;
        val[I] = 2*(.b$val - .a$val)  + .v;
        log.val[I] = .b$log.val  +
            log(-2*expm1(.a$log.val - .b$log.val) + exp(log(.v) - .b$log.val));
    } else {
        val[I] = 2*(g(c,b[I]) - g(c,a[I]));
        J = I[which(a[I]==0)];
        val[J] = val[J] + b[J]^c/c;

        J = I[which(a[I]>0)];
        lb=log(b[J]); la=log(a[J]);
        val[J] = val[J] + lb/x.over.expm1(c*lb) - la/x.over.expm1(c*la);
        log.val[I] = log(val[I]);
    }

    return(list(val=val, log.val=log.val));
} ## M.abc

phi.star.abc.par <- function(a,b,c)
### Calculate the parameters in $\phi^*_{a,b,c}(x)$ in Eq. (4.12).  Same
### requirements on dimensions of input arguments as in M.abc(a, b, c)
{
    if (length(a)==0 | length(b)==0) return(numeric(0));
    if (length(c)!=1 | c<0 | c>=1) {
        stop(bold(red("c")), " must be a scalar in (0,1)");
    }
    if (any(c(a,b)<0)) {
        stop("All entries of ", bold(red("a")), " and ", bold(red("b"))," must be nonnegative");
    }
    max.l = max(length(a), length(b));
    a = a + rep(0, max.l);
    b = b + rep(0, max.l);
    
    cc = 1+sqrt(1-c);
    xc=2*cc-c; pc = 1/cc;  ## x_c and p_c in Eq. (4.7)

    rn <- function(x,n, method=1) {
        ## r_n(x) below (4.8); ONLY r_7(4) is needed
        if (method==1) { ## OK for median x and n, definitely for x=4,n=7)
            r=expm1(x);
            if (n>0) {
                r= r-sum(x^(1:n)/factorial(1:n));
            }
            return(r*factorial(n+1)/x^(n+1));
        }
        if (method==2) {
            z=1;
            zs=z;
            while ( abs(z) > 10^-20 ) {
                z = z*x/(n+2);
                zs = c(zs, z);
                n=n+1;
            }
            return(sum(zs));
        }
        if (method==3) {
            fun <- function(t) { exp(x*t)*(1-t)^n; }
            return((n+1)*integral(fun, 0, 1, no_intervals=0));
        }
    }
    
    A = pmin(b,xc);  ## A below (4.7)

    ## w in Eq. (4.11): note the index in Eq. (4.11) starts 0 and ends 8;
    ## w is a matrix, one column per entry in a (or b)
    w = matrix(nrow = 9, ncol = max.l);
    I=which(a==0);
    w[1,I] = A[I]^c/c;
    I=which(a>0);
    log.A=log(A[I]); log.a=log(a[I]);
    w[1,I] = log.A/x.over.expm1(log.A*c) - log.a/x.over.expm1(log.a*c);
    for ( k in 1:8 ) {
        w[k+1,] = (A^(k+c) - a^(k+c))/((k+c)*factorial(k));
     }
    w[9,]=rn(4,7)*w[9,];

    return(list(pc=pc, xc=xc, w=w));
} ## phi.star.abc.par

phi.star.abc <- function(x,a,b,c)
### Calculate $\phi^*_{a,b,c}(x)$ in Eq. (4.12).  The input arguments x, a,
### and b are arrays or scalars.  All the entires in a and b must be
### nonnegative, and c must be a scalar in (0,1).  If x, a, and b have
### different lengths, the shorter ones internally are padded with their
### repeats to make their lengths equal to the longest one.  The output is a
### vector of the same length as the (padded) x, a, and b, such that its i-th
### entry is calculated for x[i], a[i], b[i], and c.
{
    if (length(x)==0 | length(a)==0 | length(b)==0) return(numeric(0));
    if (length(c)!=1 | c<=0 | c>=1) {
        stop(bold(red("c")), " must be a scalar in (0,1)");
    }
    if (any(c(a,b)<0)) {
        stop("All entries of ", bold(red("a"))," and ", bold(red("b"))," must be nonnegative");
    }
    max.l = max(length(x), length(a), length(b));
    x = x + rep(0, max.l);
    a = a + rep(0, max.l);
    b = b + rep(0, max.l);

    par=phi.star.abc.par(a,b,c);
    A = pmin(b, par$xc); B=pmax(a, par$xc); ## A and B below Eq. (4.7)

    val=rep(0, max.l);

    ## val.s stores the \phi^*_{a,b,c}(x)*|x|^(1-c).  It is used in Eq. (7.1)
    ## to calculate R(v), where c is $\delta$.
    val.s=rep(0, max.l);
    
    m = M.abc(B, b, c);

    I = which(m$val>0 & x!=0);
    if (length(I)>0) {
        y = abs(x[I])^(1/par$pc);  ## par$pc is the scalar $p_c$ in Eq. (4.7)
        ## For the next line, see below (4.12)
        s=exp( (c-1/par$pc)*log(b[I]) + b[I] - m$log.val[I]) * par$pc;
        I = I[which(b[I]^(1/par$pc)>=sign(x[I]*y))];
        .tx = truncated.exp(b[I]^(1/par$pc) - sign(x[I])*y, s);
        log.v = m$log.val[I] + log(.tx);
        val[I] = 2/par$pc * (y/abs(x[I])) * exp(log.v);
        val.s[I] = val[I] * abs(x[I])^(1-c);
    }
    
    m = M.abc(a,A,c)$val;
    I = which(m>0 & a<=x & x<=A);

    if (length(I)>0) {
        betax = c/(expm1(c*log(A[I])) - expm1(c*log(a[I])));
        for ( k in 1:8 ) {
            betax = rbind(betax, (k+c)*x[I]^k / (A[I]^(k+c) - a[I]^(k+c)));
        }
        ## par$w stores the $w_k$ in (4.11); a matrix, one column for each 
        ## pair of entries of _a_ and _b_
        .mw = 2*m[I]*
            colSums(as.matrix(par$w[,I]*betax))/colSums(as.matrix(par$w[,I]));
        val[I] = val[I] + .mw*x[I]^(c-1);
        val.s[I] = val.s[I] + .mw;
    }
    return(list(val=val, val.s=val.s));
} ## phi.star.abc

sample.phi.star.abc <- function(a,b,c)
### Algorithm 4.1.  Same requirement on the dimensions of the input
### arguments as in M.abc(a,b,c).  Both the sampled values and their
### logarithms are returned.  If a sample value is non-positive, then its
### logarithm is returned as -Inf.
{
    if (length(a)==0 | length(b)==0) {
        return(list(val=numeric(0), log.val=numeric(0)));
    }

    if (length(c)!=1 | c<=0 | c>=2) {
        stop(bold(red("c")), " must be a scalar in (0,2)");
    }
    max.l = max(length(a), length(b));
    a = a + rep(0, max.l);
    b = b + rep(0, max.l);
    if (any(a>=b) | any(a<0)) {
        stop("In sample.phi.star.abc(", bold(red("a,b,c")), "): entries of ", bold(red("a"))," must be nonnegative and less than the corresponding entries of ", bold(red("b")));
    }
    
    par=phi.star.abc.par(a,b,c);
    pc=par$pc; A = pmin(b, par$xc); B = pmax(a, par$xc);

    log.A = log(A); 
    
    M0=M.abc(a,A,c); M1=M.abc(B,b,c);

    val = rep(0, max.l); log.val = log(val);
    i = (runif(max.l)<1/(1+exp(M1$log.val - M0$log.val)));
    I = which(i);  len.I = length(I);
    if (len.I > 0) { 
        k = sample.mat(par$w[,I])-1;
        l.a = (a[I] / A[I])^(k+c);
        log.val[I] = log.A[I] + log(runif(len.I)*(1-l.a) + l.a)/(k+c);
        val[I] = exp(log.val[I]);
    }
    
    I = which(!i); len.I = length(I);
    if (len.I > 0) {
        s = pc*exp((c-1/pc)*log(b[I]) + b[I] - M1$log.val[I]);
        Y = sample.truncated.exp(s);
        
        D = b[I]^(1/pc) - Y;
        val[I] = abs(D)^pc*sign(D);
        log.val[I] = log(pmax(0, val[I]));
    }
    
    return(list(val=val, log.val=log.val, log.M=cbind(M0$log.val, M1$log.val)));
} ## sample.phi.star.abc

G.abc <- function(a,b,c)
### Evaluate the function $G(a,b,c)$ in Eq. (4.15).  Same requirements on the
### input arguments as in M.abc(a,b,c)
{
    if (length(a)==0 | length(b)==0) return(numeric(0));
    if (length(c)!=1 | c<0 | c>=2) {
        stop("In G.abc(", bold(red("a,b,c")), "): ", bold(red("c")), " must be a scalar in (0,2)");
    }
    if (any(c(a,b)<0)) {
        stop("In G.abc(", bold(red("a,b,c")), "): All entries of ", bold(red("a")), " and ", bold(red("b")), " must be nonnegative");
    }
    max.l = max(length(a), length(b));
    a = a + rep(0, max.l);
    b = b + rep(0, max.l);

    val=rep(0, length(a));
    I = which(a<b & b<=a+1);
    if (length(I)>0) {
        J = I[which(a[I]==0)];
        val[J] = b[J]^c/c;
        J = I[which(a[I]>0)];
        lb=log(b[J]); la=log(a[J]);
        val[J] = exp(-a[J])*(lb/x.over.expm1(lb*c) - la/x.over.expm1(la*c));
    }

    I = which(a+1<b & b<=a+2);
    if (length(I)>0) {
        val[I] = G.abc(a[I], a[I]+1, c) + G.abc(a[I]+1, b[I], c);
    }

    I = which(b>a+2);
    if (length(I)>0) {
        a2 = a[I]+2;
        if (c>=1) {
            val[I]= 2*(a2^(c-1)*exp(-a2) - b[I]^(c-1)*exp(-b[I]));
        } else {
            val[I]=a2^(c-1)*(exp(-a2) - exp(-b[I]));
        }
        val[I]=G.abc(a[I], a2, c) + val[I];
    }
    
    return(val);
} ## G.abc

g.star.abc <- function(x,a,b,c)
### Calculate $g^*_{a,b,c}(x)$ in Eq. (4.16).  Same requirement on the input
### arguments as phi.stat.abc(x,a,b,c)
{
    if (length(x)==0 | length(a)==0 | length(b)==0) return(numeric(0));
    if (length(c)!=1 | c<=0 | c>=2) {
        stop("In g.star.abc(", bold(red("x,a,b,c")), "): ", bold(red("c")), " must be a scalar in (0,2)");
    }
    if (any(c(a,b)<0)) {
        stop("In g.star.abc(", bold(red("x,a,b,c")), "): All entries of", bold(red("a")), " and ", bold(red("b")), " must be nonnegative");
    }
    max.l = max(length(x), length(a), length(b));
    x = x + rep(0, max.l);
    a = a + rep(0, max.l);
    b = b + rep(0, max.l);

    val = rep(0, max.l);
    if (c>=1) {
        A = pmin(b,c-1); B = pmax(a, c-1);

        I = which( A>a & A>=x);
        if (length(I)>0)  {
            gx = G.abc(a[I], A[I], c);
            s = A[I]^(c-1)*exp(-A[I])/gx;
            val[I] = 2*gx*truncated.exp(A[I] - x[I], s);
        }

        I = which( B<b & x> B);
        if (length(I)) {
            gx = G.abc(B[I], b[I], c);
            t = B[I]^(c-1) * exp(-B[I])/gx;
            val[I] = 2*gx*truncated.exp(x[I] - B[I], t);
        }
    } else {
        I = which(x>a);
        if (length(I)>0) {
            gx = c*G.abc(a[I], b[I], c);
            s = exp(-a[I])/gx;
            xc=x[I]^c;
            val = 2*gx*(xc/x[I])*truncated.exp(xc - a[I]^c, s);
        }
    }
    return(val);        
} ## g.star.abc

sample.g.star.abc <- function(a,b,c)
### Algorithm 4.2.  The input arguments a and b are arrays.  If they have
### different lengths, the shorter one internally are padded with its
### repeats to have the same length as the longer one.  The the entries of
### (padded) a must be nonnegative and strictly less than the corresponding
### ones in (padded) b.  The input argument c must be a scalar in (0,2).
### Both the sampled values and their logarithms are returned.
{
    if (length(a)==0 | length(b)==0) {
        return(list(val=numeric(0), log.val=numeric(0)));
    }
    if (length(c)!=1 | c<=0 | c>=2) {
        stop("In sample.g.star.abc(", bold(red("a,b,c")), "):", bold(red("c")), " must be a scalor in (0,2)");
    }
    max.l = max(length(a), length(b));
    a = a + rep(0, max.l);
    b = b + rep(0, max.l);

    if (any(a<0 | a>=b)) {
        stop(c("Entries of _a_ must be nonnegative and less than the",
               " corresponding\n  entries of _b_"));
    }

    val = rep(0, max.l); log.val = val;
    if (c>=1) {
        A = pmin(b,c-1); B = pmax(a,c-1);
        M0 = G.abc(a,A,c); M1=G.abc(B,b,c);
        i = (runif(max.l) < M0/(M0+M1));
        
        I = which(i); 
        if (length(I)>0) {
            s = A[I]^(c-1) * exp(-A[I])/M0[I];
            val[I] = A[I] - sample.truncated.exp(s);
        }

        I = which(!i);
        if (length(I)>0) {
            t = B[I]^(c-1) * exp(-B[I])/M1[I];
        }
        val[I] = B[I] + sample.truncated.exp(t);
        log.val = log(pmax(0, val));
    } else {
        s = exp(-a)/(c*G.abc(a,b,c));
        log.val = log(a^c + sample.truncated.exp(s))/c;
        val = exp(log.val);
    }
    return(list(val=val, log.val=log.val));
} ## sample.g.star.abc

kappa.alpha <- function(alpha)
### $\kappa_{\alpha,i}$, $i=1,2,3,4$, in Eq. (5.9) and Prop. 6.2
{
    if (length(alpha)!=1 | alpha<=0 | alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in (0,1)");
    }
    delta = 1-alpha;

    c.a=c.alpha(alpha);

    kappa.alpha = rep(0,4);
    kappa.alpha[1] = 2*c.a*(1/delta-2);
    kappa.alpha[2] = c.a*(4 + 1/exp(1));
    kappa.alpha[3] = kappa.alpha[1] + (log(2))^(-alpha)*kappa.alpha[2] + 1;
    kappa.alpha[4] = kappa.alpha[1]*(log(2))^delta + kappa.alpha[2];

    return(kappa.alpha);
} ## kappa.alpha

psi.alpha <- function(t, alpha)
### $\psi_\alpha(t)$ in Eq. (5.9).  The input argument t is an array of
### nonnegative numbers adn alpha must be a scalar in (0,1)
{
    if (length(t)==0) return(numeric(0));
    if (any(t<0)) {
        stop("All entries of _t_ cannot have negative entries");
    }
    if (length(alpha)!=1 | alpha<=0 | alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in (0,1)");
    }

    k.a = kappa.alpha(alpha);
    delta = 1-alpha;

    val = rep(1, length(t));
    val[which(is.infinite(t))]=Inf;
    
    I = which(t>0 & is.finite(t));

    ell = log1p.recip(t[I]);
    val[I] = 1 + k.a[1]*t[I]*ell^delta + k.a[2]*ell^(-alpha);
    return(val);
} ## psi.alpha

Q.alpha.z <- function(theta, z, alpha)
### $Q_{\alpha,z}$ in Eq. (5.10).  The input arguments theta and z are
### arrays.  If they have different lengths, the shorter one internally
### are padded with its repeats to have the same length as the longer one.
### All the entries of theta must be in [0, pi), all the entries of z
### must be nonnegative, and alpha must be a scalar in (0,1)
{
    nz = length(z);
    if (length(theta)==0 | nz==0) return(numeric(0));
    if (any(theta<0 | theta>=pi)) {
        stop("entries of theta must be in [0,pi)");
    }
    if (any(z<0)) {
        stop("entries of _z_ must be nonnegative");
    }
    if (length(alpha)!=1 | alpha<=0 | alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in (0,1)");
    }

    max.l = max(length(theta), nz);
    theta = theta + rep(0, max.l);
    z = z + rep(0, max.l);

    tau = rep(0, max.l);
    I = which(z>0);
    tau[I] = exp(log(z[I]) + log.H.alpha(theta[I], alpha));
    psi = psi.alpha(tau, alpha);

    val = rep(0, max.l);
    I = which(psi<Inf);
    val[I] = exp(log(psi[I]) - tau[I]);

    return(val);
} ## Q.alpha.z

get.t.n <- function(theta0, Delta)
### Compute $t_0$, ... $t_m$, $t_{m+1} = \theta_0$ in Eq. (6.2).  Note
### they are indexed by 1, ... m+2 in the output.  Both input arguments
### theta0 and Delta must be are scalars.  
{
    f <- function(theta) { ## cf. below (5.2)
        sx = sinc(theta);
        return(1 + 1/(sx*sx) - 2*cos(theta)/sx);
    }
    
    L=log1p(Delta);
    T=theta0;
    cnt=1;
    while (T[cnt]>0) {
        t = T[cnt]*(1 - L/f(T[cnt]));
        t = max(t,0);
        T=c(T,t); cnt=cnt+1
    }
    T=T[seq(cnt,1, by=-1)];

    return(list(t.n=T, theta0=theta0, Delta=Delta));
} ## get.t.n

get.theta.n <- function(alpha, alpha0, theta0, Delta)
### Compute $\theta_0$, ..., $\theta_N$, $\theta_{N+1}=\pi$ in Prop. 6.1.
### Note they are indexed by 1, ... N+2 in the output.  All the input
### arguments must be scalars.
{
    if (length(alpha)!=1 | alpha<=0 | alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in (0,1)");
    }

    f.alpha <- function(theta) {  ## f_alpha in Prop. 6.1
        A = (1+Delta/2)*theta/(pi - alpha*theta);

        C = pi*(1/(sin(alpha0*theta0))^2 - 1/(pi - alpha0*theta0)^2);
        ## cf. the second before the last formula in Prop. 2.2

        theta1 = A*pi/(1+alpha*A);
        theta2 = theta + log1p(Delta/(2+Delta))/C;
        return(min(c(theta1,theta2,pi)));
    }
    T=theta0;
    cnt=1;
    while(T[cnt]<pi) {
        T = c(T, f.alpha(T[cnt]));
        cnt = cnt+1;
    }
    return(list(theta.n=T, alpha=alpha, alpha0=alpha0, theta0=theta0,
                Delta=Delta));
} ## get.theta.n

log.J.alpha <- function(theta, theta.n.struct)
### Log of $J_\alpha(\theta)$ in Eq. (6.6).  The input argument theta is
### an array, and theta.n.struct should be generated by get.theta.n().
{
    if (length(theta)==0) return(numeric(0));
    alpha = theta.n.struct$alpha;
    theta.n = theta.n.struct$theta.n;
    delta = 1-alpha;

    Ka = rep(0, length(theta));
    for ( i in 1:(length(theta.n)-1) ) {
        I = which(theta>=theta.n[i] & theta<theta.n[i+1]);
        Ka[I] = K.alpha(theta.n[i], alpha);
    }
    
    val = rep(-Inf, length(theta));
    I = which(Ka!=0);
    val[I] = log(Ka[I]) + log1p((delta*pi/alpha)/(pi - theta[I]))/delta;

    return(val);
} ## log.J.alpha

get.theta.z <- function(z, t.n.struct, theta.n.struct)
### Compute $\theta_z$ in Eq. (6.7).  The input argument z is an array with
### all entries positive, t.n.struct should be generated by get.t.n(), and
### theta.n.struct by get.theta.n()
{
    nz = length(z);
    if (any(z<=0)) {
        stop("All entries of z must be strictly positive");
    }
    
    t.n=t.n.struct$t.n;
    theta0 = t.n.struct$theta0;
    Delta = t.n.struct$Delta;
    
    if (theta0!=theta.n.struct$theta0 | Delta != theta.n.struct$Delta) {
        stop("_t.n.struct_ and _theta.n.struct_ not compatiable");
    }
    alpha = theta.n.struct$alpha;
    alpha0 = theta.n.struct$alpha0;
    theta.n = theta.n.struct$theta.n;
    
    val = rep(0, nz);

    if (nz>0) {
        log.z = log(z);
    }
    else {
        log.z = numeric(0);
    }
    I = which(z<1 & log.z + log.H.alpha(theta0, alpha) >=0);
    if (length(I)>0) {
        for ( i in 1:length(t.n) ) {
            J = which(log.H.alpha(t.n[i], alpha) + log.z[I]>=0);
            if (length(J)>0) {
                val[I[J]] = t.n[i];
                I = I[-J];
            }
            if (length(I)==0) { break; }
        }
    }

    I = which(z<1 & log.z + log.H.alpha(theta0, alpha)<0);
    if (length(I)>0) {
        delta = 1-alpha;
        for ( i in seq(length(theta.n)-1, 1, by=-1)) {
            J = which(log.H.alpha(theta.n[i], alpha) + log.z[I] <= 0);
            if (length(J)>0) {
                x = -delta*(log.z[I[J]] + log(K.alpha(theta.n[i], alpha)));
                val[I[J]] = pmin((1-delta/alpha/expm1(x))*pi, theta.n[i+1]);
                I = I[-J];
            }
            if (length(I)==0) { break; }
        }
    }

    return(list(z=z, theta.z=val, theta.n=theta.n, alpha=alpha,
                theta0=theta0, alpha0=alpha0, Delta=Delta));
} ## get.theta.z

omega.s.z <- function(theta.z.struct)
### Compute $\omega_z$ and $s_z$ in Prop. 6.2.  The input argument
### theta.z.struct should be generated by get.theta.z().
{
    z = theta.z.struct$z;
    theta.z = theta.z.struct$theta.z;
    alpha = theta.z.struct$alpha;

    nz = length(z);
    if (nz==0) {
        return(list(omega.z=numeric(0), s.z=numeric(0)));
    }

    log.tau = log(z) + log.H.alpha(theta.z, alpha);
    tau = exp(log.tau);
    dlogH = dlog.H.alpha(theta.z, alpha);
    k.a = kappa.alpha(alpha);
    omega.z = 2*k.a[3]/(tau*dlogH);
    s.z = dlogH * exp((1+alpha)*log.tau - tau);

    return(list(omega.z=omega.z, s.z = s.z));
} ## omega.z

get.acdpi.n <- function(theta.z.struct)
### Set up $c_n$, $d_n$ in Eq. (6.13), $a_n$ in in Eq. (6.14), and $\pi_n$
### in Eqs. (6.15)-(6.16).  The input argument theta.z.struct should be
### generated by get.theta.z().  It is a list that contains the following
### arrays: z, theta.z, and theta.n.  The output is a list of four matrices
### of the same size, consisting of values of $a_n$, $c_n$, $d_n$, $\pi_n$,
### repsectively.  The length of each column is the length of theta.n minus
### 1.  The i-th column in each matrix corresponds to z[i].  If z[i] has
### $\theta_z<=\theta_n$, then all the i-th cocolumns of the matrices are 0.
### [This is because the corresponding interval $I_n$ is empty; cf.
### comment below Eq. (6.11)].  NOTE: theta.n and theta.z store $\theta_n$
### and $\theta_z$ defined in Prop. 6.1 and Eq. (6.7); $b_n$ in Eq. (6.14)
### is only used to define $\pi_n$, so is not output.
{
    val = theta.z.struct;

    z = theta.z.struct$z;
    nz = length(z);
    if (nz==0) {
        val$a.n=numeric(0); val$c.n=numeric(0);
        val$d.n=numeric(0); val$pi.n=numeric(0);
        return(val);
    }

    alpha = theta.z.struct$alpha;
    theta.z = theta.z.struct$theta.z;
    theta.n = theta.z.struct$theta.n;
    
    delta = 1 - alpha;

    ## theta.n stores $\theta_0$, ... $\theta_N$, $\theta_{N+1}=\pi$ in
    ## Prop. 6.1, so it has N+2 entries
    N = length(theta.n)-1;  ## This N is one plus the N in Prop. 6.1

    a.n = matrix(rep(0, N*nz), nr=N);
    c.n = matrix(rep(0, N*nz), nr=N);
    d.n = matrix(rep(0, N*nz), nr=N);
    pi.n = matrix(rep(0, N*nz), nr=N);
    
    log.K = log(K.alpha(theta.n, alpha));
    log.z = log(z);

    b.n = delta*pi/(pi - alpha*theta.n);
    l.n = log1p(delta*pi/(alpha*(pi - theta.n)));
    for ( n in 1:N) {
        I = which(theta.z>theta.n[n]); ## cf. between Eqs. (6.11) and (6.12)
        if (length(I)>0) {
            J = I[which(theta.z[I] >= theta.n[n+1])];
            a.n[n,J] = exp(log.z[J] + log.K[n] + l.n[n+1]/delta);

            J = I[which(theta.z[I] < theta.n[n+1])];
            a.n[n,J] = 1;
            
            c.n[n,I] = log1p.recip(a.n[n,I]);
            d.n[n,I] = log1p.recip(exp(log.z[I] + log.K[n] + l.n[n]/delta));
  
            if (theta.n[n]<=(1-delta/alpha)*pi) {
                pi.n[n,I] = (1+a.n[n,I])^2*alpha*(pi - theta.n[n])^2/pi*
                    exp(-delta*(log.z[I] + log.K[n]));
            } else {
                pi.n[n,I] = (1+a.n[n,I])^2*delta^2*pi*
                    exp(delta*(log.z[I] + log.K[n]))/(alpha*b.n[n]^2)
            }
        }
    }

    val$a.n=a.n; val$c.n=c.n; val$d.n=d.n; val$pi.n=pi.n;
    
    return(val);
} ## get.acdpi.n

integral.P.star <- function(acdpi.struct)
### The integrals of $P^*_{n,i}$, $i=1,2,3$, in Prop. 6.5.  [The formulas
### are in a subsequent paragraph.]  The input argument acdpi.struct should
### be generated by get.acdpi.n().  It is a list that contains matrices a.n,
### c.n, and d.n.  These matrices have the same size.  The output contains
### three matrices, consisting of the integrals of $P^*_{n,i}$, $i=1,2,3$.
### Each matrix has the same size as a.n.
{
    a.n = acdpi.struct$a.n;
    c.n = acdpi.struct$c.n;
    d.n = acdpi.struct$d.n;
    pi.n = acdpi.struct$pi.n;
    alpha = acdpi.struct$alpha;
    theta.n = acdpi.struct$theta.n; # last entry of _theta.n_ is \pi
    
    N = nrow(d.n); # N = length of theta.n minus 1
                    
    int.P1 = matrix(0, nr=N, nc=ncol(d.n));
    int.P2 = int.P1; int.P3 = int.P1;

    delta = 1-alpha;
    dp1 = 1+delta;
    k.a = kappa.alpha(alpha);
    I = which(theta.n[1:N] <= (1-delta/alpha)*pi); n.I=length(I);
    if ( n.I > 0 ) {
        int.P1[I,] = 2*exp(-dp1*log1p(delta)) * (1+a.n[I,]) * k.a[1] *
            matrix(G.abc(dp1*c.n[I,], dp1*d.n[I,], dp1), nr = n.I);
        int.P2[I,] = 2*k.a[2] * exp(-delta*log(delta)) *
            matrix(G.abc(delta*c.n[I,], delta*d.n[I,], delta), nr = n.I);
        int.P3[I,] = (expm1(-delta*c.n[I,])-expm1(-delta*d.n[I,]))/delta;
    }
    
    I = which(theta.n[1:N] > (1-delta/alpha)*pi); n.I = length(I);
    if (n.I > 0 ) {
        int.P1[I,] = 2 * exp(-dp1*log1p(-delta)) * k.a[1] *
            matrix(G.abc(alpha*c.n[I,], alpha*d.n[I,], dp1), nr = n.I);
        int.P2[I,] = 2 * k.a[2] * exp(-delta*log(delta)) *
            matrix(M.abc(delta*c.n[I,], delta*d.n[I,], delta)$val, nr = n.I);
        int.P3[I,] = (expm1(delta*d.n[I,]) - expm1(delta*c.n[I,]))/delta;
    }

    return(list(int.P1 = int.P1, int.P2 = int.P2, int.P3 = int.P3));
} # integral.P.star

P.star.n <- function(tx, thetax, ax, cx, dx, pix, alpha)
### Calculate $P^*_n(t)$ in Prop. 6.5.  The input arguments tx, thetax, ax,
### cx, dx, and pix are arrays and must have the same length, correponding
### to a 6-tuple $(t, \theta_n, a_n, c_n, d_n, pi_n)$ in the formulas.
### Also, from Algorithm 6.2, if an entry of tx is not between the
### corresponding entries of cx and dx, the exact value of $P^*_n(t)$ is not
### needed, so the returned value is NA.  The returned value consists of a
### 3-row matrix with the i-th row consisting of the values of
### $P^*_{n,i}(t)$, and a vector of the values of $P^*_n(t)$.  [For different
### entries of tx, the value of n in the subscript may be different but
### need not be specified.]
{
    if (length(alpha)!=1 | alpha<=0 | alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in (0,1)");
    }
    if (length(tx)==0) { return(numeric(0)); }
    l = c(length(tx), length(thetax), length(ax), length(cx), length(dx),
          length(pix));
    if (max(l) > min(l)) {
        stop("The first 6 input arrays must have the same length");
    }

    delta=1-alpha;
    dp1 = delta + 1;
    k.a = kappa.alpha(alpha);

    P = matrix(rep(NA, 3*length(tx)), nr=3); 
    I = which(thetax<=(1-delta/alpha)*pi & tx>cx & tx<=dx);
    if (length(I)>0) {
        P[,I]=rbind(exp(-delta*log1p(delta))*(1 + ax[I])*k.a[1]*
                    g.star.abc(dp1*tx[I], dp1*cx[I], dp1*dx[I], dp1),
                    k.a[2]*delta^alpha*
                    g.star.abc(delta*tx[I], delta*cx[I], delta*dx[I], delta),
                    exp(-delta*tx[I]));
    }

    I = which(thetax>(1-delta/alpha)*pi & tx>cx & tx<=dx);
    if (length(I)>0) {
        P[,I]=rbind(exp(-delta*log1p(-delta))*k.a[1]*
                    g.star.abc(alpha*tx[I], alpha*cx[I], alpha*dx[I], dp1),
                    k.a[2]*delta^alpha*
                    phi.star.abc(delta*tx[I], delta*cx[I], delta*dx[I],
                                 delta)$val,exp(delta*tx[I]));
    }

    return(list(mat=P, val=pix*colSums(P)));
} # P.star.n

sample.P.star.n <- function(thetax, cx, dx, int.P.star, alpha)
### Algorithm 6.1 to sample from the normalized $P^*_n(t)$.  The input
### arguments thetax, cx, and dx are arrays of the same length,
### corresponding to the triple $(\theta_n, c_n, d_n)$ in the algorithm.
### Each entry of cx must be strictly smaller than the corresponding entry
### of dx.  The input argument int.P.star has 3 rows and its number of
### columns is the same as the length of cx, each column consisting of the
### integrals of $\int^{d_n}_{c_n} P^*_{n,i}$, $i=1,2,3$.  Each column sum
### of int.P.star must be strictly positive.
{
    if (length(cx) == 0) { return(numeric(0)); }
    if (length(alpha)!=1 | alpha<=0 | alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in (0,1)");
    }
    l = c(length(thetax), length(cx), length(dx), ncol(int.P.star));
    if (max(l) > min(l)) {
        stop("_thetax_, _cx_, _dx_, and the rows of _int.P.star_",
             " must have the same length");
    }
    if (any(cx>=dx | colSums(int.P.star)<=0)) {
        stop("Entries of _cx_ must be strictly less the corresponding",
             " entries of _dx_, and all column sums of _int.P.star_ must",
             " be positive");
    }
    delta=1-alpha;
    dp1 = delta + 1;
    
    X = rep(0, length(cx));
    intvl = sample.mat(int.P.star);

    ## For i=3 in Algorithm 7.1
    wx = log1p(runif(length(cx))* expm1(delta*(dx-cx)))/delta
    
    I = which(thetax<=(1-delta/alpha)*pi);
    if (length(I)>0) {
        J = I[which(intvl[I]==1)];
        if (length(J)>0) {
            X[J] = sample.g.star.abc(dp1*cx[J], dp1*dx[J], dp1)$val/dp1;
        }
        J = I[which(intvl[I]==2)];
        if (length(J)>0) {
            X[J] = sample.g.star.abc(delta*cx[J], delta*dx[J], delta)$val/delta;
        }
        J = I[which(intvl[I]==3)];
        X[J] = dx[J] - wx[J];
    }

    I = which(thetax>(1-delta/alpha)*pi);
    if (length(I)>0) {
        J = I[which(intvl[I]==1)];
        if (length(J)>0) {
            X[J] = sample.g.star.abc(alpha*cx[J], alpha*dx[J], dp1)$val/alpha;
        }
        J = I[which(intvl[I]==2)];
        if (length(J)>0) {
            X[J] = sample.phi.star.abc(delta*cx[J],
                                       delta*dx[J],delta)$val/delta;
        }
        J = I[which(intvl[I]==3)];
        X[J] = cx[J] + wx[J];
    }
    return(X);
} ## sample.P.star.n

sample.Q.alpha.z <- function(z, alpha, alpha0, theta0, Delta)
### Algorithm 6.2 to sample from the normalized $Q_{\alpha,z}$.  The input
### argument z is an array of values in (0,1).  [For z>=1, use the other
### two algorithms.]  The input argument alpha must be a scalar in
### [alpha0, 1), and alpha0, theta0, Delta correspond to $\alpha_0$,
### $\theta_0$, and $\Delta$ in Eq. (6.1).
{
    nz = length(z);
    if (nz==0) return(numeric(0));
    if (any(z<=0 | z>=1)) {
        stop("All entries of _z_ must be in (0,1)");
    }
    
    if (length(Delta)!=1 | Delta<=0 | Delta>=1) {
        stop("_Delta_ must be a scalar in (0,1)");
    }
    if (length(alpha0)!=1 | alpha0<=.5 | alpha0>=1) {
        stop("_alpha0_ must be a scalar in (0.5, 1)");
    }
    if (length(theta0)!=1 | theta0<=pi/3*(1+1/alpha0) | theta0>=pi) {
        stop("_theta0_ must be a scalar in (pi/3+pi/(3*_alpha0_), pi)");
    }
    if (length(alpha)!=1 | alpha<alpha0| alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in [_alpha0_,1)=[", alpha0, '1)');
    }
    delta = 1-alpha;
    log.z = log(z);
    
    k.a = kappa.alpha(alpha);
    t.n.struct = get.t.n(theta0, Delta);
    theta.n.struct = get.theta.n(alpha, alpha0, theta0, Delta);
    theta.z.struct = get.theta.z(z, t.n.struct, theta.n.struct);
    acdpi.struct = get.acdpi.n(theta.z.struct);

    t.n = t.n.struct$t.n;
    lh.t.n = log.H.alpha(t.n, alpha);
    
    theta.n = theta.n.struct$theta.n;
    lk.theta.n = log(K.alpha(theta.n, alpha));
    
    theta.z = theta.z.struct$theta.z;

    ## \int Q^*_E
    osz.struct = omega.s.z(theta.z.struct);
    int.Q.star = osz.struct$omega.z;

    ## \int Q^*_{D_i}
    for ( i in 1:(length(t.n)-1) ) {
        Ints = rep(0, nz);
        I = which(theta.z > t.n[i]);
        if (length(I) > 0 ) {
            tau = exp(log.z[I] + lh.t.n[i]);
            Ints[I] = (1+Delta)*psi.alpha(tau, alpha) *
                (pmin(t.n[i+1], theta.z[I]) - t.n[i]);
        }
        int.Q.star = rbind(int.Q.star, Ints);
    }

    ## \int Q^*_{I_n}
    int.P.star = integral.P.star(acdpi.struct);
    sum.int.P.star = int.P.star$int.P1+int.P.star$int.P2+int.P.star$int.P3;
    int.Q.star = rbind(int.Q.star, (1+Delta)*acdpi.struct$pi.n*sum.int.P.star);

    val = rep(0, nz); ## storage of sampled value

    ### Auxilliary; see Algorithm 6.2 for justification
    S = rep(0, nz); ## sampled interval indices, one per entry of z
    U = rep(0, nz); ## storage of samples from unif(0,1)
    r = rep(Inf, nz); ## storage of ratios
    cache = rep(0, nz); ## storage of samples from various envelopes

    I = 1:nz;  ## Indices of _z_ with no corresponding sampled values

    while (length(I)>0) {
        ## For each _z_ with no corresponding sampled values, sample from
        ## intervals {E, D_0, ... D_m, I_0, ..., I_N} (Algorithm 6.2, line 2)
        S[I] = sample.mat(int.Q.star[,I]); 
        U[I] = runif(length(I));
        r[I] = Inf;

        ## Sampled interval is E.
        J = I[which(S[I]==1)]; n.E=length(J);
        ## In E
        if (n.E>0) {
            ## Use cache to store Z (Algorithm 6.2, line 5)
            cache[J] = sample.truncated.exp(osz.struct$s.z[J]);
            val[J] =  cache[J] + theta.z[J];
            sub = J[which(val[J]<pi)];
            if (length(sub)>0) {
                r[sub] = U[sub] * osz.struct$omega.z[sub] *
                    truncated.exp(cache[sub], osz.struct$s.z[sub])/
                    Q.alpha.z(val[sub], z[sub], alpha);
            }
        }

        ## Sampled interval is D_i
        J = I[which(S[I]>=2 & S[I]<=length(t.n))]; n.D=length(J);
        if (n.D>0) {
            ## Left and right end points of D_i
            .a = t.n[S[J]-1]; .b = pmin(t.n[S[J]], theta.z[J]);
            val[J] = .a + (.b - .a)*runif(length(J));
            r[J] = U[J]*(1+Delta)*
                psi.alpha(exp(log.z[J] + lh.t.n[S[J]-1]), alpha)/
                Q.alpha.z(val[J], z[J], alpha);
        }

        ## Sampled interval is I_n
        J = I[which(S[I] > length(t.n))]; n.I=length(J);
        if (n.I>0) {
            .i = S[J] - length(t.n); ## indices (plus 1) of sampled I_n
            
            ## The next are the indices of the selected entries in the matrices
            ## in acdpi.struct, each treated as an array
            .s = (J-1)*(length(theta.n)-1) + .i;
            
            .int = rbind(int.P.star$int.P1[.s], int.P.star$int.P2[.s],
                         int.P.star$int.P3[.s]);
            ## Use cache to store t (Algorithm 6.2, line 12)
            cache[J] = sample.P.star.n(theta.n[.i], acdpi.struct$c.n[.s],
                                       acdpi.struct$d.n[.s], .int, alpha);##
            .b = which(cache[J]> acdpi.struct$c.n[.s] &
                       cache[J]<=acdpi.struct$d.n[.s]);
            if (length(.b)>0) {
                .t = cache[J[.b]];
                ## Calculate \theta=\lambda_z^-1(t) and \lambda'_z(\theta)
                ## in Eq. (6.19).
                .x = log.z[J[.b]] + lk.theta.n[.i[.b]] + .t + log1p(-exp(-.t));
                .inv = pi - (pi*delta/alpha)/expm1(-delta*.x);
                .deriv = pi*expm1(-.t)/(pi-.inv)/(pi-alpha*.inv);

                val[J[.b]] = .inv;
                .P =  P.star.n(.t, theta.n[.i[.b]],
                               acdpi.struct$a.n[.s[.b]],
                               acdpi.struct$c.n[.s[.b]],
                               acdpi.struct$d.n[.s[.b]],
                               acdpi.struct$pi.n[.s[.b]],
                               alpha)$val * abs(.deriv);
                r[J[.b]] = U[J[.b]]*(1+Delta)*.P/
                    Q.alpha.z(val[J[.b]], z[J[.b]], alpha);
            }
        }
        I = I[which(r[I]>1)];
    }

    return(list(val=val, intvl=S, n.types=c(1, length(t.n)-1,
                                            length(theta.n)-1)));
} ## sample.Q.alpha.z

sample.chi.alpha.z.5.1 <- function(z, alpha, alpha0, theta0, Delta)
### Algorithm 5.1 to sample from the normalized \chi_{\alpha,z}.  Same
### requirements on the input arguments as in sample.Q.alpha.z().
{
    cat("Entering sample.chi.alpha.z.5.1...\n");
    nz = length(z);
    if (nz==0) { return(numeric(0)); }
    if (any(z>=1 | z<=0)) {
        stop("All entries of _z_ must be in (0,1)");
    }

    delta = 1-alpha;
    c.a = c.alpha(alpha);
    
    chi.P.ratio <- function(v, tau.y, b)
    {
        ## Evaluate the r.h.s. of Eq. (7.1).  The input arguments v, tau.y,
        ## and b must be of the same length; tau.y should consist of values
        ## of $\tau*y$, with $y = \exp(v)-1$, and b should consist of values
        ### of $\ell(1/\tau)$
        if (length(v) != length(tau.y) | length(v) != length(b) ) {
            stop("_v_, _tau.y_, _b_ must have the same length");
        }
        
        R = phi.star.abc(v, 0, b, delta)$val.s;  ## R(v) in (8.1)
        D = exp(tau.y-v)*R + v^alpha/c.a; ## Denominator of the rhs
        I = which(v>=b);
        if (length(I)>0) {
            D[I] = D[I] + (v[I]/b[I])^alpha;
        }

        return(exp(-alpha)*(x.over.expm1(-delta/alpha*v))^alpha/D);        
    }
    
    y = rep(0, nz);
    log.y=log(y); ## 1/16
    theta = rep(0, nz);
    accept = rep(FALSE, nz);  ## acceptance/rejection
    log.z = log(z);
    
    I = 1:nz; ## Entries not yet sampled
    S = rep(0, nz);
    n.types=c();
    
    while (length(I)>0) {
        Q.struct=sample.Q.alpha.z(z[I], alpha, alpha0, theta0, Delta);

        theta[I]=Q.struct$val;

        tau = exp(log.z[I] + log.H.alpha(theta[I], alpha));

        t.i = which(tau<Inf);
        if (length(t.i)>0) {
            v = rep(0, length(t.i)); log.v=log(v);
            tau = tau[t.i];

            ell = log1p.recip(tau);
            elld = ell^delta; elldt = elld/(tau*ell);
            ## $I_i$, $i=0,1,2$ in Eqs. (5.7)-(5.8)
            I.i =rbind(4*elldt +(2/delta-4)*elld, elldt/exp(1), 1/(c.a*tau));
            
            ix = sample.mat(I.i)-1;
            sub = which(ix==0);  ## Case i=0 in Algorithm 5.1
            if (length(sub)>0) {
                .phi = sample.phi.star.abc(0, ell[sub], delta);
                v[sub] = .phi$val;
                log.v[sub] = .phi$log.val;
            }
            sub = which(ix>0); ## Case i>0 in Algorithm 5.1
            if (length(sub)>0) {
                v[sub] = log1p((rexp(length(sub))+2-ix[sub])/tau[sub]);
                log.v[sub]=log(v[sub]);
            }

            sub = which(log.v>-Inf);
            if (length(sub)>0) {
                J = I[t.i[sub]];
                y[J] = expm1(v[sub]);
                log.y[J]=log.v[sub] - log(x.over.expm1(v[sub]));
                cpr=chi.P.ratio(v[sub], tau[sub]*y[J], ell[sub]);
                accept[J] = (runif(length(sub)) <=cpr);
            }
        }
        I = I[which(!accept[I])];
    }
    one.minus.x = -expm1(-delta/alpha*log1p(y)); ## one minus undershoot
    
    return(list(y=y, log.y=log.y, theta=theta, undershoot.gap=one.minus.x,
                undershoot=1-one.minus.x, alpha=alpha));
} # sample.chi.alpha.z.5.1

compare.sample.algs <- function(z, alpha, alg=1:3)
### Given an array of z's, compare the cpu times of the three routines to
### sample y's and run KS test to check if the samples have different
### distributions
{
    alpha0=2/3; theta0=6*pi/7; Delta=.5;

    a = c();
    for ( i in 1:3) {
        if (any(alg==i)) { a = c(a,i);}
    }
    if (length(alg)==0) {
        cat("No algorithm specified.  Implement all 3 algorithms.\n");
        alg = 1:3;
    } else {
        alg = a;
    }
    na = length(alg);
    
    Y=c(); Log.Y=c(); Theta=c();
    S = c();
    L=c();
    if (any(alg==1)) {
        S1=system.time(X1 <- sample.chi.alpha.z.3.1(z,alpha));
        L=c(L, "1");
        S = rbind(S,S1);
        Y = rbind(Y, X1$y);
        Log.Y = rbind(Log.Y, X1$log.y);
        Theta = rbind(Theta, X1$theta);
    }
    
    if (any(alg==2)) {
        S2=system.time(X2 <- sample.chi.alpha.z.3.2(z,alpha));
        L=c(L, "2");
        S = rbind(S,S2);
        Y = rbind(Y, X2$y);
        Log.Y = rbind(Log.Y, X2$log.y);
        Theta = rbind(Theta, X2$theta);
    }

    if (any(alg==3)) {
        S3=system.time(X3 <- sample.chi.alpha.z.5.1(z,alpha,alpha0,theta0,
                                                         Delta));
        L=c(L, "3");
        S = rbind(S,S3);
        Y = rbind(Y, X3$y);
        Log.Y = rbind(Log.Y, X3$log.y);
        Theta = rbind(Theta, X3$theta);
    }

    rownames(S) = L;
    if (na<=1) {
        KS.y=numeric(0);
        KS.log.y=numeric(0);
        KS.theta=numeric(0);
    } else {
        KL=c();
        KS.y = matrix(rep(0,na*(na-1)), nc=2);
        cnt=0;
        for ( i in 1:(na-1) ) {
            for ( j in (i+1):na ) {
                cnt = cnt+1;
                kst=ks.test(Y[i,],Y[j,]);
                KS.y[cnt,]=c(kst$statistic, kst$p.value);
                KL = c(KL, sprintf("%dvs%d", alg[i], alg[j]));
            }
        } 
        colnames(KS.y) = c("D", "p-value");
        rownames(KS.y) = KL;

        KL=c();
        KS.log.y = matrix(rep(0,na*(na-1)), nc=2);
        cnt=0;
        for ( i in 1:(na-1) ) {
            for ( j in (i+1):na ) {
                cnt = cnt+1;
                kst=ks.test(Log.Y[i,],Log.Y[j,]);
                KS.log.y[cnt,]=c(kst$statistic, kst$p.value);
                KL = c(KL, sprintf("%dvs%d", alg[i], alg[j]));
            }
        } 
        colnames(KS.log.y) = c("D", "p-value");
        rownames(KS.log.y) = KL;

        KL=c();
        KS.theta = matrix(rep(0,na*(na-1)), nc=2);
        cnt=0;
        for ( i in 1:(na-1) ) {
            for ( j in (i+1):na ) {
                cnt = cnt+1;
                kst=ks.test(Theta[i,],Theta[j,]);
                KS.theta[cnt,]=c(kst$statistic, kst$p.value);
                KL = c(KL, sprintf("%dvs%d", alg[i], alg[j]));
            }
        } 
        colnames(KS.theta) = c("D", "p-value");
        rownames(KS.theta) = KL;
    }
    return(list(alpha=alpha, alg=alg, CPU=S, KS.y=KS.y, KS.log.y=KS.log.y, 
                KS.theta=KS.theta));
} ## compare.sample.algs

sample.chi.alpha.z <- function(z, alpha, alpha0=2/3, theta0=6*pi/7, Delta=.5)
### Sample from the normalized $\chi_{\alpha,z}$ by calling one of the
###             sample.chi.alpha.z.?.?()
### functions depending on the value of $(z,\alpha)$.  The input argument z
### must be an array of positive numbers, and alpha must be a scaler in (0,1).
{
    z0=1; z1=(1-alpha)*10^-30; alpha1=max(alpha0, 0.9);

    cat("Entering sample.chi.alpha.z ...\n");
    
    nz = length(z);
    if (nz==0) return(numeric(0));
    if (length(alpha)!=1 | alpha<=0 | alpha>=1) {
        stop(bold(red("alpha")), " must be a scalar in (0,1)");
    }
    if (any(z<=0)) {
        stop("All the entries of ", bold(red("z")), " must be positive");
    }

    y=rep(0,nz); log.y = rep(0,nz); theta = rep(0,nz);
    
    ## Call sample.chi.alpha.z.3.1() for entries in _z_ that are >= _z0_,
    ## regardless of the value of _alpha_
    I = which(z>=z0);
    if (length(I)>0) {
        S = sample.chi.alpha.z.3.1(z[I], alpha);
        y[I]=S$y;
        log.y[I]=S$log.y;
        theta[I]=S$theta;
    }
        
    ## If _alpha_<= _alpha1_, call sample.chi.alpha.z.3.2() for entries in
    ## _z_ that are < _z0_, otherwise, call the same function for entries
    ## in _z_ in [_z1_, _z0_]. 
    I = which(z<z0);
    if (alpha>alpha1) {
        I = I[which(z[I]>=z1)];
    }
    if (length(I)>0) {
        S = sample.chi.alpha.z.3.2(z[I], alpha);
        y[I]=S$y;
        log.y[I]=S$log.y;
        theta[I]=S$theta;
    }

    ## If _alpha_ > _alpha1_, call sample.chi.alpha.z.5.1() for entries in
    ## _z_ that are < _z1_
    I=which(z<z1);
    if (alpha>alpha1 & length(I)>0) {
        S = sample.chi.alpha.z.5.1(z[I], alpha,alpha0,theta0,Delta);
        y[I]=S$y;
        log.y[I]=S$log.y;
        theta[I]=S$theta;
    }
    
    return(list(log.y=log.y, y=y, theta=theta));
} # sample.chi.alpha.z

sample.fp <- function(n, alpha, b, diff.b, linv.B, ...)
### Algorithm 7.1.  b is the boundary function, diff.b is its derivative,
### and linv.B is the log of the inverse of B(t) = t^{-1/alpha} * b(t).
### For example, for constant boundary at M:
###   f <- function(t,a,M) { return(rep(M,length(t))); }
###   df <- function(t,a,M) { return(rep(0, length(t))); }
###   linv.F <- function(s,a,M) { a * (log(M) - log(s));}
###   sample.fp(1000, .7, f, df, linv.F, 5);
{
    z = rep(0,n);  log.z=rep(-Inf, n); 
    y = rep(0,n);  log.y=rep(-Inf, n); theta=rep(NA,n);
    
    delta = 1-alpha;

    ## Sample _z_ 
    I = 1:n;
    while(length(I)>0) {
        log.z[I] = sample.z(length(I), alpha, log=TRUE);
        z[I] = exp(log.z[I]);
        I = I[which(z[I]<=0)];
    }
    
    log.s = log(alpha) + delta/alpha*(log(delta) - log.z);
    log.t = linv.B(exp(log.s), alpha, ...);
    t = exp(log.t);
    
    db = -diff.b(t, alpha, ...);
    I = which(runif(n)>=db/(db + b(t,alpha,...)/alpha/t));

    if (length(I)>0) {
        S=sample.chi.alpha.z(z[I], alpha);
        y[I] = S$y; log.y[I] = S$log.y; theta[I] = S$theta;
    }

    return(list(log.z=log.z, z=z, log.t=log.t, t=t, log.y=log.y,
                y=y, theta=theta));    
} # sample.fp
