## $Id: ar-fp-runs.R,v 1.1 2025/01/04 06:07:15 zchi Exp zchi $

source("ar-fp.R");
library("ggplot2");

batch.compare.sample.algs <- function(zs, alphas, size.per.rep, n.rep,
                                      fname=c(), algs=1:3)
### Applies compare.sample.algs to a list of z values and alpha values.
### The output is a list of length equal to
###          length(zs) x length(alphas) x n.rep
### The result for zs[i] is in the i-th block of length 
###          length(alphas) x n.rep
### Within the block, the result for zs[i] and alphas[j] is in the j-th
### sub-block of length _n.rep_, so that the k-th rep is stored in the k-th
### entry of the sub-block.  Each entry itself is a list.  
{
    K = c();
    cnt=1;
    for ( z in zs) {
        for ( alpha in rep(alphas, each=n.rep) ) {
            cat('z=', z, ' alpha=', alpha, '...\n');
            K[[cnt]] = compare.sample.algs(rep(z, size.per.rep), alpha,
                                           alg=algs);
            K[[cnt]]$z=z;
            K[[cnt]]$size=size.per.rep;
            cnt=cnt+1;
        }
    }
    if ( is.character(fname)) {
        Batch.Result=K;
        save(zs, alphas, algs, n.rep, Batch.Result, file=fname);
        cat(sprintf("Output saved in '%s'\n", fname));
    } else {
        return(K);
    }    
}

check.batch.compare <- function(fname, print.fig=FALSE)
### Summarizes the output of a call to batch.compare.sample.algs() that
### is stored in the file with name fname.  If print.fig=TRUE, then all
### figures displayed are saved as PDF files.
{
    load(batch.compare.sample.algs.output.fname);
    n.z=length(zs); n.a = length(alphas);
    run.times=list();
    p.vals=list();
    hc.ks = list();
    for ( i in 1:n.z ) {
        rt = c();
        pv = c();
        for ( j in 1:(n.a*n.rep)) {
            k = (i-1)*n.a*n.rep + j;
            X = Batch.Result[[k]];
            rt = cbind(rt, X$CPU[,"user.self"]);

            ## _pv_ consists of _n.a_ blocks of columns, the l-th block
            ## for the l-th entry in _alphas_.  The l-th block consists of
            ## _n.rep_ columns, one column per rep, and consists of p-values
            ## of KS test on the sampled values from a pair of algorithms:
            ## A vs B, A<B.  The entries in the column are arranged according
            ## to the dictionary order of (A,B).
            pv=cbind(pv, X$KS.log.y[,"p-value"]); 
        }
        rt.mean = t(matrix(apply(matrix(t(rt), nr=n.rep), 2, mean), nr=n.a));
        if (n.rep>1) {
            rt.std = t(matrix(apply(matrix(t(rt), nr=n.rep), 2, std), nr=n.a));
        } else {
            rt.std = t(matrix(apply(matrix(t(rt*0), nr=n.rep), 2, std),nr=n.a));
        }
        run.times[[i]] = rt.mean;
        rt.mean = as.vector(t(rt.mean));
        rt.std = as.vector(t(rt.std));
        alg.lab=rep(algs, each=n.a);
        df = data.frame(alpha=alphas, rt.mean=rt.mean, rt.std=rt.std,
                        algorithm=factor(alg.lab));

        figure(1);
        ## print(ggplot(data=df, aes(x=alpha, y=rt.mean, group=algorithm,
        ##                           shape=algorithm))+ geom_point(size=2) + 
        ##       geom_line(color='gray', show.legend=TRUE) +
        ##       geom_errorbar(aes(ymin=rt.mean-rt.std, ymax=rt.mean+rt.std),
        ##                     width=0)+
        ##       theme_bw()+ylab("run time")+ scale_shape(solid=FALSE)+
        ##       ggtitle(sprintf("z=%E", zs[i])));
        print(ggplot(data=df, aes(x=alpha, y=rt.mean, group=algorithm,
                                  shape=algorithm))+
              geom_point(size=2, show.legend=FALSE) + 
              geom_line(color='gray', show.legend=FALSE) +
              geom_errorbar(aes(ymin=rt.mean-rt.std, ymax=rt.mean+rt.std),
                            width=0)+
              theme(axis.text = element_text(color="black", size=15))+
              expand_limits(y=0)+
              labs(x=NULL, y=NULL)+ scale_shape(solid=FALSE));
        if (print.fig) {
            sfname=sprintf("%sz=%E-runtime.pdf",
                           batch.compare.sample.algs.output.fname, zs[i]);
            ggsave(sfname);
            cat("Plot saved in ", sfname, "\n");
        }

        if (length(algs>1)) {
            alg.pairs=c();
            for ( a in 1:(length(algs)-1) ) {
                for ( b in (a+1):length(algs) ) {
                    alg.pairs=c(alg.pairs, sprintf("%d vs %d",algs[a],algs[b]));
                }
            }
            pv.mean = t(matrix(apply(matrix(t(pv), nr=n.rep), 2, mean),
                               nr=n.a));
            pv.std = t(matrix(apply(matrix(t(pv), nr=n.rep), 2, std),
                               nr=n.a));
            p.vals[[i]] = pv.mean;
            pv.mean = as.vector(t(pv.mean));
            pv.std = as.vector(t(pv.std));
            df2=data.frame(alpha=alphas, pv.mean = pv.mean, pv.std=pv.std,
                           alg.pairs=rep(alg.pairs,each=n.a));

            figure(2);
            ## print(ggplot(data=df2,aes(x=alpha, y=pv.mean, group=alg.pairs,
            ##                           shape=alg.pairs))+ geom_point(size=2)+
            ##       geom_line(color='gray', show.legend=TRUE) +
            ##       theme_bw() + ylab("p-value") + scale_shape(solid=FALSE)+
            ##       ggtitle(sprintf("z=%E", zs[i])));
            print(ggplot(data=df2,aes(x=alpha, y=pv.mean, group=alg.pairs,
                                      shape=alg.pairs))+
                  geom_point(size=2, show.legend=FALSE)+
                  geom_line(color='gray', show.legend=FALSE) +
                  geom_errorbar(aes(ymin=pv.mean-pv.std, ymax=pv.mean+pv.std),
                                width=0)+
                  theme(axis.text = element_text(color="black", size=15)) +
                  labs(x=NULL, y=NULL) + scale_shape(solid=FALSE) +
                  scale_y_continuous(limits=c(0,1))+
                  geom_hline(yintercept=0.05)
                  );

            if (print.fig) {
                sfname=sprintf("%sz=%E-pvalue.pdf",
                               batch.compare.sample.algs.output.fname, zs[i]);
                ggsave(sfname);
                cat("Plot saved in ", sfname, "\n");
            }
        }

        uicontinue();
    }
    return(list(run.times=run.times, p.vals=p.vals, hc.ks=hc.ks));
}

batch.sample.fp <- function(alphas, size.per.rep, n.rep, fname=c())
### Applies sample.fp to a list of alpha values.  The output is a list of
### length equal to
###          length(alphas) x n.rep
### The result for alphas[i] is in the i-th block of length _n.rep_,
### so that the k-th rep is stored in in the k-th entry of the block.
### Each entry itself is a list.    
{
    K = c();
    cnt=1;
    for ( alpha in rep(alphas, each=n.rep) ) {
        cat('alpha=', alpha, '...\n');
        K[[cnt]] = list(CPU=system.time(sample.fp(size.per.rep, alpha)),
                        alpha=alpha, size=size.per.rep);
        cnt=cnt+1;
    }
    
    if ( is.character(fname)) {
        Batch.Result=K;
        save(alphas, n.rep, Batch.Result, file=fname);
        cat(sprintf("Output saved in '%s'\n", fname));
    } else {
        return(K);
    }    
}

check.batch.sample.fp <- function(fname, print.fig=FALSE)
### Summarizes the output of a call to batch.sample.fp() that is stored
### in the file with name fname.  If print.fig=TRUE, then all figures
### displayed are saved as PDF files.
{    
    load(batch.sample.fp.output.fname);
    N = length(Batch.Result);
    alphas=c(); cputimes=c();
    for ( i in 1:N ) {
        alphas = c(alphas, Batch.Result[[i]]$alpha);
        cputimes = c(cputimes, Batch.Result[[i]]$CPU[3]);
    }
    na = length(unique(alphas));
    n.rep = N/na;
    alphas = alphas[seq(1,N, by=n.rep)];
    cputimes = matrix(cputimes, nrow = n.rep);
    cpu.mean = apply(cputimes, 2, mean);
    cpu.std = apply(cputimes, 2, std);
    size.per.rep = Batch.Result[[1]]$size;
    df = data.frame(alphas=alphas, cpu.mean=cpu.mean, cpu.std=cpu.std);
    print(ggplot(data=df, aes(x=alphas, y=cpu.mean))  +
          geom_point(size=2) + geom_line(color='gray')+
          geom_errorbar(aes(ymin=cpu.mean-cpu.std, ymax=cpu.mean+cpu.std))+
          theme(axis.text = element_text(color="black", size=15))+
          labs(x=NULL, y=NULL)+ scale_shape(solid=FALSE));
    if (print.fig) {
        sfname=sprintf("%s.pdf", batch.sample.fp.output.fname);
        ggsave(sfname);
        cat("Plot saved in ", sfname, "\n");
    }

    return(list(alphas=alphas, cpu.mean=cpu.mean, cpu.std=cpu.std));
}


record.03.13b <- function() {
    alphas=c(.995, .996, 0.997, .998, .999, .9999);
    size.per.rep=100000;
    n.rep=100;
    fname=sprintf("sample.fp.batch_%sb", substring(Sys.time(),1,10));
    batch.sample.fp(alphas, size.per.rep, n.rep, fname);
}

record.03.13 <- function() { ### 3/13/2025
    alphas=c(seq(.005, .995, by=.03), .996, 0.997, .998, .999, .9999);
    size.per.rep=10000;
    n.rep=100;
    fname=sprintf("sample.fp.batch_%s", substring(Sys.time(),1,10));
    batch.sample.fp(alphas, size.per.rep, n.rep, fname);
}

record.03.03 <- function() {  ### 3/3/2025
    zs = c(10^seq(-4, -1, by=1), 2, 3, 4)
    alphas = c(seq(.905, .995, by=.03), .999);
    size.per.rep=200000;
    n.rep=1;
    fname=sprintf("compare.sample.algs.batch_%s",substring(Sys.time(),1,10));
    alg=1:2;
    batch.compare.sample.algs(zs, alphas, size.per.rep, n.rep, fname, alg);
}


record.02.24 <- function() { ### 2/24/2025
    zs=10^seq(-60, -5, by=5);
    alphas=seq(.995, .9995, by=.0005);
    size.per.rep=100000;
    n.rep=1;
    fname=sprintf("compare.sample.algs.batch_%s",substring(Sys.time(),1,10));
    alg=2:3;
    batch.compare.sample.algs(zs, alphas, size.per.rep, n.rep, fname, alg);
}

record.02.23 <- function() {  ### 2/23/2025
    zs=10^seq(log(0.05, base=10), log(4, base=10), by=.15);
    alphas=seq(.005, .995, by=.03);
    size.per.rep=200000;
    n.rep=1;
    fname=sprintf("compare.sample.algs.batch_%s",substring(Sys.time(),1,10));
    alg=1:2;
    batch.compare.sample.algs(zs, alphas, size.per.rep, n.rep, fname, alg);
}

record.02.20 <- function() { ### run on 2/20/2025
    zs=10^seq(log(0.05, base=10), log(4, base=10), by=.15);
    alphas=seq(.005, .995, by=.03);
    size.per.rep=1000;
    n.rep=200;
    fname=sprintf("compare.sample.algs.batch_%s",substring(Sys.time(),1,10));
    alg=1:2;
    batch.compare.sample.algs(zs, alphas, size.per.rep, n.rep, fname, alg);
}
