# Simulating integrated Brownian motion

The state is ordered as $`(f',f)`$. For interval length $`\Delta`$ and
derivative diffusion variance rate $`\lambda`$, its innovation
covariance is

``` math
\lambda
\begin{pmatrix}
\Delta & \Delta^2/2\\
\Delta^2/2 & \Delta^3/3
\end{pmatrix}.
```

The analytic Cholesky factor used by Stan can be checked directly.

``` r

delta <- 0.8
lambda <- 0.4
L <- ibm_transition_cholesky(delta, lambda)
Q <- ibm_transition_covariance(delta, lambda)
stopifnot(isTRUE(all.equal(L %*% t(L), Q)))
Q
```

    ##       [,1]       [,2]
    ## [1,] 0.320 0.12800000
    ## [2,] 0.128 0.06826667

Constant `lambda` gives ordinary IBM. Supplying one value per interval
gives the locally varying diffusion used by the adaptive model.

``` r

t <- c(0, 0.1, 0.7, 2.2, 3)
ordinary <- simulate_ibm(t, lambda = 0.3, seed = 1)
adaptive <- simulate_ibm(t, lambda = c(0.02, 0.02, 1.2, 0.02), seed = 2)

plot(t, adaptive$f, type = "l", xlab = "t", ylab = "f(t)")
```

![](simulating-ibm_files/figure-html/simulate-1.png)
