## ---- echo=FALSE, cache=FALSE-------------------------------------------------
suppressPackageStartupMessages(library("sparseMVN"))
suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("mvtnorm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("forcats"))
suppressPackageStartupMessages(library("kableExtra"))

knitr::opts_chunk$set(prompt=TRUE, cache=FALSE,error=FALSE,
                      echo=FALSE, eval=TRUE,
                      comment="#", collapse=FALSE) #$
options(replace.assign=TRUE,
        digits=4, dplyr.summarise.inform=FALSE)

## ---- echo=FALSE--------------------------------------------------------------
N <- 5
k <- 2
p <- k ## dimension of mu
nv1 <- N*k+p
nels1 <- nv1^2
nnz1 <- N*k^2 + 2*p*N*k + p^2
nnz1LT <- N*k*(k+1)/2  + p*N*k + p*(p+1)/2
Q <- 1000
nv2 <- Q*k+p
nels2 <- nv2^2
nnz2 <- Q*k^2 + 2*p*Q*k + p^2
nnz2LT <- Q*k*(k+1)/2 + p*Q*k + p*(p+1)/2
options(scipen=999)

## ----blockarrow, results="asis"-----------------------------------------------
Mat1a <- as(kronecker(diag(N), matrix(1, k, k)),"sparseMatrix")
Mat1a <- rbind(Mat1a, Matrix(1, p, N*k))
Mat1a <- cbind(Mat1a, Matrix(1, k*N+p, p))

tmp_a <- as.matrix(Mat1a)
colnames(tmp_a) <- letters[1:NCOL(tmp_a)]


tmp_a %>%
  as_tibble(.name_repair="minimal") %>%
  mutate(across(colnames(tmp_a), ~ifelse(.x == 1, "|", "."))) %>%
  kbl(col.names=NULL, caption="A block-arrow pattern") %>%
  kable_minimal()


## ---- banded, echo=FALSE,results="asis"---------------------------------------
Mat1b <- kronecker(Matrix(1, k, k), diag(N))
Mat1b <- rbind(Mat1b, Matrix(1, p, N * k))
Mat1b <- cbind(Mat1b, Matrix(1, k*N+p, p))
tmp_b <- as.matrix(Mat1b)
colnames(tmp_b) <- letters[1:NCOL(tmp_b)]


tmp_b  %>%
  as_tibble(.name_repair="minimal") %>%
  mutate(across(colnames(tmp_b), ~ifelse(.x == 1, "|", "."))) %>%
  kbl(col.names=NULL, caption="A banded sparsity pattern") %>%
  kable_minimal()


## ---- echo=FALSE, results="hide"----------------------------------------------
Mat2 <- as(kronecker(diag(Q),matrix(1,k,k)),"lMatrix") %>%
    rbind(Matrix(TRUE,p,Q*k)) %>%
    cbind(Matrix(TRUE, k*Q+p, p)) %>%
    as("dgCMatrix") %>%
    as("symmetricMatrix")
A2 <- as(Mat2,"matrix")

## ---- eval=FALSE, echo=TRUE, prompt=FALSE-------------------------------------
#  rmvn.sparse(n, mu, CH, prec=TRUE, log=TRUE)
#  dmvn.sparse(x, mu, CH, prec=TRUE, log=TRUE)

## ---- results='hide'----------------------------------------------------------
D <- sparseMVN::binary.sim(N=50, k=2, T=50)
priors <- list(inv.A=diag(2), inv.Omega=diag(2))
start <- rep(c(-1,1),51)
opt <- trustOptim::trust.optim(start,
                               fn=sparseMVN::binary.f,
                               gr=sparseMVN::binary.grad,
                               hs=sparseMVN::binary.hess,
                               data=D, priors=priors,
                               method="Sparse",
                               control=list(function.scale.factor=-1))

## -----------------------------------------------------------------------------
R <- 100
pm <- opt[["solution"]]
H <- -opt[["hessian"]]
CH <- Cholesky(H)

## ---- echo=FALSE--------------------------------------------------------------
load("runtimes.Rdata")

## ---- echo=FALSE--------------------------------------------------------------
tab1 <- filter(runtimes, stat %in% c("density","rand")) %>%
  dplyr::group_by(N, k, stat, pattern, type) %>%
  dplyr::summarize(mean_ms=mean(time/1000000)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = c(pattern, type), values_from = mean_ms)

tab2 <- as_tibble(runtimes) %>%
  dplyr::filter(stat %in% c("chol","solve")) %>%
  dplyr::group_by(N, k, stat, pattern, type) %>%
  dplyr::summarize(mean_ms=mean(time/1000000)) %>%
  tidyr::pivot_wider(names_from = c(stat, pattern), values_from = mean_ms)

## ----cases, echo=FALSE, fig.cap="Cases for timing comparison. $N$ and $k$ refer, respectively, to the number of blocks in the block-arrow structure (analogous to heterogeneous units in the  binary choice example), and the size of each block.  The total number of variables is $M=(N+1)k$, and the total number of elements in the matrix is $M^2$. The three rightmost columns present the number of nonzeros in the full matrix and lower triangle, and the sparsity (proportion of matrix elements that are nonzero)."----
cases %>%
  select(k,N, nvars,nels,
         nnz, nnzLT, pct.nnz) %>%
  kable(digits=c(rep(0,6),3), format.args=list(big.mark=","),
        col.names=c('k', 'N', 'variables', 'elements', 'full', 'lower tri', 'sparsity')) %>%
  kable_minimal()  %>%
   add_header_above(c(" "=4, "nonzeros"=3))

## ---- echo=FALSE--------------------------------------------------------------
tmp <- dplyr::filter(tab1, N==min(tab1[['N']]), k==min(tab1[['k']]),
              stat=="density")
sm <- with(tmp, c(dense_cov,sparse_cov))

## ----densRand, echo=FALSE, fig.width=5, fig.cap="Mean computation time for simulating 1,000 MVN samples, and computing 1,000 MVN densities, averaged over 200 replications. Densities were computed using [dmvnorm] and [dmvn.sparse], while random samples were generated using [rmvnorm] and [rmvn.sparse]."----

DF1 <- tab1 %>%
  dplyr::mutate(stat=forcats::fct_recode(stat, random='rand')) %>%
  dplyr::rename(dense=dense_cov, sparse=sparse_cov) %>%
  tidyr::pivot_longer(cols=c(dense, sparse), names_to='pattern', values_to='value')


P1 <- ggplot(DF1, aes(x=N, y=value, color=pattern, shape=pattern, linetype=pattern)) +
  geom_line() +
  geom_point(size=2) +
  scale_x_continuous("Number of blocks (N)") +
  scale_y_continuous("Computation time (milliseconds)", labels=scales::comma) +
  scale_color_manual("Pattern",values=c(dense='red', sparse='blue')) +
  scale_shape("Pattern") +
  scale_linetype("Pattern") +
  facet_grid(stat~k, scales="free_y", labeller=label_bquote(cols = k==.(k))) +
  theme_bw() +
  theme(strip.background=element_rect(fill='white'))

print(P1)


## ----cholSolve,  echo=FALSE, fig.width=5, fig.keep='all', fig.cap="Computation time for Cholesky decomposition of sparse and dense matrices, and inversion of dense matrices."----

DF2 <- tab2 %>%
    dplyr::rename(`dense inversion`=solve_dense,
           `dense Cholesky`=chol_dense, `sparse Cholesky`=chol_sparse) %>%
  tidyr::pivot_longer(cols=c(`dense inversion`,`sparse Cholesky`,`dense Cholesky`),
               names_to = 'pattern', values_to = 'value')


P2 <- ggplot(DF2, aes(x=N, y=value, color=pattern, shape=pattern, linetype=pattern)) +
  geom_line() +
  geom_point(size=2) +
  scale_x_continuous("Number of blocks (N)") +
  scale_y_continuous("Computation time (milliseconds)", labels=scales::comma) +
  facet_grid(.~k, scales="free_y", labeller=label_bquote(cols = k==.(k))) +
  theme_bw() +
  theme(strip.background=element_rect(fill='white'))

print(P2)


