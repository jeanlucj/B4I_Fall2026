# The bivariate DGE-IGE model in matrix-variate form

[B4I_Proposal_Models.docx](B4I_Proposal_Models.docx) gives the joint
producer–associate model with the oat and pea responses **stacked into one long
vector**, each random term pre-multiplied by $I_2 \otimes Z$. This note rewrites
it with the two responses as **two columns of a matrix**. The two forms are
algebraically identical — one identity converts between them — but the
matrix-variate form is shorter, makes the trait-order convention explicit, and
reveals structure in the interaction term that the stacked form hides.

Math renders on GitHub and in any Markdown viewer with LaTeX support.

---

## 1. The model as the proposal writes it

Transcribed from the proposal's equations 31–34. Writing $\mathrm{Pr}$ for the producer
effect and $As$ for the associate effect, with $\to$ reading "effect on":

$$
\begin{pmatrix} y^{oat} \\ y^{pea} \end{pmatrix}
= \left(I_2 \otimes Z_{oat}\right)
  \begin{pmatrix} \mathrm{Pr}_i^{oat} \\ As_i^{oat \to pea} \end{pmatrix}
+ \left(I_2 \otimes Z_{p}\right)
  \begin{pmatrix} As_j^{pea \to oat} \\ \mathrm{Pr}_j^{pea} \end{pmatrix}
+ \left(I_2 \otimes Z_{s}\right)
  \begin{pmatrix} S_{ij}^{oat \times pea \to oat} \\ S_{ij}^{pea \times oat \to pea} \end{pmatrix}
+ \begin{pmatrix} \epsilon_{ijk}^{oat} \\ \epsilon_{ijk}^{pea} \end{pmatrix}
$$

with the three random terms distributed as

$$
\begin{pmatrix} \mathrm{Pr}_i^{oat} \\ As_i^{oat \to pea} \end{pmatrix} \sim
\begin{pmatrix}
  \sigma^2_{\mathrm{Pr},oat} & \sigma_{\mathrm{Pr} As,oat} \\
  \sigma_{\mathrm{Pr} As,oat} & \sigma^2_{As,oat \to pea}
\end{pmatrix} \otimes \mathbf{G}_{oat}
$$

$$
\begin{pmatrix} As_j^{pea \to oat} \\ \mathrm{Pr}_j^{pea} \end{pmatrix} \sim
\begin{pmatrix}
  \sigma^2_{As,pea \to oat} & \sigma_{\mathrm{Pr} As,pea} \\
  \sigma_{\mathrm{Pr} As,pea} & \sigma^2_{\mathrm{Pr},pea}
\end{pmatrix} \otimes \mathbf{G}_{pea}
$$

$$
\begin{pmatrix} S_{ij}^{oat \times pea \to oat} \\ S_{ij}^{pea \times oat \to pea} \end{pmatrix} \sim
\begin{pmatrix}
  \sigma^2_{s,oat \times pea \to oat} & \sigma_{sop} \\
  \sigma_{sop} & \sigma^2_{s,pea \times oat \to pea}
\end{pmatrix} \otimes \mathbf{G}_{oat} \otimes \mathbf{G}_{pea}
$$

> **Two transcription notes, not corrections.** The proposal's equation 31 uses
> $\sim$ where $=$ belongs — the left side is the data, not a random variable
> with that distribution. And it omits the fixed term $X\beta$ that equations 1
> and 21 carry; the generic $X\beta + Zu$ for the experimental design is
> introduced earlier and taken as understood. Both are kept out of the way
> below by writing the fixed part explicitly.
>
> Equations 32–34 likewise write $\sim$ directly in front of a covariance
> matrix, leaving the $\mathcal{N}(0,\,\cdot)$ implicit; equations 25–27, the
> single-response versions, write it out in full. Reproduced as given above and
> written out properly in §4.
>
> The residual covariance is left open in the proposal — "there may also be
> reason and possibility to estimate error covariance components between the
> vectors $\epsilon^{oat}_{ijk}$ and $\epsilon^{pea}_{ijk}$". The
> implementation does estimate it, unstructured; see §7.

---

## 2. The one identity that converts the two forms

Everything rests on this. For any matrix $U$ with $p$ columns and any conformable
$Z$:

$$
\left(I_p \otimes Z\right)\,\mathrm{vec}(U) \;=\; \mathrm{vec}(Z\,U)
$$

**What $\mathrm{vec}$ does.** $\mathrm{vec}(U)$ stacks the *columns* of $U$ on
top of each other into one long vector. So if
$U_{oat} = \left[\,\mathrm{Pr}^{oat} \;\middle|\; As^{oat \to pea}\,\right]$ is
$n_{oat} \times 2$, then

$$
\mathrm{vec}(U_{oat}) = \begin{pmatrix} \mathrm{Pr}^{oat} \\ As^{oat \to pea} \end{pmatrix}
$$

which is *exactly* the stacked vector in the proposal's equation 31.

**Why the identity holds.** $I_p \otimes Z$ is block-diagonal with $Z$ repeated
$p$ times, so it applies the same $Z$ separately to each stacked block — which
is the same thing as applying $Z$ to each column of $U$ at once. Nothing is
assumed; it is a restatement.

So **every $\left(I_2 \otimes Z\right)$ in the proposal is a $Z\,U$ waiting to
be written as one.** That is the whole conversion.

---

## 3. The matrix normal distribution

The same move on the variance side needs one definition.

A random matrix $U$ of size $n \times p$ is **matrix normal**, written
$U \sim \mathcal{MN}_{n \times p}(M, A, B)$, when

$$
\mathrm{vec}(U) \sim \mathcal{N}\!\left(\mathrm{vec}(M),\; B \otimes A\right)
$$

- $M$ is the $n \times p$ mean.
- $A$ is $n \times n$: the covariance **among rows**. Here, among accessions —
  this is where a genomic relationship matrix goes.
- $B$ is $p \times p$: the covariance **among columns**. Here, among traits —
  this is where the producer–associate covariance goes.

The structure says the row and column covariances are **separable**: the
covariance between entries $(i,t)$ and $(i',t')$ factorises as
$A_{ii'} \times B_{tt'}$. Two accessions' effects covary through kinship, two
traits' effects covary through the $2\times2$ matrix, and the two channels
multiply.

Note the ordering: the column covariance comes **first** in the Kronecker
product, because $\mathrm{vec}$ stacks columns. That is why the proposal writes
its $2\times2$ matrices on the left of $\otimes \mathbf{G}$ — the proposal's
three variance statements are already in exactly this form, which is the clue
that the whole model is matrix-variate.

---

## 4. The model in matrix-variate form

Let $n$ be the number of plots and collect the two responses as the columns of
an $n \times 2$ matrix:

$$
\mathbf{Y} = \left[\, y^{oat} \;\middle|\; y^{pea} \,\right]
$$

Then the entire model is

$$
\boxed{\;
\mathbf{Y} \;=\; \mathbf{X}\mathbf{B}
\;+\; \mathbf{Z}_{oat}\,\mathbf{U}_{oat}
\;+\; \mathbf{Z}_{pea}\,\mathbf{U}_{pea}
\;+\; \mathbf{Z}_{s}\,\mathbf{U}_{s}
\;+\; \mathbf{E} \;}
$$

with each random term matrix normal:

$$
\begin{aligned}
\mathbf{U}_{oat} &\sim \mathcal{MN}\!\left(0,\; \mathbf{G}_{oat},\; \boldsymbol{\Sigma}_{oat}\right) \\[2pt]
\mathbf{U}_{pea} &\sim \mathcal{MN}\!\left(0,\; \mathbf{G}_{pea},\; \boldsymbol{\Sigma}_{pea}\right) \\[2pt]
\mathbf{U}_{s}   &\sim \mathcal{MN}\!\left(0,\; \mathbf{G}_{oat} \otimes \mathbf{G}_{pea},\; \boldsymbol{\Sigma}_{s}\right) \\[2pt]
\mathbf{E}       &\sim \mathcal{MN}\!\left(0,\; \mathbf{I}_n,\; \boldsymbol{\Sigma}_{\epsilon}\right)
\end{aligned}
$$

where the effect matrices hold the two roles of each species side by side,

$$
\mathbf{U}_{oat} = \left[\, \mathrm{Pr}^{oat} \;\middle|\; As^{oat \to pea} \,\right],\qquad
\mathbf{U}_{pea} = \left[\, As^{pea \to oat} \;\middle|\; \mathrm{Pr}^{pea} \,\right],
$$

$$
\mathbf{U}_{s} = \left[\, S^{oat \times pea \to oat} \;\middle|\; S^{pea \times oat \to pea} \,\right],
$$

and the $2\times2$ column covariances are the proposal's, unchanged:

$$
\boldsymbol{\Sigma}_{oat} =
\begin{pmatrix}
  \sigma^2_{\mathrm{Pr},oat} & \sigma_{\mathrm{Pr} As,oat} \\
  \sigma_{\mathrm{Pr} As,oat} & \sigma^2_{As,oat \to pea}
\end{pmatrix},\qquad
\boldsymbol{\Sigma}_{pea} =
\begin{pmatrix}
  \sigma^2_{As,pea \to oat} & \sigma_{\mathrm{Pr} As,pea} \\
  \sigma_{\mathrm{Pr} As,pea} & \sigma^2_{\mathrm{Pr},pea}
\end{pmatrix},
$$

$$
\boldsymbol{\Sigma}_{s} =
\begin{pmatrix}
  \sigma^2_{s,oat \times pea \to oat} & \sigma_{sop} \\
  \sigma_{sop} & \sigma^2_{s,pea \times oat \to pea}
\end{pmatrix},\qquad
\boldsymbol{\Sigma}_{\epsilon} =
\begin{pmatrix}
  \sigma^2_{\epsilon,oat} & \sigma_{\epsilon,oat\,pea} \\
  \sigma_{\epsilon,oat\,pea} & \sigma^2_{\epsilon,pea}
\end{pmatrix}.
$$

**Reading the columns.** Column 1 of every term contributes to oat yield, column
2 to pea yield. So $\mathbf{U}_{oat}$'s first column is the oat's effect on
itself (producer) and its second is the oat's effect on its partner (associate);
$\mathbf{U}_{pea}$ is the mirror, which is why its columns appear in the other
order. **The column position *is* the role.** In the stacked form that
information lives in the order of the sub-vectors, where it is easy to lose.

### Dimensions

| object | size | what one row is |
|---|---|---|
| $\mathbf{Y}$ | $n \times 2$ | a plot: its oat yield and its pea yield |
| $\mathbf{X}$, $\mathbf{B}$ | $n \times p$, $p \times 2$ | design and fixed effects, per trait |
| $\mathbf{Z}_{oat}$, $\mathbf{U}_{oat}$ | $n \times n_{oat}$, $n_{oat} \times 2$ | an oat accession: its producer and associate effects |
| $\mathbf{Z}_{pea}$, $\mathbf{U}_{pea}$ | $n \times n_{pea}$, $n_{pea} \times 2$ | a pea accession, same |
| $\mathbf{Z}_{s}$, $\mathbf{U}_{s}$ | $n \times n_{oat}n_{pea}$, $n_{oat}n_{pea} \times 2$ | an oat × pea **combination** |
| $\mathbf{E}$ | $n \times 2$ | a plot's two residuals |
| $\boldsymbol{\Sigma}_\bullet$ | $2 \times 2$ | among-trait covariance |

$\mathbf{Z}_{oat}$ and $\mathbf{Z}_{pea}$ are plot-to-accession incidence
matrices; $\mathbf{Z}_{s}$ is plot-to-combination. That last one is where the
dimensionality problem lives: its $n_{oat}n_{pea}$ columns are 186,966 for the
real panel, which is why
[interaction_decomposition.md](interaction_decomposition.md) §1.3–1.5 replaces
it with a low-rank basis.

### Checking the equivalence, term by term

Apply $\mathrm{vec}$ to the boxed equation. Using the identity from §2 on each
term, and $\mathrm{vec}(\mathbf{Y}) = (y^{oat\prime}, y^{pea\prime})'$:

$$
\mathrm{vec}(\mathbf{Z}_{oat}\mathbf{U}_{oat})
= \left(I_2 \otimes \mathbf{Z}_{oat}\right)\mathrm{vec}(\mathbf{U}_{oat})
= \left(I_2 \otimes \mathbf{Z}_{oat}\right)\begin{pmatrix} \mathrm{Pr}^{oat} \\ As^{oat \to pea} \end{pmatrix}
$$

which is the proposal's first term exactly, and likewise for the other two. On
the variance side, $\mathbf{U}_{oat} \sim \mathcal{MN}(0, \mathbf{G}_{oat},
\boldsymbol{\Sigma}_{oat})$ means
$\mathrm{vec}(\mathbf{U}_{oat}) \sim \mathcal{N}(0, \boldsymbol{\Sigma}_{oat}
\otimes \mathbf{G}_{oat})$, which is the proposal's variance statement exactly.
The two forms are the same model.

---

## 5. What the rewrite makes visible: the interaction is a three-way array

This is the part the stacked form hides, and it is not cosmetic.

$\mathbf{U}_{s}$ has $n_{oat}n_{pea}$ rows, indexed by a **pair** $(i,j)$. So
each of its two columns is itself a matrix waiting to be unfolded. Write
$\mathbf{I}^{oat}$ for the $n_{oat} \times n_{pea}$ matrix obtained by reshaping
column 1, and $\mathbf{I}^{pea}$ likewise:

$$
\mathbf{U}_{s} = \left[\, \mathrm{vec}\big(\mathbf{I}^{oat\prime}\big) \;\middle|\; \mathrm{vec}\big(\mathbf{I}^{pea\prime}\big) \,\right]
$$

> **The transpose is load-bearing.** Which flattening you use is fixed by the
> order of the Kronecker product, and getting it wrong is the silent failure
> `docs/specific-combination_kronecker.md` was written about. In
> $\mathbf{G}_{oat} \otimes \mathbf{G}_{pea}$ the **first** factor indexes the
> slower-varying position, so row $(i-1)n_{pea} + j$ is combination $(i,j)$:
> oat slow, pea fast. That is the **row-major** flattening of
> $\mathbf{I}$ — equivalently $\mathrm{vec}$ of its transpose, since
> $\mathrm{vec}$ stacks columns. Verified:
> $\left(\mathbf{G}_{oat} \otimes \mathbf{G}_{pea}\right)_{(i-1)n_{pea}+j,\;(i'-1)n_{pea}+j'}
> = \mathbf{G}_{oat,ii'}\,\mathbf{G}_{pea,jj'}$, and the column-major
> alternative fails.
>
> This is the same fact as `byrow = TRUE` in `beta_from_vector()` and in the
> reshape at `code/dge_ige_functions.R:395`.

Stack those two matrices and the interaction is a **three-way array**
$\mathcal{I}$ of size $n_{oat} \times n_{pea} \times 2$ — oat by pea by trait —
whose covariance is separable in all three directions:

$$
\mathrm{Cov}\big(\mathcal{I}\big) \;=\; \boldsymbol{\Sigma}_{s} \otimes \mathbf{G}_{oat} \otimes \mathbf{G}_{pea}
$$

That is a **tensor normal** distribution: relatedness among oats in mode 1,
relatedness among peas in mode 2, the producer–associate covariance in mode 3.
Written as one long vector with $\mathbf{G}_{oat} \otimes \mathbf{G}_{pea}$
treated as a single opaque kernel, the two-mode structure is invisible. Written
as an array, it is the first thing you see.

**And that structure is exactly what the decomposition exploits.** For a single
trait, $\mathbf{I}^{oat}$ is an $n_{oat} \times n_{pea}$ matrix and its SVD gives
oat scores and pea loadings — the whole subject of
[interaction_decomposition.md](interaction_decomposition.md). That document's
central object, `I_hat = kron_A %*% Beta %*% t(kron_B)`, is the mode-1 and
mode-2 structure of $\mathcal{I}$ made explicit. Seeing the interaction as an
array rather than a vector is what makes "score for an oat, loading for a pea"
the obvious thing to look for.

---

## 6. GMA and SMA in the same notation

The proposal's equations 28–30 relate the producer–associate parameterisation to
the general/specific mixing-ability one. In matrix form they become single
statements about row sums:

$$
GMA^{oat} = \mathbf{U}_{oat}\,\mathbf{1}_2
= \mathrm{Pr}^{oat} + As^{oat \to pea},\qquad
GMA^{pea} = \mathbf{U}_{pea}\,\mathbf{1}_2
= As^{pea \to oat} + \mathrm{Pr}^{pea}
$$

$$
SMA = \mathbf{U}_{s}\,\mathbf{1}_2
= S^{oat \times pea \to oat} + S^{pea \times oat \to pea}
$$

So **GMA is the row sum of the effect matrix** — each accession's total
contribution across both species — and post-multiplying by $\mathbf{1}_2$ is all
that distinguishes the two parameterisations. Its variance follows immediately
from the matrix-normal form:

$$
\mathrm{Var}\big(GMA^{oat}\big) = \left(\mathbf{1}_2' \boldsymbol{\Sigma}_{oat} \mathbf{1}_2\right) \mathbf{G}_{oat}
= \left(\sigma^2_{\mathrm{Pr},oat} + 2\sigma_{\mathrm{Pr} As,oat} + \sigma^2_{As,oat \to pea}\right)\mathbf{G}_{oat}
$$

which shows why the producer–associate covariance matters so much in practice:
it enters GMA's variance **doubled**. A strongly negative
$\sigma_{\mathrm{Pr} As,oat}$ can leave almost no GMA variance even when both component
variances are large — the tension the proposal raises and
[B4I_followups.md](B4I_followups.md) §5 reports credible intervals for.

`accession_effects()` in `code/dge_ige_functions.R` computes exactly this:
`GMA = PrEff + AsEff`.

---

## 7. How this maps onto what the code actually fits

The matrix-variate form is not just tidier — **it is the form `BGLR::Multitrait`
takes.** `fit_producer_associate()` in `code/dge_ige_functions.R` passes a
two-column response matrix and gets a $2\times2$ covariance per term:

| matrix-variate object | in the code |
|---|---|
| $\mathbf{Y}$ | `Y`, a two-column matrix passed as `y =` |
| $\mathbf{Z}_{oat}\,\mathbf{U}_{oat}$ | `ETA$G_oat$X = Z_oat %*% L_oat`, `model = "BRR"` |
| $\boldsymbol{\Sigma}_{oat}$ | `fit$ETA$G_oat$Cov$Omega` |
| $\boldsymbol{\Sigma}_{s}$ | `fit$ETA$G_mix$Cov$Omega` |
| $\boldsymbol{\Sigma}_{\epsilon}$ | `fit$resCov$R`, fitted `type = "UN"` |
| $\mathbf{G}_{oat} \otimes \mathbf{G}_{pea}$ | the row-wise Kronecker basis, truncated at `kron_rank` |

The kinship enters through a factorisation rather than directly, because BGLR
takes a design matrix and not a covariance. With $\mathbf{L}$ such that
$\mathbf{L}\mathbf{L}' = \mathbf{G}_{oat}$ and coefficients
$\mathbf{B}_{oat}$ ($q \times 2$) given $\mathrm{vec}(\mathbf{B}_{oat}) \sim
\mathcal{N}(0, \boldsymbol{\Omega} \otimes \mathbf{I}_q)$:

$$
\mathbf{U}_{oat} = \mathbf{L}\,\mathbf{B}_{oat}
\;\sim\; \mathcal{MN}\!\left(0,\; \mathbf{L}\mathbf{L}',\; \boldsymbol{\Omega}\right)
= \mathcal{MN}\!\left(0,\; \mathbf{G}_{oat},\; \boldsymbol{\Omega}\right)
$$

So `Omega` *is* $\boldsymbol{\Sigma}_{oat}$, and a Bayesian ridge regression on
the factored design matrix *is* the matrix-normal random effect.
[interaction_decomposition.md](interaction_decomposition.md) §1.1 derives the
$\mathbf{L} = \mathbf{U}\boldsymbol{\Lambda}^{1/2}$ construction and why the
$\sqrt{\lambda}$ scaling is what makes this work.

> ### The trait order in the code is the reverse of the proposal's
>
> This document, following the proposal, puts **oat in column 1**. The code
> does not:
>
> ```r
> DGE_IGE_TRAITS <- c("peaYield", "oatYield")
> ```
>
> So the fitted `Omega` has pea first, and $\boldsymbol{\Sigma}_{oat}$ as printed
> by the code is this document's $\boldsymbol{\Sigma}_{oat}$ with its rows and
> columns swapped — the diagonal entries exchange meaning, with
> `Omega[1,1]` the oat's *associate* variance and `Omega[2,2]` its *producer*
> variance.
>
> The code never relies on position. `DGE_IGE_ROLES` maps `Pr` and `As` to trait
> *names* per term and `covariance_components()` looks them up, which is why
> `dge_ige_functions.R` calls the trait ordering "the single highest-consequence
> convention in the project — swapping it turns every producer into an associate
> with no error."
>
> This is the clearest argument for the matrix-variate notation. In the stacked
> form the convention is buried in the order of sub-vectors inside
> $I_2 \otimes Z$ and has to be tracked by hand. In matrix form it is a column
> label on an object you can print.

---

## 8. Why bother

| | stacked | matrix-variate |
|---|---|---|
| Response | one $2n$-vector | one $n \times 2$ matrix |
| Each random term | $\left(I_2 \otimes Z\right)\mathrm{vec}(U)$ | $Z\,U$ |
| Trait covariance | inside a Kronecker product | a named $2\times2$ matrix |
| Role of an effect | position in a stacked vector | column of a matrix |
| Interaction structure | an opaque $n_{oat}n_{pea}$ kernel | a three-way separable array |
| Relation to the code | needs translating | is what BGLR fits |

The stacked form has one genuine advantage, which is why it appears in the
proposal: it is the form in which the model is a *standard* mixed model, so it
can be handed to any mixed-model software that accepts a single response vector
and user-specified covariance structures — which is exactly what the proposal
says Sub-objective 1.4 will test ("asreml, BGLR, BLUPF90"). The matrix-variate
form assumes the software knows about multiple traits.

Both are worth having written down. This one is better for reasoning about the
interaction, and for reading the code.
