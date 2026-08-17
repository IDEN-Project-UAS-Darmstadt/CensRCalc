# Full Details and Usage

## Introduction

This vignette provides an overview of the functionalities and usage of
the
[`help(CensRCalc)`](https://iden-project-uas-darmstadt.github.io/CensRCalc/reference/CensRCalc.md)
package. This package is used to calculate the parameter of the random
censoring distribution so that the expected censoring proportion in a
simulated dataset is equal to a given target value. This package only
covers right-censoring.

The theory behind this package is based on the work by Wan
([2017](#ref-Wan2017)).

## Mathematical Background

In this chapter we describe the mathematical background of the censoring
parameter calculation.

### Censoring without accrual

First, we consider the special case where patients are accrued at the
same time, i.e. at the start of a clinical trial. We assume independent
censoring, i.e. the event time conditional on a \\d\\-dimensional
covariate vector \\\mathbf{X}\\, that is \\T^{\*}\mid\mathbf{X}\sim
F\_{T^\*\mid\mathbf{X}}\\, is independent of the random censoring time
\\C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}\sim
F\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}\\ where
\\\lambda\_{C\_{\text{rnd}}}\\ is the parameter of the random censoring
distribution:

\\\begin{equation\*} T^{\*}\mid\mathbf{X} \perp\\\\\\\perp
C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}} \end{equation\*}\\

The administrative censoring time \\\tau\_{\text{adm}}\\ is fixed. Thus,
the actual censoring time \\C\\ is the minimum of the random censoring
time \\C\_{\text{rnd}}\\ and the administrative censoring time
\\\tau\_{\text{adm}}\\:

\\\begin{equation\*} C := \min(C\_{\text{rnd}}, \tau\_{\text{adm}})
\end{equation\*}\\

Then the observed time \\T\\ is the minimum of the event time \\T^{\*}\\
and the censoring time \\C\\:

\\\begin{equation\*} T := \min(T^{\*}, C) \end{equation\*}\\

The event indicator is given by: \\\begin{equation\*} \delta := I(T^{\*}
\leq C) = \begin{cases} 1 & \text{, } T^{\*} \leq C \\ 0 & \text{,
otherwise} \end{cases} \end{equation\*}\\

Let \\p\_{C}\\ be the target value of the expected censoring proportion.
The goal is to calculate the parameter \\\lambda\_{C\_{\text{rnd}}}\\ of
the random censoring distribution \\F\_{C\_{\text{rnd}}}\\ so that the
expected censoring proportion in a simulated dataset is \\p\_{C}\\. For
this, we derive a formula for \\\lambda\_{C\_{\text{rnd}}}\\ depending
on the expected censoring proportion \\p\_{C}\\ and the parameters of
the event time distribution \\F\_{T^{\*}}\\ as well as the
administrative censoring time \\\tau\_{\text{adm}}\\. This leads to an
equation that can be solved for \\\lambda\_{C\_{\text{rnd}}}\\ as a root
finding problem.

We start with the probability of censoring for a given subject,
\\P(\delta=0\mid\mathbf{X})\\. Using the assumption of independent
censoring we know that

\\\begin{align\*} P(\delta=0\mid\mathbf{X}, \lambda\_{C\_{\text{rnd}}})
&= P(\\C\_{\text{rnd}} \leq T^{\*}, C\_{\text{rnd}} \leq
\tau\_{\text{adm}}\\ \vee \\\tau\_{\text{adm}} \< T^{\*},
\tau\_{\text{adm}} \< C\_{\text{rnd}}\\\mid\mathbf{X},
\lambda\_{C\_{\text{rnd}}}) \\ &= P(C\_{\text{rnd}} \leq T^{\*},
C\_{\text{rnd}} \leq \tau\_{\text{adm}}\mid\mathbf{X},
\lambda\_{C\_{\text{rnd}}}) + P(\tau\_{\text{adm}} \< T^{\*},
\tau\_{\text{adm}} \< C\_{\text{rnd}}\mid\mathbf{X},
\lambda\_{C\_{\text{rnd}}}) \\ &= P(C\_{\text{rnd}} \leq \min(T^{\*},
\tau\_{\text{adm}})\mid\mathbf{X}, \lambda\_{C\_{\text{rnd}}}) +
P(\tau\_{\text{adm}} \< T^{\*}\mid\mathbf{X})P(\tau\_{\text{adm}} \<
C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}) \\ &=
\\\begin{aligned}\[t\] &\int\_{0}^{\tau\_{\text{adm}}} P(c \leq
T^{\*}\mid\mathbf{X})dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c) +
\\ &P(\tau\_{\text{adm}} \< T^{\*}\mid\mathbf{X})P(\tau\_{\text{adm}} \<
C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}) \end{aligned} \\ &=
\int\_{0}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c)dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c) +
S\_{T^{\*}\mid\mathbf{X}}(\tau\_{\text{adm}})S\_{C\_{\text{rnd}}\mid
\lambda\_{C\_{\text{rnd}}}}(\tau\_{\text{adm}}) \\ &=
\int\_{0}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c)f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c)dc +
S\_{T^{\*}\mid\mathbf{X}}(\tau\_{\text{adm}})S\_{C\_{\text{rnd}}\mid
\lambda\_{C\_{\text{rnd}}}}(\tau\_{\text{adm}}) \end{align\*}\\ where
\\S\_{T^{\*}\mid\mathbf{X}}\\ and
\\S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}\\ are the survival
functions of the event time and random censoring time conditional on the
covariate vector \\\mathbf{X}\\ and the random censoring parameter
\\\lambda\_{C\_{\text{rnd}}}\\, respectively.

The first term (the integral) describes the probability of random
censoring before the administrative censoring. The second term describes
the probability that the administrative censoring happens before the
event or random censoring. Note that \\P(\delta=0\mid\mathbf{X},
\lambda\_{C\_{\text{rnd}}})\\ is a random variable which takes on values
depending on the covariates \\\mathbf{X}\\ of a specific subject. Now we
have to integrate \\P(\delta=0\mid\mathbf{X},
\lambda\_{C\_{\text{rnd}}})\\ over the covariate distribution \\F\_{X}\\
to get the probability of censoring in the simulated dataset,
independent of the covariates of a specific subject:

\\\begin{equation\*} P(\delta = 0\mid\lambda\_{C\_{\text{rnd}}}) =
E\_{X}(P(\delta=0\mid\mathbf{X}, \lambda\_{C\_{\text{rnd}}})) =
\int\_{D\_{X}} P(\delta=0\mid\mathbf{X}=x, \lambda\_{C\_{\text{rnd}}})
dF\_{X} = \int\_{D\_{X}} P(\delta=0\mid\mathbf{X}=x,
\lambda\_{C\_{\text{rnd}}})f\_{X}(x) dx, \end{equation\*}\\

where \\D\_{X}\\ is the domain of the covariate vector \\X\\. This
d-dimensional integral must usually be computed numerically. The
packages `help(calculus)` and `help(cubature)` are used for this
purpose.

Now we can solve for \\\lambda\_{C\_{\text{rnd}}}\\ using a root finding
algorithm. The function to be solved is given by:

\\\begin{equation\*} P(\delta = 0\mid\lambda\_{C\_{\text{rnd}}}) =
p\_{C} \Longleftrightarrow P(\delta = 0\mid\lambda\_{C\_{\text{rnd}}}) -
p\_{C} = 0. \end{equation\*}\\

To solve this root finding problem we use the
[`stats::uniroot`](https://rdrr.io/r/stats/uniroot.html) function that
implements a bracketing root finding algorithm.

### Censoring with accrual

In the general case, patients are accrued within a time span \\\[0,
\tau\_{\text{acc}}\]\\. We consider uniform distributed accrual times
and thus, the random accrual time is \\R \sim U(\[0,
\tau\_{\text{acc}}\])\\. The other parameters and variables are
identical to those in the special case. Stochastically, the random
process of accrual can be equivalently modelled by making the
administrative censoring time a random variable, \\C\_{\text{adm}}\\,
given by the following definition:

\\\begin{equation\*} C\_{\text{adm}} := \tau\_{\text{adm}} - R.
\end{equation\*}\\ In other words, after accrual a patient can only be
observed until the administrative censoring time is reached.
Consequently, the time a patient can be observed after accrual, the
follow-up time, is given by the difference between the administrative
censoring time and the actual accrual time of that patient.

It follows that \\C\_{\text{adm}} \sim U(\[\tau\_{\text{adm}} -
\tau\_{\text{acc}}, \tau\_{\text{adm}}\])\\. Furthermore, we assume
independence between \\T^{\*}\mid\mathbf{X}\\ and
\\C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}\\ as well as
\\C\_{\text{adm}}\\:

\\\begin{equation\*} T^{\*}\mid\mathbf{X} \perp\\\\\\\perp
C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}} \quad \wedge \quad
T^{\*}\mid\mathbf{X} \perp\\\\\\\perp C\_{\text{adm}} \end{equation\*}\\

Again, we start with the probability of censoring for a given subject:
\\\begin{align\*} P(\delta=0\mid\mathbf{X}, \lambda\_{C\_{\text{rnd}}})
&= P(\\C\_{\text{rnd}} \< T^{\*}, C\_{\text{rnd}} \< C\_{\text{adm}}\\
\vee \\C\_{\text{adm}} \< T^{\*}, C\_{\text{adm}} \<
C\_{\text{rnd}}\\\mid\mathbf{X}, \lambda\_{C\_{\text{rnd}}}) \\ &=
P(C\_{\text{rnd}} \< T^{\*}, C\_{\text{rnd}} \<
C\_{\text{adm}}\mid\mathbf{X}, \lambda\_{C\_{\text{rnd}}}) +
P(C\_{\text{adm}} \< T^{\*}, C\_{\text{adm}} \<
C\_{\text{rnd}}\mid\mathbf{X}, \lambda\_{C\_{\text{rnd}}}) \\ &=
\\\begin{aligned}\[t\] &\int\_{0}^{\tau\_{\text{adm}}} P(c\_{\text{rnd}}
\< T^{\*}, c\_{\text{rnd}} \<
C\_{\text{adm}}\mid\mathbf{X})dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}}) +
\\ &\int\_{\tau\_{\text{adm}} - \tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
P(c\_{\text{adm}} \< T^{\*}, c\_{\text{adm}} \<
C\_{\text{rnd}}\mid\mathbf{X},
\lambda\_{C\_{\text{rnd}}})dF\_{C\_{\text{adm}}}(c\_{\text{adm}})
\end{aligned} \\ &= \\\begin{aligned}\[t\]
&\int\_{0}^{\tau\_{\text{adm}}} P(c\_{\text{rnd}} \<
T^{\*}\mid\mathbf{X})P(c\_{\text{rnd}} \< c\_{\text{adm}})
dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}}) +
\\ &\int\_{\tau\_{\text{adm}} - \tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
P(c\_{\text{adm}} \< T^{\*}\mid\mathbf{X})P(c\_{\text{adm}} \<
C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}})
dF\_{C\_{\text{adm}}}(c\_{\text{adm}}) \end{aligned} \\ &=
\\\begin{aligned}\[t\] &\int\_{0}^{\tau\_{\text{adm}}} P(c\_{\text{rnd}}
\< T^{\*}\mid\mathbf{X})P(c\_{\text{rnd}} \< c\_{\text{adm}})
dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}}) +
\\ &\int\_{\tau\_{\text{adm}} - \tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
P(C\_{\text{adm}} \< T^{\*}\mid\mathbf{X})P(C\_{\text{adm}} \<
C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}})
dF\_{C\_{\text{adm}}}(c\_{\text{adm}}) \end{aligned} \\ &=
\\\begin{aligned}\[t\] &\int\_{0}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})S\_{C\_{\text{adm}}}(c\_{\text{rnd}})dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}}) +
\\ &\int\_{\tau\_{\text{adm}} - \tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{adm}})S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{adm}})
dF\_{C\_{\text{adm}}}(c\_{\text{adm}}) \end{aligned} \\ &=:
I\_{\text{rnd}} + I\_{\text{adm}} \end{align\*}\\ where
\\I\_{\text{rnd}}\\ describes the random censoring and
\\I\_{\text{adm}}\\ the administrative censoring.

Since \\C\_{\text{adm}} \in \[\tau\_{\text{adm}} - \tau\_{\text{acc}},
\tau\_{\text{adm}}\]\\ we have \\P(c\_{\text{rnd}} \< C\_{\text{adm}}) =
S\_{C\_{\text{adm}}}(c\_{\text{rnd}}) = 1\\ for \\c\_{\text{rnd}} \in
(0, \tau\_{\text{adm}} - \tau\_{\text{acc}})\\. So we can rewrite
\\I\_{\text{rnd}}\\: \\\begin{equation\*} I\_{\text{rnd}} =
\int\_{0}^{\tau\_{\text{adm}} - \tau\_{\text{acc}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}}) +
\int\_{\tau\_{\text{adm}} - \tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})S\_{C\_{\text{adm}}}(c\_{\text{rnd}})dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})
\end{equation\*}\\

which gives us: \\\begin{align\*} P(\delta=0\mid\mathbf{X},
\lambda\_{C\_{\text{rnd}}}) &= \\\begin{aligned}\[t\]
&\int\_{0}^{\tau\_{\text{adm}} - \tau\_{\text{acc}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}}) +
\\ &\int\_{\tau\_{\text{adm}} - \tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})S\_{C\_{\text{adm}}}(c\_{\text{rnd}})dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}}) +
\\ &\int\_{\tau\_{\text{adm}} - \tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{adm}})S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{adm}})
dF\_{C\_{\text{adm}}}(c\_{\text{adm}}) \end{aligned} \\ &=
\\\begin{aligned}\[t\] &\int\_{0}^{\tau\_{\text{adm}} -
\tau\_{\text{acc}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})dc\_{\text{rnd}} +
\\ &\int\_{\tau\_{\text{adm}} - \tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})S\_{C\_{\text{adm}}}(c\_{\text{rnd}})f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})dc\_{\text{rnd}} +
\\ &\int\_{\tau\_{\text{adm}} - \tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{adm}})S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{adm}})f\_{C\_{\text{adm}}}(c\_{\text{adm}})
dc\_{\text{adm}} \end{aligned} \\ &= I\_{\text{rnd,nadm}} +
I\_{\text{rnd,adm}} + I\_{\text{adm}} \end{align\*}\\

The first integral \\I\_{\text{rnd,nadm}}\\ describes the random
censoring when no administrative censoring can occur. The second
integral \\I\_{\text{rnd,adm}}\\ describes the random censoring when
also administrative censoring can occur and the third one
\\I\_{\text{adm}}\\ describes the administrative censoring.

Like in the case without accrual time, to get the probability of
censoring in the simulated dataset, independent of the covariates of a
specific subject, we integrate \\P(\delta=0\mid\mathbf{X},
\lambda\_{C\_{\text{rnd}}})\\ over the covariate distribution \\F\_{X}\\
and then the parameter \\\lambda\_{C\_{\text{rnd}}}\\ can be computed as
a root finding problem.

### Edge cases

Lastly, we consider two edge cases, when \\\tau\_{\text{acc}} = 0\\ and
when \\\tau\_{\text{adm}} = \infty\\. In the first case, the accrual
time span is zero, i.e. all patients are accrued at or before the
beginning of the study. So the general formula should collapse to the
special formula we showed before. If \\\tau\_{\text{acc}} = 0\\ it
follows that \\C\_{\text{adm}} \sim U(\[\tau\_{\text{adm}},
\tau\_{\text{adm}}\]) = \delta\_{\tau\_{\text{adm}}}\\, where
\\\delta\_{\tau\_{\text{adm}}}\\ is the Dirac distribution at
\\\tau\_{\text{adm}}\\. It follows:

\\\begin{align\*} I\_{\text{adm}} &= \int\_{\tau\_{\text{adm}} -
\tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{adm}})S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{adm}})
dF\_{C\_{\text{adm}}}(c\_{\text{adm}}) \\ &=
\int\_{\tau\_{\text{adm}}}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{adm}})S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{adm}})
d\delta\_{\tau\_{\text{adm}}}(c\_{\text{adm}}) \\ &=
\int\_{\\\tau\_{\text{adm}}\\}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{adm}})S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{adm}})
d\delta\_{\tau\_{\text{adm}}}(c\_{\text{adm}}) \\ &=
S\_{T^{\*}\mid\mathbf{X}}(\tau\_{\text{adm}})S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(\tau\_{\text{adm}})
\end{align\*}\\

For \\I\_{\text{rnd,adm}}\\ we have: \\\begin{align\*}
I\_{\text{rnd,adm}} = \int\_{\tau\_{\text{adm}}}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c)S\_{C\_{\text{adm}}}(c)
dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c) =
\int\_{\\\tau\_{\text{adm}}\\}
S\_{T^{\*}\mid\mathbf{X}}(c)S\_{C\_{\text{adm}}}(c)
dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c) = 0
\end{align\*}\\ because
\\F\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}\\ is a continuous
distribution function (hence, a continuous measure), so the integral
over one point is zero.

So when \\\tau\_{\text{acc}} = 0\\ we get: \\\begin{equation\*}
P(\delta=0\mid\mathbf{X}, \lambda\_{C\_{\text{rnd}}}) =
\int\_{0}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c)f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c)dc +
S\_{T^{\*}\mid\mathbf{X}}(\tau\_{\text{adm}})S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(\tau\_{\text{adm}})
\end{equation\*}\\ which is the special case we saw in the [Censoring
without accrual](#sec:no_accrual).

In the second case, the administrative censoring time is infinite, i.e.
patients are observed until the event or random censoring occurs. This
results in three improper integrals in the general formula. We can solve
these by considering the limit of these integrals as
\\\tau\_{\text{adm}} \to \infty\\. First we take a look at
\\I\_{\text{rnd}}\\. If \\\tau\_{\text{adm}} \to \infty\\, then
\\\tau\_{\text{adm}}-\tau\_{\text{acc}} \to \infty\\. So we have:

\\\begin{align\*} I\_{\text{rnd,adm}} &= \int\_{\infty}^{\infty}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})S\_{C\_{\text{adm}}}(c\_{\text{rnd}})dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})
\\ &= \lim\_{a \to \infty} \int\_{a}^{a}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})S\_{C\_{\text{adm}}}(c\_{\text{rnd}})dF\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})
\\ &= \lim\_{a \to \infty} 0 = 0. \end{align\*}\\

Next, we consider \\I\_{\text{adm}}\\. If \\\tau\_{\text{adm}} \to
\infty\\, then \\C\_{\text{adm}} \sim U(\[\infty, \infty\]) =
\delta\_{\infty} = \lim\_{a \to \infty} \delta\_{a}\\. It follows:

\\\begin{align\*} I\_{\text{adm}} &= \int\_{\infty}^{\infty}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{adm}})S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{adm}})
dF\_{C\_{\text{adm}}}(c\_{\text{adm}}) \\ &= \lim\_{a \to \infty}
\int\_{a}^{a}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{adm}})S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{adm}})
d\delta\_{a}(c\_{\text{adm}}) \\ &= \lim\_{a \to \infty}
S\_{T^{\*}\mid\mathbf{X}}(a)S\_{
C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(a) = 0. \end{align\*}\\

Note that \\\lim\_{a \to \infty} S\_{X}(a) = 0\\ for every continuous
survival function \\S\_{X}\\.

So, when \\\tau\_{\text{adm}} \to \infty\\ the general formula collapses
to: \\\begin{equation\*} P(\delta=0\mid\mathbf{X},
\lambda\_{C\_{\text{rnd}}}) = \int\_{0}^{\infty}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})dc\_{\text{rnd}}
\end{equation\*}\\

We see when \\\tau\_{\text{adm}} \to \infty\\, the accrual time does not
affect the calculation of the censoring parameter. In other words,
accrual time becomes irrelevant because each patient is observed until
the event or random censoring.

## Basic Usage

``` r

library(CensRCalc)
library(progressr) # For progress bar during root finding
library(logger) # For logging messages

tau_acc <- 1
tau_adm <- 3


event_surv <- function(t, x) 1 - punif(t, min = 0, max = x)
cens_surv <- function(t, lambda_c) 1 - punif(t, min = 0, max = lambda_c)
cens_dens <- function(t, lambda_c) dunif(t, min = 0, max = lambda_c)
cov_dens <- function(x) dunif(x, min = 1, max = 2)
cov_bounds <- list(x = c(1, 2))
bounds <- c(0.1, 10)

log_threshold(DEBUG, namespace = "CensRCalc") # Set logger to DEBUG level
estimator <- estimate_cens_prop(event_surv, cov_dens, cov_bounds,
  time_accrual = tau_acc,
  time_admin_cens = tau_adm,
  cens_surv = cens_surv,
  cens_dens = cens_dens,
  target_bounds = bounds
)
#> DEBUG [2026-08-17 11:33:09] Density total integral calculated as: 0.999559635073658
#> DEBUG [2026-08-17 11:33:09] Density total integral calculated as: 0.999559635073658
#> DEBUG [2026-08-17 11:33:09] Density total integral calculated as: 0.999999903558124
#> DEBUG [2026-08-17 11:33:09] Density total integral calculated as: 0.999999903558124
#> DEBUG [2026-08-17 11:33:09] Density total integral calculated as: 0.999999998533039
#> DEBUG [2026-08-17 11:33:09] Density total integral calculated as: 0.999999998533039
#> DEBUG [2026-08-17 11:33:09] Density total integral calculated as: 0.999999977140023
#> DEBUG [2026-08-17 11:33:09] Density total integral calculated as: 0.999999977140023
#> DEBUG [2026-08-17 11:33:09] Density total integral calculated as: 1.00000017787208
#> DEBUG [2026-08-17 11:33:09] Density total integral calculated as: 1.00000017787208
log_threshold(INFO, namespace = "CensRCalc") # Reset logger to INFO level

lambda_c_to_plot <- seq(bounds[1], bounds[2], length.out = 100)
cens_prop <- sapply(lambda_c_to_plot, estimator)
plot(lambda_c_to_plot, cens_prop,
  type = "l", xlab = expression(lambda[C]),
  ylab = "Censoring Proportion"
)
```

![](full_details_files/figure-html/basic_estimator-1.png)

``` r

p_c <- 0.3 # Target censoring proportion

with_progress({
  result <- find_cens_param(p_c, event_surv, cov_dens, cov_bounds,
    time_accrual = tau_acc,
    time_admin_cens = tau_adm,
    cens_surv = cens_surv,
    cens_dens = cens_dens,
    target_bounds = bounds
  )
  print(result)
})
#> $parameter
#> [1] 2.499999
#> 
#> $cens_prop
#> [1] 0.3

lambda_c <- result$parameter
```

``` r

set.seed(123)

n <- 10000
x <- runif(n, min = 1, max = 2) # Covariate
t_true <- runif(n, min = 0, max = x) # Event times uniform (0, X)
c_rnd <- runif(n, min = 0, max = lambda_c) # Random censoring times uniform
r <- runif(n, min = 0, max = tau_acc) # Recruitment times uniform (0, tau_acc)
c_adm <- tau_adm - r # Administrative follow-up after recruitment
c <- pmin(c_rnd, c_adm)
t_obs <- pmin(t_true, c)
delta <- as.integer(t_true <= c)

mean(delta == 0) # Empirical censoring proportion
#> [1] 0.2989
```

## Analytical Background of the Example

For this, we use the simplest distributions for the event and censoring
times, i.e. uniform distributions: \\\begin{align\*}
T^{\*}\mid\mathbf{X} &\sim \text{U}(0, \lambda\_{T^{\*}}) \\
C\_{\text{rnd}} &\sim \text{U}(0, \lambda\_{C\_{\text{rnd}}})
\end{align\*}\\ where \\\lambda\_{T^{\*}} = \beta_0 + \beta_1X\\. Let
\\\beta_0 = 0\\, \\\beta_1 = 1\\ and \\X \sim \text{U}(1,2)\\.

Let \\\tau\_{\text{acc}} = 1\\, i.e. patients are accrued within a one
year time span. The special case is not shown because there is nothing
to be aware of when setting \\\tau\_{\text{acc}} = 0\\. Furthermore,
\\\tau\_{\text{adm}} = 3\\, i.e. patients are observed for three years.
That means \\C\_{\text{adm}} \sim U(\[2, 3\])\\. The target value of the
expected censoring proportion is \\p\_{C} = 0.3\\. The goal is to
calculate the parameter \\\lambda\_{C\_{\text{rnd}}}\\ of the random
censoring distribution so that the expected censoring proportion in a
simulated dataset is \\p\_{C}\\.

For calculating the censoring probability for a specific subject, we
need the survival functions of the event time, random censoring time and
the administrative censoring time as well as the density functions of
the random and administrative censoring time. The survival functions are
given by:

\\\begin{align\*} S\_{T^{\*}\mid\mathbf{X}}(t) &= 1 - F\_{T^{\*}}(t) =
1 - \frac{t-0}{X-0} = \frac{X-t}{X} \\
S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(t) &= 1 -
F\_{C\_{\text{rnd}}}(t) = 1 - \frac{t-0}{\lambda\_{C\_{rnd}}-0} =
\frac{\lambda\_{C\_{rnd}}-t}{\lambda\_{C\_{rnd}}} \\
S\_{C\_{\text{adm}}}(t) &= 1 - F\_{C\_{\text{adm}}}(t) = 1 - \frac{t -
(\tau\_{\text{adm}} - \tau\_{\text{acc}})}{\tau\_{\text{adm}} -
(\tau\_{\text{adm}} - \tau\_{\text{acc}})} = 1 - \frac{t - 2}{1} = 3 - t
\end{align\*}\\

The density functions are given by: \\\begin{align\*}
f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(t) &=
\frac{1}{\lambda\_{C\_{rnd}}-0} = \frac{1}{\lambda\_{C\_{rnd}}} \\
f\_{C\_{\text{adm}}}(t) &= \frac{1}{\tau\_{\text{adm}} -
(\tau\_{\text{adm}} - \tau\_{\text{acc}})} = 1 \end{align\*}\\

The probability is then calculated as follows: \\\begin{align\*}
P(\delta=0\mid\mathbf{X}, \lambda\_{C\_{\text{rnd}}}) &=
\\\begin{aligned}\[t\] &\int\_{0}^{2}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})dc\_{\text{rnd}} +
\\ &\int\_{2}^{3}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})S\_{C\_{\text{adm}}}(c\_{\text{rnd}})f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})dc\_{\text{rnd}} +
\\ &\int\_{2}^{3}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{adm}})S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{adm}})f\_{C\_{\text{adm}}}(c\_{\text{adm}})
dc\_{\text{adm}} \end{aligned} \end{align\*}\\

Since the support of \\T^\*\\ is \\\[0,X\]\\ with \\X \in \[1,2\]\\, the
survival function \\S\_{T^{\*}\mid\mathbf{X}}(t) = 0\\ for all \\t \geq
2\\. Thus, the second and third integrals, which integrate over
\\\[2,3\]\\, are zero. The first integral has to be split up into two
integrals:

\\\begin{align\*} P(\delta=0\mid\mathbf{X}, \lambda\_{C\_{\text{rnd}}})
&= \int\_{0}^{2}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})dc\_{\text{rnd}}
\\ &= \int\_{0}^{X}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})dc\_{\text{rnd}} +
\int\_{X}^{2}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})dc\_{\text{rnd}}
\\ \end{align\*}\\

Using the same reason as above, the second integral over \\\[X,2\]\\, is
also zero. Thus, we have:

\\\begin{align\*} P(\delta=0\mid\mathbf{X}, \lambda\_{C\_{\text{rnd}}})
&= \int\_{0}^{X}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})dc\_{\text{rnd}}
\\ &= \int\_{0}^{X} \frac{X - c\_{\text{rnd}}}{X} \cdot
\frac{1}{\lambda\_{C\_{rnd}}} dc\_{\text{rnd}} \\ &=
\frac{1}{\lambda\_{C\_{rnd}} X} \int\_{0}^{X} (X - c\_{\text{rnd}})
dc\_{\text{rnd}} \\ &= \frac{1}{\lambda\_{C\_{rnd}} X} \left\[ X
c\_{\text{rnd}} - \frac{c\_{\text{rnd}}^2}{2} \right\]\_{0}^{X} \\ &=
\frac{1}{\lambda\_{C\_{rnd}} X} \left( X^2 - \frac{X^2}{2} \right) \\ &=
\frac{X}{2 \lambda\_{C\_{rnd}}} \end{align\*}\\

Next we calculate the expected censoring proportion in the simulated
dataset, independent of the covariates of a specific subject:

\\\begin{align\*} P(\delta=0\mid\lambda\_{C}) &=
E_X\left\[P(\delta=0\mid\mathbf{X}, \lambda\_{C})\right\] \\ &=
\int\_{D\_{\mathbf{X}}} P(\delta=0\mid\mathbf{X}=x,
\lambda\_{C\_{rnd}})f\_{X}(x)dx \\ &= \int\_{1}^{2} \frac{X}{2
\lambda\_{C\_{rnd}}}f\_{X}(x)dx \\ &= \int\_{1}^{2} \frac{X}{2
\lambda\_{C\_{rnd}}} \cdot 1 dx \\ &= \frac{1}{2 \lambda\_{C\_{rnd}}}
\int\_{1}^{2}X dx \\ &= \frac{1}{2 \lambda\_{C\_{rnd}}}
\left\[\frac{X^2}{2} \right\]\_{1}^{2} \\ &= \frac{1}{4
\lambda\_{C\_{rnd}}} \left(4 - 1\right) \\ &= \frac{3}{4
\lambda\_{C\_{rnd}}} \end{align\*}\\

In the last step we can calculate the parameter \\\lambda\_{C\_{rnd}}\\:
\\\begin{align\*} &P(\delta=0\mid\lambda\_{C\_{\text{rnd}}})
\overset{!}{=} 0.3 \\ &\leftrightarrow \frac{3}{4 \lambda\_{C\_{rnd}}} =
0.3 \\ &\leftrightarrow \frac{3}{4 \cdot 0.3} = \lambda\_{C\_{rnd}} \\
&\leftrightarrow \lambda\_{C\_{rnd}} = \frac{30}{12} = 2.5
\end{align\*}\\

## Session Info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] logger_0.4.2    progressr_1.0.0 CensRCalc_0.1.1
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6         knitr_1.51        rlang_1.3.0       xfun_0.60        
#>  [5] otel_0.2.0        textshaping_1.0.5 calculus_1.1.0    jsonlite_2.0.0   
#>  [9] glue_1.8.1        backports_1.5.1   htmltools_0.5.9   ragg_1.5.2       
#> [13] sass_0.4.10       rmarkdown_2.31    evaluate_1.0.5    jquerylib_0.1.4  
#> [17] fastmap_1.2.0     yaml_2.3.12       lifecycle_1.0.5   compiler_4.6.1   
#> [21] mathjaxr_2.0-0    fs_2.1.0          Rcpp_1.1.2        htmlwidgets_1.6.4
#> [25] systemfonts_1.3.2 digest_0.6.39     R6_2.6.1          cubature_2.1.4-1 
#> [29] Rdpack_2.6.6      rbibutils_2.4.1   checkmate_2.3.4   bslib_0.12.0     
#> [33] tools_4.6.1       pkgdown_2.2.1     cachem_1.1.0      desc_1.4.3
```

## References

Wan, Fei. 2017. “Simulating Survival Data with Predefined Censoring
Rates for Proportional Hazards Models.” *Statistics in Medicine* 36
(February): 838–54. <https://doi.org/10.1002/sim.7178>.
