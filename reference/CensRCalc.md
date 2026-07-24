# CensRCalc: Analytical Censoring Parameter Determination

Calculations for the expected censoring proportion under independent
censoring to determine a censoring distribution parameter that achieves
a target censoring proportion in simulation studies. The expected
censoring proportion is evaluated by numerical integration over time and
covariates, with support for administrative censoring and accrual
periods.

## Details

Let the event time be \\T^\*\mid\mathbf{X}\sim F\_{T^\*\mid\mathbf{X}}\\
with survival function
\\S\_{T^\*\mid\mathbf{X}}=1-F\_{T^\*\mid\mathbf{X}}\\ and the random
censoring time be \\C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}\sim
F\_{C\_{\text{rnd}} \mid\lambda\_{C\_{\text{rnd}}}}\\ with survival
function
\\S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}=1-F\_{C\_{\text{rnd}}
\mid\lambda\_{C\_{\text{rnd}}}}\\, where \\\lambda\_{C\_{\text{rnd}}}\\
is the random censoring parameter of interest. Administrative censoring
occurs at time \\\tau\_{\text{adm}}\\. Accrual takes place over
\\\[0,\tau\_{\text{acc}}\]\\ where the random accrual time is \\R\sim
U(0,\tau\_{\text{acc}})\\. The equivalent random administrative
censoring time after accrual is
\\C\_{\text{adm}}=\tau\_{\text{adm}}-R\\. Hence, \\C\_{\text{adm}}\sim
U(\tau\_{\text{adm}}-\tau\_{\text{acc}}, \tau\_{\text{adm}})\\.

Under these assumptions, the probability of censoring dependent on
\\\lambda\_{C\_{\text{rnd}}}\\ and the covariate vector \\\mathbf{X}\\
is

\\ \begin{aligned} P(\delta=0\mid \mathbf{X},\lambda\_{C\_{\text{rnd}}})
= &\int\_{0}^{\tau\_{\text{adm}} - \tau\_{\text{acc}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})
f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})dc\_{\text{rnd}} +
\\ &\int\_{\tau\_{\text{adm}} - \tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(c\_{\text{rnd}})S\_{C\_{\text{adm}}}(c\_{\text{rnd}})
f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(c\_{\text{rnd}})dc\_{\text{rnd}} +
\\ &\int\_{\tau\_{\text{adm}} - \tau\_{\text{acc}}}^{\tau\_{\text{adm}}}
S\_{T^{\*}\mid\mathbf{X}}(C\_{\text{adm}})
S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(C\_{\text{adm}})
f\_{C\_{\text{adm}}}(C\_{\text{adm}})dc\_{\text{adm}} \end{aligned} \\

Next we take the expectation over the covariate distribution to get the
expected censoring proportion:
\\P(\delta=0\mid\lambda\_{C\_{\text{rnd}}})=E_X\left\[P(\delta=0\mid
\mathbf{X},\lambda\_{C\_{\text{rnd}}})\right\].\\ This is implemented in
[`estimate_cens_prop()`](https://iden-project-uas-darmstadt.github.io/CensRCalc/reference/estimate_cens_prop.md).
To obtain the censoring parameter \\\lambda\_{C\_{\text{rnd}}}\\ that
achieves a target censoring proportion \\p_C\\, we solve the root
finding problem \\P(\delta=0\mid\lambda\_{C\_{\text{rnd}}}) - p_C = 0\\

Special cases include no accrual time (\\\tau\_{\text{acc}}=0\\) and no
administrative censoring (\\\tau\_{\text{adm}}=\infty\\). See the
package vignette
[`vignette("full_details", package = "CensRCalc")`](https://iden-project-uas-darmstadt.github.io/CensRCalc/articles/full_details.md)
for details.

Numerically, the integrals are evaluated using adaptive cubature
([`calculus::integral()`](https://calculus.eguidotti.com/reference/integral.html)),
which was extended to handle discrete and mixed discrete/continuous
covariates.

## See also

Useful links:

- <https://iden-project-uas-darmstadt.github.io/CensRCalc>

- Report bugs at
  <https://github.com/iden-project-uas-darmstadt/CensRCalc/issues>

## Author

**Maintainer**: Lukas Klein <lukas.klein@h-da.de>

Authors:

- Lukas Klein <lukas.klein@h-da.de>

- Henrik Stahl <henrik.stahl@h-da.de>

Other contributors:

- Antje Jahn-Eimermacher \[contributor\]

- Gunter Grieser \[contributor\]
