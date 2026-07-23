#' CensRCalc: Analytical Censoring Parameter Determination
#'
#' @description
#' \loadmathjax
#' Provides calculations for the expected censoring
#' proportion under independent censoring to determine a
#' censoring distribution parameter that achieves a target censoring
#' proportion in simulation studies. The expected censoring proportion is
#' evaluated by numerical integration over time and covariates, with support for
#' administrative censoring and accrual periods.
#'
#' @details
#' Let the event time be
#' \mjeqn{T^*\mid\mathbf{X}\sim F_{T^*\mid\mathbf{X}}}{T*~F_T*|X} with survival
#' function \mjeqn{S_{T^*\mid\mathbf{X}}=1-F_{T^*\mid \mathbf{X}}
#' }{S_T*|X = 1-F_T*|X} and the random censoring time be \mjeqn{
#' C_{\text{rnd}}\sim F_{C_{\text{rnd}}\mid \lambda_{C_{\text{rnd}}}}
#' }{C_rnd|lambda_C_rnd~F_C_rnd|lambda_C_rnd} with survival
#' function \mjeqn{S_{C_{\text{rnd}}\mid\lambda_{C_{\text{rnd}}}}=
#' 1-F_{C_{\text{rnd}}\mid\lambda
#' {C_{\text{rnd}}}}}{S_C(t|lambda_C_rnd)=1-F_C(t|lambda_C_rnd)}, where
#' \mjeqn{\lambda_{C_{\text{rnd}}}}{lambda_C_rnd}
#' is the random censoring parameter of interest. Administrative censoring
#' occurs at time \mjeqn{\tau_{adm}}{tau_adm}. Accrual takes place over
#' \mjeqn{[0,\tau_{\text{acc}}]}{[0,tau_acc]} where the random accrual time is
#' \mjeqn{R\sim U(0,\tau_{\text{acc}})}{R~U(0,tau_acc)}. The equivalent random
#' administrative censoring time after accrual is
#' \mjeqn{C_{\text{adm}}=\tau_{adm}-R}{C_adm=tau_adm-R}. Hence,
#' \mjeqn{C_{\text{adm}}\sim U(\tau_{adm}-\tau_{\text{acc}},
#' \tau_{adm})}{C_adm~U(tau_adm-tau_acc, tau_adm)}.
#'
#' Under these assumptions, the probability of censoring dependent on
#' \mjeqn{\lambda_{C_{\text{rnd}}}}{lambda_C_rnd} and the covariate vector
#' \mjeqn{\mathbf{X}}{X} is
#'
#' \mjdeqn{
#' \begin{aligned}
#' P(\delta=0\mid \mathbf{X},\lambda_{C_{\text{rnd}}}) =
#' &\int_{0}^{\tau_{\text{adm}} - \tau_{\text{acc}}}
#' S_{T^{*}|\mathbf{X}}(c_{\text{rnd}})
#' f_{C_{\text{rnd}}|\lambda_{C_{\text{rnd}}}}(c_{\text{rnd}})dc_{\text{rnd}}
#' + \\\\\\
#' &\int_{\tau_{\text{adm}} - \tau_{\text{acc}}}^{\tau_{\text{adm}}}
#' S_{T^{*}|\mathbf{X}}(c_{\text{rnd}})S_{C_{\text{adm}}}(c_{\text{rnd}})
#' f_{C_{\text{rnd}}|\lambda_{C_{\text{rnd}}}}(c_{\text{rnd}})dc_{\text{rnd}}
#' + \\\\\\
#' &\int_{\tau_{\text{adm}} - \tau_{\text{acc}}}^{\tau_{\text{adm}}}
#' S_{T^{*}|\mathbf{X}}(C_{\text{adm}})
#' S_{C_{\text{rnd}}| \lambda_{C_{\text{rnd}}}}(C_{\text{adm}})
#' f_{C_{\text{adm}}}(C_{\text{adm}}) dc_{\text{adm}}
#' \end{aligned}
#' }{P(delta=0|X,lambda_C_rnd) =
#' \int_0^(tau_adm - tau_acc)
#' S_T*|X(c_rnd)f_C_rnd|lambda_C_rnd(c_rnd)dc_rnd +
#' \int_(tau_adm - tau_acc)^tau_adm
#' S_T*|X(c_rnd)S_C_adm(c_rnd)f_C_rnd|lambda_C_rnd(c_rnd)dc_rnd +
#' \int_(tau_adm - tau_acc)^tau_adm
#' S_T*|X(c_adm)S_C_rnd| lambda_C_rnd(c_adm)f_C_adm(c_adm) dc_adm}
#'
#' Next we take the expectation over the covariate distribution to get the
#' expected censoring proportion:
#' \mjdeqn{P(\delta=0\mid\lambda_{C_{\text{rnd}}})=E_X\Bigg[P(\delta=0\mid
#' \mathbf{X},\lambda_{C_{\text{rnd}}})\Bigg]}{P(delta=0|lambda_C_rnd)=
#' E_X[P(delta=0|X,lambda_C_rnd)]}. This is implemented in
#' [estimate_cens_prop()]. To obtain the censoring parameter
#' \mjeqn{\lambda_{C_{\text{rnd}}}}{lambda_C_rnd}
#' that achieves a target censoring proportion \mjeqn{p_C}{p_C}, we solve the
#' root finding problem
#' \mjeqn{\lambda_{C_{\text{rnd}}}}{lambda_C_rnd}.
#' This is implemented in [find_cens_param()].
#'
#' Special cases include no accrual time
#' (\mjeqn{\tau_{\text{acc}}=0}{tau_acc=0}) and no administrative
#' censoring (\mjeqn{\tau_{adm}=\infty}{tau_adm=Inf}). See the package vignette
#' \code{vignette("full_details", package = "CensRCalc")} for details and
#' derivations.
#'
#' Numerically, the integrals are evaluated using adaptive cubature
#' ([calculus::integral()]), which was extended to handle discrete and mixed
#' discrete/continuous covariates.
#'
#' @import mathjaxr
#' @importFrom Rdpack reprompt
#' @rdname CensRCalc
#'
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
