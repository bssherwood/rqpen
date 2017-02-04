\name{rqPen-package}
\alias{rqPen-package}
\docType{package}
\title{
Penalized quantile regression for LASSO, SCAD, and MCP penalty functions including group penalties
}
\description{
This package provides functions to find solutions to penalized quantile regression problems.  Throughout this package, the estimated coefficients are the minimizers of the penalized quantile regression objective function:
\deqn{\beta = \frac{1}{n}\sum_{i=1}^{n} \rho_{\tau}(y_i - x_i^T \beta) + \sum_{j=1}^{p} p_{\lambda}(|\beta_{j}|)},
where \deqn{\rho_{\tau}(u) = u(\tau - I(u < 0))}.  This package can handle three different penalty functions with \eqn{\lambda > 0}:

LASSO: \deqn{p_{\lambda}(|\beta_j|) = \lambda|\beta_j|}

SCAD: \deqn{p_{\lambda}(|\beta_j|) = \lambda|\beta_j|I(0\leq |\beta_j| < \lambda) + \frac{a\lambda|\beta_j|-(\beta_j^2+\lambda^2)/2}{a-1}I(\lambda \leq |\beta_j| \leq a\lambda) + \frac{(a+1)\lambda^2}{2}I(|\beta_j| > a\lambda),}
for \eqn{a > 2}

MCP: \deqn{p_{\lambda}(|\beta_j|) = \lambda(|\beta_j|-\frac{\beta_{j}^2}{2a\lambda})I(0 \leq |\beta_j| \leq a\lambda) + \frac{a\lambda^2}{2}I(|\beta_j| > a\lambda),}
for \eqn{a > 1}.
}

