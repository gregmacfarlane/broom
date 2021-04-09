#' Tidying methods for Apollo model outputs
#' 
#' These functions tidy the coefficients of abstract discrete choice models
#' estimated with the `apollo` package.
#' 
#' @param x an object returned from [apollo::apollo_estimate()]
#' @param nest.params A character vector with the names of nesting coefficient 
#'   parameters. For these parameters the null hypothesis is that the 
#'   value is 1 rather than the traditional 0; the t-statistic for this
#'   parameter will be calculated appropriately if given.
#' @template param_confint
#' @template param_unused_dots
#' 
#' @details As of now, the objects returned from apollo are instances of
#' class `maxLik`.  This requires use of the method call, which might disrupt 
#' some downstream operations. It might be appropriate to manually change the
#' model object class to `apollo`.
#' 
#' @evalRd return_tidy(regression = TRUE)
#' 
#' @examples
#' \dontrun{
#'   apollo_mnl
#'   class(apollo_mnl) <- c("apollo", class(apollo_mnl))
#'   tidy(apollo_mnl) 
#'   tidy(apollo_mnl, se.type = "robust") 
#'   tidy(apollo_mnl, conf.int = TRUE) 
#' 
#'   glance(apollo_mnl)
#' }
#' 
#' @aliases apollo_tidiers
#' @export
#' @family apollo tidiers
#' 
#' @seealso [tidy()], [apollo::apollo_estimate()]
#' 
tidy.apollo <- function(x, conf.int = FALSE, se.type = c("default", "robust"), 
                        nest.params = NULL,
                        conf.level = 0.95, ...){
  
  # get estimated variables 
  varnames <- rownames(x$varcov)
  
  se.type <- match.arg(se.type)
  if(se.type == "default"){
    se <- sqrt(diag(x$varcov))
    robust <- FALSE
  } else {
    se <- sqrt(diag(x$robvarcov))
    robust <- TRUE
  }
  
  # vector of null hypothesis values
  t_null <- rep(0, length(varnames))
  names(t_null) <- varnames
  if(!is.null(nest.params)){
    t_null[nest.params] <- 1
  }
  
  # construct parameter table
  ret <- tibble::tibble(
    term = varnames,
    estimate = x$estimate[varnames],
    std.error = se,
    statistic = (estimate - t_null) / std.error,
    p.val = 2 * pnorm(statistic)
  ) 
  
  if (conf.int) {
    ret <- ret %>%
      dplyr::mutate(
        conf.low = estimate - conf.level * std.error,
        conf.high = estimate + conf.level * std.error
      )
  }
  
  ret
}

#' @templateVar class apollo
#' @template title_desc_glance
#' 
#' @inherit tidy.apollo params examples
#' 
#' @evalRd return_glance(
#'   "logLik", 
#'   "rho2",
#'   "rho20",
#'   "AIC", 
#'   "BIC",
#'   "nobs"
#' )
#' @export
#' @family mlogit tidiers
#' @seealso [glance()], [apollo::apollo_estimate()]
#' 
#' 
#'
glance.apollo  <- function(x, ...) {
  # compute mcfadden r2
  # model log likelihood
  llM <- as.numeric(logLik(x))
  # null model: equal odds for all alternatives
  ll0 <- as.numeric(x$LL0)
  # market shares model: odds equal to chosen proportions
  # cannot see how to compute this from apollo output
  # llC <- sum(x$freq*log(prop.table(x$freq))) 
  
  res <- as_glance_tibble(
    logLik = llM,
    rho20 = 1 - llM / ll0,
    AIC = stats::AIC(x),
    nobs =  x$nObs,
    na_types = "rrri"
  )
  
  res
}

