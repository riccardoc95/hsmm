# Cloglog/Gompertz duration model with covariates.
DurationModel.Gompertz <- R6::R6Class(
  "DurationModel.Gompertz",
  inherit = DurationModel,
  public = list(
    initialize = function(params) {
      super$initialize(params)
    }
  )
)
