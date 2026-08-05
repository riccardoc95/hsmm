# Duration model without external dwell-time covariates.
DurationModel.Geometric <- R6::R6Class(
  "DurationModel.Geometric",
  inherit = DurationModel,
  public = list(
    initialize = function(params) {
      super$initialize(params)
    }
  )
)
