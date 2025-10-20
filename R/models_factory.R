DurationModelFactory <- R6Class("DurationModelFactory", public = list(
  create = function(params) {
    if (is.null(params$covariates_q)){
      type <- "geometric"
    } else{
      type <- "gompertz"
    }
    
    switch(
      type,
      "geometric" = DurationModel.Geometric$new(params),
      "gompertz" = DurationModel.Gompertz$new(params),
      stop("Unknown transition model type: ", type)
    )
  }
))

TransitionModelFactory <- R6Class("TransitionModelFactory", public = list(
  create = function(params, p.array) {
    if (is.null(params$covariates_omega)){
      type <- "homogeneous"
    } else{
      type <- "timevaring"
    }
    switch(
      type,
      "homogeneous" = TransitionModel.Homogeneous$new(params, p.array),
      "timevaring" = TransitionModel.TimeVarying$new(params, p.array),
      stop("Unknown transition model type: ", type)
    )
  }
))

EmissionModelFactory <- R6Class("EmissionModelFactory", public = list(
  create = function(params, Pi, post.pi) {
    type <- params$model_type
    switch(
      type,
      "torus" = TorusModel$new(params, Pi, post.pi),
      "gaussian" = GaussianModel$new(params, Pi, post.pi),
      "poisson" = PoissonModel$new(params, Pi, post.pi),
      "exponential" = ExponentialModel$new(params, Pi, post.pi),
      "gamma" = GammaModel$new(params, Pi, post.pi),
      "beta" = BetaModel$new(params, Pi, post.pi),
      "student" = StudentTModel$new(params, Pi, post.pi),
      stop("Unknown emission model type: ", type)
    )
  }
))