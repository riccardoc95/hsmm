DurationModelFactory <- R6::R6Class(
  "DurationModelFactory",
  public = list(
    create = function(params) {
      if (is.null(params$covariates_q) || ncol(params$covariates_q) == 0) {
        DurationModel.Geometric$new(params)
      } else {
        DurationModel.Gompertz$new(params)
      }
    }
  )
)

TransitionModelFactory <- R6::R6Class(
  "TransitionModelFactory",
  public = list(
    create = function(params, p.array) {
      if (is.null(params$covariates_omega) || ncol(params$covariates_omega) == 0) {
        TransitionModel.Homogeneous$new(params, p.array)
      } else {
        TransitionModel.TimeVarying$new(params, p.array)
      }
    }
  )
)

EmissionModelFactory <- R6::R6Class(
  "EmissionModelFactory",
  public = list(
    create = function(params, Pi, post.pi) {
      type <- tolower(params$model_type)
      switch(
        type,
        torus = TorusModel$new(params, Pi, post.pi),
        gaussian = GaussianModel$new(params, Pi, post.pi),
        poisson = PoissonModel$new(params, Pi, post.pi),
        exponential = ExponentialModel$new(params, Pi, post.pi),
        gamma = GammaModel$new(params, Pi, post.pi),
        beta = BetaModel$new(params, Pi, post.pi),
        student = StudentTModel$new(params, Pi, post.pi),
        studentt = StudentTModel$new(params, Pi, post.pi),
        tstudent = StudentTModel$new(params, Pi, post.pi),
        var = VarGaussianModel$new(params, Pi, post.pi),
        var_gaussian = VarGaussianModel$new(params, Pi, post.pi),
        vargaussian = VarGaussianModel$new(params, Pi, post.pi),
        stop("Unknown emission model type: ", type)
      )
    }
  )
)
