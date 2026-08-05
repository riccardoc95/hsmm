# Small common interface used by all emission models.
EmissionModel <- R6::R6Class(
  "EmissionModel",
  public = list(
    density = NULL,

    initialize = function(params, Pi, post.pi) {
      invisible(self)
    },

    compute_density = function(params) {
      invisible(self)
    },

    update = function(params, Pi, post.pi) {
      invisible(self)
    }
  )
)
