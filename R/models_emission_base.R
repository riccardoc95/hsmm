library(R6)


EmissionModel <- R6Class(
  "EmissionModel",
  public = list(
    density = NULL,
    initialize = function(params, Pi, post.pi) {
      invisible(self)
    },
    compute_density = function(params, p.array) {
      invisible(self)
    },
    update = function(params, Pi, post.pi) {
      invisible(self)
    }
  )
)
