library(hsmm)

options(hsmm.make_plots = FALSE)

files <- list.files(
  "inst/tests/manual",
  pattern = "^[0-9][0-9]_.*[.]R$",
  full.names = TRUE
)

results <- data.frame(
  file = basename(files),
  passed = FALSE,
  error = "",
  stringsAsFactors = FALSE
)

for (i in seq_along(files)) {
  cat("\n\nRunning", basename(files[i]), "\n")

  tryCatch({
    source(files[i], local = new.env())
    results$passed[i] <- TRUE
  }, error = function(e) {
    results$error[i] <- conditionMessage(e)
    cat("ERROR:", conditionMessage(e), "\n")
  })
}

cat("\n\nFINAL SUMMARY\n")
print(results, row.names = FALSE)

if (!all(results$passed)) {
  stop("At least one manual test did not pass")
}

cat("\nAll manual tests passed.\n")
