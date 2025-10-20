# Set the folder containing your R scripts
folder_path <- "R"  # <-- change this

# Get all .R files in the folder (sorted alphabetically)
r_files <- list.files(folder_path, pattern = "\\.R$", full.names = TRUE)
r_files <- sort(r_files)

# Output file path
output_file <- file.path(folder_path, "combined_scripts.R")

# Combine all scripts into one file
cat("# Combined R script generated on", Sys.time(), "\n\n", file = output_file)

for (f in r_files) {
  cat("# ---- Start of", basename(f), "----\n\n", file = output_file, append = TRUE)
  lines <- readLines(f, warn = FALSE)
  cat(paste(lines, collapse = "\n"), file = output_file, append = TRUE)
  cat("\n\n# ---- End of", basename(f), "----\n\n", file = output_file, append = TRUE)
}

cat("✅ Combined", length(r_files), "files into:", output_file, "\n")