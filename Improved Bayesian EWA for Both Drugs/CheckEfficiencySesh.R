# Load required libraries
#install.packages("microbenchmark")
library(microbenchmark)

# Function to run a file and measure its execution time
run_file <- function(file_path) {
  env <- new.env()
  system.time(source(file_path, local = env))
}

# Run each file multiple times and compare
benchmark_results <- microbenchmark(
  original = run_file("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/Model_Parameters_BothDrugs.R"),
  new = run_file("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/Model_Params_Improved.R"),
  times = 10  # Number of times to run each file
)

# Print results
print(benchmark_results)

# Plot results
plot(benchmark_results)

# Calculate and print the speedup factor
original_mean <- mean(benchmark_results$time[benchmark_results$expr == "original"])
new_mean <- mean(benchmark_results$time[benchmark_results$expr == "new"])
speedup <- original_mean / new_mean

cat("Speedup factor:", speedup, "\n")
cat("The new version is", ifelse(speedup > 1, "faster", "slower"), "than the original version.\n")