## read.R ##

# Read matrix, meta-data (TO DO: add sleection of the chosen organism)

expression.matrix <- as.matrix(read.csv(
  "example_data/expression_matrix_preprocessed.csv",
  sep = ",", row.names = 1
))

meta <- data.frame(
  srr = colnames(expression.matrix),
  timepoint = rep(c("0h", "12h", "36h"), each = 2)
)
