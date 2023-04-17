alg1 <- "ii_srz_insert_first.txt"
alg2 <- "ii_srz_insert_best.txt"
alpha <- 0.05
print(getwd())
best.known <- read.csv("bestSolutions.txt")

a.cost <- read.csv(alg1)$Solution
a.cost <- 100 * (a.cost - best.known$BS) / best.known$BS # Relative Percentage Deviation
b.cost <- read.csv(alg2)$Solution
b.cost <- 100 * (b.cost - best.known$BS) / best.known$BS # Relative Percentage Deviation

p <- wilcox.test (a.cost, b.cost, paired=T)$p.value
print(p)
if (p<alpha) {
  print("The two sets are statistically different")
} else {
  print("The two sets are NOT statistically different.")
}
