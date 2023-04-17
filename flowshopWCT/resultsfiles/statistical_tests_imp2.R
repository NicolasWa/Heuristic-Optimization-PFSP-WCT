alg1 <- "ILS_srz_insert_firstNEW.txt" #for ILS
alg2 <- "RII_srz_insert_firstNEW.txt" #for RII
alpha <- 0.05
print(getwd())

RPD_ILS <- read.csv(alg1)$Avg_RPD
RPD_ILS_100jobs <- RPD_ILS[1:30]
RPD_ILS_50jobs <- RPD_ILS[31:60]
RPD_RII <- read.csv(alg2)$Avg_RPD
RPD_RII_100jobs <- RPD_RII[1:30]
RPD_RII_50jobs <- RPD_RII[31:60]

print("WILCOXON TEST FOR 50 JOBS INSTANCES")
p <- wilcox.test (RPD_ILS_50jobs, RPD_RII_50jobs, paired=T)$p.value
print(p)
if (p<alpha) {
  print("The two sets (50 jobs instances) are statistically different")
} else {
  print("The two sets (50jobs instances) are NOT statistically different.")
}

print("WILCOXON TEST FOR 100 JOBS INSTANCES")
p <- wilcox.test (RPD_ILS_100jobs, RPD_RII_100jobs, paired=T)$p.value
print(p)
if (p<alpha) {
  print("The two sets (100 jobs instances) are statistically different")
} else {
  print("The two sets (100jobs instances) are NOT statistically different.")
}

print("Avg RPD ILS (all instances) : ")
mean(RPD_ILS)
print("Avg RPD ILS (50jobs) : ")
mean(RPD_ILS_50jobs)
print("Avg RPD ILS (100jobs) : ")
mean(RPD_ILS_100jobs)


print("Avg RPD RII (all instances) : ")
mean(RPD_RII)
print("Avg RPD RII (50jobs) : ")
mean(RPD_RII_50jobs)
print("Avg RPD RII (100jobs) : ")
mean(RPD_RII_100jobs)


# Creating the plot 50 jobs
plot(RPD_ILS_50jobs, RPD_RII_50jobs, pch = 19, col = "lightblue", main = "Correlation plot for instances of 50 jobs")
# Regression line
abline(lm(RPD_RII_50jobs ~ RPD_ILS_50jobs), col = "red", lwd = 3)
#paste("Correlation:", round(cor(RPD_ILS_100jobs, RPD_RII_100jobs), 4))
corr_50_jobs <- cor.test(x=RPD_ILS_50jobs, y=RPD_RII_50jobs, method = 'spearman')
corr_50_jobs



# Creating the plot 100 jobs
plot(RPD_ILS_100jobs, RPD_RII_100jobs, pch = 19, col = "lightblue", main = "Correlation plot for instances of 100 jobs")
# Regression line
abline(lm(RPD_RII_100jobs ~ RPD_ILS_100jobs), col = "red", lwd = 3)
#paste("Correlation:", round(cor(RPD_ILS_100jobs, RPD_RII_100jobs), 4))
corr_100_jobs <- cor.test(x=RPD_ILS_100jobs, y=RPD_RII_100jobs, method = 'spearman')
corr_100_jobs

all_performances <- read.csv("all_performances.txt")
RPD_RII_tab <- read.csv(alg2)
RPD_ILS_tab <- read.csv(alg1)
