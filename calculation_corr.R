#calculate correlations

#read in data
data<-read.delim("amphibian_dormancy/sample_data.txt", header=T)
str(dat)

##write a function to do the calculation (SPearman correlation)
#parts to change: 
#data- whatever your data frame is named, 
#the 'fresh_CFU' and 'Frozen_counts' to whatever column you're naming
#under part 2, change the species to whatever you're comparing

calculate_correlation <- function(factor_level, data) {
  subset_data <- subset(data, Species == factor_level)
  correlation_test <- cor.test(subset_data$Frozen_counts, subset_data$Fresh_CFU, method='spearman')
  return(correlation_test)
}

#Calculate correlation with p-value for each level of the factor
factor_levels <- unique(data$Species)
correlation_results <- lapply(factor_levels, function(level) {
  calculate_correlation(level, data)
})

#view results
correlation_results
