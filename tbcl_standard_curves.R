##TbCl strandard curve
library(ggplot2)
tbcl<-read.delim('~/Documents/GitHub/amphibian_dormancy/tbcl_std_curves.txt', header=T)

ggplot(tbcl, aes(No_spores, Fluorescens))+
  facet_wrap(~Name)+
  geom_point()+
  theme_bw()+
  ylab("Fluorescens")+
  xlab('# Spores/ml')+
  geom_smooth(method='lm')+
  scale_x_log10()

#get lm for all data
linear_model <- lm(Fluorescens ~ No_spores, data = tbcl)

# 3. View the results (Slope and Intercept)
summary(linear_model)

# 4. Extract the equation components
intercept <- coef(linear_model)[1]
slope <- coef(linear_model)[2]

cat("The linear equation is: y =", slope, "* x +", intercept)
#The linear equation is: y = 0.002476504 * x + 71978.03
