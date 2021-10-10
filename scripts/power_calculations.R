###power calculations
library(pwr)
library(ggplot2)
library(pairs)
library(plotmatrix)
library(data.table)
library(tidyverse)


#power.prop.test(p1=0.15, p2=0.30, power=0.85, sig.level=0.05) ### test, how many number of patients do you need to power0.85


#How to get the direction of effect size
#n1=138; n2=138
p1=0.31; p2=0.15
h = abs(2*asin(sqrt(p1))-2*asin(sqrt(p2))); h ## non directional
h2 <- ES.h(p1, p2);h2
desired_p <- 0.05/30

### what power, do you get if you unequal numbers in each group 
# test.result.p <- pwr.2p2n.test(h, n1=n1, n2=n2, sig.level=0.05)
# test.result.p$power
# plot(test.result.p)
# test.result.p$sig.level
###real worrk
grp1.n <- seq(30,5000,5)  ##group 1, N
grp2.n <- seq(30,5000,5)-15  ## group 2, N - 15

total_n = seq(500,5000,100)
grp1x.n <- 5
grp2x.n <- total_n-grp1x.n

grp1.length <- length(grp1.n)
grp2.length <- length(grp2.n)

power.holder <- array(numeric(grp1.length*grp2.length), dim=c(grp1.length,grp2.length),dimnames=list(grp1.n,grp2.n))
p.holder <- array(numeric(grp1.length*grp2.length), dim=c(grp1.length,grp2.length),dimnames=list(grp1.n,grp2.n))
#p.holder <- array(numeric(grp1.length*grp2.length), dim=c(grp1.length,grp2.length),dimnames=list(grp1.n,grp2.n))
#power.holder <- data.frame(grp1=NA,stringsAsFactors=FALSE)

for (i in 1:length(grp1.n)){
  for (j in 1:length(grp2.n)){
    result.pwr.2p2n.test <- pwr.2p2n.test(h, n1=grp1.n[i], n2=grp2.n[j], sig.level=desired_p,alternative = "greater")
    power.holder[i,j] <- result.pwr.2p2n.test$power
    p.holder[i,j] <- result.pwr.2p2n.test$sig.level
    #return(result.pwr.2p2n.test)
  }
}



df.power.holder <- data.table::melt(power.holder,variable.name = c("group_1(N)", "group_1(N-15)"),value.name = "power")
df.p.holder<- data.table::melt(p.holder,variable.name = c("group_1(N)", "group_1(N-15)"),value.name = "sig.level")

head(df.power.holder)
filter(df.power.holder,power<1)
table(df.power.holder$power)

table(df.p.holder$sig.level)
filter(df.power.holder,power<1)





dim(df.power.holder)
length.df.power.holder <- seq(1,NROW(df.power.holder),1)
df.power.holder$varx <- length.df.power.holder
#df.power.holder <- df.power.holder %>% mutate(grp_of_interest=ifelse(df.power.holder$Var1 == df.power.holder$Var2, as.numeric('1'),0))

df.power.holder

#distinct(df.power.holder, Var1,.keep_all = TRUE)
#filter(df.power.holder,grp_of_interest==1)
#df.power.holder %>% filter(Var1 %in% grp1.n,Var2%in%grp2)
#df.power.holder <- arrange(df.power.holder, Var1, Var2)

df.power.holder.subset <- df.power.holder[1:30,]
ggplot(df.power.holder.subset,aes(x=varx,y=power)) + geom_point()
##sanity checks
#df.power.holder.distinct <- distinct_at(df.power.holder, vars(c(Var1, Var2)))
#dim(df.power.holder.distinct)
#filter(df.power.holder,Var1==25 & Var2==5)

###use a list
list.grp1.n <- as.list(grp1.n)
list.grp2.n <- as.list(grp2.n)
names_sample_groups <- paste0("N=",grp1.n," * ",grp2.n)
df_sample_sizes <- data.frame(stringsAsFactors=FALSE)

for (i in 1:length(list.grp1.n)) {
  result.pwr.2p2n.test <- pwr.2p2n.test(h, n1=list.grp1.n[[i]], n2=list.grp2.n[[i]], sig.level=0.05)$power
  #result.pwr.2p2n.test <- ceiling(result.pwr.2p2n.test$power)
  df_sample_sizes <- as.data.frame(c(df_sample_sizes,result.pwr.2p2n.test))
}
names(df_sample_sizes) <- names_sample_groups


as.vector(df_sample_sizes)
str(df_sample_sizes)
head(df_sample_sizes)
dim(df_sample_sizes)

# for (i in 1:length(g1.n)){
#     result.pwr.2p2n.test <- pwr.2p2n.test(h, n1=i, n2=i-15, sig.level=0.05)
#     #samsize[j,i] <- ceiling(result.pwr.2p2n.test$n)
#     #return(result.pwr.2p2n.test)
# }
# 
# 
# for (k in 1:g1.n.length) {
#   for (l in 1:g2.n.length) {
#     print(paste("k =", k, "l= ",l))
#   }
# }
# 
# # range of correlations
# r <- seq(.1,.5,.01)
# nr <- length(r)
# 
# # power values
# p <- seq(.4,.9,.1)
# np <- length(p)
# 
# # obtain sample sizes
# samsize <- array(numeric(nr*np), dim=c(nr,np))
# 
# samsize <- array(numeric(nr*np), dim=c(nr,np))
# for (i in 1:np){
#   for (j in 1:nr){
#     result <- pwr.r.test(n = NULL, r = r[j],
#                          sig.level = .05, power = p[i],
#                          alternative = "two.sided")
#     samsize[j,i] <- ceiling(result$n)
#   }
# }
