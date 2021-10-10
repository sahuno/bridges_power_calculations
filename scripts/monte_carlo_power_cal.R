###samuel ahuno
##date: Jan 17th 2021

library(pwr)
library(ggplot2)
library(data.table)
library(tidyverse)
library(Rmpfr)
library(bigmemory)

##generate sample sizes per group
total_n = seq(7,100000,1)
grp1x.n <- rep(2,length(total_n))
grp2x.n <- total_n-grp1x.n

head(total_n)
head(grp1x.n)
head(grp2x.n)
##get effect sizes
p1=1; p2=0.098
h = abs(2*asin(sqrt(p1))-2*asin(sqrt(p2))); h ## non directional
h2 <- ES.h(p1, p2);h2
desired_p <- 0.05/30


#function()
#{
#}
##the real work
df.results.pwr.2p2n.test <- data.frame(power=NA,stringsAsFactors=FALSE)
for (i in 1:length(grp1x.n)){
  df.results.2p2n <- pwr.2p2n.test(h, n1=grp1x.n[i], n2=grp2x.n[i], sig.level=desired_p,alternative = "greater")$power
  df.results.pwr.2p2n.test <- rbind(df.results.pwr.2p2n.test,df.results.2p2n)
  #p.holder[i,j] <- result.pwr.2p2n.test$sig.level
}

##get results ready
##remove first row  (NA)
df.results.pwr.2p2n.test <- df.results.pwr.2p2n.test[-1,]
df.results.pwr.2p2n.test <- as.data.frame(df.results.pwr.2p2n.test)
names(df.results.pwr.2p2n.test) <- "power"

head(df.results.pwr.2p2n.test)
str(df.results.pwr.2p2n.test)


df_sample <- data.frame(grp1x.n,grp2x.n,key_sample=paste0("n",grp1x.n,"_",grp2x.n))
df.results.pwr.2p2n.test.key <- cbind(df_sample,df.results.pwr.2p2n.test)

summary(df.results.pwr.2p2n.test.key$power)

df.power.above80 <- df.results.pwr.2p2n.test.key %>% filter(power >= 0.80)

dim(df.power.above80)
head(df.power.above80)

df.power.above80 <- df.power.above80[1:100,]
length_df <- NROW(df.power.above80)
#plot.power80 <- ggplot(df.power.above80,aes(x=reorder(as.factor(key_sample),power),y=power)) + geom_point() + theme(axis.text.x = element_text(size=8,angle = 90, vjust = 1,hjust=1)) + 
# scale_y_continuous(labels = function(x) paste0(x, "%"),breaks=c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1))

#ggsave(filename = paste0("figures/power_above_80_percent_first",length_df,"points.pdf"), plot = plot.power80, width = 34, height = 15)



#grp1x.length <- length(grp1x.n)
#grp2x.length <- length(grp2x.n)
#power.holder <- array(numeric(grp1x.length*grp2x.length), dim=c(grp1x.length,grp2x.length),dimnames=list(grp1x.n,grp2x.n))
#p.holder <- array(numeric(grp1x.length*grp2x.length), dim=c(grp1x.length,grp2x.length),dimnames=list(grp1x.n,grp2x.n))
#mat <- mpfrArray(grp1x.length*grp2x.length, dim = c(grp1x.length,grp2x.length), dimnames=list(grp1x.n,grp2x.n))



# for (i in 1:length(grp1x.n)){
#   for (j in 1:length(grp2x.n)){
#     result.pwr.2p2n.test <- pwr.2p2n.test(h, n1=grp1x.n[i], n2=grp2x.n[j], sig.level=desired_p,alternative = "greater")
#     power.holder[i,j] <- result.pwr.2p2n.test$power
#     #p.holder[i,j] <- result.pwr.2p2n.test$sig.level
#   }
# }


############### end  of work ###################



# 
# ###prepare power table for merging
# df.power.holder <- data.table::melt(power.holder,variable.name = c("group_1", "group_2"),value.name = "power")
# names(df.power.holder) <- c("group1", "group2", "power")
# head(df.power.holder)
# df.power.holder <- df.power.holder %>% mutate(key_power = paste0("n",df.power.holder$group1,"_",df.power.holder$group2))
# 
# ##prepare sig level table for merging
# df.p.holder<- data.table::melt(p.holder,variable.name = c("group_1", "group_1"),value.name = "sig.level")
# names(df.p.holder) <- c("group1", "group2", "sig.level")
# df.p.holder <- df.p.holder %>% mutate(key_sig_level = paste0("n",df.p.holder$group1,"_",df.p.holder$group2))
# unique(df.p.holder$key_sig_level)
# 
# ##merge sig.level and power
# df.merged.power_n_sig <- merge(df.power.holder,df.p.holder,by.x="key_power",by.y="key_sig_level")
# head(df.merged.power_n_sig)
# 
# 
# df_sample <- data.frame(grp1x.n,grp2x.n,key_sample=paste0("n",grp1x.n,"_",grp2x.n))
# 
# #unique(df.merged.power_n_sig$key_power)
# 
# ##not found 
# df.merged.power_n_sig %>% filter(key_power %in% df_sample$key_sample)
# df.merged.power_n_sig %>% filter(key_power =="n5_495")
# 
# df.merged.power_n_sig_unique <- distinct(df.merged.power_n_sig, key_power, .keep_all = TRUE)
# summary(df.merged.power_n_sig_unique$power)
# ggplot(df.merged.power_n_sig_unique,aes(x=key_power,y=power)) + geom_point() + theme(axis.text.x = element_text(size=8,angle = 90, vjust = 1,hjust=1)) + 
#   scale_y_continuous(labels = function(x) paste0(x, "%"),breaks=c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1))
# 
# 

#df.merged.power_n_sig[match(df.merged.power_n_sig$key_power,df_sample$key_sample),]
#merge(df.merged.power_n_sig,df_sample,by.y="key_sample",by.x="key_power")


################ scratch ###########

# total_test1 =3
# total_test2 =17
# test.binorm <- rbinom(1000,total_test1,0.9)
# test.binorm2 <- rbinom(1000,total_test2,0.15)
# 
# keep.p <-  data.frame(stringsAsFactors=FALSE)
# for (i in 1:length(test.binorm)) {
#   pid_5_test <- array(c(test.binorm[i],(total_test1 - test.binorm[i]),test.binorm2[i],total_test2 - test.binorm2[i]),dim=c(2,2))
#   fisher_p <- fisher.test(pid_5_test)$p.value
#   keep.p <<- as.data.frame(c(keep.p,fisher_p))
# }
# length(which(keep.p<desired_p))/length(keep.p)
##check tcga, how many case do we see, deeletion of one copy of the genes
# table(test.binorm)
# chisq.sig
# 
# test.binorm[1]
# comp1 = 5 - test.binorm[1]
# test.binorm2[1]
# comp2 <- 15 - test.binorm2[1]
fisher_p$p.value
fisher.test(pid_5)
sample(x=5,size=100,prob=0.8,replace=TRUE)