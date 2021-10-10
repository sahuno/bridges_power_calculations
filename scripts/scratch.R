
#generate random number
n1 <- seq(30,150,5) 
n2 <- seq(30,150,5)-15
p<- sample(x=seq(0.01,0.9,0.001), size=length(n1)*length(n2), replace = TRUE, prob = NULL)

p.matrix <- array(data=p,dim=c(length(n1),length(n2)),dimnames=list(n1,n2))

df.p <- data.table::melt(p.matrix,value.name = "power")


#df.p <- arrange(df.power.holder, Var1, Var2)
df.p$n_of_interest <- as.data.frame(ifelse(df.p$Var1[i] == n1[i] & df.p$Var2[i] == n2[i], as.numeric('1'),0))

for (i in 1:length(df.p)) {
  df.p$n_of_interest <- as.data.frame(ifelse(df.p$Var1[i] == n1[i] & df.p$Var2[i] == n2[i], as.numeric('1'),0))
}

i=1
j=1
as.data.frame(n1,n2)
key_df <- data.frame(n1,n2)
key_df$key_subset <- paste(key_df$n1,key_df$n2,sep="_")


df.p$key_df_p <- paste(df.p$Var1,df.p$Var2,sep="_")

df.p %>% filter(key_df_p %in% key_subset)
#which(df.p$

df.p[match(key_df$key_subset,df.p$key_df_p),]
head(df.p)


########################################## Sunday jan 17 2021
n1 <- seq(30,150,5) 
n2 <- seq(30,150,5)-15
h = abs(2*asin(sqrt(p1))-2*asin(sqrt(p2))); h ## non directional

px.holder <- data.frame(grp1=NA,stringsAsFactors=FALSE)
i=1
for (i in 1:length(df.p)) {
  results.pwr.2p2n.test <- pwr.2p2n.test(h, n1=n1[i], n2=n2[i], sig.level=desired_p,alternative = "greater")$power
}
results.pwr.2p2n.test$sig.level

results.pwr.2p2n.test <- data.frame(power=NA,stringsAsFactors=FALSE)
for (i in 1:length(n1)) {
  df <- pwr.2p2n.test(h, n1=n1[i], n2=n2[i], sig.level=desired_p,alternative = "greater")$power
  #df.sig <- pwr.2p2n.test(h, n1=n1[i], n2=n2[i], sig.level=desired_p,alternative = "greater")$sig.level
  results.pwr.2p2n.test <- rbind(results.pwr.2p2n.test,df)
}
mydf <- data.frame(px.holder)
