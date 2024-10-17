#1. Simple MCM and MMCM

mhccexecutelong <- function(nvec,apmat)
{
  k <- length(nvec)
  n <- sum(nvec)
  smatch <- as.matrix(nonbimatch(distancematrix(as.matrix(dist(apmat, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))))$matches)
  multcm <- rep(0,k)
  cs <- c(0,cumsum(nvec))
  for(i in 1:k)
  {
    for(j in (cs[i]+1):(cs[i+1]))
    {
      multcm[i]<-multcm[i]+((as.numeric(smatch[j,3])>cs[i])&&(as.numeric(smatch[j,3])<=cs[i+1]))
    }
  }
  multcm<-multcm/2
  A<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in i:k)
    {
      for(l in (cs[i]+1):(cs[i+1]))
      {
        A[i,j]<-A[i,j]+((as.numeric(smatch[l,3])>cs[j])&&(as.numeric(smatch[l,3])<=cs[j+1]))
      }
    }
  }
  av<-t(A)[lower.tri(t(A))]
  return(list(as.matrix(multcm),as.matrix(av)))
}
#nvec<- c(100,130,110,45) if there are 4 classes and they have 100, 130, 110, 45 samples, resp.
#nvec<-c(2,3,2); d<-3; X1,X2 in class1, X3,X4,X5 in class2, X6,X7 in class3. X1<- (3,2,4), X2<- (.5,2,-1), X3<-(0,0,2), X4<-(-1,2,3),...
#apmat<-[3 2 4
#          .5 2 -1
#           0 0 2
#           -1 2 3
# ...
# ]      
#r<-mhccexecutelong(nvec,apmat)
#r[[2]]
#av<-r[[2]]: (A_{1,2},A_{1,3},...,A_{1,K},A_{2,3},A_{2,4},...,A_{2,K},A_{3,4},...) (bold A)



mhcccreate<-function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  mu1<-rep(0,k)
  sig1<-matrix(0,k,k)
  for(i in 1:k)
  {
    mu1[i]<-(nvec[i]*(nvec[i]-1))/(2*(n-1))
    sig1[i,i]<-(((nvec[i]*(nvec[i]-1))/(2*(n-1)))*(1-((nvec[i]*(nvec[i]-1))/(2*(n-1)))))+((nvec[i]*(nvec[i]-1)*(nvec[i]-2)*(nvec[i]-3))/(4*(n-1)*(n-3)))
  }
  for(i in 1:k)
  {
    for(j in setdiff(1:k,i))
    {
      sig1[i,j]<-((nvec[i])*(nvec[j])*(nvec[i]-1)*(nvec[j]-1))/(2*((n-1)^2)*(n-3))
    }
  }
  mu<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in i:k)
    {
        mu[i,j]<-(((nvec[i])*(nvec[j]))/(n-1))
    }
  }
  muv<-t(mu)[lower.tri(t(mu))]
  bigsig<-matrix(0,(((k^2)-k)/2),(((k^2)-k)/2))
  for(i in 1:(k-1))
  {
    for(j in (i+1):k)
    {
      bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((i-1)*k)-((i*(i-1))/2)+j-i)]<- ((nvec[i]*nvec[j]*(nvec[i]-1)*(nvec[j]-1))/((n-1)*(n-3)))+(((nvec[i]*nvec[j])/(n-1))*(1-((nvec[i]*nvec[j])/(n-1))))
    }
  }
  if(k>=3)
  {
    for(i in 1:(k-2))
    {
      for(j in (i+1):(k-1))
      {
        for(l in (j+1):k)
        {
          bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-((nvec[i]*(nvec[i]-1)*nvec[j]*nvec[l])/((n-1)*(n-3))) - ((((nvec[i])^2)*nvec[j]*nvec[l])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-((nvec[i]*(nvec[i]-1)*nvec[j]*nvec[l])/((n-1)*(n-3))) - ((((nvec[i])^2)*nvec[j]*nvec[l])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-((nvec[j]*(nvec[j]-1)*nvec[i]*nvec[l])/((n-1)*(n-3))) - ((((nvec[j])^2)*nvec[i]*nvec[l])/((n-1)^2))
          bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-((nvec[j]*(nvec[j]-1)*nvec[i]*nvec[l])/((n-1)*(n-3))) - ((((nvec[j])^2)*nvec[i]*nvec[l])/((n-1)^2))
          bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-((nvec[l]*(nvec[l]-1)*nvec[j]*nvec[i])/((n-1)*(n-3))) - ((((nvec[l])^2)*nvec[j]*nvec[i])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-((nvec[l]*(nvec[l]-1)*nvec[j]*nvec[i])/((n-1)*(n-3))) - ((((nvec[l])^2)*nvec[j]*nvec[i])/((n-1)^2))
        }
      }
    }
  }
  if(k>=4)
  {
    for(i in 1:(k-3))
    {
      for(j in (i+1):(k-2))
      {
        for(l in (j+1):(k-1))
        {
          for(m in (l+1):k)
          {
            bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((l-1)*k)-((l*(l-1))/2)+m-l)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((j-1)*k)-((j*(j-1))/2)+m-j)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((i-1)*k)-((i*(i-1))/2)+m-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+m-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((j-1)*k)-((j*(j-1))/2)+m-j),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((l-1)*k)-((l*(l-1))/2)+m-l),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
          }
        }
      }
    }
  }
 return(list(as.matrix(mu1),sig1,as.matrix(muv),bigsig))
}

#s<-mhcccreate(nvec)
#mean1<-s[[3]], var1<-s[[4]]
#statistic<- t(av-mean1)%*%solve(var1)%*%(av-mean1)
#Reject if statistic > z_{0.05}.

#Install the multicross package.

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------


#2. Feature Selection using MCM and MMCM in the Location Setup

#final step selection using univariate test selects much less variables than  true number.
library(MASS)
library(nbpMatching)
library(taRifx)
#source("pval_mmcc_test8.R")
mhcccreate<-function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  mu1<-rep(0,k)
  sig1<-matrix(0,k,k)
  for(i in 1:k)
  {
    mu1[i]<-(nvec[i]*(nvec[i]-1))/(2*(n-1))
    sig1[i,i]<-(((nvec[i]*(nvec[i]-1))/(2*(n-1)))*(1-((nvec[i]*(nvec[i]-1))/(2*(n-1)))))+((nvec[i]*(nvec[i]-1)*(nvec[i]-2)*(nvec[i]-3))/(4*(n-1)*(n-3)))
  }
  for(i in 1:k)
  {
    for(j in setdiff(1:k,i))
    {
      sig1[i,j]<-((nvec[i])*(nvec[j])*(nvec[i]-1)*(nvec[j]-1))/(2*((n-1)^2)*(n-3))
    }
  }
  mu<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in i:k)
    {
      mu[i,j]<-(((nvec[i])*(nvec[j]))/(n-1))
    }
  }
  muv<-t(mu)[lower.tri(t(mu))]
  bigsig<-matrix(0,(((k^2)-k)/2),(((k^2)-k)/2))
  for(i in 1:(k-1))
  {
    for(j in (i+1):k)
    {
      bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((i-1)*k)-((i*(i-1))/2)+j-i)]<- ((nvec[i]*nvec[j]*(nvec[i]-1)*(nvec[j]-1))/((n-1)*(n-3)))+(((nvec[i]*nvec[j])/(n-1))*(1-((nvec[i]*nvec[j])/(n-1))))
    }
  }
  if(k>=3)
  {
    for(i in 1:(k-2))
    {
      for(j in (i+1):(k-1))
      {
        for(l in (j+1):k)
        {
          bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-((nvec[i]*(nvec[i]-1)*nvec[j]*nvec[l])/((n-1)*(n-3))) - ((((nvec[i])^2)*nvec[j]*nvec[l])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-((nvec[i]*(nvec[i]-1)*nvec[j]*nvec[l])/((n-1)*(n-3))) - ((((nvec[i])^2)*nvec[j]*nvec[l])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-((nvec[j]*(nvec[j]-1)*nvec[i]*nvec[l])/((n-1)*(n-3))) - ((((nvec[j])^2)*nvec[i]*nvec[l])/((n-1)^2))
          bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-((nvec[j]*(nvec[j]-1)*nvec[i]*nvec[l])/((n-1)*(n-3))) - ((((nvec[j])^2)*nvec[i]*nvec[l])/((n-1)^2))
          bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-((nvec[l]*(nvec[l]-1)*nvec[j]*nvec[i])/((n-1)*(n-3))) - ((((nvec[l])^2)*nvec[j]*nvec[i])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-((nvec[l]*(nvec[l]-1)*nvec[j]*nvec[i])/((n-1)*(n-3))) - ((((nvec[l])^2)*nvec[j]*nvec[i])/((n-1)^2))
        }
      }
    }
  }
  if(k>=4)
  {
    for(i in 1:(k-3))
    {
      for(j in (i+1):(k-2))
      {
        for(l in (j+1):(k-1))
        {
          for(m in (l+1):k)
          {
            bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((l-1)*k)-((l*(l-1))/2)+m-l)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((j-1)*k)-((j*(j-1))/2)+m-j)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((i-1)*k)-((i*(i-1))/2)+m-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+m-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((j-1)*k)-((j*(j-1))/2)+m-j),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((l-1)*k)-((l*(l-1))/2)+m-l),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
          }
        }
      }
    }
  }
  return(list(as.matrix(mu1),sig1,as.matrix(muv),bigsig))
}


mhccexecutelong<-function(nvec,apmat)
{
  k<-length(nvec)
  n<-sum(nvec)
  smatch<-as.matrix(nonbimatch(distancematrix(as.matrix(dist(apmat, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))))$matches)
  multcm<-rep(0,k)
  cs<-c(0,cumsum(nvec))
  for(i in 1:k)
  {
    for(j in (cs[i]+1):(cs[i+1]))
    {
      multcm[i]<-multcm[i]+((as.numeric(smatch[j,3])>cs[i])&&(as.numeric(smatch[j,3])<=cs[i+1]))
    }
  }
  multcm<-multcm/2
  A<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in i:k)
    {
      for(l in (cs[i]+1):(cs[i+1]))
      {
        A[i,j]<-A[i,j]+((as.numeric(smatch[l,3])>cs[j])&&(as.numeric(smatch[l,3])<=cs[j+1]))
      }
    }
  }
  av<-t(A)[lower.tri(t(A))]
  return(list(as.matrix(multcm),as.matrix(av)))
}
library(mvtnorm)
#rho1 = 0.5
# rho2 = 0.5
n = 1000
theta = 0.2#0.2#1
Xtotal= c()
#set.seed(28)
nvec = rep(200,5)
set.seed(1234)
l = sample(100,24)
mu = rep(0,100)
for(i in 1:length(nvec)){
  set.seed(i)
  mu[l] = theta*i
  Xtotal = rbind(Xtotal,rmvnorm(nvec[i],mean=mu))
}

d = ncol(Xtotal)#100# d is total number of variables to start with
k = length(nvec)
ll<-mhcccreate(nvec)
mu1<-ll[[1]]
sig1<-ll[[2]]
muv<-ll[[3]]
bigsig<-ll[[4]]

pval_test<-function(X){
  X = as.matrix(X)
  lll<-mhccexecutelong(nvec,X)
  av<-lll[[2]]
  stbig<-t(as.matrix(av)-as.matrix(muv))%*%solve(bigsig)%*%(as.matrix(av)-as.matrix(muv))
  # print(as.numeric(stbig))
  stbig<-(1-pchisq(stbig,df=((k^2)-k)/2))*(d/ncol(X))
  #print(stbig)
  stsmall<-(sum(av)-sum(muv))/sum(bigsig)
  stsmall<-pnorm(stsmall)*(d/ncol(X))
  return(stbig)
}



hclust_fun2<-function(X,alpha){
  cor_X <- cor(X, use = "pairwise.complete.obs", method = "pearson")
  
  hclust_resu <- hclust(as.dist(1-cor_X),"single")#"single"
  
  l = list()
  l[[1]] =  abs(hclust_resu$merge[1,])
  for (i in 2:nrow(hclust_resu$merge)){
    temp1 = hclust_resu$merge[i,1]
    temp2 = hclust_resu$merge[i,2]
    if((temp1>0) & (temp2>0)){
      l[[length(l)+1]] = c(l[[temp1]],l[[temp2]])
    }
    if((temp1>0) & (temp2<0)){
      l[[length(l)+1]] = unique(c(l[[temp1]],-temp2))
    }
    
    if((temp1<0) & (temp2>0)){
      l[[length(l)+1]] = unique(c(l[[temp2]],-temp1))
    }
    if((temp1<0) & (temp2<0)){
      l[[length(l)+1]] = c(-temp2,-temp1)
    }  
  }
  
  
  merged_matr = hclust_resu$merge
  merged_list = list()
  for(i in 1:nrow(merged_matr)){
    merged_list[[i]] = list()
    merged_list[[i]][[1]] = merged_matr[i,]
    merged_list[[i]][[2]] = l[[i]]
    
  }
  
  
  attempt2<-function(n){
    
    r = list()
    a = merged_list[[n]][[1]][1]
    b = merged_list[[n]][[1]][2]
    
    if(a>0){
      a1 = sort(merged_list[[a]][[2]])
    }else{
      a1 = -a
    }
    if(b>0){
      b1 = sort(merged_list[[b]][[2]])
    }else{
      b1 = -b
    }
    
    c = sort(unique(c(a1,b1)))
    r[[1]] = c
    r[[4]] = as.numeric(pval_test(X[,c])<alpha)
    
    if(r[[4]]==0){#& length(r[[1]])>1
      if(length(r[[1]])==1){
        r[[2]] = "empty"
        r[[3]] = "empty"
        
      }
      else{
        r[[2]] = array(-1616)
        r[[3]] = array(-2828)
      }
      
    }
    else{
      if(a<0){
        r[[2]] = list()
        r[[2]][[1]] = -a
        r[[2]][[4]] = as.numeric(pval_test(X[,-a])<alpha)
        r[[2]][[2]] = "empty"#array()
        r[[2]][[3]] = "empty"#array()
        # r[[2]] = -a
        # return(r)
      }
      else{
        r[[2]] = attempt2(a)
      }
      
      if(b<0){
        r[[3]] = list()
        r[[3]][[1]] = -b
        r[[3]][[4]] = as.numeric(pval_test(X[,-b])<alpha)
        r[[3]][[2]] = "empty"#array()
        r[[3]][[3]] = "empty"#array()
      }
      else{
        r[[3]] = attempt2(b)
      }
      
    }
    
   
    
    return(r)
  }#ends attempt2()
  
  return(attempt2(nrow(merged_matr)))
  
}
tree_s<-hclust_fun2(Xtotal,alpha=0.05)

var1 = list()
pval_mmcc<-function(l,alpha){
  
  if(class(l[[2]])=="character"){
    if(l[[4]]==1){
      print(l[[1]])
      assign("var1",list.append(var1,l[[1]]), envir = .GlobalEnv)
    }
    
  }
  else if(l[[4]]==1 & class(l[[2]])=="list"){
    if(l[[2]][[4]]==0 & l[[3]][[4]]==0){
      print(l[[1]])
      assign("var1",list.append(var1,l[[1]]), envir = .GlobalEnv)
      
    }
    else if(l[[2]][[4]]==1 & l[[3]][[4]]==0){
      pval_mmcc(l[[2]],alpha)
      
    }
    else if(l[[2]][[4]]==0 & l[[3]][[4]]==1){
      pval_mmcc(l[[3]],alpha)
      
    }
    else{
      pval_mmcc(l[[2]],alpha)     
      pval_mmcc(l[[3]],alpha)
      
      
    }
    
    
  }
  return(var1)
  
}

selected_var = sort(unlist(pval_mmcc(tree_s,0.05)))
true_var = sort(l)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------


#3. Feature Selection using MCM and MMCM in the Scale Setup

#final step selection using univariate test selects much less variables than  true number.
library(MASS)
library(nbpMatching)
library(taRifx)
#source("pval_mmcc_test8.R")
mhcccreate<-function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  mu1<-rep(0,k)
  sig1<-matrix(0,k,k)
  for(i in 1:k)
  {
    mu1[i]<-(nvec[i]*(nvec[i]-1))/(2*(n-1))
    sig1[i,i]<-(((nvec[i]*(nvec[i]-1))/(2*(n-1)))*(1-((nvec[i]*(nvec[i]-1))/(2*(n-1)))))+((nvec[i]*(nvec[i]-1)*(nvec[i]-2)*(nvec[i]-3))/(4*(n-1)*(n-3)))
  }
  for(i in 1:k)
  {
    for(j in setdiff(1:k,i))
    {
      sig1[i,j]<-((nvec[i])*(nvec[j])*(nvec[i]-1)*(nvec[j]-1))/(2*((n-1)^2)*(n-3))
    }
  }
  mu<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in i:k)
    {
      mu[i,j]<-(((nvec[i])*(nvec[j]))/(n-1))
    }
  }
  muv<-t(mu)[lower.tri(t(mu))]
  bigsig<-matrix(0,(((k^2)-k)/2),(((k^2)-k)/2))
  for(i in 1:(k-1))
  {
    for(j in (i+1):k)
    {
      bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((i-1)*k)-((i*(i-1))/2)+j-i)]<- ((nvec[i]*nvec[j]*(nvec[i]-1)*(nvec[j]-1))/((n-1)*(n-3)))+(((nvec[i]*nvec[j])/(n-1))*(1-((nvec[i]*nvec[j])/(n-1))))
    }
  }
  if(k>=3)
  {
    for(i in 1:(k-2))
    {
      for(j in (i+1):(k-1))
      {
        for(l in (j+1):k)
        {
          bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-((nvec[i]*(nvec[i]-1)*nvec[j]*nvec[l])/((n-1)*(n-3))) - ((((nvec[i])^2)*nvec[j]*nvec[l])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-((nvec[i]*(nvec[i]-1)*nvec[j]*nvec[l])/((n-1)*(n-3))) - ((((nvec[i])^2)*nvec[j]*nvec[l])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-((nvec[j]*(nvec[j]-1)*nvec[i]*nvec[l])/((n-1)*(n-3))) - ((((nvec[j])^2)*nvec[i]*nvec[l])/((n-1)^2))
          bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-((nvec[j]*(nvec[j]-1)*nvec[i]*nvec[l])/((n-1)*(n-3))) - ((((nvec[j])^2)*nvec[i]*nvec[l])/((n-1)^2))
          bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-((nvec[l]*(nvec[l]-1)*nvec[j]*nvec[i])/((n-1)*(n-3))) - ((((nvec[l])^2)*nvec[j]*nvec[i])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-((nvec[l]*(nvec[l]-1)*nvec[j]*nvec[i])/((n-1)*(n-3))) - ((((nvec[l])^2)*nvec[j]*nvec[i])/((n-1)^2))
        }
      }
    }
  }
  if(k>=4)
  {
    for(i in 1:(k-3))
    {
      for(j in (i+1):(k-2))
      {
        for(l in (j+1):(k-1))
        {
          for(m in (l+1):k)
          {
            bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((l-1)*k)-((l*(l-1))/2)+m-l)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((j-1)*k)-((j*(j-1))/2)+m-j)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((i-1)*k)-((i*(i-1))/2)+m-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+m-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((j-1)*k)-((j*(j-1))/2)+m-j),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((l-1)*k)-((l*(l-1))/2)+m-l),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
          }
        }
      }
    }
  }
  return(list(as.matrix(mu1),sig1,as.matrix(muv),bigsig))
}


mhccexecutelong<-function(nvec,apmat)
{
  k<-length(nvec)
  n<-sum(nvec)
  smatch<-as.matrix(nonbimatch(distancematrix(as.matrix(dist(apmat, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))))$matches)
  multcm<-rep(0,k)
  cs<-c(0,cumsum(nvec))
  for(i in 1:k)
  {
    for(j in (cs[i]+1):(cs[i+1]))
    {
      multcm[i]<-multcm[i]+((as.numeric(smatch[j,3])>cs[i])&&(as.numeric(smatch[j,3])<=cs[i+1]))
    }
  }
  multcm<-multcm/2
  A<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in i:k)
    {
      for(l in (cs[i]+1):(cs[i+1]))
      {
        A[i,j]<-A[i,j]+((as.numeric(smatch[l,3])>cs[j])&&(as.numeric(smatch[l,3])<=cs[j+1]))
      }
    }
  }
  av<-t(A)[lower.tri(t(A))]
  return(list(as.matrix(multcm),as.matrix(av)))
}
library(mvtnorm)
#rho1 = 0.5
# rho2 = 0.5
n = 1000
theta = 7#0.2#1
Xtotal= c()
#set.seed(28)
nvec = rep(200,5)
set.seed(1234)
l = sample(100,24)
mu = rep(0,100)
diag_sigma = rep(1,100)
for(i in 1:length(nvec)){
  set.seed(i)
  #mu[l] = theta*i
  diag_sigma[l] = 1+(theta*i)
  Xtotal = rbind(Xtotal,rmvnorm(nvec[i],mean=mu,sigma = diag(diag_sigma)))
}

d = ncol(Xtotal)#100# d is total number of variables to start with
k = length(nvec)
ll<-mhcccreate(nvec)
mu1<-ll[[1]]
sig1<-ll[[2]]
muv<-ll[[3]]
bigsig<-ll[[4]]

pval_test<-function(X){
  X = as.matrix(X)
  lll<-mhccexecutelong(nvec,X)
  av<-lll[[2]]
  stbig<-t(as.matrix(av)-as.matrix(muv))%*%solve(bigsig)%*%(as.matrix(av)-as.matrix(muv))
  # print(as.numeric(stbig))
  stbig<-(1-pchisq(stbig,df=((k^2)-k)/2))*(d/ncol(X))
  #print(stbig)
  stsmall<-(sum(av)-sum(muv))/sum(bigsig)
  stsmall<-pnorm(stsmall)*(d/ncol(X))
  return(stbig)
}



hclust_fun2<-function(X,alpha){
  cor_X <- cor(X, use = "pairwise.complete.obs", method = "pearson")
  
  hclust_resu <- hclust(as.dist(1-cor_X),"single")#"single"
  
  l = list()
  l[[1]] =  abs(hclust_resu$merge[1,])
  for (i in 2:nrow(hclust_resu$merge)){
    temp1 = hclust_resu$merge[i,1]
    temp2 = hclust_resu$merge[i,2]
    if((temp1>0) & (temp2>0)){
      l[[length(l)+1]] = c(l[[temp1]],l[[temp2]])
    }
    if((temp1>0) & (temp2<0)){
      l[[length(l)+1]] = unique(c(l[[temp1]],-temp2))
    }
    
    if((temp1<0) & (temp2>0)){
      l[[length(l)+1]] = unique(c(l[[temp2]],-temp1))
    }
    if((temp1<0) & (temp2<0)){
      l[[length(l)+1]] = c(-temp2,-temp1)
    }  
  }
  
  
  merged_matr = hclust_resu$merge
  merged_list = list()
  for(i in 1:nrow(merged_matr)){
    merged_list[[i]] = list()
    merged_list[[i]][[1]] = merged_matr[i,]
    merged_list[[i]][[2]] = l[[i]]
    
  }
  
  
  attempt2<-function(n){
    
    r = list()
    a = merged_list[[n]][[1]][1]
    b = merged_list[[n]][[1]][2]
    
    if(a>0){
      a1 = sort(merged_list[[a]][[2]])
    }else{
      a1 = -a
    }
    if(b>0){
      b1 = sort(merged_list[[b]][[2]])
    }else{
      b1 = -b
    }
    
    c = sort(unique(c(a1,b1)))
    r[[1]] = c
    r[[4]] = as.numeric(pval_test(X[,c])<alpha)
    
    if(r[[4]]==0){#& length(r[[1]])>1
      if(length(r[[1]])==1){
        r[[2]] = "empty"
        r[[3]] = "empty"
        
      }
      else{
        r[[2]] = array(-1616)
        r[[3]] = array(-2828)
      }
      
    }
    else{
      if(a<0){
        r[[2]] = list()
        r[[2]][[1]] = -a
        r[[2]][[4]] = as.numeric(pval_test(X[,-a])<alpha)
        r[[2]][[2]] = "empty"#array()
        r[[2]][[3]] = "empty"#array()
        # r[[2]] = -a
        # return(r)
      }
      else{
        r[[2]] = attempt2(a)
      }
      
      if(b<0){
        r[[3]] = list()
        r[[3]][[1]] = -b
        r[[3]][[4]] = as.numeric(pval_test(X[,-b])<alpha)
        r[[3]][[2]] = "empty"#array()
        r[[3]][[3]] = "empty"#array()
      }
      else{
        r[[3]] = attempt2(b)
      }
      
    }
    
    
    
    return(r)
  }#ends attempt2()
  
  return(attempt2(nrow(merged_matr)))
  
}
tree_s<-hclust_fun2(Xtotal,alpha=0.05)

var1 = list()
pval_mmcc<-function(l,alpha){
  
  if(class(l[[2]])=="character"){
    if(l[[4]]==1){
      print(l[[1]])
      #print(l[[1]])
      assign("var1",list.append(var1,l[[1]]), envir = .GlobalEnv)
    }
    
  }
  else if(l[[4]]==1 & class(l[[2]])=="list"){
    if(l[[2]][[4]]==0 & l[[3]][[4]]==0){
      print(l[[1]])
      assign("var1",list.append(var1,l[[1]]), envir = .GlobalEnv)
      
    }
    else if(l[[2]][[4]]==1 & l[[3]][[4]]==0){
      pval_mmcc(l[[2]],alpha)
      
    }
    else if(l[[2]][[4]]==0 & l[[3]][[4]]==1){
      pval_mmcc(l[[3]],alpha)
      
    }
    else{
      pval_mmcc(l[[2]],alpha)     
      pval_mmcc(l[[3]],alpha)
      
      
    }
    
    
  }
  return(var1)
  
}

selected_var = sort(unlist(pval_mmcc(tree_s,0.05)))
true_var = sort(l)
intersect(true_var,selected_var)
