

###############################
###    power of a matrix    ###
###  a is original matrix   ###
###  alpha is the power     ###
###  returns new matrix     ###
###############################
matpower <- function(a,alpha){
      small <- .00000001
if (length(c(a))==1) {return(a^(alpha))} else {
p1<-nrow(a)
eva<-eigen(a)$values
eve<-eigen(a)$vectors
eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
index<-(1:p1)[abs(eva)>small]
evai<-eva
evai[index]<-(eva[index])^(alpha)
ai<-eve%*%diag(evai,length(evai))%*%t(eve)
return(ai)}
}

############################################  
##  deviance between span(u) and span(v)  ##
############################################ 
dist<-function(u,v)
{
 dv<-proj(u)-proj(v)
 return(sqrt(sum((svd(dv)$d)^2)/p))
}

#######################################
##  make a basis matrix orthonormal  ##
#######################################
stv<-function(v)
{return(v%*%matpower(t(v)%*%v,-0.5))}


####################################################################
##  standardize n*p X to Z with zero mean and identity variance   ##
####################################################################
stand<-function(x){
n<-nrow(x)
p<-ncol(x)
xb <- apply(x, 2, mean)
xb <- t(matrix(xb, p, n))
x1 <- x - xb
sigma <- t(x1) %*% (x1)/(n-1)
eva <- eigen(sigma)$values
eve <- eigen(sigma)$vectors
sigmamrt <- eve %*% diag((1/sqrt(abs(eva)))*I(eva>1e-05)) %*% t(eve)
z <- sigmamrt %*% t(x1)
return(t(z))
}

#################################################################
##   projection matrix of v w.r.t. inner product (v,u)=v'S u   ## 
#################################################################
proj<-function(v,S)
{
 v<-as.matrix(v)
 p<-nrow(v)
 if (missing(S)) {S<-diag(1,p)}
 return(v%*%matpower(t(v)%*%S%*%v,-1)%*%t(v)%*%S)
}

#####################################
### discretizing y into H slices  ###
#####################################
slicing<-function(y,H) 
{
 n<-length(y)
 if (length(levels(as.factor(y)))>H)
  {
   ytilde<-rep(0,H+1)
   ytilde[1]<-min(y)
   for (h in 1:(H-1))
   {
    ytilde[h+1]<-quantile(y,h/H)
    ytilde[h+1]<-max(y[y<=ytilde[h+1]])
    if (ytilde[h+1]==ytilde[h]) ytilde[h+1]=min(y[y>ytilde[h]])
   }  
  }
 if (length(levels(as.factor(y)))<=H)
  {
   H <- length(levels(as.factor(y)))
   ytilde<-rep(0,H+1)
   ytilde[1]=min(y)
   for (h in 1:(H-1))
   {
    ytilde[h+1]<-min(y[y>ytilde[h]])
   }
  } 
 ytilde[H+1]=max(y)+1
 prop<-rep(1,H)
 for (i in 1:H)
  {
   prop[i] = sum((y >= ytilde[i])&(y < ytilde[i+1]))/n
  }
res<-list()
res$H<-H
res$ytilde<-ytilde
res$prop<-prop
return(res)
}


#################################################
##  kernel matrix for directional regression   ##
##  X is n*p matrix, Y is n-dim vector         ##
##  H is the number of slices                  ## 
##  fixed.vx = TRUE if cov(x) is treated known ##
##  vx is cov(x) and useful if fixed.vx=TRUE   ##
#################################################
mdr<-function(x,y,H,fixed.vx,vx){
  n<-nrow(x)
  p<-ncol(x)
  x<- t(t(x)-apply(x,2,mean))
  if (missing(fixed.vx)) {fixed.vx=FALSE}
  if (missing(vx)) {vx<-diag(1,p)}
 
  if (fixed.vx==TRUE) {mx<-vx} else
                     {mx<-cov(x)}
  imx<-matpower(mx,-1)

  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  ind<-matrix(0,n,H)

  xbar<-matrix(0,p,H)
   for (j in 1:(H-1))
   {
     ind[,j]<-((y >= ytilde[j])&(y < ytilde[j+1]))
     xbar[,j]<- (t(x)%*%(ind[,j]))/sum(ind[,j])
   }
  ind[,H]<-(y >= ytilde[H])
  xbar[,H]<- (t(x)%*%(ind[,H]))/sum(ind[,H])

  A<-matrix(0,p,p)
  B<-matrix(0,p,p)
  C<-0
  for (q in 1:H)
   {
    xq<-(t(x))[,ind[,q]==1]-xbar[,q]  
    A<-A + prop[q]*((xq%*%t(xq)/(sum(ind[,q])-1)+
           xbar[,q]%*%t(xbar[,q]))%*%imx%*%(xq%*%t(xq)/(sum(ind[,q])-1)+
           xbar[,q]%*%t(xbar[,q])) - mx)
    B<-B + sqrt(prop[j])*(xbar[,q]%*%t(xbar[,q]))
    C<-C + sqrt(prop[j])*(t(xbar[,q])%*%imx%*%xbar[,q])
   }
  C<-as.vector(C)
  M<-imx%*%(2*A + 2*(B%*%imx%*%B) + 2*B*C)%*%imx
  return(M)
}

library(optimx)

###############################################################
###  The following R codes for PWSVM are downloaded from    ###
### https://github.com/liyf1988/PrincipalWeightedSVM/blob/master/PrincipalWeightedSVM.r
### posted by Dr. Yufeng Liu.
###############################################################
#svm function that return optimized coefficients
#Input: x = covariates (continuous), y = outcome (binary, {-1,1})
#Output: optimized coefficients in svm
svm.f = function(x,y,weight,lambda) {
  n = length(y)
  x_bar = colMeans(x)
  x=t(t(x)-x_bar)
  sigma_x = t(x)%*%x/n;      # weight = weight; lambda = lambda
  svm.beta = function(beta, gs=NULL) {
    tmp = sum(sapply(1:n, function(i) ((y[i]==1)-y[i]*weight)*
              max(1-y[i]*(beta[1]+beta[-1]%*%(x[i,]-x_bar)),0)))
    Lambda = t(beta[-1])%*%sigma_x%*%beta[-1] + (lambda/n)*tmp
    return(Lambda[1,1])
  }
  svm.g = function(beta, gs=NULL) {
    tmp = matrix(NA,nrow=n,ncol=dim(x)[2]+1)
    for (i in 1:n) {
      tmp[i,1] = ((y[i]==1)-y[i]*weight)*(abs(1-y[i]*(beta[1]+
                        beta[-1]%*%(x[i,]-x_bar)))>0)*y[i]
      tmp[i,2:(dim(x)[2]+1)] = ((y[i]==1)-y[i]*weight)*
               (abs(1-y[i]*(beta[1]+beta[-1]%*%(x[i,]-x_bar)))>0)%*%
               (x[i,]-x_bar)*y[i]
    }
    Lambda_g = c(-(lambda/n)*sum(tmp[,1]),2*sigma_x%*%beta[-1]-
                 (lambda/n)*colSums(tmp[,2:(dim(x)[2]+1)]))
    return(Lambda_g)
  }
  svm.h = function(beta, gs=NULL) {
    return(rbind(rep(0,dim(x)[2]+1),cbind(rep(0,dim(x)[2]),2*sigma_x)))
  }
  svm_opt = optimx(rep(0,dim(x)[2]+1),fn=svm.beta)#,gr=svm.g,hess=svm.h)
  return(svm_opt)
}

#Perform svm using a vector of weights, and do pca on the coefficients
#Input: x = covariates (continuous), y = outcome (binary, {-1,1})
#output: principal component analysis of the coefficients in svm, 
#which contain the central space of y|x;
#svm coefficients estimated along a vector of weights from 0 to 1
pwsvm = function(x,y){
  p<-dim(x)[2]
  beta_vector = matrix(NA,nrow=20,ncol=p+1)
  w_vector = seq(from=0.3,to=0.7,length=10)
  M = matrix(0,dim(x)[2],dim(x)[2])
  for (i in 1:10) {
    coef = svm.f(x,y,w_vector[i],1)[2,1:(p+1)]
    beta_vector[i,] = unlist(coef)
    M = M + unlist(coef)[-1]%*%t(unlist(coef)[-1])
  }
  pca_M = prcomp(M,scale.=F,center=F)
  return(list("pc"=pca_M,"M"=M,"beta"=beta_vector))
}


########################################################
##              dual inverse regression               ##
##                                                    ##
##  X is n*p matrix, Y and T are n-dim vectors        ##
##  Hy is the number of slices for M_{Y|X}            ##
##  Ht is the number of slices for M_{T|X}            ##
##  d1 is the dimension of \beta_{Y1,T}               ##
##  d0 is the dimension of \beta_{Y0,T}               ##
##  methody is the inverse regression method for Y|X  ##
##  methodt is the inverse regression method for T|X  ##
##      -- options include SIR, SAVE, dirtnl. reg.    ##
##  penalty is by default FALSE,                      ##
##      -- interpretation is better if it is TRUE     ##
##  q1y is rank(M_{Y1|X}), specified if known         ##
##  q0y is rank(M_{Y0|X}), specified if known         ##
##  q1t & q0t are rank(M_{T|X}), specified if known   ##
##      -- allowed different for ease of computation  ##
##  fixe.vx means cov(X) is known (specified by vx)   ##
##  output:                                           ##
##    by1t is \beta_{Y1,T}, bty1 is \beta_{T,Y1}      ##
##    b0 is the better of by1t and bty1               ## 
##    the others follow likewise                      ##   
######################################################## 
dirm<-function(X=X,Y=Y,T=T,Hy=Hy,Ht=Ht,d1=d1,d0=d0,
               methody=methody,methodt=methodt,penalty=penalty,
               q1y=q1y,q0y=q0y,q1t=q1t,q0t=q0t,
               fixed.vx=fixed.vx,vx=vx)
{
 p<-ncol(X)
 mx<-cov(X)

 if ((missing(Hy))) {Hy <- p+1}
 if ((missing(Ht))) {Ht <- p+1}
 if (missing(penalty)) {penalty=FALSE}
 if (missing(q1y)) {q1y=p+1}
 if (missing(q0y)) {q0y=p+1}
 if (missing(q1t)) {q1t=p+1}
 if (missing(q0t)) {q0t=p+1}
 if (missing(fixed.vx)) {fixed.vx=FALSE}
 if (missing(vx)) {vx=mx}

 X1<-X[T==1,]
 Y1<-Y[T==1]
 X0<-X[T==0,]
 Y0<-Y[T==0]
  
 if ((methody=="sir")|(methody=="save"))
   {
    library(dr)
    if (d1>0) {my1<-dr(Y1~X1,method=methody,nslices=Hy)$M} else 
              {my1<-diag(0,p)}
    if (d0>0) {my0<-dr(Y0~X0,method=methody,nslices=Hy)$M} else
              {my0<-diag(0,p)}
   }
 if (methody=="dr")
   {
    if (d1>0) {my1<-mdr(X1,Y1,Hy,fixed.vx=fixed.vx,vx=vx)} else {
               my1<-diag(0,p)}
    if (d0>0) {my0<-mdr(X0,Y0,Hy,fixed.vx=fixed.vx,vx=vx)} else {
               my0<-diag(0,p)}
   }
 if (methody=="pwsvm")
   {
    S1<-2*Y1-1
    S0<-2*Y0-1
    if (d1>0) {my1<-pwsvm(X1,S1)$M} else {my1<-diag(0,p)}
    if (d0>0) {my0<-pwsvm(X0,S0)$M} else {my0<-diag(0,p)}
   }

 if ((methodt=="sir")|(methodt=="save"))
   {
    library(dr)
    if ((d1>0)|(d0>0)) {mt<-dr(T~X,method=methodt,nslices=Ht)$M} else
                       {mt<-diag(0,p)}
   }
 if (methodt=="dr")
   {  
    if ((d1>0)|(d0>0)) {mt<-mdr(X,T,Ht,fixed.vx=fixed.vx,vx=vx)} else {
                        mt<-diag(0,p)}
   }
 if (methodt=="pwsvm")
   {  
    S<-2*T-1
    if ((d1>0)|(d0>0)) {mt<-pwsvm(X,S)$M} else {mt<-diag(0,p)}
   }

 my1t<-my1%*%vx%*%mt
 my0t<-my0%*%vx%*%mt
 U1<-svd(my1)$u
 U0<-svd(my0)$u
 Ut<-svd(mt)$u

 if (methody=="sir") {py<- Hy} else {py<-p+1}
 if (methody=="sir") {pt<- Ht} else {pt<-p+1}

 a1<-sqrt(sum((svd(my1t)$d)^2))
 a0<-sqrt(sum((svd(my0t)$d)^2))
 f1y<-rep(0,py)
 if (d1>0)
 {
  f1y[1]<-a1
  for (q in 1:(py-1))
   {
    b1qy<-svd(proj(U1[,1:q],mx)%*%mt)$u[,1:d1]
    f1y[q+1]<-sqrt(sum((svd(my1%*%(mx-mx%*%proj(b1qy,mx))%*%mt)$d)^2))
   }
  if (penalty==TRUE) {f1y = f1y[-1]*st(f1y)} else {f1y<-f1y[-1]}
  if (q1y>py-1) {q1y<-which.min(f1y)}
  by1t<-svd(proj(U1[,1:q1y],mx)%*%mt)$u[,1:d1]
 } else {by1t<-rep(0,p)}
 
 f0y<-rep(0,py)
 if (d0>0)
 {
  f0y[1]<-a0
  for (q in 1:(py-1))
   {
    b0qy<-svd(proj(U0[,1:q],mx)%*%mt)$u[,1:d0]
    f0y[q+1]<-sqrt(sum((svd(my0%*%(mx-mx%*%proj(b0qy,mx))%*%mt)$d)^2))
   }
  if (penalty==TRUE) {f0y = f0y[-1]*st(f0y)} else {f0y<-f0y[-1]}
  if (q0y>py-1) {q0y<-which.min(f0y)}
  by0t<-svd(proj(U0[,1:q0y],mx)%*%mt)$u[,1:d0]
 } else {by0t<-rep(0,p)}

 f1t<-rep(0,pt)
 if (d1>0)
 {
  f1t[1]<-a1
  for (q in 1:(pt-1))
   {
    b1qt<-svd(proj(Ut[,1:q],mx)%*%my1)$u[,1:d1]
    f1t[q+1]<-sqrt(sum((svd(my1%*%(mx-mx%*%proj(b1qt,mx))%*%mt)$d)^2))
   }
  if (penalty==TRUE) {f1t = f1t[-1]*st(f1t)} else {f1t<-f1t[-1]}
  if (q1t>pt-1) {q1t<-which.min(f1t)}
  bty1<-svd(proj(Ut[,1:q1t],mx)%*%my1)$u[,1:d1]
 } else {bty1<-rep(0,p)}

 f0t<-rep(0,pt)
 if (d0>0)
 {
  f0t[1]<-a0
  for (q in 1:(pt-1))
   {
    b0qt<-svd(proj(Ut[,1:q],mx)%*%my0)$u[,1:d0]
    f0t[q+1]<-sqrt(sum((svd(my0%*%(mx-mx%*%proj(b0qt,mx))%*%mt)$d)^2))
   }
  if (penalty==TRUE) {f0t = f0t[-1]*st(f0t)} else {f0t<-f0t[-1]}
  if (q0t>pt-1) {q0t<-which.min(f0t)}
  bty0<-svd(proj(Ut[,1:q0t],mx)%*%my0)$u[,1:d0]
 } else {bty0<-rep(0,p)}

 if (f1y[q1y]<f1t[q1t]) {b1<-by1t} else {b1<-bty1}
 if (f0y[q0y]<f0t[q0t]) {b0<-by0t} else {b0<-bty0}

 res<-list()

 res$my1t<-my1t
 res$my0t<-my0t

 res$q1y <- q1y
 res$q0y <- q0y
 res$q1t <- q1t
 res$q0t <- q0t

 res$f1y<-f1y
 res$f0y<-f0y
 res$f1t<-f1t
 res$f0t<-f0t

 res$objy1t<-f1y[q1y+1]
 res$objy0t<-f0y[q0y+1]
 res$objty1<-f1t[q1t+1]
 res$objty0<-f0t[q0t+1]

 res$by1t <- by1t
 res$by0t <- by0t
 res$bty1 <- bty1
 res$bty0 <- bty0

 res$b1<-b1
 res$b0<-b0

 return(res)
}

################################################################
## ladle estimator for order determination of dual inv. reg.  ##
## X = n*p matrix, Y and T are n-dimensional vectors          ##
## Hy = number of slices of M_{Y|X}, Ht the same for M_{T|X}  ##
## nboot = number of bootstrap resamplings                    ##
## methody = method for M_{Y|X}, methodt the same for M_{T|X} ##
## output:                                                    ##  
##  d1 is the estimated dimensions for S_{Y1,T}               ##
##  d0 is the estimated dimensions for S_{Y0,T}               ##
################################################################ 
ladle.dirm<-function(X=X,Y=Y,T=T,Hy=Hy,Ht=Ht,nboot=nboot,
                     methody=methody,methodt=methodt,
                     stand.X=stand.X,fixed.vx=fixed.vx,vx=vx)
{
 p<-ncol(X)
 if (missing(stand.X)) {stand.X=TRUE}
 if (stand.X==TRUE) {X<-stand(X)}
 mx<-cov(X)
 imx<-matpower(mx,-1)
 mx1<-cov(X[T==1,])
 mx0<-cov(X[T==0,])
 
   if (p > 10) {r<-as.integer(p/log(p))} else {r <- p-1}
   obs<-cbind(Y,T,X)
   if ((missing(Hy))) {Hy <- p+1}
   if ((missing(Ht))) {Ht <- p+1}
   if (missing(nboot)) {nboot <- as.integer(n/2)}
   if (missing(fixed.vx)) {fixed.vx=FALSE}
   if (missing(vx)) {vx<-cov(X)}

 m<-dirm(X=X,Y=Y,T=T,Hy=Hy,Ht=Ht,d1=1,d0=1,
               methody=methody,methodt=methodt,q1y=1,q0y=1,q1t=1,q0t=1,
               fixed.vx=fixed.vx,vx=vx)

 my1t<-m$my1t
 my0t<-m$my0t
 Bhat1<-svd(my1t)$u
 Bhat0<-svd(my0t)$u
 lam1<-(svd(my1t)$d[1:(r+1)])^2
 lam0<-(svd(my0t)$d[1:(r+1)])^2

 obs<-cbind(Y,T,X)
 fn1<-rep(0,r+1)
 fn0<-rep(0,r+1)
  for (j in 1:nboot)
   {
    u<-round(runif(n,min=-0.5,max=n+0.5))
    bs<-obs[u,]    
    Xs<-bs[,-c(1,2)]
    Ys<-bs[,1]
    Ts<-bs[,2]

    ms<-dirm(X=Xs,Y=Ys,T=Ts,Hy=Hy,Ht=Ht,d1=1,d0=1,
                 methody=methody,methodt=methodt,q1y=1,q0y=1,q1t=1,q0t=1,
                 fixed.vx=fixed.vx,vx=mx)
    my1ts<-ms$my1t
    my0ts<-ms$my0t
    Bstar1<-svd(my1ts)$u
    Bstar0<-svd(my0ts)$u
    if (methody=="sir") {r<-min(r,Hy-1)}
    if (methodt=="sir") {r<-min(r,Ht-1)}
     for (i in 1:r)
     {
      fn1[i+1]<-fn1[i+1]+1-abs(det(t(Bstar1[,1:i])%*%Bhat1[,1:i]))
      fn0[i+1]<-fn0[i+1]+1-abs(det(t(Bstar0[,1:i])%*%Bhat0[,1:i]))
     }
  }
 fn1<-fn1/nboot
 fn0<-fn0/nboot
 
 res<-list()

 res$fn1 <- fn1
 res$fn0 <- fn0

 res$f1 <-fn1 / (1+sum(fn1))
 res$f0 <-fn0 / (1+sum(fn0))

 res$lam1 <- lam1
 res$lam0 <- lam0

 res$phi1 <- lam1/(1+sum(lam1))
 res$phi0 <- lam0/(1+sum(lam0))

 res$g1 <- res$f1 + res$phi1
 res$g0 <- res$f0 + res$phi0

 res$d1 <- which.min(res$g1) - 1
 res$d0 <- which.min(res$g0) - 1

return(res)
}

st<-function(v)
{
 k<-length(v)
 u<-v[-k]
 for (i in 1:(k-1)) {u[i]<-v[i+1]/min(v[1:i])}
 return(u)
}

#####################################################################
###  Ladle estimator for uniqueness of local minimum SDR space    ###
###  the working method in the paper to detect the global space   ###
###                                                               ###
###  X = n*p matrix, Y and T are n-dimensional vectors            ###
###  Hy = number of slices of M_{Y|X}, Ht the same for M_{T|X}    ###
###  nboot = number of bootstrap resamplings                      ###
###  d1 = dimension of S_{Y1,T}, d0 = dimension of S_{Y0,T}       ###
###  methody = method for M_{Y|X}, methodt the same for M_{T|X}   ###
###  q1y is rank(M_{Y1|X}), specified if known                    ###
###  q0y is rank(M_{Y0|X}), specified if known                    ###
###  q1t and q0t are rank(M_{T|X}), specified if known            ###
###     --- allowed different for convenience of implementation   ###
###  output:                                                      ###
###    d1 is the estimated dimension for P(by1t)-P(bty1)          ###
###     --- d1=0 means the global space exists                    ###
###    d0 follows likewsie                                        ### 
#####################################################################
ladle.glb<-function(X=X,Y=Y,T=T,Hy=Hy,Ht=Ht,nboot=nboot,d1=d1,d0=d0,
                methody=methody,methodt=methodt,
                q1y=q1y,q0y=q0y,q1t=q1t,q0t=q0t,
                fixed.vx=fixed.vx,vx=vx)
{
 p<-ncol(X)
 n<-nrow(X)
  if ((missing(Hy))) {Hy <- p+1}
  if ((missing(Ht))) {Ht <- p+1}
  if (missing(nboot)) {nboot <- as.integer(n/2)}
  if (missing(q1y)) 
       {q1y<-p+1}
  if (missing(q0y)) 
       {q0y<-p+1}
  if (missing(q1t)) 
       {q1t<-p+1}
  if (missing(q0t)) {q0t<-q1t}
  if (missing(fixed.vx)) {fixed.vx=FALSE}
  if (missing(vx)) {vx=cov(X)}

 b<-dirm(X=X,Y=Y,T=T,Hy=Hy,Ht=Ht,d1=d1,d0=d0,
                  methody=methody,methodt=methodt,penalty=TRUE, 
                  q1y=q1y,q0y=q0y,q1t=q1t,q0t=q0t,fixed.vx=fixed.vx,vx=vx)
 by1t<-b$by1t
 by0t<-b$by0t
 bty1<-b$bty1
 bty0<-b$bty0
 q1y<-b$q1y
 q0y<-b$q0y
 q1t<-b$q1t
 q0t<-b$q0t

 B1<-proj(bty1)-proj(by1t)
 B0<-proj(bty0)-proj(by0t)

 Bhat1<-svd(B1)$u
 Bhat0<-svd(B0)$u

 q1=2*d1+1
 q0=2*d0+1
 lam1<-(svd(B1)$d[1:q1])^2
 lam0<-(svd(B0)$d[1:q0])^2

 if ((lam1[1]<01e-06)&(lam1[1]<01e-06))
  {
    res<-list()   
    res$d1 <- 0
    res$d0 <- 0
  } else {

obs<-cbind(Y,T,X)
f1<-rep(0,q1)
f0<-rep(0,q0)
 for (j in 1:nboot)
  {
   u<-round(runif(n,min=-0.5,max=n+0.5))
   bs<-obs[u,]    
   Xs<-bs[,-c(1,2)]
   Ys<-bs[,1]
   Ts<-bs[,2]

   bs<-dirm(X=Xs,Y=Ys,T=Ts,Hy=Hy,Ht=Ht,d1=d1,d0=d0,
                  methody=methody,methodt=methodt,penalty=TRUE, 
                  q1y=q1y,q0y=q0y,q1t=q1t,q0t=q0t,fixed.vx=fixed.vx,vx=vx)
   by1ts<-bs$by1t
   by0ts<-bs$by0t
   bty1s<-bs$bty1
   bty0s<-bs$bty0

   B1s<-proj(bty1s)-proj(by1ts)
   B0s<-proj(bty0s)-proj(by0ts)

  Bstar1<-svd(B1s)$u
  Bstar0<-svd(B0s)$u
    for (i in 1:(q1-1))
     {
      f1[i+1]<-f1[i+1]+1-abs(det(t(Bstar1[,1:i])%*%Bhat1[,1:i]))
     }
    for (i in 1:(q0-1))
     {
      f0[i+1]<-f0[i+1]+1-abs(det(t(Bstar0[,1:i])%*%Bhat0[,1:i]))
     }
  }
 f1<-f1/nboot
 f0<-f0/nboot
 
 res<-list()

 res$f1 <- f1
 res$f0 <- f0

 res$ff1 <-f1 / (1+sum(f1))
 res$ff0 <-f0 / (1+sum(f0))

 res$lam1 <- lam1
 res$lam0 <- lam0

 res$phi1 <- lam1/(1+sum(lam1))
 res$phi0 <- lam0/(1+sum(lam0))

 res$g1 <- res$ff1 + res$phi1
 res$g0 <- res$ff0 + res$phi0

 res$d1 <- which.min(res$g1) - 1
 res$d0 <- which.min(res$g0) - 1
 }
 
return(res)
}



