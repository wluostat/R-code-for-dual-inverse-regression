source("ladle-prepare.R")
library(abind)
library(scatterD3)
library(plot3D)
library(energy)
library(dr)

################################################
###  simplified projection for orthonormal basis
################################################
proj<-function(v)
{return(v%*%matpower(t(v)%*%v,-1)%*%t(v))}

#### coefficients Ai and Bi
qs<-function(as,bs,uo,vo)
{
 if (vo>=1-0.0001) {return(diag(1,2))} else {
 c<-(1-uo)*vo/((as*bs-vo^2)*(uo))
 if (c==0) {return(diag(1,2))} else {
 aa<-max(min(2*c*sqrt(as*bs)/(1-2*c*vo),1),-1)
 a<-sin(asin(aa)/2)*sqrt((1-2*c*vo)/as)
 b<-c/a
 return(rbind(c(b,a),c(a,b))) }}
}

#### combine two matrices by rows
matdb<-function(u,v)
{
 u<-as.matrix(u)
 v<-as.matrix(v)
 m<-nrow(u)
 n<-nrow(v)
 rec<-numeric(0)
 for (i in 1:m)
 for (j in 1:n)
  {rec<-rbind(rec,c(u[i,],v[j,]))}
 return(rec)
}

#### list all the orthogonal matrices
mas<-function(s,ep)
{
 if (s==1)
  {
   w=2
   rec<-list()
   rec[[1]]=matrix(1,1,1)
   rec[[2]]=matrix(-1,1,1)
  } else {
 som<-sin(c(seq(from=0,to=pi/2,by=ep),pi/2)) #sqrt(seq(from=0,to=1,by=ep))
 m<-4*length(som)
 com<-c(sqrt(1-som^2))
 som<-c(rep(som,2),-rep(som,2))
 com<-rep(c(com,-com),2)
 ss<-s*(s-1)/2
  if (ss==1) {ms<- matrix(1:m,,1)} else {
              ms<- matdb((1:m),(1:m))}
  if (ss>2)
   {
    for (i in 3:ss) {ms<-matdb(ms,1:m)}
   }
 rec<-list()
 w<-1
  for (i in 1:(m^ss))
   {
    reci<-diag(1,s)
    for (j in 1:(s-1))
    for (k in (j+1):s)
     {
      rjk<-diag(1,s)
      jk<-(k-1)*(k-2)/2+j
      rjk[j,j]<-com[ms[i,jk]]
      rjk[k,k]<-com[ms[i,jk]]
      rjk[j,k]<-som[ms[i,jk]]
      rjk[k,j]<- -som[ms[i,jk]]
      rec[[w]]<-reci%*%rjk
      w<-w+1
     }
   }
  w<- w-1
  }
 res<-list()
 res$beta<-rec
 res$num<-w
 return(res)
}

#### list all the semi-orthogonal matrices (s>d)
masi<-function(q,ep)
{
  if (q==1)
   {
    beta<-list()
    beta[[1]]<-matrix(1,1,1)
    beta[[2]]<-matrix(-1,1,1) 
    res<-list()
    res$beta<-beta
    res$num<-2 
   } else {
    se<-c(seq(from=0,to=2*pi,by=ep),2*pi)
    mse<-matrix(se,,1)
     if (q>2)
      {
       for (i in 2:(q-1)) {mse<-matdb(mse,se)}
      }
    beta<-list()
     for (i in 1:nrow(mse))
      {
       bi<-1
        for (j in 1:(q-1))
         {bi<-bi*c(rep(cos(mse[i,j]),j),sin(mse[i,j]),rep(1,q-j-1))}
       beta[[i]]<-matrix(bi,,1)
      }
    res<-list()
    res$beta<-beta
    res$num<-nrow(mse) 
   }
 return(res)
}

mass<-function(s,d,ep)
{
 if (ep>=2*pi)
  {
   res<-list()
   beta<-list()
   beta[[1]]<-matrix(diag(1,s)[,1:d],s,d)
   beta[[2]]<-matrix(diag(-1,s)[,1:d],s,d)
   res$beta<-beta
   res$num<-2
  } 
 if (s==1)
  {
   res<-list()
   beta<-list()
   beta[[1]]<-matrix(1,1,1)
   beta[[2]]<-matrix(-1,1,1)
   res$beta<-beta
   res$num<-2
  } 
 if ((ep<2*pi)&(s>1))
  {
   if (d==1)
    {
     mm<-masi(s,ep)
     beta<-mm$beta
     w<-mm$num
    } else {
     mm<-masi(s,ep)
     beta<-mm$beta
     w<-mm$num
     for (i in 2:d)
      {
       mi<-masi(s-i+1,ep)
       nb<-list()
       nw<-0
        for (j in 1:w)
        for (k in 1:mi$num)
         {
          nw<-nw+1
          nb[[nw]]<-cbind(beta[[j]],
     matrix(svd(diag(1,s)-beta[[j]]%*%t(beta[[j]]))$u[,1:(d-i+2)],s,)%*%
     mi$beta[[k]])
         }
       beta<-nb
       w<-nw
      }
    }
   res<-list()
   res$beta<-beta
   res$num<-w
  }
 return(res)
}

#### tranfer into binary code
tentwo<-function(a,b)
{
 r<-rep(0,b)
 rs<-a
 for (i in 1:b)
  {
   r[i]<- rs-as.integer(rs/2)*2
   rs<-(rs-r[i])/2
  }
return(r)
}

cu<-function(v)
{
p<-length(v)
u<-v
for (i in 1:p)
{u[i]<-sum(v[1:i])}
return(u)
}

ryt<-function(Y,T,X,r,d,methody,methodt,slicey,slicet)
{
 Z<-stand(X)
 if ((methody=="sir")|(methody=="save"))
  {
   sy<-dr(Y~Z,method=methody,nslices=slicey)
   my<-sy$evectors%*%diag(sy$evalues)%*%t(sy$evectors)
  }
 if (methody=="dr") {my<-mhat_dr(Z,Y,slicey)}
 if ((methodt=="sir")|(methodt=="save"))
  {
   st<-dr(T~Z,method=methodt,nslices=slicet)
   mt<-st$evectors%*%diag(st$evalues)%*%t(st$evectors)
  }
 if (methodt=="dr") {mt<-mhat_dr(Z,T,slicet)}

 myt<-t(my)%*%mt
 s<-svd(myt)
 uy<-svd(my%*%s$u[,1:d])$u[,1:d]
 ut<-svd(mt%*%s$v[,1:d])$u[,1:d]
 syt<-svd(t(uy)%*%ut)
 if (r==0) {ryt<-t(uy)%*%ut} else {
 ryt<-syt$d[r]*diag(1,d-r+1)-
      (syt$u[,r:d])%*%diag(syt$d[r:d])%*%t(syt$v[,r:d]) }
 return(ryt%*%t(ryt))
}

#########  detecting multiplicity of singular values of Uy'Ut
agv_dir<-function(Y,T,X,r,d,k,rep,methody,methodt,slicey,slicet)
{
 Z<-stand(X)
 mhat<-ryt(Y,T,Z,r,d,methody,methodt,slicey,slicet)
 eva<-c(eigen(mhat)$values,0)
 ds<-ncol(mhat)
 n<-length(Y)
 p<-ncol(X)
 eve<-numeric(0)
 for (i in 1:rep)
 {
  Zs<-cbind(Z,matrix(rnorm(n*k),,k))
  ms<-ryt(Y,T,Zs,r,d,methody,methodt,slicey,slicet)
  eve<-eve+c(0,apply((eigen(ms)$vectors[-c(1:p),1:p])^2,1,mean)) 
 }
 eve<-eve/rep

 res<-list(0)
 res$vec <- eve
 res$val <- eva
 res$obj <- cu(eve) + eva/(1+cu(eva))
 res$d <- which.min(res$obj) - 1
 return(res)
}

#########  detecting q 
agv_q<-function(Y,T,X,k,rep,methody,methodt,slicey,slicet)
{
 Z<-stand(X)
 if ((methody=="sir")|(methody=="save"))
  {
   sy<-dr(Y~Z,method=methody,nslices=slicey)
   my<-sy$evectors%*%diag(sy$evalues)%*%t(sy$evectors)
  }
 if (methody=="dr") {my<-mhat_dr(Z,Y,slicey)}
 if ((methodt=="sir")|(methodt=="save"))
  {
   st<-dr(T~Z,method=methodt,nslices=slicet)
   mt<-st$evectors%*%diag(st$evalues)%*%t(st$evectors)
  }
 if (methodt=="dr") {mt<-mhat_dr(Z,T,slicet)}
 eva<-c(svd(cbind(my,mt))$d,0)

 n<-length(Y)
 p<-ncol(X)
 eve<-rep(0,p+1)
 for (i in 1:rep)
 {
  Zs<-cbind(Z,matrix(rnorm(n*k),,k))
  if ((methody=="sir")|(methody=="save"))
  {
   sy<-dr(Y~Zs,method=methody,nslices=slicey)
   my<-sy$evectors%*%diag(sy$evalues)%*%t(sy$evectors)
  }
  if (methody=="dr") {my<-mhat_dr(Zs,Y,slicey)}
  if ((methodt=="sir")|(methodt=="save"))
  {
   st<-dr(T~Zs,method=methodt,nslices=slicet)
   mt<-st$evectors%*%diag(st$evalues)%*%t(st$evectors)
  }
  if (methodt=="dr") {mt<-mhat_dr(Zs,T,slicet)}
  rea<-matrix((svd(cbind(my,mt))$u[-c(1:p),1:p])^2,k,)
  eve<-eve+c(0,apply(rea,2,mean)) 
 }
 eve<-eve/rep

 res<-list(0)
 res$vec <- eve
 res$val <- eva
 res$obj <- cu(eve) + eva/(1+cu(eva))
 res$d <- which.min(res$obj) - 1
 return(res)
}

##### switch columns
swc<-function(m,p,d)
{
if (is.vector(m)==TRUE) {m<-matrix(m,1,)}
sm<-m
return(m[,c(t(matrix(1:(p*d),p,d)))])
}
##### vec to mat by rows






mae<-function(q,ep)
{
 m<-list()
 w<-2
 m[[1]]<-matrix(0,q,q)
 m[[2]]<-diag(1,q)
 if (q>1)
 {
  for (k in 1:(q-1))
   {
    ak<-mass(q,k,ep)
    for (i in 1:(ak$num))
    {
     w<-w+1
     aki<-ak$beta[[i]]
     m[[w]]<-aki%*%t(aki)
    }
   }
 }
 res<-list()
 res$m<-m
 res$num<-w
 return(res)
}

############################################
## exhaustively list all possible values of 
## a row below in a b_{-A}(b_A)\inv
############################################
es<-function(q,ep)
{
 rec<-matrix(0,1,q)
 a<-as.integer(1/ep)
 for (i in 1:a)
  {
   ri<-numeric(0)
   for (j in 1:q)
    {
     aj<-rec+rep(ep,nrow(rec))%*%t(diag(1,q)[,j])
     ri<-rbind(ri,aj)
    }
   rec<-rbind(rec,ri)
  }
 res<-matrix(rep(c(-1,1),each=2^(q-1)),,1)
 if (q>1)
  {for (j in 2:q) {res<-cbind(res,rep(rep(c(-1,1),2^{j-1}),each=2^(q-j)))}}
 re<-numeric(0)
 for (i in 1:nrow(rec))
 for (j in 1:nrow(res))
  {re<-rbind(re,sqrt(rec[i,])*res[j,])}
 return(re)
}

#############################
## list all D's in case(iii)
#############################
Ds<-function(q,ep)
{
 s<-seq(from=0,to=0.5,by=ep/2)
 fs<-c(0.5-sqrt(0.25-s^2),0.5+sqrt(0.25-s^2),seq(from=0,to=1,by=ep))
 rec<-matrix(fs[order(fs)],,1)
 #matrix(seq(from=0,to=1,by=ep),,1)
 if (q>1)
  {
   for (i in 2:q)
    {
     rec<-matdb(rec,matrix(seq(from=0,to=1,by=ep),,1))
     rec<-rec[rec[,i-1]>=rec[,i],]
    }
  }
 return(rec)
}

##################################################################
###  main function: my=M_Y, mt=M_T, dy=rank(my), dt=rank(mt),  ###
###                 q,d,r same as in the paper, ep=epsilon     ###
################################################################## 
spnew<-function(my,mt,dy,dt,q,d,r,ep)
{
 p<-nrow(my)
 omy<-my
 omt<-mt
# s<-svd(t(my)%*%mt)
# Ut<-matrix(svd(mt%*%s$v[,1:d])$u,p,d)
# Uy<-matrix(svd(my%*%s$u[,1:d])$u,p,d)

 s<-svd(t(svd(my)$u[,1:dy])%*%svd(mt)$u[,1:dt])
 Ut<-matrix(svd(matrix(svd(mt)$u[,1:dt],,dt)%*%matrix(s$v[,1:d],,d))$u,p,d)
 Uy<-matrix(svd(matrix(svd(my)$u[,1:dy],,dy)%*%matrix(s$u[,1:d],,d))$u,p,d)
 U<-matrix(svd(cbind(my,mt))$u[,1:q],p,q)
 if (q<p) {V<-matrix(svd(diag(1,p)-U%*%t(U))$u[,1:(p-q)],p,p-q)} else {
           V<-matrix(0,p,0)}  
 if (r>0) {Uc<-matrix(Uy%*%svd(t(Uy)%*%Ut)$u[,1:r],p,r)} else {
           Uc<-matrix(0,p,0)}
 if (d>r)
    {
     if (r>0)
      {
       uy<-matrix(Uy%*%svd(t(Uy)%*%Ut)$u[,-(1:r)],p,d-r)
       Ut<-matrix(Ut%*%svd(t(Uy)%*%Ut)$v[,-(1:r)],p,d-r)
       Uy<-uy
      } 
    } else {
     bf<-list()
     bf$betas<-Uc
     w<-1
    }

  if (d>r)
   {
    nc<-cbind(Ut,svd(Uy-proj(Ut)%*%Uy)$u[,1:(d-r)])   
    B<-t(nc)%*%Ut
    Bs<-t(nc)%*%(svd(Uy-proj(Ut)%*%Uy)$u[,1:(d-r)])
    A<-t(nc)%*%Uy
    A1<-A[1:(d-r),]
    A2<-A[-(1:(d-r)),]
   }

 if (p==q)
  {
   mw<-mae(d-r,ep)
   bf<-list()
   w<-0
   for (i in 1:(mw$num))
     {
      w<-w+1
      bi<-rbind(A1,(mw$m[[i]])%*%A2)
      bf[[w]]<-svd(cbind(Uc,nc%*%bi))$u[,1:d]
     }
   }

 if (q<p)
  {
   if ((p-q)>=(d-r))
   {
    w<-0
    mw<-mas(d-r,ep)
    Dw<-Ds(d-r,ep)
    if ((p-q)==(d-r)) {mv<-mas(p-q,ep)}
    if ((p-q)>(d-r))  {mv<-mass(p-q,d-r,ep)}
    mu<-numeric(0)
    bf<-list()
    for (i in 1:mw$num)
    for (j in 1:nrow(Dw))
    for (k in 1:mv$num)
     {
      w<-w+1
      muij<-(mw$beta[[i]])%*%diag(Dw[j,],d-r)%*%t(mw$beta[[i]])
      Dsj<-diag(sqrt(abs(Dw[j,]-(Dw[j,])^2)),d-r)
      bijk<-cbind(Uc,nc%*%(B%*%A1 + Bs%*%muij%*%A2) + 
            V%*%(mv$beta[[k]])%*%Dsj%*%t(mw$beta[[i]])%*%A2)
      bf[[w]]<-svd(bijk)$u[,1:d]
     }
   } else {
    w<-0
    mw<-mas(d-r,ep)
    Dw<-Ds(p-q,ep)
    for (l in 1:(d-r-(p-q))) 
     {Dw<-matdb(Dw,rep(rep(c(0,1),each=2^(d-r-(p-q)-l)),2^(l-1)))}
    mv<-mas(p-q,ep)
    bf<-list()
    for (i in 1:mw$num)
    for (j in 1:nrow(Dw))
    for (k in 1:mv$num)
     {
      w<-w+1
      muij<-(mw$beta[[i]])%*%diag(Dw[j,],d-r)%*%t(mw$beta[[i]])
      Dsj<-diag(sqrt(abs(Dw[j,]-(Dw[j,])^2)),d-r)
      bijk<-cbind(Uc,nc%*%(B%*%A1 + Bs%*%muij%*%A2) + 
            V%*%cbind(mv$beta[[k]],matrix(0,p-q,d-r-(p-q)))%*%
            Dsj%*%t(mw$beta[[i]])%*%A2)
      bf[[w]]<-svd(bijk)$u[,1:d]
     }
   } 
 }
 rec<-rep(0,w)
  for (i in 1:w) 
    {rec[i]<-sqrt(sum((t(my)%*%bf[[i]]%*%t(bf[[i]])%*%mt-t(my)%*%mt)^2))/
             sqrt(sum((t(my)%*%mt)^2))}
 res<-list()
 res$beta<-bf
 res$num<-w
 res$score<-rec
 return(res)
}

############################################# ##########
###  approximation of each estimate to solving (5)   ###
###  assuming the true value of M_Y and M_T known    ###
######################################################## 
score<-function(rec,tmy,tmt,sig)
{
 w<-rec$num
 bf<-rec$beta
 res<-rep(0,w)
 for (i in 1:w) 
   {
    pi<-sig%*%bf[[i]]%*%matpower(t(bf[[i]])%*%sig%*%bf[[i]],-1)%*%t(bf[[i]])
    res[i]<-sqrt(sum((t(tmy)%*%pi%*%sig%*%tmt-t(tmy)%*%sig%*%tmt)^2))/
             sqrt(sum((t(tmy)%*%sig%*%tmt)^2))
   }
 return(res)
}



