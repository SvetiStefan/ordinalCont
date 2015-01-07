#z=NULL
#for (i in 1:84){
#for (j in 1:3){
#z=rbind(z, c(i,j,y[i,j],laser[i]))
#}}
#z=as.data.frame(z)
#names(z)=c("ID","time","vas","laser")
#save(z,file="data/pain.RData")



load('data/pain.RData')
source('R/function.R')

beta=c(-1.90, 1.74, 1.01)
M = 0.63
B = 0.43
T = 0.85
pars=c(beta, M, B, T)

v=z$vas/10
v[v==0]=0.001
v[v==1]=0.999
#covs=c(2,3)  #### fix this! Has to be laser*time + other covs
#x=z[covs]
xx=z$time*z$laser
x1=sapply(xx==1,ifelse,1,0)
x2=sapply(xx==2,ifelse,1,0)
x3=sapply(xx==3,ifelse,1,0)
x=cbind(x1,x2,x3)

