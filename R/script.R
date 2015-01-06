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
pars=c(beta,M, B, T)

v=z$vas
covs=c(2,3)  #### fix this! Has to be laser*time + other covs
x=z[covs]
optim(pars,contOrdEst, v=v, x=x)
