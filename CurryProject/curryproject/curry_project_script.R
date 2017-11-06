setwd("/Users/stephenau/Documents/CurryProject") 
install.packages("jpeg")
library("jpeg")
par(mai=c(0,0,0,0))
graphics.off()
shots<-readJPEG("NBACOURT.jpg",native=TRUE)

  res = dim(shots)[1:2] # get the resolution
    plot(1,1,xlim=c(1,res[1]),ylim=c(1,res[2]),asp=1,type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  rasterImage(shots,1,1,res[1],res[2])
  

stephcurry<-locator()
kobebryant<-locator()
kylekorver<-locator()
klaythompson<-locator()
text(stephcurry,label=1:length(stephcurry$x),col="gold1")
text(kobebryant,label=1:length(kobebryant$x),col="purple1")
text(kylekorver,label=1:length(kylekorver$x),col="red1")
text(klaythompson,label=1:length(klaythompson$x),col="blue1")

stephthrees<-data.frame(stephcurry)
kobethrees<-data.frame(kobebryant)
korverthrees<-data.frame(kylekorver)
thompsonthrees<-data.frame(klaythompson)

library('spatstat')
source('quadrat_functions.r')


graphics.off()

curry<-read.csv("stephsthrees.csv")
kobe<-read.csv("kobesthrees.csv")
korver<-read.csv("korversthrees.csv")
thompson<-read.csv("thompsonsthrees.csv")

pp_curry<-data.frame(x=curry$x,y=curry$y)
quadrat(pp_curry,20,20,TRUE,'CURRY')#the negative binomial p-value is the only process that produces a value above 0 so there is small evidence that Stephen Currys shot chart is clustered(0.1%) under this model.
pp_kobe<-data.frame(x=kobe$x,y=kobe$y)
quadrat(pp_kobe,20,20,TRUE,'KOBE')#p-values show we cna reject all these processes assuming a cut-off value at 0.05. The poisson process produces a p-value greater than 0 so there is small chance of CSR under this model
pp_korver<-data.frame(x=korver$x,y=korver$y)
quadrat(pp_korver,20,20,TRUE,'KORVER')#negative binomial p-value is 0.159, greater than 0.05 so there is evidence of clustering for korver's shot chart
pp_thompson<-data.frame(x=thompson$x,y=thompson$y)
quadrat(pp_thompson,20,20,TRUE,'THOMPSON')#p-values reveal we can reject all three processes under the quadrat analysis. the negative binomial is the only process that produces a p-value greaterthan 0
#these lines serve to perform quadrat analysis against statistical models that suggest clustering, random, and normal processes.

bndry<-owin(xrange=c(0,599),yrange=c(0,800))
help(ppp)
curry.pp<-ppp(curry$x,curry$y,window=bndry)
kobe.pp<-ppp(kobe$x,kobe$y,window=bndry)
korver.pp<-ppp(korver$x,korver$y,window=bndry)
thompson.pp<-ppp(thompson$x,thompson$y,window=bndry)
#This stores point patterns processes into variables that are later used to plot

layout(matrix(1:4,2,2))
par(mai=rep(0.5,4))
plot(curry.pp,cex=0.5,main="CURRY")
plot(kobe.pp,meansize=2,main="KOBE")
plot(korver.pp,cex=0.5,main="KORVER")
plot(thompson.pp,cex=0.5,main="THOMPSON")


layout(matrix(1:4,2,2))
par(mai=rep(0.5,4))
plot(density.ppp(curry.pp,50))
plot(density.ppp(kobe.pp,50))
plot(density.ppp(korver.pp,50))
plot(density.ppp(thompson.pp,50))

GSC<-Gest(curry.pp)
KSC<-Kest(curry.pp)

GKB<-Gest(kobe.pp)
KKB<-Kest(kobe.pp)

GKK<-Gest(korver.pp)
KKK<-Kest(korver.pp)

GKT<-Gest(thompson.pp)
KKT<-Kest(thompson.pp)

layout(matrix(c(1,2),1,2))
plot(GSC,cbind(theo,rs)~r,col=c("black","red"),ylab="G(r)",main="StephenCurry shotchart",lty=1,legendargs=list(cex=0.75),cex.main=0.8)
plot(KSC,sqrt(cbind(theo,border)/pi)-r~r,xlim=c(0,20),col=c("black","red"),ylab="L(r)",main="StephenCurry shotchart",lty=1,legendargs=list(cex=0.75),cex.main=0.8)

plot(GKB,cbind(theo,rs)~r,col=c("black","red"),ylab="G(r)",main="KOBE shotchart",lty=1,legendargs=list(cex=0.75),cex.main=0.8)
plot(KKB,sqrt(cbind(theo,border)/pi)-r~r,xlim=c(0,20),col=c("black","red"),ylab="L(r)",main="KOBE shotchart",lty=1,legendargs=list(cex=0.75),cex.main=0.8)

plot(GKK,cbind(theo,rs)~r,col=c("black","red"),ylab="G(r)",main="Kylekorver shotchart",lty=1,legendargs=list(cex=0.75),cex.main=0.8)
plot(KKK,sqrt(cbind(theo,border)/pi)-r~r,xlim=c(0,20),col=c("black","red"),ylab="L(r)",main="Kylekorver shotchart",lty=1,legendargs=list(cex=0.75),cex.main=0.8)

plot(GKT,cbind(theo,rs)~r,col=c("black","red"),ylab="G(r)",main="KlayThompson shotchart",lty=1,legendargs=list(cex=0.75),cex.main=0.8)
plot(KKT,sqrt(cbind(theo,border)/pi)-r~r,xlim=c(0,20),col=c("black","red"),ylab="L(r)",main="KlayThompson shotchart",lty=1,legendargs=list(cex=0.75),cex.main=0.8)

myKfun <- function(S, ..., lam) { Kinhom(S, lambda=lam[S], ...) }
GSCE<-envelope(curry.pp,Gest)
KSCE<-envelope(curry.pp,Kest)
layout(matrix(c(1,2),1,2))
plot(GSCE,cbind(theo,obs,lo,hi)~r,col=c("black","red","black","black"),lty=c(1,1,2,2),ylab="G(r)",main="Curryshots: G(r)",legendargs=list(cex=0.75),cex.main=0.8)
plot(KSCE,sqrt(cbind(theo,obs,lo,hi)/pi)-r~r,col=c("black","red","black","black"),lty=c(1,1,2,2),ylab="L(K(r))",main="Curryshots: K(r)",legendargs=list(cex=0.75),cex.main=0.8)

curryVS<-read.csv('stephcurrycomposite.csv')
attach(curryVS)
names(curryVS)
pairs(cbind(x,y,yKlay),main="StephvsKlay")
cor(cbind(x,y,xKlay,yKlay))
