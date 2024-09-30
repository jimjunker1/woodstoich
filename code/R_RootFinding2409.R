print('file: R_RootFinding2409.R')
print(date())
# purpose:
# created by: Dieter.Wolf-Gladrow@awi.de
#    9/2024 version 1.0
# This software is provided 'as is' without warranty of
# any kind. But it's mine, so you can't sell it.
#
# ---------------------------------------------------------

# ,xaxs='i'
# ,yaxs='i' # apply exact plot limits

flag = 2

# 2 = DEB (9/2024)
# set.seed(1953) # set seed for random number generators
# M = 1e3 # number of Monte Carlo runs
# par(mfrow=c(1,2))
# install.packages('latex2exp') # do only once on your computer
# library(latex2exp)
# text(0,0.2,TeX('$\\mu = 0$'),col='black')
# title(ylab=TeX('$p_{\\mu_1,\\mu_2}(x)$'),line=2.5)
# text(0.3,100,bquote(~delta == .(delta)),col='black',pos=4,cex=1.5)
# install.packages('rootSolve')
library(rootSolve)
# ----------------------------------------------------------------------------------
if (flag == 0) {
  print('')
  # print(' ---------------------------------------------------')
  # print('')
  # print('Results: ')
  # print(c(round(,),''))
  # print(c(,''))
  # print(c('',))
}

# ----------------------------------------------------------------------------------
if (flag == 0) {
  print('')
  #  Q = read.table('.txt',header=TRUE) # Q = read.csv('.csv',header=TRUE)
  # library(latex2exp)
  # sflag = 1
  # if (sflag == 1) {
  # png('.png',width=16,height=16,units='cm',res=300)
  # plot(x,y,type='p',lwd=4,col='blue',xlab='x',ylab='y',las=1,cex=0.6,cex.lab=1.5)
  # dev.off()
  # }
  # print(' ---------------------------------------------------')
  # print('')
  # print('Results: ')
  # print(c(round(,),''))
  # print(c(,''))
  # print(c('',))
}
# ----------------------------------------------------------------
# Remarks:
# ----------------------------------------------------------------
# ----------------------------------------------------------------------------------
if (flag == 2) {
  print('2 = DEB (9/2024)')
  model <- function(x) {
    # x[1] = rdot; x[2] = JECC; x[3] = JENC
  kg = 0.8
  MV = 0.007
  kM = 0.25
  YVEc = 0.99
  YVEn = 0.01
  kappa = 0.58
  sigma = 0  # <- simpilfied for this case
  JAEcmax = 0.0006
  JAEnmax = 0.08
  MEnmax = 0.13
  MEcmax = 0.385
  MEc = 0.1
  MEn = 0.04
  # Please check especially equation 1 (F1)
  #          1.  2. 3.    4.          5.      5           43
  F1 <- 1/MV*(kg/(1+(YVEc*(kappa*x[2]-(1+sigma)*kM*MV/YVEc))^(-1)+
            #    6.              6
                 (YVEn*kappa*x[3])^(-1)+
            # 7.    8           9.      9           8
              (YVEc*(kappa*x[2]-(1+sigma)*kM*MV/YVEc)+
                     YVEn*kappa*x[3])))^(-1) - x[1] # = 0
            #                       721
  F2 <- MEc*(JAEcmax/MEcmax)-MV*x[1] - x[2]
  F3 <- MEn*(JAEnmax/MEnmax)-MV*x[1] - x[3]
  q = c(F1 = F1, F2 = F2, F3 = F3)
  return(q)}
  rdotSt = 0.5 # start values (guess)
  JEcGSt =  0.0008
  JEnGSt =  0.0004
  (ss <- multiroot(f = model, start = c(rdotSt,JEcGSt,JEnGSt)))
  #  Q = read.table('.txt',header=TRUE) # Q = read.csv('.csv',header=TRUE)
  # library(latex2exp)
  # sflag = 1
  # if (sflag == 1) {
  # png('.png',width=16,height=16,units='cm',res=300)
  # plot(x,y,type='p',lwd=4,col='blue',xlab='x',ylab='y',las=1,cex=0.6,cex.lab=1.5)
  # dev.off()
  # }
  # print(' ---------------------------------------------------')
  # print('')
  # print('Results: ')
  # print(c(round(,),''))
  # print(c(,''))
  # print(c('',))
}
# ----------------------------------------------------------------
# Remarks:
# ----------------------------------------------------------------
# ----------------------------------------------------------------------------------
if (flag == 1) {
  print('1 = example from package rootSolve')
  model <- function(x) {
      F1 <- x[1] + x[2] + x[3]^2 -12
      F2 <- x[1]^2 - x[2] + x[3] -2
      F3 <- 2*x[1] - x[2]^2 + x[3] -1
    q = c(F1 = F1, F2 = F2, F3 = F3)
    return(q)}
 #   + }
  # first solution
  # (ss <- multiroot(f = model, start = c(1, 1, 1)))
  # $root
  # [1] 1 2 3
  #
  # $f.root
  # F1            F2            F3
  # 2.902585e-10  4.794444e-09 -8.682278e-09
  #
  # $iter
  # [1] 6
  #
  # $estim.precis
  # [1] 4.588994e-09

  #  Q = read.table('.txt',header=TRUE) # Q = read.csv('.csv',header=TRUE)
  # library(latex2exp)
  # sflag = 1
  # if (sflag == 1) {
  # png('.png',width=16,height=16,units='cm',res=300)
  # plot(x,y,type='p',lwd=4,col='blue',xlab='x',ylab='y',las=1,cex=0.6,cex.lab=1.5)
  # dev.off()
  # }
  # print(' ---------------------------------------------------')
  # print('')
  # print('Results: ')
  # print(c(round(,),''))
  # print(c(,''))
  # print(c('',))
}
# ----------------------------------------------------------------
# Remarks:
# ----------------------------------------------------------------
# -------------------------------------------------------------
# sflag = 1
# if (sflag == 1) {
# png('.png',width=16,height=12,units='cm',res=300)
# plot(x,y,type='p',lwd=4,col='blue',xlab='x',ylab='y',las=1,cex=0.6,cex.lab=1.5)
# dev.off()
# }

# ylab = NA
# title(ylab=expression(paste('Vitamin ',B[12],' (pmol ',L^{-1},')')),line=2)


# jpeg('1710.jpeg',width=16,height=12,units='cm',quality=100,res=300)

# par(mfrow=c(2,2))  # allow for 2 x 2 panels in one plot

# title(ylab=expression(paste(NO[3],' (',mu,'mol ',L^-1,')')),line=2)

# pch=20: bullet,21: filled circle,22: filled square,
# 23: filled diamond,24: filled triangle point-up, 25: filled triangle point down.

# system('pwd')  # Mac print working directory
# getwd()        # PC

# bty='n',col.axis='white',xaxt='n',yaxt='n'  # = Visible off

# text(-2,0.35,expression(paste(hat(mu),' = 0')))
# text(x99,0.1,expression(x[99]),col='red')   # index
# text(xt,55,paste('= ',as.character(round(meansdest,3))),col='blue',pos=4)
# pos=4 -# text to the right of
# \u00B1  # +- symbol
# text(4, 8.4, 'expression(hat(beta) == (X^t * X)^{-1} * X^t * y)',cex = .8)

# This does not work:
# text(xt,55,expression(paste('= ',as.character(round(meansdest,3)))),col='blue')

# xcL = qnorm(alpha/2,mean=mu0,sd=sigma)
# xcU = qnorm(1-alpha/2,mean=mu0,sd=sigma)
# x1=xcU; x2=xmax; y1=dnorm(x1,mean=mu0,sd=sigma); y2=dnorm(x2,mean=mu0,sd=sigma);
# dx=(x2-x1)/50;
# xn=seq(x1,x2,dx); yn=dnorm(xn,mean=mu0,sd=sigma); xf=c(x2,x1,xn); yf=c(0,0,yn)
# quartz(title='I love my Mac',5,10)
# par(mfrow=c(3,1))
# plot(x,yA,type='l',lwd=4,col='black',xlab='(a)',ylab='Y',las=1,cex=0.6,cex.lab=1.5)
# polygon(xf,yf,col='blue')

# write.table(X1,file='xyz.txt')
# out = read.table('xyz.txt',header=TRUE)
# X1 = out$x

# install.packages('latex2exp') # do only once on your computer
# library(latex2exp)
# title(ylab=TeX('$p_{\\mu_1,\\mu_2}(x)$'))
# text(-3,0.36,bquote(~gamma[1] == .(gamma1) %+-% .(qr)),col='blue',cex=1.5,pos=4)
