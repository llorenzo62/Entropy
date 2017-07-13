source('ebc.R')



#estimation one sample-100 all meth.

par(mfrow=c(3,3))
sample=get_sample(N=100,dist='normal',Sh=3.5,okgraph=T)

Sh=c()
for (met in c('ML','MM','Jeffreys','Laplace','SG','minimax','CS','shrink')){
    Sh=append(Sh,ebc_sample(sample,method=met, bins=set_bins('dyadic',1e5),okplot=T,npts=6))
}
summary(Sh)

#big simulation

global=data.frame()
#entropy estimation methods to test
methods= c('ML','MM','Jeffreys','Laplace','SG','minimax','CS','shrink')
#distributions to test
distributions=c('normal','exp','uniform')
#dyadic bins for algorithm
bins=set_bins('dyadic')
#sample size to test
ext=c(30,50,75,100,1000)
#entropies to test
Sh_ref=c(-0.88,0.01,1.577, 2.153, 3.3682, 3.7181)

par(mfrow=c(2,3))
par(oma=c(0,0,3,0))


for (dst in distributions){
    for (Sh in rep(Sh_ref,1)){
        
        for (N in rep(ext,1)){
            
            sample=get_sample(N,dist=dst,Sh=Sh)
            r=c()
            for (met in methods){
                r=append(r,ebc_sample(sample,bins=bins,method=met,okplot=F)[1])
            }
            results=data.frame(dist=dst,N=N,H0=Sh,rbind(r))
            global=rbind(global,results)
            
        }
    }
    boxplot(global[,4:ncol(global)]-global[,3],ylab='dSh')
}


names(global)=c('dist','N','H0',methods)
global_1=global

perf=global_1
v_met=c(4:(ncol(perf)))
perf[,v_met]=perf[,v_met]-perf$H0
perf[abs(perf$H0)>0.5,v_met]=perf[abs(perf$H0)>0.5,v_met]/perf[abs(perf$H0)>0.5,3]

table(perf$N,perf$H0)
table(perf[perf$dist=='normal',]$N,perf[perf$dist=='normal',]$H0)

#BIG repetition
# for (repetition in c(1:5)){
#     global=data.frame()
#     #entropy estimation methods to test
#     methods= c('ML','MM','Jeffreys','Laplace','SG','minimax','CS','shrink')
#     #distributions to test
#     distributions=c('normal','exp','uniform')
#     #dyadic bins for algorithm
#     bins=set_bins('dyadic')
#     #sample size to test
#     ext=c(30,50,75,100,1000)
#     #entropies to test
#     Sh_ref=c(-0.88,0.01,1.577, 2.153, 3.3682, 3.7181)
#     
#     par(mfrow=c(2,3))
#     par(oma=c(0,0,3,0))
#     #Simulation
#     
#     for (dst in distributions){
#         for (Sh in rep(Sh_ref,4)){
#             
#             for (N in rep(ext,4)){
#                 
#                 sample=get_sample(N,dist=dst,Sh=Sh)
#                 r=c()
#                 for (met in methods){
#                     r=append(r,ebc_sample(sample,bins=bins,method=met,okplot=F)[1])
#                 }
#                 results=data.frame(dist=dst,N=N,H0=Sh,rbind(r))
#                 global=rbind(global,results)
#                 #mtext(paste(dst,' H0=',Sh),outer=T)
#                 
#                 #abline(h=0)
#                 
#             }
#         }
#         boxplot(global[,4:ncol(global)]-global[,3],ylab='dSh')
#     }
#     
#     names(global)=c('dist','N','H0',methods)
#     global_1=rbind(global_1,global)
#     
#     perf=global_1
#     v_met=c(4:(ncol(perf)))
#     perf[,v_met]=perf[,v_met]-perf$H0
#     perf[abs(perf$H0)>0.5,v_met]=perf[abs(perf$H0)>0.5,v_met]/perf[abs(perf$H0)>0.5,3]
#     
#     table(perf$N,perf$H0)
#     table(perf[perf$dist=='normal',]$N,perf[perf$dist=='normal',]$H0)
# }
# table(perf$N,perf$H0)
# table(perf[perf$dist=='normal',]$N,perf[perf$dist=='normal',]$H0)

#END of BIG repetition


par(mfrow=c(1,1))

#names(perf)[v_met]=methods
plot(perf[,v_met])
boxplot(perf[,v_met],ylab='dSh')
abline(h=c(0.2,0.0,-0.2))

par(mfrow=c(2,2))
for (i in names(table(perf$dist))){
    draw=perf[perf$dist==i,]
    boxplot(draw[,v_met],ylab='dSh',main=paste('dist: ',i))
    #boxplot(draw$dSh/draw$Sh0~log10(draw$N),xlab='log10(N)',ylab='dSh/Sh0',main=paste('method: ',i))
    abline(h=c(0.2,0.0,-0.2))
}

par(mfrow=c(3,3))
for (i in v_met){
    
    #plot(log10(draw$N),draw$dSh/draw$Sh0,xlab='log10(N)',ylab='dSh/Sh0',main=paste('method: ',i))
    boxplot(perf[,i]~perf$H0,xlab='H0',ylab='dSh',main=paste('method: ',names(perf)[i]))
    abline(h=c(0.2,0.0,-0.2))
    #boxplot(perf[,i]~log10(perf$N),xlab='log10(N)',ylab='dSh',main=paste('method: ',names(perf)[i]))
    #abline(h=c(0.2,0.0,-0.2))
}
par(mfrow=c(3,3))
for (i in v_met){
    
    #plot(log10(draw$N),draw$dSh/draw$Sh0,xlab='log10(N)',ylab='dSh/Sh0',main=paste('method: ',i))
    #boxplot(perf[,i]~perf$H0,xlab='H0',ylab='dSh',main=paste('method: ',names(perf)[i]))
    #abline(h=c(0.2,0.0,-0.2))
    boxplot(perf[,i]~log10(perf$N),xlab='log10(N)',ylab='dSh',main=paste('method: ',names(perf)[i]))
    abline(h=c(0.2,0.0,-0.2))
}

v_min=c()
for (i in c(1:nrow(perf))) {
    a=min(abs(perf[i,v_met]))
    m=v_met[abs(perf[i,v_met])==a]
    v_min=append(v_min,m)
}    
table(v_min)
table(v_min)/sum(table(v_min))
select=v_met


v_perf=data.frame()
for (i in as.numeric(names(table(perf$N)))){
    draw=perf[perf$N==i,]
    v_min=c()
    for (j in c(1:nrow(draw))) {
        a=min(abs(draw[j,v_met]))
        m=v_met[abs(draw[j,v_met])==a]
        v_min=append(v_min,m)
    }    
    a=table(v_min)/sum(table(v_min))
    v_aux=c(rep(0,length(v_met)))
    for (k in c(as.numeric(names(a)))) {
        v_aux[(k-v_met[1]+1)]=a[names(a)==k]
    }
    aux=c(i,rbind(v_aux))
    v_perf=rbind(v_perf,aux)
}
names(v_perf)=c('N',names(perf)[v_met])
v_perf

par(mfrow=c(1,1))
for (i in (select-2)){ #ncol(v_perf))){
    
    if (i==(select[1]-2)) {
        
        plot(log10(v_perf$N),100*v_perf[,i],ylim=c(0,50),main='% minima',xlab='log10(N)',ylab='% minima',col=i-1,type='l')
    }
    else {
        
        lines(100*v_perf[,i]~log10(v_perf$N),col=i-1)
    }
}
legend(x='topright',legend=names(perf)[select],col=(select-3),lty=1)

#End big simulation ****************************************************************************************

#Mixing two samples iid


methods=c('MM','CS','shrink')
par(mfrow=c(4,3))
H0=3.5
dist1='uniform'
dist2=dist1
p1=50
base=1e5
factor=set_bins('fib',100)
res=data.frame()
sa=get_sample(base,dist1=dist,Sh=H0)
sb=get_sample(base,dist2=dist,Sh=H0)
for (f in append(0,factor)){
    #sa=get_sample(base,dist=dist,p1=(-f),Sh=H0)
    #sb=get_sample(base,dist=dist,p1=f,Sh=H0)
    sc=append((sa-f),(sa+f))
    #H=(log((f+1)/f)+H0)
    hist(sc,breaks=30,main=paste('factor=',f))
    H=H0+log(2)
    print(paste('H_esp=',H))
    v=c()
    for (met in methods){
        a=ebc_sample(sc,method=met,okplot=F)
        v=append(v,a[1])
    }
    res=rbind(res,data.frame(factor=f,H0=H,rbind(v)))
    
}
names(res)[3:5]=methods
res[,3:5]=res[,3:5]-res$H0
par(mfrow=c(1,1))
plot(res$factor,res$factor,type='l',col='white',ylim=c(min(res[,3:5]),(max(res[,3:5]))),xlab='factor',ylab='Sh-H0')
lines((res$MM)~res$factor,col=1)
lines((res$CS)~res$factor,col=2)
lines((res$shrink)~res$factor,col=3)
legend(x='bottomright',legend=names(res)[3:5],col=c(1,2,3),lty=1)
summary(res[,3:5])
res


methods=c('MM','CS')
par(mfrow=c(4,3))
H1=2.7
H2=2.7
dist1='uniform'
dist2='normal'
p11=5
p12=-10
base=1e5
tau=1
factor=set_bins('fib',100)
res=data.frame()
sa=get_sample((base*tau),dist=dist1,Sh=H1,p1=p11)
sb=get_sample(base,dist=dist2,Sh=H2,p1=p12)
theta=tau/(tau+1)
H=(tau*H1+H2)/(tau+1)+log((tau+1)/(tau^theta))
for (f in append(0,factor)){
    #sa=get_sample(base,dist=dist,p1=(-f),Sh=H0)
    #sb=get_sample(base,dist=dist,p1=f,Sh=H0)
    sc=append((sa-f),(sb+f))
    #H=(log((f+1)/f)+H0)
    hist(sc,breaks=30,main=paste('factor=',f))
    
    
    print(paste('H_esp=',H))
    v=c()
    for (met in methods){
        a=ebc_sample(sc,method=met,okplot=F)
        v=append(v,a[1])
    }
    res=rbind(res,data.frame(factor=f,H0=H,rbind(v)))
    
}
names(res)[3:ncol(res)]=methods
res[,3:ncol(res)]=res[,3:ncol(res)]-res$H0
par(mfrow=c(1,1))
plot(res$factor,res$factor,type='l',col='white',ylim=c(min(res[,3:ncol(res)]),(max(res[,3:ncol(res)]))),xlab='factor',ylab='Sh-H0')
#lines((res$MM)~res$factor,col=1)
lines((res$CS)~res$factor,col=2)
abline(h=0)
#lines((res$shrink)~res$factor,col=3)
legend(x='bottomright',legend=names(res)[3:ncol(res)],col=c(1,2,3),lty=1)
summary(res[,3:ncol(res)])
res


par(mfrow=c(3,2))
method='shrink'
sa=get_sample(1e2,dist='gamma',Sh=3.47)
ebc_sample(sa,okplot=T,method=method)
for (factor in c(1,2,10,100,1000)){
    print(ebc_sample(as.integer(factor*sa)/factor,okplot=T,method=method))
}



sn=get_sample(N=1000,dist='normal',Sh= 2.42639,p1=15)
sb=rbinom(1000,size=30,prob=0.5)
par(mfrow=c(2,1))
hist(sb, main='Binomial')
hist(sn, main='Normal')
par(mfrow=c(1,2))

met='ML'
par(mfrow=c(1,2))
ebc_sample(sb,method=met,okplot=T)
abline(v=0,h=log(sqrt(2*pi*exp(1)*30*0.5*0.5)))
mtext('Binomial(size=30,prob=0.5)')

ebc_sample(sn,method=met,okplot=T)
abline(v=0,h=log(sqrt(2*pi*exp(1)*30*0.5*0.5)))
abline(h=log(1000),lty=2 )
mtext('isoentropic Normal')



par(mfrow=c(2,2))
for (i in c(1,10,100,1000)){
    ebc_sample(ceiling(sn*i)/i,method=met,okplot=T)
    abline(v=0,h=log(sqrt(2*pi*exp(1)*30*0.5*0.5)))
    abline(h=log(1000),lty=2 )
    mtext(paste('ceiling(normal·',i,')/',i)) 
}

par(mfrow=c(1,1))
res=ebc_points(sb,method=met,okplot=T)
names(res)[2]='Binomial'
res$normal=ebc_points(sn,method=met,okplot=T)$ebc
for (i in c(1,10,100,1000)){
    res[,ncol(res)+1]=ebc_points(ceiling(sn*i)/i,method=met,okplot=T)$ebc
    names(res)[ncol(res)]=paste('ceiling(',i,')')
}
plot(a$logN,a$logN,type='l',ylim=c(min(res[,2:ncol(res)]),max(res[,2:ncol(res)])),col='white',
     xlab='-ln(bin size)',ylab='ebc estimate')
for ( i in c(2:3)){
    lines(res[,i]~res$logN,col=i)
}
for ( i in c(4:ncol(res))){
    points(res[,i]~res$logN,pch=i)
}
abline(v=0,h=log(sqrt(2*pi*exp(1)*30*0.5*0.5)))
#abline(h=log(1000),lty=2 )
legend(x='topleft',legend=names(res)[4:ncol(res)],pch=c(4:ncol(res)))
legend(x='bottomright',legend=names(res)[2:3],lty=1,col=c(2:3))








plot(res$logN,res$normal,type='p',pch=1,xlab='-ln(d_n)',ylab=paste(met,'estimate'))
abline(v=0,h=log(sqrt(2*pi*exp(1)*30*0.5*0.5)))
points(res$Binomial~res$logN,pch=2)
points(res[,4]~res$logN,pch=3)
points(res[,5]~res$logN,pch=4)
legend(x='topleft',legend=names(res)[2:5],pch=c(1:4))
########### END MIXING 2 SAMPLES IID*****************************************************************************

# BIG SIMULATION 2D VERSION**************************************************************************************

N=c(100,500,1000,3000,5000)
H_ref=c(0.1,1.0,2.0,3.0)
distributions=c('normal','exp','uniform')
methods=c('MM','CS','shrink')
met='CS'
res=data.frame()
par(mfrow=c(3,3))
for (big_i in c(1:300)){
    d1=distributions[round(runif(1,1,length(distributions)))]
    d2=distributions[round(runif(1,1,length(distributions)))]
    H1=H_ref[round(runif(1,1,length(H_ref)))]
    H2=H_ref[round(runif(1,1,length(H_ref)))]
    H0=H1+H2
    
    
    for (size in N){
        s1=get_sample(size,dist=d1,Sh=H1)
        s2=get_sample(size,dist=d2,Sh=H2)
        v=c(size)
        v=append(v,ebc_sample(s1,method=met))
        v=append(v,ebc_sample(s2,method=met))
        v=append(v,ebc_sample2d(s1,s2,method=met,okplot=F)-H0)
        res=rbind(res,rbind(v))
        print(v)
        
    }
}


res[abs(res$H1+res$H2)>0.5,4:ncol(res)]=res[abs(res$H1+res$H2)>0.5,4:ncol(res)]/(res$H1+res$H2)

names(res)=c('N','H1','H2','CS')
par(mfrow=c(1,1))
boxplot(res[4:ncol(res)],ylab='relative error estimate')
abline(h=c(-0.2,0.0,0.2))

library(rgl)
library(akima)
au=(1+sqrt(5))/2
N=10000
H_ref=c(-10:10)
method='CS'
res=data.frame()
par(mfrow=c(1,1))
for (Shs in H_ref){
    s=get_sample(N,dist='uniform',Sh=Shs)
    Hs=ebc_sample(s,method=method)
    H=c()
    for (She in H_ref){
        #y=cos(cos((s+au)*(s-au)*(s^2-s+au)))+get_sample(N,dist='normal',Sh=She)
        y=exp(cos((au*pi*s)))+get_sample(N,dist='normal',Sh=She)
        Hy=ebc_sample(y,method=method)
        Hsy=ebc_sample2d(y,s,method=method,okplot=F)
        print(paste('Hs=',Hs,'  Hy=',Hy,' Hsy=',Hsy))
        H=append(H,(Hs+Hy-Hsy))
        
        a=summary(lm(y~(s)))  
        
        res=rbind(res,c(a$adj.r.squared,Shs,She,Hs,Hy,Hsy))
        
    }
}
names(res)=c('R2','H0x','H0e','Hx','Hy','Hxy')
dH=-(res$Hxy-res$Hx-res$Hy)


co=round(3*abs(res$R2)+1)
plot3d(cbind(res[,c(2,3)],(dH)),type='s',radius=0.5,col=co,main='y=cos((x²+x+pi))+Error')
b=interp(res$H0x,res$H0e,res$R2)
a=interp(res$H0x,res$H0e,(dH))
c=interp(res$H0x,res$H0e,rep(0.2,nrow(res)))
d=interp(res$H0x,res$H0e,rep(-0.2,nrow(res)))
surface3d(c$x,c$y,c$z,col='grey',alpha=0.75)
surface3d(d$x,d$y,d$z,col='grey',alpha=0.75)
surface3d(a$x,a$y,a$z,col=round(abs(b$z)*3+1))

plot(dH~res$R2, ylab='Hx+Hy-Hxy', xlab='adjusted R²')

abline(h=c(-0.2,0.2))

