#Functions

library(entropy)


get_sample<-function(N,dist='normal',Sh=NaN,p1=NaN,p2=NaN,okgraph=FALSE){
    problem=TRUE
    if (okgraph) {x=c(0:1e4)/100}
    if (dist=='normal'){
        if (okgraph) {x=x-50}
        if (is.nan(p1)) {p1=0}
        if (! is.nan(Sh)) {
            p2=exp(Sh-log(sqrt(2*pi*exp(1))))
            
            problem=FALSE
        }
        else {
            if (! is.nan(p2)) {
                Sh=log(sqrt(2*pi*exp(1))*p2)
                
                problem=FALSE
            }
        }
        if (!problem) {
            print(paste('N(mean=',p1,',sd=',p2,'); Sh0=',Sh))
            sample=rnorm(N,mean=p1,sd=p2)
            if (okgraph) {
                pdf=dnorm(x,mean=p1,sd=p2)
                tit=paste('Normal(mean=',p1,'; sd=',p2,')')
            }
            
            #return(sample)
        }
        
    }
    if (dist=='exp') {
        if (! is.nan(Sh)) {
            p1=exp(1-Sh)
            
            problem=FALSE
        }
        else {
            if (! is.nan(p1)) {
                Sh=1-log(p1)
                
                problem=FALSE
            }
        }
        if (!problem) {
            print(paste('Exp(rate=',p1,'); Sh0=',Sh))
            sample=rexp(N,rate=p1)
            if (okgraph) {
                pdf=dexp(x,rate=p1)
                tit=paste('Exp(rate=',p1,')')
            }
            #return(sample)
        }
    }
    if (dist=='uniform') {
        if (! is.nan(Sh)) {
            if (is.nan(p1)) {p1=0}
            p2=exp(Sh)+p1
            
            problem=FALSE
        }
        else {
            if (! (is.nan(p1) & is.nan(p2))) {
                if (is.nan(p1)) {p1=0}
                if (is.nan(p2)) {p2=0}
                Sh=log(p2-p1)
                
                problem=FALSE
            }
        }
        if (!problem) {
            print(paste('Uniform(min=',p1,',max=',p2,'); Sh0=',Sh))
            sample=runif(N,min=p1,max=p2)
            if (okgraph) {
                pdf=dexp(x,mean=p1,sd=p2)
                tit=paste('Uniform(min=',p1,';max=',p2,')')
            }
            #return(sample)
        }
    }    
    if (dist=='gamma') {
        if (!is.nan(Sh)) {
            if (is.nan(p1)) {
                p1=1
                p2=exp(1-Sh)
            }
            else{
                aux=p1+log(gamma(p1))+(1-p1)*digamma(p1)
                p2=exp(aux-Sh)
            }
            
            
            
            problem=FALSE
        }
        else {
            if (!is.nan(p1)) {
                if (is.nan(p2)) {p2=1}
                Sh=p1-log(p2)+log(gamma(p1))+(1-p1)*digamma(p1)
                
                problem=FALSE
            }
        }
        
        if (!problem) {
            print(paste('Gamma(shape=',p1,',rate=',p2,'); Sh0=',Sh))
            sample=rgamma(N,shape=p1,rate=p2)
            if (okgraph) {
                pdf=dgamma(x,shape=p1,rate=p2)
                tit=paste('Gamma(shape=',p1,';rate=',p2,')')
            }
            #return(sample)
        }
    }
    if (problem) {
        print('Some kind of problem with:')
        print(paste('Distr: ',dist))
        print(paste('N=',N))
        print(paste('Sh=',Sh))
        print(paste('p1=',p1))
        print(paste('p2=',p2))
        return(NaN)
    }
    else {
        if (okgraph) {
            #plot(x,pdf,type='l',main=tit)
            hist(sample,breaks=10,main=tit,freq=F)
            #lines(pdf~x, col = 2, add = TRUE)
        }
        return(sample)
    }
}


set_bins <- function(base=dyadic,limit=1e4){
    if (base=='dyadic') {
        exponent=ceiling(log(limit)/log(2))
        sample=2^c(1:exponent)
    }
    else {
        if (base=='fib') {
            a=1
            b=1
            c=a+b
            sample=c(2)
            while (c<limit){
                a=b
                b=c
                c=a+b
                sample=append(sample,c)
            }
            
        }
        else {
            if (base=='factorial'){
                cont=2
                f=2
                sample=c(1,2)
                while (f<limit){
                    cont=cont+1
                    f=f*cont
                    sample=append(sample,f)
                }
            }
            else {
                exponent=ceiling(log(limit)/log(base))
                sample=base^c(1:exponent)
                sample=as.integer(sample)
                sample=as.integer(names(table(sample)))
            }
        }
    }
    return(sample)
}


plateau<-function(sec){
    #Finds a plateau at the end of the vector sec
    n=length(sec)
    lim=(max(sec)-min(sec))/n
    p=n
    for (i in (n:2)){
        if (abs(sec[i]-sec[i-1])<lim && abs(sec[i]-sec[n])<lim){p=i-1}
        else {break}
    }
    return(p)
}



evaluate<- function(ML,size,npts=5,D1=1,okplot=FALSE,title='') {
    ok=TRUE
    bottom=1
    if (okplot) { 
        plot(size,ML,xlab='-ln(bin size)',ylab='Sh estimate',main=title)
        abline(v=0)
        }
    
    plt=plateau(ML)
    if (length(ML)-plt>3) {
        hplt=mean(ML[plt:length(ML)])
        print(paste('plateau estimation',hplt))
        #print(length(ML)-plt+1)
        #print(length(size[plt:length(size)]))
        if (okplot) {lines(rep(hplt,length(ML)-plt+1)~size[plt:length(size)],col='green')}
        cero=min(abs(size))
        ptcero=1
        #print(paste(cero,'>',size[ptcero]))
        while (abs(size[ptcero])>cero && ptcero<length(size)) {
            #print(paste(cero,'>',size[ptcero]))
            ptcero=ptcero+1}
        #print(paste('Pt-0=',ptcero,'fin_pl=',plt))
        if (plt-ptcero<2) {
            print('plateau estimation')
            return(hplt)
            
        }
        else {
            ML=ML[1:plt]
            size=size[1:plt]
        }
    }
    top=length(ML)
    ntps_max=min(c(12,top-1))
    if (npts>ntps_max){npts=ntps_max}
    minmax=c(1,0,0)
    inf.a=NaN
    #print(ML)
    #print(size)
    for (n in c(ntps_max:npts)){
        #print((length(ML)-n))
        for (bottom in c((length(ML)-n):1)){
            top=bottom+n
            model=lm(ML[bottom:top]~size[bottom:top])
            
            a=summary(model)
            if (! is.data.frame(inf.a)){
                inf.a=data.frame(H0=a$coefficients[1],slope=a$coefficients[2],
                                 R2=a$adj.r.squared,F=a$fstatistic[1],
                                 bottom=bottom,top=top)
                minmax=c(abs(D1-inf.a$slope),inf.a$F,inf.a$R2)
            }
            if(!(is.nan(a$fstatistic[1]) | is.nan(a$coefficients[2]) | is.nan(a$adj.r.squared))){
                poll=0
                if (abs(D1-a$coefficients[2])<minmax[1]) {poll=poll+0.3}
                if (a$fstatistic[1]>minmax[2]) {poll=poll+0.2}
                if (abs(a$adj.r.squared)>minmax[3]) {poll=poll+0.2}
                if (poll>=0.5){
                    
                    
                    inf.a=data.frame(H0=a$coefficients[1],slope=a$coefficients[2],
                                     R2=a$adj.r.squared,F=a$fstatistic[1],
                                     bottom=bottom,top=top)
                    minmax=c(abs(D1-inf.a$slope),inf.a$F,inf.a$R2)
                    #print(inf.a)
                }
            }
            #v_par=rbind(v_par,inf.a)
            if (okplot){
                points(size,ML)
                lines(model$fitted.values~size[bottom:top], col='grey')
            }
            
        }
        
    }
    if (okplot) {
        bottom=inf.a$bottom
        top=inf.a$top
        model=lm(ML[bottom:top]~size[bottom:top])
        lines(model$fitted.values~size[bottom:top], col='red')
    } 
    print(inf.a)
    a=as.numeric(inf.a)
    return(a)
}



ebc_sample<- function(sample,method='MM',bins=set_bins('dyadic',1e4),okplot=FALSE,npts=5){
    
    #sample NaN values: these are converted to 0
    sample[is.na(sample)]=0.0
    
    size=-log((max(sample)-min(sample))/bins)
    v=c()
    for (i in bins) {
        tries=discretize(sample,i)
        v=append(v,entropy(tries,method=method,verbose=F))
        
    }
    tit=paste(method,' (N=',length(sample),')')
    
    a=evaluate(v,size,npts=npts,okplot=okplot,title=tit)
   return(a[1])
    
    
    
}


ebc_sample2d<- function(sx,sy,method='MM',bins=set_bins('dyadic',1e3),okplot=FALSE,npts=5){
    
    #sample NaN values: these are converted to 0
    sx[is.na(sx)]=0.0
    sy[is.na(sy)]=0.0
    R1=max(sx)-min(sx)
    R2=max(sy)-min(sy)
    
    size=-log((R1*R2)/(bins^2))
    v=c()
    for (i in bins) {
        tries=discretize2d(sx,sy,i,i)
        v=append(v,entropy(tries,method=method,verbose=F))
        
    }
    tit=paste(method,' (N=',length(sx),')')
    
    a=evaluate(v,size,npts=npts,D1=2,okplot=okplot,title=tit)
    return(a[1])
    
}

boot_sample<-function(sample,method='ML',bins=set_bins('dyadic',1e4),okplot=FALSE,npts=5,
                      try.size=NaN,try.number=100,try.replace=TRUE){
    if (is.nan(try.size)) {try.size=floor(0.75*length(sample))}
    if (okplot) {par(mfrow=c(3,3))}
    res=c()
    for (nt in c(1:try.number)){
        s=sample(sample,size=try.size,replace=try.replace)
        a=ebc_sample(s,method=method,bins=bins,okplot=okplot,npts=npts)
        #print(a)
        res=append(res,a[1])
    }
    if (okplot) {par(mfrow=c(1,1))}
    return(res)
}

position<-function(bool) {return(c(1:length(bool))[bool])}


ebc_points<- function(sample,method='MM',bins=set_bins('dyadic',1e4),okplot=FALSE,npts=5){
    
    
    size=-log((max(sample)-min(sample))/bins)
    v=c()
    for (i in bins) {
        tries=discretize(sample,i)
        v=append(v,entropy(tries,method=method,verbose=F))
        
    }
    tit=paste(method,' (N=',length(sample),')')

    #a=evaluate(v,size,npts=npts,plot=okplot,title=tit)
    return(data.frame(logN=size,ebc=v))
    
    
    
}


ebc_points2d<- function(sx,sy,method='MM',bins=set_bins('dyadic',5e2),okplot=FALSE,npts=5){
    
    R1=max(sx)-min(sx)
    R2=max(sy)-min(sy)
    
    size=-log((R1*R2)/(bins^2))
    v=c()
    for (i in bins) {
        tries=discretize2d(sx,sy,i,i)
        v=append(v,entropy(tries,method=method,verbose=F))
        
    }
    tit=paste(method,' (N=',length(sample),')')
    
    #a=evaluate(v,size,npts=npts,plot=okplot,title=tit)
    return(data.frame(logN=size,ebc=v))
    
    
    
}


explore_I<-function(func=expression(s),N=1e3,H_ref=c(-5,5),npts=20,okplot=FALSE
                    ,method='CS',ref_pl=c(-0.2,0.2),dist_s='normal',dist_e='normal')
    {
    library('rgl')
    library('akima')
    au=(1+sqrt(5))/2
    H0_ref=H_ref
    He_ref=H_ref
    dh=(H0_ref[2]-H0_ref[1])/npts
    H0_ref=c(0:npts)*dh+H0_ref[1]
    dh=(He_ref[2]-He_ref[1])/npts
    He_ref=c(0:npts)*dh+He_ref[1]
    #H_ref=c(-200:200)/20
    res=data.frame()
    par(mfrow=c(1,1))
    for (Shs in H0_ref){
        s=get_sample(N,dist=dist_s,Sh=Shs)
        Hs=ebc_sample(s,method=method)
        if (Shs<(-19)){
            print(Shs)
        }
        for (She in He_ref){
            #y=cos(cos((s+au)*(s-au)*(s^2-s+au)))+get_sample(N,dist='normal',Sh=She)
            if (She<(-19) & (abs(Shs) <=2 )){
                print(paste(Shs,She))
            }
            y=eval(func)
            y=y+get_sample(N,dist=dist_e,Sh=She)
            Hy=ebc_sample(y,method=method)
            Hsy=ebc_sample2d(y,s,method=method,okplot=F)
            print(paste('Hx=',Hs,'  Hy=',Hy,' Hxy=',Hsy))
  
            a=summary(lm(as.formula(paste('y~',as.character(func)))))  
            
            res=rbind(res,c(a$adj.r.squared,Shs,She,Hs,Hy,Hsy))
            
        }
    }
    names(res)=c('R2','H0x','H0e','Hx','Hy','Hxy')
    res$I=res$Hx+res$Hy-res$Hxy
    
    if (okplot){
        
        #Surface
        cols=c('blue','green','red','black','red','green','blue')
        co=cols[round(3*(res$R2)+4)]
        plot3d(cbind(res$H0x,res$H0e,res$I),type='s',radius=0.01,col=co,
               main=paste('y=',as.character(func),'+Error'),
               xlab = paste('H0x (',dist_s,')'),
               ylab=paste('H0e (',dist_e,')'),
               zlab='Ixy')
        b=interp(res$H0x,res$H0e,res$R2)
        a=interp(res$H0x,res$H0e,res$I)
        c=interp(res$H0x,res$H0e,rep(ref_pl[1],nrow(res)))
        d=interp(res$H0x,res$H0e,rep(ref_pl[2],nrow(res)))
        surface3d(c$x,c$y,c$z,col='grey',alpha=0.75)
        surface3d(d$x,d$y,d$z,col='grey',alpha=0.75)
        surface3d(a$x,a$y,a$z,col=round(abs(b$z)*3+1))
               
        #plot(res$I~res$R2, ylab='Hx+Hy-Hxy', xlab='adjusted R²')
        #abline(h=ref_pl)
        #a=(res$H0x-res$H0e)
        #plot(res$I~a,xlab='H0s-H0e',ylab='I',col=co,pch=18,
        #     main=paste('y=',as.character(func),'+Error'))
        #abline(h=ref_pl,v=0.0)
        #a=(res$Hy-res$Hx)
        #plot(res$I~a,xlab='Hy-Hx',ylab='I',col=co,pch=18,
        #     main=paste('y=',as.character(func),'+Error'))
        #abline(h=ref_pl,v=0.0)
        
        
        
        #Contour plot H0
        
        cH=as.vector(res$I)
        
        c=cbind(res$H0x,res$H0e)
        dv=(ceiling(max(c))-floor(min(c)))/npts
        v=min(c)+0:(npts-1)*dv
        minh=min(c)
        
        
        cx=round((c[,1]-minh)/dv)
        cy=round((c[,2]-minh)/dv)
        m=matrix(0.0,npts,npts)
        for (i in c(1:nrow(c))){
            m[cx[i],cy[i]]=cH[i]
        }
        k=min(m)+0:20*(max(m)-min(m))/20
        filled.contour(v,v,m,levels=k)
        title(main=paste('Hx+Hy-Hxy  (y=',as.character(func),'+error)'),
              xlab = paste('H0x (',dist_s,')'),
              ylab=paste('H0e (',dist_e,')'))
        #print(paste(max(m),min(m)))
        
        
    }
    return(res)
    
}

#END FUNCTIONS *****************************************************************************************************

