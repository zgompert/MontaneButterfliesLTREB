## R script to update prediction plots with observed data
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales)
library(fields)

## helper functions
stand<-function(x=NA){
	x<-(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
	return(x)
}

standNew<-function(x=NA,newx){
	newx<-(newx-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
	return(newx)
}



## read in the data
dat<-read.csv("foristerCastData2022.csv")
ndat<-read.csv("foristerCastData2022_April.csv")


#  Castle Peak   Donner Pass Lang Crossing Sierra Valley    Washington 
#        41890        103477        102832        101514        108523

sites<-unique(dat$site_name)
myargs<-commandArgs(trailingOnly=TRUE)
j<-as.numeric(myargs[1])

load(paste("forecast2022_",sites[j],".rdat",sep=""))
cat("working on site",sites[j],"\n")
sub_dat<-dat[dat$site_name==sites[j],]
sub_ndat<-ndat[ndat$site_name==sites[j],]
spKeep<-names(which(tapply(sub_dat$pa,INDEX=sub_dat$genus_species,sum) > 10))
sub_dat<-sub_dat[ (sub_dat$genus_species %in% spKeep),]
sp<-unique(sub_dat$genus_species)
bb<-as.matrix(cbind(stand(sub_dat$PrevFall_tmn),stand(sub_dat$PrevFall_ppt),stand(sub_dat$Win_tmn),stand(sub_dat$pck),stand(sub_dat$ordDate),stand(sub_dat$ordDate)^2,stand(sub_dat$Year)))
bbi<-as.matrix(cbind(bb,bb[,4]*bb[,5],bb[,4]*bb[,6]))
dmin<-min(sub_dat$ordDat)
dmax<-max(sub_dat$ordDat)
ndays<-dmin:dmax
ndays<-standNew(sub_dat$ordDat,ndays)
nD<-length(ndays)
nSp<-length(sp)
b1<-standNew(sub_dat$PrevFall_tmn,sub_ndat$PrevFall_tmn)
b2<-standNew(sub_dat$PrevFall_ppt,sub_ndat$PrevFall_ppt)
b3<-standNew(sub_dat$Win_tmn,sub_ndat$Win_tmn)
b4<-standNew(sub_dat$pck,sub_ndat$pck)
b7<-standNew(sub_dat$Year,sub_ndat$Year)
nbb<-as.matrix(cbind(rep(b1,nD*nSp),rep(b2,nD*nSp),rep(b3,nD*nSp),rep(b4,nD*nSp),rep(ndays,nSp),rep(ndays^2,nSp),rep(b7,nD*nSp)))
nbbi<-as.matrix(cbind(nbb,nbb[,4]*nbb[,5],nbb[,4]*nbb[,6]))
newSp<-rep(sp,each=nD)
D<-list(X=bbi,N=dim(bbi)[1],K=9,L=length(sp),y=sub_dat$pa,ll=as.numeric(as.factor(sub_dat$genus_species)),newX=nbbi,newN=dim(nbbi)[1],newll=as.numeric(as.factor(newSp)))

	
spid<-as.numeric(as.factor(sub_dat$genus_species))
nspid<-as.numeric(as.factor(newSp))

pac<-tapply(INDEX=sub_dat$genus_species,X=sub_dat$pa,sum)

mn<-mean(sub_dat$ordDate)
sdd<-sd(sub_dat$ordDate)
od<-(nbbi[,5]*sdd)+mn
odd<-as.Date(sprintf("%05d", 22000+od), format = "%y%j")


## plots of posteriors
ppA<-extract(fit,par="pp")[[1]]
ppB<-extract(fit,par="beta")[[1]]
nppA<-extract(fit,par="ppNew")[[1]]
med<-apply(ppB,c(2,3),median)
lb<-apply(ppB,c(2,3),quantile,probs=0.05)
ub<-apply(ppB,c(2,3),quantile,probs=0.95)
brks<-quantile(bb[,4],probs=seq(0,1,1/18))
cs<-rev(heat.colors(n=18))
css<-rep(cs[1],length(spid))
for(i in 2:18){
        css[bb[,4] > brks[i]]<-cs[i]
}

pdf(paste("bayHglm_2022_fit_",sites[j],".pdf",sep=""),width=9,height=12)
par(mfrow=c(4,3))
par(mar=c(4.5,5.5,2.5,1.5))
med_PA<-apply(ppA,2,median)
for(i in 1:length(sp)){
        a<-which(spid==i)
       	plot(bbi[a,5],med_PA[a],pch=19,xlab="Day",ylab="Prob. present",col=alpha(css[a],0.9),cex.lab=1.4)
	        title(main=sp[i],cex.main=1.4)
}
dev.off()
	
pdf(paste("bayHglm_2022_forecast_",sites[j],".pdf",sep=""),width=9,height=12)
par(mfrow=c(4,3))
par(mar=c(4.5,5.5,2.5,1.5))
med_PA<-apply(nppA,2,median)
lb_PA<-apply(nppA,2,quantile,probs=0.05)
ub_PA<-apply(nppA,2,quantile,probs=0.95)
for(i in 1:length(sp)){
        a<-which(nspid==i)
       	plot(od[a],med_PA[a],pch=19,xlab="Day",ylab="Prob. present",col=alpha("darkgray",0.9),ylim=c(0,1),cex.lab=1.4,type='n')
	polygon(c(od[a],rev(od[a])),c(lb_PA[a],rev(ub_PA[a])),col=alpha("darkgray",.4),border=NA)
	lines(od[a],med_PA[a],col=alpha("darkgray",0.9),lwd=1.8)
        title(main=sp[i],cex.main=1.4)
}
dev.off()

null_PA<-matrix(NA,nrow=length(unique(odd)),ncol=8000)
for(i in 1:8000){ ## 8000 HMC samples
	sam_PA<-rbinom(n=length(med_PA),size=1,nppA[i,])
	null_PA[,i]<-tapply(X=sam_PA,INDEX=odd,sum)
}
lb<-apply(null_PA,1,quantile,probs=0.05)
ub<-apply(null_PA,1,quantile,probs=0.95)

pdf(paste("bayHglm_2022_forecast_spec_",sites[j],".pdf",sep=""),width=4.5,height=4.5)
par(mar=c(4.5,4.5,.5,.5))
plot(sort(unique(odd)),tapply(X=med_PA,INDEX=odd,sum),type='n',xlab="Date",ylab="Number of species",col=alpha("darkgray",0.9),cex.lab=1.2,ylim=c(0,1.02*max(ub)))
polygon(c(sort(unique(odd)),rev(sort(unique(odd)))),c(lb,rev(ub)),col=alpha("darkgray",.4),border=NA)
lines(sort(unique(odd)),tapply(X=med_PA,INDEX=odd,sum),lwd=1.8,col=alpha("darkgray",0.9))
dev.off()

ppM<-matrix(NA,nrow=nSp,ncol=nD)
for(i in 1:nSp){
	a<-which(nspid==i)
        ppM[i,]<-med_PA[a][order(od[a])]
}
colTab<-rev(heat.colors(10))
pdf(paste("bayHglm_2022_forecast_all_",sites[j],".pdf",sep=""),width=6,height=18)
par(mar=c(4.5,9,.5,.5))
image.plot(t(ppM[rev(order(apply(ppM,1,which.max))),]),xlab="Date",ylab="",axes=FALSE,cex.lab=1.2,col=colTab)
axis(2,at=(1:nSp-1)/(nSp-1),sp[rev(order(apply(ppM,1,which.max)))],las=2,cex.axis=.7)
sdy<-seq(1,nD,7)
axis(1,at=sdy/(nD-1),sort(unique(format(odd,format="%m-%d")))[sdy],las=2,cex.axis=.7)
box()
dev.off()

med<-apply(ppB,c(2,3),median)
lb<-apply(ppB,c(2,3),quantile,probs=0.05)
ub<-apply(ppB,c(2,3),quantile,probs=0.95)

pdf(paste("effects_2022_",sites[j],".pdf",sep=""),width=9,height=9)
par(mar=c(4.5,5.5,2.5,1.5))
covnms<-c("Fall_tmn","Fall_ppt","Win_tmn","Pack","Day","Day_2","Year","Pack-Day","Pack-Day_2")
par(mfrow=c(3,3))
for(i in 1:9){
  	plot(med[,i],xlab="Species",ylab="Coefficient",pch=19,col=alpha("gray20",.8),
        ylim=c(min(lb[,i]),max(ub[,i])),cex.lab=1.4)
	abline(h=0,lty=2)
	segments(1:nSp,lb[,i],1:nSp,ub[,i])
	title(covnms[i],cex.main=1.4)
}
dev.off()

pdf(paste("hier_effects_2022_",sites[j],".pdf",sep=""),width=5,height=6)
par(mar=c(4.5,5.5,2.5,1.5))
plot(fit,pars="mu",cex.axis=1.2)
plot(fit,pars="sigma",cex.axis=1.2)
dev.off()

