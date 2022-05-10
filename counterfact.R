library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(fields)
library(scales)

## get workspaces
ff<-c(list.files(pattern=glob2rx("forecast2022_*rdat")))

yrs<-c(2002,2011,2015)

for(j in 1:length(ff)){
	load(ff[j])
	
	ppBeta<-extract(fit,par="beta")[[1]]
	ppAlpha<-extract(fit,par="alpha")[[1]]
	nsams<-dim(ppAlpha)[1]
	pest<-vector("list",length(yrs))
	for(y in 1:length(yrs)){
		bs<-unique(bbi[which(sub_dat$Year==yrs[y]),c(1:4)]) ## covars from yrs[y]
		b7<-standNew(sub_dat$Year,sub_ndat$Year)## current year
		## i.e., assumes past climate but this year
		pest[[y]]<-vector("list",nSp)
		for(i in 1:nSp){
			pest[[y]][[i]]<-matrix(NA,nrow=length(ndays),ncol=nsams)
			for(k in 1:nsams){
				est<-ppAlpha[k,i] + bs[1]*ppBeta[k,i,1] + bs[2]*ppBeta[k,i,2] + 
					bs[3]*ppBeta[k,i,3] + bs[4]*ppBeta[k,i,4] +
					ndays*ppBeta[k,i,5] + ndays^2*ppBeta[k,i,6] + b7*ppBeta[k,i,7] 
					bs[4]*ndays*ppBeta[k,i,8] + bs[4]*ndays^2*ppBeta[k,i,9] 
				pest[[y]][[i]][,k]<-1/(1+exp(-1*est))
			}
			

			#plot(med,ylim=c(0,1),type='l')
			#lines(lb)
			#lines(ub)
		}
	}
	mn<-mean(sub_dat$ordDate)
	sdd<-sd(sub_dat$ordDate)
	od<-(nbbi[,5]*sdd)+mn
	odd<-as.Date(sprintf("%05d", 22000+od), format = "%y%j")
	pid<-as.numeric(as.factor(sub_dat$genus_species))
	nspid<-as.numeric(as.factor(newSp))


	nppA<-extract(fit,par="ppNew")[[1]]
	med_PA<-apply(nppA,2,median)
	lb_PA<-apply(nppA,2,quantile,probs=0.05)
	ub_PA<-apply(nppA,2,quantile,probs=0.95)
	cs<-c("green","blue","red")
	pdf(paste("bayHglm_2022_cofrcast_",sites[j],".pdf",sep=""),width=9,height=12)
	par(mfrow=c(4,3))
	par(mar=c(4.5,5.5,2.5,1.5))

	for(i in 1:nSp){
	        a<-which(nspid==i)
	        plot(od[a],med_PA[a],pch=19,xlab="Day",ylab="Prob. present",col=alpha("darkgray",0.9),ylim=c(0,1),cex.lab=1.4,type='n')
        	polygon(c(od[a],rev(od[a])),c(lb_PA[a],rev(ub_PA[a])),col=alpha("darkgray",.4),border=NA)
       		lines(od[a],med_PA[a],col=alpha("darkgray",0.9),lwd=1.8)
        	title(main=sp[i],cex.main=1.4)
		for(y in 1:length(yrs)){
			med<-apply(pest[[y]][[i]],1,median)
			lb<-apply(pest[[y]][[i]],1,quantile,probs=0.05)
			ub<-apply(pest[[y]][[i]],1,quantile,probs=0.95)
			polygon(c(od[a],rev(od[a])),c(lb,rev(ub)),col=alpha(cs[y],.2),border=NA)
			lines(od[a],med,col=alpha(cs[y],.5),lwd=1.5)				
		}
	}
	dev.off()
	
	## species number
	null_PA<-matrix(NA,nrow=length(unique(odd)),ncol=8000)
	for(i in 1:8000){ ## 8000 HMC samples
       		sam_PA<-rbinom(n=length(med_PA),size=1,nppA[i,])
        	null_PA[,i]<-tapply(X=sam_PA,INDEX=odd,sum)
	}
	lb<-apply(null_PA,1,quantile,probs=0.05)
	ub<-apply(null_PA,1,quantile,probs=0.95)

	pdf(paste("bayHglm_2022_cofrcast_spec_",sites[j],".pdf",sep=""),width=4.5,height=4.5)
	par(mar=c(4.5,4.5,.5,.5))
	plot(sort(unique(odd)),tapply(X=med_PA,INDEX=odd,sum),type='n',xlab="Date",ylab="Number of species",col=alpha("darkgray",0.9),cex.lab=1.2,ylim=c(0,1.2*max(ub)))
	polygon(c(sort(unique(odd)),rev(sort(unique(odd)))),c(lb,rev(ub)),col=alpha("darkgray",.4),border=NA)
	lines(sort(unique(odd)),tapply(X=med_PA,INDEX=odd,sum),lwd=1.8,col=alpha("darkgray",0.9))
	for(y in 1:length(yrs)){
		sam_PA<-vector("list",nSp)
		med<-matrix(NA,nrow=nSp,ncol=length(ndays))
		for(i in 1:nSp){ 
			NN<-prod(dim(pest[[y]][[i]]))
       			sam_PA[[i]]<-matrix(rbinom(n=NN,size=1,pest[[y]][[i]]),nrow=dim(pest[[y]][[i]])[1],
				ncol=dim(pest[[y]][[i]])[2])
			med[i,]<-apply(pest[[y]][[i]],1,median)
		}
		null_PA<-Reduce("+",sam_PA)
		lb<-apply(null_PA,1,quantile,probs=0.05)
		ub<-apply(null_PA,1,quantile,probs=0.95)
		med<-apply(med,2,sum)
		polygon(c(sort(unique(odd)),rev(sort(unique(odd)))),c(lb,rev(ub)),col=alpha(cs[y],.4),border=NA)
		lines(sort(unique(odd)),med,lwd=1.8,col=alpha(cs[y],0.9))
	}


	dev.off()
	
	## heat maps
	ppMref<-matrix(NA,nrow=nSp,ncol=nD)
	for(i in 1:nSp){
        	a<-which(nspid==i)
        	ppMref[i,]<-med_PA[a][order(od[a])]
	}

	for(y in 1:length(yrs)){

		ppM<-matrix(NA,nrow=nSp,ncol=nD)
		for(i in 1:nSp){
		        ppM[i,]<-apply(pest[[y]][[i]],1,median)
		}
		colTab<-rev(heat.colors(10))
		pdf(paste("bayHglm_2022_forecast_all_cfyr_",yrs[y],"_",sites[j],".pdf",sep=""),width=6,height=4.5)
		par(mar=c(4.5,4.5,.5,.5))
		image.plot(t(ppM[rev(order(apply(ppMref,1,which.max))),]),xlab="Date",ylab="",axes=FALSE,cex.lab=1.2,col=colTab)
		axis(2,at=((1:nSp)-1)/nSp,sp[rev(order(apply(ppMref,1,which.max)))],las=2,cex.axis=.2)
		sdy<-seq(1,nD,7)
		axis(1,at=sdy/nD,sort(unique(format(odd,format="%m-%d")))[sdy],las=2,cex.axis=.7)
		box()
		dev.off()
	}

}
#dim(ppBeta)
#[1] 8000   69    9
#dim(ppAlpha)
#[1] 8000   69

