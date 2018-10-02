###################################################################################################
########## THIS FILE CONTAINS ALL CODE USED TO PRODUCE OUTCOMES IN SIMULATION SECTION #############
###################################################################################################

###################################################################################################
########## File meant to be used for public viewing, and uses datasets not      ###################
########## used in the accompanied manuscript and thus will not produce the      ##################
########## same results. The dataset used in the accompnied manuscript            #################
########## is not available, to protect patient confidentiality                 ###################
###################################################################################################

###################################################################################################
############################## Code written in R version 3.4.3 ####################################
###################################################################################################


load("forPublic/RSCdataObject.Rdata")
sapply( list.files("Package/", full.names=TRUE), source )

## The files in the "Package" folder, along with a properly formatted reference hospital,
## are sufficient to run the software

## Everything in this file exists to assess the performance of the software, not just to
## run the software itself

## Structure of full.out:
## 1 - Joint probability vectors (rows) for the predictors in 151 hospitals (columns)
## 2 - Prevalence of outcome in 151 hospitals
## 3 - Table referring the names of levels in joint probability vector to the marginal levels they represent
##     Used to extract marginal distributions from [[1]]
## 4 - Expected prevalence of outcome in each joint level in predictors
## 5 - Sample size of each hospital


## This paper assumes the predictor has non-zero prevalence in all levels

jsum<-full.out[[1]]
jsum<-jsum+1e-10
jsum<-apply(jsum,2,function(x) x/sum(x))

## Like with the example data, we separate the hospitals into some "indices" and "references"

py<-full.out[[4]]

cond<-full.out[[3]][,-1]
IR.cutoff<-floor(ncol(jsum)/4)
indices<-jsum[,c(1:IR.cutoff)]
ref<-jsum[,c((IR.cutoff+1):ncol(jsum))]

## But this time, instead of viewing the outcome incidence as observed, we view them as the "true values"
## used to simulate "observed hospitals"

true.overs<-full.out[[2]]
true.smrs<-(full.out[[2]][1:IR.cutoff])/c(py%*%indices)



##############################################################
#### We rake everything once again just list Example Data ####
##############################################################

rakeref<-function(k){                    
  print(k)
  index<-cbind(cond,jsum[,k])
  e<-lapply(1:ncol(cond),function(i) tapply(jsum[,k],index[,i],sum))         ## extract marginal distributions
  
  rakes<-lapply(c(1:ncol(jsum))[-k],function(i) RSCpoint(e=e,pt=jsum[,i]))   ## rake each fake index with the other hospitals
  rjsum<-sapply(rakes,function(x) x$estimate)                                ## put everything in one matrix
  
  return(rjsum)
}


### Estimated Run-time of the following bit is 250 secs,
###   so the code itself has been commented out and
###   the results have been saved in "RSCrakedObject"
### Feel free to run it yourself

#   all.rakerefs<-lapply(c(1:ncol(jsum)),rakeref)                   
#   save(all.rakerefs,file="forPublic/RSCrakedObject.Rdata")  


load(file="forPublic/RSCrakedObject.Rdata")

nlist.ref<-full.out[[5]][c((IR.cutoff+1):(length(full.out[[5]])))]
raked.refs<-lapply(c((IR.cutoff+1):ncol(jsum)),function(x) all.rakerefs[[x]][,c((IR.cutoff+1):(ncol(jsum)-1))])

nlist.ind<-full.out[[5]][1:IR.cutoff]
raked.indices<-lapply(c(1:IR.cutoff),function(x) all.rakerefs[[x]][,c((IR.cutoff):(ncol(jsum)-1))])


######################################
######################################
######################################




get.outputs<-function(tick){

  set.seed(tick)

  ## In every simulated scenario, we keep the original predictor distributions, 
  ## but we simulate the outcome via a Poisson distribution
  
  full.out[[2]]<-rpois(n=length(true.overs),lambda=true.overs*full.out[[5]])/full.out[[5]]
  

  true.values.ref<-full.out[[2]][c((IR.cutoff+1):(length(full.out[[2]])))]
  true.values.n<-(true.values.ref*nlist.ref)
  
  shape.list<-RSCconcentration(0.95,py,ref,nlist.ref=true.values.n,true.values.ref,raked.refs,Dexpect=NULL,verbose=TRUE,id=TRUE)
  
  shape<-shape.list[1]
  good.shape.id<-shape.list[2]
  
  print(paste("shape found",good.shape.id,tick))
  
  #################################################################
  ###### Apply to "index hospitals" with concentration fixed ######
  #################################################################

  true.values.ind<-full.out[[2]][1:IR.cutoff]
  nlist.ind<-full.out[[5]][1:IR.cutoff]
  raked.indices<-lapply(c(1:IR.cutoff),function(x) all.rakerefs[[x]][,c((IR.cutoff):(ncol(jsum)-1))])
  
  all.bds<-sapply(1:length(raked.indices),function(i){
    print(i)
    out<-RSC(n=nlist.ind[i],op=true.values.ind[i],py=py,raked.index=raked.indices[[i]],s=shape)
    return(out)
  })
  
  Ep<-c(py%*%indices)
  Op<-true.values.ind
  real.mean<-Op/Ep        ### SIR if we had full information
  real.sd<-sqrt(Op*nlist.ind)/(Ep*nlist.ind)
  
  all.bds.gamma<-apply(all.bds,1,function(x) x*nlist.ind*Ep)
  
  coverage<-pgamma(q=all.bds.gamma[,3],shape=Op*nlist.ind+1/3,rate=1)-
    pgamma(q=all.bds.gamma[,1],shape=Op*nlist.ind+1/3,rate=1)
  
  test.outputs<-rbind(lower=all.bds[1,],upper=all.bds[3,],tmean=real.mean,tsd=real.sd,
                      coverage=coverage,emean=all.bds[2,],tobs=Op,tdenom=Ep,shape=shape,good.shape.id=good.shape.id)
  
  return(test.outputs)

}



### Estimated Run-time of the following bit is 250-300 secs PER simulation,
###   so the code itself has been commented out and
###   all 1000 results have been saved in "RSCsimObject"
### Feel free to run it yourself
### NOTE - RSCconcentrationBase uses random number generation 
###        to assess concentration parameters. Your results may vary slightly. 


  #bunch.of.outputs<-list()

##  for(tick in c(1:1000)) bunch.of.outputs[[tick]]<-get.outputs(tick)
##  save(bunch.of.outputs,file="forPublic/RSCsimObject.Rdata") 


load("forPublic/RSCsimObject.Rdata")

  




####################################################
############ NON-RANDOM DENOMINATOR ################
####################################################

#### This section uses the overall reference dataset as a single "reference hospital"
#### instead of assessing the variance between hospitals within the reference dataset

bare.ref.sub<-apply(sapply(1:ncol(ref),function(i) ref[,i]*nlist.ref[i]),1,sum)
bare.ref<-bare.ref.sub/sum(bare.ref.sub)            ### There's only one reference here

rakeref.bare<-function(k){             #### FUNCTION TO RAKE ONE HOSPITAL AGAINST FULL DATA
  print(k)
  index<-cbind(cond,jsum[,k])    
  e<-lapply(1:ncol(cond),function(i) tapply(jsum[,k],index[,i],sum)) #### index marginals
  
  rakes<-RSCpoint(e=e,pt=bare.ref)             #### raked refereces
  rjsum<-rakes$estimate                        #### raked reference joints
  
  return(rjsum)
}


all.bare.rakerefs<-sapply(1:IR.cutoff,rakeref.bare)   ##### Rake each index
bare.adj<-c(py%*%all.bare.rakerefs)*nlist.ind         ##### Adjusted SIR denominator of each index


genned.obs<-sapply(1:length(bunch.of.outputs),function(j) bunch.of.outputs[[j]][3,]*c(py%*%indices))
                                                      ##### All observed outcome incidences in the simulated indices


bunch.of.bare.outputs.func<-function(j){              ##### Get my SIRs for each set of simulated indices
  
  print(j)
  
  true.obs<-genned.obs[,j]*nlist.ind
  
  bare.mean<-true.obs/bare.adj
  bare.sd<-sqrt(true.obs)/bare.adj
  
  bare.upper<-bare.mean+qnorm(0.975)*bare.sd
  bare.lower<-bare.mean-qnorm(0.975)*bare.sd
  
  return(rbind(bare.lower,bare.upper,bare.mean))
}



bunch.of.bare.outputs<-lapply(1:length(bunch.of.outputs),bunch.of.bare.outputs.func)




################################################################
############## All values presented in paper ###################
################################################################



### How often is the theoretical coverage rate NOT 95%?

mean(sapply(bunch.of.outputs,function(x) x["good.shape.id",1])==0.95)



#######################################
### Our method vs Fixed Denominator ###
#######################################

is.in<-sapply(c(1:length(bunch.of.outputs)),function(i) true.smrs<bunch.of.outputs[[i]][2,] & 
                true.smrs>bunch.of.outputs[[i]][1,])
coverage.per.hospital<-apply(is.in,1,mean)

is.in.bare<-sapply(c(1:length(bunch.of.bare.outputs)),function(i) true.smrs<bunch.of.bare.outputs[[i]][2,] & 
                     true.smrs>bunch.of.bare.outputs[[i]][1,])
coverage.per.hospital.bare<-apply(is.in.bare,1,mean)


summary(coverage.per.hospital)  ### Our method coverage rates
IQR(coverage.per.hospital)

summary(coverage.per.hospital.bare)  ### Fixed coverage rates
IQR(coverage.per.hospital.bare)

#### Figure 2

png(file="sumplotSIRsim1.png")
plot(true.smrs,coverage.per.hospital*100,xlab="True Standardized Incidence Ratio",ylab="Prediction Rate",log="x",yaxt="n",ylim=c(80,100))
points(true.smrs,coverage.per.hospital.bare*100,pch=4)
axis(2, at=c(80,85,90,95,100), lab=paste(c(80,85,90,95,100),"%"))
abline(h=95,lty=2)
legend("bottomright",pch=c(1,4),legend=c("Our Proposed Method","Fixed Denominator Method"))
dev.off()



#######################################
############ Percent Biases ###########
#######################################

b.plot.y<-unlist(lapply(bunch.of.outputs,function(x) (true.smrs-x[3,])/true.smrs))
b.plot.x<-rep(true.smrs,times=length(bunch.of.outputs))
b.plot.F<-apply(sapply(bunch.of.outputs,function(x) (true.smrs-x[3,])/true.smrs),1,mean)

summary(b.plot.F)   #### within hospital summaries

summary(b.plot.y)   #### overall bias summaries
IQR(b.plot.y)

max(b.plot.x)       #### Highest SIR in dataset


#### Fixed Denominator Method (not shown in paper but we did mention they were similar)

b.plot.y.bare<-unlist(lapply(bunch.of.bare.outputs,function(x) (true.smrs-x[3,])/true.smrs))
b.plot.F.bare<-apply(sapply(bunch.of.bare.outputs,function(x) (true.smrs-x[3,])/true.smrs),1,mean)

summary(b.plot.F.bare)
summary(b.plot.y.bare)


#### Figure 3

png(file="sumplotSIRsim2.png")
plot(b.plot.x,b.plot.y*100,xlab="True Standardized Incidence Ratio",ylab="Percent Bias",
     pch=c(16,4)[as.numeric(unlist(is.in)==FALSE)+1],col=c("grey","black")[as.numeric(unlist(is.in)==FALSE)+1],
     log="x",yaxt="n")
axis(2, at=c(-200,-100,0,100), lab=paste(c(-200,-100,0,100),"%"))
abline(h=0,lty=2)
dev.off()


