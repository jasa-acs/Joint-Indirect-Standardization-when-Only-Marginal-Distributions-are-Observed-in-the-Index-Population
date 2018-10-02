###################################################################################################
########## THIS FILE CONTAINS ALL CODE USED TO PRODUCE OUTCOMES IN EXAMPLE DATA SECTION ###########
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

#################################################

cond<-full.out[[3]][,-1]
IR.cutoff<-floor(ncol(jsum)/4)           ## We set 1/4 of the hospitals as "index hospitals" ...
indices<-jsum[,c(1:IR.cutoff)]
ref<-jsum[,c((IR.cutoff+1):ncol(jsum))]  ## ... and the rest as "references"

rakeref<-function(k){                    ## This function rakes an index (real or fake) with the remaining hospitals
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
 

###########################################################
##########PLAYING WITH CONCENTRATION PARAMETERS############
###########################################################

## In a scenario where someone's actually using the software,
##  the "concentration parameter" is pre-computed and associated 
##  with the reference dataset itself

## The code in this section is where the "pre-computing" happens,
##  but running it is the domain of the software developper,
##  not the user



## We now compute the concentration parameter for the 

py<-full.out[[4]]
nlist.ref<-full.out[[5]][c((IR.cutoff+1):(length(full.out[[5]])))]
true.values.ref<-full.out[[2]][c((IR.cutoff+1):(length(full.out[[2]])))]
raked.refs<-lapply(c((IR.cutoff+1):ncol(jsum)),function(x) all.rakerefs[[x]][,c((IR.cutoff+1):(ncol(jsum)-1))])

### Estimated Run-time of the following bit is 250 secs,
###   so the code itself has been commented out and
###   the results have been saved in "RSCshapeObject"
### Feel free to run it yourself
### NOTE - RSCconcentrationBase uses random number generation 
###        to assess concentration parameters. Your results may vary slightly. 

# shape<-RSCconcentration(0.95,py,ref,nlist.ref,true.values.ref,raked.refs,Dexpect=NULL,verbose=TRUE)                   
# save(shape,file="forPublic/RSCshapeObject.Rdata")   

load(file="forPublic/RSCshapeObject.Rdata")








######################################################################################################
######## With the concentration parameter fixed, we apply algorithm to the "index hopsitals" #########
######################################################################################################

### The code in this section mostly reflects the code for RSCconcentrationBase
### Such complex code would not be necessary for users of software produced by this algorithm
### as they would only have to apply to a single index hospital

true.values.ind<-full.out[[2]][1:IR.cutoff]
nlist.ind<-full.out[[5]][1:IR.cutoff]
raked.indices<-lapply(c(1:IR.cutoff),function(x) all.rakerefs[[x]][,c((IR.cutoff):(ncol(jsum)-1))])


#### remainig computation times are trivial

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
          coverage=coverage,emean=all.bds[2,],tobs=Op,tdenom=Ep)

save(test.outputs,file="forPublic/RSCexampleObject.Rdata")      
load(file="forPublic/RSCexampleObject.Rdata")
    
    




##########################################
##### NUMBERS QUOTED IN THE PAPER ########
##########################################

## Full sample size
  sum(full.out[[5]])

## Range of sample sizes within each hospital
  min(full.out[[5]])
  max(full.out[[5]])
  mean(full.out[[5]])
  sd(full.out[[5]])

## Prevalence of high dose per hospital
  summary(full.out[[2]])
  
## Range of prevalence of outcome within each joint predictor level
  min(py)
  max(py)

## Shape
  shape

## Our widths vs. full knowledge widths
  ours<-test.outputs[2,]-test.outputs[1,]
  full<-2*qnorm(0.975)*test.outputs[4,]
  mean(ours/full)
  max(ours/full)

## Shared widths vs. interval widths
  shared<-sapply(1:length(ours),function(i){
    max(0,
        min(test.outputs[3,i]+qnorm(0.975)*test.outputs[4,i],test.outputs[2,i])-
         max(test.outputs[3,i]-qnorm(0.975)*test.outputs[4,i],test.outputs[1,i]))
  })
  
  mean(shared/ours)
  mean(shared/full)
  
  

##############################
###### FIGURE ONE ############
##############################

xlim<-c(0,2.5)

draw.predict.II<-function(test.outputs,title="No Title"){
  
  coverage.rate<-mean(test.outputs[5,])
  
  plot(seq(xlim[1],xlim[2],length=100),seq(xlim[1],xlim[2],length=100),type="l",
       xlab="Observed Standardized Incidence Ratio",
       ylab="Standardized Incidence Ratio Uncertainty Interval",
       main="")
  
  for(j in 1:ncol(test.outputs)){
    
    rect(xleft=test.outputs[3,j]-0.01,xright=test.outputs[3,j]+0.01,ytop=test.outputs[2,j],ybottom=test.outputs[1,j],
         col=NA,border="black")
    
    rect(xleft=test.outputs[3,j]-0.01,xright=test.outputs[3,j]+0.01,
         ytop=test.outputs[3,j]-qnorm(0.975)*test.outputs[4,j],ybottom=test.outputs[3,j]+qnorm(0.975)*test.outputs[4,j],
         col=rgb(0.25,0.25,0.25,alpha=0.5),border=NA)
  }
  
  legend("topleft",legend=c("Using Full Data", "Using Marginal Data"),pch=c(15,0),col=c(rgb(0.25,0.25,0.25,alpha=0.5),"black"),bty="n")
  abline(v=1,lty=2)
  abline(h=1,lty=2)
  
}



png(file="sumplotSIR.png")
draw.predict.II(test.outputs,"High DLP")
dev.off()


