
RSCconcentrationBase<-function(s,py,ref,nlist.ref,true.values.ref,raked.refs,Dexpect=NULL,verbose=FALSE){
                                          
  set.seed(1234)
  
  ########## This function is meant to be used only by software developers
  ########## to assess performance of a candidate concentration parameter
  
  ########## s - current estimate for the concentration parameter
  ########## py - probability of outcome in each level of joint predictor
  ########## ref - true probability vector of predictor in each reference hospital
  ########## nlist.ref - sample size of each reference hospital
  ########## true.values.ref - true probability of outcome in each reference hospital
  ########## raked.refs - raked probability vectors of the predictor in each reference hospital
  ##########              with each other reference hospital, while pretending the former
  ##########              is an index hospital with only marginals available
  ##########              compute this with repeated "RSCpoint"
  ########## Dexpect - mechanism for changing the expected value of the Dirichlet random variable
  ##########           for purposes of this paper, it's set to a vector of equal values for now
  ########## verbose - if true, prints out coverage rate immediately upon completion
  
  if(is.null(Dexpect)) Dexpect<-rep(1/(ncol(ref)-1),ncol(ref)-1)
  
  ########## The "concentration parameter" is to be used in a Dirichlet random variable
  ########## the primary established means of generating a Dirichlet random variable
  ########## is to generate a bunch of gamma random variables, then standardize to add up to one
  
  ########## The following function randomly generates an SIR
  
  gen.one.val<-function(ref.num){
    
    smr.up<-rgamma(n=1,shape=(true.values.ref*nlist.ref)[ref.num]+1/3,
                   rate=1)                                               ## generate SIR numerator
    
    w.list<-rgamma(n=length(Dexpect),shape=Dexpect*s)                    ## generate many gamma rv
    w<-w.list/sum(w.list)                                            ## standardize to make Dirichlet
    
    smr.down<-(py%*%(raked.refs[[ref.num]])%*%w)*nlist.ref[ref.num]  ## compute SIR denominator
    
    return(smr.up/smr.down)
  }
  
  ########## Existing software for generating gamma random variables break and end up
  ########## difficult to standardize when shape parameters are too small
  ########## But we know that, theoretically, Dirichlet random variable converge
  ########## to multinomial when concentration small
  
  gen.one.val.II<-function(ref.num){
    
    smr.up<-rgamma(n=1,shape=(true.values.ref*nlist.ref)[ref.num]+1/3,
                   rate=1)
    
    w<-c(rmultinom(n=1,size=1,prob=Dexpect))
    smr.down<-(py%*%(raked.refs[[ref.num]])%*%w)*nlist.ref[ref.num]
          
    return(smr.up/smr.down)
  }
  
  ######### Generate some random SIRs, which is a slightly different process
  ######### depending on how small the concentration parameter being tested is
  
  if(s>0.001) {
    gen.multi.val<-function(ref.num) replicate(n=1000,expr=gen.one.val(ref.num))
  }else{
    gen.multi.val<-function(ref.num) replicate(n=1000,expr=gen.one.val.II(ref.num))
  }
  
  all.gen<-sapply(1:ncol(ref),gen.multi.val)
  all.bds<-apply(all.gen,2,function(x) quantile(x[x<Inf],probs=c(0.025,0.975),na.rm=T))   
                  ####### Prediction interval of SIRs given this concentration parameter
  
  Ep<-c(py%*%ref)
  all.bds.gamma<-apply(all.bds,1,function(x) x*(nlist.ref*Ep))
  
  coverage<-pgamma(q=all.bds.gamma[,2],shape=true.values.ref*nlist.ref+1/3,rate=1)-
    pgamma(q=all.bds.gamma[,1],shape=true.values.ref*nlist.ref+1/3,rate=1)
  coverage.rate<-mean(coverage,na.rm=T)
                  ####### Compare distribution of real predictor distribution
                  ####### with results from pretending each reference only have marginals
  
  if(verbose==TRUE) print(paste("Trying concentration",s,"coverage",coverage.rate))
  
  return(coverage.rate)
  
}



