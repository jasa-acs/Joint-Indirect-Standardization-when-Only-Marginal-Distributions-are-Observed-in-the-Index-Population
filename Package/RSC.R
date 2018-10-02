
RSC<-function(n, op, py,raked.index=NULL,e=NULL,ref=NULL,s,Dexpect=NULL){
  
  set.seed(1234)
  
  ########## specify EITHER "raked.index"
  ##########   OR specify e and ref (required RSCpoint)
  
  ########## n - sample size of the index
  ########## op - observed prevalence of outcome in index hospital
  ########## py - probability of outcome in each level of joint predictor
  
  ########## raked.index - reference hospitals raked to match marginals of the index hospital
  
  ########## e - marginals of the index hospital
  ########## ref - true probability vector of predictor in each reference hospital
  
  ########## s - concentration parameter of the reference dataset
  ########## Dexpect - mechanism for changing the expected value of the Dirichlet random variable
  ##########           for purposes of this paper, it's set to a vector of equal values for now
  
  if(is.null(Dexpect) & !is.null(ref)) Dexpect<-rep(1/(ncol(ref)),ncol(ref))
  if(is.null(Dexpect) & !is.null(raked.index)) Dexpect<-rep(1/(ncol(raked.index)),ncol(raked.index))
  
  if(is.null(raked.index)){
    raked.index<-lapply(1:ncol(ref),function(i) RSCpoint(e=e,pt=ref[,i]))}   ## rake the index with each reference hospital
  
  
  ### Generate SIRs
  ### Dirichlet random variables are generated differently depending on value of shape
  
  gen.one.val<-function(){
    
    smr.up<-rgamma(n=1,shape=n*op+1/3,rate=1)                        ## generate SIR numerator
    
    w.list<-rgamma(n=length(Dexpect),shape=Dexpect*s)                    ## generate many gamma rv
    w<-w.list/sum(w.list)                                            ## standardize to make Dirichlet
    
    smr.down<-(py%*%(raked.index%*%w))*n  ## compute SIR denominator
    
    return(smr.up/smr.down)
  }
  
  gen.one.val.II<-function(){
    
    smr.up<-rgamma(n=1,shape=n*op+1/3,rate=1)                     
    
    w<-c(rmultinom(n=1,size=1,prob=Dexpect))
    smr.down<-(py%*%(raked.index%*%w))*n
    
    return(smr.up/smr.down)
  }
  
  if(s>0.001) {
    gen.multi.val<-function() replicate(n=1000,expr=gen.one.val())
  }else{
    gen.multi.val<-function() replicate(n=1000,expr=gen.one.val.II())
  }
  
  all.gen<-gen.multi.val()
  
  ######################
  ######################
  
  out<-quantile(all.gen,c(0.025,0.5,0.975),na.rm=T)
  return(out)
}
