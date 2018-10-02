
RSCpoint<-function(e,pt,A=NULL,n=NULL,verbose=FALSE){
  
  ########## require - RSCform, RSCcond
  
  ########## e - list of marginal vectors
  ########## pt - initial joint vector estimate
  ########## A - matrix transforming pt into induced marginals
  ########## n - values to restrict marginals to
  ########## verbose - keep track of execution process
  
  ########## Step I - automatic formulations when bounds, A, n, left null
  
  if(class(e)=='list' && is.null(A)){                                   ## create A and n from e
    form<-RSCform(e)
    A<-form$A
    n<-form$n  }
  
  ########### Step II - Raking Procedure
  
  cur.g<-pt              ## initial gh estimate
  t<-1                   ## loop iteration number
  ding<-FALSE            ## indication of whether Ap=n has been met
  
  while(ding==FALSE){
    u.cell<-(t-1)%%length(n)+1                                  ## which A-condition is being considered
    v.cells<-A[u.cell,]==1                                      ## which joint vector elements are being considered
    
    adjust<-n[u.cell]/(A[u.cell,]%*%cur.g)                      ## actual adjustment
      if(A[u.cell,]%*%cur.g==0) adjust=1
    
    new.g<-cur.g*(v.cells*adjust+!v.cells)                      ## make adjustment
    new.error<-sum((n-(A%*%new.g))^2)                           ## have overall restrictions (Ap=n) been met
    error.tolerance<-exp(-20+0.01*t/nrow(A))                    ## how much error in Ap=n will I tolerate
    
    if(abs(new.error)<error.tolerance) {ding<-TRUE              ## if Ap=n is met, stop loop
      }else {                                                   ## otherwise
        cur.g<-new.g                                            ## set guess and restart loop
        if(verbose==TRUE) print(t)
        t<-t+1
    }                          
  }
  
  return(list(estimate=new.g,iterations=t,error=new.error))     ## result obtained once loop ends
}
