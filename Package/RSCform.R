

RSCform<-function(e){
  e<-lapply(e,function(x) x/sum(x)) 
  cond<-as.data.frame(RSCcond(ncatv=sapply(e,length)))
  Araw<-Reduce(cbind,lapply(1:ncol(cond),
                            function(i) sapply(levels(cond[,i]),function(x) as.numeric(cond[,i]==x))))
  A<-t(sapply(1:ncol(Araw),function(i) Araw[,i]))
  n<-unlist(e)
  
  return(list(A=A,n=n))
}

