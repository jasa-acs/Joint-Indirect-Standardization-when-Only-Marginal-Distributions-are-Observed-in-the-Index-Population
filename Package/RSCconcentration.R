
RSCconcentration<-function(desired.coverage,py,ref,nlist.ref,true.values.ref,raked.refs,Dexpect=NULL,verbose=FALSE,id=FALSE){
  
  ######## require RSCconcentrationBase
  ######## Gets you to concentration parameter that produces a desired coverage rate
  
  ######## desired.coverage - desired coverage rate
  ######## id - whether to give information on expected coverage
  ######## all other arguments passed to RSCconcentrationBase
  
  RSCconcentrationT<-function(s) RSCconcentrationBase(s,py=py,ref=ref,
                nlist.ref=nlist.ref,true.values.ref=true.values.ref,raked.refs=raked.refs,Dexpect=Dexpect,
                verbose=verbose)-desired.coverage
  
  maxval.shape<-RSCconcentrationT(10e-4)
  minval.shape<-RSCconcentrationT(10e4)
  
  if(maxval.shape<0) { shape<-10e-4 ; shape.id<-maxval.shape+desired.coverage
   } else if(minval.shape>0) { shape<-10e4  ; shape.id<-minval.shape+desired.coverage
   } else {
     shape<-uniroot(f=RSCconcentrationT, interval=c(10e-4,10e4), tol=0.001)$root
     shape.id<-0.95}
  
  if(id==FALSE) return(shape)
  if(id==TRUE) return(c(shape,shape.id))
}



