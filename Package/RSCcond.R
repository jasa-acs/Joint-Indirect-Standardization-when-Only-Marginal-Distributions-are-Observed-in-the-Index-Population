
RSCcond<-function(ncatv,npred=NULL){   ######MATRIX OF A B C D E ... COMBINATIONS
  
        ######## ncatv - total number of categories among each predictor
        ######## npred - number of predictors
  
  if(is.null(npred)) npred=length(ncatv)
  
  total.length<-prod(ncatv)  ####### total number of rows
  
  get.vector<-function(i){
    combo.per.cat<-prod(ncatv[(i+1):npred])   ##### number of combinations per category
    return(as.factor(rep(rep((1:ncatv[i])-1,each=combo.per.cat),length=total.length)))}
  
  return(sapply(1:npred,get.vector))

}
