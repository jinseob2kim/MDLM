#' This is data to be included in MDLM
#'
#' @author Jinseob Kim \email{jinseob2kim@gmail.com}
#' @references \url{https://github.com/jinseob2kim/MDLM}
"exMDLM"


mdlm1=function(ybar, xbar, data, cor.mat=NULL, ini.beta=NULL){
  data= data[complete.cases(data[,c(ybar,xbar)]),]

  if (is.null(ini.beta)){
    ini.beta=unlist(IniBetaMake(ybar,xbar,data))
  }
  if (is.null(cor.mat)){
    cor.mat=LowerToMatrix(rep(0,length(xbar)*(length(xbar)-1)/2))
  }
  f= function(x){CostFun(bs=x, "y",xbar=c("x1","x2"),data,cor.mat)}
  res.optim= optim(ini.beta,f,NULL,hessian=T)
  return(summary.optim2(res.optim,f,nrow(data)))
}




#exMDLM=read.csv("/home/secondmath/Dropbox/Survival_new/exdata2.csv")
#save(exMDLM, file="data/exMDLM.rda")
#dd=data(exMDLM)
#mdlm1(ybar="y",xbar=c("x1","x2"),data=a)


