#' This is data to be included in MDLM
#'
#' @author Jinseob Kim \email{jinseob2kim@gmail.com}
#' @references \url{https://github.com/jinseob2kim/MDLM}
"exMDLM"



#' MDLM of fixed vector space
#'
#'
#'
#' @param ybar Dependent variable.
#' @param xbar Independent variable.
#' @param data Data.
#' @param cor.mat Dependency matrix of vector space
#' @param ini.beta Initial beta.
#'
#' @return beta, standard error, z value, p value, rmse
#' @examples
#' # ADD_EXAMPLES_HERE
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
  return(list(summary=summary.optim2(res.optim,f,nrow(data)), rmse=sqrt(1/nrow(data)*res.optim$value)))
}




#exMDLM=read.csv("/home/secondmath/Dropbox/Survival_new/exdata2.csv")
#save(exMDLM, file="data/exMDLM.rda")
#dd=data(exMDLM)
#mdlm1(ybar="y",xbar=c("x1","x2"),data=exMDLM)
#b.start=c(0,0,0,0)
#r.start=0.1

#f3=function(x){CostFun(bs=x[1:length(b.start)],  "y",xbar=c("x1","x2") ,a,LowerToMatrix(x[-(1:length(b.start))]))}
#mdlm3 = constrOptim(c(b.start,r.start),f3, NULL,ui=ui,ci=ci)
#rmse[3]=sqrt(1/nsample * mdlm3$value)
#summary.optim2(mdlm3,f3,nrow(a))

