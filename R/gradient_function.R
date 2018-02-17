

#' Generate initial beta
#'
#' Get initial beta from conventional linear model
#'
#' @param ybar Dependent variable.
#' @param xbar Independent variable.
#' @param data Data.
#'
#' @return Initial beta
#' @examples
#' # ADD_EXAMPLES_HERE
IniBetaMake=function(ybar,xbar,data){
  coef.list= coef(lm(as.formula(paste(ybar,"~",paste(xbar,collapse="+"))),data=data))
  p=length(xbar)
  return(list(beta.ini=coef.list[-1],beta0.ini=rep(coef.list[1]/p,p)))
}



#' Square sum for CostFun
#'
#' Square sum for CostFun
#'
#' @param bs Intial beta
#' @param ybar Dependent variable.
#' @param xbar Independent variable.
#' @param data Data.
#' @param rmat Correlation matrix.
#'
#' @return Squared sum value
#' @examples
#'
SqrtFun=function(bs=b.start, ybar="y",xbar,data,rmat=r){
  y=data[,ybar]
  x=data[,xbar]
  p=ncol(x)
  n=nrow(data)
  #beta=bs[1:p]
  beta=matrix(rep(bs[1:p],n),ncol=p,byrow=T)
  #beta0=bs[-(1:p)]
  beta0=matrix(rep(bs[-(1:p)],n),ncol=p,byrow=T)
  lx= beta*x+beta0

  jsum=function(i){
    ri = matrix(rep(rmat[i,],n),ncol=p,byrow=T)
    return(rowSums(lx[,i]*lx*ri))
  }
  return(rowSums(sapply(1:p,jsum)))
}



#' Cost function for MDLM
#'
#' Least square estimation for MDLM
#'
#' @param bs Intial beta
#' @param ybar Dependent variable.
#' @param xbar Independent variable.
#' @param data Data.
#' @param rmat Correlation matrix.
#'
#' @return Sum of squared error for MDLM
#' @examples
#' # ADD_EXAMPLES_HERE
CostFun= function(bs=b.start, ybar="y",xbar,data,rmat=r){
  y=data[,ybar]
  x=data[,xbar]
  p=ncol(x)
  n=nrow(data)
  beta=matrix(rep(bs[1:p],n),ncol=p,byrow=T)
  beta0=matrix(rep(bs[-(1:p)],n),ncol=p,byrow=T)
  return(sum( (y-sqrt(SqrtFun(bs, ybar,xbar,data,rmat)))^2))
}



#' Convert lower triangle to Correlation matrix
#'
#' Convert lower triangle to Correlation matrix
#'
#' @param lower Lower triangle. The length is must be x(x-1)/2 form
#'
#' @return Matrix
#' @examples
#' # ADD_EXAMPLES_HERE
LowerToMatrix=function(lower=c(4,7,8)){
  n=(1+sqrt(1+8*length(lower)))/2
  dmat=diag(n)
  dmat[lower.tri(dmat)] = lower
  dmat[upper.tri(dmat)] = lower
  return(dmat)
}
