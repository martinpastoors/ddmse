ddWt<-function(wt,biomass,bref,delta=-0.2,min=0.5,max=1.5){
  adj=(biomass%/%bref)^delta
  
  wt%*%qmin(qmax(adj,min),max)}

ddMat<-function(wt,x) 
  1/(1+exp(-x[1]%*%(wt%-%x[2])))

ddM<-function(wt,par){
    par["m1"]%*%(wt%^%par["m2"])}
  
ddFn<-function(year,x,par,massFlag=TRUE,matFlag=TRUE,mFlag=TRUE){
  
  if (massFlag){
    lsr=landings.wt(x)[,year]/stock.wt(x)[,year]
    dsr=discards.wt(x)[,year]/stock.wt(x)[,year]
    dsr[!is.finite(dsr)]=0
    
    stock.wt(   x)[,year]=ddWt(stock.wt(x)[,year],biomass(x)[,ac(an(year)-1)],par["bref"],par["delta"])
    landings.wt(x)[,year]=stock.wt(x)[,year]*lsr
    discards.wt(x)[,year]=stock.wt(x)[,year]*dsr}
    
  if (matFlag) 
    mat(x)[,year]=ddMat(stock.wt(x)[,year],par[c("matk","w50")])

  if (mFlag)   
    m(  x)[,year]=par["m1"]%*%(stock.wt(x)[,year]%^%par["m2"])
  
  x}
