ddWt<-function(wt,biomass,bmsy,alpha=-0.2,min=0.5,max=1.5){
  adj=(biomass%/%bmsy)^alpha
  
  wt%*%qmin(qmax(adj,min),max)}


ddMat<-function(wt,x) 
  1/(1+exp(-x[1]%*%(wt%-%x[2])))

ddM<-function(wt, F, gamma=2.0753+0.9914*F, lambda=-0.0276-0.1088*F, refWt=310, refM=0.10, M1=0.05){
  
  ref=refM/(gamma%*%exp(lambda*(refWt/0.0003)^(1/3)))
  
  exp(((wt/0.0001)^(1/3))%*%lambda)%*%gamma%*%ref+M1}

ddFn<-function(year,x,par){
  lsr=landings.wt(x)[,year]/stock.wt(x)[,year]
  dsr=discards.wt(x)[,year]/stock.wt(x)[,year]
  dsr[!is.finite(dsr)]=0
  
  stock.wt(   x)[,year]=ddWt(stock.wt(x)[,year],biomass(x)[,ac(an(year)-1)],par["bmsy"],alpha=-0.25)
  landings.wt(x)[,year]=stock.wt(x)[,year]*lsr
  discards.wt(x)[,year]=stock.wt(x)[,year]*dsr
  
  mat(x)[,year]=ddMat(stock.wt(x)[,year],par[c("k","w50")])
  m(  x)[,year]=par["m1"]%*%(stock.wt(x)[,year]%^%par["m2"])
  
  x}

ddFn2<-function(x,par){
  lsr=landings.wt(x)/stock.wt(x)
  dsr=discards.wt(x)/stock.wt(x)
  dsr[!is.finite(dsr)]=0
  
  stock.wt(   x)=ddWt(stock.wt(x),biomass(x),par["bmsy"],alpha=-0.25)
  landings.wt(x)=stock.wt(x)*lsr
  discards.wt(x)=stock.wt(x)*dsr
  
  mat(x)=ddMat(stock.wt(x),par[c("k","w50")])
  m(  x)=par["m1"]%*%(stock.wt(x)%^%par["m2"])
  
  x}
