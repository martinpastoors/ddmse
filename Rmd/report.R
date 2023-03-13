# ------------------------------------------------------------------------------
# report.R
#
# Generating DDMSE results to be included into Rmd files
#
# Laurie Kell and Martin Pastoors
#
# 09/02/2023 First coding
# 13/02/2023 Implemented all codes from cond, DD and OM
# ------------------------------------------------------------------------------

rm(list=ls())

library(FLCore)    # install.packages("FLCore", repos="http://flr-project.org/R")
library(FLBRP)     # install.packages("FLBRP", repos="http://flr-project.org/R")
library(FLasher)   # install.packages("FLasher", repos="http://flr-project.org/R")
library(FLAssess)  # install.packages("FLAssess", repos="http://flr-project.org/R")
library(FLife)     # install.packages("FLife", repos="http://flr-project.org/R")
library(ggplotFL)
library(kobe)
library(grid)

library(popbio)

library(GGally)
library(ggpubr)
library(ggrepel)   # labelling points on a graph

library(rjson) # install.packages("rjson"); used for dropbox function
library(RJSONIO)

library(plyr)
library(dplyr)
library(reshape)

library(readxl)
  
library(mgcv)
library(mgcViz)

library(tidyverse)

theme_set(theme_bw(16))

source("R/dd.R")
source("R/get_dropbox.r")
source("R/pm.R")
source("R/kobePhaseMar4.R")

dropboxdir<-try(file.path(get_dropbox(), "DDMSE"))

if ("try-error" %in% is(dropboxdir)) 
  dropboxdir="~/Dropbox/DDMSE"

myparams <- data.frame(
  stock    = c("mac"                       , "whb"            ),
  mystkname=c("Northeast Atlantic mackerel", "Blue whiting"   ),
  mylatin = c("Scomber scombrus"           ,"Micromesistius poutassou"),
  m1txt   = c("0.15/(300^-0.288)"          , "0.2/(150^-0.288)"),
  m1      = c( 0.15/(300^-0.288)           ,  0.2/(150^-0.288) ), # 300 g MAC, 150 g WHB 
  m2      = c(-0.288                       , -0.288          ),
  m3      = c(NA                           , NA               ),
  m       = c("lorenzen"                   , "lorenzen"       ),
  blim    = c(2000000                      , 1500000          ),
  btrig   = c(2580000                      , 2250000          ),
  bmsy    = c(3500000                      , 4000000          ), # guesstimates
  alpha   = c(-0.2                         , -0.2             ),
  delta   = c(-0.2                         , -0.2             ),
  matk    = c(0.2                           , 0.1              ), # source: conditioning
  wt1     = c(0.01                         , NA               ), # for resetting weight at age 1
  minage  = c(0                            , 1                ),
  maxage  = c(12                           , 10               ),
  minyear = c(1991                         , 2000             ),
  maxyear = c(2050                         , 2050             ),
  m1scaler= c(8                            , 10               ),
  w50     = c(0.166                        , 0.08             ),
  matk2   = c(28                           , 63               ),
  steepness=c(0.8                          , 0.5               ))

# save(myparams, file=file.path(dropboxdir, "data", "inputs", "myparams.RData"))

mystk     <- "mac";  
# mystk     <- "whb";  
# for (mystk in c("mac","whb")) {
  
  mystkname  <- myparams[myparams$stock==mystk,"mystkname"]
  mylatin    <- myparams[myparams$stock==mystk,"mylatin"]
  maxyear    <- myparams[myparams$stock==mystk,"maxyear"]
  minage     <- myparams[myparams$stock==mystk,"minage"]
  maxage     <- myparams[myparams$stock==mystk,"maxage"]
  nages      <- maxage - minage + 1  
  
  figuresdir <- file.path(dropboxdir, "results", mystk, "figures")
  tablesdir  <- file.path(dropboxdir, "results", mystk, "tables")

  # load(file = file.path(dropboxdir, "results", mystk, paste(mystk,"section0.RData", sep="_")))
  
  
  # Get FLStock and FLBRP object  ----------------------------------------------
  
  load(file.path(dropboxdir,"data","om",paste0(mystk,".RData"))) # only used for FLBRP part; comes from 
  stk   = get(mystk)
  stkR  = get(paste0(mystk,"R"))
  
  if (mystk == "mac") {
    # ts    = as.data.frame(readxl::read_excel(file.path(dropboxdir,"data/inputs/dat.xlsx"),sheet="ts"))
    ts=model.frame(FLQuants(stk,tb=stock,ssb=ssb),drop=T)
    
    stock.wt(   stk)[1] = myparams[myparams$stock==mystk,"wt1"]
    catch.wt(   stk)[1] = myparams[myparams$stock==mystk,"wt1"]
    landings.wt(stk)[1] = myparams[myparams$stock==mystk,"wt1"]
    discards.wt(stk)[1] = myparams[myparams$stock==mystk,"wt1"]
    mat(        stk)[1] = 0
  } else if (mystk == "whb") {
    
    ts=model.frame(FLQuants(stk,tb=stock,ssb=ssb),drop=T)
    ts$tb=ts$tb/1e6
  }
  
  save(.,
       file = file.path(dropboxdir, "results", mystk, paste(mystk,"section0.RData", sep="_")))
  
  # ------------------------------------------------------------------------------
  # 1. Introduction
  # ------------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ------------------------------------------------------------------------------
  # 2. Conditioning
  # ------------------------------------------------------------------------------
  load(file = file.path(dropboxdir, "results", mystk, paste(mystk,"section0.RData", sep="_")))
  
  section <- "02"
  
  wt     = transform(as.data.frame(stock.wt(stk),drop=TRUE),cohort=year-as.numeric(age))
  wt     = merge(wt,ts,by="year")
  wt     = ddply(wt,.(age), transform, rsdl =data/mean(data))
  wt     = ddply(wt,.(age), transform, rsdl_=rsdl/var(data)^0.5)
  wt$age = factor(wt$age,levels=sort(unique(wt$age)))
  wt     = arrange(wt, year, age)
  
  mat    = transform(as.data.frame(mat(stk),drop=TRUE),cohort=year-as.numeric(age))
  mat    = merge(mat,ts)
  mat    = ddply(mat,.(age), transform, rsdl =data/mean(data))
  mat    = ddply(mat,.(age), transform, rsdl_=rsdl/var(data)^0.5) 
  mat$age= factor(mat$age,levels=sort(unique(mat$age)))
  mat    = arrange(mat, year, age)
  
  m      = transform(as.data.frame(m(stk),drop=TRUE),cohort=year-as.numeric(age))  
  m      = merge(m,ts)
  m      = ddply(m,.(age), transform, rsdl =data/mean(data))
  m      = ddply(m,.(age), transform, rsdl_=rsdl/var(data)^0.5) 
  m$age  = factor(m$age,levels=sort(unique(m$age))) 
  m      = arrange(m, year, age)
  
  # plot weight at age ---------------------------------------------------------
  p <-
    ggplot(wt)+  
    geom_line(aes(year,data,group=cohort),linetype=3,colour="grey10")+
    geom_line(aes(year,data,group=age,col=age))+
    geom_point(aes(year,data,col=age))+
    ggrepel::geom_text_repel(data = tail(wt, length(unique(wt$age))), 
                             aes(year, data,label=age, colour=age),
                             hjust = "outward", vjust=0.5, direction = "y", min.segment.length = 100) +
    xlab("Year")+ylab("Mass-at-age")+
    theme(legend.position="none")+
    # guides(colour = guide_legend(nrow = 1)) +
    scale_color_manual("Age",values=rainbow(length(unique(wt$age))))
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "waa.jpg", sep="_")),
       width=10, height=6, units="in", res=300)
  print(p) 
  dev.off()

  # plot maturity at age -------------------------------------------------------
  p <-
    ggplot(mat)+  
    geom_line(aes(year,data,group=cohort),linetype=3,colour="grey10")+
    geom_line(aes(year,data,group=age,col=age))+
    geom_point(aes(year,data,col=age))+
    ggrepel::geom_text_repel(data = tail(mat, length(unique(mat$age))), 
                             aes(year, data,label=age, colour=age),
                             hjust = "outward", vjust=0.5, direction = "y", min.segment.length = 100) +
    xlab("Year")+ylab("Maturity-at-age")+
    theme(legend.position="none")+
    # guides(colour = guide_legend(nrow = 1)) +
    scale_color_manual("Age",values=rainbow(length(unique(mat$age))))
  

  jpeg(filename=file.path(figuresdir, paste(section,mystk, "mat.jpg", sep="_")),
       width=10, height=6, units="in", res=300)
  print(p) 
  dev.off()
  
  
  # plot m at age --------------------------------------------------------------
  p <-
    ggplot(m)+  
    geom_line(aes(year,data,group=cohort),linetype=3,colour="grey10")+
    geom_line(aes(year,data,group=age,col=age))+
    geom_point(aes(year,data,col=age))+
    ggrepel::geom_text_repel(data = tail(m, length(unique(m$age))), 
                             aes(year, data,label=age, colour=age),
                             hjust = "outward", vjust=0.5, direction = "y", min.segment.length = 100) +
    xlab("Year")+ylab("M-at-age")+
    theme(legend.position="none")+
    # guides(colour = guide_legend(nrow = 1)) +
    scale_color_manual("Age",values=rainbow(length(unique(m$age))))
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "m.jpg", sep="_")),
       width=10, height=6, units="in", res=300)
  print(p) 
  dev.off()
  
  # Density dependence in weights at age ---------------------------------------
  
  p1=ggplot(wt)+
    geom_boxplot(aes(age,data))+
    ylab("Mass-at-age")+xlab("Age")
  p2=subset(wt,rsdl<3&rsdl>0) %>% 
    mutate(tb = tb/1000) %>% 
    ggplot()+
    geom_point(aes(tb,rsdl))+
    geom_smooth(aes(tb,rsdl),method="lm")+
    xlab("Total biomass (1000 tonnes)")+ylab("Residual") +
    scale_x_continuous(breaks = scales::pretty_breaks(n=2)) +
    facet_wrap(age~.)
  p <-
    ggarrange(p1,p2 ,label.y=1,  
            widths = c(4,6), heights = c(1,1,1),
            nrow   = 1,          ncol    = 2)+
    theme_bw(16)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())

  jpeg(filename=file.path(figuresdir, paste(section,mystk, "dd_waa.jpg", sep="_")),
       width=10, height=8, units="in", res=300)
  print(p) 
  dev.off()
  
  
  # GAMM model and diagnostics of weight vs Total biomass -----------------------
  
  # drop age 0 if existing
  wt2 = subset(wt,age!=0)
  nn  = length(unique(wt2$age))
  
  # run gammV  
  gmr=mgcViz::gammV(log(data) ~ s(tb/mean(tb), bs="tp") + age, 
                    data=wt2)
  
  # plot gammV
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "gamm.jpg", sep="_")),
       width=10, height=8, units="in", res=300)
  plot(gmr) + labs(x="Total biomass", y="Weight deviation")
  dev.off()  
  
  # summary table
  fileConn <-file(file.path(tablesdir, paste(section,mystk, "gamm summary.txt", sep="_")))
  summary(gmr) %>% capture.output() %>% writeLines(., con=fileConn)
  close(fileConn)
  
  # summary coefficients
  fileConn <-file(file.path(tablesdir, paste(section,mystk, "gamm coefficients.txt", sep="_")))
  as.data.frame(coefficients(gmr)) %>% capture.output() %>% writeLines(., con=fileConn)
  close(fileConn)
  
  # plot(exp(coefficients(gmr)[2:nn]))  # NOT NEEDED?
  
  # connect prediction to wt dataset  
  wt2=cbind(wt2,hat=predict(gmr))
  wt2=merge(wt2,data.frame(age=1:nn,
                           mean=exp(c(0,coefficients(gmr)[2:nn]))))
  
  # Linear   
  p <-
    ggplot(wt2)+
    geom_line(aes(tb,exp(hat)/mean,col=age))+
    geom_point(aes(tb,data/mean,col=age))+
    xlab("Total biomass")+ylab("Mass-Mean") 
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "lm_waa.jpg", sep="_")),
       width=10, height=8, units="in", res=300)
  print(p) 
  dev.off()
  
  # Linear model summary  
  fileConn <-file(file.path(tablesdir, paste(section,mystk, "lm summary.txt", sep="_")))
  summary(lm(hat~tb,data=transform(wt2,tb=tb/mean(tb)))) %>% capture.output() %>% writeLines(., con=fileConn)
  close(fileConn)
  
  # Density dependence in maturity at age --------------------------------------
  
  p1=ggplot(mat)+
    geom_boxplot(aes(age,data))+
    xlab("Maturity-at-age")+ylab("Age")
  p2=ggplot(mat)+
    geom_point(aes(tb,rsdl))+
    geom_smooth(aes(tb,rsdl),method="lm")+
    facet_wrap(age~.)+
    xlab("Biomass")+ylab("Residual")
  
  p <-
    ggarrange(p1,p2 ,label.y=1,  
            widths = c(4,6), heights = c(1,1,1),
            nrow   = 1,          ncol    = 2)+
    theme_bw(16)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "dd_mat.jpg", sep="_")),
       width=10, height=8, units="in", res=300)
  print(p) 
  dev.off()
  
  # Estimating maturity ogive --------------------------------------------------
  
  matFn <-function(wt,x) 1/(1+exp(-x[1]*(wt-x[2])))
  fn    <-function(x,dat) sum((dat$mat-matFn(dat$wt,x))^2)
  
  dat=model.frame(FLQuants(stk,"wt"=stock.wt,"mat"=function(x) mat(x)))
  dat$age=factor(dat$age)

  matPar=optim(par=c(k=0.5,w50=.2), fn, dat=dat)
  save(matPar,file=file.path(dropboxdir,paste0("data/om/",mystk,"matPar.RData")))
 
  p <-
    ggplot(dat,aes(wt,mat))+
    geom_point(aes(wt,mat,col=age))+
    geom_smooth(se=F,span=0.2)+
    xlab("Mass-at-age")+ylab("Matutrity-at-age")+
    geom_line(aes(wt,mat),
              col="red",
              data=data.frame(mat=matFn(dat$wt,matPar$par),wt=dat$wt))

  # print(paste("k=", matPar$par[1]))
  # print(paste("W50=", matPar$par[2]))

  jpeg(filename=file.path(figuresdir, paste(section,mystk, "mat_ogive.jpg", sep="_")),
       width=10, height=8, units="in", res=300)
  print(p) 
  dev.off()

  # clean up and save
  rm(dat, gmr, m, mat, p, p1, p2, ts, wt, wt2)
  save(list=ls(),
       file = file.path(dropboxdir, "results", mystk, paste(mystk,"section2.RData", sep="_")))


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ----------------------------------------------------------------------------
  # 3. Density dependence and theoretical OM
  # ----------------------------------------------------------------------------
  
  load(file = file.path(dropboxdir, "results", mystk, paste(mystk,"section2.RData", sep="_")))
  
  section <- "03"
  
  # Get theoretical values from literature -------------------------------------
  par=FLife::teleost[c("linf","k","t0","l50","a","b"),dimnames(teleost)$iter==mylatin]
  
  # Uses life history theory to derive parameters for biological relationships
  par=FLife::lhPar(par,
                   m1=myparams[myparams$stock==mystk,"m1"],
                   m2=myparams[myparams$stock==mystk,"m2"],
                   m3=myparams[myparams$stock==mystk,"m3"])
  rm(param.)
  
  # WHAT IS THE ASSUMED SELECTIVITY AND WHERE DOES IT COME FROM?
  
  # derive a FLBRP equilibrium object with Lorenzen natural mortality
  eq =FLife::lhEql(par,m=myparams[myparams$stock==mystk,"m"])

  # plot the BRP object --------------------------------------------------------
  p <-ggplotFL::plot(eq) + 
           theme_bw() +
           labs(title=paste("FLBRP object"))
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "brp.jpg", sep="_")),
       width=10, height=8, units="in", res=300)
  print(p) 
  dev.off()
  
  # convert FLBRP to FLStock and project
  om=as(eq,"FLStock")
  om=fwd(om,fbar=fbar(om)[,-1],sr=eq)

  # update the par object
  par=rbind(par,
            FLPar(c(bref =c(refpts(eq)["msy","biomass"]),
                    delta= myparams[myparams$stock==mystk,"delta"],
                    matk = myparams[myparams$stock==mystk,"matk"],
                    w50  = par["a"]*par["l50"]^par["b"],
                    alpha=myparams[myparams$stock==mystk,"alpha"])))

  fileConn <-file(file.path(tablesdir, paste(section,mystk, "par.txt", sep="_")))
  as.data.frame(par) %>% dplyr::select(-iter) %>% pander::pandoc.table(style="simple") %>% capture.output() %>% writeLines(., con=fileConn)
  close(fileConn)

  # save par file
  # save(par,file=file.path(dropboxdir,paste0("data/om/par",mystk,".RData")))
  
  # Simulation of density dependence in mass-at-age ----------------------------
  # for different levels of biomass relative to Bmsy
  
  dat=ddWt(stock.wt(eq),biomass(eq),refpts(eq)["msy","biomass"]) 
  dat=merge(as.data.frame(dat,drop=T),
            as.data.frame(biomass(eq)%/%refpts(eq)["msy","biomass"],drop=T),
            by="year")
  names(dat)[3:4]=c("mass","biomass")
  
  p1 <- ggplot(dat)+    
           theme_bw() + theme(legend.position = "none") +
           geom_line(aes(age,mass,col=biomass,group=year))+
           scale_x_continuous(limits=c(0,maxage), breaks=(seq(0,maxage,1)))+
           # scale_y_continuous(limits=c(0,800))+
           geom_line(aes(age,data),data=as.data.frame(stock.wt(eq)),col="red",size=2)+
           labs(x="Age", y="Mass-at-age", colour="Biomass/Bmsy") 
  
  # jpeg(filename=file.path(figuresdir, paste(section,mystk, "sim_dd_mass.jpg", sep="_")),
  #      width=10, height=10, units="in", res=300)
  # print(p) 
  # dev.off()
  
  # Simulation of density dependence in maturity-at-age ------------------------ 
  # for different levels of biomass relative to Bmsy
  
  dat=transform(dat,mat=1/(1+exp(-par["matk"]*(mass-par["w50"]))))
  
  p2 <- 
    ggplot(dat)+           
    theme_bw() + theme(legend.position = "none") +
    geom_line(aes(age,mat,col=biomass,group=year))+
    scale_x_continuous(limits=c(0,6), breaks=(seq(0,maxage,1)))+
    geom_line(aes(age,data),data=as.data.frame(mat(eq)),col="red",size=2)+
    labs(x="Age", y="Maturity-at-age", colour="Biomass/Bmsy")
  
  # jpeg(filename=file.path(figuresdir, paste(section,mystk, "sim_dd_maturity.jpg", sep="_")),
  #      width=10, height=10, units="in", res=300)
  # print(p) 
  # dev.off()
  
  # Simulation of density dependence in M-at-age -------------------------------
  # for different levels of biomass relative to Bmsy
  
  dat=transform(dat,m=c(par["m1"])*(mass^c(par["m2"])))

  p3 <-
    ggplot(dat)+           
    theme_bw() + theme(legend.position = c(0.8,0.8)) +
    geom_line(aes(age,m,col=biomass,group=year))+
    scale_x_continuous(limits=c(0,maxage), breaks=(seq(0,maxage,1)))+
    geom_line(aes(age,data),data=as.data.frame(m(eq)),col="red",size=2)+
    labs(x="Age", y="M-at-age", colour="Biomass/Bmsy")

  p <-
    ggarrange(p1,p2, p3, ncol = 3)+
    theme_bw(16)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "sim_dd3.jpg", sep="_")),
       width=10, height=6, units="in", res=300)
  print(p) 
  dev.off()
  
  # Comparison of equilibrium yields curves without DD and with DD mass --------
  
  x         = propagate(eq,length(dimnames(eq)$year))
  fbar(x)   = fbar(x)[,1]
  fbar(x)[] = c(fbar(eq))*2
  
  # landing weight at age and discard weight at age relative to stock weight
  lsr       = landings.wt(eq)/stock.wt(eq)
  dsr       = discards.wt(eq)/stock.wt(eq)
  
  # iterate
  for (i in seq(10)){
    stock.wt(   x)=ddWt(stock.wt(eq),
                        biomass(x),
                        refpts(eq)["msy","biomass"],
                        delta=myparams[myparams$stock==mystk,"alpha"])
    landings.wt(x)=stock.wt(x)*lsr
    discards.wt(x)=stock.wt(x)*dsr
    x=brp(x)}
  
  df <- 
    bind_rows(
      model.frame(FLQuants(eq,ssb    =function(x) ssb(x),
                              biomass=function(x) biomass(x), 
                              catch  =function(x) catch(x),
                              f      =function(x) fbar(x)),
                           drop=T) %>% mutate(scen="no-DD"),
      model.frame(FLQuants(x, ssb    =function(x) ssb(x),
                              biomass=function(x) biomass(x), 
                              catch  =function(x) catch(x),
                              f      =function(x) fbar(x)),
                  drop=T) %>% mutate(scen="DD mass")
    )

  # Comparison of equilibrium yields curves with DD mass+maturity
  
  x=propagate(eq,length(dimnames(eq)$year))
  fbar(x)[,1]=c(seq(0,max(fbar(eq)),length.out=51),seq(1,50,length.out=51)*max(fbar(eq)))[-52]
  fbar(x)=fbar(x)[,1]
  
  lsr=landings.wt(eq)/stock.wt(eq)
  dsr=discards.wt(eq)/stock.wt(eq)
  # iterate
  for (i in seq(10)){
    stock.wt(   x)=ddWt(stock.wt(eq),
                        biomass(x),
                        refpts(eq)["msy","biomass"],
                        delta=myparams[myparams$stock==mystk,"alpha"])
    landings.wt(x)=stock.wt(x)*lsr
    discards.wt(x)=stock.wt(x)*dsr
    
    mat(x) =1/(1+exp(-par["matk"]*(stock.wt(x)-par["w50"])))
    x=brp(x)
  }
  
  df <- 
    bind_rows(
      df,
      model.frame(FLQuants(x,  ssb    =function(x) ssb(x),
                               biomass=function(x) biomass(x), 
                               catch  =function(x) catch(x),
                               f      =function(x) fbar(x)),
                  drop=T) %>% mutate(scen="DD mass+mat")
    )
 
  # Comparison of equilibrium yields curves with DD mass+mat+M 
  
  x=propagate(eq,length(dimnames(eq)$year))
  fbar(x)[,1]=c(seq(0,max(fbar(eq)),length.out=51),seq(1,50,length.out=51)*max(fbar(eq)))[-52]
  fbar(x)=fbar(x)[,1]
  
  lsr=landings.wt(eq)/stock.wt(eq)
  dsr=discards.wt(eq)/stock.wt(eq)
  # iterate
  for (i in seq(10)){
    stock.wt(   x)=ddWt(stock.wt(eq),biomass(x),refpts(eq)["msy","biomass"],delta=myparams[myparams$stock==mystk,"alpha"])
    landings.wt(x)=stock.wt(x)*lsr
    discards.wt(x)=stock.wt(x)*dsr
    
    mat(x) = 1/(1+exp(-par["matk"]*(stock.wt(x)-par["w50"])))
    m(x)   = par["m1"]%*%(stock.wt(x)%^%par["m2"])
    x      = brp(x)
  }
  
  df <- 
    bind_rows(
      df,
      model.frame(FLQuants(x, ssb    =function(x) ssb(x),
                           biomass=function(x) biomass(x), 
                           catch  =function(x) catch(x),
                           f      =function(x) fbar(x)),
                  drop=T) %>% mutate(scen="DD mass+mat+M")
    )
  save(df,file=file.path(dropboxdir,"data", "om", paste(mystk, "Eq.RData", sep="_")))

  p <-
    df %>% 
    mutate(scen = factor(scen, levels=c("no-DD", "DD mass", "DD mass+mat", "DD mass+mat+M"))) %>% 
    ggplot(aes(biomass, catch))+
    theme_bw() +
    geom_line(aes(colour=scen)) +
    labs(x="Total biomass", y="Catch", colour="") +
    scale_colour_manual(values =c("no-DD"         = "black",
                                  "DD mass"       = "red", 
                                  "DD mass+mat"   = "blue",
                                  "DD mass+mat+M" = "darkgreen")) 
  
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "eq_yield.jpg", sep="_")),
       width=10, height=6, units="in", res=300)
  print(p) 
  dev.off()
  
  # Check 1 --------------------------------------------------------------------

  # ptm <- proc.time()
  # 
  # # generate FLquants with different F values
  # f = FLQuant(c(rep(seq(0, 1,length.out=51)    ,each=101),
  #               rep(seq(1,25,length.out=61)[-1],each=101)),
  #             dimnames=list(year=dimnames(om)$year,
  #                           iter=seq(101)))  %*%
  #   refpts(eq)["msy","harvest"]
  # 
  # # Add F's to eq object
  # fbar(eq)=f
  # eq=brp(eq)
  # 
  # om=as(eq,"FLStock")
  # om=fwd(om,f=f[,-1],sr=eq)
  # 
  # # simulate with DD for 100 years into the future (THIS TAKES ABOUT 10 MINUTES)
  # om2=om
  # 
  # 
  # for (year in dimnames(om2)$year[-1]){
  #   om2=ddFn(year,om2,par)
  #   om2=fwd(om2,f=f[,year],sr=eq)
  # }
  # catch(om2)=computeCatch(om2,slot="all")
  # 
  # # Stop the clock
  # proc.time() - ptm
  # 
  # # plot(window(FLStocks("OM"=om,DD=om2),start=40))  # DIFFICULT TO INTERPRET
  # 
  # t <-
  #   bind_rows(
  #     model.frame(FLQuants(om2[,dim(m(om2))[2]],"biomass"=biomass,"catch"=catch)) %>% mutate(scen="DD sim"),
  #     model.frame(FLQuants(om[,dim(m(om))[2]],"biomass"=biomass,"catch"=catch)) %>% mutate(scen="no-DD sim"),
  #     filter(df, scen=="DD mass+mat+M") %>% mutate(scen="DD equilibrium")
  #   )
  # 
  # p <-
  #   ggplot(data=t) +
  #   geom_line(aes(biomass,catch, col=scen))+
  #   labs(title="Check on equilibrium and simulation results") +
  #   scale_colour_manual(values =c("no-DD sim"         = "black",
  #                                 "DD sim"       = "red",
  #                                 "DD equilibrium"   = "blue"))


  # clean up and save
  rm(dat, df, p, p1, p2, p3, t, x)
  save(list=ls(),
       file = file.path(dropboxdir, "results", mystk, paste(mystk,"section3.RData", sep="_")))
  # load(file = file.path(dropboxdir, "results", mystk, paste(mystk,"section3.RData", sep="_")))
  


  
  
  
  
  
  
  
  
  
  
  
  
  
  # ------------------------------------------------------------------------------
  # 4. OM based on ICES assessment
  # ------------------------------------------------------------------------------
  # load(file.path(dropboxdir, "data/om/parMac.RData"))
  # par_laurie <- par
  
  # rm(list=ls()[!grepl("dropboxdir|figuresdir|tablesdir|mystk",ls())])
  # load(file = file.path(dropboxdir, "results", mystk, paste(mystk,"section0.RData", sep="_")))
  load(file = file.path(dropboxdir, "results", mystk, paste(mystk,"section3.RData", sep="_")))

  # bind_rows(
  #   as.data.frame(par_laurie) %>% mutate(par="par_laurie"),
  #   as.data.frame(par) %>% mutate(par="par")
  # )  

  rm(eq, control)
  
  section <- "04"
  
  ices=window(stk,start=myparams[myparams$stock==mystk,"minyear"])
  stock(ices) =computeStock(ices)
  name(ices) <- mystkname

  assess_df <-
    model.frame(FLQuants(ices, "Biomass"=function(x) biomass(x),
                               "SSB"=function(x) ssb(x),
                               "F"=function(x) fbar(x),
                               "Yield"=function(x) catch(x),
                               "Rec"=function(x) rec(x)),drop=T) %>% mutate(scen="sam")
  
  # Start with the SAM assessment in forward projection
  #@# years changed and single iter to get rid of noise in 2020 @#@#@#@#@#@#@#@# 
  sam_sr  = fmle(as.FLSR(window(iter(ices,1), start=2000, end=2019),model="bevholtSV"),
                 fixed=list(s=myparams[myparams$stock==mystk,"steepness"], 
                            spr0=mean(spr0(ices))), 
                 control=list(silent=TRUE))
  sam_eq  = FLBRP(ices,ab(sam_sr))
  sam     = fwdWindow(ices,sam_eq,end=2050) 

  # plot(sam_sr)
  
  F  =propagate(window(fbar(sam),start=2020),101)
  F[]=rep(c(seq(0,                                   c(refpts(sam_eq)["msy",  "harvest"]),length.out=51),
            seq(c(refpts(sam_eq)["msy",  "harvest"]),c(refpts(sam_eq)["crash","harvest"])*1.2,length.out=51)[-1]),each=dim(F)[2])
  # F[]=rep(c(seq(0,                               c(refpts(eq)["crash","harvest"])*1.2,length.out=101)),each=dim(F)[2])
  
  sam_control=as(FLQuants("f"=F),"fwdControl")
  sam        =fwd(sam,
                  control=sam_control,
                  sr=sam_eq)
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  
  df_om  =model.frame(FLQuants(sam[,"2050"], "Biomass"=function(x) biomass(x),
                                    "SSB"=function(x) ssb(x),
                                    "F"=function(x) fbar(x),
                                    "Yield"=function(x) catch(x)),drop=T) %>% mutate(scen="sam")
                
  
  # VPA assessment based on SAM assessment -------------------------------------
  
  # run VPA
  vpa =ices+FLAssess:::VPA(ices)
  stock(vpa) = computeStock(vpa)
  name(vpa)  = paste(mystkname, "VPA")
  
  # add to assess_df
  assess_df <-
    bind_rows(
      assess_df,
      model.frame(FLQuants(vpa, "Biomass"=function(x) biomass(x),
                                "SSB"=function(x) ssb(x),
                                "F"=function(x) fbar(x),
                                "Yield"=function(x) catch(x),
                                "Rec"=function(x) rec(x)),drop=T) %>% mutate(scen="vpa") )
  
  #@# years changed and single iter to get rid of noise in 2020 @#@#@#@#@#@#@#@# 
  vpa_sr  =fmle(as.FLSR(window(iter(vpa,1),start=2000, end=2019),
                        model="bevholtSV"),
            fixed=list(s   =myparams[myparams$stock==mystk,"steepness"], spr0=mean(spr0(vpa))), 
            control=list(silent=TRUE))
  vpa_eq  = FLBRP(vpa,ab(vpa_sr))
  vpa     = fwdWindow(vpa,vpa_eq,end=2050) 
  
  # plot(vpa_sr)
  
  F  =propagate(window(fbar(vpa),start=2020),101)
  F[]=rep(c(seq(0,                                   c(refpts(vpa_eq)["msy",  "harvest"]),length.out=51),
            seq(c(refpts(vpa_eq)["msy",  "harvest"]),c(refpts(vpa_eq)["crash","harvest"])*1.2,length.out=51)[-1]),each=dim(F)[2])
  # F[]=rep(c(seq(0,                               c(refpts(eq)["crash","harvest"])*1.2,length.out=101)),each=dim(F)[2])
  
  vpa_control = as(FLQuants("f"=F),"fwdControl")
  vpa         = fwd(vpa,
                    control=vpa_control,
                    sr=vpa_eq)
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  
  df_om  =bind_rows(
    df_om,
    model.frame(FLQuants(vpa[,"2050"], "Biomass"=function(x) biomass(x),
                                       "SSB"=function(x) ssb(x),
                                       "F"=function(x) fbar(x),
                                       "Yield"=function(x) catch(x)),drop=T) %>% mutate(scen="vpa"))

  # VPA assessment with age varying M -----------------------------------------------------


  # RUN NEXT SECTION ONLY ONCE !!!!!! =====================================================
  # load(file.path(dropboxdir, paste0("data/om/par",mystk,".RData")))
  vpaM_par <- par
  vpaM_par["m1"]  =par["m1"] /myparams[myparams$stock==mystk,"m1scaler"]    
  vpaM_par["w50"] =myparams[myparams$stock==mystk,"w50"]
  vpaM_par["matk"]=myparams[myparams$stock==mystk,"matk2"]               
  # save(par, file=file.path(dropboxdir, paste0("data/om/par",mystk,".RData")))
  # =======================================================================================
  

  fileConn <-file(file.path(tablesdir, paste(section,mystk, "par_om.txt", sep="_")))
  as.data.frame(vpaM_par) %>% dplyr::select(-iter) %>% pander::pandoc.table(style="simple") %>% capture.output() %>% writeLines(., con=fileConn)
  close(fileConn)
  
  vpaM        = ices
  m(vpaM)     = vpaM_par["m1"]%*%(stock.wt(vpaM)%^%vpaM_par["m2"])
  vpaM        = vpaM+FLAssess:::VPA(vpaM)
  stock(vpaM) = computeStock(vpaM)
  name(vpaM)  = paste(mystkname, "VPA-M")
  
  vpaM_hist   = vpaM

  # add to assess df
  assess_df <-
    bind_rows(
      assess_df,
      model.frame(FLQuants(vpaM, "Biomass"=function(x) biomass(x),
                                 "SSB"=function(x) ssb(x),
                                 "F"=function(x) fbar(x),
                                 "Yield"=function(x) catch(x),
                                 "Rec"=function(x) rec(x)),drop=T) %>% mutate(scen="vpaM") )

  #@# years changed and single iter to get rid of noise in 2020 @#@#@#@#@#@#@#@# 
  # stock recruitment estimation  
  vpaM_sr  = fmle(as.FLSR(window(iter(vpaM,1), start=2000, end=2019),
                         model="bevholtSV"),
                fixed=list(s=myparams[myparams$stock==mystk,"steepness"], spr0=mean(spr0(window(iter(vpaM,1),start=2000,end=2019)))), 
                control=list(silent=TRUE))
  
  vpaM_eq  = FLBRP(window(iter(vpaM,1),start=2000,end=2019),sr=ab(vpaM_sr))
  vpaM     = fwdWindow(window(vpaM,end=2020),vpaM_eq,end=2050) 

  F  =propagate(window(fbar(vpaM),start=2020),101)
  F[]=rep(c(seq(0,                                    c(refpts(vpaM_eq)["msy",  "harvest"]),length.out=51),
            seq(c(refpts(vpaM_eq)["msy",  "harvest"]),c(refpts(vpaM_eq)["crash","harvest"])*1.5,length.out=51)[-1]),each=dim(F)[2])

  control_Frange = as(FLQuants("f"=F),"fwdControl")
  vpaM           = fwd(vpaM,
                     control=control_Frange,
                     sr=vpaM_eq)
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  
  d1=model.frame(FLQuants(vpaM[,"2050"],
                          "Biomass"=function(x) biomass(x),
                          "SSB"=function(x) ssb(x),
                          "F"=function(x) fbar(x),
                          "Yield"=function(x) catch(x)),drop=T)
  
  
  # plot of SRR
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "om_srr.jpg", sep="_")),
       width=10, height=12, units="in", res=300)
  print(plot(vpaM_sr)) 
  dev.off()
  
  # plot if BRP
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "om_brp.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  print(plot(vpaM_eq)) 
  dev.off()
  
  
  df_helper <-
    df_om %>% 
    group_by(scen) %>% 
    summarise(Yield = max(Yield, na.rm=TRUE)) %>% 
    left_join(df_om)
  
  # Plot of yield vs biomass
  p1<-
    df_om %>% 
    filter(scen %in% c("sam","vpa","vpa-M")) %>% 
    
    ggplot(aes(x=Biomass, y=Yield)) +
    theme_bw() +
    theme(legend.position="none") +
    geom_line(aes(colour=scen)) +
    geom_segment(data=df_helper %>% filter(scen %in% c("sam","vpa","vpa-M")),
                 aes(x=Biomass, xend=Biomass, y=0, yend=Yield, colour=scen)) +
    ggrepel::geom_text_repel(data=df_helper %>% filter(scen %in% c("sam","vpa","vpa-M")),
                             aes(y=0, x=Biomass, label=format(Biomass/1000000,digits=2,nsmall=1), colour=scen),
                             min.segment.length = 0, show.legend=FALSE) +
    
    geom_segment(data=df_helper %>% filter(scen %in% c("sam","vpa","vpa-M")),
                 aes(x=0, xend=Biomass, y=Yield, yend=Yield, colour=scen)) +
    ggrepel::geom_text_repel(data=df_helper %>% filter(scen %in% c("sam","vpa","vpa-M")),
                             aes(x=0, y=Yield, label=format(Yield/1000000,digits=2,nsmall=1), colour=scen),
                             min.segment.length = 0, show.legend=FALSE)

  # jpeg(filename=file.path(figuresdir, paste(section,mystk, "om_biomass_yield.jpg", sep="_")),
  #      width=10, height=10, units="in", res=300)
  # print(p) 
  # dev.off()

  # Plot of F vs biomass
  p2 <-
    df_om %>% 
    filter(F < 1.5) %>% 
    filter(scen %in% c("sam","vpa","vpa-M")) %>% 
    
    ggplot(aes(x=Biomass, y=F)) +
    theme_bw() +
    theme(legend.position=c(0.8, 0.8)) +
    geom_line(aes(colour=scen)) +
    geom_segment(data=df_helper %>% filter(scen %in% c("sam","vpa","vpa-M")),
                 aes(x=Biomass, xend=Biomass, y=0, yend=F, colour=scen)) +
    ggrepel::geom_text_repel(data=df_helper %>% filter(scen %in% c("sam","vpa","vpa-M")),
                             aes(y=0, x=Biomass, label=format(Biomass/1000000,digits=2,nsmall=1), colour=scen),
                             min.segment.length = 0, show.legend=FALSE) +
    
    geom_segment(data=df_helper %>% filter(scen %in% c("sam","vpa","vpa-M")),
                 aes(x=0, xend=Biomass, y=F, yend=F, colour=scen)) +
    ggrepel::geom_text_repel(data=df_helper %>% filter(scen %in% c("sam","vpa","vpa-M")),
                             aes(x=0, y=F, label=format(F,digits=2,nsmall=1), colour=scen),
                             min.segment.length = 0, show.legend=FALSE)
  
  p <- ggarrange(p1,p2, ncol=2, nrow=1)
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "om_biomass_yield_F.jpg", sep="_")),
       width=10, height=6, units="in", res=300)
  print(p) 
  dev.off()
  
  # Comparing assessments
  p <-
    assess_df %>%
    tidyr::pivot_longer(names_to = "variable", values_to = "data", Biomass:Rec) %>% 
    dplyr::mutate(variable = factor(variable, levels=c("Yield","Rec","F","SSB","Biomass"))) %>% 
    ggplot(aes(x=year, y=data, group=scen)) +
    theme_bw() +
    geom_line(aes(colour=scen)) +
    expand_limits(y=0) +
    facet_wrap(~variable, ncol=2, scales="free_y")
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "comparing_assessments.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  print(p)
  dev.off()
  
  df_om  =bind_rows(
    df_om,
    d1 %>% mutate(scen="vpa-M"))
  
  df_helper <-
    df_om %>% 
    group_by(scen) %>% 
    summarise(Yield = max(Yield, na.rm=TRUE)) %>% 
    left_join(df_om)
  

  ##############################################################################
  # Setting Bref at BMSY!!
  
  vpaM_par["bref"]=subset(d1, Yield==max(Yield))[,"Biomass"]
  # vpaM_par["bref"]=refpts(vpaM_eq)["msy","biomass"]
  # save(vpaM_par, file=file.path(dropboxdir, paste0("data/om/par",mystk,".RData")))
  vpaM_par[c("m1","m2","delta","bref","matk","w50")]
  ##############################################################################

  # OM with DD in mass ---------------------------------------------------------
  
  # Start the clock
  ptm <- proc.time()
  
  vpaDDM=vpaM
  # mat(vpaDDMM)[1] <-0

  i <- ac(2020)  
  for (i in ac(2020:2050)) {
    vpaDDM =ddFn(i,vpaDDM,vpaM_par,massFlag=TRUE,matFlag=FALSE,mFlag=FALSE)
    control=as(FLQuants("f"=F[,i]),"fwdControl")
    vpaDDM =fwd(vpaDDM,
                control=control,
                sr=vpaM_eq)
  }
  
  # make data.frame
  d2=model.frame(FLQuants(vpaDDM[,"2050"],
                          "Biomass"=function(x) biomass(x),
                          "SSB"=function(x) ssb(x),
                          "F"=function(x) fbar(x),
                          "Yield"=function(x) catch(x)),drop=T)
  
  # Four panel plot
  p1=ggplot(stock.wt(vpaDDM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(stock.wt(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mass-at-age")
  p2=ggplot(mat(vpaDDM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(mat(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mat-at-age")
  p3=ggplot(m(vpaDDM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(m(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("M-at-age")
  p4=ggplot(d2)+
    geom_line(aes(SSB,Yield),col="black")+
    geom_line(aes(SSB,Yield),col="red",data=d1)
  
  p <- ggarrange(p1,p2,p3,p4, ncol=2, nrow=2)
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "VPADDM_4panels.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  print(p) 
  dev.off()
  
  # add to OM data.frame
  df_om  =bind_rows(
    df_om,
    d2 %>% mutate(scen="vpa-DDM"))
  
  # Stop the clock
  proc.time() - ptm
  
  # OM with DD in mass and maturity --------------------------------------------
  
  # Start the clock
  ptm <- proc.time()
  
  vpaDDMM =vpaM
  
  for (i in ac(2020:2050)) {
    vpaDDMM =ddFn(i,vpaDDMM,vpaM_par,massFlag=TRUE,matFlag=TRUE,mFlag=FALSE)
    control=as(FLQuants("f"=F[,i]),"fwdControl")
    vpaDDMM =fwd(vpaDDMM,
                 control=control,
                 sr=vpaM_eq)
  }
  
  # make data.frame
  d3=model.frame(FLQuants(vpaDDMM[,"2050"],
                          "Biomass"=function(x) biomass(x),
                          "SSB"=function(x) ssb(x),
                          "F"=function(x) fbar(x),
                          "Yield"=function(x) catch(x)),drop=T)
  
  # Four panel plot
  p1=ggplot(stock.wt(vpaDDMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(stock.wt(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mass-at-age")
  p2=ggplot(mat(vpaDDMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(mat(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mat-at-age")
  p3=ggplot(m(vpaDDMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(m(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("M-at-age")
  p4=ggplot(d3)+
    geom_line(aes(SSB,Yield),col="black")+
    geom_line(aes(SSB,Yield),col="red",data=d1)
  
  p <- ggarrange(p1,p2,p3,p4, ncol=2, nrow=2)
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "VPADDMM_4panels.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  print(p) 
  dev.off()
  
  # add to OM data.frame
  df_om  =bind_rows(
    df_om,
    d3 %>% mutate(scen="vpa-DDMM"))
  
  # Stop the clock
  proc.time() - ptm
  
  # plot(mat(vpaDDMM)["2"])
  # plot(stock.wt(vpaDDMM)["2"])
  
  # OM with DD in mass and maturity and M --------------------------------------
  
  # Start the clock
  ptm <- proc.time()
  
  vpaDDMMM=vpaM
  
  for (i in ac(2020:2050)) {
    vpaDDMMM =ddFn(i,vpaDDMMM,vpaM_par,massFlag=TRUE,matFlag=TRUE,mFlag=TRUE)
    control=as(FLQuants("f"=F[,i]),"fwdControl")
    vpaDDMMM =fwd(vpaDDMMM,control=control,sr=vpaM_eq)
  }

  # Make data.frame
  d4=model.frame(FLQuants(vpaDDMMM[,"2050"],
                          "Biomass"=function(x) biomass(x),
                          "SSB"=function(x) ssb(x),
                          "F"=function(x) fbar(x),
                          "Yield"=function(x) catch(x)),drop=T)
  
  # Four panel plot
  p1=ggplot(stock.wt(vpaDDMMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(stock.wt(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mass-at-age")
  p2=ggplot(mat(vpaDDMMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(mat(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mat-at-age")
  p3=ggplot(m(vpaDDMMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(m(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("M-at-age")
  p4=ggplot(d4)+
    geom_line(aes(SSB,Yield),col="black")+
    geom_line(aes(SSB,Yield),col="red",data=d1)
  
  p <- ggarrange(p1,p2,p3,p4, ncol=2, nrow=2)
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "VPADDMMM_4panels.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  print(p) 
  dev.off()
  
  # add to OM data.frame
  df_om  =bind_rows(
    df_om,
    d4 %>% mutate(scen="vpa-DDMMM"))

  # approximate reference points
  df_helper <-
    df_om %>% 
    group_by(scen) %>% 
    summarise(Yield = max(Yield, na.rm=TRUE)) %>% 
    left_join(df_om)
  
  
  # Stop the clock
  proc.time() - ptm
  
  # Plotting overview scenarios ========================================================
  
  # scenarios
  scenarios = c("vpa-M", "vpa-DDM","vpa-DDMM","vpa-DDMMM") 
  
  # Plot of yield vs biomass
  p1<-
    df_om %>% 
    filter(scen %in% scenarios) %>% 
    mutate(scen = factor(scen, levels=scenarios)) %>% 
    
    ggplot(aes(x=Biomass, y=Yield)) +
    theme_bw() +
    theme(legend.position="none") +
    geom_line(aes(colour=scen)) +
    geom_segment(data=df_helper %>% filter(scen %in% scenarios) %>% mutate(scen = factor(scen, levels=scenarios)),
                 aes(x=Biomass, xend=Biomass, y=0, yend=Yield, colour=scen)) +
    ggrepel::geom_text_repel(data=df_helper %>% filter(scen %in% scenarios),
                             aes(y=0, x=Biomass, label=format(Biomass/1000000,digits=2,nsmall=1), colour=scen),
                             min.segment.length = 0, show.legend=FALSE) +
    
    geom_segment(data=df_helper %>% filter(scen %in% scenarios) %>% mutate(scen = factor(scen, levels=scenarios)),
                 aes(x=0, xend=Biomass, y=Yield, yend=Yield, colour=scen)) +
    ggrepel::geom_text_repel(data=df_helper %>% filter(scen %in% scenarios),
                             aes(x=0, y=Yield, label=format(Yield/1000000,digits=2,nsmall=1), colour=scen),
                             min.segment.length = 0, show.legend=FALSE)
  
  # jpeg(filename=file.path(figuresdir, paste(section,mystk, "om_biomass_yield.jpg", sep="_")),
  #      width=10, height=10, units="in", res=300)
  # print(p) 
  # dev.off()
  
  # Plot of F vs biomass
  p2 <-
    df_om %>% 
    filter(F < 1.5) %>% 
    filter(scen %in% scenarios) %>% 
    mutate(scen = factor(scen, levels=scenarios)) %>% 
    
    ggplot(aes(x=Biomass, y=F)) +
    theme_bw() +
    theme(legend.position=c(0.8, 0.8)) +
    geom_line(aes(colour=scen)) +
    geom_segment(data=df_helper %>% 
                   filter(scen %in% scenarios) %>% 
                   mutate(scen = factor(scen, levels=scenarios)),
                 aes(x=Biomass, xend=Biomass, y=0, yend=F, colour=scen)) +
    ggrepel::geom_text_repel(data=df_helper %>% 
                               filter(scen %in% scenarios) %>% 
                               mutate(scen = factor(scen, levels=scenarios)),
                             aes(y=0, x=Biomass, label=format(Biomass/1000000,digits=2,nsmall=1), colour=scen),
                             min.segment.length = 0, show.legend=FALSE) +
    
    geom_segment(data=df_helper %>% 
                   filter(scen %in% scenarios) %>% 
                   mutate(scen = factor(scen, levels=scenarios)),
                 aes(x=0, xend=Biomass, y=F, yend=F, colour=scen)) +
    ggrepel::geom_text_repel(data=df_helper %>% 
                               filter(scen %in% scenarios) %>% 
                               mutate(scen = factor(scen, levels=scenarios)),
                             aes(x=0, y=F, label=format(F,digits=2,nsmall=1), colour=scen),
                             min.segment.length = 0, show.legend=FALSE)
  
  p <- ggarrange(p1,p2, ncol=2, nrow=1)
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "om_biomass_yield_F.jpg", sep="_")),
       width=10, height=6, units="in", res=300)
  print(p) 
  dev.off()
  
  oms <- FLStocks()
  oms[["Base"]] =vpaM
  oms[["M"]]    =vpaDDM
  oms[["MM"]]   =vpaDDMM
  oms[["MMM"]]  =vpaDDMMM
  
  rfs=ldply(oms, function(x) { 
    model.frame(FLQuants(x[,"2050"],
                         biomass=function(x) stock(x),
                         ssb    =function(x) ssb(  x),
                         catch  =function(x) catch(x),
                         f      =function(x) fbar( x)),drop=TRUE)[,-1]})
  
  # reference points (at maximum catch)
  rfpts=ddply(rfs,.(.id), with, {  
    flag=catch==max(catch)
    data.frame(Bmsy  =biomass[flag],
               SSBmsy=ssb[flag],
               MSY   =catch[flag],
               Fmsy  =f[flag],
               Virgin=max(ssb),
               B0    =max(biomass))}) %>% 
    dplyr::rename(scen=.id)
  
  # Make table
  fileConn <-file(file.path(tablesdir, paste(section,mystk, "refpoints.txt", sep="_")))
  rfpts %>% 
    pander::pandoc.table(style="simple", big.mark=",", justify="left", split.tables=400) %>% 
    capture.output() %>% writeLines(., con=fileConn)
  close(fileConn)
  

  # clean up and save
  rm(dat, df, p, p1, p2, p3, p4, t, x, d1, d2, d3, d4, x)
  save(list=ls(),
       file = file.path(dropboxdir, "results", mystk, paste(mystk,"section4.RData", sep="_")))
#  load(file = file.path(dropboxdir, "results", mystk, paste(mystk,"section4.RData", sep="_")))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ------------------------------------------------------------------------------
  # 5. MSE
  # ------------------------------------------------------------------------------
  
  load(file = file.path(dropboxdir, "results", mystk, paste(mystk,"section4.RData", sep="_")))
  section <- "05"
  

  om=fwdWindow(window(vpaM,end=2020),end=2050,vpaM_eq)
  F    =rep(c(seq(0,                                    c(refpts(vpaM_eq)["msy",  "harvest"]),  length.out=51),
              seq(c(refpts(vpaM_eq)["msy",  "harvest"]),c(refpts(vpaM_eq)["crash","harvest"])*1.5,length.out=51)[-1]))
  F    =FLQuant(rep(F,each=31),dimnames=list(year=2020:2050,iter=seq(101)))
  om=fwd(om,fbar=F,sr=vpaM_eq)
  
  
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  #@# Projection for range of F                                              @#@
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  oms=FLStocks("Base"     =fwd(propagate(fwdWindow(window(iter(om,1),end=2020),end=2050,vpaM_eq),101),fbar=F,sr=vpaM_eq))
  oms[["DD Mass"]]        =    propagate(fwdWindow(window(iter(om,1),end=2020),end=2050,vpaM_eq),101)
  oms[["DD Mass, Mat"]]   =    propagate(fwdWindow(window(iter(om,1),end=2020),end=2050,vpaM_eq),101)
  oms[["DD Mass, Mat, M"]]=    propagate(fwdWindow(window(iter(om,1),end=2020),end=2050,vpaM_eq),101)
  for (iYear in ac(2020:2050)){
    oms[["DD Mass"]]=ddFn(iYear,oms[["DD Mass"]],vpaM_par,TRUE,FALSE,FALSE)
    oms[["DD Mass"]]=fwd(oms[["DD Mass"]],fbar=F[,iYear],sr=vpaM_eq)
    
    oms[["DD Mass, Mat"]]=ddFn(iYear,oms[["DD Mass, Mat"]],vpaM_par,TRUE,TRUE,FALSE)
    oms[["DD Mass, Mat"]]=fwd(oms[["DD Mass, Mat"]],fbar=F[,iYear],sr=vpaM_eq)
    
    oms[["DD Mass, Mat, M"]]=ddFn(iYear,oms[["DD Mass, Mat, M"]],vpaM_par,TRUE,TRUE,TRUE)
    oms[["DD Mass, Mat, M"]]=fwd(oms[["DD Mass, Mat, M"]],fbar=F[,iYear],sr=vpaM_eq)
  } 
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  # reference points data.frame (equilibrium)
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  rfs=ldply(oms, function(x) { 
    model.frame(FLQuants(x[,"2050"],
                         biomass=function(x) stock(x),
                         ssb    =function(x) ssb(  x),
                         catch  =function(x) catch(x),
                         f      =function(x) fbar( x)),drop=TRUE)[,-1]})
  
  # reference points (at maximum catch)
  rfpts=ddply(rfs,.(.id), with, {  
    flag=catch==max(catch)
    data.frame(Bmsy  =biomass[flag],
               SSBmsy=ssb[flag],
               MSY   =catch[flag],
               Fmsy  =f[flag],
               Virgin=max(ssb),
               B0    =max(biomass))}) %>% 
    dplyr::rename(scen=.id)
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  #@# projection with different Fs                                           @#@
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  prj=FLStocks("Base" =fwdWindow(window(iter(oms[[1]],1:4),end=2020),end=2050,vpaM_eq),
               "M"    =fwdWindow(window(iter(oms[[2]],1:4),end=2020),end=2050,vpaM_eq),
               "MM"   =fwdWindow(window(iter(oms[[3]],1:4),end=2020),end=2050,vpaM_eq),
               "MMM"  =fwdWindow(window(iter(oms[[4]],1:4),end=2020),end=2050,vpaM_eq))
  
  F=fbar(prj[[1]][,ac(2020:2050)])%=%rep(rfpts$Fmsy,each=31)
  
  control=as(FLQuants("f"=F),"fwdControl")
  prj[[1]] =fwd(prj[[1]],
                control=control,
                sr=vpaM_eq)
  
  for (i in ac(2020:2050)) {
    
    control=as(FLQuants("f"=F[,i]),"fwdControl")
    
    prj[[2]] =ddFn(i,prj[[2]],vpaM_par,massFlag=TRUE,matFlag=FALSE,mFlag=FALSE)
    prj[[2]] =fwd(prj[[2]],
                  control=control,
                  sr=vpaM_eq)
    
    prj[[3]] =ddFn(i,prj[[3]],vpaM_par,massFlag=TRUE,matFlag=TRUE,mFlag=FALSE)
    prj[[3]] =fwd(prj[[3]],
                  control=control,
                  sr=vpaM_eq)
    
    prj[[4]] =ddFn(i,prj[[4]],vpaM_par,massFlag=TRUE,matFlag=TRUE,mFlag=TRUE)
    prj[[4]] =fwd(prj[[4]],
                  control=control,
                  sr=vpaM_eq)}
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  #@# Compare                                                                @#@
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  
  # summary data.frame projections at different Fmsy  
  prj_df <-
    ldply(prj, function(x) {
      model.frame(FLQuants(x,
                           rec    =function(x) rec(x),
                           biomass=function(x) stock(x),
                           ssb    =function(x) ssb(  x),
                           catch  =function(x) catch(x),
                           f      =function(x) fbar( x)),
                  drop=TRUE)}
    )
  
  nyears = length(unique(prj_df$year))
  nscen  = length(unique(prj_df$.id))
  
  prj_df <- prj_df %>% 
    bind_cols(Fmsy    = rep(rep(rfpts$Fmsy,     each=(nyears)), nscen)) %>% 
    bind_cols(Fmsy_id = rep(rep(unique(prj_df$.id), each=(nyears)), nscen)) %>% 
    mutate(Fmsy = format(Fmsy, digits=2)) 

  # summary of weight, maturity and M data.frame projections at different Fmsy  
  prj_df2 <-
    ldply(prj, function(x) {
      model.frame(FLQuants(x,
                           stock.wt = function(x) stock.wt(x),
                           mat      = function(x) mat(x),
                           m        = function(x) m(x)),
                  drop=TRUE)}
    )
  
  nyears = length(unique(prj_df2$year))
  nscen  = length(unique(prj_df2$.id))
  nages  = length(unique(prj_df2$age))
  
  prj_df2 <- 
    prj_df2 %>% 
    bind_cols(Fmsy    = rep(rep(rfpts$Fmsy,     each=(nyears*nages)), nscen)) %>% 
    bind_cols(Fmsy_id = rep(rep(unique(prj_df2$.id), each=(nyears*nages)), nscen)) %>% 
    mutate(Fmsy = format(Fmsy, digits=2)) 
  

  # plot of deterministic projections
  p <-
    prj_df %>% 
    filter(.id == Fmsy_id) %>% 
    pivot_longer(names_to = "variable", values_to = "data", rec:f) %>% 
    mutate(variable = factor(variable, levels=c("rec","f","catch","ssb","biomass"))) %>% 
    ggplot(aes(x=year, y=data, group=Fmsy)) +
    theme_bw() +
    geom_line(aes(colour=Fmsy)) +
    expand_limits(y=0) +
    facet_grid(variable ~ .id, scales="free_y")
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "deterministic_projections_at_Fmsy.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  print(p) 
  dev.off()

  
  # plot of stock weights etc. in deterministic projections
  p <-
    prj_df2 %>% 
    filter(.id == Fmsy_id) %>% 
    filter(age %in% c(1:6)) %>% 
    pivot_longer(names_to = "variable", values_to = "data", stock.wt:m) %>% 
    mutate(variable = factor(variable, levels=c("stock.wt","mat","m"))) %>% 
    ggplot(aes(x=year, y=data, group=.id)) +
    theme_bw() +
    geom_line(aes(colour=.id)) +
    expand_limits(y=0) +
    labs(x="", y="", colour="scen") +
    facet_grid(variable ~ age, scales="free_y")
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "deterministic_projections_of_MMM.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  print(p) 
  dev.off()

  # check whether values are on the equilbrium curves (= prf_df)
  
  eqCurves=ldply(oms, function(x) {
    model.frame(FLQuants(x[,"2050"],
                         biomass=function(x) stock(x),
                         ssb    =function(x) ssb(  x),
                         catch  =function(x) catch(x),
                         f      =function(x) fbar( x)),drop=TRUE)[,-1]})
  
  ts=ldply(prj, function(x)
    model.frame(FLQuants(x,
                         ssb  =function(x) ssb(  x),
                         f    =function(x) fbar(x),
                         catch=function(x) catch(x)),drop=T))
  ts=transform(ts,What=.id)
  ts=transform(ts,.id=unique(eqCurves$.id)[an(iter)])
  

  ggplot()+
    geom_line(data=eqCurves,
              aes(ssb,catch,col=.id))+
    geom_point(data=subset(ts,year==2050),
               aes(ssb,catch,col=.id))

  ggplot()+
    geom_line(data=subset(ts,year==2050),
              aes(ssb,catch,col=.id))+
    geom_point(data=subset(prj_df,year==2050),
               aes(ssb,catch,col=.id))
  
  
  
  
  # =====================================================================
    
  # Next running the stochastic loop
  devRec=rlnorm(100,rec(prj[[1]])[,ac(2020:2050),,,,1]%=%0,0.3)
  base=fwdWindow(window(propagate(iter(prj[[1]],1),100),end=2020),end=2050,vpaM_eq)
  
  prjs=mlply(data.frame(F=rfpts$Fmsy), 
             .progress = "tk", 
             function(F) {

    ## FMSY estimate by OM ######################################
    F=FLQuant(F,dimnames=list(year=2020:2050))
    
    ## No DD ####################################################
    control=as(FLQuants("f"=F),"fwdControl")
    x1     =fwd(base,
                control=control,
                sr=vpaM_eq,
                residuals=devRec)
    
    ## With DD ####################################################
    x2=x1
    x3=x1
    x4=x1
    print("")
    for (iYr in ac(2020:2050)) {
      
      cat(paste(iYr," "))
      
      control=as(FLQuants("f"=F[,iYr]),"fwdControl")
      
      ## DD M ##################################################         
      x2 =ddFn(iYr,x2,vpaM_par,massFlag=TRUE,matFlag=FALSE,mFlag=FALSE)
      x2 =fwd(x2,
              control=control,
              sr=vpaM_eq,
              residuals=devRec)
      
      ## DD MM##################################################
      x3 =ddFn(iYr,x3,vpaM_par,massFlag=TRUE,matFlag=TRUE,mFlag=FALSE)
      x3 =fwd(x3,
              control=control,
              sr=vpaM_eq,
              residuals=devRec)
      
      ## DD MMM ################################################
      x4 =ddFn(iYr,x4,vpaM_par,massFlag=TRUE,matFlag=TRUE,mFlag=TRUE)
      x4 =fwd(x4,
              control=control,
              sr=vpaM_eq,
              residuals=devRec)
    }  ## year loop
    
    FLStocks(list("Base"=x1,"M"=x2,"MM"=x3,"MMM"=x4))        
    
  }) ## F loop

  # generating stochastic projection data frame
  prjs_df=ldply(prjs, function(x) ldply(x, function(x) { 
    model.frame(FLQuants(x,
                         rec    =function(x) rec(x),
                         biomass=function(x) stock(x),
                         ssb    =function(x) ssb(  x),
                         catch  =function(x) catch(x),
                         f      =function(x) fbar( x)),drop=TRUE)})) %>% 
    mutate(Fmsy=factor(F, levels=rfpts$Fmsy,labels=c("Base","M","MM","MMM"))) %>% 
    as_tibble()

  prjs_df2 <-
    prjs_df %>% 
    filter(year >= 2020) %>% 
    group_by(.id, iter, F, Fmsy) %>% 
    mutate(cumcatch = cumsum(catch)) %>% 
    ungroup() %>% 
    dplyr::select(.id, iter, F, Fmsy, year, cumcatch)


  probs <- c(0.025, 0.5, 0.975)
  # p_names <- map_chr(p, ~paste0(.x*100, "%"))
  probs_names <- c("lower","median","upper")
  probs_funs <- map(probs, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
    set_names(nm = probs_names)  
  
  prjs_df_summ <-
    prjs_df %>% 
    pivot_longer(names_to = "variable", values_to = "data", rec:f) %>% 
    group_by(.id, F, Fmsy, variable, year) %>% 
    summarize_at(vars(data), tibble::lst(!!!probs_funs)) %>% 
    mutate(F = format(F, digits=2)) 
    
  prjs_df2_summ <-
    prjs_df2 %>% 
    pivot_longer(names_to = "variable", values_to = "data", cumcatch) %>% 
    group_by(.id, F, Fmsy, variable, year) %>% 
    summarize_at(vars(data), tibble::lst(!!!probs_funs)) %>% 
    mutate(F = format(F, digits=2)) 
  

  
  # only with the right Fmsy
  p <-
    prjs_df_summ %>% 
    filter(.id == Fmsy) %>% 
    mutate(variable = factor(variable, levels=c("rec","f","catch","ssb","biomass"))) %>% 
    ggplot(aes(x=year, y=median)) +
    theme_bw() +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=Fmsy), alpha=0.3) +
    geom_line(aes(colour=Fmsy)) +
    labs(y="", x="", colour="Fmsy", fill="Fmsy") +
    expand_limits(y=0) +
    facet_grid(variable ~ .id, scales="free_y")

  jpeg(filename=file.path(figuresdir, paste(section,mystk, "stochastic_projections_by_scenario.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  print(p) 
  dev.off()

  # stochastic projections at 'right' Fmsy  
  p <-
    prjs_df_summ %>% 
    filter(.id == Fmsy) %>% 
    mutate(variable = factor(variable, levels=c("rec","f","catch","ssb","biomass"))) %>% 
    ggplot(aes(x=year, y=median)) +
    theme_bw() +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=.id), alpha=0.3) +
    geom_line(aes(colour=.id)) +
    labs(y="", x="", colour="scen", fill="scen") +
    expand_limits(y=0) +
    facet_grid(variable ~ ., scales="free_y")
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "stochastic_projections_at_Fmsy.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  print(p) 
  dev.off()
  
  p <-
    prjs_df_summ %>% 
    mutate(variable = factor(variable, levels=c("rec","f","catch","ssb","biomass"))) %>% 
    ggplot(aes(x=year, y=median)) +
    theme_bw() +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=F), alpha=0.3) +
    geom_line(aes(colour=F)) +
    labs(y="", x="", colour="Fmsy", fill="Fmsy") +
    expand_limits(y=0) +
    facet_grid(variable ~ .id, scales="free_y")

  # prjs_df_summ %>% filter(year==2050) %>% 
  #   mutate(specified=ifelse(.id==Fmsy, "correct","misspecified")) %>% 
  #   ggplot(aes(x=.id, y=median)) +
  #   theme_bw() + theme(legend.position="bottom") +
  #   geom_point(aes(colour=specified),
  #              position=position_dodge(.4)) +
  #   geom_errorbar(aes(ymin=lower, ymax=upper, colour=specified), 
  #                 width=0.2, position=position_dodge(.4)) +
  #   geom_text(aes(label=paste(Fmsy, F), colour=specified), 
  #             hjust=0, position=position_dodge(.4)) +
  #   expand_limits(y=0) +
  #   facet_wrap(~variable, scales="free_y")
  
    
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "stochastic_projections_at_different_Fmsys.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  print(p) 
  dev.off()

  prjs_df2050 <-
    prjs_df %>% 
    filter(year == 2050) %>% 
    filter(.id == Fmsy) %>% 
    pivot_longer(names_to = "variable", values_to = "data", rec:f) %>% 
    mutate(variable = factor(variable, levels=c("rec","f","catch","ssb","biomass")))
  
  r <-
    prjs_df2050 %>% 
    filter(variable != "f") %>% 
    group_by(year, variable) %>% 
    summarise(
      min=min(data, na.rm=TRUE), 
      max=max(data, na.rm=TRUE)) %>% 
    mutate(binwidth = (max-min)/10)

  m <-
    prjs_df2050 %>% 
    filter(variable != "f") %>% 
    group_by(year, .id, variable) %>% 
    summarise(median =median(data, na.rm=TRUE)) 
  
  # cumulative catch
  p <-
    ggplot() +
    theme_bw() +
    geom_histogram(data= subset(prjs_df2050, variable == "rec"),
                   aes(x=data, fill=.id), 
                   position=position_dodge(), binwidth = filter(r, variable=="rec")$binwidth, alpha=0.5) +
    geom_histogram(data= subset(prjs_df2050, variable == "catch"),
                   aes(x=data, fill=.id), 
                   position=position_dodge(), binwidth = filter(r, variable=="catch")$binwidth, alpha=0.5) +
    geom_histogram(data= subset(prjs_df2050, variable == "ssb"),
                   aes(x=data, fill=.id), 
                   position=position_dodge(), binwidth = filter(r, variable=="ssb")$binwidth, alpha=0.5) +
    geom_histogram(data= subset(prjs_df2050, variable == "biomass"),
                   aes(x=data, fill=.id), 
                   position=position_dodge(), binwidth = filter(r, variable=="biomass")$binwidth, alpha=0.5) +
    # geom_freqpoly(aes(x=cumcatch, colour=.id), stat="identity") +
    geom_vline(data=m, aes(xintercept=median, colour=.id), show.legend = FALSE) +
    labs(y="", x="", colour="Fmsy", fill="scen") +
    facet_grid(.id ~ variable, scales="free")
  
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "histograms-2050.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  print(p) 
  dev.off()

#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
#@#  Equilibrium Check                                                       @#@
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  
  eqCurves=ldply(oms, function(x) { 
    model.frame(FLQuants(x[,dim(oms[[1]])[2]],
                         biomass=function(x) stock(x),
                         ssb    =function(x) ssb(  x),
                         catch  =function(x) catch(x),
                         f      =function(x) fbar( x)),drop=TRUE)[,-1]})
  
  ts=ldply(prj, function(x) model.frame(FLQuants(x, ssb  =function(x) ssb(  x), 
                                                 f    =function(x) fbar(x), 
                                                 catch=function(x) catch(x)),drop=T))
  ts=transform(ts,What=.id)
  ts=transform(ts,.id=unique(eqCurves$.id)[an(iter)])
  
  p <-
    ggplot()+
    geom_line( aes(ssb,catch,col=.id),data=eqCurves)+
    geom_point(aes(ssb,catch,col=.id),data=subset(ts,year==2050))+
    theme_bw()+theme(legend.position="bottom")+
    xlab("SSB")+ylab("Yield")
  
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "equilibrium_check.jpg", sep="_")),
       width=10, height=6, units="in", res=300)
  print(p) 
  dev.off()
  
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  
  
  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
  
  eqCurves=ldply(oms, function(x) { 
    model.frame(FLQuants(x[,dim(oms[[1]])[2]],
                         biomass=function(x) stock(x),
                         ssb    =function(x) ssb(  x),
                         catch  =function(x) catch(x),
                         f      =function(x) fbar( x)),drop=TRUE)[,-1]})
  
  fmsy=transform(ddply(eqCurves,.(.id), with, 
                       data.frame(Fmsy=f[catch==max(catch)])),lower=Fmsy*0.8,upper=Fmsy*1.2)
  msy=ddply(merge(eqCurves,fmsy,by=".id"),.(.id), with, 
            data.frame(lcatch=catch[min((lower-f)^2)==(lower-f)^2],
                       msy   =catch[catch==max(catch)],
                       ucatch=catch[min((upper-f)^2)==(upper-f)^2]))
  bmsy=ddply(merge(eqCurves,fmsy,by=".id"),.(.id), with, 
             data.frame(lssb   =ssb[min((lower-f)^2)==(lower-f)^2],
                        bmsy   =ssb[catch==max(catch)],
                        ussb   =ssb[min((upper-f)^2)==(upper-f)^2]))
  
  sch=ldply(prjs, function(x) ldply(x, function(x) { 
    model.frame(FLQuants(x,
                         biomass=function(x) stock(x),
                         ssb    =function(x) ssb(  x),
                         catch  =function(x) catch(x),
                         f      =function(x) fbar( x)),drop=TRUE)}))
  sch=transform(sch,.id=factor(.id,labels=c("Base","DD Mass","DD Mass, Mat","DD Mass, Mat, M")))
  
  sch=subset(sch,year==2050)
  sch=merge(sch,bmsy,by=".id")
  sch=merge(sch,msy,by=".id") 
  sch=merge(sch,fmsy,by=".id") 
  
  # Base
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "kobe base.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  kobe:::kobePhaseMar2(subset(transmute(subset(sch,.id=="Base"),
                                        stock=ssb/bmsy,
                                        harvest=catch/msy,
                                        run=ac(f))),
                       xlab=expression(B/B[MSY]),
                       ylab=expression(Catch/MSY),
                       col=c("red","grey","grey","grey")) 
  dev.off()

  # M
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "kobe M.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  kobe:::kobePhaseMar2(subset(transmute(subset(sch,.id=="DD Mass"),
                                        stock=ssb/bmsy,
                                        harvest=catch/msy,
                                        run=ac(f))),
                       xlab=expression(B/B[MSY]),ylab=expression(Catch/MSY),col=c("grey","red","grey","grey")) 
  dev.off()

  # MM
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "kobe MM.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  kobe:::kobePhaseMar2(subset(transmute(subset(sch,.id=="DD Mass, Mat"),stock=ssb/bmsy,harvest=catch/msy,run=ac(f))),
                       xlab=expression(B/B[MSY]),ylab=expression(Catch/MSY),col=c("grey","grey","red","grey")) 
  dev.off()

  # MMM
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "kobe MMM.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  kobe:::kobePhaseMar2(subset(transmute(subset(sch,.id=="DD Mass, Mat, M"),stock=ssb/bmsy,harvest=catch/msy,run=ac(f))),
                       xlab=expression(B/B[MSY]),ylab=expression(Catch/MSY),col=c("grey","grey","grey","red")) 
  dev.off()
  
  ## Yield & F #################################################################
  kobePhase=kobe:::kobePhase

  # Base
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "kobe FY base.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
       kobePhaseMar4(subset(transmute(subset(sch,.id=="Base"),
                                      stock  =f/Fmsy,
                                      harvest=catch/msy,
                                      run    =ac(f))),
                quadcol=c("yellow","yellow","green","red"),
                xlab=expression(F/F[MSY]),ylab=expression(Catch/MSY),col=c("red","grey","grey","grey"),
                xlim=3) 
  dev.off()
  
  # M
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "kobe FY M.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
       kobePhaseMar4(subset(transmute(subset(sch,.id=="DD Mass, Mat"),
                                      stock  =f/Fmsy,
                                      harvest=catch/msy,
                                      run    =ac(f))),
                quadcol=c("yellow","yellow","green","red"),
                xlab=expression(F/F[MSY]),ylab=expression(Catch/MSY),col=c("grey","grey","red","grey"),
                xlim=3) 
  dev.off()
  
  # MM
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "kobe FY MM.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
       kobePhaseMar4(subset(transmute(subset(sch,.id=="DD Mass, Mat, M"),
                                 stock  =f/Fmsy,
                                 harvest=catch/msy,
                                 run    =ac(f))),
                quadcol=c("yellow","yellow","green","red"),
                xlab=expression(F/F[MSY]),ylab=expression(Catch/MSY),col=c("grey","grey","grey","red"),
                xlim=3) 
  dev.off()
  
  # MMM
  jpeg(filename=file.path(figuresdir, paste(section,mystk, "kobe FY MMM.jpg", sep="_")),
       width=10, height=10, units="in", res=300)
  dev.off()
  
  # ggdensity(transmute(sch,stock=ssb/bmsy,harvest=catch/msy,run=ac(signif(f,3)),.id=.id),x="stock",fill=".id")+
  #   geom_vline(aes(xintercept=1),col="red")+
  #   facet_grid(run~.) 

  
  # ggdensity(transmute(sch,stock=ssb/bmsy,catch=catch/msy,run=ac(signif(f,3)),.id=.id),x="catch",fill=".id")+
  #   geom_vline(aes(xintercept=1),col="red")+
  #   facet_grid(run~.) 

  
  rm(p, p1, p2, p3, p4)
  save(list=ls(),
       file = file.path(dropboxdir, "results", mystk, paste(mystk,"section5.RData", sep="_")))
  load(file = file.path(dropboxdir, "results", mystk, paste(mystk,"section5.RData", sep="_")))
