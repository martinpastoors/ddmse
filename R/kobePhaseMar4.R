

kobePhaseMar4=function(pts,trks=NULL,mns=FALSE,size=1,
                       xlab=expression(B:B[MSY]),
                       ylab=expression(F:F[MSY]),
                       xlim=2,ylim=xlim,
                       quadcol=c("red","green","yellow","yellow"),
                       col =colorRampPalette(c("orange","blue"),space="Lab"),
                       shade=.5,col2=grey(shade),col3=grey(shade*1.1),
                       layer=NULL,
                       bref=1,
                       fref=1){
  
  if (!("run" %in% names(pts)))
    pts=cbind(pts,run=factor(1))
  if (!is.null(trks) & !("run" %in% names(trks)))
    trks=cbind(trks,run=factor(1))
  
  if ("function" %in% is(col))
    col=col(length(unique(pts$run)))
  
  if (length(size)==1) size=rep(size,2)
  
  ##### Density plots   #############################################################################################
  
  # second density plot, oriented vertically (hence the 'coord_flip()' at the end
  dH<-ggplot(pts) + 
    #geom_density(aes(x = harvest, y =  ..count.., group=run), fill=col2, col=col3, position = "stack") + 
    geom_density(aes(x = harvest, y = ..count..,               fill=run, alpha=0.4)) + 
    geom_vline(xintercept=fref,col="red",data=data.frame(fref=fref))+
    coord_cartesian(xlim=c(0,ylim))   +
    scale_fill_manual(values=col)          +
    xlab(xlab) + ylab(ylab)                +
    theme(legend.position = "none", 
          axis.title.x = element_text(colour ='NA'), 
          axis.text.x  = element_text(colour ="NA"), 
          axis.ticks.x = element_line(colour ="NA"),
          axis.ticks =   element_line(colour ="NA"),
          
          axis.title.y = element_blank(), 
          axis.text.y  = element_blank(), 
          axis.ticks.y = element_blank(), 
          
          plot.margin = unit(c(0, 0, 1, 0), "lines"),
          panel.background = element_rect(fill   ="NA", colour ="NA"), 
          panel.border     = element_rect(fill   ="NA", colour ="NA"), 
          panel.grid.major = element_line(colour ="NA"), 
          panel.grid.minor = element_line(colour ="NA")                    
    )
  
  # kobe phase plot
  kC=kobePhase(pts,quadcol=quadcol,bref=bref,fref=fref) +
    geom_point(aes(stock,harvest,group=run),col="black",size=size[1]) +  
    geom_point(aes(stock,harvest,col=run,group=run),size=size[1]*.5) +  
    coord_cartesian(xlim=c(0,xlim),ylim=c(0,ylim)) +
    scale_colour_manual(values=col)      +
    xlab(xlab) + ylab(ylab)              +
    theme(legend.position = "none",
          axis.text.y  = element_text(colour="black", angle=90), 
          plot.margin = unit(c(0, 0, 1, 1), "lines"))
  
  if ("LayerInstance"%in%is(layer))
    kC=kC+layer
  else if ("list"%in%is(layer))
    kC=kC+layer[[1]]+layer[[2]]
  
  
  if (mns)
    kC=kC+geom_point(aes(stock,harvest,col=run,group=run),size=6.0*size[1], colour="black",  data=ddply(pts,.(run),function(x) data.frame(stock=median(x$stock),harvest=median(x$harvest)))) +
    geom_point(aes(stock,harvest,col=run,group=run),size=4.5*size[1], colour="cyan",   data=ddply(pts,.(run),function(x) data.frame(stock=median(x$stock),harvest=median(x$harvest))))
  if (!is.null(trks))
    kC=kC+geom_path(aes(stock,harvest, col=run,group=run),size=1*size[2], data=trks)   
  
  fnVP=function(dH,kC){
    vplayout <- function(x, y)
      viewport(layout.pos.row = x, layout.pos.col = y)
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(4, 5)))  # 5 by 5 grid
    print(dH +coord_flip(xlim=c(0,ylim)), vp=vplayout(1:4,5))         # 2nd to the left +opts(legend.position = c(0,1.05)) + opts(legend.text = theme_text(colour = "black")) 
    print(kC, vp=vplayout(1:4,1:4))                     # the main x/y plot will instead spread across most of the grid
  }
  
  fnVP(dH,kC)
  
  invisible(list(harvest=dH,phase=kC))}
