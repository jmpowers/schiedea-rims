plot.lmList <-
  function(object,ord.var=1,...) {
    ## Ord.var indicates which coefficient to 
    ##     sort the data frame of coefficients by.
    require(reshape)
    ## create a data.frame of coefs from list of glm fits.
    cL <- coef(object) 
    ## Adds a column for group ID --
    ## rownames (the variable upon which the fits were conditioned)
    ## -- here genotype. 
    cL$grp  <- rownames(cL) 
    if (is.numeric(ord.var) & length(ord.var)==1) {
      ## identify the ordering coefficient by name
      ord.var <- names(cL)[ord.var+1] 
      if (!ord.var %in% names(cL)) stop("unknown ordering variable")
      ## create a new order of factor levels for group ID
      cL$grp <- reorder(cL$grp,cL[[ord.var]])
    } else  
    ##otherwise just use the original ordering
    cL$grp <- reorder(cL$grp,ord.var) 
    ##"melt" or stack the data frame to 
    ##   allow a dotplot of each coefficient.
    dotplot(grp~value|variable,data=melt(cL),...) 
  }


qqmath.lmList <- function(object,...) {
  require(reshape)
  qqmath(~value|variable,data=melt(coef(object)),
         prepanel = prepanel.qqmathline,
         panel = function(x, ...) {
           panel.qqmathline(x, ...)
           panel.qqmath(x, ...)
         },
         scale=list(y=list(relation="free")))
}

overdisp <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
#sligthty different version that works
# If the p-value is < 0.05, the data are overdispersed from the model.
overdisp1 <- function(model) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp2 <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  (rdf <- nrow(model@frame)-model.df) #model@frame?
  rp <- residuals(model)              #missing pearson
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE,log.p=TRUE) #log.p=T
  c(chisq=Pearson.chisq,ratio=prat,p=exp(pval))
}

comp_od <- function(mod) {
#for handling overdispersion - see bolker GLMM FAQ
cc <- coef(summary(mod))
phi <- overdisp(mod)["ratio"]
cc <- within(as.data.frame(cc),
             {   `Std. Error` <- `Std. Error`*sqrt(phi)
             `z value` <- Estimate/`Std. Error`
             `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
             })
printCoefmat(cc,digits=3)
}

sprint <- function(m) lme4:::printMer(m,correlation=FALSE)

locscaleplot <- function(model,col="black", ...) {
  f <- fitted(model)
  r <- sqrt(abs(residuals(model)))  ## had better be Pearson resids
  plot(f,r,col=col, ...) 
  L1 <- loess(r~f)
  fvec = seq(min(f),max(f),length.out=150)
  lines(fvec,predict(L1,fvec),col=3, lwd=4)
}

dfun <- function(x) {
  x$AIC <- x$AIC-min(x$AIC)
  names(x)[names(x)=="AIC"] <- "dAIC"
  x
}

badcorr <- function(x) {
  cc <- attr(x,"correlation")
  diag(cc) <- 0
  any(abs((abs(cc)-1))<1e-5)
}
anybadcorr <- function(x) {
  any(sapply(VarCorr(x),badcorr))
}

locscaleplot <- function(model,col="black") {
  f <- fitted(model)
  r <- abs(residuals(model))
  plot(f,r,col=col) 
  L1 <- loess(r~f)
  fvec = seq(min(f),max(f),length.out=150)
  lines(fvec,predict(L1,fvec),col=2)
}

printvc <- function(m,digits=2,ctol=1e-3) {
  v <- VarCorr(m)
  prtfun <- function(x) {
    cc <- attr(x,"correlation")
    diag(cc) <- 0
    corstr <- ifelse(abs(abs(cc)-1)<ctol,"="," ")
    ss <- format(x,digits=digits) ## sprintf(fmt,x)
    ss <- paste(ss,corstr,sep="")
    m <- matrix(ss,nrow=nrow(x))
    m[upper.tri(m)] <- ""
    dimnames(m) <- dimnames(x)
    ## writeLines(apply(m,1,paste,collapse=" "))
    print(m,quote=FALSE)
  }
  for (i in seq_along(v)) {
    cat(names(v)[i],":\n",sep="")
    prtfun(v[[i]])
    cat("\n")
  }
}

plot.ICtab <- function(x,sort=TRUE,within) {
  z <- with(x,data.frame(n=attr(x,"row.names"),dAIC))
  if (sort) z <- transform(z,n=reorder(n,dAIC))
  if (!missing(within)) z$dAIC[z$dAIC>within] <- NA
  dotplot(n~dAIC,data=z)
}
