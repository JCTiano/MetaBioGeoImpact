## Functions that are used in the meta-analysis --> Clean version
  
# 0. Average per study

avgstud <- function(df){
  
  dfn <- NULL
  
  for(i in (unique(df$respvar))){
    
    subset <- df[df$respvar == i,]
    
    for(s in (unique(subset$article_type_id))){
      
      subset2  <- subset[subset$article_type_id == s, ]
      newn     <- (subset2$n1 + subset2$n2) / 2
      avglnrr  <- sum(subset2$lnRR*newn)/sum(newn)
      avgvarln <- sum(subset2$ogVarLnRR*newn)/sum(newn)
      
      dfn <- rbind(dfn, data.frame(i, s, avglnrr, avgvarln, sum(newn)))
    }
    
  }
  return(dfn)
  
}


avgstud2 <- function(df){
  
  dfn <- NULL
  
  for(i in (unique(df$respvar))){
    
    subset <- df[df$respvar == i,]
    
    for(us in (unique(subset$article_type_id))){
      
      subset2  <- subset[subset$article_type_id == us, ]
      
      if(nrow(subset2) == 1){
        
        v <- i
        M <- subset2$lnRR
        Vm <- subset2$ogVarLnRR
        SD <- Vm^2
        SE <- sqrt(Vm)
        LL <- round(M - 1.96 * SE, digits = 2)
        UL <- round(M + 1.96 * SE, digits = 2)
        Z  <- M/SE
        pval <- 2*pnorm(q=Z, lower.tail=T)  
        n    <- length(subset2$respvar)
        s    <- length(unique(subset2$article_id))
        dfn <- rbind(dfn, data.frame(v, M, Vm, SD, SE, LL, UL, Z, pval, n, s, us))
        
      }else{
        DF <- sum(complete.cases(subset2$lnRR), na.rm = T) - 1 # degrees of freedom for observations that can be used
        
        # Calculations for log response ratios (lnRR)
        Q  <- sum(subset2$WY2, na.rm = T) - (sum(subset2$WY, na.rm = T))^2/sum(subset2$W, na.rm = T)
        C  <- sum(subset2$W, na.rm = T) - sum(subset2$W^2, na.rm = T)/sum(subset2$W, na.rm = T)
        T2 <- pmax(0, (Q - DF)/C)         # tau squared (between studies variance) never negative
        subset2$WRand  <- 1/(subset2$ogVarLnRR + T2)        # Weight = 1/variance within + variance between studies
        subset2$WYRand <- subset2$WRand * subset2$lnRR
        
        v <- i
        M  <- round((sum(subset2$WYRand, na.rm = T)/sum(subset2$WRand, na.rm = T)), digits = 2) # summary EF calculation
        Vm <- 1/(sum(subset2$WRand, na.rm = T))     # variance
        SD <- sd(subset2$lnRR, na.rm = T)
        SE <- sqrt(Vm)
        LL <- round(M - 1.96 * SE, digits = 2)
        UL <- round(M + 1.96 * SE, digits = 2)
        Z  <- M/SE
        pval <- 2*pnorm(q=Z, lower.tail=T)     # p value for 2 tailed test
        n    <- length(subset2$respvar)
        s    <- length(unique(subset2$article_id))
        
        dfn <- rbind(dfn, data.frame(v, M, Vm, SD, SE, LL, UL, Z, pval, n, s, us))
      }
      
    }
    
  }
  rem <- which(dfn$Vm == 0)
  dfn <- dfn[-rem,]
  return(dfn)
  
}

# 1. Summarizer

# Function to calculate average and tailed distributions for random effects assumption.
# Calculates the mean and variance around true response ratio assuming a random effects model.
# Average effect size +- confidence interval is calculated per variable --> Not taking into account slice depths.
# JT changed 'iv' to 'i' (typo?)

summarizer <- function(df) {
  
  dfn <- NULL
  
  for(i in (unique(df$respvar))){
    
    tt <- subset(df, df$respvar == i)
    DF <- sum(complete.cases(tt$lnRR), na.rm = T) - 1 # degrees of freedom for observations that can be used
    
    # Calculations for log response ratios (lnRR)
    Q  <- sum(tt$WY2, na.rm = T) - (sum(tt$WY, na.rm = T))^2/sum(tt$W, na.rm = T)
    C  <- sum(tt$W, na.rm = T) - sum(tt$W^2, na.rm = T)/sum(tt$W, na.rm = T)
    T2 <- pmax(0, (Q - DF)/C)         # tau squared (between studies variance) never negative
    tt$WRand  <- 1/(tt$VarLnRR + T2)        # Weight = 1/variance within + variance between studies
    tt$WYRand <- tt$WRand * tt$lnRR
    
    v <- i
    M  <- round((sum(tt$WYRand, na.rm = T)/sum(tt$WRand, na.rm = T)), digits = 2) # summary EF calculation
    Vm <- 1/(sum(tt$WRand, na.rm = T))     # variance
    SD <- sd(tt$lnRR, na.rm = T)
    SE <- sqrt(Vm)
    LL <- round(M - 1.96 * SE, digits = 2)
    UL <- round(M + 1.96 * SE, digits = 2)
    Z  <- M/SE
    pval <- 2*pnorm(q=Z, lower.tail=T)     # p value for 2 tailed test
    n    <- length(tt$respvar)
    s    <- length(unique(tt$article_id))
    
    dfn <- rbind(dfn, data.frame(v, M, Vm, SD, SE, LL, UL, Z, pval, n, s))
    
  }
  return(dfn)
}

# 1b_review. Summarizer

# Function to calculate average and tailed distributions for random effects assumption.
# Calculates the mean and variance around true response ratio assuming a random effects model.
# Average effect size +- confidence interval is calculated per variable --> Not taking into account slice depths.
# JT changed 'iv' to 'i' (typo?)
# 
# df <- dfshort
# i <- unique(df$respvar)[1]

summarizerb <- function(df) {
  
  dfn <- NULL
  
  for(i in (unique(df$respvar))){

    print(i)
    tt <- subset(df, df$respvar == i)
    sel <- which(tt$VarLnRR != 0)
    tt <- tt[sel,]
    
    if(length(tt$respvar) < 3){
      print(paste0("for ", i, " not enough rows."))
    } else{
      
    tt <- escalc(data = tt, yi = lnRR, vi = VarLnRR) # Just transforming in escalc object
    varcov <- vcalc(vi = vi, cluster = article_type_id, data = tt, nearpd = TRUE) # cluster variances in variance-covariance matrix per study type for use in aggregate.
    tt <- aggregate(x = tt, cluster = article_type_id, V = varcov, struct = "CS", weighted = FALSE, addk = TRUE) # aggregate effect sizes in single study effect size.
      
    out <- rma(data = tt, yi = yi, vi = vi)
    
    M <- round(out$b, 2)
    SD <- out$se/sqrt(length(tt$respvar))
    Vm <- SD^2
    SE <- out$se
    LL <- round(out$ci.lb, 2)
    UL <- round(out$ci.ub, 2)
    Z <- out$zval
    pval <- round(out$pval, 4)
    n    <- length(tt$respvar)
    s    <- length(unique(tt$article_id))
    v <- i
    
    dfn <- rbind(dfn, data.frame(v, M, Vm, SD, SE, LL, UL, Z, pval, n, s))
    }
  }
  return(dfn)
}


# JT !!! corrected the weighting. Was erronously weighted by 1 / the within study weight( this was the error) + between study variance ! 
# it is now corrected to have the random weights equal to 1 / within + between study VARIANCE

summarizerJT <- function(df) {
  
  dfn <- NULL
  
  for(i in (unique(df$respvar))){
    
    tt <- subset(df, df$respvar == i)
    DF <- sum(complete.cases(tt$lnRR), na.rm = T) - 1 # degrees of freedom for observations that can be used
    
    # Calculations for log response ratios (lnRR) ## JT ## modified to use metafor method for T2 calculation
    #ll = rma(yi = tt$lnRR, vi = tt$VarLnRR)         # Use rma function to get tau^2 and other values
    #T2 = ll$tau2
    Q  <- sum(tt$WY2, na.rm = T) - (sum(tt$WY, na.rm = T))^2/sum(tt$W, na.rm = T)
    C  <- sum(tt$W, na.rm = T) - sum(tt$W^2, na.rm = T)/sum(tt$W, na.rm = T)
    T2 <- pmax(0, (Q - DF)/C)         # tau squared (between studies variance) never negative
    tt$WRand  <- 1/(tt$VarLnRR + T2)  # Weight = 1/variance within + variance between studies
    tt$WYRand <- tt$WRand * tt$lnRR
    
    v <- i
    M  <- round((sum(tt$WYRand, na.rm = T)/sum(tt$WRand, na.rm = T)), digits = 2) # summary EF calculation
    Vm <- 1/(sum(tt$WRand, na.rm = T))     # variance
    SD <- sd(tt$lnRR, na.rm = T)
    SE <- sqrt(Vm)
    LL <- round(M - 1.96 * SE, digits = 2)
    UL <- round(M + 1.96 * SE, digits = 2)
    Z  <- M/SE
    pval <- 2*pnorm(q=Z, lower.tail=T)     # p value for 2 tailed test
    n    <- length(tt$respvar)
    s    <- length(unique(tt$article_id))
    
    dfn <- rbind(dfn, data.frame(v, M, Vm, SD, SE, LL, UL, Z, pval, n, s))
    
  }
  return(dfn)
}

# 2. Summarizer2


# Function to calculate average and tailed distributions for random effects assumption.
## With study types separated

summarizer2 <- function(df) {
  
  dfn <- NULL
  
  for(i in (unique(df$respvar))){
    
    tt <- subset(df, df$respvar == i)
    
    for(j in (unique(tt$hstype))){
      
      ts <- subset(tt, tt$hstype == j)
      
      DF <- sum(complete.cases(ts$lnRR), na.rm = T) - 1 # degrees of freedom for observations that can be used
      
      # Calculations for log response ratios (lnRR)
      Q  <- sum(ts$WY2, na.rm = T) - (sum(ts$WY, na.rm = T))^2/sum(ts$W, na.rm = T)
      C  <- sum(ts$W, na.rm = T) - sum(ts$W^2, na.rm = T)/sum(ts$W, na.rm = T)
      T2 <- pmax(0, (Q - DF)/C)       # tau squared (between studies variance) never negative
      ts$WRand  <- 1/(ts$VarLnRR + T2)        # Weight = 1/variance within + variance between studies
      ts$WYRand <- ts$WRand * ts$lnRR
      
      v <- i
      M  <- round((sum(ts$WYRand, na.rm = T)/sum(ts$WRand, na.rm = T)), digits = 2) # summary EF calculation
      Vm <- 1/(sum(ts$WRand, na.rm = T))                              # variance
      SD <- sd(ts$lnRR, na.rm = T)
      SE <- sqrt(Vm)
      LL <- round(M - 1.96 * SE, digits = 2)
      UL <- round(M + 1.96 * SE, digits = 2)
      Z  <- M/SE
      pval <- 2*pnorm(q=Z, lower.tail=T) # p value for 2 tailed test
      n    <- length(ts$respvar)
      s    <- length(unique(ts$article_id))
      stype <- j
      
      dfn <- rbind(dfn, data.frame(v, M, Vm, SD, SE, LL, UL, Z, pval, n, s, stype))
      
    }
    
  }
  return(dfn)
}

# 2b. Summarizer2b


# Function to calculate average and tailed distributions for random effects assumption.
## With study types separated

summarizer2b <- function(df) {
  
  dfn <- NULL
  
  for(i in (unique(df$respvar))){
    
    tt <- subset(df, df$respvar == i)
    
    for(j in (unique(tt$hstype))){
      
      ts <- subset(tt, tt$hstype == j)
      
      if(length(ts$respvar) < 2){
        
        print(paste0("for ", i, " not enough rows."))
        
      } else{
        
        ts <- escalc(data = ts, yi = lnRR, vi = VarLnRR) # Just transforming in escalc object
        varcov <- vcalc(vi = vi, cluster = article_type_id, data = ts, nearpd = TRUE) # cluster variances in variance-covariance matrix per study type for use in aggregate.
        ts <- aggregate(x = ts, cluster = article_type_id, V = varcov, struct = "CS", weighted = FALSE) # aggregate effect sizes in single study effect size.
        
        out <- rma(data = ts, yi = yi, vi = vi)
        
        M <- round(out$b, 2)
        SD <- out$se/sqrt(length(tt$respvar))
        Vm <- SD^2
        SE <- out$se
        LL <- round(out$ci.lb, 2)
        UL <- round(out$ci.ub, 2)
        Z <- out$zval
        pval <- round(out$pval, 4)
        n    <- length(ts$respvar)
        s    <- length(unique(ts$article_id))
        v <- i
        stype <- j
      
        dfn <- rbind(dfn, data.frame(v, M, Vm, SD, SE, LL, UL, Z, pval, n, s, stype))
      }
    }
  }
  return(dfn)
}

# 3. Summarizer - varslice

summarizer_varslice <- function(df, varz) {
  
  dfn <- NA
  
  for(i in (unique(varz))){
    
    tt <- subset(df, df$respvar == i)
    
    for(j in (unique(tt$depths))){
      
      dd <- subset(tt, tt$depths == j)
      dd <- escalc(data = dd, yi = lnRR, vi = VarLnRR) # Just transforming in escalc object
      varcov <- vcalc(vi = vi, cluster = article_type_id, data = dd, nearpd = TRUE) # cluster variances in variance-covariance matrix per study type for use in aggregate.
      avg <- aggregate(x = dd, cluster = article_type_id, V = varcov, struct = "CS", weighted = FALSE) # aggregate effect sizes in single study effect size.
      
      dfn <- rbind(dfn, avg)
      
    }
  }
  return(dfn)
}

# 4. Plotdepthslices
# This function plots the study-averaged response ratios per depth slice, calculated with summarizer_varslice, and adds one central mean which is calculated with a model of the form:
# rma.uni(yi, vi, data = out, mods = ~ depths - 1, subset = (respvar == "TOC"), method = "REML")
# as the true means.

plotdepthslices <- function(var, model, avgdf, maintitle, miny, maxy, textpos){
  
  df <- with(model, cbind(b, se, zval, pval, ci.lb, ci.ub))
  df <- data.frame(cbind(df, row.names(df)))
  df[,7] <- gsub("depths", "", df[,7])
  colnames(df) <- c("avg", "se", "zval", "pval", "ci.lb", "ci.ub", "depths")
  df$lab <- paste0(round(as.numeric(df$avg),2), " [", round(as.numeric(df$ci.lb),2), " - ", round(as.numeric(df$ci.ub),2), "]")
  
  ggplot(data = subset(avgdf, avgdf$respvar == var), aes(x = depths)) +
    coord_flip() +
    theme_classic() +  
    geom_pointrange(aes(y = yi, ymin = yi-vi, ymax = yi+vi), position = position_dodge2(width = 0.8), col = "azure3", size = 0.4, shape = 18) +
    geom_pointrange(data = df, aes(y = as.numeric(avg), ymin = as.numeric(ci.lb), ymax = as.numeric(ci.ub)), size = 0.8, col = "navyblue", shape = 18) + 
    geom_text(data = df, aes(x = depths, y = textpos, label = lab), size = 3) + 
    geom_hline(yintercept = as.numeric(0), linetype="dashed") +
    # labs(title = maintitle, subtitle = "est. [95% CI]", x = "Depth", y = "Log response ratio") + 
    labs(title = maintitle, x = "Depth", y = "Log response ratio") + 
    scale_x_discrete(limits = rev(levels(avgdf$depths))) + 
    theme_classic() +
    ylim(miny, maxy) + 
    theme(plot.title = element_text(hjust = 0.5),
          #      plot.subtitle = element_text(hjust = textpos),
          plot.margin = margin(r = 1),
          axis.title.x = element_text(vjust = 1.5))
}


# 5. Contvarplot
# This function returns the fitted models, for each of the depth strata.
# The fitted models are of the form: 
# rma.mv(yi = lnRR, V = VarLnRR, mods = ~depths*varname, random = ~ 1|article_type_id, data = subsetperstratum, method = "REML")

ct <- function(x){return(length(which(x > 0)))}

contvarplot <- function(resp, depvar, df, xlabz, log = FALSE) {
  
  # Subset of our variable
  sub <- subset(df, respvar == resp)
  t <- table(sub$article_type_id, sub$depths)
  
  # filters out when too few individual studies.
  lengths <- apply(t, 2, ct)
  lengths <- lengths[lengths >= 3] 
  sub2    <- subset(sub, depths %in% names(lengths))
  
  # make model with interactions
  if(log == TRUE){
    name <- log(sub2[,depvar]+1)
  }else{
    name <- sub2[,depvar]    
  }
  
  res <- rma.mv(yi = lnRR, V = VarLnRR, mods = ~depths*name, random = ~ 1|article_type_id, data = sub2, method = "REML")
  
  outdf <- NULL
  
  for(i in 1:length(names(lengths))){
    
    sub3  <- subset(sub2, depths == names(lengths)[i])
    
    if(log == TRUE){
      name2 <- log(sub3[,depvar]+1)
    }else{
      name2 <- sub3[,depvar]    
    }
    
    mod   <- rma.mv(yi = lnRR, V = VarLnRR, mods = ~name2, random = ~ 1|article_type_id, data = sub3, method = "REML", control=list(rel.tol=1e-8))
    # 
    #     pint  <- mod$pval[1]
    #     pcont <- mod$pval[2]
    #     coeff <- mod$b[2]
    # 
    #     regplot(x = mod,
    #             shade = mycol2, bg = mycol,
    #             xlab = xlabz, ylab = "lnRR", main = names(lengths)[i],
    #             cex.axis = 0.7)
    #     abline(h = 0, col = "darkgrey")
    # 
    #     mtext(side = 3, text = paste0("Sig. intercept = ", round(pint, 4)), line = -1.5, adj = 0.05, cex = 0.7)
    #     mtext(side = 3, text = paste0("Sig. var. = ", round(pcont, 4)), line = -2.5, adj = 0.05, cex = 0.7)
    #     mtext(side = 3, text = paste0("Coeff. var. = ", round(coeff, 2)), line = -3.5, adj = 0.05, cex = 0.7)
    
    newmodinf <- data.frame("Var" = resp,
                            "Covar" = xlabz,
                            "Strata" = names(lengths)[i],
                            "QE" = mod$QE,
                            "QM" = mod$QM,
                            "QEp" = mod$QEp,
                            "QMp" = mod$QMp,
                            "int_est" = mod$b[1],
                            "int_sig" = mod$pval[1],
                            "coeff_est" = mod$b[2],
                            "coeff_sig" = mod$pval[2]
    )
    
    outdf <- rbind(outdf, newmodinf)
    
    
  }
  
  return(outdf)
  
}

ct <- function(x){return(length(which(x > 0)))}

plotreg <- function(resp, depvar, df, xlabz, log = FALSE, sample="Surface", title = "") {
  
  # Subset of our variable
  sub <- subset(df, respvar == resp)
  t <- table(sub$article_type_id, sub$depths)
  
  # filters out when too few individual studies.
  lengths <- apply(t, 2, ct)
  lengths <- lengths[lengths >= 3] 
  sub2    <- subset(sub, depths %in% names(lengths))
  
  # make model with interactions
  if(log == TRUE){
    name <- log(sub2[,depvar]+1)
  }else{
    name <- sub2[,depvar]    
  }
  
  res <- rma.mv(yi = lnRR, V = VarLnRR, mods = ~depths*name, random = ~ 1|article_type_id, data = sub2, method = "REML")
  
  outdf <- NULL
  
  #for(i in 1:length(names(lengths))){
    
    sub3  <- subset(sub2, depths == sample)
    
    if(log == TRUE){
      name2 <- log(sub3[,depvar]+1)
    }else{
      name2 <- sub3[,depvar]    
    }
    
    mod   <- rma.mv(yi = lnRR, V = VarLnRR, mods = ~name2, random = ~ 1|article_type_id, data = sub3, method = "REML", control=list(rel.tol=1e-8))
    # 
    pint  <- mod$pval[1]
    pcont <- mod$pval[2]
    coeff <- mod$b[2]
    # 
    Col = adjustcolor("darkgreen", 0.5)
    Col2 = adjustcolor("darkgreen", 0.5)
    
    #regplot(x = mod,
    #        shade = mycol2, bg = mycol,
    #        xlab = xlabz, ylab = "lnRR", main = names(lengths)[i],
    #        cex.axis = 0.7, las = 1)
    regplot(x = mod, bg = Col, col = Col2, shade = adjustcolor("black", 0.3), grid = T, refline = 0, 
            ylab = "ln(RR)", main = title, log = "x",
            xlab = xlabz, las = 1, lcol = c("white","grey25","grey1","grey50"), 
            lty = c(1,0,1,1), lwd = c(2,1,1,2), cex.axis = 0.6, cex.lab = 0.75, cex.main = 0.8)
    abline(h = 0, col = "darkgrey")
    # 
    # mtext(side = 3, text = paste0("Sig. intercept = ", round(pint, 4)), line = -1.5, adj = 0.05, cex = 0.7)
    # mtext(side = 3, text = paste0("Sig. var. = ", round(pcont, 4)), line = -2.5, adj = 0.05, cex = 0.7)
    # mtext(side = 3, text = paste0("Coeff. var. = ", round(coeff, 2)), line = -3.5, adj = 0.05, cex = 0.7)
    # 
    # newmodinf <- data.frame("Var" = resp,
    #                         "Covar" = xlabz,
    #                         "Strata" = names(lengths)[i],
    #                         "QE" = mod$QE,
    #                         "QM" = mod$QM,
    #                         "QEp" = mod$QEp,
    #                         "QMp" = mod$QMp,
    #                         "int_est" = mod$b[1],
    #                         "int_sig" = mod$pval[1],
    #                         "coeff_est" = mod$b[2],
    #                         "coeff_sig" = mod$pval[2]
    # )
    # 
    # outdf <- rbind(outdf, newmodinf)
    
    
  }
  
  #return(outdf)
  
#}


plotreg2 <- function(resp, depvar, df, xlabz, logvar = FALSE, sample="Surface", title = "") {
  
  # Subset of our variable
  sub <- subset(df, respvar == resp)
  t <- table(sub$article_type_id, sub$depths)
  
  # filters out when too few individual studies.
  lengths <- apply(t, 2, ct)
  lengths <- lengths[lengths >= 3] 
  sub2    <- subset(sub, depths %in% names(lengths))
  
  # make model with interactions
  if(logvar == TRUE){
    name <- log(sub2[,depvar]+1)
  }else{
    name <- sub2[,depvar]+0.0001    
  }
  
  res <- rma.mv(yi = lnRR, V = VarLnRR, mods = ~depths*name, random = ~ 1|article_type_id, data = sub2, method = "REML")
  
  outdf <- NULL
  
  #for(i in 1:length(names(lengths))){
  
  sub3  <- subset(sub2, depths == sample)
  
  if(logvar == TRUE){
    name2 <- log(sub3[,depvar]+1)
  }else{
    name2 <- sub3[,depvar]+0.0001    
  }
  
  mod   <- rma.mv(yi = lnRR, V = VarLnRR, mods = ~name2, random = ~ 1|article_type_id, data = sub3, method = "REML", control=list(rel.tol=1e-8))
  # 
  pint  <- mod$pval[1]
  pcont <- mod$pval[2]
  coeff <- mod$b[2]
  # 
  Col = adjustcolor("darkgreen", 0.5)
  Col2 = adjustcolor("darkgreen", 0.5)
  
  #regplot(x = mod,
  #        shade = mycol2, bg = mycol,
  #        xlab = xlabz, ylab = "lnRR", main = names(lengths)[i],
  #        cex.axis = 0.7, las = 1)
  regplot(x = mod, bg = Col, col = Col2, shade = adjustcolor("black", 0.3), grid = T, refline = 0, 
          ylab = "ln(RR)", main = title, log = "x",
          xlab = xlabz, las = 1, lcol = c("white","grey25","grey1","grey50"), 
          lty = c(1,0,1,1), lwd = c(2,1,1,2), cex.axis = 0.6, cex.lab = 0.75, cex.main = 0.8)
  abline(h = 0, col = "darkgrey")
  # 
  # mtext(side = 3, text = paste0("Sig. intercept = ", round(pint, 4)), line = -1.5, adj = 0.05, cex = 0.7)
  # mtext(side = 3, text = paste0("Sig. var. = ", round(pcont, 4)), line = -2.5, adj = 0.05, cex = 0.7)
  # mtext(side = 3, text = paste0("Coeff. var. = ", round(coeff, 2)), line = -3.5, adj = 0.05, cex = 0.7)
  # 
  # newmodinf <- data.frame("Var" = resp,
  #                         "Covar" = xlabz,
  #                         "Strata" = names(lengths)[i],
  #                         "QE" = mod$QE,
  #                         "QM" = mod$QM,
  #                         "QEp" = mod$QEp,
  #                         "QMp" = mod$QMp,
  #                         "int_est" = mod$b[1],
  #                         "int_sig" = mod$pval[1],
  #                         "coeff_est" = mod$b[2],
  #                         "coeff_sig" = mod$pval[2]
  # )
  # 
  # outdf <- rbind(outdf, newmodinf)
  
  
}

#return(outdf)

#}

# 6. ModTab
# Extract info from means-per-slice models
# 11. Model fit - table
## Fits all models and returns a table

modTab <- function(model, var) {
  
  estims <- model[["beta"]]
  se <- model[["se"]]
  pval <- model[["pval"]]
  confintlb <- model[["ci.lb"]]
  confintub <- model[["ci.ub"]]
  names <- dimnames(model[["beta"]])[[1]]
  
  dfnew <- data.frame("Variable" = var                ,
                      "Strata"   = names              ,
                      "Estimate" = round(estims, 3)   ,
                      "SE"       = round(se, 3)       ,
                      "p-value"  = round(pval, 4)     ,
                      "CI-lb"    = round(confintlb, 3),
                      "CI-ub"    = round(confintub, 3))
  return(dfnew)
}

# 7. Plot HM
# Plot a heatmap of either the coefficients or the intercepts.

plothmcoeff <- function(model){
  
  mod <- model
  
  mod$Strata <- factor(mod$Strata, levels = c("Surface", "Subsurface", "Deep", "Very deep", "Full sample"))
  mod$Covar  <- factor(mod$Covar, levels = rev(c("log (timesincetrawl)", "trawling effort", "Current velocity", "Water depth", "Primary productivity", "habitat", "Season")))
  # Define the manual break levels
  break_levels <- c(-0.0001, 0.001001, 0.01001, 0.05001, 0.1, 1)
  
  # Discretize p-values into classes based on the manual breaks
  mod$pvalue_class <- cut(mod$coeff_sig, breaks = break_levels, labels = c("<0.001", "0.001-0.01", "0.01-0.05", "0.05-0.1", ">0.1"))
  
  # # Discretize p-values into classes based on the manual breaks
  # mod$pvalue_class <- cut(mod$int_sig, breaks = break_levels, labels = c("<0.001", "0.001-0.01", "0.01-0.05", "0.05-0.1", ">0.1"))
  
  # Define custom class names
  class_names <- c("<0.001", "0.001-0.01", "0.01-0.05", "0.05-0.1", ">0.1")
  
  mod$pvalue_class <- factor(mod$pvalue_class, levels = class_names)
  
  # Define custom colors
  custom_colors <- c("#993404", "#d95f0e","#fe9929", "#fed98e","#ffffd4")
  
  # Create a heatmap with discretized p-values using custom colors
  p <- ggplot(mod, aes(x = Strata, y = Covar, fill = pvalue_class, label = round(coeff_est, 2))) +
    geom_tile() +
    geom_text(color = "black", size = 4.5) +
    scale_fill_manual(
      values = custom_colors,
      labels = class_names,
      drop = FALSE
    ) +
    labs(
      title = paste(mod$Var),
      x = "",
      y = ""
    ) + 
    theme_minimal()
  p
  return(p)
}


plothmint <- function(model){
  
  mod <- model
  
  mod$Strata <- factor(mod$Strata, levels = c("Surface", "Subsurface", "Deep", "Very deep", "Full sample"))
  mod$Covar  <- factor(mod$Covar, levels = rev(c("log (timesincetrawl)", "trawling effort", "Current velocity", "Water depth", "Primary productivity", "habitat", "Season")))
  # Define the manual break levels
  break_levels <- c(-0.0001, 0.001001, 0.01001, 0.05001, 0.1, 1)
  
  # Discretize p-values into classes based on the manual breaks
  # mod$pvalue_class <- cut(mod$coeff_sig, breaks = break_levels, labels = c("<0.001", "0.001-0.01", "0.01-0.05", "0.05-0.1", ">0.1"))
  
  # # Discretize p-values into classes based on the manual breaks
  mod$pvalue_class <- cut(mod$int_sig, breaks = break_levels, labels = c("<0.001", "0.001-0.01", "0.01-0.05", "0.05-0.1", ">0.1"))
  
  # Define custom class names
  class_names <- c("<0.001", "0.001-0.01", "0.01-0.05", "0.05-0.1", ">0.1")
  
  mod$pvalue_class <- factor(mod$pvalue_class, levels = class_names)
  
  # Define custom colors
  custom_colors <- c("#993404", "#d95f0e","#fe9929", "#fed98e","#ffffd4")
  
  # Create a heatmap with discretized p-values using custom colors
  p <- ggplot(mod, aes(x = Strata, y = Covar, fill = pvalue_class, label = round(int_est, 2))) +
    geom_tile() +
    geom_text(color = "black", size = 4.5) +
    scale_fill_manual(
      values = custom_colors,
      labels = class_names,
      drop = FALSE
    ) +
    labs(
      title = paste(mod$Var),
      x = "",
      y = ""
    ) + 
    theme_minimal()
  p
  return(p)
}

# 7b. Plot HM with less covariates
# Plot a heatmap of either the coefficients or the intercepts.

plothmcoeff2 <- function(model, labz = TRUE){

mod <- model

mod$Strata <- factor(mod$Strata, levels = c("Surface", "Subsurface", "Deep", "Very deep", "Full sample"))
mod$Covar  <- factor(mod$Covar, levels = rev(c("log (timesincetrawl)", "Current velocity", "Water depth", "Primary productivity")))
# Define the manual break levels
break_levels <- c(-0.0001, 0.001001, 0.01001, 0.05001, 0.1, 1)

# Discretize p-values into classes based on the manual breaks
mod$pvalue_class <- cut(mod$coeff_sig, breaks = break_levels, labels = c("<0.001", "0.001-0.01", "0.01-0.05", "0.05-0.1", ">0.1"))

# # Discretize p-values into classes based on the manual breaks
# mod$pvalue_class <- cut(mod$int_sig, breaks = break_levels, labels = c("<0.001", "0.001-0.01", "0.01-0.05", "0.05-0.1", ">0.1"))

# Define custom class names
class_names <- c("<0.001", "0.001-0.01", "0.01-0.05", "0.05-0.1", ">0.1")

mod$pvalue_class <- factor(mod$pvalue_class, levels = class_names)

# Define custom colors
custom_colors <- c("#993404", "#d95f0e","#fe9929", "#fed98e","#ffffd4")

if (labz == TRUE) {

# Create a heatmap with discretized p-values using custom colors
p <- ggplot(mod, aes(x = Strata, y = Covar, fill = pvalue_class, label = round(coeff_est, 2))) +
  geom_tile() +
  geom_text(color = "black", size = 4.5) +
  scale_fill_manual(
    values = custom_colors,
    labels = class_names,
    drop = FALSE
  ) +
  labs(
    title = paste(mod$Var),
    x = "",
    y = ""
  ) + 
  theme_minimal()

} else {
  
p <- ggplot(mod, aes(x = Strata, y = Covar, fill = pvalue_class, label = round(coeff_est, 2))) +
    geom_tile() +
    geom_text(color = "black", size = 4.5) +
    scale_fill_manual(
      values = custom_colors,
      labels = class_names,
      drop = FALSE
    ) +
    labs(
      title = paste(mod$Var),
      x = "",
      y = ""
    ) + 
    theme_minimal() + 
    theme(axis.text.y = element_blank())

}
p
return(p)
}

