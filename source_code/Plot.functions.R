#-------------------------------------------------------------------------------
#        Summary effect size functions for plotting and data retrieval                        
#                                Justin C. Tiano                                 
#-------------------------------------------------------------------------------
#
# Various functions used to calculate and create summary effect size plots using in the 2022 meta-analysis 
# on trawling and biogeochemistry. The functions plot 2 types of summary effect sizes: 1) log response ratios
# (LnRR), and 2) standardized mean difference (Cohen's d). 'GrPlot' plots the standardized mean difference
# from converted Fishers Z effect sizes calculated for trawling gradient studies. "getdats" is used to retrieve 
# LnRR or Cohen's d effect sizes for a selected variable.
#
# CONTENTS:
# 1. function PlotSum
# 2. function GrPlot
# 3. function PlotsD
# 4. function PlotEO
# 5. function getdats
#-------------------------------------------------------------------------------

# 1. PlotSum
# Function for main summary effect size plot. Plots summary effect sizes for Cohen's d or LnRR.
# ### DEBUG ###
# df <- short
# parm <- "Chl-a"
# ### ENDDEBUG ###

PlotSum = function(df, parm = "Chl-a", start = F, row = 1, textpos = -1, 
                   xlow = -1, xhigh = 1, ylow = 0, yhigh = 2, d = FALSE, 
                   textpos2 = 2, subsetname = "") {
  
  tt = subset(df, df$respvar == parm )             # specify response variable 
  DF = sum(complete.cases(tt$lnRR), na.rm = T) - 1 # degrees of freedom for observations that can be used
  
  if (d == TRUE) {
    # Calculations for standardized mean difference (Cohen's d)
    ll = rma(yi = tt$d, vi = tt$Vd)                           # Use rma function to get tau^2 and other values
    T2 = ll$tau2
    tt$WYRand = tt$WRand * tt$d    # Weight * Effect size
    
    M  = (sum(tt$WYRand, na.rm = T)/sum(tt$WRand, na.rm = T))  # summary EF calculation
    Vm = 1/(sum(tt$WRand, na.rm = T))                          # variance
    SD = sd(tt$g, na.rm = T)                                   # Standard deviation of Cohen's d in subset
    SE = sqrt(Vm)
    LL = M - 1.96 * SE
    UL = M + 1.96 * SE
    
  } else {
    # Calculations for log response ratios (lnRR)
    Q = sum(tt$WY2, na.rm = T) - (sum(tt$WYsquared, na.rm = T))/sum(tt$W, na.rm = T)
    C = sum(tt$W, na.rm = T) - sum(tt$W^2, na.rm = T)/sum(tt$W, na.rm = T)
    T2  = pmax(0, (Q - DF)/C)       # tau squared (between studies variance) never negative
    tt$WRand = 1/(tt$W + T2)        # Weight = 1/variance within + variance between studies
    tt$WYRand = tt$WRand * tt$lnRR
    
    M  = round((sum(tt$WYRand, na.rm = T)/sum(tt$WRand, na.rm = T)), digits = 2) # summary EF calculation
    Vm = 1/(sum(tt$WRand, na.rm = T))                              # variance
    SD = sd(tt$lnRR, na.rm = T)                                    # Standard deviation of Cohen's d in subset
    SE = sqrt(Vm)
    LL = round(M - 1.96 * SE, digits = 2)
    UL = round(M + 1.96 * SE, digits = 2)
    Z  = M/SE
    pval = 2*pnorm(q=Z, lower.tail=T) # p value for 2 tailed test
  }
  if (start == TRUE) {
    plot(x = c(LL,UL), y = c(row,row), type = "l", xlim = c(xlow, xhigh), ylim = c(ylow, yhigh), yaxt = "n",
         main = "", xlab = "Summary ln(RR)", ylab = "", axes = F, cex.main = 0.8, cex.lab = 0.8)
    axis(1, cex.axis = 0.8)
    points(x = M, y = row, pch = 15, cex = 0.8)
    abline(h = c(row+0.5, 0), lty = 1)
    abline(v = 0, lty = 1, col = adjustcolor("grey40", 0.5), lwd = 2)
    text(textpos, row, labels = paste(parm,  "[",DF + 1,"]"), cex = 0.8)
    arrows(x0 = LL, x1 = UL, y0 = c(row,row), length = 0.03, angle = -90, code = 3)
    text(textpos2, row, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.8) # summary EF and 95% CI's
    mtext("Parameter [independent observations]", side = 3, at = xlow+2, font = 3, cex = 0.75)
    mtext("ln(RR) [95% CI]", side = 3, at = xhigh-2, font = 3, cex = 0.75)
    mtext("Log Response Ratio", side = 3, line = 1, font = 2, cex = 0.8)
    mtext(subsetname, side = 3, line = 0, font = 2, cex = 0.8)
    
  } else {
    lines(x = c(LL,UL), y = c(row,row))
    points(x = c(M), y =  c(row), pch = 15, cex = 0.8)
    text(textpos, row, labels = paste(parm,  "[",DF + 1,"]"), cex = 0.8)
    abline(h = row+0.5, lty = 3, col = adjustcolor("black", 0.3))
    arrows(x0 = LL, x1 = UL, y0 = c(row,row), length = 0.03, angle = -90, code = 3)
    text(textpos2, row, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.8) # summary EF and 95% CI's
  }
}

# 2. GrPlot
# Function for plot that compares gradient studies and comparative studies (input for gradient studies). Plots summary effect sizes Cohen's d (or LnRR) 

GrPlot = function(parm="TOC", start = T, row = 1, textpos = -1, xlow = -1, xhigh = 1, ylow = 0, yhigh = 2,
                  textpos2 = 10) {
  
  tt = subset(GradInfo, GradInfo$Response.variable==parm )   # specify response variable 
  DF = nrow(tt) - 1     # degrees of freedom for 'studies'(comparisons) that can be used
  
  # Calculations for standardized mean difference (Cohen's d)
  ll = rma(yi = tt$d, vi = tt$Vd)                           # Use rma function to get tau^2 and other values
  T2 = ll$tau2
  tt$WRand = 1/tt$W + T2      # Weight = 1/variance within + variance between studies
  tt$WYRand = tt$WRand * tt$d
  
  M  = round((sum(tt$WYRand, na.rm = T)/sum(tt$WRand, na.rm = T)), digits = 2) # summary EF calculation
  Vm = 1/(sum(tt$WRand, na.rm = T))                          # variance
  SD = sd(tt$g, na.rm = T)                                   # Standard deviation of Cohen's d in subset
  SE = sqrt(Vm)
  LL = round(M - 1.96 * SE, digits = 2)
  UL = round(M + 1.96 * SE, digits = 2)
  
  move = -0.2
  Col = adjustcolor("red", 0.7)
  if (start == TRUE) {
    plot(x = c(LL,UL), y = c(row+move,row+move), type = "l", xlim = c(xlow, xhigh), ylim = c(ylow, yhigh), yaxt = "n", 
         axes = F, cex.main = 0.8, cex.lab = 0.8,
         main = "", xlab = "Summary Cohen's d", ylab = "", col = Col)
    axis(1, cex.axis = 0.8)
    abline(h = c(row+0.65, -0.1), lty = 1)
    abline(v = 0, lty = 1, col = adjustcolor("grey40", 0.5), lwd = 2)
    points(x = M, y = row+move, pch = 15, col = Col, cex = 0.8)
    text(textpos-2, row, labels = paste("[",DF + 1,"]"), cex = 0.8, col = Col)

    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3, col = Col)
    mtext("Parameter[independent observations]", side = 3, at = xlow+8, font = 3, cex = 0.75)
    mtext("Standardized mean difference", side = 3, line = 1, font = 2, cex = 0.8)
    mtext("Cohen's d [95% CI]", side = 3, at = xhigh-5, font = 3, cex = 0.75)
    text(textpos2, row-0.2, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.8, col = Col) # summary EF and 95% CI's
  } else {
    lines(x = c(LL,UL), y = c(row+move,row+move), col = Col)
    points(x = c(M), y =  c(row+move), pch = 15, col = Col, cex = 0.8)
    text(textpos-2, row, labels = paste("[",DF + 1,"]"), cex = 0.8, col = Col)
    abline(h = row+0.5, lty = 3, col = adjustcolor("black", 0.3))
    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3, col = Col)
    text(textpos2, row-0.1, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.8, col = Col) # summary EF and 95% CI's
  }
}
# 3. PlotD
# Function for plot that compares gradient studies and comparative studies (input for comparative studies). Plots summary effect sizes Cohen's d

PlotD = function(df, parm = "Chl-a", start = T, row = 1, textpos = -1, xlow = -1, xhigh = 1, ylow = 0, yhigh = 2, d = FALSE, mv = 0.2,
                 textpos2 = 10) {
  tt = subset(df, df$respvar==parm )            # specify response variable 
  DF = sum(complete.cases(tt$lnRR), na.rm = T) - 1       # degrees of freedom for 'studies'(comparisons) that can be used
  move = mv

    # Calculations for standardized mean difference (Cohen's d)
    ll = rma(yi = tt$d, vi = tt$Vd) # Use rma function to get tau^2 and other values
    T2 = ll$tau2
    tt$WRand = 1/(tt$Wd + T2)       # Weight = 1/variance within + variance between studies
    tt$WYRand = tt$WRand * tt$d     # Weight * Effect size
    
    M  = round((sum(tt$WYRand, na.rm = T)/sum(tt$WRand, na.rm = T)), digits = 2) # summary EF calculation
    Vm = 1/(sum(tt$WRand, na.rm = T))                          # variance
    SD = sd(tt$g, na.rm = T)                                   # Standard deviation of Cohen's d in subset
    SE = sqrt(Vm)
    LL = round(M - 1.96 * SE, digits = 2)
    UL = round(M + 1.96 * SE, digits = 2)
    
  if (start == TRUE) {
    plot(x = c(LL,UL), y = c(row+move,row+move), type = "l", xlim = c(xlow, xhigh), ylim = c(ylow, yhigh), yaxt = "n",
         main = "Parameter [Independent Observations]", xlab = "Summary effect size", ylab = "")
    points(x = M, y = row+move, pch = 15, cex = 0.8)
    abline(v = 0, lty = 2)
    text(textpos, row, labels = paste(parm,  "[",DF + 1," ]"), cex = 0.75)
    abline(h = row+0.5, lty = 3, col = adjustcolor("black", 0.3))
    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3)
    text(textpos2, row+0.1, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.75) # summary EF and 95% CI's
    
  } else {
    lines(x = c(LL,UL), y = c(row+move,row+move))
    points(x = c(M), y =  c(row+move), pch = 15, cex = 0.8)
    text(textpos, row, labels = paste(parm,  "[",DF + 1,"]"), cex = 0.75)
    abline(h = row+0.5, lty = 3, col = adjustcolor("black", 0.3))
    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3)
    text(textpos2, row+0.1, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.75) # summary EF and 95% CI's
  }
}

# 4. PlotEO
# Function for plot that compares experimental and observational studies. Plots summary effect sizes Cohen's d (or LnRR) 

# data = Ob
# parm = "Mean grain size"

PlotEO = function(data = df1, parm="Chl-a", start = T, row = 1, textpos = -1, xlow = -1, xhigh = 1, ylow = 0, yhigh = 2, d = FALSE, 
                  mv = 0.3, Col = "black", lab = TRUE, textpos2 = 10) {
  tt = subset(data, data$respvar == parm )            # specify response variable 
  DF = sum(complete.cases(tt$lnRR), na.rm = T) - 1         # degrees of freedom for observations that can be used
  move = mv
  
  if (d == TRUE) {
    # Calculations for standardized mean difference (Cohen's d)
    ll = rma(yi = tt$d, vi = tt$Vd)                           # Use rma function to get tau^2 and other values
    T2 = ll$tau2
    tt$WRand = 1/(tt$Wd + T2)      # Weight = 1/variance within + variance between studies
    tt$WYRand = tt$WRand * tt$d    # Weight * Effect size
    
    M  = round((sum(tt$WYRand, na.rm = T)/sum(tt$WRand, na.rm = T)), digits = 2) # summary EF calculation
    Vm = 1/(sum(tt$WRand, na.rm = T))                          # variance
    SD = sd(tt$g, na.rm = T)                                   # Standard deviation of Cohen's d in subset
    SE = sqrt(Vm)
    LL = round(M - 1.96 * SE, digits = 2)
    UL = round(M + 1.96 * SE, digits = 2)
    
  } else {
    # Calculations for log response ratios (lnRR)
    Q = sum(tt$WY2, na.rm = T) - (sum(tt$WYsquared, na.rm = T))/sum(tt$W, na.rm = T)
    C = sum(tt$W, na.rm = T) - sum(tt$W^2, na.rm = T)/sum(tt$W, na.rm = T)
    T2  = pmax(0, (Q - DF)/C)     # tau squared (between studies variance) never negative
    tt$WRand = 1/(tt$W + T2)      # Weight = 1/variance within + variance between studies
    tt$WYRand = tt$WRand * tt$lnRR
    
    M  = round((sum(tt$WYRand, na.rm = T)/sum(tt$WRand, na.rm = T)), digits = 2) # summary EF calculation
    Vm = 1/(sum(tt$WRand, na.rm = T))                          # variance
    SD = sd(tt$g, na.rm = T)                                   # Standard deviation of Cohen's d in subset
    SE = sqrt(Vm)
    LL = round(M - 1.96 * SE, digits = 2)
    UL = round(M + 1.96 * SE, digits = 2)

  }
  if (lab == TRUE) {
    
  if (start == TRUE) {
    plot(x = c(LL,UL), y = c(row+move,row+move), type = "l", xlim = c(xlow, xhigh), ylim = c(ylow, yhigh), yaxt = "n",
         main = "", xlab = "Summary effect size", ylab = "", col = Col, cex.main = 0.8, cex.lab = 0.8, axes = F)
    points(x = M, y = row+move, pch = 15, cex = 0.8)
    abline(v = 0, lty = 1, col = adjustcolor("grey40", 0.5), lwd = 2)
    text(textpos, row, labels = paste(parm,  "[",DF + 1,"]"), cex = 0.8, col = Col)
    abline(h = c(row+0.65, -0.5), lty = 1)
    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3, col = Col)
    axis(1, cex.axis = 0.8)
    mtext("Parameter [independent observations]", side = 3, at = xlow+5, font = 3, cex = 0.75)
    mtext("Log Response Ratio", side = 3, line = 1, font = 2, cex = 0.8)
    text(textpos2, row+0.25, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.8) # summary EF and 95% CI's
    mtext("ln(RR) [95% CI]", side = 3, at = xhigh-2, font = 3, cex = 0.75)
    
  } else {
    lines(x = c(LL,UL), y = c(row+move,row+move), col = Col)
    points(x = c(M), y =  c(row+move), pch = 15, cex = 0.8)
    text(textpos, row, labels = paste(parm,  "[",DF + 1,"]"), cex = 0.8, col = Col)
    abline(h = row+0.5, lty = 3, col = adjustcolor("black", 0.3))
    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3, col = Col)
    text(textpos2, row+0.25, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.8) # summary EF and 95% CI's
  }
  } else {
    points(x = c(M), y =  c(row+move), pch = 15, col = Col, cex = 0.8)
    text(textpos, row, labels = paste("[",DF + 1,"]"), cex = 0.8, col = Col)
    abline(h = row+0.5, lty = 3, col = adjustcolor("black", 0.3))
    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3, col = Col)
    text(textpos2, row-0.25, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.8, col = Col) # summary EF and 95% CI's
  }
}

# 5. getdats
# Function to retrieve data for a certain variable while calculating effect size. Output is a dataframe used for statistical analysis for 
# Chl-a and OC data

getdats = function(dat = df1, parm="Chl-a", d = FALSE) {
  tt = subset(dat, dat$respvar == parm )            # specify response variable 
  DF = sum(complete.cases(tt$Spooled), na.rm = T) - 1     # degrees of freedom for 'studies'(comparisons) that can be used
  
  if (d == TRUE) {
    # Calculations for standardized mean difference (Cohen's d)
    ll = rma.uni(yi = tt$d, vi = tt$Vd)                           # Use rma function to get tau^2 and other values
    T2 = ll$tau2
    
    tt$WRand = 1/(tt$Wd + T2)      # Weight = 1/variance within + variance between studies
    tt$WYRand = tt$WRand * tt$d    # Weight * Effect size
    M  = (sum(tt$WYRand, na.rm = T)/sum(tt$WRand, na.rm = T))  # summary EF calculation
    Vm = 1/(sum(tt$WRand, na.rm = T))                          # variance
    SD = sd(tt$g, na.rm = T)                                   # Standard deviation of Cohen's d in subset
    SE = sqrt(Vm)
    LL = M - 1.96 * SE
    UL = M + 1.96 * SE
    
    return(tt)
    
  } else {
    # Calculations for log response ratios (lnRR)
    ll = rma.uni(yi = tt$lnRR, vi = tt$VarLnRR)
    T2 = ll$tau2
    tt$WRand = 1/(tt$Wd + T2)      # Weight = 1/variance within + variance between studies
    tt$WYRand = tt$WRand * tt$lnRR
    M  = (sum(tt$WYRand, na.rm = T)/sum(tt$WRand, na.rm = T)) # summary EF calculation
    Vm = 1/(sum(tt$WRand, na.rm = T))                          # variance
    SD = sd(tt$lnRR, na.rm = T)                                    # Standard deviation of Cohen's d in subset
    SE = sqrt(Vm)
    LL = M - 1.96 * SE
    UL = M + 1.96 * SE
    
    return(tt)
  }
}

# 6. Bplot
# univariate boxplot of lnrr.

bplot <- function(df, var, mtitle, xtitle = "", droplev = FALSE, bwex = 0.8, ptscex = 0.6, ...) {
  
  mycol <- rgb(95, 158, 160, max = 255, alpha = 125, names = "Pine Green")
  
  par(mar = c(3,3,3,12))
  boxplot(df$lnRR~df[,var], 
          main = mtitle, ylab = "LnRR", xlab = xtitle,
          col = mycol,
          droplevels = droplev,
          boxwex = bwex,
          outcol = "white",
          cex = 0.8,
          cex.axis = 0.6,
          ...)
  abline(h = 0)
  
  stripchart(df$lnRR~df[,var],              # Data
             method = "jitter", # Random noise
             pch = 19,          # Pch symbols
             cex = ptscex,
             col = "black",     # Color of the symbol
             vertical = TRUE,   # Vertical mode
             add = TRUE) 

  m <- rma.uni(data = df, yi = lnRR, vi = VarLnRR, method = "REML", mods = ~df[,var])
  rsq   <- m$R2
  pval  <- m$QMp
  
  mtext(side = 4, text = paste0("Linear model:"), line = 1, padj = -4, las = 2)  
  mtext(side = 4, text = paste0("Rsq. = ", round(rsq, 5), " %"), line = 1, padj = 0, las = 2)
  mtext(side = 4, text = paste0("pval. = ", round(pval, 5)), line = 1, padj = -2, las = 2) 
}

# 7. Plot continuous variable

    cplot <- function(m, scol, bgcol, xlabz, ylabz = "lnRR", mtitle = "", updown = 0, ...) {
  
  regplot(x = m, 
          shade = mycol2, bg = mycol,
          xlab = xlabz, ylab = "lnRR", main = mtitle,
          cex.axis = 0.7, ...)
  
  r2    <- m$R2
  pint  <- m$pval[1]
  pcont <- m$pval[2]
  coeff <- m$b[2]
  
  mtext(side = 3, text = paste0("Rsq = ", round(r2, 2), " %"), line = -1.5, adj = 0.05, cex = 0.7, padj = updown)
  mtext(side = 3, text = paste0("Sig. intercept = ", round(pint, 4)), line = -2.5, adj = 0.05, cex = 0.7, padj = updown)
  mtext(side = 3, text = paste0("Sig. var. = ", round(pcont, 4)), line = -3.5, adj = 0.05, cex = 0.7, padj = updown)
  mtext(side = 3, text = paste0("Coeff. var. = ", round(coeff, 2)), line = -4.5, adj = 0.05, cex = 0.7, padj = updown)
  
    }
    

# 8. Plot good data + model output
var = "Mean grain size"
model = modMgsslice
avgdf = out
maintitle = "Mean grain size"
miny = -0.5
maxy = 0.7
textpos = 0.5
    
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

# 9. Contvarplot + helper

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
  
  for(i in 1:length(names(lengths))){
    
    sub3  <- subset(sub2, depths == names(lengths)[i])
    
    if(log == TRUE){
      name2 <- log(sub3[,depvar]+1)
    }else{
      name2 <- sub3[,depvar]    
    }
    
    mod   <- rma.mv(yi = lnRR, V = VarLnRR, mods = ~name2, random = ~ 1|article_type_id, data = sub3, method = "REML")
    
    pint  <- mod$pval[1]
    pcont <- mod$pval[2]
    coeff <- mod$b[2]
    
    regplot(x = mod, 
            shade = mycol2, bg = mycol,
            xlab = xlabz, ylab = "lnRR", main = names(lengths)[i],
            cex.axis = 0.7)
    abline(h = 0, col = "darkgrey")
    
    mtext(side = 3, text = paste0("Sig. intercept = ", round(pint, 4)), line = -1.5, adj = 0.05, cex = 0.7)
    mtext(side = 3, text = paste0("Sig. var. = ", round(pcont, 4)), line = -2.5, adj = 0.05, cex = 0.7)
    mtext(side = 3, text = paste0("Coeff. var. = ", round(coeff, 2)), line = -3.5, adj = 0.05, cex = 0.7)
    
  }
  
  return(res)
  
}

# 10. Theme for the violin plots

theme_agile <- function(base_size = 11, base_family = "Arial", lines_lwd = 0.50, plot_grid = TRUE, axis_font = base_family, title_size = base_size*1.2, legend_size = base_size,
                        bg_col = "white",title_font = base_family , base_col  = "black", axis_lines = TRUE,
                        minor_grid = ifelse(plot_grid, TRUE, FALSE), vert_grid = ifelse(plot_grid, TRUE, FALSE), ticks_type = "outer", horz_grid = ifelse(plot_grid, TRUE, FALSE), alpha_leg = 0.1, bord_size = 0,
                        legend_bg = "white", strip_bg = "white", grid_thick = 1,
                        grid_type = "solid", ticks_xy  = "xy", grid_cols = c("grey50", "grey70")){
  theme_bw()+
    ggplot2::theme(
      plot.margin = grid::unit(c(1, 1, .5, .7), "cm"),
      text = ggplot2::element_text(family = base_family, size = base_size),
      axis.line =  element_line(size = ifelse(axis_lines, grid::unit(lines_lwd, "mm"),0), color = "black"),
      axis.ticks.length = grid::unit(ifelse(ticks_type == "outer", 0.15, -0.15), "cm"),
      axis.ticks.x =  element_line(size = ifelse(stringr::str_detect(ticks_xy, "x"), grid::unit(lines_lwd, "cm"),0), color = "black"),
      axis.ticks.y =  element_line(size = ifelse(stringr::str_detect(ticks_xy, "y"), grid::unit(lines_lwd, "cm") ,0), color = "black"),
      axis.text.x = ggplot2::element_text(size = base_size, colour = base_col , family = axis_font,margin=margin(ifelse(ticks_type == "inner", 11, 5),5,10,5,"pt")),
      axis.text.y = ggplot2::element_text(size = base_size, colour = base_col , family = axis_font, margin=margin(5,ifelse(ticks_type == "inner", 11, 5),10,5,"pt")),
      axis.title.y = ggplot2::element_text(size =  base_size, colour = base_col , vjust = 1.5, family = axis_font),
      axis.title.x = ggplot2::element_text(size = base_size,colour = base_col ,vjust = -.5, family = axis_font),
      panel.background = ggplot2::element_rect(fill = bg_col),
      plot.background = ggplot2::element_rect(fill = bg_col),
      panel.border = ggplot2::element_rect(colour = "black", fill=NA, size = bord_size),
      panel.grid.major.x = ggplot2::element_line(linetype = grid_type,colour = ifelse(vert_grid, grid_cols[1],bg_col), size = ifelse(vert_grid,0.25 * grid_thick, 0)),
      panel.grid.minor.x = ggplot2::element_line(linetype = grid_type,colour = ifelse(vert_grid, ifelse(minor_grid, grid_cols[2 - (length(grid_cols) == 1)   ],bg_col),bg_col), size = ifelse(vert_grid,0.15* grid_thick, 0)),
      panel.grid.major.y = ggplot2::element_line(linetype = grid_type,colour = ifelse(horz_grid, grid_cols[1],bg_col), size = ifelse(horz_grid,0.25* grid_thick, 0)),
      panel.grid.minor.y = ggplot2::element_line(linetype = grid_type,colour = ifelse(horz_grid, ifelse(minor_grid, grid_cols[2 - (length(grid_cols) == 1)  ],bg_col),bg_col), size = ifelse(horz_grid,0.15* grid_thick, 0)),
      plot.title = ggplot2::element_text(face="bold",vjust = 2, colour = base_col , size = title_size, family = title_font),
      legend.background = ggplot2::element_rect(fill = scales::alpha(legend_bg, alpha_leg)), legend.key = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend_size, family = base_family),
      legend.title = element_blank(),
      strip.background =  ggplot2::element_rect(colour = strip_bg, fill = strip_bg),
      strip.text.x = ggplot2::element_text(size = base_size + 1),
      strip.text.y = ggplot2::element_text(size = base_size + 1)
    )
}

theme_agile2 <- function(base_size = 11, base_family = "Arial", lines_lwd = 0.50, plot_grid = TRUE, axis_font = base_family, title_size = base_size*1.2, legend_size = base_size,
                        bg_col = "white",title_font = base_family , base_col  = "black", axis_lines = TRUE,
                        minor_grid = ifelse(plot_grid, TRUE, FALSE), vert_grid = ifelse(plot_grid, TRUE, FALSE), ticks_type = "outer", horz_grid = ifelse(plot_grid, TRUE, FALSE), alpha_leg = 0.1, bord_size = 0,
                        legend_bg = "white", strip_bg = "white", grid_thick = 1,
                        grid_type = "solid", ticks_xy  = "xy", grid_cols = c("grey50", "grey70")){
  theme_bw()+
    ggplot2::theme(
      plot.margin = grid::unit(c(1, 1, .5, .7), "cm"),
      text = ggplot2::element_text(family = base_family, size = base_size),
      axis.line =  element_line(size = ifelse(axis_lines, grid::unit(lines_lwd, "mm"),0), color = "black"),
      axis.ticks.length = grid::unit(ifelse(ticks_type == "outer", 0.15, -0.15), "cm"),
      axis.ticks.x =  element_line(size = ifelse(stringr::str_detect(ticks_xy, "x"), grid::unit(lines_lwd, "cm"),0), color = "black"),
      axis.ticks.y =  element_line(size = ifelse(stringr::str_detect(ticks_xy, "y"), grid::unit(lines_lwd, "cm") ,0), color = "black"),
      axis.text.x = ggplot2::element_text(size = base_size, colour = base_col , family = axis_font,margin=margin(ifelse(ticks_type == "inner", 11, 5),5,10,5,"pt")),
      axis.text.y = ggplot2::element_text(size = base_size, colour = base_col , family = axis_font, margin=margin(5,ifelse(ticks_type == "inner", 11, 5),10,5,"pt")),
      axis.title.y = ggplot2::element_text(size =  base_size, colour = base_col , vjust = 1.5, family = axis_font),
      axis.title.x = ggplot2::element_text(size = base_size,colour = base_col ,vjust = 1, family = axis_font),
      panel.background = ggplot2::element_rect(fill = bg_col),
      plot.background = ggplot2::element_rect(fill = bg_col),
      panel.border = ggplot2::element_rect(colour = "black", fill=NA, size = bord_size),
      panel.grid.major.x = ggplot2::element_line(linetype = grid_type,colour = ifelse(vert_grid, grid_cols[1],bg_col), size = ifelse(vert_grid,0.25 * grid_thick, 0)),
      panel.grid.minor.x = ggplot2::element_line(linetype = grid_type,colour = ifelse(vert_grid, ifelse(minor_grid, grid_cols[2 - (length(grid_cols) == 1)   ],bg_col),bg_col), size = ifelse(vert_grid,0.15* grid_thick, 0)),
      panel.grid.major.y = ggplot2::element_line(linetype = grid_type,colour = ifelse(horz_grid, grid_cols[1],bg_col), size = ifelse(horz_grid,0.25* grid_thick, 0)),
      panel.grid.minor.y = ggplot2::element_line(linetype = grid_type,colour = ifelse(horz_grid, ifelse(minor_grid, grid_cols[2 - (length(grid_cols) == 1)  ],bg_col),bg_col), size = ifelse(horz_grid,0.15* grid_thick, 0)),
      plot.title = ggplot2::element_text(face="bold",vjust = 2, colour = base_col , size = title_size, family = title_font, hjust = 0.5),
      legend.background = ggplot2::element_rect(fill = scales::alpha(legend_bg, alpha_leg)), legend.key = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend_size, family = base_family),
      legend.title = element_blank(),
      strip.background =  ggplot2::element_rect(colour = strip_bg, fill = strip_bg),
      strip.text.x = ggplot2::element_text(size = base_size + 1),
      strip.text.y = ggplot2::element_text(size = base_size + 1)
    )
}
