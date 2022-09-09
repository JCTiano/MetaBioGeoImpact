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

PlotSum = function(parm="Chl-a", start = F, row = 1, textpos = -1, xlow = -1, xhigh = 1, ylow = 0, yhigh = 2, d = FALSE, textpos2 = 2) {
  tt = subset(df1, df1$Response.variable==parm )            # specify response variable 
  DF = sum(complete.cases(tt$lnRR), na.rm = T) - 1       # degrees of freedom for observations that can be used
  
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
    text(textpos, row, labels = paste(parm,  "[",DF + 1," ]"), cex = 0.75)
    arrows(x0 = LL, x1 = UL, y0 = c(row,row), length = 0.03, angle = -90, code = 3)
    text(textpos2, row, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.75) # summary EF and 95% CI's
    mtext("Parameter[independent observations]", side = 3, at = xlow+6, font = 3, cex = 0.75)
    mtext("ln(RR)[95% CI]", side = 3, at = xhigh-2, font = 3, cex = 0.75)
    mtext("Log Response Ratio", side = 3, line = 1, font = 2, cex = 0.8)
    
  } else {
    lines(x = c(LL,UL), y = c(row,row))
    points(x = c(M), y =  c(row), pch = 15, cex = 0.8)
    text(textpos, row, labels = paste(parm,  "[",DF + 1,"]"), cex = 0.75)
    abline(h = row+0.5, lty = 3, col = adjustcolor("black", 0.3))
    arrows(x0 = LL, x1 = UL, y0 = c(row,row), length = 0.03, angle = -90, code = 3)
    text(textpos2, row, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.75) # summary EF and 95% CI's
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
    text(textpos, row, labels = paste("[",DF + 1," ]"), cex = 0.75, col = Col)

    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3, col = Col)
    mtext("Parameter[independent observations]", side = 3, at = xlow+8, font = 3, cex = 0.75)
    mtext("Standardized mean difference", side = 3, line = 1, font = 2, cex = 0.8)
    mtext("Cohen's d [95% CI]", side = 3, at = xhigh-5, font = 3, cex = 0.75)
    text(textpos2, row-0.1, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.75, col = Col) # summary EF and 95% CI's
  } else {
    lines(x = c(LL,UL), y = c(row+move,row+move), col = Col)
    points(x = c(M), y =  c(row+move), pch = 15, col = Col, cex = 0.8)
    text(textpos, row, labels = paste("[",DF + 1,"]"), cex = 0.75, col = Col)
    abline(h = row+0.5, lty = 3, col = adjustcolor("black", 0.3))
    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3, col = Col)
    text(textpos2, row-0.1, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.75, col = Col) # summary EF and 95% CI's
  }
}
# 3. PlotD
# Function for plot that compares gradient studies and comparative studies (input for comparative studies). Plots summary effect sizes Cohen's d

PlotD = function(parm="Chl-a", start = T, row = 1, textpos = -1, xlow = -1, xhigh = 1, ylow = 0, yhigh = 2, d = FALSE, mv = 0.2,
                 textpos2 = 10) {
  tt = subset(df1, df1$Response.variable==parm )            # specify response variable 
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

PlotEO = function(data = df1, parm="Chl-a", start = T, row = 1, textpos = -1, xlow = -1, xhigh = 1, ylow = 0, yhigh = 2, d = FALSE, 
                  mv = 0.3, Col = "black", lab = TRUE, textpos2 = 10) {
  tt = subset(data, data$Response.variable==parm )            # specify response variable 
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
    text(textpos, row, labels = paste(parm,  "[",DF + 1," ]"), cex = 0.75, col = Col)
    abline(h = c(row+0.65, -0.5), lty = 1)
    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3, col = Col)
    axis(1, cex.axis = 0.8)
    mtext("Parameter[independent observations]", side = 3, at = xlow+5, font = 3, cex = 0.75)
    mtext("Log Response Ratio", side = 3, line = 1, font = 2, cex = 0.8)
    text(textpos2, row+0.25, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.75) # summary EF and 95% CI's
    mtext("ln(RR)[95% CI]", side = 3, at = xhigh-2, font = 3, cex = 0.75)
    
  } else {
    lines(x = c(LL,UL), y = c(row+move,row+move), col = Col)
    points(x = c(M), y =  c(row+move), pch = 15, cex = 0.8)
    text(textpos, row, labels = paste(parm,  "[",DF + 1,"]"), cex = 0.75, col = Col)
    abline(h = row+0.5, lty = 3, col = adjustcolor("black", 0.3))
    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3, col = Col)
    text(textpos2, row+0.25, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.75) # summary EF and 95% CI's
  }
  } else {
    points(x = c(M), y =  c(row+move), pch = 15, col = Col, cex = 0.8)
    text(textpos, row, labels = paste("[",DF + 1,"]"), cex = 0.75, col = Col)
    abline(h = row+0.5, lty = 3, col = adjustcolor("black", 0.3))
    arrows(x0 = LL, x1 = UL, y0 = c(row+move,row+move), length = 0.03, angle = -90, code = 3, col = Col)
    text(textpos2, row-0.25, labels = paste(M,"[",LL,",",UL,"]"), cex = 0.75, col = Col) # summary EF and 95% CI's
  }
}

# 5. getdats
# Function to retrieve data for a certain variable while calculating effect size. Output is a dataframe used for statistical analysis for 
# Chl-a and OC data

getdats = function(dat = df1, parm="Chl-a", d = FALSE) {
  tt = subset(dat, dat$Response.variable==parm )            # specify response variable 
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