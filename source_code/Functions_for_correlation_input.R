#-------------------------------------------------------------------------------
#        Functions for inputting data and calculating Fishers Z effect sizes (r) 
#                       for trawling gradient studies
#                                Justin Tiano                                 
#-------------------------------------------------------------------------------
#
# For trawling gradient studies, log response ratio, and directly calculating standardized 
# mean differences is not possible. A type of correlational effect size (Fishers Z or 'r'), 
# here was used to assess potential patterns between in the studies between bottom trawl intensity
# and biogeochemical parameters. These values are then converted to to standardized mean 
# differences (Cohn's d) to be able to compare with observational and experimental studies.
#
# CONTENTS:
# 1. function Input
# 2. function Input2
# 3. function Input3
# 4. function plotR
# 5. function plotR2
# 6. function getdat
#-------------------------------------------------------------------------------

# Function to input data

Input = function(ID = NULL, var1 = NULL, var2 = NULL, var3 = NULL, var4 = NULL, var5 = NULL, 
                var6 = NULL, var7 = NULL, var8 = NULL, var9 = NULL, var10 = NULL, var11 = NULL,
                var12 = NULL, var13  = NULL, var14 = NULL, var15 = NULL, var16 = NULL, var17 = NULL,
                subset = FALSE, subname = NULL){
  
  if (subset == TRUE) {
  
    study = subname
    S = study
    S$r        = NA
    S$r.pvalue = NA
    S$FishersZ = NA
    S$SEz      = NA
    
    # 1st variable 
    S$r[S$Response.variable==var1$var]        = var1[2]
    S$r.pvalue[S$Response.variable==var1$var] = var1[3]
    S$FishersZ[S$Response.variable==var1$var] = var1[4]
    S$SEz[S$Response.variable==var1$var]      = var1[5]
    # 2nd variable
    S$r[S$Response.variable==var2$var]        = var2[2]
    S$r.pvalue[S$Response.variable==var2$var] = var2[3]
    S$FishersZ[S$Response.variable==var2$var] = var2[4]
    S$SEz[S$Response.variable==var2$var]      = var2[5]
    # 3rd variable
    S$r[S$Response.variable==var3$var]        = var3[2]
    S$r.pvalue[S$Response.variable==var3$var] = var3[3]
    S$FishersZ[S$Response.variable==var3$var] = var3[4]
    S$SEz[S$Response.variable==var3$var]      = var3[5]
    # 3rd variable
    S$r[S$Response.variable==var4$var]        = var4[2]
    S$r.pvalue[S$Response.variable==var4$var] = var4[3]
    S$FishersZ[S$Response.variable==var4$var] = var4[4]
    S$SEz[S$Response.variable==var4$var]      = var4[5]
    # 5th variable
    S$r[S$Response.variable==var5$var]        = var5[2]
    S$r.pvalue[S$Response.variable==var5$var] = var5[3]
    S$FishersZ[S$Response.variable==var5$var] = var5[4]
    S$SEz[S$Response.variable==var5$var]      = var5[5]
    # 6th variable
    S$r[S$Response.variable==var6$var]        = var6[2]
    S$r.pvalue[S$Response.variable==var6$var] = var6[3]
    S$FishersZ[S$Response.variable==var6$var] = var6[4]
    S$SEz[S$Response.variable==var6$var]      = var6[5]
    # 7th variable
    S$r[S$Response.variable==var7$var]        = var7[2]
    S$r.pvalue[S$Response.variable==var7$var] = var7[3]
    S$FishersZ[S$Response.variable==var7$var] = var7[4]
    S$SEz[S$Response.variable==var7$var]      = var7[5]
    # 8th variable
    S$r[S$Response.variable==var8$var]        = var8[2]
    S$r.pvalue[S$Response.variable==var8$var] = var8[3]
    S$FishersZ[S$Response.variable==var8$var] = var8[4]
    S$SEz[S$Response.variable==var8$var]      = var8[5]
    # 8th variable
    S$r[S$Response.variable==var9$var]        = var9[2]
    S$r.pvalue[S$Response.variable==var9$var] = var9[3]
    S$FishersZ[S$Response.variable==var9$var] = var9[4]
    S$SEz[S$Response.variable==var9$var]      = var9[5]
    # 10th variable
    S$r[S$Response.variable==var10$var]        = var10[2]
    S$r.pvalue[S$Response.variable==var10$var] = var10[3]
    S$FishersZ[S$Response.variable==var10$var] = var10[4]
    S$SEz[S$Response.variable==var10$var]      = var10[5]
    # 11th variable
    S$r[S$Response.variable==var11$var]        = var11[2]
    S$r.pvalue[S$Response.variable==var11$var] = var11[3]
    S$FishersZ[S$Response.variable==var11$var] = var11[4]
    S$SEz[S$Response.variable==var11$var]      = var11[5]
    # 12th variable
    S$r[S$Response.variable==var12$var]        = var12[2]
    S$r.pvalue[S$Response.variable==var12$var] = var12[3]
    S$FishersZ[S$Response.variable==var12$var] = var12[4]
    S$SEz[S$Response.variable==var12$var]      = var12[5]
    # 13th variable
    S$r[S$Response.variable==var13$var]        = var13[2]
    S$r.pvalue[S$Response.variable==var13$var] = var13[3]
    S$FishersZ[S$Response.variable==var13$var] = var13[4]
    S$SEz[S$Response.variable==var13$var]      = var13[5]
    # 14th variable
    S$r[S$Response.variable==var14$var]        = var14[2]
    S$r.pvalue[S$Response.variable==var14$var] = var14[3]
    S$FishersZ[S$Response.variable==var14$var] = var14[4]
    S$SEz[S$Response.variable==var14$var]      = var14[5]
    # 15th variable
    S$r[S$Response.variable==var15$var]        = var15[2]
    S$r.pvalue[S$Response.variable==var15$var] = var15[3]
    S$FishersZ[S$Response.variable==var15$var] = var15[4]
    S$SEz[S$Response.variable==var15$var]      = var15[5]
    # 16th variable
    S$r[S$Response.variable==var15$var]        = var16[2]
    S$r.pvalue[S$Response.variable==var15$var] = var16[3]
    S$FishersZ[S$Response.variable==var15$var] = var16[4]
    S$SEz[S$Response.variable==var15$var]      = var16[5]
    # 17th variable
    S$r[S$Response.variable==var15$var]        = var17[2]
    S$r.pvalue[S$Response.variable==var15$var] = var17[3]
    S$FishersZ[S$Response.variable==var15$var] = var17[4]
    S$SEz[S$Response.variable==var15$var]      = var17[5]
    
    return(S)
    
  } else {
  
    study = subset(df, df$article_id == ID)
    S = study
    S$r        = NA
    S$r.pvalue = NA
    S$FishersZ = NA
    S$SEz      = NA
    
    # 1st variable 
    S$r[S$Response.variable==var1$var]        = var1[2]
    S$r.pvalue[S$Response.variable==var1$var] = var1[3]
    S$FishersZ[S$Response.variable==var1$var] = var1[4]
    S$SEz[S$Response.variable==var1$var]      = var1[5]
    # 2nd variable
    S$r[S$Response.variable==var2$var]        = var2[2]
    S$r.pvalue[S$Response.variable==var2$var] = var2[3]
    S$FishersZ[S$Response.variable==var2$var] = var2[4]
    S$SEz[S$Response.variable==var2$var]      = var2[5]
    # 3rd variable
    S$r[S$Response.variable==var3$var]        = var3[2]
    S$r.pvalue[S$Response.variable==var3$var] = var3[3]
    S$FishersZ[S$Response.variable==var3$var] = var3[4]
    S$SEz[S$Response.variable==var3$var]      = var3[5]
    # 3rd variable
    S$r[S$Response.variable==var4$var]        = var4[2]
    S$r.pvalue[S$Response.variable==var4$var] = var4[3]
    S$FishersZ[S$Response.variable==var4$var] = var4[4]
    S$SEz[S$Response.variable==var4$var]      = var4[5]
    # 5th variable
    S$r[S$Response.variable==var5$var]        = var5[2]
    S$r.pvalue[S$Response.variable==var5$var] = var5[3]
    S$FishersZ[S$Response.variable==var5$var] = var5[4]
    S$SEz[S$Response.variable==var5$var]      = var5[5]
    # 6th variable
    S$r[S$Response.variable==var6$var]        = var6[2]
    S$r.pvalue[S$Response.variable==var6$var] = var6[3]
    S$FishersZ[S$Response.variable==var6$var] = var6[4]
    S$SEz[S$Response.variable==var6$var]      = var6[5]
    # 7th variable
    S$r[S$Response.variable==var7$var]        = var7[2]
    S$r.pvalue[S$Response.variable==var7$var] = var7[3]
    S$FishersZ[S$Response.variable==var7$var] = var7[4]
    S$SEz[S$Response.variable==var7$var]      = var7[5]
    # 8th variable
    S$r[S$Response.variable==var8$var]        = var8[2]
    S$r.pvalue[S$Response.variable==var8$var] = var8[3]
    S$FishersZ[S$Response.variable==var8$var] = var8[4]
    S$SEz[S$Response.variable==var8$var]      = var8[5]
    # 8th variable
    S$r[S$Response.variable==var9$var]        = var9[2]
    S$r.pvalue[S$Response.variable==var9$var] = var9[3]
    S$FishersZ[S$Response.variable==var9$var] = var9[4]
    S$SEz[S$Response.variable==var9$var]      = var9[5]
    # 10th variable
    S$r[S$Response.variable==var10$var]        = var10[2]
    S$r.pvalue[S$Response.variable==var10$var] = var10[3]
    S$FishersZ[S$Response.variable==var10$var] = var10[4]
    S$SEz[S$Response.variable==var10$var]      = var10[5]
    # 11th variable
    S$r[S$Response.variable==var11$var]        = var11[2]
    S$r.pvalue[S$Response.variable==var11$var] = var11[3]
    S$FishersZ[S$Response.variable==var11$var] = var11[4]
    S$SEz[S$Response.variable==var11$var]      = var11[5]
    # 12th variable
    S$r[S$Response.variable==var12$var]        = var12[2]
    S$r.pvalue[S$Response.variable==var12$var] = var12[3]
    S$FishersZ[S$Response.variable==var12$var] = var12[4]
    S$SEz[S$Response.variable==var12$var]      = var12[5]
    # 13th variable
    S$r[S$Response.variable==var13$var]        = var13[2]
    S$r.pvalue[S$Response.variable==var13$var] = var13[3]
    S$FishersZ[S$Response.variable==var13$var] = var13[4]
    S$SEz[S$Response.variable==var13$var]      = var13[5]
    # 14th variable
    S$r[S$Response.variable==var14$var]        = var14[2]
    S$r.pvalue[S$Response.variable==var14$var] = var14[3]
    S$FishersZ[S$Response.variable==var14$var] = var14[4]
    S$SEz[S$Response.variable==var14$var]      = var14[5]
    # 15th variable
    S$r[S$Response.variable==var15$var]        = var15[2]
    S$r.pvalue[S$Response.variable==var15$var] = var15[3]
    S$FishersZ[S$Response.variable==var15$var] = var15[4]
    S$SEz[S$Response.variable==var15$var]      = var15[5]
    # 16th variable
    S$r[S$Response.variable==var15$var]        = var16[2]
    S$r.pvalue[S$Response.variable==var15$var] = var16[3]
    S$FishersZ[S$Response.variable==var15$var] = var16[4]
    S$SEz[S$Response.variable==var15$var]      = var16[5]
    # 17th variable
    S$r[S$Response.variable==var15$var]        = var17[2]
    S$r.pvalue[S$Response.variable==var15$var] = var17[3]
    S$FishersZ[S$Response.variable==var15$var] = var17[4]
    S$SEz[S$Response.variable==var15$var]      = var17[5]
    
    return(S)
  
  }
}

Input2 = function(ID = study,rows = 17){
  
  S <- data.frame(matrix(NA,    # Create empty data frame
                         nrow = rows,
                         ncol = 9))
  names(S) = names(ID[c(1,2,11,13,25,26,30,24,34)])
  S[1:rows,] = ID[1,c(1,2,11,13,25,26,30,24,34)] # fill in data
  
  return(S)
}

Input3 = function(S, var1 = NULL, var2 = NULL, var3 = NULL, var4 = NULL, var5 = NULL, 
                  var6 = NULL, var7 = NULL, var8 = NULL, var9 = NULL, var10 = NULL, var11 = NULL,
                  var12 = NULL, var13  = NULL, var14 = NULL, var15 = NULL, var16 = NULL, var17 = NULL){
  # 1st variable
  S$Response.variable[1] = var1[["var"]]
  S$r[1]                 = var1[[2]]
  S$n[1]                 = var1[[6]]
  S$V[1]                 = var1[[7]]
  S$SE[1]                = var1[[5]]
  # 2nd variable
  S$Response.variable[2] = var2[["var"]]
  S$r[2]                 = var2[["r"]]
  S$n[2]                 = var2[["n"]]
  S$V[2]                 = var2[["V"]]
  S$SE[2]                = var2[["SEz"]]
  # 3rd variable
  S$Response.variable[3] = var3$var
  S$r[3]                 = var3$r
  S$n[3]                 = var3$n
  S$V[3]                 = var3$V
  S$SE[3]                = var3$SEz
  # 4th variable 
  S$Response.variable[4] = var4$var
  S$r[4]                 = var4$r
  S$n[4]                 = var4$n
  S$V[4]                 = var4$V
  S$SE[4]                = var4$SEz
  
  S = S[1:length(rows),]
  return(S)
  
}      

# Function to calculate the correlational effect size
plotR = function(data , Funit = "hours/km2", var = "TOC", plot = TRUE){
  c = subset(data, data$covariates.Unit==Funit)
  X = as.numeric(c$covariates.Value[c$Response.variable==var])
  Y = as.numeric(c$Response.value[c$Response.variable==var])
  a   = cor.test(X,Y, method = "pearson")
  r   = a$estimate                                   # correlation value
  n   = sum(c$replicates[c$Response.variable==var])  # number of replicates
  z   = 0.5 * log((1 + r)/(1 - r))                   # Fishers Z transformation
  Vz  = 1/(n-3)                                      # variance of z
  SEz = sqrt(Vz)                                     # standard error of z
  r   = (exp(2 * z) - 1)/ (exp(2 * z) + 1)           # weighted r value 
  dig = 3
  
  if (plot == TRUE) {
    
    plot(X, Y, main = var)
    abline(lm(Y~X))
    legend("topleft", bty = "n", c(paste("r=",round(r, digits = dig)),
                                   paste("p=",round(a$p.value,digits = dig)),
                                   paste("z=",round((z),digits = dig)),
                                   paste("SEz=",round((SEz),digits = dig),
                                         paste("n=",n))))
    print(paste("r=",round(r, digits = dig)))   # print values
    print(paste("p=",round(a$p.value, digits = dig)))
    print(paste("z=",round(z, digits = dig)))
    print(paste("SEz=",round(SEz, digits = dig)))
    print(paste("n=",n))
    print(var)
    
    return(data.frame(var = var, r = r, p.value = a$p.value, z =  z, SEz = SEz, n = n, V = Vz))
    
  } else {
    
    return(data.frame(var = var, r = r, p.value = a$p.value, z =  z, SEz = SEz, n = n, V = Vz))
    
  }

}   

# 2nd effect size function for special cases
plotR2 = function(data , Funit = "hours/km2", var = "TOC", plot = TRUE){
  c = subset(data, data$FEunit==Funit)
  X = as.numeric(c$covariates.Value[c$Response.variable==var])
  Y = as.numeric(c$Response.value[c$Response.variable==var])
  a   = cor.test(X,Y, method = "pearson")
  r   = a$estimate                                   # correlation value
  n   = sum(c$replicates[c$Response.variable==var])  # number of replicates
  z   = 0.5 * log((1 + r)/(1 - r))                     # Fishers Z transformation
  Vz  = 1/(n-3)                                      # variance of z
  SEz = sqrt(Vz)                                     # standard error of z
  r   = (exp(2 * z) - 1)/ (exp(2 * z) + 1)           # weighted r value 
  dig = 3
  
  if (plot == TRUE) {
    
    plot(X, Y, main = var)
    abline(lm(Y~X))
    legend("topleft", bty = "n", c(paste("r=",round(r, digits = dig)),
                                   paste("p=",round(a$p.value,digits = dig)),
                                   paste("z=",round((z),digits = dig)),
                                   paste("SEz=",round((SEz),digits = dig),
                                         paste("n=",n))))
    print(paste("r=",round(r, digits = dig)))   # print values
    print(paste("p=",round(a$p.value, digits = dig)))
    print(paste("z=",round(z, digits = dig)))
    print(paste("SEz=",round(SEz, digits = dig)))
    print(paste("n=",n))
    print(var)
    
    return(data.frame(var = var, r = r, p.value = a$p.value, z =  z, SEz = SEz, n = n, V = Vz))
    
  } else {
    
    return(data.frame(var = var, r = r, p.value = a$p.value, z =  z, SEz = SEz, n = n, V = Vz))
    
  }
  
}

# simple function to calculate r effect size
getdat = function(data = study, Funit = "Teffort", var = "TOC"){
  c = subset(data, data$FEunit==Funit)
  a = cor.test(as.numeric(c$Feffort), c$Response.value, method = "pearson")
  r   = a$estimate                                   # correlation value
  n   = sum(c$replicates[c$Response.variable==var])  # number of replicates
  z   = 0.5 * log((1 + r)/(1 - r))                     # Fishers Z transformation
  Vz  = 1/(n-3)                                      # variance of z
  SEz = sqrt(Vz)                                     # standard error of z
  r   = (exp(2 * z) - 1)/ (exp(2 * z) + 1)           # weighted r value 
  return(data.frame(var = var, r = r, p.value = a$p.value, z =  z, SEz = SEz, n = n, V = Vz))
}