# The following function is from (http://menugget.blogspot.com/2012/10/dineof-data-interpolating-empirical.html)

#This function is based on the DINEOF (Data Interpolating Empirical Orthogonal Functions)
#procedure described by Beckers and Rixon (2003).
#
#The arguments are:
#Xo - a gappy data field
#nu - a maximum number of EOFs to iterate (leave equalling "NULL" if algorithm shold proceed until convergence)
#ref.pos - a vector of non-gap reference positions by which errors will be assessed via root mean squared error ("RMS"). 
#If ref.pos = NULL, then either 30 or 1% of the non-gap values (which ever is larger) will be sampled at random.
#delta.rms - is the threshold for RMS convergence.
#
#The results object includes:
#Xa - the data field with interpolated values (via EOF reconstruction) included.
#n.eof - the number of EOFs used in the final solution
#RMS - a vector of the RMS values from the iteration
#NEOF - a vector of the number of EOFs used at each iteration
#
#Beckers, Jean-Marie, and M. Rixen. "EOF Calculations and Data Filling from Incomplete Oceanographic Datasets."
#Journal of Atmospheric and Oceanic Technology 20.12 (2003): 1839-1856.
#
#4.0: Incorporates irlba algorithm and is only the dineof procedure
#
dineof <- function(Xo, n.max=NULL, ref.pos=NULL, delta.rms=1e-5){
  
  library(irlba)
  
  if(is.null(n.max)){
    n.max=dim(Xo)[2]
  } 
  
  na.true <- which(is.na(Xo))
  na.false <- which(!is.na(Xo))
  if(is.null(ref.pos)) ref.pos <- sample(na.false, max(30, 0.01*length(na.false)))
  
  Xa <- replace(Xo, c(ref.pos, na.true), 0)
  rms.prev <- Inf
  rms.now <- sqrt(mean((Xa[ref.pos] - Xo[ref.pos])^2))
  n.eof <- 1
  RMS <- rms.now
  NEOF <- n.eof
  Xa.best <- Xa
  n.eof.best <- n.eof 
  while(rms.prev - rms.now > delta.rms & n.max > n.eof){ #loop for increasing number of eofs
    while(rms.prev - rms.now > delta.rms){ #loop for replacement
      rms.prev <- rms.now
      SVDi <- irlba(Xa, nu=n.eof, nv=n.eof) 
      RECi <- as.matrix(SVDi$u[,seq(n.eof)]) %*% as.matrix(diag(SVDi$d[seq(n.eof)], n.eof, n.eof)) %*% t(as.matrix(SVDi$v[,seq(n.eof)]))
      Xa[c(ref.pos, na.true)] <- RECi[c(ref.pos, na.true)]
      rms.now <- sqrt(mean((Xa[ref.pos] - Xo[ref.pos])^2))
      print(paste(n.eof, "EOF", "; RMS =", round(rms.now, 8)))
      RMS <- c(RMS, rms.now)
      NEOF <- c(NEOF, n.eof)
      gc()
      if(rms.now == min(RMS)) {
        Xa.best <- Xa
        n.eof.best <- n.eof
      }
    }
    n.eof <- n.eof + 1
    rms.prev <- rms.now
    SVDi <- irlba(Xa, nu=n.eof, nv=n.eof) 
    RECi <- as.matrix(SVDi$u[,seq(n.eof)]) %*% as.matrix(diag(SVDi$d[seq(n.eof)], n.eof, n.eof)) %*% t(as.matrix(SVDi$v[,seq(n.eof)]))
    Xa[c(ref.pos, na.true)] <- RECi[c(ref.pos, na.true)]
    rms.now <- sqrt(mean((Xa[ref.pos] - Xo[ref.pos])^2))
    print(paste(n.eof, "EOF", "; RMS =", round(rms.now, 8)))
    RMS <- c(RMS, rms.now)
    NEOF <- c(NEOF, n.eof)
    gc()
    if(rms.now == min(RMS)) {
      Xa.best <- Xa
      n.eof.best <- n.eof
    }
  }
  
  Xa <- Xa.best
  n.eof <- n.eof.best
  rm(list=c("Xa.best", "n.eof.best", "SVDi", "RECi"))
  
  Xa[ref.pos] <- Xo[ref.pos]
  
  RESULT <- list(
    Xa=Xa, n.eof=n.eof, RMS=RMS, NEOF=NEOF
  )
  
  RESULT
}
