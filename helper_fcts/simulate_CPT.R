pow <- function(b, e) return(b^e)

simulate_CPT_individ <- function(alpha, beta, gamma.gain, gamma.loss, lambda, sens, n=NULL) {
  if (is.null(n)) n <- 60
  v.x.a <- v.y.a <-w.x.a <-w.y.a <- z.a <- Vf.a <- numeric(60) 
  v.x.b <- v.y.b <-w.x.b <-w.y.b <- z.b <- Vf.b <- numeric(60)
  # Item-Loop, mixed gambles,gamble A
    v.x.a <- pow(mixed_prospects.a[,3],alpha)       
    v.y.a <- (-1 * lambda) * pow((-1 * mixed_prospects.a[,1]),beta)       
    
    z.a   <- pow(mixed_prospects.a[,4],gamma.gain) + pow(mixed_prospects.a[,2],gamma.loss) 
    w.x.a <- pow(mixed_prospects.a[,4],gamma.gain) / pow(z.a,(1/gamma.gain)) 
    w.y.a <- pow(mixed_prospects.a[,2],gamma.loss) / pow(z.a,(1/gamma.loss)) 
    
    Vf.a  <- w.x.a * v.x.a + w.y.a * v.y.a
  
    # Item-Loop, mixed gambles,gamble B
  
    v.x.b <- pow(mixed_prospects.b[,3],alpha)       
    v.y.b <- (-1 * lambda) * pow((-1 * mixed_prospects.b[,1]),beta)       
    
    z.b   <- pow(mixed_prospects.b[,4],gamma.gain) + pow(mixed_prospects.b[,2],gamma.loss) 
    w.x.b <- pow(mixed_prospects.b[,4],gamma.gain) / pow(z.b,(1/gamma.gain)) 
    w.y.b <- pow(mixed_prospects.b[,2],gamma.loss) / pow(z.b,(1/gamma.loss)) 
    
    Vf.b  <- w.x.b * v.x.b + w.y.b * v.y.b
    
    binval <- (1)/(1+exp((-1*sens)*(Vf.b-Vf.a)))
    rawdata = as.numeric(runif(n)<binval)
    return(rawdata)
}
