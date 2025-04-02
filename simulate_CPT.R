pow <- function(b, e) return(b^e)

simulate_CPT_individ <- function(alpha, beta, gamma, delta, lambda, luce, n=NULL) {
  if (is.null(n)) n <- 60
  v.x.a <- v.y.a <-w.x.a <-w.y.a <- z.a <- Vf.a <- numeric(60) 
  v.x.b <- v.y.b <-w.x.b <-w.y.b <- z.b <- Vf.b <- numeric(60)
  # Item-Loop, mixed gambles,gamble A
    v.x.a <- pow(prospects.a[,3],alpha)       
    v.y.a <- (-1 * lambda) * pow((-1 * prospects.a[,1]),beta)       
    
    z.a   <- pow(prospects.a[,4],gamma) + pow(prospects.a[,2],delta) 
    w.x.a <- pow(prospects.a[,4],gamma) / pow(z.a,(1/gamma)) 
    w.y.a <- pow(prospects.a[,2],delta) / pow(z.a,(1/delta)) 
    
    Vf.a  <- w.x.a * v.x.a + w.y.a * v.y.a
  
    # Item-Loop, mixed gambles,gamble B
  
    v.x.b <- pow(prospects.b[,3],alpha)       
    v.y.b <- (-1 * lambda) * pow((-1 * prospects.b[,1]),beta)       
    
    z.b   <- pow(prospects.b[,4],gamma) + pow(prospects.b[,2],delta) 
    w.x.b <- pow(prospects.b[,4],gamma) / pow(z.b,(1/gamma)) 
    w.y.b <- pow(prospects.b[,2],delta) / pow(z.b,(1/delta)) 
    
    Vf.b  <- w.x.b * v.x.b + w.y.b * v.y.b
    
    binval <- (1)/(1+exp((-1*luce)*(Vf.b-Vf.a)))
    rawdata = as.numeric(runif(n)<binval)
    return(rawdata)
}
