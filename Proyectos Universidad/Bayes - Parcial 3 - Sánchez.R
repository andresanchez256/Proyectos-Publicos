#####--- Parcial No. 3 , Estadistica Bayesiana ---####

# Andrés Felipe Sánchez Felipe #

## Librerias utilizadas ##

library(coda)
library(xtable)

#### Valores y Datos Iniciales ####
x <- c(33, 14, 27, 90, 12, 17)   # Poblacion de las 6 ciudades, en miles
y <- c(1, 3, 2, 12, 1, 1)        # Casos de la condicion no genetica 
n = length(y)

#### MODELO 1 ####
## Hiperparámetros ##
aalpha <- 1
balpha <- 1
abeta <- 10
bbeta <- 1 
ac<-0

## Algoritmo de MCMC del modelo 1 ##
S = 55000
THETA1 <- matrix(NA, nrow = S, ncol = n)
AB <- matrix(NA, nrow = S, ncol = 2)
LLm1 <- NA
ncat <- floor(S/10) 

# Valores iniciales de la cadena
set.seed(16)
delta <- 1.4
alpha<-rgamma(1,aalpha,balpha)
beta<-rgamma(1,abeta,bbeta)
theta<- NULL

#### Implementacion manual del MCMC del M1 ####
set.seed(2409)  # Establecer una semilla
for (s in 1:S) {
  # Actualizacion de los Thetas #
  for (i in 1:n) {
    theta[i] <- rgamma(1,shape = alpha + y[i], rate = beta + x[i])
  }
  
  # Actualizacion de los Alpha - Algoritmo de Metropolis #
  
  ### Propuesta ###
  alpha.prop <- abs(runif(1,alpha - delta,alpha + delta))
  
  ### Tasa de aceptacion en escala log ### 
  log_r <- sum(dgamma(theta, shape = alpha.prop, rate = beta, log = T)) +
    dgamma(alpha.prop, shape = aalpha, rate = balpha, log = T) - 
    sum(dgamma(theta, shape = alpha, rate = beta, log = T)) - 
    dgamma(alpha, shape = aalpha, rate = balpha, log = T)
  
  ### --- Actualizacion --- ###
  if (log(runif(1)) < log_r) { 
    alpha <- alpha.prop
    ac <- ac + 1
  }
  
  # Actualizacion de los Betas #
  beta <- rgamma(1, shape = abeta + n*alpha, rate = bbeta + sum(theta))
  
  # Almacenamiento de los parametros #
  THETA1[s,] <- theta
  AB[s,] <- c(alpha,beta)
  
  # Log-verosimilitud #
  LLm1[s]<-sum(dpois(x = y, lambda = theta*x, log = T))
  
  # Contador para el algoritmo #
  if (s%%ncat == 0) cat("El MCMC del 1er Modelo está ", 100*round(s/S, 1), "% completado \n", sep = "")
}

(Tasa <- ac/S) #Tasa de aceptación para alpha
AB<-AB[-c(1:5000),]
THETA1<-THETA1[-c(1:5000),]

# Tamaños efectivos de muestra
(ESSTHETA1 <- coda::effectiveSize(THETA1))
(ESSAB <- coda::effectiveSize(AB))

# Coeficientes de Variacion MC
(SETHETA1 <- apply(THETA1,2,sd)/sqrt(ESSTHETA1))
(SEAB <- apply(AB,2,sd)/sqrt(ESSAB))
(CVTHETA1 <- SETHETA1/colMeans(THETA1))
(CVAB <- SEAB/colMeans(AB))

ESS1<-c(ESSAB,ESSTHETA1)
CV1<-c(CVAB,CVTHETA1)*100
Tab1 <- t(rbind(ESS1,CV1))
rownames(Tab1) <- c('alpha','beta','theta_1','theta_2','theta_3','theta_4','theta_5','theta_6')
colnames(Tab1) <- c('T.E. de Muestra', 'C.V. en%')
Tab1

xtable(Tab1)
# Diagramas de seguimiento
windows()
par(mfrow=c(4,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(AB[,1],type = 'l',ylab = expression(alpha))
plot(AB[,2],type = 'l',ylab = expression(beta))
plot(THETA1[,1],type = 'l',ylab = expression(theta[1]))
plot(THETA1[,2],type = 'l',ylab = expression(theta[2]))
plot(THETA1[,3],type = 'l',ylab = expression(theta[3]))
plot(THETA1[,4],type = 'l',ylab = expression(theta[4]))
plot(THETA1[,5],type = 'l',ylab = expression(theta[5]))
plot(THETA1[,6],type = 'l',ylab = expression(theta[6]))

### Inferencia Posterior ###

Theta1Post <- matrix(NA, nrow = n,ncol = 5)
for (i in 1:n) {
  Theta1Post[i,] <- c(mean(THETA1[,i]),sd(THETA1[,i]),quantile(THETA1[,i],c(0.025,0.5,0.975)))
}
ABPost <- matrix(NA, nrow = 2,ncol = 5)
for (j in 1:2) {
  ABPost[j,] <- c(mean(AB[,j]),sd(AB[,j]),quantile(AB[,j],c(0.025,0.5,0.975)))
}
Post1 <- rbind(ABPost,Theta1Post)
rownames(Post1) <- c('alpha','beta','theta_1','theta_2','theta_3','theta_4','theta_5','theta_6')
colnames(Post1) <- c('Media', 'Desviacion', '2.5%','50%', '97.5%')
Post1

xtable(Post1)

## Graficas de Theta_2 vs Theta_j ##
windows()
par(mfrow=c(2,3),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(x=THETA1[,2], y=THETA1[,1], xlab = expression(theta[2]), 
     ylab = expression(theta[1]), main = expression(paste(theta[2], 'vs', theta[1]))
     , pch = 20,xlim = c(0,0.78),ylim = c(0,0.78))
abline(b=1, a=0, col='blue',lwd=2)
plot(x=THETA1[,2], y=THETA1[,3], xlab = expression(theta[2]), 
     ylab = expression(theta[3]), main = expression(paste(theta[2], 'vs', theta[3]))
     , pch = 20,xlim = c(0,0.78),ylim = c(0,0.78))
abline(b=1, a=0, col='blue',lwd=2)
plot(x=THETA1[,2], y=THETA1[,4], xlab = expression(theta[2]), 
     ylab = expression(theta[4]), main = expression(paste(theta[2], 'vs', theta[4]))
     , pch = 20,xlim = c(0,0.78),ylim = c(0,0.78))
abline(b=1, a=0, col='blue',lwd=2)
plot(x=THETA1[,2], y=THETA1[,5], xlab = expression(theta[2]), 
     ylab = expression(theta[5]), main = expression(paste(theta[2], 'vs', theta[5]))
     , pch = 20,xlim = c(0,0.78),ylim = c(0,0.78))
abline(b=1, a=0, col='blue',lwd=2)
plot(x=THETA1[,2], y=THETA1[,6], xlab = expression(theta[2]), 
     ylab = expression(theta[6]), main = expression(paste(theta[2], 'vs', theta[6]))
     , pch = 20,xlim = c(0,0.78),ylim = c(0,0.78))
abline(b=1, a=0, col='blue',lwd=2)

# P(theta_2 > theta_j) #  

mean(THETA1[,2]>THETA1[,1])
mean(THETA1[,2]>THETA1[,3])
mean(THETA1[,2]>THETA1[,4])
mean(THETA1[,2]>THETA1[,5])
mean(THETA1[,2]>THETA1[,6])

maxTHETA1<-apply(THETA1, 1, max)
mean(THETA1[,2]==maxTHETA1) #p(theta_2 = max(theta))

y/x
## DIC del 1er modelo ##

lhat<-sum(dpois(x = y, lambda = colMeans(THETA1)*x,log = T))
pDIC1 <- 2*(lhat - mean(LLm1))
(DIC1 <- -2*lhat + 2*pDIC1 )

## Valores ppp del 1er modelo ##

B <- dim(THETA1)[1]
PP1 <- matrix(NA,nrow = B, ncol = 2)

set.seed(2409)
for (b in 1:B) {
  yp1 <- rpois(n = n, lambda = THETA1[b,]*x)  
  PP1[b,] <- c(mean(yp1),sd(yp1))
}

(PPmean <- mean(PP1[,1] > mean(y)))
(PPsd <- mean(PP1[,2] > sd(y)))

#### MODELO 2 ####
## Hiperparámetros ##
aalpha <- 1
balpha <- 1
abeta <- 10
bbeta <- 1 
alambda <- 0
blambda <- 10
acth<-rep(0,6)
aca<-0
acl<-0

## Algoritmo de MCMC del modelo 2 ##
S = 55000
THETA2 <- matrix(NA, nrow = S, ncol = n)
ABL <- matrix(NA, nrow = S, ncol = 3)
LLm2 <- NA
ncat <- floor(S/10) 

# Valores iniciales de la cadena
set.seed(16)
alpha<-rgamma(1,aalpha,balpha)
beta<-rgamma(1,abeta,bbeta)
lambda<-runif(1,alambda,blambda)
theta<- rgamma(n,alpha,beta)
deltath<-c(0.3, 0.4, 0.4, 0.3, 0.4, 0.4) 
deltaalp<-1.2 
deltalam<-4.5

#### Implementacion manual del MCMC del M2 ####

set.seed(2409)  # Establecer una semilla
for (s in 1:S) {
  # Actualizacion de los Thetas - Algoritmo de Metropolis #
  
  for (i in 1:n) {
    
    # Propuesta #
    
    theta.prop <- abs(runif(1,theta[i] - deltath[i], theta[i] + deltath[i]))
    theta.prop <- min(theta.prop, 2 - theta.prop)
    
    # Tasa de aceptacion en escala log
    
    log_r <- dnbinom(x= y[i], size = (theta.prop*x[i])/lambda, mu = theta.prop*x[i], log = T) +
             dgamma(theta.prop, shape = alpha, rate = beta, log=T)  -
             dnbinom(x= y[i], size = (theta[i]*x[i])/lambda, mu = theta[i]*x[i], log = T) -
             dgamma(theta[i], shape = alpha, rate = beta, log=T)
    # Actualizar thetas #
    
    if (log(runif(1)) < log_r) { 
      theta[i] <- theta.prop
      acth[i] <- acth[i] + 1
    }
    
  }
  
  # Actualizacion de los Alpha - Algoritmo de Metropolis #
  
  ### Propuesta ###
  alpha.prop <- abs(runif(1,alpha - deltaalp,alpha + deltaalp))
  
  ### Tasa de aceptacion en escala log ### 
  log_r <- sum(dgamma(theta, shape = alpha.prop, rate = beta, log = T)) +
    dgamma(alpha.prop, shape = aalpha, rate = balpha, log = T) - 
    sum(dgamma(theta, shape = alpha, rate = beta, log = T)) - 
    dgamma(alpha, shape = aalpha, rate = balpha, log = T)
  
  ### Actualizar alpha ###
  if (log(runif(1)) < log_r) { 
    alpha <- alpha.prop
    aca <- aca + 1
  }
  
  # Actualizacion de los Betas #
  beta <- rgamma(1, shape = abeta + n*alpha, rate = bbeta + sum(theta))
  
  # Actualizacion de los Lambda - Algoritmo de Metropolis #
  
  ### Propuesta ###
  lambda.prop <- abs(runif(1,lambda - deltalam,lambda + deltalam))
  
  ### Tasa de aceptacion en escala log ### 
  log_r <- sum(dnbinom(x= y, size = (theta*x)/lambda.prop, mu = theta*x, log = T)) +
           dunif(lambda.prop,alambda,blambda,log = T)-
           sum(dnbinom(x= y, size = (theta*x)/lambda, mu = theta*x, log = T)) -
           dunif(lambda,alambda,blambda,log = T)
  
  ### Actualizar lambda ###
  if (log(runif(1)) < log_r) { 
    lambda <- lambda.prop
    acl <- acl + 1
  }
  # Almacenamiento de los parametros #
  THETA2[s,] <- theta
  ABL[s,] <- c(alpha,beta,lambda)
  
  # Log-verosimilitud #
  LLm2[s]<-sum( dnbinom(x= y, size = (theta*x)/lambda, mu = theta*x, log = T))
  
  # Contador para el algoritmo #
  if (s%%ncat == 0) cat("El MCMC del 2do Modelo está ", 100*round(s/S, 1), "% completado \n", sep = "")
}

(Tasathetas <- acth/S) #Tasas de aceptación para los thetas
(Tasaalpha <- aca/S) #Tasas de aceptación para los thetas
(Tasalambda <- acl/S) #Tasas de aceptación para los thetas
ABL<-ABL[-c(1:5000),]
THETA2<-THETA2[-c(1:5000),]
# Tamaños efectivos de muestra
(ESSTHETA2 <- apply(THETA2,2,effectiveSize))
(ESSABL <- apply(ABL,2,effectiveSize))

# Coeficientes de Variacion MC
(SETHETA2 <- apply(THETA2,2,sd)/sqrt(ESSTHETA2))
(SEABL <- apply(ABL,2,sd)/sqrt(ESSABL))
(CVTHETA2 <- SETHETA2/colMeans(THETA2))
(CVABL <- SEABL/colMeans(ABL))

ESS2<-c(ESSABL,ESSTHETA2)
CV2<-c(CVABL,CVTHETA2)*100
Tab2 <- t(rbind(ESS2,CV2))
rownames(Tab2) <- c('alpha','beta','lambda','theta_1','theta_2','theta_3','theta_4','theta_5','theta_6')
colnames(Tab2) <- c('T.E. de Muestra', 'C.V. en%')
Tab2

xtable(Tab2)
# Diagramas de seguimiento
windows()
par(mfrow=c(3,3),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(ABL[,1],type = 'l',ylab = expression(alpha))
plot(ABL[,2],type = 'l',ylab = expression(beta))
plot(ABL[,3],type = 'l',ylab = expression(lambda))
plot(THETA2[,1],type = 'l',ylab = expression(theta[1]))
plot(THETA2[,2],type = 'l',ylab = expression(theta[2]))
plot(THETA2[,3],type = 'l',ylab = expression(theta[3]))
plot(THETA2[,4],type = 'l',ylab = expression(theta[4]))
plot(THETA2[,5],type = 'l',ylab = expression(theta[5]))
plot(THETA2[,6],type = 'l',ylab = expression(theta[6]))

### Inferencia Posterior ###

Theta2Post <- matrix(NA, nrow = n,ncol = 5)
for (i in 1:n) {
  Theta2Post[i,] <- c(mean(THETA2[,i]),sd(THETA2[,i]),quantile(THETA2[,i],c(0.025,0.5,0.975)))
}
ABLPost <- matrix(NA, nrow = 3,ncol = 5)
for (j in 1:3) {
  ABLPost[j,] <- c(mean(ABL[,j]),sd(ABL[,j]),quantile(ABL[,j],c(0.025,0.5,0.975)))
}
Post2 <- rbind(ABLPost,Theta2Post)
rownames(Post2) <- c('alpha','beta','lambda','theta_1','theta_2','theta_3','theta_4','theta_5','theta_6')
colnames(Post2) <- c('Media', 'Desviacion', '2.5%','50%', '97.5%')
Post2

xtable(Post2)

## Graficas de Theta_2 vs Theta_j ##
windows()
par(mfrow=c(2,3),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(x=THETA2[,2], y=THETA2[,1], xlab = expression(theta[2]), 
     ylab = expression(theta[1]), main = expression(paste(theta[2], 'vs', theta[1]))
     , pch = 20,xlim = c(0,0.78),ylim = c(0,0.78))
abline(b=1, a=0, col='red',lwd=2)
plot(x=THETA2[,2], y=THETA2[,3], xlab = expression(theta[2]), 
     ylab = expression(theta[3]), main = expression(paste(theta[2], 'vs', theta[3]))
     , pch = 20,xlim = c(0,0.78),ylim = c(0,0.78))
abline(b=1, a=0, col='red',lwd=2)
plot(x=THETA2[,2], y=THETA2[,4], xlab = expression(theta[2]), 
     ylab = expression(theta[4]), main = expression(paste(theta[2], 'vs', theta[4]))
     , pch = 20,xlim = c(0,0.78),ylim = c(0,0.78))
abline(b=1, a=0, col='red',lwd=2)
plot(x=THETA2[,2], y=THETA2[,5], xlab = expression(theta[2]), 
     ylab = expression(theta[5]), main = expression(paste(theta[2], 'vs', theta[5]))
     , pch = 20,xlim = c(0,0.78),ylim = c(0,0.78))
abline(b=1, a=0, col='red',lwd=2)
plot(x=THETA2[,2], y=THETA2[,6], xlab = expression(theta[2]), 
     ylab = expression(theta[6]), main = expression(paste(theta[2], 'vs', theta[6]))
     , pch = 20,xlim = c(0,0.78),ylim = c(0,0.78))
abline(b=1, a=0, col='red',lwd=2)

# P(theta_2 > theta_j) #  

mean(THETA2[,2]>THETA2[,1])
mean(THETA2[,2]>THETA2[,3])
mean(THETA2[,2]>THETA2[,4])
mean(THETA2[,2]>THETA2[,5])
mean(THETA2[,2]>THETA2[,6])

maxTHETA2<-apply(THETA2, 1, max)
mean(THETA2[,2]==maxTHETA2) #p(theta_2 = max(theta))

y/x
## DIC del 1er modelo ##
lambdahat<-mean(ABL[,3])
muhat<-colMeans(THETA2)*x
lhat2<-sum(dnbinom(x = y, size = muhat/lambdahat,mu= muhat,log = T))
pDIC2 <- 2*(lhat2 - mean(LLm2))
(DIC2 <- -2*lhat2 + 2*pDIC2 )

## Valores ppp del 1er modelo ##

B <- dim(THETA2)[1]
PP2 <- matrix(NA,nrow = B, ncol = 2)

set.seed(2409)
for (b in 1:B) {
  m <-THETA2[b,]*x
  s <-m/ABL[b,3]
  yp2 <- rnbinom(n = n, size = s, mu = m)  
  PP2[b,] <- c(mean(yp2),sd(yp2))
}

(PP2mean <- mean(PP2[,1] > mean(y)))
(PP2sd <- mean(PP2[,2] > sd(y)))
