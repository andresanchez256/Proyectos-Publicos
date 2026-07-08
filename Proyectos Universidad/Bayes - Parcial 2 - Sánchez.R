#####--- Parcial No. 2 , Estadistica Bayesiana ---####

# Andrés Felipe Sánchez Felipe #

## Librerias utilizadas ##

library(coda)
library(xtable)
library(sf)
library(dplyr)
library(ggplot2)
library(stringr)

###        CARGA DE DATOS       ###
## Cargar datos de Entrenamiento ##

datosE <- read.csv2("C:/Users/ANDRES/Trabajos Universidad/Estadistica/Bayesiana/SB11_1.txt")
  
## Cargar datos de Prueba ##
  
datosP <- read.csv2("C:/Users/ANDRES/Trabajos Universidad/Estadistica/Bayesiana/SB11_2.txt")

##### Punto 1 #####

# Medias #

mediasE <- aggregate(punt_mate~cole_dept,datosE,mean)
mediasE$cole_dept[mediasEnt$cole_dept=='5']<-'05'
mediasE$cole_dept[mediasEnt$cole_dept=='8']<-'08'

## Mapa departamentos

DepMAP<- sf::st_read("C:/Users/ANDRES/Trabajos Universidad/Estadistica/Bayesiana/MGN_DPTO_POLITICO.shp",quiet=TRUE)

MAPCOL <- DepMAP %>% left_join(mediasE,by=c("DPTO_CCDGO"='cole_dept'))

windows()
ggplot(MAPCOL) +
  geom_sf(data=MAPCOL,aes(fill=punt_mate),col='white',linetype='solid') +
  coord_sf(xlim=c(-78.9,-66.4),ylim=c(-3.8,11.9)) + 
  labs(x='',y='',fill="Puntaje Promedio") + 
  ggtitle('PUNTAJES PROMEDIO EN MATEMATICAS POR DEPARTAMENTO') + 
  scale_fill_gradient(low="darkgreen",high="cyan",breaks=seq(0,100,5)) + 
  theme(legend.position="right",
        title=element_text(family = 'italic',color = 'royalblue',face = 'bold',size=26,vjust = 1,
                           hjust = 0.5),
        legend.text=element_text(face="mono",size=15,hjust=1,vjust=0.5,
                                 color="black",angle=0),
        legend.direction="vertical") + theme_bw()

##### Punto 3 #####
# MODELO 1 : MODELO NORMAL SEMI-CONJUGADO #

y <- datosE$punt_mate

# estadisticos
N      <- length(y)
mean.y <- mean(y)
var.y  <- var(y)
sum.y  <- N*mean(y)
ss.y   <- (N-1)*var(y)

# hiperparametros

mu0 <- 50
t20 <- 25
s20 <- 100
nu0 <- 1

# setup

S    <- 25000                            # numero de muestras 
MS_1   <- matrix(data=NA, nrow=S, ncol=2)  # matriz para almacenar las muestras
LP1  <- matrix(data=NA, nrow=S, ncol=1)
ncat <- floor(S/10)                      # mostrar anuncios cada 10%

# inicializar la cadena

mu   <- mean.y
sig2 <- var.y

# Algoritmo MCMC
set.seed(2409)
for(s in 1:S) {
  # actualizar el valor de theta
  t2n <- 1/(1/t20 + N/sig2)
  mun <- t2n*(mu0/t20 + sum.y/sig2)
  mu  <- rnorm(1, mun, sqrt(t2n))
  # actualizar el valor de sigma^2
  an  <- (nu0 + N)/2
  bn  <- (nu0*s20 + ss.y + N*(mean.y - mu)^2)/2
  sig2 <- 1/rgamma(1, an, bn)
  # almacenar
  MS_1[s,] <- c(mu, sig2)
  # log-verosimilitud
  LP1[s] <- sum(dnorm(x = y, mean = mu, sd = sqrt(sig2), log = T))
  # progreso
  if (s%%ncat == 0) cat("Algoritmo MCMC modelo 1, ", 100*round(s/S, 1), "% completado \n", sep = "")
}
### Diagnostico del Modelo ###

plot(LP1,type = 'l')
effectiveSize(MS_1)

### Retirar muestras de calentamiento ###
MS_1<-MS_1[5001:S,]
LP1<-LP1[5001:S]

### Tamaño efectivo de muestra ###

effectiveSize(MS_1)

### Tabla de reporte ###

Theta   <- cbind(mean(MS_1[,1]),sd(MS_1[,1]),quantile(MS_1[,1],probs = 0.025),quantile(MS_1[,1],probs = 0.975))
Sigma   <- cbind(mean(sqrt(MS_1[,2])),sd(sqrt(MS_1[,2])),quantile(sqrt(MS_1[,2]),probs = 0.025),quantile(sqrt(MS_1[,2]),probs = 0.975))
tablaM1 <- rbind(Theta,Sigma)
colnames(tablaM1) <- c('Media','Desv. Estandar', 'Perc. 2.5', 'Perc. 97.5')
rownames(tablaM1) <- c(expression(theta),expression(sigma))

xtable(tablaM1) # Exportar a formato Latex

##### Punto 4 #####
# MODELO 2 : MODELO JERARQUICO PARA MEDIA #

## Grupos ##

g<-factor(datosE$cole_dept)
niveles<-levels(g)
m<-length(niveles)

### Organizar los Y por departamentoss###

Y<-list()
for (j in 1:m) {
  Y[[j]]<-datosE[datosE[,'cole_dept']==niveles[j],'punt_mate']
}

## Estadísticos ##
n <- NULL
ybar <- NULL
s2 <- NULL

for (j in 1:m){
  n[j]<-length(Y[[j]])
  ybar[j]<-mean(Y[[j]])
  s2[j]<-var(Y[[j]])
}


### hiperparámetros para el modelo 2 ###

mu0  <- 50 ;g20  <- 25
eta0 <- 1  ;t20  <- 100
nu0  <- 1  ;s20  <- 100

# setup
S     <- 25000                            # numero de iteraciones
ncat  <- floor(S/10)                      # progreso
THETA_2 <- matrix(data=NA, nrow=S, ncol=m)  # almacenar thetas 
MST_2   <- matrix(data=NA, nrow=S, ncol=3)  # almacenar mu, sigma^2, tau^2
LP2   <- matrix(data=NA, nrow=S, ncol=1)

### valores inciales ###
theta  <- ybar
sigma2 <- mean(s2)
mu     <- mean(theta)
tau2   <- var(theta)

### Algoritmo MCMC ###
set.seed(2409)
for (s in 1:S) {
  # actualizar theta
  for(j in 1:m) {
    vtheta   <- 1/(n[j]/sigma2+1/tau2)
    etheta   <- vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
    theta[j] <- rnorm(1,etheta,sqrt(vtheta))
  }
  # actualizar sigma^2
  nun <- nu0+sum(n)
  ss  <- nu0*s20
  for (j in 1:m) {
    ss <- ss+sum((Y[[j]]-theta[j])^2)
  }
  sigma2 <- 1/rgamma(1,nun/2,ss/2)
  # actualizar mu
  vmu <- 1/(m/tau2+1/g20)
  emu <- vmu*(m*mean(theta)/tau2 + mu0/g20)
  mu  <- rnorm(1,emu,sqrt(vmu)) 
  # actualizar tau2
  etam <- eta0+m
  ss   <- eta0*t20 + sum( (theta-mu)^2 )
  tau2 <- 1/rgamma(1,etam/2,ss/2)
  # almacenar valores
  THETA_2[s,] <- theta
  MST_2[s,]   <- c(mu,sigma2,tau2)
  # log-verosimilitud (opcion 2)
  LP2[s] <- sum(dnorm(x = y, mean = rep(theta, n), sd = sqrt(sigma2), log = T))
  # progreso
  if (s%%ncat == 0) cat("Algoritmo MCMC modelo 2, ", 100*round(s/S, 1), "% completado \n", sep = "")
} 

### Diagnostico del Modelo ###

plot(LP2,type = 'l')
effectiveSize(MST_2)
effectiveSize(THETA_2)


### Retirar muestras de calentamiento ###

MST_2<-MST_2[5001:S,]
THETA_2<-THETA_2[5001:S,]
LP2<-LP2[5001:S]

### Tamaños efectivos de muestra ###

effectiveSize(MST_2)
effectiveSize(THETA_2)

### Tabla de reporte ###

Mu    <- cbind(mean(MST_2[,1]),sd(MST_2[,1]),quantile(MST_2[,1],probs = 0.025),quantile(MST_2[,1],probs = 0.975))
Tau   <- cbind(mean(sqrt(MST_2[,3])),sd(sqrt(MST_2[,3])),quantile(sqrt(MST_2[,3]),probs = 0.025),quantile(sqrt(MST_2[,3]),probs = 0.975))
Sigma <- cbind(mean(sqrt(MST_2[,2])),sd(sqrt(MST_2[,2])),quantile(sqrt(MST_2[,2]),probs = 0.025),quantile(sqrt(MST_2[,2]),probs = 0.975))
tablaM2 <- rbind(Mu,Tau,Sigma)
colnames(tablaM2) <- c('Media','Desv. Estandar', 'Perc. 2.5', 'Perc. 97.5')
rownames(tablaM2) <- c(expression(mu),expression(tau),expression(sigma))

xtable(tablaM2) # Exportar a formato Latex

##### Punto 5 #####
# MODELO 3 : MODELO JERARQUICO PARA MEDIA Y VARIANZA#

### hiperparámetros para el modelo 2 ###

mu0  <- 50 ; g20  <- 25
eta0 <- 1  ; t20  <- 100
a0  <- 1   ; b0  <- 1/100; kappa0 <- 1  

# setup
S     <- 25000                             # numero de iteraciones
ncat  <- floor(S/10)                       # progreso
THETA_3  <- matrix(data=NA, nrow=S, ncol=m)
SIGMA2_3 <- matrix(data=NA, nrow=S, ncol=m)
MTSN_3   <- matrix(data=NA, nrow=S, ncol=4)  # mu, tau2, sigma02, nu0
LP3    <- matrix(data=NA, nrow=S, ncol=1)
nu0s   <- 1:100                            # rango de valores para muestrear en p(nu_0 | rest)

### valores inciales ###
theta  <- ybar
sigma2 <- s2 
mu     <- mean(theta)
tau2   <- var(theta)
s20    <- 100
nu0    <- 1

### Algoritmo MCMC ###
set.seed(2409)
for(s in 1:S) {
  # actualizar thetas
  for(j in 1:m) {
    vtheta   <- 1/(n[j]/sigma2[j]+1/tau2)
    etheta   <- vtheta*(ybar[j]*n[j]/sigma2[j]+mu/tau2)
    theta[j] <- rnorm(1,etheta,sqrt(vtheta))
  }
  # actualizar sigma2s
  for(j in 1:m) { 
    nun       <- nu0+n[j]
    ss        <- nu0*s20+ sum((Y[[j]]-theta[j])^2)
    sigma2[j] <- 1/rgamma(1,nun/2,ss/2)
  }
  # actualizar s20
  s20 <- rgamma(1,a0+m*nu0/2,b0+nu0*sum(1/sigma2)/2)
  # actualizar nu0
  lpnu0 <- .5*nu0s*m*log(s20*nu0s/2)-m*lgamma(nu0s/2)+(nu0s/2)*sum(log(1/sigma2)) - nu0s*s20*sum(1/sigma2)/2 - kappa0*nu0s
  nu0   <- sample(x = nu0s, size = 1, prob=exp( lpnu0-max(lpnu0) ))
  # actualizar mu
  vmu <- 1/(m/tau2+1/g20)
  emu <- vmu*(m*mean(theta)/tau2 + mu0/g20)
  mu  <- rnorm(1,emu,sqrt(vmu))
  # actualizar tau2
  etam <-eta0+m
  ss   <- eta0*t20 + sum( (theta-mu)^2 )
  tau2 <-1/rgamma(1,etam/2,ss/2)
  # almacenar
  THETA_3[s,]  <- theta
  SIGMA2_3[s,] <- sigma2
  MTSN_3[s,]   <- c(mu,tau2,s20,nu0)
  # log-verosimilitud (opcion 2)
  LP3[s] <- sum(dnorm(x = y, mean = rep(theta, n), sd = sqrt(rep(sigma2, n)), log = T))
  ### progreso
  if (s%%ncat == 0) cat("Algoritmo MCMC modelo 3, ", 100*round(s/S, 1), "% completado \n", sep = "")
}

### Diagnostico del Modelo ###

plot(LP3,type = 'l')
effectiveSize(MTSN_3)
effectiveSize(THETA_3)
effectiveSize(SIGMA2_3)


### Retirar muestras de calentamiento ###

MTSN_3<-MTSN_3[5001:S,]
THETA_3<-THETA_3[5001:S,]
SIGMA2_3<-SIGMA2_3[5001:S,]
LP3<-LP3[5001:S]

### Tamaños efectivos de muestra ###

effectiveSize(MTSN_3)
effectiveSize(THETA_3)
effectiveSize(SIGMA2_3)

### Tabla de reporte ###

Mu    <- cbind(mean(MTSN_3[,1]),sd(MTSN_3[,1]),quantile(MTSN_3[,1],probs = 0.025),quantile(MTSN_3[,1],probs = 0.975))
Tau   <- cbind(mean(sqrt(MTSN_3[,2])),sd(sqrt(MTSN_3[,2])),quantile(sqrt(MTSN_3[,2]),probs = 0.025),quantile(sqrt(MTSN_3[,2]),probs = 0.975))
Sigma <- cbind(mean(sqrt(MTSN_3[,3])),sd(sqrt(MTSN_3[,3])),quantile(sqrt(MTSN_3[,3]),probs = 0.025),quantile(sqrt(MTSN_3[,3]),probs = 0.975))
tablaM3 <- rbind(Mu,Tau,Sigma)
colnames(tablaM3) <- c('Media','Desv. Estandar', 'Perc. 2.5', 'Perc. 97.5')
rownames(tablaM3) <- c(expression(mu),expression(tau),expression(sigma))

xtable(tablaM3) # Exportar a formato Latex

# Maximo de nu posterior#
freq <- data.frame(table(MTSN_3[,4]))
(max.post_nu <- freq[which.max(freq$Freq),1])

##### Punto 6 #####
### DIC ###
# DIC del Modelo 1 #
mu_hat   <- mean(MS_1[,1])
sig2_hat <- mean(MS_1[,2])
lp_hat.m1 <- sum(dnorm(x = y, mean = mu_hat, sd = sqrt(sig2_hat), log = T))
pDIC.m1  <- 2*(lp_hat.m1 - mean(LP1))
(DIC.m1   <- -2*lp_hat.m1 + 2*pDIC.m1)

# DIC del Modelo 2 #
theta_hat  <- colMeans(THETA_2)
sigma2_hat <- colMeans(MST_2)[2]
lp_hat.m2 <- sum(dnorm(x = y, mean = rep(theta_hat, n), sd = sqrt(sigma2_hat), log = T))
pDIC.m2  <- 2*(lp_hat.m2 - mean(LP2))
(DIC.m2   <- -2*lp_hat.m2 + 2*pDIC.m2)

# DIC del Modelo 3 #
theta_hat  <- colMeans(THETA_3)
sigma2_hat <- colMeans(SIGMA2_3)
lp_hat.m3 <- sum(dnorm(x = y, mean = rep(theta_hat, n), sd = sqrt(rep(sigma2_hat, n)), log = T))
pDIC.m3  <- 2*(lp_hat.m3 - mean(LP3))
(DIC.m3   <- -2*lp_hat.m3 + 2*pDIC.m3)

### WAIC ###
# WAIC del Modelo 1 #
lppd.m1  <- 0
pWAIC.m1 <- 0
for (i in 1:N) {
  # lppd
  tmp1    <- dnorm(x = y[i], mean = MS_1[,1], sd = sqrt(MS_1[,2]))
  lppd.m1 <- lppd.m1 + log(mean(tmp1))
  # pWAIC
  tmp2 <- dnorm(x = y[i], mean = MS_1[,1], sd = sqrt(MS_1[,2]), log = T)
  pWAIC.m1 <- pWAIC.m1 + 2*( log(mean(tmp1)) - mean(tmp2) )
}
(WAIC.m1 <- -2*lppd.m1 + 2*pWAIC.m1)

# WAIC del Modelo 2 #
lppd.m2  <- 0
pWAIC.m2 <- 0
for (i in 1:N) {
  # lppd
  tmp1    <- dnorm(x = y[i], mean = THETA_2[,g[i]], sd = sqrt(MST_2[,2]))
  lppd.m2 <- lppd.m2 + log(mean(tmp1))
  # pWAIC
  tmp2 <- dnorm(x = y[i], mean = THETA_2[,g[i]], sd = sqrt(MST_2[,2]), log = T)
  pWAIC.m2 <- pWAIC.m2 + 2*( log(mean(tmp1)) - mean(tmp2) )
}
(WAIC.m2 <- -2*lppd.m2 + 2*pWAIC.m2)

# WAIC del Modelo 3 #
lppd.m3  <- 0
pWAIC.m3 <- 0
for (i in 1:N) {
  # lppd
  tmp1    <- dnorm(x = y[i], mean = THETA_3[,g[i]], sd = sqrt(SIGMA2_3[,g[i]]))
  lppd.m3 <- lppd.m3 + log(mean(tmp1))
  # pWAIC
  tmp2 <- dnorm(x = y[i], mean = THETA_3[,g[i]], sd = sqrt(SIGMA2_3[,g[i]]), log = T)
  pWAIC.m3 <- pWAIC.m3 + 2*( log(mean(tmp1)) - mean(tmp2) )
}
(WAIC.m3 <- -2*lppd.m3 + 2*pWAIC.m3)

### MSE ###
### Datos de Prueba ###
# Se ordenan los datos segun el factor #

y_test <- datosP$punt_mate
N_test <- length(y_test)
n_test <- NULL
Y.test <- list()
for (i in 1:m) { 
  Y.test[[i]] <-datosP[datosP[,'cole_dept']==niveles[i],'punt_mate']
  n_test[i]<-length(Y.test[[i]])
  
}

# MSE del Modelo 1 #
sum.YP1 <- c(rep(0,N_test))
set.seed(2409)
for (b in 1:B) {
  yrep <- rnorm(n = N_test, mean = MS_1[b,1], sd = sqrt(MS_1[b,2]))
  sum.YP1 <- yrep + sum.YP1
  # progreso
  if (b%%ncat == 0) cat("MSE del modelo 1, ", 100*round(b/B, 1), "% completado \n", sep = "")
}
mean.YP1<-sum.YP1/B
(MSE_m1 <- mean((y_test - mean.YP1)^2))

# MSE del Modelo 2 #
sum.YP2 <- c(rep(0,N_test))
set.seed(2409)
for (b in 1:B) {
  theta   <- THETA_2[b,]
  sigma2  <- MST_2[b,2]
  yrep    <- rnorm(n = N_test, mean = rep(theta, n_test), sd = sqrt(sigma2))
  sum.YP2 <- yrep + sum.YP2
  # progreso
  if (b%%ncat == 0) cat("MSE del modelo 2, ", 100*round(b/B, 1), "% completado \n", sep = "")
}
mean.YP2<-sum.YP2/B
(MSE_m2 <- mean((y_test - mean.YP2)^2))

# MSE del Modelo 3 #
sum.YP3 <- c(rep(0,N_test))
set.seed(2409)
for (b in 1:B) {
  theta   <- THETA_3[b,]
  sigma2  <- SIGMA2_3[b,]
  yrep    <- rnorm(n = N_test, mean = rep(theta, n_test), sd = sqrt(rep(sigma2, n_test)))
  sum.YP3 <- yrep + sum.YP3
  # progreso
  if (b%%ncat == 0) cat("MSE del modelo 3, ", 100*round(b/B, 1), "% completado \n", sep = "")
}
mean.YP3<-sum.YP3/B
(MSE_m3 <- mean((y_test - mean.YP3)^2))

### Tabla de reporte ###

DIC <- c(DIC.m1,DIC.m2,DIC.m3)
WAIC <-c(WAIC.m1,WAIC.m2,WAIC.m3)
MSE <-c(MSE_m1,MSE_m2,MSE_m3)
CompMod <- rbind(DIC,WAIC,MSE)
colnames(CompMod) <- c('Modelo 1','Modelo 2', 'Modelo 3')
rownames(CompMod) <- c('DIC','WAIC','MSE')

xtable(CompMod)

##### Punto 7 #####
### Intervalos de Credibilidad ###
that <- colMeans(THETA_3)
ic1  <- apply(X = THETA_3, MARGIN = 2, FUN = function(x) quantile(x, c(0.025,0.975)))
ic2  <- apply(X = THETA_3, MARGIN = 2, FUN = function(x) quantile(x, c(0.005,0.995)))
departamentos <- c('ANTIOQUIA','ATLANTICO','BOGOTA','BOLIVAR','BOYACA','CALDAS','CAQUETA','CAUCA','CESAR',
                   'CORDOBA','CUNDINAMARCA','HUILA','LA GUAJIRA','MAGDALENA','META','NARIÑO','N. DE SANTANDER',
                   'QUINDIO','RISARALDA','SANTANDER','TOLIMA','V. DEL CAUCA','ARAUCA','CASANARE','PUTUMAYO')

ranking <- order(that, decreasing = T)
that <- that[ranking]
ic1  <- ic1[,ranking]
ic2  <- ic2[,ranking]
departamentos <- departamentos[ranking] 

k    <- m
that <- that[1:k]
ic1  <- ic1[,1:k]
ic2  <- ic2[,1:k]
departamentos <- departamentos[1:k] 

colo <- c("green","darkgreen")[as.numeric((ic2[1,] < 50) & (50 < ic2[2,]))+1]

# grafico
windows()
par(mfrow=c(3,1),mar=c(3,3,1.5,1),mgp=c(1.75,.75,0))
plot(NA, NA, xlab = "Departamento", ylab = "Puntaje", main = paste0("Intervalos de credibilidad para el puntaje promedio de los departamentos"), 
     xlim = c(1,8), ylim = c(49.9,75), cex.axis = 0.75, xaxt = "n")
axis(side = 1, at = 1:8, labels = departamentos[1:8], cex.axis = 0.75)
abline(h = 50, col = "gray", lwd = 2)
for (j in 1:8) {
  segments(x0 = j, y0 = ic1[1,j], x1 = j, y1 = ic1[2,j], lwd = 3, col = colo[j])
  segments(x0 = j, y0 = ic2[1,j], x1 = j, y1 = ic2[2,j], lwd = 1, col = colo[j])
  lines(x = j, y = that[j], type = "p", pch = 16, cex = 0.8, col = colo[j])
}
plot(NA, NA, xlab = "Departamento", ylab = "Puntaje", 
     xlim = c(9,16), ylim = c(48,62.5), cex.axis = 0.75, xaxt = "n")
axis(side = 1, at = 9:16, labels = departamentos[9:16], cex.axis = 0.75)
abline(h = 50, col = "gray", lwd = 2)
for (j in 9:16) {
  segments(x0 = j, y0 = ic1[1,j], x1 = j, y1 = ic1[2,j], lwd = 3, col = colo[j])
  segments(x0 = j, y0 = ic2[1,j], x1 = j, y1 = ic2[2,j], lwd = 1, col = colo[j])
  lines(x = j, y = that[j], type = "p", pch = 16, cex = 0.8, col = colo[j])
}
plot(NA, NA, xlab = "Departamento", ylab = "Puntaje", 
     xlim = c(17,k), ylim = c(30,60), cex.axis = 0.75, xaxt = "n")
axis(side = 1, at = 17:k, labels = departamentos[17:k], cex.axis = 0.75)
abline(h = 50, col = "gray", lwd = 2)
for (j in 17:k) {
  segments(x0 = j, y0 = ic1[1,j], x1 = j, y1 = ic1[2,j], lwd = 3, col = colo[j])
  segments(x0 = j, y0 = ic2[1,j], x1 = j, y1 = ic2[2,j], lwd = 1, col = colo[j])
  lines(x = j, y = that[j], type = "p", pch = 16, cex = 0.8, col = colo[j])
}

##### Punto 8 #####
### Puntajes de matematicas en Casanare ###
Casanare <- datosE[datosE[,'cole_dept']==85,2]

# estadisticos de prueba

(ts_obs <- mean(Casanare))
n_Casn <- length(Casanare)

B = 20000

TS_Cas <- matrix(data=NA, nrow=B, ncol=3)
ncat <- floor(B/10)

## Se establece una semilla ##
set.seed(2409)  
for(b in 1:B){
  ## Predictiva posterior del modelo 1 ##
  
  theta.m1   <- MS_1[b,1]
  sigma2.m1  <- MS_1[b,2]
  yp_cas.m1  <- rnorm(n = n_Casn, mean = theta.m1, sd = sqrt(sigma2.m1))
  
  ## Predictiva posterior del modelo 2 ##
  
  theta.m2   <- THETA_2[b,24]
  sigma2.m2  <- MST_2[b,2]
  yp_cas.m2  <- rnorm(n = n_Casn, mean = theta.m2, sd = sqrt(sigma2.m2))
  
  ## Predictiva posterior del modelo 3 ##
  
  theta.m3   <- THETA_3[b,24]
  sigma2.m3  <- SIGMA2_3[b,24]
  yp_cas.m3  <- rnorm(n = n_Casn, mean = theta.m3, sd = sqrt(sigma2.m3))
  TS_Cas[b,] <- as.numeric(c(mean(yp_cas.m1),mean(yp_cas.m2),mean(yp_cas.m3)))
  # progreso
  if (b%%ncat == 0) cat("PP para el departamento del Casanare, ", 100*round(b/B, 1), "% completado \n", sep = "")
}
### ppp (posterior predictive p-values)

(ppp_cas.m1<-mean(TS_Cas[,1]>ts_obs))
(ppp_cas.m2<-mean(TS_Cas[,2]>ts_obs))
(ppp_cas.m3<-mean(TS_Cas[,3]>ts_obs))

### Tabla de reporte ###

ppp_cas <- c(ppp_cas.m1,ppp_cas.m2,ppp_cas.m3)
PPPCAS <- rbind(ppp_cas)
colnames(PPPCAS) <- c('Modelo 1','Modelo 2', 'Modelo 3')
rownames(PPPCAS) <- c('Valor ppp')

xtable(PPPCAS)
# Gráficos #
par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(x = NA, y = NA, ylab = "Densidad", xlab = "Media predictiva posterior", 
     cex.axis = 0.7, xlim = c(22,77), ylim = c(0.0031, 0.079)) 
lines(density(TS_Cas[,1], adjust = 2.2), col = "purple")
lines(density(TS_Cas[,2], adjust = 2.2), col = "blue")
lines(density(TS_Cas[,3], adjust = 2.2), col = "green")
abline(v = ts_obs, lty = 2, col = "red")
legend("topleft", legend = c("Media muestral","Modelo 1","Modelo 2","Modelo 3"), col = c("red","purple","blue","green"), lwd = 2, bty = "n")







