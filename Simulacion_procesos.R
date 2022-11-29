# -----------------------------------------------------------------------------#
# Programa para la simulación de procesos estocásticos lineales ----------------
## Arturo Yesid Gonzalez ----
# -----------------------------------------------------------------------------#

# limpieza del ambiente ----
rm(list = ls())


library(pacman)
p_load(car, forecast, TSstudio, timeSeries, TSA, ggplot2, glue)

set.seed(1234)

salida_proc_ar <- "Procesos_AR_p/"
salida_proc_ma <- "Procesos_MA_q/"

# Proceso AR(1) ----
# -----------------------------------------------------------------------------#
getwd()


# multiples graficos de proceso AR(1)
phi_sec = seq(from = -0.94, to=0.94, 0.05)

# Ciclo para generar graficos
for (i in 1:length(phi_sec)) {
  nombre <- paste("Grafico_Proc_AR(1)","_phi=",phi_sec[i],".png", sep ="")
  
  png(glue(getwd(),"/","{salida_proc_ar}{nombre}")) # Generacion de nombre archivo con apuntador a carpeta
  AR_test <- arima.sim(list(order(1,0,0), ar=phi_sec[i]), sd = sqrt(0.5) 
                       , n=500) # ARIMA(p, 0, 0) = AR(p) -- AR(1) sigma = 0.5
  par(mfrow=c(1,3), las=1)
  plot(as.ts(AR_test), main="serie")
  acf(AR_test,lag.max = 18, main="FAC simple", ylab="", xlab="Rezago")
  pacf(AR_test, lag.max = 18, main="FACP simple", ylab="", xlab="Rezago")
  dev.off()
  
}

AR2 <- arima.sim(list(order(2,0,0), ar= c(0.5,0.3)), sd= sqrt(5), n=5000)

grafico1 <- autoplot(as.ts(AR2)) + ylab("Serie") + xlab("tiempo")
grafico2 <- Acf(as.ts(AR2), plot = F) %>% autoplot() + ylab("ACF - FAC") + xlab("rezago")
grafico3 <- Pacf(as.ts(AR2), plot = F) %>% autoplot() + ylab("PACF - FACP") + xlab("rezago")
windows()
gridExtra::grid.arrange(grafico1,grafico2,grafico3, nrow=1)


AR3 <- arima.sim(list(order(3,0,0), ar= c(0.5,-0.2,0.4)), sd= sqrt(5), n=5000)

grafico1 <- autoplot(as.ts(AR3)) + ylab("Serie") + xlab("tiempo")
grafico2 <- Acf(as.ts(AR3), plot = F) %>% autoplot() + ylab("ACF - FAC") + xlab("rezago")
grafico3 <- Pacf(as.ts(AR3), plot = F) %>% autoplot() + ylab("PACF - FACP") + xlab("rezago")
windows()
gridExtra::grid.arrange(grafico1,grafico2,grafico3, nrow=1)


# Mensaje Proc. AR(p) ----
# La caracterización de procesos AR(p): cortan la PACF en el nivel p y FAC
# presenta un comportamiento de caida exponencial
# ---------------------------------------------------------------------------#


# Proceso MA(1) ----
# -----------------------------------------------------------------------------#

# multiples graficos de proceso MA(1)
theta_sec = seq(from = -0.94, to=0.94, 0.05)

# Ciclo for

for (i in 1:length(theta_sec)) {
  nombre <- paste("Grafico_Proc_MA(1)","_theta=",theta_sec[i],".png", sep ="")
  
  png(glue(getwd(),"/","{salida_proc_ma}{nombre}"))
  MA_test <- arima.sim(list(order(0,0,1), ma=theta_sec[i]), sd = sqrt(0.5) 
                       , n=500) # ARIMA(0, 0, q) = MA(q) -- MA(1) sigma = 0.5
  par(mfrow=c(1,3), las=1)
  plot(as.ts(MA_test), main="serie")
  acf(MA_test,lag.max = 18, main="FAC simple", ylab="", xlab="Rezago")
  pacf(MA_test, lag.max = 18, main="FACP", ylab="", xlab="Rezago")
  dev.off()
  
}

MA2 <- arima.sim(list(order(0,0,2), ma= c(0.4,0.3)), sd= sqrt(5), n=5000)

grafico1 <- autoplot(as.ts(MA2)) + ylab("Serie") + xlab("tiempo")
grafico2 <- Acf(as.ts(MA2), plot = F) %>% autoplot() + ylab("ACF - FAC") + xlab("rezago")
grafico3 <- Pacf(as.ts(MA2), plot = F) %>% autoplot() + ylab("PACF - FACP") + xlab("rezago")
windows()
gridExtra::grid.arrange(grafico1,grafico2,grafico3, nrow=1)


MA3 <- arima.sim(list(order(0,0,3), ma= c(-0.5,0.2,-0.2)), sd= sqrt(5), n=5000)

grafico1 <- autoplot(as.ts(MA3)) + ylab("Serie") + xlab("tiempo")
grafico2 <- Acf(as.ts(MA3), plot = F) %>% autoplot() + ylab("ACF - FAC") + xlab("rezago")
grafico3 <- Pacf(as.ts(MA3), plot = F) %>% autoplot() + ylab("PACF - FACP") + xlab("rezago")
windows()
gridExtra::grid.arrange(grafico1,grafico2,grafico3, nrow=1)



# Mensaje Proc. MA(q) ----
# La caracterización de procesos MA(q): cortan la ACF en el nivel q y PACF
# presenta un comportamiento aleatorio, en casi todos los casos
# ---------------------------------------------------------------------------#



# Simulación procesos ARMA

ARMA11 <- arima.sim(list(order(1,0,1), ar=c(-0.7), ma= c(0.4)), sd= sqrt(5), n=5000)

grafico1 <- autoplot(as.ts(ARMA11)) + ylab("Serie") + xlab("tiempo")
grafico2 <- Acf(as.ts(ARMA11), plot = F) %>% autoplot() + ylab("ACF - FAC") + xlab("rezago")
grafico3 <- Pacf(as.ts(ARMA11), plot = F) %>% autoplot() + ylab("PACF - FACP") + xlab("rezago")
windows()
gridExtra::grid.arrange(grafico1,grafico2,grafico3, nrow=1)


ARMA22 <- arima.sim(list(order(2,0,2), ar=c(-0.7,0.2), ma= c(0.4,0.3)), sd= sqrt(5), n=5000)

grafico1 <- autoplot(as.ts(ARMA22)) + ylab("Serie") + xlab("tiempo")
grafico2 <- Acf(as.ts(ARMA22), plot = F) %>% autoplot() + ylab("ACF - FAC") + xlab("rezago")
grafico3 <- Pacf(as.ts(ARMA22), plot = F) %>% autoplot() + ylab("PACF - FACP") + xlab("rezago")
windows()
gridExtra::grid.arrange(grafico1,grafico2,grafico3, nrow=1)


# Mensaje Proc. ARMA(p,q) ----
# La caracterización de procesos ARMA(p,q) a partir de los gráficos ACF y PACF
# no es intuitivo
# ---------------------------------------------------------------------------#





# ---------------------------------------------------------------------------#
# Simulación caminatas aleatorias ----
# ---------------------------------------------------------------------------#

# Creación de variable tiempo de simulación
tiempo <- 1000

# inicialización de los procesos a simular
y <- a <- rnorm(n=tiempo, sd=1, mean = 0)


# Ciclo for para la creación de la serie " proc. estocástico"

for (t in 2:tiempo) {
  y[t] <- y[t-1] + a[t]
}

grafico1 <- autoplot(as.ts(y), main("serie")) + ylab("Serie") + xlab("tiempo")
grafico2 <- Acf(as.ts(y), plot = F) %>% autoplot() + ylab("ACF - FAC") + xlab("rezago")
grafico3 <- Pacf(as.ts(y), plot = F) %>% autoplot() + ylab("PACF - FACP") + xlab("rezago")
windows()
gridExtra::grid.arrange(grafico1,grafico2,grafico3, nrow=1)



# ---------------------------------------------------------------------------
# Simulación caminata aleatoria con deriva "drift"
# ---------------------------------------------------------------------------
nu <- -0.3
beta <- 0.0005

for (t in 2:tiempo) {
  y[t] <- nu + beta*t + y[t-1] + a[t]
}

grafico1 <- autoplot(as.ts(y)) + ylab("") + xlab("tiempo")
grafico2 <- Acf(as.ts(y),plot = F) %>% autoplot() + ylab("") + xlab("rezago")
grafico3 <- Pacf(as.ts(y),plot = F) %>% autoplot() + ylab("") + xlab("rezago")

windows()
gridExtra::grid.arrange(grafico1, grafico2, grafico3, nrow=1)

