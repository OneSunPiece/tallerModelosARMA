# Visuals
library(ggplot2)
# Stats
library(stats)
# Periodograma
library(TSA)
# forecast
library(forecast)
# dplyr
library(dplyr)

#rm(list = ls(all = TRUE))

source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Descomp.Loess.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Filtro.lineal.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Loess.Optimo.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SelectModel.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SuavizamientoEstacional.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-estimar.recursiva.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexpo.ErrorARMA.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencial.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-BP.LB.test-pruebaDW1.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")

data_dir <- "./data/anex-EMC-SeriesIndiceEmpalmados-vtas nominales lineas mercancia-may2023.csv" 
data <- read.table(file.choose(), header = TRUE, sep = ";",
    skip = 10, dec = ",", colClasses = c(rep("NULL", 14), "numeric", rep("NULL", 7)))

serie <- ts(data, freq = 12, start = c(2013, 1))
serie_log <- ts(log(data), freq = 12, start = c(2013, 1))

# Punto 1
# Graficamos la serie en escala normal y lognormal
par(mfrow = c(1, 2))
plot(serie, ylab = "PIB Nominal", main = "Escala Normal")
plot(serie_log, ylab = "PIB Nominal", main = "Escala Lognormal")


# Punto 2 A
# Graficamos diferentes componentes de la serie para su análisis
descom_serie_log <- decompose(serie_log, type = "additive")
# 1. Graficando la tendencia
plot(descom.serie_log$trend)
# 2. Graficando la varianza
x_no_trend <- diff(serie_log)
plot(x_no_trend, ylab = expression(log(Y[t]) - log(Y[t - 1])))
# 3. Graficando boxplots para estacionalidad
boxplot(serie_log ~ cycle(serie_log), main = "Gráfico de boxplots")
abline(h = mean(x_no_trend))
# 4. Graficando periodograma para estacionalidad
periodogram(x_no_trend, main = "Periodograma")
# 5. Graficamos la ACF
acf(as.numeric(serie), lag.max = 36, ci.type = "ma",
    col = 4, ci.col = 2, main = "ACF de la serie")
acf(as.numeric(serie_log), lag.max = 36, ci.type = "ma",
    col = 4, ci.col = 2, main = "ACF de la serie logaritmica")

# Punto 2 B
# Ecuacion Teorica \hat{Y}=e^{\beta_0+\beta_1t+\beta_2t^2+\beta_3t^3+\sum_{j=1}^{k} [ \alpha_j sin(2\pi F_jt) + \gamma_j cos(2 \pi F_j t)] }+E_t \\ E_t \sim idd N(0,\sigma^2)
# Validación cruzada
#   Inicializo variables para el modelo de ajuste
m <- 12
n <- length(serie) - m

t <- 1:n
t2 <- t^2
t3 <- t^3

yt <- ts(serie[t], freq = m, start = c(2013, 1)) # Escala normal
yt_log <- ts(serie_log[t], freq = m, start = c(2013, 1)) # Escala lognormal

# Funciones trigonometricas en las frecuencias j/12, j=2,3,4,5
sen2 <- sin(pi * t / 3)
cos2 <- cos(pi * t / 3)
sen3 <- sin(pi * t / 2)
cos3 <- cos(pi * t / 2)
sen4 <- sin(2 * pi * t / 3)
cos4 <- cos(2 * pi * t / 3)
sen5 <- sin(5 * pi * t / 6)
cos5 <- cos(5 * pi * t / 6)

# Matriz de diseño
x <- data.frame(t, t2, t3, sen2, cos2, sen3, cos3, sen4, cos4, sen5, cos5)
head(x)

# Ajuste del modelo
modelo <- lm(yt ~ ., data = x) # DISCUSS: "x" se hace con los valores 1:113, y escala normal
summary(modelo)
#Calculo valores ajustados del modelo exponencial cubico con estacionalidad trigonometrica
varianza <- summary(modelo)$sigma^2
yt_hat_modelo_ajustado <- ts(
    exp(fitted(modelo)) * exp(varianza / 2),
    freq = m,
    start = start(yt) # DISCUSS: yt es escala normal
)

# Inicializo variables para el pronostico
t_pronostico <- c((n + 1):(n + m))
t2_pronostico <- t_pronostico^2
t3_pronostico <- t_pronostico^3

yt_pronostico <- ts(serie_log[t_pronostico], freq = m, start = c(2013, 1))

#Funciones trigonometricas en los periodos de pronostico
sen1_pronostico <- sin(pi * t_pronostico / 6)
cos1_pronostico <- cos(pi * t_pronostico / 6)
sen2_pronostico <- sin(pi * t_pronostico / 3)
cos2_pronostico <- cos(pi * t_pronostico / 3)
sen3_pronostico <- sin(pi * t_pronostico / 2)
cos3_pronostico <- cos(pi * t_pronostico / 2)
sen4_pronostico <- sin(2 * pi * t_pronostico / 3)
cos4_pronostico <- cos(2 * pi * t_pronostico / 3)
sen5_pronostico <- sin(5 * pi * t_pronostico / 6)
cos5_pronostico <- cos(5 * pi * t_pronostico / 6)
cos6_pronostico <- cos(pi * t_pronostico)

#Matriz de diseño en los pronosticos
x_pronosticos <- data.frame(
    t = t_pronostico,
    t2 = t2_pronostico,
    t3 = t3_pronostico,
    sen2 = sen2_pronostico,
    cos2 = cos2_pronostico,
    sen3 = sen3_pronostico,
    cos3 = cos3_pronostico,
    sen4 = sen4_pronostico,
    cos4 = cos4_pronostico,
    sen5 = sen5_pronostico,
    cos5 = cos5_pronostico)

head(x_pronosticos)

# Creo una lista con el nombre de los parametros
parameter_names <- c(paste0("beta", 0:3), "alfa2", "gamma2",
    "alfa3", "gamma3", "alfa4",
    "gamma4", "alfa5", "gamma5")
modelo_exponencial <- regexponencial(respuesta = yt, # DISCUSS: yt es escala normal
    data = x, # DISCUSS: x son los datos 1:113
    names.param = parameter_names)

summary(modelo_exponencial)

# Calculo valores ajustados del modelo 2
yt_hat_modelo_exponencial <- ts(fitted(modelo_exponencial),freq=m,start=start(yt))
# Grafico de los valores ajustados
plot(serie,ylab="PIB Nominal de productos de aseo",main="Modelo Exponencial Multiplicativo con estacionalidad trigonometrica")
lines(yt_hat_modelo_exponencial, col = 2, lwd = 2)
legend("topleft", legend = c("Original", "Modelo Exponencial2"), lty = 1, col = c(1, 2))

#Calculo de los criterios AIC y BIC en el modelo
nparmod2=length(coef(modelo_exponencial)[coef(modelo_exponencial)!=0])
Criterios2=exp.crit.inf.resid(residuales=residuals(modelo_exponencial),n.par=nparmod2)

#Pronosticos del modelo exponencial en la escala original, solo son de tipo puntual por ser modelo no lineal
pronostico_exponencial <- predict(modelo_exponencial, newdata = x_pronosticos, interval = "prediction", level = 0.95)
yt_pronostico_exponencial <- ts(pronostico_exponencial, freq=m, start = start(yt_pronostico))
#yt_pronostico_exponencial
#precision pronosticos puntuales modelo exponencial
accuracy(yt_pronostico_exponencial,yt_pronostico)

modelo1 <- Descomp.Loess(serie.ajuste=yt,h=m,tipo.descomp="multiplicative",grado=1,criterio="aicc")

#plot(modelo1$St,ylab=expression(hat(S)[t]),lwd=2)

#Serie de pronósticos puntuales 
#modelo_exponencial$tablapron