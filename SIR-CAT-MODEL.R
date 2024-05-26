# Coded by Santiago Caballero Rosas

############################################
#       |\      _,,,---,,_                 #
# ZZZzz /,`.-'`'    -.  ;-;;,_             #
#      |,4-  ) )-,_. ,\ (  `'-'            #
#     '---''(_/--'  `-'\_)  Modelo SIR-CAT # 
############################################

# Instalar la paquetería que necesitamos
library(deSolve)
library(rootSolve)
library(phaseR)

# Declarar los valores de parámetros que permanecen constantes
beta_HG = 0.2  # Tasa de transmisión de gatos a humanos
beta_GG = 0.3  # Tasa de transmisión entre gatos
gamma_H = 0.1  # Tasa de recuperación en humanos
gamma_G = 0.1  # Tasa de recuperación en gatos

# Definir las ecuaciones diferenciales del modelo de Toxoplasma
toxoplasma_model <- function(t, y, parms) {
  SH <- y[1]  # Humanos sanos
  IH <- y[2]  # Humanos infectados
  TH <- y[3]  # Humanos tratados
  SG <- y[4]  # Gatos sanos
  IG <- y[5]  # Gatos infectados
  
  dSH <- -beta_HG * IG * SH
  dIH <- beta_HG * IG * SH - gamma_H * IH
  dTH <- gamma_H * IH
  dSG <- -beta_GG * IG * SG
  dIG <- beta_GG * IG * SG
  
  return(list(c(dSH, dIH, dTH, dSG, dIG)))
}

# Condiciones iniciales
initial_conditions <- c(SH = 9800, IH = 200, TH = 0, SG = 300, IG = 200)

# Tiempo de integración
tspan <- seq(from = 0, to = 100, by = 0.01)

# Parámetros
parms <- c(beta_HG = beta_HG, beta_GG = beta_GG, gamma_H = gamma_H, gamma_G = gamma_G)

# Integrar el modelo
out <- ode(y = initial_conditions, times = tspan, func = toxoplasma_model, parms = parms)

# Ajustar los márgenes de la figura
par(mar = c(5, 5, 2, 2))  # Márgenes: inferior, izquierdo, superior, derecho

# Graficar los resultados de la dinámica de humanos
matplot(out[,1], out[,2:4], type = "l", xlab = "Tiempo", ylab = "Población", 
        col = c("yellow", "red", "green"), lty = 1, lwd = 2, 
        main = "Dinámica de la población humana")
legend("topright", legend = c("SH (Humanos sanos)", "IH (Humanos infectados)", "TH (Humanos tratados)"), 
       col = c("yellow", "red", "green"), lty = 1, lwd = 2)

# Graficar los resultados de la dinámica de gatos
matplot(out[,1], out[,5:6], type = "l", xlab = "Tiempo", ylab = "Población", 
        col = c("blue", "purple"), lty = 1, lwd = 2, 
        main = "Dinámica de la población de gatos")
legend("topright", legend = c("SG (Gatos sanos)", "IG (Gatos infectados)"), 
       col = c("blue", "purple"), lty = 1, lwd = 2)

# Graficar los resultados combinados
matplot(out[,1], out[,2:6], type = "l", xlab = "Tiempo", ylab = "Población", 
        col = c("yellow", "red", "green", "blue", "purple"), lty = 1, lwd = 2, 
        main = "Modelo SIR-CAT combinado")
legend("topright", legend = c("SH (Humanos sanos)", "IH (Humanos infectados)", "TH (Humanos tratados)", "SG (Gatos sanos)", "IG (Gatos infectados)"), 
       col = c("yellow", "red", "green", "blue", "purple"), lty = 1, lwd = 2)

################################################################
#PUNTOS DE EQUILIBRIO
################################################################

# Función para encontrar los puntos de equilibrio
find_equilibrium <- function(y) {
  SH <- y[1]
  IH <- y[2]
  TH <- y[3]
  SG <- y[4]
  IG <- y[5]
  
  return(c(
    -beta_HG * IG * SH,
    beta_HG * IG * SH - gamma_H * IH,
    gamma_H * IH,
    -beta_GG * IG * SG,
    beta_GG * IG * SG
  ))
}

# Parámetros iniciales para encontrar el equilibrio
initial_guess <- c(SH = 9800, IH = 200, TH = 0, SG = 300, IG = 200)

# Encontrar los puntos de equilibrio
equilibrium <- multiroot(f = find_equilibrium, start = initial_guess)$root
print(paste("Puntos de equilibrio: SH =", equilibrium[1], "IH =", equilibrium[2], 
            "TH =", equilibrium[3], "SG =", equilibrium[4], "IG =", equilibrium[5]))

# Calcular la matriz jacobiana en el punto de equilibrio
jacobian <- function(y) {
  SH <- y[1]
  IH <- y[2]
  TH <- y[3]
  SG <- y[4]
  IG <- y[5]
  
  J <- matrix(0, nrow = 5, ncol = 5)
  J[1, 1] <- -beta_HG * IG
  J[1, 5] <- -beta_HG * SH
  J[2, 1] <- beta_HG * IG
  J[2, 2] <- -gamma_H
  J[2, 5] <- beta_HG * SH
  J[3, 2] <- gamma_H
  J[4, 4] <- -beta_GG * IG
  J[4, 5] <- -beta_GG * SG
  J[5, 4] <- beta_GG * IG
  J[5, 5] <- beta_GG * SG
  
  return(J)
}

# Evaluar la estabilidad de los puntos de equilibrio
J_eq <- jacobian(equilibrium)
eigenvalues <- eigen(J_eq)$values
print("Autovalores de la matriz jacobiana en el punto de equilibrio:")
print(eigenvalues)

# Evaluar estabilidad
if (all(Re(eigenvalues) < 0)) {
  print("El punto de equilibrio es un atractor (estable).")
} else if (any(Re(eigenvalues) > 0)) {
  print("El punto de equilibrio es un repulsor (inestable).")
} else {
  print("El punto de equilibrio es un punto silla (ni atractor ni repulsor).")
}
