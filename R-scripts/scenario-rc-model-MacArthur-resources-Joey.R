### PJK Apr 30 2020
library(tidyverse)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot())



### Arrhenius function to model the temperature dependence
arrhenius_function <- function(Temp, E, b1, ref_temp = 1) {
  k <- 8.62e-05 #Boltzmann's constant
  E <- E # 0.6 # activation energy (eV)
  T <- Temp+273.15 #range of temp in K
  Tc <- ref_temp+273.15 #reference temperature
  
  metabolism<-(b1*exp(1)^(E*(1/(k*Tc)-1/(k*T))))
  return(metabolism)
}



### Function to calculate MCT components under different temperatures
temp_dependences_MacArthur <- function(T = 25, ref_temp2 = 1,
                                       r_Ea1 = 0.5, r_Ea2 = 0.5, 
                                       c_Ea11 = 0.8, c_Ea12 = 0.8, 
                                       c_Ea21 = 0.8, c_Ea22 = 0.8, 
                                       K_Ea1 = -0.3, K_Ea2 = -0.3, 
                                       v_Ea1 = 0.0, v_Ea2 = 0.0, 
                                       m_Ea1 = 0.6, m_Ea2 = 0.6){
  
  # resource growth rates
  r1 = arrhenius_function(Temp = T, E = r_Ea1, b1 = 0.05) 
  r2 = arrhenius_function(Temp = T, E = r_Ea2, b1 = 0.05) 
  
  # resource carrying capacity
  K1 = arrhenius_function(Temp = T, E = K_Ea1, b1 = 2) 
  K2 = arrhenius_function(Temp = T, E = K_Ea2, b1 = 2) 
  
  # cij = per capita consumption of comsumer i on resource j
  c11 = arrhenius_function(Temp = T, E = c_Ea11, b1 = 0.2)
  c12 = arrhenius_function(Temp = T, E = c_Ea12, b1 = 0.4)
  c21 = arrhenius_function(Temp = T, E = c_Ea21, b1 = 0.4)
  c22 = arrhenius_function(Temp = T, E = c_Ea22, b1 = 0.2) 
  
  # vij = conversion factor that converts resource j into biomass of consumer i
  v11 = arrhenius_function(Temp = T, E = v_Ea1, b1 = 0.2) 
  v12 = arrhenius_function(Temp = T, E = v_Ea1, b1 = 0.4) 
  v21 = arrhenius_function(Temp = T, E = v_Ea2, b1 = 0.4) 
  v22 = arrhenius_function(Temp = T, E = v_Ea2, b1 = 0.2)
  
  # mortality rates
  m1 = arrhenius_function(Temp = T, E = m_Ea1, b1 = 0.01)
  m2 = arrhenius_function(Temp = T, E = m_Ea2, b1 = 0.01)
  
  # Absolute competition coefficients
  beta11 = v11 * c11 * (K1/r1) * c11 + v12 * c12 * (K2/r2) * c12
  beta21 = v21 * c21 * (K1/r1) * c11 + v22 * c22 * (K2/r2) * c12
  beta22 = v21 * c21 * (K1/r1) * c21 + v22 * c22 * (K2/r2) * c22
  beta12 = v11 * c11 * (K1/r1) * c21 + v12 * c12 * (K2/r2) * c22
  g1 = v11 * c11 * K1 + v12 * c12 * K2 - m1
  g2 = v21 * c21 * K1 + v22 * c22 * K2 - m2
  
  # Relative competition coefficients 
  a11 = beta11 / g1
  a21 = beta21 / g2
  a22 = beta22 / g2
  a12 = beta12 / g1
  
  # MCT components
  rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
  stabil_potential <- 1 - rho #stabilizing potential
  fit_ratio <- sqrt((a11*a12)/(a22*a21))  #fitness ratio
  coexist <- rho < fit_ratio &  fit_ratio < 1/rho
  
  # report results
  data.frame(T = T, 
             a11 = a11, a12 = a12, a22 = a22, a21 = a21, g1 = g1, g2 = g2,
             stabil_potential = stabil_potential, fit_ratio = fit_ratio, rho = rho, coexist = coexist, 
             m1 = m1, m2 = m2, r1 = r1, r2 = r2, K1 = K1, K2 = K2,
             c11 = c11,  c12 = c12,  c21 = c21, c22 = c22)
}



### Function to plot: parameter values, consumption vector slopes, ND & FD
plot_MacArthur <- function(Data.parameter, Data.temperature){
  
  # panel a: parameter value ~ temperature | species
  panel.a = 
    ggplot(Data.parameter, 
           aes(x = T, y = value, col = parameter)) +
    geom_point() +
    xlab("Temperature") +
    ylab(label = "Parameter value") + 
    ggtitle("(a)") +
    theme(legend.position = "bottom")
  
  # panel b: consumption vector slope ~ temperature | species
  Temp = data.frame(T = Data.temperature$T, 
                    slope1 = Data.temperature$c12 / Data.temperature$c11, 
                    slope2 = Data.temperature$c22 / Data.temperature$c21)
  Temp = Temp %>% gather(value=value, key=parameter, -T)
  panel.b = 
    ggplot(Temp, 
           aes(x = T, y = value, col = parameter)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype=5) +
    xlab("Temperature") +
    ylab("Parameter value") + 
    ggtitle("(b)") +
    theme(legend.position = "bottom")
  
  # panel c: ND ~ temperature
  panel.c = 
    Data.temperature %>%
    filter((g1 > 0) & (g2 > 0)) %>%
    ggplot(aes(x = T, y = stabil_potential, col = T)) +
    geom_point() +
    scale_color_viridis_c("Temperature", end = 0.8) +
    xlab("Temperature") +
    ylab(expression(paste("Stabilization potential (1-", rho, ")"))) + 
    ggtitle("(c)") +
    theme(legend.position = "bottom")
  
  # panel d: FD ~ temperature
  panel.d = 
    Data.temperature %>%
    filter((g1 > 0) & (g2 > 0)) %>%
    ggplot(aes(x = T, y = fit_ratio, col = T)) +
    geom_point() +
    scale_color_viridis_c("Temperature", end = 0.8) +
    xlab("Temperature") +
    ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
    ggtitle("(d)") +
    theme(legend.position = "bottom")
  
  # panel e: FD ~ ND
  panel.e = 
    Data.temperature %>%
    filter((g1 > 0) & (g2 > 0)) %>%
    ggplot(aes(x = stabil_potential, y = fit_ratio, col = T)) +
    geom_point() +
    geom_ribbon(data = data.frame(x = seq(min(Data.temperature$stabil_potential)*0.99, max(Data.temperature$stabil_potential)*1.01, 0.001)),
                aes(x = x,
                    y = NULL, 
                    ymin = 1-x,
                    ymax = 1/(1-x)),
                fill = "grey", color = "black", alpha = 0.3) +
    geom_hline(yintercept = 1, linetype=5) +
    scale_x_continuous(expand=c(0, 0)) + 
    scale_color_viridis_c("Temperature", end = 0.8) +
    xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
    ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
    ggtitle("(e)") +
    theme(legend.position = "bottom")
  
  # Combine panels
  return(panel.a + panel.b + panel.c + panel.d + panel.e + plot_layout(nrow=1))
  
}



### Function to plot: alpha, competitive impact, sensitivity to competition
plot_MacArthur_alpha <- function(Data.temperature){
  
  # panel a: alpha ~ temperature | species
  Data.alpha = Data.temperature[, c("T", "a11", "a12", "a22", "a21")] %>% gather(value=value, key=parameter, -T)
  panel.a = 
    ggplot(Data.alpha, 
           aes(x = T, y = value, col = parameter)) +
    geom_point() +
    scale_color_manual(values=c("a11"="#08519c", "a21"="#6baed6", "a22"="#a50f15", "a12"="#fb6a4a")) + 
    xlab("Temperature") +
    ylab(expression(alpha[ij])) + 
    ggtitle("(a)") +
    theme(legend.position = "bottom")
  
  # panel b: relative impact ~ temperature | species
  Temp = data.frame(T = Data.temperature$T, 
                    rel.impact1 = Data.temperature$a21 / Data.temperature$a11, 
                    rel.impact2 = Data.temperature$a12 / Data.temperature$a22)
  Temp = Temp %>% gather(value=value, key=parameter, -T)
  panel.b = 
    ggplot(Temp, 
           aes(x = T, y = value, col = parameter)) +
    geom_point() +
    scale_color_manual(values=c("rel.impact1"="#08519c", "rel.impact2"="#a50f15")) + 
    xlab("Temperature") +
    ylab(expression(paste("Relative impact (", alpha[ji], "/", alpha[ii], ")"))) + 
    ggtitle("(b)") +
    theme(legend.position = "bottom")
  
  # panel c: sensitivity to competition ~ temperature | species
  Temp2 = data.frame(T = Data.temperature$T, 
                     sen1 = Data.temperature$a12 * Data.temperature$a11, 
                     sen2 = Data.temperature$a21 * Data.temperature$a22)
  Temp2 = Temp2 %>% gather(value=value, key=parameter, -T)
  panel.c = 
    ggplot(Temp2, 
           aes(x = T, y = value, col = parameter)) +
    geom_point() +
    scale_color_manual(values=c("sen1"="#08519c", "sen2"="#a50f15")) + 
    xlab("Temperature") +
    ylab(expression(paste("Sensitivity (", alpha[ii], "*", alpha[ij], ")"))) + 
    ggtitle("(c)") +
    theme(legend.position = "bottom")
  
  # Combine panels
  return(panel.a + panel.b + panel.c + plot_layout(nrow=1))
  
}



### SIMULATE !
### r varies with temperature & species have different activation energy for r (benefitting Sp1: r_Ea2 > r_Ea1)
# ### Other parameters vary with temperature (same AE between species for other parameters) 
# Data.temperature.r = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
#                                                 r_Ea1 = 0.3, r_Ea2 = 0.5, 
#                                                 K_Ea1 = -0.3, K_Ea2 = -0.3, 
#                                                 c_Ea11 = 0.8, c_Ea12 = 0.8, 
#                                                 c_Ea21 = 0.8, c_Ea22 = 0.8, 
#                                                 v_Ea1 = 0.0, v_Ea2 = 0.0, 
#                                                 m_Ea1 = 0.6, m_Ea2 = 0.6)

### Other parameters do not vary with temperature (AE=0) to disentangle underlying mechanism 
### Patterns are similar with previous simulations
Data.temperature.r = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
                                                r_Ea1 = 0.3, r_Ea2 = 0.5, 
                                                K_Ea1 = -0.0, K_Ea2 = -0.0, 
                                                c_Ea11 = 0.0, c_Ea12 = 0.0, 
                                                c_Ea21 = 0.0, c_Ea22 = 0.0, 
                                                v_Ea1 = 0.0, v_Ea2 = 0.0, 
                                                m_Ea1 = 0.0, m_Ea2 = 0.0)
Data.parameter.r = Data.temperature.r[, c("T", "r1", "r2")] %>% gather(value=value, key=parameter, -T)
Plot.r = plot_MacArthur(Data.parameter.r, Data.temperature.r)
Plot.alpha.r = plot_MacArthur_alpha(Data.temperature.r)
ggsave(file="../figures/MacArthur-scenario-r.png", plot=Plot.r, device = "png", width=17.5, height=4)
ggsave(file="../figures/MacArthur-scenario-r-alpha.png", plot=Plot.alpha.r, device = "png", width=10.5, height=4)



### K varies with temperature & species have different activation energy for K (benefitting Sp1: K_Ea2 > K_Ea1)
# ### Other parameters vary with temperature (same AE between species for other parameters) 
# Data.temperature.K = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
#                                                 r_Ea1 = 0.5, r_Ea2 = 0.5, 
#                                                 K_Ea1 = -0.3, K_Ea2 = -0.1, 
#                                                 c_Ea11 = 0.8, c_Ea12 = 0.8, 
#                                                 c_Ea21 = 0.8, c_Ea22 = 0.8, 
#                                                 v_Ea1 = 0.0, v_Ea2 = 0.0, 
#                                                 m_Ea1 = 0.6, m_Ea2 = 0.6)

### Other parameters do not vary with temperature (AE=0) to disentangle underlying mechanism
### Patterns are similar with previous simulations
Data.temperature.K = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
                                                r_Ea1 = 0.0, r_Ea2 = 0.0, 
                                                K_Ea1 = -0.3, K_Ea2 = -0.1, 
                                                c_Ea11 = 0.0, c_Ea12 = 0.0, 
                                                c_Ea21 = 0.0, c_Ea22 = 0.0, 
                                                v_Ea1 = 0.0, v_Ea2 = 0.0, 
                                                m_Ea1 = 0.0, m_Ea2 = 0.0)
Data.parameter.K = Data.temperature.K[, c("T", "K1", "K2")] %>% gather(value=value, key=parameter, -T)
Plot.K = plot_MacArthur(Data.parameter.K, Data.temperature.K)
Plot.alpha.K = plot_MacArthur_alpha(Data.temperature.K)
ggsave(file="../figures/MacArthur-scenario-K.png", plot=Plot.K, device = "png", width=17.5, height=4)
ggsave(file="../figures/MacArthur-scenario-K-alpha.png", plot=Plot.alpha.K, device = "png", width=10.5, height=4)



### m varies with temperature & species have different activation energy for m (benefitting Sp1: m_Ea1 < m_Ea2)
# ### Other parameters vary with temperature (same AE between species for other parameters) 
# Data.temperature.m = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
#                                                 r_Ea1 = 0.5, r_Ea2 = 0.5, 
#                                                 K_Ea1 = -0.3, K_Ea2 = -0.3, 
#                                                 c_Ea11 = 0.8, c_Ea12 = 0.8, 
#                                                 c_Ea21 = 0.8, c_Ea22 = 0.8, 
#                                                 v_Ea1 = 0.0, v_Ea2 = 0.0, 
#                                                 m_Ea1 = 0.4, m_Ea2 = 0.6)

### Other parameters do not vary with temperature (AE=0) to disentangle underlying mechanism
### Patterns are similar with previous simulations
Data.temperature.m = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
                                                r_Ea1 = 0.0, r_Ea2 = 0.0, 
                                                K_Ea1 = -0.0, K_Ea2 = -0.0, 
                                                c_Ea11 = 0.0, c_Ea12 = 0.0, 
                                                c_Ea21 = 0.0, c_Ea22 = 0.0, 
                                                v_Ea1 = 0.0, v_Ea2 = 0.0, 
                                                m_Ea1 = 0.4, m_Ea2 = 0.6)
Data.parameter.m = Data.temperature.m[, c("T", "m1", "m2")] %>% gather(value=value, key=parameter, -T)
Plot.m = plot_MacArthur(Data.parameter.m, Data.temperature.m)
Plot.alpha.m = plot_MacArthur_alpha(Data.temperature.m)
ggsave(file="../figures/MacArthur-scenario-m.png", plot=Plot.m, device = "png", width=17.5, height=4)
ggsave(file="../figures/MacArthur-scenario-m-alpha.png", plot=Plot.alpha.m, device = "png", width=10.5, height=4)



### c varies with temperature & species have different activation energy for c (benefitting Sp1: c_Ea1i > c_Ea2i)
# ### Other parameters vary with temperature (same AE between species for other parameters) 
# Data.temperature.c = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1),
#                                                 r_Ea1 = 0.5, r_Ea2 = 0.5,
#                                                 K_Ea1 = -0.3, K_Ea2 = -0.3,
#                                                 c_Ea11 = 0.8, c_Ea12 = 0.8,
#                                                 c_Ea21 = 0.6, c_Ea22 = 0.6,
#                                                 v_Ea1 = 0.0, v_Ea2 = 0.0,
#                                                 m_Ea1 = 0.6, m_Ea2 = 0.6)
# Data.parameter.c = Data.temperature.c[, c("T", "c11", "c12", "c21", "c22")] %>% gather(value=value, key=parameter, -T)
# Plot.c = plot_MacArthur(Data.parameter.c, Data.temperature.c)
# Plot.alpha.c = plot_MacArthur_alpha(Data.temperature.c)
# ggsave(file="../figures/MacArthur-scenario-c-2.png", plot=Plot.c, device = "png", width=17.5, height=4)
# ggsave(file="../figures/MacArthur-scenario-c-alpha-2.png", plot=Plot.alpha.c, device = "png", width=10.5, height=4)

### Other parameters do not vary with temperature (AE=0) to disentangle underlying mechanism
### Patterns are different with previous simulations, indicating interactive effect among parameters
Data.temperature.c = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
                                                r_Ea1 = 0.0, r_Ea2 = 0.0, 
                                                K_Ea1 = -0.0, K_Ea2 = -0.0, 
                                                c_Ea11 = 0.8, c_Ea12 = 0.8, 
                                                c_Ea21 = 0.6, c_Ea22 = 0.6, 
                                                v_Ea1 = 0.0, v_Ea2 = 0.0, 
                                                m_Ea1 = 0.0, m_Ea2 = 0.0)
Data.parameter.c = Data.temperature.c[, c("T", "c11", "c12", "c21", "c22")] %>% gather(value=value, key=parameter, -T)
Plot.c = plot_MacArthur(Data.parameter.c, Data.temperature.c)
Plot.alpha.c = plot_MacArthur_alpha(Data.temperature.c)
ggsave(file="../figures/MacArthur-scenario-c.png", plot=Plot.c, device = "png", width=17.5, height=4)
ggsave(file="../figures/MacArthur-scenario-c-alpha.png", plot=Plot.alpha.c, device = "png", width=10.5, height=4)



### preference varies with temperature & species have different activation energy for c (c_Ea11 > c_Ea12 & c_Ea22 > c_Ea21)
# ### Other parameters vary with temperature (same AE between species for other parameters) 
# Data.temperature.p = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1),
#                                                 r_Ea1 = 0.5, r_Ea2 = 0.5,
#                                                 K_Ea1 = -0.3, K_Ea2 = -0.3,
#                                                 c_Ea11 = 0.95, c_Ea12 = 0.65,
#                                                 c_Ea21 = 0.6, c_Ea22 = 0.8,
#                                                 v_Ea1 = 0.0, v_Ea2 = 0.0,
#                                                 m_Ea1 = 0.6, m_Ea2 = 0.6)

### Other parameters do not vary with temperature (AE=0) to disentangle underlying mechanism
### Patterns are similar with previous simulations
Data.temperature.p = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
                                                r_Ea1 = 0.0, r_Ea2 = 0.0, 
                                                K_Ea1 = -0.0, K_Ea2 = -0.0, 
                                                c_Ea11 = 0.95, c_Ea12 = 0.65, 
                                                c_Ea21 = 0.6, c_Ea22 = 0.8, 
                                                v_Ea1 = 0.0, v_Ea2 = 0.0, 
                                                m_Ea1 = 0.0, m_Ea2 = 0.0)
Data.parameter.p = Data.temperature.p[, c("T", "c11", "c12", "c21", "c22")] %>% gather(value=value, key=parameter, -T)
Plot.p = plot_MacArthur(Data.parameter.p, Data.temperature.p)
Plot.alpha.p = plot_MacArthur_alpha(Data.temperature.p)
ggsave(file="../figures/MacArthur-scenario-p.png", plot=Plot.p, device = "png", width=14, height=4)
ggsave(file="../figures/MacArthur-scenario-p-alpha.png", plot=Plot.alpha.p, device = "png", width=10.5, height=4)




