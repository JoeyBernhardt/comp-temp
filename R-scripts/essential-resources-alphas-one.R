
library(tidyverse)
all_rstars <- read_csv("data/all-rstars.csv")

#### goal here is to find the alphas depending on what zone we are in

### We use Arrhenius function to model the temperature dependence
arrhenius_function <- function(Temp, E, b1, ref_temp = 0) {
	k<-8.62e-05 #Boltzmann's constant
	E <- E # 0.6 # activation energy (eV)
	T<-Temp+273.15 #range of temp in K
	Tc <- ref_temp+273.15 #reference temperature
	
	metabolism<-(b1*exp(1)^(E*(1/(k*Tc)-1/(k*T))))
	return(metabolism)
}

T <- 30

# temp_dependences_MacArthur <- function(T = 25, r_Ea = 0.6, c_Ea = 0.6, K_Ea = 0, v_Ea = 0.0, m_Ea = 0.2, ref_temp2 = 0)

find_alphas <- function(T, r_Ea = 0.6, c_Ea_1P = 0.9, c_Ea_1N = 0.7, c_Ea_2P = 0.1, c_Ea_2N = 0.3, k_Ea = 0, m_Ea = 0.2, ref_temp = 0){

	c1P_b <- 0.1
	c2P_b <- 0.5 ## slope for species 2 is 2.5, 
	c1N_b <- 0.3
	c2N_b <- 0.2
	
	R1P_b <- 0.1
	R2P_b <- 0.5
	R1N_b <- 5
	R2N_b <- 1.6
	
	
	
	
	k1N_b <- 33
	k2N_b <- 17
	k1P_b <- 1.2
	k2P_b <- 2
	
	D <- 0.1
	
	## how can we define the consumption vector lines?
	
	SN <- 1000/200
	SP <- 140/350

	
	r1_b <- 1.2
	r2_b <- 1.25

	m1_b <- 0.1
	m2_b <- 0.1
	
	c1P = arrhenius_function(Temp = T, E = c_Ea_1P, b1 = c1P_b)
	c1N = arrhenius_function(Temp = T, E = c_Ea_1N, b1 = c1N_b)
	c2P = arrhenius_function(Temp = T, E = c_Ea_2P, b1 = c2P_b)
	c2N = arrhenius_function(Temp = T, E = c_Ea_2N, b1 = c2N_b) ### cij = per capita consumption of comsumer i on resource j
	
	k1N = arrhenius_function(Temp = T, E = k_Ea, b1 = k1N_b)
	k2P = arrhenius_function(Temp = T, E = k_Ea, b1= k2P_b) #half saturation constant for N resource consumption
	k1P = arrhenius_function(Temp = T, E = k_Ea, b1 = k1P_b)
	k2N = arrhenius_function(Temp = T, E = k_Ea, b1= k2N_b)
	
	r1 = arrhenius_function(Temp = T, E = r_Ea, b1 = r1_b)
	r2 = arrhenius_function(Temp = T, E = r_Ea, b1= r2_b) #population growth rates
	m1 = arrhenius_function(Temp = T, E = m_Ea, b1 = m1_b)
	m2 = arrhenius_function(Temp = T, E = m_Ea, b1 = m2_b) # mortality rates
	
	R1N <- (m1 * k1N)/ (r1 - m1)
	R1P <- (m1 * k1P)/ (r1 - m1)
	R2P <- (m2 * k2P)/ (r2 - m2)
	R2N <- (m2 * k2N)/ (r2 - m2)
	
	cons_vec1_intercept <- max(R1P, R2P) -(c1P/c1N)*max(R1N, R2N)
	cons_vec2_intercept <- max(R1P, R2P) -(c2P/c2N)*max(R1N, R2N)
	supply_vec <- SP/SN
	
	cons_vec1_fun <- function(x){
		y <- (c1P/c1N)*x + cons_vec1_intercept
		return(y)
	}
	
	cons_vec2_fun <- function(x){
		y <- (c2P/c2N)*x + cons_vec2_intercept
		return(y)
	}
	
	supply_vec_fun <- function(x){
		y <- (SP/SN)*x
		return(y)
	}
	
	zone_middle <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN) | cons_vec1_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec2_fun(SN)
	zone_bottom <- cons_vec2_fun(SN) >= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN)
	zone_top <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) >= cons_vec1_fun(SN)
	
	zones <- c("zone_middle", "zone_bottom", "zone_top")
	zone <- zones[c(zone_middle, zone_bottom, zone_top)]
	zone <- "zone_middle"
	alphas <- function(zone) {
		if (zone == "zone_middle") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_top") {
			a11 <- c1N / (D * (SN - R1N))
			a12 <- c2N / (D * (SN - R1N))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_bottom") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1P / (D * (SP - R2P))
			a22 <- c2P / (D * (SP - R2P))
		}
		alphas1 <- data.frame(a11 = a11, a12 = a12, a21 = a21, a22 = a22)
		return(alphas1)
	}
	
	alphas_calc <- alphas(zone)
	r_s <- data.frame(r1 = r1, r2 = r2)
	comps <- data.frame(D = D, zone = zone)
	
	
	
	rho <- sqrt((alphas_calc$a12*alphas_calc$a21)/(alphas_calc$a11*alphas_calc$a22)) #niche overlap
	stabil_potential <- 1 - rho #stabilizing potential
	fit_ratio <- sqrt((alphas_calc$a11*alphas_calc$a12)/(alphas_calc$a22*alphas_calc$a21)) #fitness ratio 
	coexist <- rho < fit_ratio &  fit_ratio < 1/rho
	
	
	rs <- data.frame(R1N = R1N, R1P = R1P, R2N = R2N, R2P= R2P,
			   K1 = (r1)/alphas_calc$a11, K2 = (r2)/alphas_calc$a22, T = T, m1 = m1, m2 = m2, c1P = c1P,  c1N = c1N,
			   c2P = c2P, c2N = c2N, stabil_potential = stabil_potential,
			   fit_ratio = fit_ratio, rho = rho, coexist = coexist,
					 a11 = alphas_calc$a11, a12 = alphas_calc$a12, a21 = alphas_calc$a21, a22 = alphas_calc$a22)
	
	alphas_calc2 <- bind_cols(r_s, comps, rs)
	return(alphas_calc2)
}

temps <- seq(1,50, by = 0.5)
results_1 <- temps %>% 
	map_df(find_alphas) 

results_1 %>% 
	ggplot(aes(x = T, y = fit_ratio)) + geom_point()

results_1 %>% 
	ggplot(aes(x = T, y = round(stabil_potential, digits = 2))) + geom_point() +
	ylab("Stabilization potential")

results_1 %>% 
	ggplot(aes(y = fit_ratio, x = round(stabil_potential, digits = 2))) + geom_point() +
	ylab("Fitness difference")

sqrt((alphas_calc$a12*alphas_calc$a21)/(alphas_calc$a11*alphas_calc$a22))

results_1 %>% 
	mutate(rho2 = sqrt((a12*a21)/(a11*a22))) %>%
	mutate(denom = a12*a21) %>% 
	mutate(numer = a11*a22) %>%
	mutate(thinger = denom/numer) %>% View

### Explore different temperature dependences
results_tilman <- data.frame()
for(i in 1:6){
	for(j in 1:6){
		r_Ea <- 0.3
		c_Ea <- 0.1*(i-1)
		m_Ea <- 0.1*(j-1)
		# hold <- temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), r_Ea = r_Ea , c_Ea = c_Ea, m_Ea = m_Ea, ref_temp2 = 1)
		hold <- find_alphas(T = seq(30, 30, by = 0.1), r_Ea = r_Ea , c_Ea = c_Ea, m_Ea = m_Ea, ref_temp = 1)
		hold$r_Ea <- r_Ea
		hold$c_Ea <- c_Ea
		hold$m_Ea <- m_Ea
		results_tilman <- bind_rows(results_tilman, hold)
	}
}

temps <- seq(1,30, by = 0.1)

hold <- find_alphas(T = T, r_Ea = r_Ea , c_Ea = c_Ea, m_Ea = m_Ea, ref_temp = 1)




results <- data.frame()
for(i in 1:350){
	for(j in 1:6){
		for(z in 1:6){
		r_Ea <- 0.3
	 	c_Ea <- 0.1*(z-1)
		# c_Ea <- 0.5
		m_Ea <- 0.1*(j-1)
		T <- 0.1*(i-1)
		hold <- find_alphas(T = T, r_Ea = r_Ea , c_Ea = c_Ea, m_Ea = m_Ea, ref_temp = 0)
		# hold <- find_alphas(T = seq(0, 35, by = 0.1), r_Ea = r_Ea , c_Ea = c_Ea, m_Ea = m_Ea, ref_temp = 1)
		hold$r_Ea <- r_Ea
		hold$c_Ea <- c_Ea
		hold$m_Ea <- m_Ea
		results <- bind_rows(results, hold)
	}}}


unique(results$rho)
length(unique(results$a11))
max(results$fit_ratio)
min(results$fit_ratio)

results %>% 
	# filter(fit_ratio < 2) %>% 
	# filter(T %in% c(0, 10, 20, 30, 40)) %>% 
	ggplot(aes(x = R1N, y = R1P, color = T)) + geom_point() +
	geom_segment(aes(x = R1N, y = R1P, xend = R1N + c1N/20, yend = R1P + c1P/20, color = T),
				 size = 0.5, linetype = 1, alpha = 1) +
	geom_segment(aes(x = R1N, y = R1P, xend = R1N, yend = 0.2, color = T),
				 size = 0.5, linetype = 1, alpha = 1) +
	geom_segment(aes(x = R1N, y = R1P, xend = 5, yend = R1P, color = T),
				 size = 0.5, linetype = 1, alpha = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	scale_colour_viridis_c() +
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	coord_cartesian() +
	# xlim(0, 0.15) + ylim(0, 0.15) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 
ggsave("figures/essential-resources-tilman-plots.png", width = 14, height = 8)


x <- seq(0.5,1.5,length = 100)
lines.df <- data.frame(x = x, y = 1/(x), y2 = x)

results %>% 
	# filter(fit_ratio < 2) %>% 
ggplot(aes(x = stabil_potential, y = fit_ratio, color = T, shape = coexist))+
	geom_hline(yintercept = 1, color = "grey", linetype = 2)+
	geom_vline(xintercept = 0, color = "grey", linetype = 2)+
	geom_line(data = lines.df, aes(x = 1-x, y = y, color = NULL, shape = NULL))+
	geom_line(data = lines.df, aes(x = 1-x, y = y2, color = NULL, shape = NULL))+
	geom_point()+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	scale_color_viridis_c(option = "magma", end = 0.8)+
	scale_shape_manual(values = c(0,19)) +
	ylab("Fitness ratio") + xlab("Niche differences") 
ggsave("figures/essential-resources-chesson-plots-one.png", width = 14, height = 8)


results %>% 
	ggplot(aes(x = T, y = fit_ratio)) + geom_point()

unique(results$c_Ea)
results %>% 
	dplyr::filter(m_Ea == "0") %>% 
	dplyr::filter(c_Ea == "0.1") %>% 
	ggplot(aes(x = T, y = stabil_potential)) + geom_line() +
	geom_hline(yintercept = 0) +
	ylab("Niche difference")

unique(results$c_Ea)
results %>% 
	# dplyr::filter(m_Ea == "0") %>% 
	# dplyr::filter(c_Ea == "0") %>% 
	ggplot(aes(x = T, y = fit_ratio)) + geom_line(color = "blue") +
	geom_hline(yintercept = 0) +
	ylab("Fitness ratio") +
	facet_grid(c_Ea ~ m_Ea) + xlab("Temperature (°C)")
ggsave("figures/fitness-ratio-temperature-two.png", width = 14, height = 12)

results %>% 
	# dplyr::filter(m_Ea == "0") %>% 
	# dplyr::filter(c_Ea == "0") %>% 
	ggplot(aes(x = T, y = stabil_potential)) + geom_line(color = "blue") +
	geom_hline(yintercept = 0) +
	ylab("Niche difference") +
	facet_grid(c_Ea ~ m_Ea) + xlab("Temperature (°C)")
ggsave("figures/niche-difference-temperature-two.png", width = 14, height = 12)

results %>% 
	# filter(r1 > 0) %>%
	ggplot(aes(x = T, y = R2N, group = m_Ea, color = m_Ea))+
	geom_line()+
	scale_color_viridis_c(end = 0.8)+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	ylab(label = "N* (species 2)") + xlab("Temperature (°C)")

results %>% 
	# filter(r1 > 0) %>%
	ggplot(aes(x = T, y = R2P, group = m_Ea, color = m_Ea))+
	geom_line()+
	scale_color_viridis_c(end = 0.8)+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	ylab(label = "P* (species 1)") + xlab("Temperature (°C)")
ggsave("figures/essential-resources-P-star-plots.png", width = 14, height = 8)
