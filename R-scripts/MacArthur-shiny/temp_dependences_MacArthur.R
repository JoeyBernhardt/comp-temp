temp_dependences_MacArthur <- function(T = 25, ref_temp2 = 1,
									   r_EaN = 0.5, r_EaP = 0.5, 
									   c_Ea1N = 0.8, c_Ea1P = 0.8, 
									   c_Ea2N = 0.8, c_Ea2P = 0.8, 
									   K_EaN = -0.3, K_EaP = -0.3, 
									   v_EaN = 0.0, v_EaP = 0.0, 
									   m_Ea1 = 0.6, m_Ea2 = 0.6,
									   c1N_b = 0.2, c2P_b = 0.2,
									   c1P_b = 0.4, c2N_b = 0.4,
									   v1N_b = 0.2, v2N_b = 0.4, 
									   v1P_b = 0.4, v2P_b = 0.2
){
	
	# resource growth rates
	rN = arrhenius_function(Temp = T, E = r_EaN, b1 = 0.05) 
	rP = arrhenius_function(Temp = T, E = r_EaP, b1 = 0.05) 
	
	# resource carrying capacity
	KN = arrhenius_function(Temp = T, E = K_EaN, b1 = 2) 
	KP = arrhenius_function(Temp = T, E = K_EaP, b1 = 2) 
	
	# cij = per capita consumption of comsumer i on resource j
	c1N = arrhenius_function(Temp = T, E = c_Ea1N, b1 = c1N_b)
	c1P = arrhenius_function(Temp = T, E = c_Ea1P, b1 = c1P_b) ## species 1 consumes more P than N
	c2N = arrhenius_function(Temp = T, E = c_Ea2N, b1 = c2N_b) ## species 2 consumes more N than P
	c2P = arrhenius_function(Temp = T, E = c_Ea2P, b1 = c2P_b) 
	
	# vij = conversion factor that converts resource j into biomass of consumer i
	v1N = arrhenius_function(Temp = T, E = v_EaN, b1 = v1N_b)
	v2N = arrhenius_function(Temp = T, E = v_EaN, b1 = v2N_b) ## species 2 converts N more efficiently 
	v1P = arrhenius_function(Temp = T, E = v_EaP, b1 = v1P_b) ## species 1 converts P more efficiently 
	v2P = arrhenius_function(Temp = T, E = v_EaP, b1 = v2P_b)
	
	# mortality rates
	m1 = arrhenius_function(Temp = T, E = m_Ea1, b1 = 0.01)
	m2 = arrhenius_function(Temp = T, E = m_Ea2, b1 = 0.01)
	
	# Absolute competition coefficients
	beta11 = v1N * c1N * (KN/rN) * c1N + v1P * c1P * (KP/rP) * c1P ### intra
	beta12 = v1N * c1N * (KN/rN) * c2N + v1P * c1P * (KP/rP) * c2P ### inter
	beta22 = v2N * c2N * (KN/rN) * c2N + v2P * c2P * (KP/rP) * c2P ### intra
	beta21 = v2N * c2N * (KN/rN) * c1N + v2P * c2P * (KP/rP) * c1P ### inter
	
	
	g1 = v1N * c1N * KN + v1P * c1P * KP - m1 ### growth rate of the consumer 1
	g2 = v2N * c2N * KN + v2P * c2P * KP - m2 ### growth rate of consumer 2
	
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
			   m1 = m1, m2 = m2, rN = rN, rP = rP, KN = KN, KP = KP,
			   c1N = c1N,  c1P = c1P,  c2N = c2N, c2P = c2P)
}

