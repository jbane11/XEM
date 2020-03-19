	real function slac_c12_fit(x)

	real alpha,A,C

	A = 12

	alpha = -0.070+2.189*x - 24.667*x**2 + 145.291*x**3
     >        -497.237*x**4 + 1013.129*x**5 - 1208.393*x**6
     >        +775.767*x**7 - 205.872*x**8

	C = exp( 0.017 + 0.018*log(x) + 0.005*log(x)**2)

	
	slac_c12_fit = C*A**alpha

	return 

	end
