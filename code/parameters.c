int set_params() {

	//cosmology-specific parameters
	rho_cosm = 2.78E+11; //[h^2 M_sol Mpc^-3]
	
	G = 6.673e-11; //[N m^2 kg^-2]
	
	//HOD parameters at Mmax -19.5
	//http://arxiv.org/pdf/1005.2413v2.pdf
	log_Mmin = 11.57;
	sigma_logm = 0.17;
	log_M0 = 12.23;
	log_M1 = 12.87;
	M0 = pow(10.0, log_M0); //[M_sol]
	M1 = pow(10.0, log_M1); //[M_sol]
	alpha = 0.99;
	
	//critical linear overdensity for top-hat collapse
	delta_c = 1.686;
	
	//mass-concentration-redshift relation for CDM haloes
	c0 = 11.0;
	beta = -0.13;
	
	//virialisation overdensity treshold
	//only EdS??
	//del_nl = 0.004;
	del_nl = 200.0;
	//del_nl = 0.2;
	
	//from  seminal paper by Bullock et Bullock (not really, its from Bullock at al where they give the cdm = 1/(1+z) (m/mstar)^a
	//double omega_m_z = omega_m*pow(a_cur, -3.0)*pow(H0/H, 2.0);
	//https://arxiv.org/pdf/astro-ph/9908159v3.pdf
	//del_nl = del_nl(z) in functions.c 		

	//for D+ calculation
	lambda_F4 = 1.0/0.042; //[h Mpc^-1]	
	lambda_F5 = 1.0/0.133;
	lambda_F6 = 1.0/0.419;
	lambda_GR = 0.0;	
	
	
	
	return 0;	
}
