double z_to_t_int(double z, void * params) {
	Parameters *parameters = (Parameters*) params;	
	return 1.0/((1.0 + z)*H_a(z_to_a(z), parameters));
}

double select_lambda(int gravity) {
	double lambda_zero;
	switch (gravity) {
		case F4:
			lambda_zero = lambda_F4;
			break;
		case F5:
			lambda_zero = lambda_F5;
			break;
		case F6:
			lambda_zero = lambda_F6;
			break;
		case GR:
			lambda_zero = lambda_GR;
			break;
		default:
			printf("Gravity can only be one of the 4 (integers): F4, F5, F6, GR\n");
			exit(0);	
	}		
	return lambda_zero;
}

double z_to_t(double z, Parameters *params) {
	double error, result;
	gsl_function F;	
	F.function = &z_to_t_int;
	F.params = params;  	
	gsl_integration_qagiu(&F, z, 0, 1e-7, z_to_t_workspace_size, z_to_t_workspace, &result, &error);	
	return result;	
}

double t_to_z(double t, Spline *redshift_reverse) {
	return splint_generic(redshift_reverse, t);
}

double a_to_t(double a, Parameters *params) {
	return z_to_t(a_to_z(a), params);
}

Spline* prep_redshift_spline(BinInfo *z_binInfo, Parameters *params) {

	double z_cur, t_cur;
	Spline* toReturn;
	double* zs = malloc((size_t)z_binInfo->bins * sizeof(*zs));
	double* ts = malloc((size_t)z_binInfo->bins * sizeof(*ts));	
	
	for (int i = 0; i < z_binInfo->bins; i++) {
		z_cur = bin_to_x(z_binInfo, i);
		t_cur = z_to_t(z_cur, params);
		zs[i] = z_cur;
		ts[i] = t_cur;
	}		

	toReturn = new_spline(z_binInfo->bins, ts, zs);
	double minval = splint_generic(toReturn, toReturn->xmin);
	double maxval = splint_generic(toReturn, toReturn->xmax);
	
	if (!equalsDouble(z_binInfo->xmax, minval)) {
		printf("redshift spline didn't go well (Dplus.c prep_redshift_spline()) \n");
		printf("%lf \t %lf \n", minval, z_binInfo->xmax);
		exit(0);
	}
	
	if (!equalsDouble(z_binInfo->xmin, maxval)) {
		printf("redshift spline didn't go well (Dplus.c prep_redshift_spline()) \n");
		printf("%lf \t %lf \n", maxval, z_binInfo->xmin);
		exit(0);
	}
	
	return toReturn;
}

/*
double t_to_a_analytic(double t) {
	double t_lambda = 2.0/(3.0*H0*sqrt(omega_v));
	return pow(omega_m/omega_v, 1.0/3.0)*pow(sinh(t/t_lambda),2.0/3.0);
}
*/

//http://arxiv.org/pdf/1412.5195v3.pdf (19)
//f[0] - Dplus dot, f[1] - Dplus double dot
//y[0] - Dplus evaluated at some value of a, y[1] - Dplus dot, evaluated at some a
int growth_function(double t, const double Dplus_vals[], double dDplus_dt[], void *params) {
	double lambda_zero, lambda, k, omega_m_z_cur, bracket_factor, a_cur, z_cur, t0, H;
	double R_bar, R0_bar;
	
	Minim_Params *parameters = (Minim_Params*) params;
	lambda_zero = parameters->lambda0;
	k = parameters->k;
	t0 = parameters->t0;
	Spline* redshift_reverse = parameters->redshift_reverse;
	Parameters *cosmo_params = parameters->cosmo_params;
	
	//generic
	z_cur = t_to_z(t, redshift_reverse);	
	a_cur = z_to_a(z_cur);
		
	//omega_m + omega_v
	//a_cur = t_to_a_analytic(t);
	//z_cur = a_to_z(a_cur);
	
	//EdS
	//a_cur = pow(t/t0, 2.0/3.0);
	
	H = H_a(a_cur, cosmo_params);
	omega_m_z_cur = omega_m_z(z_cur, cosmo_params);	
	
	
	//(7), (17)
	R0_bar = 3.0*pow(cosmo_params->H0, 2.0)*(cosmo_params->omega_m_0 + 4.0*cosmo_params->omega_v_0); 
	R_bar = 3.0*pow(cosmo_params->H0, 2.0)*(cosmo_params->omega_m_0*pow(a_cur, -3.0) + 4.0*cosmo_params->omega_v_0);
	lambda = sqrt(pow(lambda_zero,2.0)*pow(R0_bar/R_bar,3.0));	
	bracket_factor = lambda*lambda*k*k/(a_cur*a_cur);
	
	dDplus_dt[0] = Dplus_vals[1]; //the value of the first derivative 
	dDplus_dt[1] = (3.0/2.0)*H*H*omega_m_z_cur*(1.0 + (1.0/3.0)*(bracket_factor)/(1.0 + bracket_factor))*Dplus_vals[0] - 2.0*H*Dplus_vals[1]; //the value of the second derivative
	return GSL_SUCCESS;
}

//only used when the max redshift specified is not in a matter-dominated period, need to find initial conditions
double *prep_Dplus_calc(int gravity, double k, BinInfo *z_binInfo, Parameters *params) {
	double z_EdS = 4.0;
	double t_start, t_end, Dplus_a_start, Dplusdot_a_start, dt, t_curr, t_next;
	double z_start, z_end, t0, lambda_zero, z_cur;
	int status;
	
	z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	Spline* redshift_reverse = prep_redshift_spline(z_binInfo, params);
	
	z_start = z_EdS;
	z_end = z_binInfo->xmax;
	t_start = z_to_t(z_start, params);	
	t_end = z_to_t(z_end, params);
	t0 = t_start;	
	dt = (t_end - t_start)/(double)(Dplus_bins - 1);
	
	//for EdS (matter dominated, flat) analytic
	Dplus_a_start = 1.0;
	Dplusdot_a_start = 2.0/(3.0*t0);
	
	lambda_zero = select_lambda(gravity);	
	
	double current_Dplus_values[2] = {Dplus_a_start, Dplusdot_a_start};

	Minim_Params *parameters = malloc(sizeof(*parameters));
	parameters->lambda0 = lambda_zero;
	parameters->k = k;
	parameters->t0 = t0;
	parameters->redshift_reverse = redshift_reverse;
	parameters->cosmo_params = params;
	
	gsl_odeiv2_system sys = {growth_function, NULL, 2, parameters}; //rk8pd
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, dt, 1e-7, 0.0);		

	t_curr = t_start;		
	for (int i = 1; i < Dplus_bins; i++) {	
		//the next iteration goes from t_curr to t_next in some number of steps chosen by the driver
		t_next = t_curr + dt;			
		status = gsl_odeiv2_driver_apply(d, &t_curr, t_next, current_Dplus_values);		
		
		if (status != GSL_SUCCESS) break;	
	}	
	gsl_odeiv2_driver_free(d);
	gsl_integration_workspace_free(z_to_t_workspace);
	free(redshift_reverse);
	free(parameters);	
	double *toReturn = malloc(2*sizeof(*toReturn));
	toReturn[0] = current_Dplus_values[0];
	toReturn[1] = current_Dplus_values[1];
	return toReturn;
}



//int Dplus_calc(int gravity, double k, double* zs, double* Dpluses)
int Dplus_calc(int gravity, double k, double* zs, double* Dpluses, BinInfo *z_binInfo, Parameters *params) {
	double t_start, t_end, Dplus_a_start, Dplusdot_a_start, dt, t_curr, t_next;
	double z_start, z_end, t0, lambda_zero, z_cur;
	double z_EdS = 4.0;
	int status;	
	
	if (z_binInfo->xmin > z_binInfo->xmax) {
		printf("check z binning, swap zmax, zmin order for Dplus calc \n");
		exit(0);
	}
	
	//redshift-time relation reversed:	
	z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	Spline* redshift_reverse = prep_redshift_spline(z_binInfo, params);
	
	//can be done better - since now z_start has to be quite high. If z_start < 3, Do one run from z~3 to z_start, figure 
	//out dplus, dplus dot at the redshift wanted, then use them again to get the dpluses correctly
	z_start = z_binInfo->xmax;
	z_end = z_binInfo->xmin;
	t_start = z_to_t(z_start, params);	
	t_end = z_to_t(z_end, params);
	t0 = t_start;	
	dt = (t_end - t_start)/(double)(Dplus_bins - 1);
	
	if (z_binInfo->xmax < z_EdS) {
		double *init_values = prep_Dplus_calc(gravity, k, z_binInfo, params);
		Dplus_a_start = init_values[0];
		Dplusdot_a_start = init_values[1];
	} else {
	//starting at a redshift where EdS applies, Dplus = (t/t0)^(2/3)
		Dplus_a_start = 1.0;
		Dplusdot_a_start = 2.0/(3.0*t0);
	}
	
	lambda_zero = select_lambda(gravity);	
	
	double current_Dplus_values[2] = {Dplus_a_start, Dplusdot_a_start};

	Minim_Params *parameters = malloc(sizeof(*parameters));
	parameters->lambda0 = lambda_zero;
	parameters->k = k;
	parameters->t0 = t0;
	parameters->redshift_reverse = redshift_reverse;
	parameters->cosmo_params = params;
		
	gsl_odeiv2_system sys = {growth_function, NULL, 2, parameters}; //rk8pd
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, dt, 1e-7, 0.0);		

	t_curr = t_start;
	zs[0] = z_start;
	Dpluses[0] = current_Dplus_values[0];			
	for (int i = 1; i < Dplus_bins; i++) {	
		//the next iteration goes from t_curr to t_next in some number of steps chosen by the driver
		t_next = t_curr + dt;			
		status = gsl_odeiv2_driver_apply(d, &t_curr, t_next, current_Dplus_values);
		z_cur = t_to_z(t_curr, redshift_reverse);	
		zs[i] = z_cur;
		Dpluses[i] = current_Dplus_values[0];
		
		if (status != GSL_SUCCESS) {
			printf ("error, return value=%d\n", status);
			break;
		}	
	}	
	gsl_odeiv2_driver_free(d);
	gsl_integration_workspace_free(z_to_t_workspace);
	free(redshift_reverse);
	free(parameters);	
	return 0;
}

Spline *Dplus_spline(int gravity, double k, BinInfo *z_binInfo, Parameters *params) {

	double *zs = malloc((size_t)Dplus_bins * sizeof(*zs));
	double *Dpluses = malloc((size_t)Dplus_bins * sizeof(*Dpluses));

	Dplus_calc(gravity, k, zs, Dpluses, z_binInfo, params);
	//starts from high z, goes to z = 0
	for (int i = 0; i < Dplus_bins; i++) {
	
		if (Dpluses[i] <= 0.0) {
			printf("Dplus %lf at index %d \n", Dpluses[i], i);
			exit(0);
		}
		
		//printf("z: %lf, Dplus: %lf, test: %lf \n", zs[i], Dpluses[i], params->H0);
		
		//normalise to 1 at z = 0
		Dpluses[i] /= Dpluses[Dplus_bins - 1];
		
		if (zs[i] == 0.0 && Dpluses[i] != 1.0) {
			printf("Dplus not normalized correctly. ind: %d, z: %lf, Dplus: %lf \n", i, zs[i], Dpluses[i]);
			exit(0);
		}
	}
	Spline *toReturn = new_spline(Dplus_bins, zs, Dpluses);	
	free(zs);
	free(Dpluses);
	return toReturn;
}

/*
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
	double lambda_zero, k, omega_m_z, bracket_factor, a_cur, z_cur;
	double R_bar, R0_bar, lambda, fR0;
	double* parameters = (double*) params;
	lambda_zero = parameters[0];
	k = parameters[1];	
	
	//z_cur = splint_generic(redshift_spline_reversed, t);
	//a_cur = z_to_a(z_cur);
	a_cur = t_to_a_analytic(t);
	z_cur = a_to_z(a_cur);
	
	//(7), (17)
	R0_bar = 3.0*H0*H0*(omega_m + 4.0*omega_v);
	R_bar = 3.0*H0*H0*(omega_m*pow(a_cur, -3.0) + 4.0*omega_v);
	lambda = sqrt(pow(lambda_zero,2.0)*pow(R0_bar/R_bar,3.0));
	
	omega_m_z = omega_m*pow(a_cur, -3.0)*pow(H0/H, 2.0);
	omega_m_z = (1.0 - pow(H/H0, 2.0))/(1.0 - pow(a_cur, -3.0));
	H = H_z(z_cur);
	bracket_factor = lambda*lambda*k*k/(a_cur*a_cur);
	//printf("bracket: %le, H: %le, omega_m_Z: %le\n", bracket_factor, H, omega_m_z);
	
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
	gsl_matrix * m = &dfdy_mat.matrix; 
	gsl_matrix_set(m, 0, 0, 0.0);
	gsl_matrix_set(m, 0, 1, 1.0);
	gsl_matrix_set(m, 1, 0, (3.0/2.0)*H*H*omega_m_z*(1.0 + (1.0/3.0)*(bracket_factor)/(1.0 + bracket_factor)));
	gsl_matrix_set(m, 1, 1, -2.0*H);
	dfdt[0] = 0.0;
    dfdt[1] = 0.0;
	return GSL_SUCCESS;
}
*/




