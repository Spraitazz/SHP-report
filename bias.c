//http://arxiv.org/pdf/1308.5183v3.pdf
double f_v(double v) {
	double A, q, p;
	A = 0.216;
	q = 0.707;
	p = 0.3;
	return A*(1.0 + 1.0/pow(q*v*v, p))*exp(-q*v*v/2.0);
}

double f_v_over_m(double v, void *params) {
	double A, q, p, fv, m;
	Bias_Params *bias_params = (Bias_Params*) params;
	A = 0.216;
	q = 0.707;
	p = 0.3;
	fv = A*(1.0 + 1.0/pow(q*v*v, p))*exp(-q*v*v/2.0);
	m = v_to_m(v, bias_params);
	if (m < 1.0) printf("m: %le in fv over m \n", m);
	return fv / m;
}

double ln_fv_deriv(double v) {
	double q = 0.707;
	double p = 0.3;
	double ans = (1.0/(1.0 + 1.0/pow(q*v*v, p)))*(1.0/pow(q,p))*(-2.0*p)*pow(v, -2.0*p - 1.0) - q*v;
	//double ans2 = (-2.0*p/v - q*v*(pow(q, p) * pow(v, 2.0*p) + 1.0)) / (pow(q, p) * pow(v, 2.0*p) + 1.0);
	//THIS HAPPENED, CHECK WHAT IS UP
	if ((-1e-6 < ans) && (ans < 1e-6)) {
		printf("ln fv derivative ?? in bias.c; ans: %le, v: %le \n", ans, v);
		exit(0); 
	} 
	return ans;
}

double m_to_v(double mass, Parameters *params, Spline* variance_spline) {
	double R = m_to_R(mass, params);
	double sigma_R = sqrt(splint_generic(variance_spline, R));
	if (sigma_R < 1e-6 || sigma_R > 1e6 || isnan(sigma_R) || (sigma_R != sigma_R)) {
		printf("m: %le, R: %lf, sigma: %lf in m to v (bias.c)\n", mass, R, sigma_R);
		exit(0);
	}	
	double v = delta_c/sigma_R;
	return v;
}

double v_to_m(double v, Bias_Params *bias_params) {
	double sigma_R = delta_c / v;
	double R = splint_generic(bias_params->variance_reversed, sigma_R*sigma_R);
	double mass = R_to_m(R, bias_params->parameters);
	if (R < 1e-6 || R > 1e6 || isnan(R) || (R != R)) {
		printf("m: %le, R: %lf, sigma: %lf in v to m\n", mass, R, sigma_R);
		exit(0);
	}	
	return mass;
}

//http://arxiv.org/pdf/1308.5183v3.pdf
double b_v(double v) {
	return 1.0 - 1.0/delta_c - (v/delta_c)*ln_fv_deriv(v);
}

double b_m(double m, Parameters *params, Spline *variance_spline) {
	return b_v(m_to_v(m, params, variance_spline));
}

double bf_over_m(double v, void *params) {
	Bias_Params *parameters = (Bias_Params*) params;
	return b_v(v)*f_v(v)/v_to_m(v, parameters);
}

double f_over_m(double v, void *params) {
	Bias_Params *parameters = (Bias_Params*) params;
	return f_v(v)/v_to_m(v, parameters);
}

double fv_dv_int(double vmin, double vmax, Parameters *params, Spline *variance) {

	double result, error;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	
	Bias_Params *send = malloc(sizeof(*send));
	send->parameters = params;
	send->variance = variance;
	send->variance_reversed = reversed_spline(variance);
	
	gsl_function F;	
	F.function = &f_v_over_m;
	F.params = send;        

    gsl_integration_qags(&F, vmin, vmax, 0, 1e-7, 1000, w, &result, &error);	
	gsl_integration_workspace_free(w);
	return result;
}

double calc_b_eff(Parameters *params, Spline* variance_spline, double Mmin, double Mmax) {
		
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	double error_top, result_top, error_bot, result_bot, result_fv, error_fv;	
	double vmin = m_to_v(Mmin, params, variance_spline);
	double vmax = m_to_v(Mmax, params, variance_spline);		
		
	Bias_Params *send = malloc(sizeof(*send));
	send->parameters = params;
	send->variance = variance_spline;
	send->variance_reversed = reversed_spline(variance_spline);	

	gsl_function F;	
	F.function = &bf_over_m;
	F.params = send;        

    gsl_integration_qags(&F, vmin, vmax, 0, 1e-7, 1000, w, &result_top, &error_top);		
	
	F.function = &f_over_m;
	
	gsl_integration_qags(&F, vmin, vmax, 0, 1e-7, 1000, w, &result_bot, &error_bot);	

	double b_eff = result_top / result_bot;
	
	free(send);
	gsl_integration_workspace_free(w);	
	return b_eff;
}

/*
double b_eff(Pk_Spline *pk_spline, Parameters *params, double *mass_lims) {
	double R1 = m_to_R(mass_lims[0], params);
	double R2 = m_to_R(mass_lims[1], params);
	BinInfo *R_binInfo_beff = prep_bins(R1, R2, 300, LOG_BIN);
	Spline *variance_beff = prep_variances(R_binInfo_beff, pk_spline);	
	double toReturn = calc_b_eff(params, variance_beff, mass_lims[0], mass_lims[1]);		
	free(R_binInfo_beff);
	free_spline(variance_beff);
	return toReturn;
}
*/

