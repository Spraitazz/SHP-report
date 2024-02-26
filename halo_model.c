

double ukm_nfw_profile(double k, double rs, double cdm){
    // For u(k|m) = (1/M)\int_vol \rho(x|m) exp^{-ikx}
    // Spherically symmetric profile, truncated at the virial radius. 
    // u(k|m) = \int_0^{r_{vir}} dr 4pi r^2 [sin(kr)/kr] rho(r|m)/m
    
    double Interim = sin(k*rs)*(gsl_sf_Si((1.0 + cdm)*k*rs) - gsl_sf_Si(k*rs)) - sin(cdm*k*rs)*pow((1.0 + cdm)*k*rs, -1.0) + cos(k*rs)*(gsl_sf_Ci((1.0 + cdm)*k*rs) - gsl_sf_Ci(k*rs));
    
    //double Interim = sin(k)*(gsl_sf_Si((1.0 + cdm)*k) - gsl_sf_Si(k)) - sin(cdm_test*k)*pow((1.0 + cdm)*k*rs, -1.0) + cos(k*rs)*(gsl_sf_Ci((1.0 + cdm)*k*rs) - gsl_sf_Ci(k*rs));
    
    // for the NFW profile, defined by rhos, rs and c, mass is not independent. 
    Interim /= log(1.0 + cdm) - cdm/(1.0 + cdm);
    if (Interim > 1.0 || Interim < 0.0) {
    	printf("NORMALISATION? u(k|M): %lf, k: %lf, rs: %lf, cdm: %lf\n", Interim, k, rs, cdm);
    	exit(0);
    }
    return Interim;
}

double twohalo_term(Spline *Pk_spline, double k, double twohalo_prefactor, double R_nl_sq) {
	return splint_Pk(Pk_spline, k) * volume * exp(-k*k*R_nl_sq) * twohalo_prefactor;
}

int haloModel_realSpace_onlySats(char out[], BinInfo *k_bins, double shot_noise, double twohalo_prefactor, Spline *Pk_spline, double scaling_factor, double R_nl, HOD_Parameters *HOD_params, int avg_sat_no) {

	FILE *f = fopen(out, "w");
	if (out == NULL) {
		printf("wrong file in haloModel_out: %s\n", out);
		exit(0);
	}	
	
	double k, model_Pk, R_nl_sq;
	R_nl_sq = R_nl * R_nl;
	
	for (int i = 0; i < k_bins->bins; i++) {
		k = bin_to_x(k_bins, i);
		model_Pk = shot_noise*pow(ukm_nfw_profile(k, HOD_params->rs, HOD_params->cdm), 2.0) + (volume/(double)avg_sat_no)*(1.0 - pow(ukm_nfw_profile(k, HOD_params->rs, HOD_params->cdm), 2.0));// + twohalo_term(Pk_spline, k*scaling_factor, twohalo_prefactor, R_nl_sq);
		
		fprintf(f, "%le \t %le \t %le \t %le\n", k, model_Pk, model_Pk, model_Pk);
	}
	
	fclose(f);	
	return 0;
}

int haloModel_RSD_onlySats(char out[], BinInfo *k_bins, double shot_noise, double twohalo_prefactor, Spline *Pk_spline, double scaling_factor, double R_nl, HOD_Parameters *HOD_params, int avg_sat_no, Parameters *par) {

	FILE *f = fopen(out, "w");
	if (out == NULL) {
		printf("wrong file in haloModel_out: %s\n", out);
		exit(0);
	}	
	
	double k, model_Pk, R_nl_sq, mono, quadro, one_halo, two_halo, veldisp;
	R_nl_sq = R_nl * R_nl;
	double sigma = 5.0;//sqrt(2.0);
	double fg = pow(par->omega_m_0, par->gamma);
	
	for (int i = 0; i < k_bins->bins; i++) {
		k = bin_to_x(k_bins, i);
		veldisp = (sqrt(pi)/2.0)*gsl_sf_erf(k*sigma)/(k*sigma);
		model_Pk = shot_noise*pow(ukm_nfw_profile(k, HOD_params->rs, HOD_params->cdm), 2.0) + (volume/(double)avg_sat_no)*(1.0 - pow(ukm_nfw_profile(k, HOD_params->rs, HOD_params->cdm), 2.0));// + twohalo_term(Pk_spline, k*scaling_factor, twohalo_prefactor, R_nl_sq);
		mono = model_Pk * veldisp * Kaiser_Monofactor(fg);//kaiserGauss_Monofactor(k*sigma, fg);
		quadro = model_Pk * veldisp * Kaiser_Quadfactor(fg);//kaiserGauss_Quadfactor(k*sigma, fg);
		
		fprintf(f, "%le \t %le \t %le \t %le\n", k, mono, quadro, model_Pk);
	}
	
	fclose(f);	
	return 0;
}

int haloModel_realSpace(char out[], BinInfo *k_bins, double shot_noise, double twohalo_prefactor, Spline *Pk_spline, double scaling_factor, double R_nl, HOD_Parameters *HOD_params) {

	FILE *f = fopen(out, "w");
	if (out == NULL) {
		printf("wrong file in haloModel_out: %s\n", out);
		exit(0);
	} else {
		printf("Printing real-space model to %s\n", out);
	}	
	
	double k = 0.0;
	double model_Pk = 0.0;
	double R_nl_sq = R_nl * R_nl;				
	
	for (int i = 0; i < k_bins->bins; i++) {
		k = bin_to_x(k_bins, i);
		//model_Pk = shot_noise*pow(ukm_nfw_profile(k, HOD_params->rs, HOD_params->cdm), 2.0) + (volume/400000.0);//*(1.0 - pow(ukm_nfw_profile(k, HOD_params->rs, HOD_params->cdm), 2.0));// + twohalo_term(Pk_spline, k*scaling_factor, twohalo_prefactor, R_nl_sq);
		model_Pk = shot_noise + twohalo_term(Pk_spline, k*scaling_factor, twohalo_prefactor, R_nl_sq);
		fprintf(f, "%le \t %le", k, model_Pk);
		fprintf(f, "\n");
	}
	
	fclose(f);	
	return 0;
}

//delete me
int haloModel_realSpace_HODsats(char out[], BinInfo *k_bins, double shot_noise, double twohalo_prefactor, Spline *Pk_spline, double scaling_factor, double R_nl, HOD_Parameters *HOD_params) {

	FILE *f = fopen(out, "w");
	if (out == NULL) {
		printf("wrong file in haloModel_out: %s\n", out);
		exit(0);
	} else {
		printf("Printing real-space model to %s\n", out);
	}	
	
	double k = 0.0;
	double model_Pk = 0.0;
	double R_nl_sq = R_nl * R_nl;				
	
	for (int i = 0; i < k_bins->bins; i++) {
		k = bin_to_x(k_bins, i);
		//model_Pk = shot_noise*pow(ukm_nfw_profile(k, HOD_params->rs, HOD_params->cdm), 2.0) + (volume/400000.0);//*(1.0 - pow(ukm_nfw_profile(k, HOD_params->rs, HOD_params->cdm), 2.0));// + twohalo_term(Pk_spline, k*scaling_factor, twohalo_prefactor, R_nl_sq);
		model_Pk = shot_noise + twohalo_term(Pk_spline, k*scaling_factor, twohalo_prefactor, R_nl_sq);
		fprintf(f, "%le \t %le", k, model_Pk);
		fprintf(f, "\n");
	}
	
	fclose(f);	
	return 0;
}

double f_v_ukm_ST(double v, void *params) {
	return f_v(v) * b_v(v); 

}

//EXPERIMENTAL INTEGRAL BEWARE. FOR b(v) TO GET THE shotnoise*ukm profile term (one-halo)
double oneHalo_int(double vmin, double vmax, Spline *variance, Parameters *params) {
	double result, error;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	
	Bias_Params *send = malloc(sizeof(*send));
	send->parameters = params;
	send->variance = variance;
	send->variance_reversed = reverse_spline(variance);
	
	gsl_function F;	
	F.function = &f_v_ukm_ST;
	F.params = send;        

    gsl_integration_qags(&F, vmin, vmax, 0, 1e-7, 1000, w, &result, &error);	
	gsl_integration_workspace_free(w);


	return result;
}
/*
//just outputs the halo model power spectrum to a given file
int haloModel_out(char out[], BinInfo *k_bins, double shot_noise, double twohalo_prefactor, int redshift_space, Spline *Pk_spline, Parameters *par, double scaling_factor, double mass_bias, double R_nl, HOD_Parameters *HOD_params) {
	
	FILE *f = fopen(out, "w");
	if (out == NULL) {
		printf("wrong file in haloModel_out: %s \n", out);
		exit(0);
	}	
	
	double k, model_Pk, model_mono, model_quad, model_hex, fg, R_nl_sq, twohalo_Pk, beta;
	fg = pow(par->omega_m_0, par->gamma);
	beta = fg/mass_bias;	
	R_nl_sq = R_nl * R_nl;
	
	if (redshift_space == REAL_SPACE) {
	
		for (int i = 0; i < k_bins->bins; i++) {
			k = bin_to_x(k_bins, i);
			model_Pk = shot_noise*pow(ukm_nfw_profile(k, rs, cdm),2.0) + twohalo_term(Pk_spline, k*scaling_factor, twohalo_prefactor, R_nl_sq);
			fprintf(f, "%le \t %le \t %le \t %le \n", k, model_Pk, model_Pk, model_Pk);
		}	
					
	} else if (redshift_space == REDSHIFT_SPACE) {
	
		for (int i = 0; i < k_bins->bins; i++) {
			k = bin_to_x(k_bins, i);
			twohalo_Pk = twohalo_term(Pk_spline, k*scaling_factor, twohalo_prefactor, R_nl_sq);
			model_mono = shot_noise + twohalo_Pk * KaiserLorentz_Monofactor(0.0, beta);
			model_quad = twohalo_Pk * KaiserLorentz_Quadfactor(0.0, beta);
			model_hex = twohalo_Pk * KaiserLorentz_Hexfactor(0.0, beta);
			//fprintf(f, "%le \t %le \t %le \n", k, model_mono, model_quad);
			fprintf(f, "%le \t %le \t %le \t %le", k, model_mono, model_quad, model_hex);
			fprintf(f, "\n");
		}
				
	} else {
		printf("can only have REDSHIFT_SPACE or REAL_SPACE in halo model\n");
	}	
	
	fclose(f);	
	return 0;
}
*/

double I_integral_HOD(double t, void *params) {
	params = params;
	double f_t = log(1.0 + t) - t/(1.0 + t);
	return f_t/(pow(t,3.0)*pow(1.0 + t,2.0));
}

/*
int print_cdm(char outFile[]) {

	Parameters *cur_params_z0 = malloc(sizeof(*cur_params_z0));
	cur_params_z0->H0 = 67.8;
	cur_params_z0->omega_m_0 = 0.25;
	cur_params_z0->omega_v_0 = 0.75;
	cur_params_z0->omega_r_0 = 0.0;
	cur_params_z0->omega = 0.0;
	cur_params_z0->gamma = 0.545;
	cur_params_z0->z = 0.0;	
	
	Parameters *cur_params_z06 = malloc(sizeof(*cur_params_z06));
	cur_params_z06->H0 = 67.8;
	cur_params_z06->omega_m_0 = 0.25;
	cur_params_z06->omega_v_0 = 0.75;
	cur_params_z06->omega_r_0 = 0.0;
	cur_params_z06->omega = 0.0;
	cur_params_z06->gamma = 0.545;
	cur_params_z06->z = 0.6;

	char model1[100];
	sprintf(model1, "%s/data/models/matter_pk_om_m_025_z_0_sigma8_06.dat", home_directory);
	SplineInfo* Pk_model1_z06 = input_spline_file(model1, Pk_model_format, NORMAL);
	
	char model2[100];
	sprintf(model2, "%s/data/models/matter_pk_om_m_015_z_0_sigma8_1.dat", home_directory);
	SplineInfo* Pk_model2_z06 = input_spline_file(model2, Pk_model_format, NORMAL);
	
	Parameters *targ_params_z0 = malloc(sizeof(*targ_params_z0));
	targ_params_z0->H0 = 67.8;
	targ_params_z0->omega_m_0 = 0.15;
	targ_params_z0->omega_v_0 = 0.85;
	targ_params_z0->omega_r_0 = 0.0;
	targ_params_z0->omega = 0.0;
	targ_params_z0->gamma = 0.545;
	targ_params_z0->z = 0.0;	
	
	Parameters *targ_params_z06 = malloc(sizeof(*targ_params_z0));
	targ_params_z06->H0 = 67.8;
	targ_params_z06->omega_m_0 = 0.15;
	targ_params_z06->omega_v_0 = 0.85;
	targ_params_z06->omega_r_0 = 0.0;
	targ_params_z06->omega = 0.0;
	targ_params_z06->gamma = 0.545;
	targ_params_z06->z = 0.6;	
	
	
	//binning
	k_bins = 500;	
	Dplus_bins = 500;
	z_bins = 500;
	k_min = 0.01;
	k_max = 1.0;	
	z_min = 0.0;
	z_max = 3.0;
	double R_min = 0.01;
	double R_max = 5.0;
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	rescaling_R_bin_info = prep_bins(R_min, R_max, 500, LOG_BIN);
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);

	z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	OLV_workspace = gsl_integration_workspace_alloc(OLV_workspace_size);

	prep_redshift_spline(rescaling_z_bin_info);

	//at z = 1.2
	SplineInfo *Pk_model1_z12 = prep_Pk_constz(GR, 0.0, 0.6, Pk_model1_z06, rescaling_k_bin_info);	

	SplineInfo_Extended *Pk_model1_z06_ext = extend_spline_model(Pk_model1_z06);
	SplineInfo_Extended *Pk_model1_z12_ext = extend_spline_model(Pk_model1_z12);
	
	SplineInfo *Pk_model2_z12 = prep_Pk_constz(GR, 0.0, 0.6, Pk_model2_z06, rescaling_k_bin_info);	
	SplineInfo_Extended *Pk_model2_z06_ext = extend_spline_model(Pk_model2_z06);
	SplineInfo_Extended *Pk_model2_z12_ext = extend_spline_model(Pk_model2_z12);	

	
	SplineInfo *variance_model1_z06 = prep_variances(rescaling_R_bin_info, Pk_model1_z06_ext);
	SplineInfo *variance_model1_z12 = prep_variances(rescaling_R_bin_info, Pk_model1_z12_ext);
	
	SplineInfo *variance_model2_z06 = prep_variances(rescaling_R_bin_info, Pk_model2_z06_ext);
	SplineInfo *variance_model2_z12 = prep_variances(rescaling_R_bin_info, Pk_model2_z12_ext);
	//print_spline(variance_model2_z12, rescaling_R_bin_info);	
		
	Mmin = 5e12;
	Mmax = 1e15;
	int bins = 200;
	FILE *f = fopen(outFile, "w");
	double dm = (Mmax - Mmin)/(double)(bins-1);
	//model1 0.6, model2 0.6, model1 1.2, model2 1.2
	double m, cdm1, cdm2, cdm3, cdm4;
	for (int i = 0; i < bins; i++) {
		m = Mmin + (double)i * dm;

		set_halo_params(m, variance_model1_z06, cur_params_z0);
		cdm1 = cdm_global;		

		set_halo_params(m, variance_model2_z06, targ_params_z0);
		cdm2 = cdm_global;
		
		set_halo_params(m, variance_model1_z12, cur_params_z06);
		cdm3 = cdm_global;		

		set_halo_params(m, variance_model2_z12, targ_params_z06);
		cdm4 = cdm_global;
		
		fprintf(f, "%le \t %le \t %le \t %le \t %le \n", m, cdm1, cdm2, cdm3, cdm4);
	
	}
	fclose(f);
	return 0;

}
*/


