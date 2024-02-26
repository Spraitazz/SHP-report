double twohalo_term(Spline *Pk_spline, double k, double twohalo_prefactor, double R_nl_sq) {
	return splint_Pk(Pk_spline, k) * volume * exp(-k*k*R_nl_sq) * twohalo_prefactor;
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
