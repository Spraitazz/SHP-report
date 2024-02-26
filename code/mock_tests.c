int ZA_mocks_test(char outPath[], int mock_no, int particle_no, Spline *Pk_model, Spline *variance_HOD, Parameters *params, double mass_lims[]) {
	
	printf("\npreparing ZA MOCKS \n");	
	
	double Mmin = mass_lims[0];
	double Mmax = mass_lims[1];
	double R1 = m_to_R(Mmin, params);
	double R2 = m_to_R(Mmax, params);
	BinInfo *R_binInfo = prep_bins(R1, R2, 300, LOG_BIN);
	Spline *variance = prep_variances(R_binInfo, Pk_model);
	
	
	//prepare cumsums from fvdv/m integral
	double m, v;
	int massBins = 300;
	BinInfo *mass_binInfo = prep_bins(Mmin, Mmax, massBins, LOG_BIN);
	
	double vmin = m_to_v(Mmin, params, variance);
	double *ms = malloc((size_t)massBins * sizeof(*ms));
	double *cumsums = malloc((size_t)massBins * sizeof(*cumsums));
	for (int i = 0; i < massBins; i++) {
		m = bin_to_x(mass_binInfo, i);		
		v = m_to_v(m, params, variance);
		ms[i] = m;
		cumsums[i] = fv_dv_int(vmin, v, params, variance);
		
	}
	
	for (int i = 0; i < massBins; i++) {
		cumsums[i] /= cumsums[massBins-1];	
	}
	
	//reversed spline needed
	Spline *fv_cumsums = input_spline_values(massBins, ms, cumsums, MEASURED);
	Spline *fv_cumsums_reverse = reverse_spline(fv_cumsums);
	free(ms);
	free(cumsums);
	
	//effective bias, should match the one given
	double b_eff = calc_b_eff(params, variance, Mmin, Mmax);	
	printf("b eff from integral (in mocks test): %lf \n\n", b_eff);		

	//displacement storage, catalogue
	double **disps_init = NULL;
	calloc2D_double(&disps_init, particle_no, 3);
	
	Catalogue *current_catalogue = NULL;
	Catalogue *satellites = NULL;
	Catalogue *centrals = NULL;
	Catalogue *final_catalogue = NULL;
	char cur_outPath[200];
	char Pk_current_out[200], Pk_current_4fold_out[200];
	char Pk_current_centrals_out[200];
	char Pk_current_sats_out[200], Pk_current_sats_4fold_out[200];

	double rand_cumsum, cur_mass;
	double shotnoise_sats = 0.0;
	double shotnoise_all = 0.0;	
	for (int i = 0; i < mock_no; i++) {				

		current_catalogue = populateHaloes_oneMass(particle_no, 5e13);		
		
		//set ST masses
		//catalogue_setMasses(current_catalogue, fv_cumsums_reverse);		
				
		//generate ZA disp field
		ZA_displacements(current_catalogue, &disps_init, Pk_model, NULL, false, 0.0);
		
		//large scale velocity field not biased tcdm pg. 3
		apply_ZA_velocities(current_catalogue, &disps_init, params);
		
		apply_ZA_displacements(current_catalogue, &disps_init);	
		
		centrals = HOD_centrals(current_catalogue);
		
		//sprintf(Pk_current_centrals_out, "%s/data/rescaling/model_to_model/HOD_runs/Pk_current_centrals_%d.dat", home_directory, i);
		//haloes_measure_Pk(centrals, Pk_current_centrals_out, 0.0, 1.0, NGP, 1.0);	
		
		satellites = HOD_satellites(current_catalogue, variance_HOD, params);
		shotnoise_sats = 0.0;//volume/(double)satellites->particle_no;
		
		sprintf(Pk_current_sats_out, "%s/data/rescaling/model_to_model/HOD3_simple_vels/Pk_current_sats_velDisp_%d.dat", home_directory, i);
		sprintf(Pk_current_sats_4fold_out, "%s/data/rescaling/model_to_model/HOD3_simple_vels/Pk_current_sats_velDisp_4fold_%d.dat", home_directory, i);
		toRedshift(satellites, params);
		haloes_measure_Pk(satellites, Pk_current_sats_out, shotnoise_sats, 1.0, NGP, 1.0);
		haloes_measure_Pk(satellites, Pk_current_sats_4fold_out, shotnoise_sats, 4.0, NGP, 1.0);	
		
		free_catalogue(current_catalogue);		
		//final_catalogue = combine_catalogues(centrals, satellites);		
		free_catalogue(centrals);
		free_catalogue(satellites);

		//current_catalogue = final_catalogue;
		//shotnoise_all = 0.0;//volume/(double)current_catalogue->particle_no;
		
		//sprintf(Pk_current_out, "%s/data/rescaling/model_to_model/HOD_runs/Pk_current_%d.dat", home_directory, i);
		//sprintf(Pk_current_4fold_out, "%s/data/rescaling/model_to_model/HOD_runs/Pk_current_4fold_%d.dat", home_directory, i);
		//haloes_measure_Pk(current_catalogue, Pk_current_out, shotnoise_all, 1.0, NGP, 1.0);
		//haloes_measure_Pk(current_catalogue, Pk_current_4fold_out, shotnoise_all, 4.0, NGP, 1.0);
		
		//sprintf(cur_outPath, outPath, i);	
		//catalogue_to_file(current_catalogue, cur_outPath);
		//free_catalogue(current_catalogue);
	}
	
	//free2D_double(&disps_init, particle_no);	
	return 0;
}

int ZA_mocks_HOD_test() {

	char Pk_current_path[100];
	sprintf(Pk_current_path, "%s/data/models/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	Spline *Pk_current_model = input_spline_file(Pk_current_path, Pk_model_format, MODEL);	
	//this will go from redshift 0.75 to redshift 3.5 later	
		
	Parameters *cur_params = malloc(sizeof(*cur_params));		
	cur_params->H0 = 100.0*0.678;
	cur_params->omega_m_0 = 0.24;
	cur_params->omega_v_0 = 0.76;
	cur_params->omega_r_0 = 0.0;
	cur_params->omega = 1.0;	
	cur_params->gamma = 0.545;
	cur_params->z = 3.5; //actually 0.75 to start off with, check filename	
	
	int particle_no = 100000;
	double Mmin_cur = 5e11;
	double Mmax_cur = 1e15;
	double R1 = m_to_R(Mmin_cur, cur_params);
	double R2 = m_to_R(Mmax_cur, cur_params);
	printf("\nM min of current: %le M_sun, M max of current: %le M_sun\n", Mmin_cur,  Mmax_cur);	
	printf("(current catalogue) R1: %lf Mpc H^-1, R2: %lf Mpc H^-1\n", R1, R2);	
	
	//--------------------------------------------------------
	
	//find optimal s, z		
	bool vary_z_current = true;	
	int k_bins = 300;
	int z_bins = 300;
	Dplus_bins = 300;
	int R_bins = 300;
	double k_min = 0.05; //k_fund
	double k_max = 3.0; //k_nyq/4
	double z_min = 0.0;
	double z_max = 4.0;	
	int current_gravity = GR;
	int target_gravity = GR;	
	

	BinInfo *R_binInfo = prep_bins(R1, R2, R_bins, LOG_BIN);
	BinInfo *k_binInfo = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	BinInfo *z_binInfo = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);	

		
	Spline *Dplus_z_cur = Dplus_spline(GR, 0.0, z_binInfo, cur_params);
	
		
	//just changing current model from z = 0.75 to z = 3.5. //Prep Pk function keeps spline model
	BinInfo *k_binInfo_curModel = prep_bins(Pk_current_model->splineInfo->xmin, Pk_current_model->splineInfo->xmax, Pk_current_model->splineInfo->lines, LOG_BIN);
	Spline *Pk_current_model_z0 = prep_Pk_constz(GR, Dplus_z_cur, 0.75, 0.0, Pk_current_model, k_binInfo_curModel, z_binInfo, cur_params);
	Pk_current_model = prep_Pk_constz(GR, Dplus_z_cur, 0.75, cur_params->z, Pk_current_model, k_binInfo_curModel, z_binInfo, cur_params);	
	free(k_binInfo_curModel);
	
	Spline *variance_HOD = prep_variances(R_binInfo, Pk_current_model_z0);
	
	

	
	
	
	//exit(0);
	
	double mass_lims_mocks[] = {Mmin_cur, Mmax_cur};
	char mocks_out[100];
	sprintf(mocks_out, "%s/data/ZA_mocks_HOD/ZA_mock_catalogue_%s.dat", home_directory, "%d");	
	ZA_mocks_test(mocks_out, 10, particle_no, Pk_current_model, variance_HOD, cur_params, mass_lims_mocks);
	//exit(0);	
	
	double weighted_SN = 194.035569;
	int final_parts = 639362;

	//HOD4_ST_realspace_sats
	char Pk_sats_model[200];
	sprintf(Pk_sats_model, "%s/data/rescaling/model_to_model/HOD3_simple_vels/Pk_sats_model.dat", home_directory);
	//haloModel_RSD_onlySats(Pk_sats_model, k_binInfo, weighted_SN, 1.0, Pk_current_model, 1.0, 0.0, HOD_parameters, final_parts, cur_params);
	
		
	//weighted SN - volume/sat no
	sat_no_tot /= 10;
	weighted_shotnoise_tot /= 10.0;

	printf("average no of sats: %d, average weighted shotnoise: %lf \n", sat_no_tot, weighted_shotnoise_tot);
	
	
	
	
	exit(0);



	return 0;
}

