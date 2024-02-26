int make_mocks() {

	char Pk_current_path[200];
	sprintf(Pk_current_path, "%s/data/models/linear_matter_pk_om_m_0.31_z_0.0_sigma8_0.83_h_0.7.dat", home_directory);
	Pk_Spline *Pk_current = input_Pk_spline_file(Pk_current_path, Pk_model_format, 0.0);
	
	
	Parameters *cur_params = malloc(sizeof(*cur_params));		
	cur_params->H0 = 100.0*0.7;
	cur_params->omega_m_0 = 0.31;
	cur_params->omega_v_0 = 0.69;
	cur_params->omega_r_0 = 0.0;
	cur_params->omega = 1.0;	
	cur_params->gamma = 0.545;
	cur_params->z = 2.5; //actually 0.0 to start off with, check filename
	
	int particle_no = 300000;
	double shotnoise = volume/particle_no;
	double Mmin_cur = 5e11;
	double Mmax_cur = 5e14;
	double *mass_lims = malloc(2 * sizeof(*mass_lims));	
	mass_lims[0] = Mmin_cur;
	mass_lims[1] = Mmax_cur;
	double R1 = m_to_R(Mmin_cur, cur_params);
	double R2 = m_to_R(Mmax_cur, cur_params);
	printf("\nM min of current: %le M_sun, M max of current: %le M_sun\n", Mmin_cur,  Mmax_cur);	
	printf("(current catalogue) R1: %lf Mpc H^-1, R2: %lf Mpc H^-1\n", R1, R2);	
	

	int z_bins = 300;
	Dplus_bins = 300;	
	double z_min = 0.0;
	double z_max = 4.0;	
	BinInfo *z_binInfo = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);	
	Spline *Dplus_z = Dplus_spline(GR, 0.0, z_binInfo, cur_params);
	free(z_binInfo);	

	Pk_set_redshift_GR(Pk_current, Dplus_z, 2.5);	
	BinInfo *R_binInfo = prep_bins(R1, R2, 300, LOG_BIN);
	Spline *variance = prep_variances(R_binInfo, Pk_current);	

	
	double b_eff = calc_b_eff(cur_params, variance, Mmin_cur, Mmax_cur);
	printf("beff in make_mocks(): %lf\n", b_eff);
	
	//with disp bias + shotnoise
	Pk_Spline *Pk_current_model = copy_Pk_spline(Pk_current);
	
	scale_spline(Pk_current_model->spline, b_eff*b_eff);
	spline_add_const(Pk_current_model->spline, shotnoise);
	//print_Pk_spline(Pk_current_model, NULL);
	exit(0);
	
	BinInfo *Pk_model_bins = prep_bins(0.05, 0.5, 200, LOG_BIN);
	char model_out[200];
	sprintf(model_out, "%s/data/mocks_new/Pk/models/Pk_model_realspace_z2.5.dat", home_directory);
	print_Pk_spline_file(Pk_current_model, Pk_model_bins, model_out);
	
	//extend_spline(Pk_current);
	char mocks_out[200];
	sprintf(mocks_out, "%s/data/mocks_new/mocks/ZA_only_mock_z2.5_%s.dat", home_directory, "%d");
	ZA_mocks(mocks_out, 10, particle_no, Pk_current, NULL, cur_params, mass_lims);

	return 0;
}


int measure_Pk() {

	char mockPath[200], Pk_out[200];
	Catalogue *current_catalogue = NULL;
	
	for (int i = 0; i < 10; i++) {
	
		sprintf(mockPath, "%s/data/mocks_new/mocks/ZA_only_mock_z2.5_%d.dat", home_directory, i);
		current_catalogue = input_catalogue_file(mockPath, 0, mocks_format_withIds);
		
		sprintf(Pk_out, "%s/data/mocks_new/Pk/measured/Pk_measured_realspace_z2.5_%d.dat", home_directory, i);
		haloes_measure_Pk(current_catalogue, Pk_out, 0.0, 1.0, CIC, 1.0);	
		
		//toRedshift(current_catalogue, cur_params);
		free_catalogue(current_catalogue);
	}

	return 0;
}


