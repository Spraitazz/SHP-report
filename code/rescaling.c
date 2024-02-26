int final_rescaling() {
	
	char Pk_current_path[200];
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
	
	int particle_no_cur = 300000;
	double shotnoise_cur = volume/particle_no_cur;
	double Mmin_cur = 5e11;
	double Mmax_cur = 5e14;
	double *mass_lims_cur = malloc(2 * sizeof(*mass_lims_cur));	
	mass_lims_cur[0] = Mmin_cur;
	mass_lims_cur[1] = Mmax_cur;
	double R1 = m_to_R(Mmin_cur, cur_params);
	double R2 = m_to_R(Mmax_cur, cur_params);
	printf("\nM min of current: %le M_sun, M max of current: %le M_sun\n", Mmin_cur,  Mmax_cur);	
	printf("(current catalogue) R1: %lf Mpc H^-1, R2: %lf Mpc H^-1\n", R1, R2);	
	
	//--------------------------------------------------------
	
	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/models/matter_pk_om_m_0.15_z_1.0_sigma8_0.69.dat", home_directory);
	Spline *Pk_target_model = input_spline_file(Pk_target_path, Pk_model_format, MODEL);
		
	Parameters *targ_params = malloc(sizeof(*targ_params));
	targ_params->H0 = 100.0*0.678;
	targ_params->omega_m_0 = 0.15;
	targ_params->omega_v_0 = 0.85;
	targ_params->omega_r_0 = 0.0;
	targ_params->omega = 1.0;	
	targ_params->gamma = 0.545;
	targ_params->z = 1.0;	
	
	int particle_no_targ = 300000;
	double shotnoise_targ = volume/particle_no_targ;	
	double R_nl = 0.0; //from linear Pk
	double Mmin_targ = 1e12;
	double Mmax_targ = 1e15;
	double *mass_lims_targ = malloc(2 * sizeof(*mass_lims_targ));	
	mass_lims_targ[0] = Mmin_targ;
	mass_lims_targ[1] = Mmax_targ;
	double R1_primed = m_to_R(Mmin_targ, targ_params);
	double R2_primed = m_to_R(Mmax_targ, targ_params);
	double *R_lims = malloc(2 * sizeof(*R_lims));
	R_lims[0] = R1_primed;
	R_lims[1] = R2_primed;
	printf("\nM min of target: %le M_sun, M max of target: %le M_sun\n", Mmin_targ,  Mmax_targ);
	printf("(target catalogue) R1': %lf Mpc H^-1, R2': %lf Mpc H^-1\n", R1_primed, R2_primed);	
	
	//---------------------------------------------------------------------

	//find optimal s, z		
	bool vary_z_current = true;	
	int k_bins = 300;
	int z_bins = 300;
	Dplus_bins = 300;
	int R_bins = 300;
	double k_min = 0.05; //k_fund
	double k_max = 0.5; //k_nyq/4
	double z_min = 0.0;
	double z_max = 4.0;	
	int current_gravity = GR;
	int target_gravity = GR;		

	BinInfo *k_binInfo = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	BinInfo *z_binInfo = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	BinInfo *R_binInfo = prep_bins(0.4*R1_primed, 2.5*R2_primed, R_bins, LOG_BIN);		
	
	//D+(z), GR only
	Spline *Dplus_z_cur;
	Spline *Dplus_z_targ;
	//MAKE PREP_PK_CONSTZ USE ARRAY OF SPLINES, in case of modified grav varying k
	if (current_gravity == GR && target_gravity == GR) {
		Dplus_z_cur = Dplus_spline(GR, 0.0, z_binInfo, cur_params);
		Dplus_z_targ = Dplus_spline(GR, 0.0, z_binInfo, targ_params);	
	}
	
	//just changing current model from z = 0.75 to z = 3.5. //Prep Pk function keeps spline model
	BinInfo *k_binInfo_curModel = prep_bins(Pk_current_model->splineInfo->xmin, Pk_current_model->splineInfo->xmax, Pk_current_model->splineInfo->lines, LOG_BIN);
	Pk_current_model = prep_Pk_constz(GR, Dplus_z_cur, 0.75, cur_params->z, Pk_current_model, k_binInfo_curModel, z_binInfo, cur_params);	
	
	double mass_lims_mocks[] = {Mmin_cur, Mmax_cur};
	char mocks_out[200];
	sprintf(mocks_out, "%s/data/final/mocks_before/DM/original_mock_%s.dat", home_directory, "%d");
	//ZA_mocks(mocks_out, 10, particle_no, Pk_current_model, NULL, cur_params, mass_lims_mocks);	
	
	double b_eff_cur = b_eff(Pk_current_model, cur_params, mass_lims_cur);		
	printf("b eff of current catalogue: %lf \n", b_eff_cur);	
	double b_eff_targ = b_eff(Pk_target_model, targ_params, mass_lims_targ);	
	printf("b eff of target catalogue: %lf \n", b_eff_targ);
	
	Spline *Pk_current_model_beff = scale_spline(Pk_current_model, b_eff_cur*b_eff_cur);	
	Spline *Pk_current_model_mocks = add_const(Pk_current_model_beff, shotnoise_cur);
	
	Spline *Pk_target_model_beff = scale_spline(Pk_target_model, b_eff_targ*b_eff_targ);	
	Spline *Pk_target_model_mocks = add_const(Pk_target_model_beff, shotnoise_targ);

		
	/*
	
	// variance target with k-limits for minimisation
	Spline *Pk_target_model_klims = prep_Pk_constz(target_gravity, Dplus_z_cur, targ_params->z, targ_params->z, Pk_target_model_mocks, k_binInfo, z_binInfo, targ_params);
	Pk_target_model_klims->model = MEASURED;
	Spline *variance_target = prep_variances(R_binInfo, Pk_target_model_klims);
	

	
	//for current simulation, varying z:		
	Spline **Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	Spline **variance_splines_zBins = malloc((size_t)z_bins * sizeof(*variance_splines_zBins));	

	//finding a redshift to start minimising from
	double z_to = 0.0;
	double z_guess = 0.0;
	double sq_diff = 0.0;
	double k_chk = 0.0;
	double min_sq_diff = DBL_MAX;	
	BinInfo *zmin_k_bins = prep_bins(k_min, k_max, 20, LOG_BIN);
	for (int i = 0; i < z_bins; i++) {	
		z_to = bin_to_x(z_binInfo, i);
		printf("preparing Pk at z: %lf, z max: %lf \n", z_to, z_max);
	
		//current Pk at different redshifts at the right k limits
		Pk_splines_zBins[i] = prep_Pk_constz(current_gravity, Dplus_z_cur, cur_params->z, z_to, Pk_current_model_mocks, k_binInfo, z_binInfo, cur_params);
		Pk_splines_zBins[i]->model = MEASURED;	
		variance_splines_zBins[i] = prep_variances(R_binInfo, Pk_splines_zBins[i]);
	
		sq_diff = 0.0;
		for (int j = 0; j < zmin_k_bins->bins; j++) {
			k_chk = bin_to_x(zmin_k_bins, j);
			sq_diff += pow(splint_Pk(Pk_splines_zBins[i], k_chk) - splint_Pk(Pk_target_model, k_chk), 2.0);
			if (sq_diff < min_sq_diff) {
				min_sq_diff = sq_diff;
				z_guess = z_to;
			}	
		}		
	}	
	printf("z guess: %lf \n", z_guess);	

	//minimization
	dsq_multimin(vary_z_current, z_guess, variance_target, variance_splines_zBins, z_binInfo, R_lims);	
	exit(0);
	*/
	
	
	//did once, found those, its enough LINEAR ONLY
	//s = 1.091418;
	//z_rescaled = 0.662463; //dsq at s, z: 2.254383e-05
	
	//linear*beff*beff + SN
	s = 1.123937;
	z_rescaled = 1.848116; //dsq at s, z: 2.066408e-05
	
	
	Parameters *cur_params_rescaled = malloc(sizeof(*cur_params));		
	cur_params_rescaled->H0 = 100.0*0.678;
	cur_params_rescaled->omega_m_0 = 0.15;
	cur_params_rescaled->omega_v_0 = 0.85;
	cur_params_rescaled->omega_r_0 = 0.0;
	cur_params_rescaled->omega = 1.0;	
	cur_params_rescaled->gamma = 0.545;
	cur_params_rescaled->z = z_rescaled;
	
	//prefactor to change redshift
	double Dplus_cur = splint_generic(Dplus_z_cur, cur_params->z);
	double Dplus_after = splint_generic(Dplus_z_cur, z_rescaled);
	double z_scalefac = Dplus_after/Dplus_cur;	
	
	
	char model_out[200];
	sprintf(model_out, "%s/data/final/mocks_before/Pk/current_model.dat", home_directory);
	print_spline_file(Pk_current_model_mocks, k_binInfo, model_out);
		
	Spline *Pk_current_model_s_scaled = scale_x(Pk_current_model_mocks, s);	
	Spline *Pk_current_model_zs_scaled = scale_spline(Pk_current_model_s_scaled, pow(s, 3.0)*pow(z_scalefac, 2.0));
	Pk_current_model_zs_scaled->model = MEASURED;
	//Spline *Pk_current_model_scaled_beff = scale_spline(Pk_current_model_zs_scaled, pow(b_eff_cur, 2.0));
	//Spline *Pk_current_model_scaled_mocks = add_const(Pk_current_model_scaled_beff, shotnoise_cur);
	BinInfo *R_bins_rescaled =  prep_bins(m_to_R(0.859981*Mmin_cur, cur_params_rescaled), m_to_R(0.859981*Mmax_cur, cur_params_rescaled), 300, LOG_BIN);
	Spline *variance_rescaled = prep_variances(R_bins_rescaled, Pk_current_model_zs_scaled);
	
	char scaled_model_out[200];
	sprintf(scaled_model_out, "%s/data/final/mocks_before/Pk/current_rescaled_model.dat", home_directory);
	print_spline_file(Pk_current_model_zs_scaled, k_binInfo, scaled_model_out);
	
	char target_model_out[200];
	sprintf(target_model_out, "%s/data/final/mocks_before/Pk/target_model.dat", home_directory);
	print_spline_file(Pk_target_model_mocks, k_binInfo, target_model_out);
	
	//exit(0);
		
	//ZA mock	
	char mockPath[200], Pk_out[200], rescaled_out[200];
	Catalogue *current_catalogue = NULL;
	
	double **frac_disps = NULL;
	calloc2D_double(&frac_disps, particle_no_cur, 3);
	
	char Pk_rescaled_measured_out[200], Pk_rescaled_measured_folded_out[200];
	char Pk_rescaled_fracDisps[200];
	
	
	for (int i = 0; i < 10; i++) {
	
		sprintf(mockPath, "%s/data/final/mocks_before/DM/original_mock_%d.dat", home_directory, i);
		current_catalogue = input_catalogue_file(mockPath, 0, mocks_format_withIds);
		
		sprintf(Pk_out, "%s/data/final/mocks_before/Pk/Pk_measured_%d.dat", home_directory, i);
		//haloes_measure_Pk(current_catalogue, Pk_out, 0.0, 1.0, CIC, 1.0);
		/*
		sprintf(Pk_realspace_out, "%s/data/rescaling/model_to_model/report_fig1/Pk_current_realspace_%d.dat", home_directory, i);
		haloes_measure_Pk(current_catalogue, Pk_realspace_out, volume/(double)particle_no, 1.0, CIC, 1.0);	
		
		toRedshift(current_catalogue, cur_params);
		
		sprintf(Pk_redshiftspace_out, "%s/data/rescaling/model_to_model/report_fig1/Pk_current_%d.dat", home_directory, i);
		haloes_measure_Pk(current_catalogue, Pk_redshiftspace_out, volume/(double)particle_no, 1.0, CIC, 1.0);				
		*/
		
		//apply box scaling
		scale_box(s);
		move_particles(current_catalogue, s);
		scale_velocities(s, current_catalogue, cur_params, targ_params);
		scale_masses(s, current_catalogue, cur_params, targ_params);
		
		sprintf(Pk_rescaled_measured_out, "%s/data/final/mocks_before/Pk/Pk_rescaled_measured_%d.dat", home_directory, i);
		sprintf(Pk_rescaled_measured_folded_out, "%s/data/final/mocks_before/Pk/Pk_rescaled_measured_4fold_%d.dat", home_directory, i);
		//haloes_measure_Pk(current_catalogue, Pk_rescaled_measured_out, 0.0, 1.0, CIC, z_scalefac);
		//haloes_measure_Pk(current_catalogue, Pk_rescaled_measured_folded_out, 0.0, 4.0, CIC, z_scalefac);
		//Spline *Pk_rescaled_measured = input_spline_file(Pk_rescaled_measured_out, k_Pkmono_format, MEASURED);
		//print_spline(Pk_rescaled_measured, Pk_testBins);
		
		//exit(0);
		
		//double *mass_lims_scaled = mass_min_max(current_catalogue);
		//double b_eff_rescaled = b_eff(Pk_current_model, cur_params_rescaled, mass_lims_scaled);		
		//printf("b eff of current catalogue: %lf \n", b_eff_cur);	
		
		//fractional disps and vels
		ZA_displacements(current_catalogue, &frac_disps, Pk_current_model_zs_scaled, Pk_target_model_mocks, true, R_nl);		
		apply_ZA_velocities(current_catalogue, &frac_disps, cur_params_rescaled);
		apply_ZA_displacements(current_catalogue, &frac_disps);// cur_params_rescaled, variance_rescaled);
		sprintf(Pk_rescaled_fracDisps, "%s/data/final/mocks_before/Pk/Pk_rescaled_measured_fracDisps_%d.dat", home_directory, i);
		haloes_measure_Pk(current_catalogue, Pk_rescaled_fracDisps, 0.0, 1.0, CIC, z_scalefac);
			
	
		//restore limits for next catalogue
		volume_limits[0] = 256.0;
		volume_limits[1] = 256.0;
		volume_limits[2] = 256.0;
		volume = volume_limits[0]*volume_limits[1]*volume_limits[2];
		
		free_catalogue(current_catalogue);

	}



	return 0;
}


int rescale_randoms_model_to_model_redshiftspace() {	

	char Pk_current_path[200];
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
	cur_params->z = 1.5; //actually 0.75 to start off with, check filename		
	
	/*
	char mocks_in[200];
	char cov_out[200];	
	sprintf(mocks_in, "/home/jonas/Testing_GR/Mocks/W1/mock_%s_zlim_0.6_0.9_Jf_1.dat", "%d");
	sprintf(cov_out, "%s/data/cov/mocks_before/305_VIPERS_mocks_cov.dat", home_directory);	
	generate_covariance_oneFile(mocks_in, cov_out, 305, 1, "%*e \t %le \t %le \t %*d");
	exit(0);	
	sprintf(mocks_in, "%s/data/ZA_mocks_report/before/HOD/ZA_mock_catalogue_HOD_%s.dat", home_directory, "%d");
	sprintf(cov_out, "%s/data/cov/mocks_before/300_mocks_cov.dat", home_directory);
	covariance_catalogue(mocks_in, 300, 0, cur_params, cov_out);	
	exit(0);
	*/
	
	int particle_no = 300000;
	double shotnoise_current = volume/particle_no;
	double Mmin_cur = 5e11;
	double Mmax_cur = 5e14;
	double R1 = m_to_R(Mmin_cur, cur_params);
	double R2 = m_to_R(Mmax_cur, cur_params);
	printf("\nM min of current: %le M_sun, M max of current: %le M_sun\n", Mmin_cur,  Mmax_cur);	
	printf("(current catalogue) R1: %lf Mpc H^-1, R2: %lf Mpc H^-1\n", R1, R2);	
	
	//--------------------------------------------------------
	
	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/models/matter_pk_om_m_0.15_z_1.0_sigma8_0.69.dat", home_directory);
	Spline *Pk_target_model = input_spline_file(Pk_target_path, Pk_model_format, MODEL);
		
	Parameters *targ_params = malloc(sizeof(*targ_params));
	targ_params->H0 = 100.0*0.678;
	targ_params->omega_m_0 = 0.15;
	targ_params->omega_v_0 = 0.85;
	targ_params->omega_r_0 = 0.0;
	targ_params->omega = 1.0;	
	targ_params->gamma = 0.545;
	targ_params->z = 1.0;		
	
	double R_nl = 0.0; //from linear Pk
	double Mmin_targ = 1e12;
	double Mmax_targ = 1e15;	
	double R1_primed = m_to_R(Mmin_targ, targ_params);
	double R2_primed = m_to_R(Mmax_targ, targ_params);
	double *R_lims = malloc(2 * sizeof(*R_lims));
	R_lims[0] = R1_primed;
	R_lims[1] = R2_primed;
	printf("\nM min of target: %le M_sun, M max of target: %le M_sun\n", Mmin_targ,  Mmax_targ);
	printf("(target catalogue) R1': %lf Mpc H^-1, R2': %lf Mpc H^-1\n", R1_primed, R2_primed);	
	
	//---------------------------------------------------------------------
	
	//find optimal s, z		
	bool vary_z_current = true;	
	int k_bins = 300;
	int z_bins = 300;
	Dplus_bins = 300;
	int R_bins = 300;
	double k_min = 0.05; //k_fund
	double k_max = 0.5; //k_nyq/4
	double z_min = 0.0;
	double z_max = 4.0;	
	int current_gravity = GR;
	int target_gravity = GR;		

	BinInfo *k_binInfo = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	BinInfo *z_binInfo = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	BinInfo *R_binInfo = prep_bins(0.4*R1_primed, 2.5*R2_primed, R_bins, LOG_BIN);		
	
	//D+(z), GR only
	Spline *Dplus_z_cur;
	Spline *Dplus_z_targ;
	//MAKE PREP_PK_CONSTZ USE ARRAY OF SPLINES, in case of modified grav varying k
	if (current_gravity == GR && target_gravity == GR) {
		Dplus_z_cur = Dplus_spline(GR, 0.0, z_binInfo, cur_params);
		Dplus_z_targ = Dplus_spline(GR, 0.0, z_binInfo, targ_params);	
	}
	

	//just changing current model from z = 0.75 to z = 3.5. //Prep Pk function keeps spline model
	BinInfo *k_binInfo_curModel = prep_bins(Pk_current_model->splineInfo->xmin, Pk_current_model->splineInfo->xmax, Pk_current_model->splineInfo->lines, LOG_BIN);
	//Spline *Pk_current_model_z0 = prep_Pk_constz(GR, Dplus_z_cur, 0.75, 0.0, Pk_current_model, k_binInfo_curModel, z_binInfo, cur_params);
	Pk_current_model = prep_Pk_constz(GR, Dplus_z_cur, 0.75, cur_params->z, Pk_current_model, k_binInfo_curModel, z_binInfo, cur_params);	
	free(k_binInfo_curModel);
	
	//Spline *variance_HOD = prep_variances(R_binInfo, Pk_current_model_z0);	
	
	
	double mass_lims_mocks[] = {Mmin_cur, Mmax_cur};
	char mocks_out[200];
	sprintf(mocks_out, "%s/data/final/mocks_before/DM/original_mock_%s.dat", home_directory, "%d");
	ZA_mocks(mocks_out, 10, particle_no, Pk_current_model, NULL, cur_params, mass_lims_mocks);
	exit(0);	
	
	BinInfo *R_binInfo_beff = prep_bins(R1, R2, 300, LOG_BIN);
	Spline *variance_current_beff = prep_variances(R_binInfo_beff, Pk_current_model);	
	double b_eff = calc_b_eff(cur_params, variance_current_beff, Mmin_cur, Mmax_cur);		
	printf("b eff in rescaling: %lf \n", b_eff);	
	
	
	
	// variance target with k-limits for minimisation
	Spline *Pk_target_model_klims = prep_Pk_constz(target_gravity, Dplus_z_cur, targ_params->z, targ_params->z, Pk_target_model, k_binInfo, z_binInfo, targ_params);
	Pk_target_model_klims->model = MEASURED;
	Spline *variance_target = prep_variances(R_binInfo, Pk_target_model_klims);
	

	
	//for current simulation, varying z:		
	Spline **Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	Spline **variance_splines_zBins = malloc((size_t)z_bins * sizeof(*variance_splines_zBins));	

	//finding a redshift to start minimising from
	double z_to = 0.0;
	double z_guess = 0.0;
	double sq_diff = 0.0;
	double k_chk = 0.0;
	double min_sq_diff = DBL_MAX;	
	BinInfo *zmin_k_bins = prep_bins(k_min, k_max, 20, LOG_BIN);
	for (int i = 0; i < z_bins; i++) {	
		z_to = bin_to_x(z_binInfo, i);
		printf("preparing Pk at z: %lf, z max: %lf \n", z_to, z_max);
	
		//current Pk at different redshifts at the right k limits
		Pk_splines_zBins[i] = prep_Pk_constz(current_gravity, Dplus_z_cur, cur_params->z, z_to, Pk_current_model, k_binInfo, z_binInfo, cur_params);
		Pk_splines_zBins[i]->model = MEASURED;	
		variance_splines_zBins[i] = prep_variances(R_binInfo, Pk_splines_zBins[i]);
	
		sq_diff = 0.0;
		for (int j = 0; j < zmin_k_bins->bins; j++) {
			k_chk = bin_to_x(zmin_k_bins, j);
			sq_diff += pow(splint_Pk(Pk_splines_zBins[i], k_chk) - splint_Pk(Pk_target_model, k_chk), 2.0);
			if (sq_diff < min_sq_diff) {
				min_sq_diff = sq_diff;
				z_guess = z_to;
			}	
		}		
	}	
	printf("z guess: %lf \n", z_guess);	

	//minimization
	dsq_multimin(vary_z_current, z_guess, variance_target, variance_splines_zBins, z_binInfo, R_lims);	
	exit(0);
	
	
	//did once, found those, its enough
	s = 1.086164;
	z_rescaled = 0.623107; //dsq at s, z: 7.544857e-06	
	
	Parameters *cur_params_rescaled = malloc(sizeof(*cur_params));		
	cur_params->H0 = 100.0*0.678;
	cur_params->omega_m_0 = 0.15;
	cur_params->omega_v_0 = 0.85;
	cur_params->omega_r_0 = 0.0;
	cur_params->omega = 1.0;	
	cur_params->gamma = 0.545;
	cur_params->z = z_rescaled;	

	
	k_binInfo_curModel = prep_bins(Pk_current_model->splineInfo->xmin, Pk_current_model->splineInfo->xmax, Pk_current_model->splineInfo->lines, LOG_BIN);
	Spline *Pk_current_model_zs_scaled = prep_Pk_constz_rescaled(current_gravity, Dplus_z_cur, cur_params->z, z_rescaled, Pk_current_model, k_binInfo_curModel, z_binInfo, cur_params, s);
	free(k_binInfo_curModel);
	


	//prefactor to change redshift
	double Dplus_cur = splint_generic(Dplus_z_cur, cur_params->z);
	double Dplus_after = splint_generic(Dplus_z_cur, z_rescaled);	
	
	//ZA mock	
	char mockPath[200], rescaled_out[200];
	Catalogue *current_catalogue = NULL;
	
	double **frac_disps = NULL;
	calloc2D_double(&frac_disps, particle_no, 3);
	
	char Pk_rescaled_measured_out[200], Pk_realspace_out[100], Pk_redshiftspace_out[100];
	Catalogue *satellites = NULL;
	Catalogue *centrals = NULL;
	
	for (int i = 0; i < 1; i++) {
	
		sprintf(mockPath, "%s/data/ZA_mocks_report/before/ZA_mock_catalogue_%d.dat", home_directory, i);
		current_catalogue = input_catalogue_file(mockPath, 0, mocks_format);
		/*
		sprintf(Pk_realspace_out, "%s/data/rescaling/model_to_model/report_fig1/Pk_current_realspace_%d.dat", home_directory, i);
		haloes_measure_Pk(current_catalogue, Pk_realspace_out, volume/(double)particle_no, 1.0, CIC, 1.0);	
		
		toRedshift(current_catalogue, cur_params);
		
		sprintf(Pk_redshiftspace_out, "%s/data/rescaling/model_to_model/report_fig1/Pk_current_%d.dat", home_directory, i);
		haloes_measure_Pk(current_catalogue, Pk_redshiftspace_out, volume/(double)particle_no, 1.0, CIC, 1.0);				
		*/
		
		//apply box scaling
		scale_box(s);
		move_particles(current_catalogue, s);
		scale_velocities(s, current_catalogue, cur_params, targ_params);
		scale_masses(s, current_catalogue, cur_params, targ_params);
		
		//fractional disps and vels
		ZA_displacements(current_catalogue, &frac_disps, Pk_current_model_zs_scaled, Pk_target_model, true, R_nl);		
		apply_ZA_velocities(current_catalogue, &frac_disps, targ_params);
		apply_ZA_displacements(current_catalogue, &frac_disps);
		
		//AT REDSHIFT 0 FOR HOD
		sprintf(Pk_rescaled_measured_out, "%s/data/tmp/Pk_rescaled.dat", home_directory);	
		haloes_measure_Pk(current_catalogue, Pk_rescaled_measured_out, 0.0, 1.0, CIC, 1.0/Dplus_cur);
		Spline *Pk_current_rescaled_measured_z0 = input_spline_file(Pk_rescaled_measured_out, k_Pkmono_format, MEASURED);
		
		Spline *variance_HOD_resc;
		
		Catalogue *centrals = HOD_centrals(current_catalogue);		
		Catalogue *satellites = HOD_satellites(current_catalogue, variance_HOD_resc, cur_params_rescaled);
		
		sprintf(rescaled_out, "%s/data/ZA_mocks_report/after/ZA_mock_catalogue_rescaled_%d.dat", home_directory, i);
		catalogue_to_file(current_catalogue, rescaled_out);
		
		//toRedshift(current_catalogue, targ_params);
		
		//char Pk_current_zs_scaled_out_RS_fracDisps[100];
		//sprintf(Pk_current_zs_scaled_out_RS_fracDisps, "%s/data/rescaling/model_to_model/fractional_runs/Pk_current_zs_scaled_redshift_fracDisps_%d.dat", home_directory, i);
		//haloes_measure_Pk(current_catalogue, Pk_current_zs_scaled_out_RS_fracDisps, 0.0, 1.0, CIC, Dplus_after/Dplus_cur);
		
	
		//Spline *Pk_current_rescaled = input_spline_file(Pk_current_zs_scaled_out, k_Pkmono_format, NORMAL, MEASURED);		
	
		//restore limits for next catalogue
		volume_limits[0] = 256.0;
		volume_limits[1] = 256.0;
		volume_limits[2] = 256.0;
		volume = volume_limits[0]*volume_limits[1]*volume_limits[2];
		
		free_catalogue(current_catalogue);

	}

	free2D_double(&frac_disps, particle_no);
	exit(0);
	

	
	//exit(0);	
		
	
	//RESCALED PARAMS HERE
	//scale velocities. Current z is z after rescaling? Parameters are now those of the target cosmology ? (H specifically)
	
	//scale_velocities(s, current_catalogue, cur_params, targ_params, z_current, z_rescaled);
	//toRedshift(current_catalogue, targ_params);
	
	//double** individual_displacements;
	//calloc2D_double(&individual_displacements, particle_no, 3);
	
	//dont forget to change to fractional displacements for ZA
	//ZA_displacements(current_catalogue, &individual_displacements, Pk_current_rescaled, Pk_target, true);
	//mass_bias_displacements(current_catalogue, &individual_displacements);
	//apply_ZA_displacements(current_catalogue, &individual_displacements, REAL_SPACE);
	//ZA_velocities(current_catalogue, &individual_displacements, targ_params);
	
	
	free(Pk_splines_zBins);
	free(variance_splines_zBins);
	free(k_binInfo);
	free(R_binInfo);
	free(z_binInfo);
	//gsl_integration_workspace_free(OLV_workspace);	
	return 0;		
}

/*
//REAL SPACE ONLY,WORKS
int rescale_randoms_model_to_model_realspace() {	

	char Pk_current_path[100];
	sprintf(Pk_current_path, "%s/data/models/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	Spline *Pk_current_model_wrongRedshift = input_spline_file(Pk_current_path, Pk_model_format, NORMAL, MODEL);	
	//this will go from redshift 0.75 to redshift 3.5 later	

	Parameters *cur_params = malloc(sizeof(*cur_params));		
	cur_params->H0 = 100.0*0.678;
	cur_params->omega_m_0 = 0.24;
	cur_params->omega_v_0 = 0.76;
	cur_params->omega_r_0 = 0.0;
	cur_params->omega = 1.0;	
	cur_params->gamma = 0.545;
	cur_params->z = 3.5;	
	
	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/models/matter_pk_om_m_0.15_z_1.0_sigma8_0.69.dat", home_directory);
	Spline *Pk_target_model = input_spline_file(Pk_target_path, Pk_model_format, NORMAL, MODEL);	
	
	Parameters *targ_params = malloc(sizeof(*targ_params));
	targ_params->H0 = 100.0*0.678;
	targ_params->omega_m_0 = 0.15;
	targ_params->omega_v_0 = 0.85;
	targ_params->omega_r_0 = 0.0;
	targ_params->omega = 1.0;	
	targ_params->gamma = 0.545;
	targ_params->z = 1.0;
	
	Mmin = 2e12;
	Mmax = 2e15;
	printf("mmin: %le M_sun,  mmax: %le M_sun\n", Mmin,  Mmax);
	R1_primed = m_to_R(Mmin, targ_params);
	R2_primed = m_to_R(Mmax, targ_params);
	printf("R1': %lf Mpc H^-1, R2': %lf Mpc H^-1\n", R1_primed, R2_primed);	
	
	//find optimal s, z		
	bool vary_z_current = true;	
	int k_bins = 300;
	int z_bins = 300;
	Dplus_bins = 300;
	int R_bins = 300;
	double k_min = 0.06; //kfund
	double k_max = 0.4; //knyq/4
	double z_min = 0.0;
	double z_max = 4.0;	
	int current_gravity = GR;
	int target_gravity = GR;	
	
	BinInfo *k_binInfo = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	BinInfo *z_binInfo = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	BinInfo *R_binInfo = prep_bins(0.5*R1_primed, 2.0*R2_primed, R_bins, LOG_BIN);
	
	Spline *Dplus_z_cur;
	Spline *Dplus_z_targ;
	//MAKE PREP_PK_CONSTZ USE ARRAY OF SPLINES, in case of modified grav varying k
	if (current_gravity == GR && target_gravity == GR) {
		Dplus_z_cur = Dplus_spline(GR, 0.0, z_binInfo, cur_params);
		Dplus_z_targ = Dplus_spline(GR, 0.0, z_binInfo, targ_params);	
	}
	
	//just changing current model from z = 0.75 to z = 3.5. //Prep Pk function makes spline nonmodel
	Spline *Pk_current_model = prep_Pk_constz(current_gravity, Dplus_z_cur, 0.75, cur_params->z, Pk_current_model_wrongRedshift, k_binInfo, z_binInfo, cur_params);
	//printf("Pk at k 0.4: %le \n", splint_Pk(Pk_current_model, 0.402705));
	//Pk_current_model = extend_spline_model(Pk_current_model->splineInfo);
	
	Spline *variance_current = prep_variances(R_binInfo, Pk_current_model);
	Spline *variance_target = prep_variances(R_binInfo, Pk_target_model);
	R_nl = 0.0;
	
	char Pk_current_model_out[100];
	sprintf(Pk_current_model_out, "%s/data/rescaling/model_to_model/realspace/Pk_current_model.dat", home_directory);
	haloModel_out(Pk_current_model_out, k_binInfo, 0.0, 1.0, REAL_SPACE, Pk_current_model, cur_params, 1.0);
	
	char Pk_target_model_out[100];
	sprintf(Pk_target_model_out, "%s/data/rescaling/model_to_model/realspace/Pk_target_model.dat", home_directory);
	haloModel_out(Pk_target_model_out, k_binInfo, 0.0, 1.0, REAL_SPACE, Pk_target_model, targ_params, 1.0);	

	//for current simulation, varying z:		
	Spline **Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	Spline **variance_splines_zBins = malloc((size_t)z_bins * sizeof(*variance_splines_zBins));	

	//finding a redshift to start minimising from
	double z_to, z_guess, sq_diff, min_sq_diff, k_chk;
	min_sq_diff = DBL_MAX;	
	BinInfo *zmin_k_bins = prep_bins(k_min, k_max, 50, LOG_BIN);
	for (int i = 0; i < z_bins; i++) {	
		z_to = bin_to_x(z_binInfo, i);
		printf("preparing Pk at z: %lf, z max: %lf \n", z_to, z_max);
	
		Pk_splines_zBins[i] = prep_Pk_constz(current_gravity, Dplus_z_cur, cur_params->z, z_to, Pk_current_model, k_binInfo, z_binInfo, cur_params);	
		variance_splines_zBins[i] = prep_variances(R_binInfo, Pk_splines_zBins[i]);
	
		sq_diff = 0.0;
		for (int j = 0; j < zmin_k_bins->bins; j++) {
			k_chk = bin_to_x(zmin_k_bins, j);
			sq_diff += pow(splint_Pk(Pk_splines_zBins[i], k_chk) - splint_Pk(Pk_target_model, k_chk), 2.0);
			if (sq_diff < min_sq_diff) {
				min_sq_diff = sq_diff;
				z_guess = z_to;
			}	
		}		
	}	
	printf("z guess: %lf \n", z_guess);
	//z_guess = 1.0;	


	//minimization
	dsq_multimin(vary_z_current, z_guess, variance_target, variance_splines_zBins, z_binInfo);	
	
	//s = 2.500123;
	//z_rescaled = 1.992451;
	
	int z_bin = x_to_bin(z_binInfo, z_rescaled);
	BinInfo *rescaled_k_binInfo = prep_bins(k_min*s, k_max*s, k_bins, LOG_BIN);	
		
	//print model after z scaling	
	Spline *Pk_current_model_z_scaled = prep_Pk_constz(current_gravity, Dplus_z_cur, cur_params->z, z_rescaled, Pk_current_model, k_binInfo, z_binInfo, cur_params);
	char Pk_current_model_z_scaled_out[100];
	sprintf(Pk_current_model_z_scaled_out, "%s/data/rescaling/model_to_model/realspace/Pk_current_model_z_scaled.dat", home_directory);
	haloModel_out(Pk_current_model_z_scaled_out, k_binInfo, 0.0, 1.0, REAL_SPACE, Pk_current_model_z_scaled, cur_params, 1.0);	
	
	//same z-scaled spline, but at k-values specific for s-scaling
	Spline *Pk_current_model_zs_scaled = prep_Pk_constz(current_gravity, Dplus_z_cur, cur_params->z, z_rescaled, Pk_current_model, rescaled_k_binInfo, z_binInfo, cur_params);	
	//zs scaling, model
	char Pk_current_model_zs_scaled_out[100];
	sprintf(Pk_current_model_zs_scaled_out, "%s/data/rescaling/model_to_model/realspace/Pk_current_model_zs_scaled.dat", home_directory);
	haloModel_out(Pk_current_model_zs_scaled_out, k_binInfo, 0.0, s*s*s, REAL_SPACE, Pk_current_model_zs_scaled, cur_params, s);
	//haloModel_out(Pk_current_model_zs_scaled_out, k_binInfo, 0.0, 1.0, REAL_SPACE, Pk_current_model_z_scaled, cur_params_rescaled, s);
	
	exit(0);

	
	//after rescaling, redshift changed, omega m and omega v of target now
	Parameters *cur_params_rescaled = malloc(sizeof(*cur_params));		
	cur_params_rescaled->H0 = 100.0*0.678;
	cur_params_rescaled->omega_m_0 = 0.15;
	cur_params_rescaled->omega_v_0 = 0.85;
	cur_params_rescaled->omega_r_0 = 0.0;
	cur_params_rescaled->omega = 1.0;	
	cur_params_rescaled->gamma = 0.545;
	cur_params_rescaled->z = z_rescaled;

	//prefactor to change redshift
	double Dplus_cur = splint_generic(Dplus_z_cur, cur_params->z);
	double Dplus_after = splint_generic(Dplus_z_cur, z_rescaled);	
	
	//ZA mock	
	char mockPath[100];
	Catalogue *current_catalogue;
	
	for (int i = 0; i < 20; i++) {
	
		sprintf(mockPath, "%s/data/ZA_mocks/ZA_mock_catalogue_%d.dat", home_directory, i);
		current_catalogue = input_catalogue_file(mockPath, 0, mocks_format);	

		//measure Pk
		char Pk_current_out[100];	
		sprintf(Pk_current_out, "%s/data/rescaling/model_to_model/realspace_runs/Pk_current_%d.dat", home_directory, i);
		haloes_measure_Pk(current_catalogue, Pk_current_out, 1.0, CIC, 1.0);

		//Proceeding with the measured Pk	
		//Spline *Pk_current = input_spline_file(Pk_current_out, k_Pkmono_format, NORMAL, MEASURED);	


		//measured after z scaling (just multiplied by the right prefactor)
		char Pk_current_z_scaled_out[100];
		sprintf(Pk_current_z_scaled_out, "%s/data/rescaling/model_to_model/realspace_runs/Pk_current_z_scaled_%d.dat", home_directory, i);
		haloes_measure_Pk(current_catalogue, Pk_current_z_scaled_out, 1.0, CIC, Dplus_after/Dplus_cur);	

		//apply box scaling
		scale_box(s);
		move_particles(current_catalogue, s);			

		//zs scaling, measured
		char Pk_current_zs_scaled_out[100];
		sprintf(Pk_current_zs_scaled_out, "%s/data/rescaling/model_to_model/realspace_runs/Pk_current_zs_scaled_%d.dat", home_directory, i);
		haloes_measure_Pk(current_catalogue, Pk_current_zs_scaled_out, 1.0, CIC, Dplus_after/Dplus_cur);
	
		//Spline *Pk_current_rescaled = input_spline_file(Pk_current_zs_scaled_out, k_Pkmono_format, NORMAL, MEASURED);		
	
	
		volume_limits[0] /= s;
		volume_limits[1] /= s;
		volume_limits[2] /= s;
		volume = volume_limits[0]*volume_limits[1]*volume_limits[2];
		
		free(current_catalogue->particles);
		free(current_catalogue);

	}


	
	

	
	exit(0);	
		
	
	//RESCALED PARAMS HERE
	//scale velocities. Current z is z after rescaling? Parameters are now those of the target cosmology ? (H specifically)
	
	//scale_velocities(s, current_catalogue, cur_params, targ_params, z_current, z_rescaled);
	//toRedshift(current_catalogue, targ_params);
	
	//double** individual_displacements;
	//calloc2D_double(&individual_displacements, particle_no, 3);
	
	//dont forget to change to fractional displacements for ZA
	//ZA_displacements(current_catalogue, &individual_displacements, Pk_current_rescaled, Pk_target, true);
	//mass_bias_displacements(current_catalogue, &individual_displacements);
	//apply_ZA_displacements(current_catalogue, &individual_displacements, REAL_SPACE);
	//ZA_velocities(current_catalogue, &individual_displacements, targ_params);
	
	
	free(Pk_splines_zBins);
	free(variance_splines_zBins);
	free(k_binInfo);
	free(R_binInfo);
	free(z_binInfo);
	//gsl_integration_workspace_free(OLV_workspace);	
	return 0;		
}
*/


//GR -> F5
//scale box -> ZA -> HOD?
/*
int rescale_testRun() {	

	double s_test = 0.9;
	double z_test = 0.3;
	double z_tmp;
	
	char Pk_current_path[100];
	sprintf(Pk_current_path, "%s/data/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	
	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/linear_matter_pk_0.76.dat", home_directory);	
		
	//Pk spline for target, model
	Pk_current = input_spline_file(Pk_current_path, Pk_model_format, NORMAL);
	//Pk_target = input_spline_file(Pk_target_path, Pk_model_format, false);	
	Pk_target = input_spline_file(Pk_current_path, Pk_model_format, NORMAL);		
	
	z_current = 0.75;
	z_target = 0.75;	
	
	//find optimal s, z		
	R1_primed = 0.05;
	R2_primed = 5.0;
	vary_z_current = true;
	//global variables
	k_bins = 300;
	z_bins = 200;
	Dplus_bins = 300;
	R_bins = 300;
	k_min = 0.01;
	k_max = 1.0;
	z_min = 0.0;
	z_max = 3.0;
	current_gravity = GR;
	target_gravity = GR;
	
	//binning info
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);	
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	rescaling_R_bin_info = prep_bins(0.6*R1_primed, 1.5*R2_primed, R_bins, LOG_BIN); 
	
	//integral workspaces
	z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	OLV_workspace = gsl_integration_workspace_alloc(OLV_workspace_size);
	
	//redshift-time relation reversed:	
	prep_redshift_spline(rescaling_z_bin_info);
	
	//for target simulation:
	Pk_target = prep_Pk_constz(target_gravity, z_target, z_test, Pk_target, rescaling_k_bin_info);
	Pk_target_extended = extend_spline_model(Pk_target);	
	variance_spline_target = prep_variances_test(rescaling_R_bin_info, Pk_target_extended, s_test);
	
	//for current simulation:		
	Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	Pk_splines_zBins_extended = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins_extended));
	variance_splines_zBins = malloc((size_t)z_bins * sizeof(*variance_splines_zBins));
	
	double min_sq_diff = DBL_MAX;
	double sq_diff;

	for (int i = 0; i < z_bins; i++) {
		z_tmp = bin_to_x(rescaling_z_bin_info, i);
		printf("Preparing OLV at different redshifts. zbin %d, z: %lf \n", i, z_tmp);
		
		Pk_splines_zBins[i] = prep_Pk_constz(current_gravity, z_current, z_tmp, Pk_current, rescaling_k_bin_info);
		
		sq_diff = pow(splint_generic(Pk_splines_zBins[i], 0.1) - splint_generic(Pk_target, 0.1), 2.0);
		if (sq_diff < min_sq_diff) {
			min_sq_diff = sq_diff;
			z_guess = z_tmp;
		}
		
		Pk_splines_zBins_extended[i] = extend_spline_model(Pk_splines_zBins[i]);

		variance_splines_zBins[i] = prep_variances(rescaling_R_bin_info, Pk_splines_zBins_extended[i]);
	}
	
	printf("z guess: %lf \n", z_guess);
	
		
	
	dsq_multimin(vary_z_current, z_guess, variance_spline_target, variance_splines_zBins);

	free(Pk_splines_zBins);
	free(Pk_splines_zBins_extended);
	free(variance_splines_zBins);
	gsl_integration_workspace_free(z_to_t_workspace);	
	gsl_integration_workspace_free(OLV_workspace);	
	return 0;		
}
*/

/*
int rescale_model() {	
	
	char Pk_current_path[100];
	sprintf(Pk_current_path, "%s/data/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	
	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/linear_matter_pk_0.76.dat", home_directory);	
	
	//Pk spline for target, model
	Pk_current = input_spline_file(Pk_current_path, Pk_model_format, NORMAL);
	Pk_target = input_spline_file(Pk_target_path, Pk_model_format, NORMAL);	
	
	z_current = 0.75;
	z_target = 0.76;	
	
	//find optimal s, z		
	//R1_primed = m_to_R(1e12, z_target);
	//R2_primed = m_to_R(5e15, z_target);
	
	printf("R1': %le, R2': %le \n", R1_primed, R2_primed);
	
	//varying z in current sim, box scaling always for current sim
	vary_z_current = true;

	k_bins = 300;
	z_bins = 200;
	Dplus_bins = 300;
	R_bins = 300;
	k_min = 0.01;
	k_max = 1.0;
	z_min = 0.0;
	z_max = 5.0;
	current_gravity = GR;
	target_gravity = GR;	
	
	//binning info
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);	
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	rescaling_R_bin_info = prep_bins(0.6*R1_primed, 1.5*R2_primed, R_bins, LOG_BIN); 
	
	//integral workspaces
	z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	OLV_workspace = gsl_integration_workspace_alloc(OLV_workspace_size);
	
	//redshift-time relation reversed:	
	prep_redshift_spline(rescaling_z_bin_info);
	
	//for target simulation:	
	Pk_target_extended = extend_spline_model(Pk_target);	
	variance_spline_target = prep_variances(rescaling_R_bin_info, Pk_target_extended);
	
	//for current simulation:		
	Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	Pk_splines_zBins_extended = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins_extended));
	variance_splines_zBins = malloc((size_t)z_bins * sizeof(*variance_splines_zBins));
	
	
	
	double min_sq_diff = DBL_MAX;
	double sq_diff;

	double z_cur;
	for (int i = 0; i < z_bins; i++) {
		z_cur = bin_to_x(rescaling_z_bin_info, i);
		printf("zbin %d, z: %lf \n", i, z_cur);

		Pk_splines_zBins[i] = prep_Pk_constz(current_gravity, z_current, z_cur, Pk_current, rescaling_k_bin_info);	
		
		sq_diff = pow(splint_generic(Pk_splines_zBins[i], 0.1) - splint_generic(Pk_target, 0.1), 2.0);
		if (sq_diff < min_sq_diff) {
			min_sq_diff = sq_diff;
			z_guess = z_cur;
		}

		//Pk_splines_zBins_extended[i] = extend_spline_model(Pk_splines_zBins[i]);
		
		//variance_splines_zBins[i] = prep_variances(rescaling_R_bin_info, Pk_splines_zBins_extended[i]);
	}
	
	printf("z guess: %lf \n", z_guess);
	//R1_primed = m_to_R(1e12, z_target);
	//R2_primed = m_to_R(5e15, z_target);
		
	
	for (int i = 0; i < z_bins; i++) {
		z_cur = bin_to_x(rescaling_z_bin_info, i);
		printf("zbin %d, z: %lf \n", i, z_cur);		

		Pk_splines_zBins_extended[i] = extend_spline_model(Pk_splines_zBins[i]);
		
		variance_splines_zBins[i] = prep_variances(rescaling_R_bin_info, Pk_splines_zBins_extended[i]);
	}
	

	//minimization to find s, z_rescaled (global vars)
	dsq_multimin(vary_z_current, z_guess, variance_spline_target, variance_splines_zBins);
	
	double dsqmin = dsq_minimized_value;
	for (int i = 0; i < 10; i++) {
		printf("\n RETRYING MINIMIZATION \n");
		dsq_multimin(vary_z_current, variance_spline_target, variance_splines_zBins);
		if (dsq_minimized_value < dsqmin) {
			dsqmin = dsq_minimized_value;
		}
	}
	

	int z_ind_final = x_to_bin(rescaling_z_bin_info, z_rescaled);
	
	Pk_current_extended = extend_spline_model(Pk_current);
	
	//current 
	char* current_outPath = prep_path("/data/rescaling/delsq_current.dat");
	//print_delsq(Pk_current_extended, current_outPath, 0.01, 0.1, 200, false);
	free(current_outPath);
	
	//target
	char* target_outPath = prep_path("/data/rescaling/delsq_target.dat");
	//print_delsq(Pk_target_extended, target_outPath, 0.01, 0.1, 200, false);
	free(target_outPath);
	
	//current with redshift scale
	char* current_outPath_z = prep_path("/data/rescaling/delsq_current_z.dat");
	//print_delsq(Pk_splines_zBins_extended[z_ind_final], current_outPath_z, 0.01, 0.1, 200, false);
	free(current_outPath_z);	
	
	//current after redshift AND box scale
	char current_outPath_zs[100];
	sprintf(current_outPath_zs, "%s/data/rescaling/delsq_current_zs.dat", home_directory);
	//print_delsq(Pk_splines_zBins_extended[z_ind_final], current_outPath_zs, 0.01, 0.1, 200, true); 

	
	generate_sigma_plots();	
	
	free(Pk_splines_zBins);
	free(Pk_splines_zBins_extended);
	free(variance_splines_zBins);
	free(rescaling_k_bin_info);
	free(rescaling_R_bin_info);
	free(rescaling_z_bin_info);
	gsl_integration_workspace_free(z_to_t_workspace);	
	gsl_integration_workspace_free(OLV_workspace);	
	return 0;		
}
*/




//CATALOGUE TO MODEL
/*
int rescale_catalogue_to_model() {	

	Parameters cur_params = {
		h = 0.678, //dimensionless hubble parameter
		H0 = 100.0*0.678, //[kms^-1 Mpc^-1]
		omega_m_0 = 0.24,
		omega_v_0 = 0.76,
		omega_r_0 = 0.0,
		omega = 1.0,	
		fg = pow(0.24, 0.545)	
	};	
	
	Parameters targ_params = {
		h = 0.678, //dimensionless hubble parameter
		H0 = 100.0*0.678, //[kms^-1 Mpc^-1]
		omega_m_0 = 0.24,
		omega_v_0 = 0.76,
		omega_r_0 = 0.0,
		omega = 1.0,	
		fg = pow(0.24, 0.4)	
	};	

	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);

	char sim_current[100];	
	sprintf(sim_current, "%s/%s", mock_directory, "AllHalo_GR_1500Mpc_a0.7_1.txt");	
	volume_limits[0] = 1500.0;  // h^-1 Mpc 
	volume_limits[1] = 1500.0;
	volume_limits[2] = 1500.0;
	volume = volume_limits[0]*volume_limits[1]*volume_limits[2];
	
	char outFile_scaled[100];
	sprintf(outFile_scaled, "%s/data/rescaling/cat_to_model/measured_pk_scaled.dat", home_directory);
	char outFile[100];
	sprintf(outFile, "%s/data/rescaling/cat_to_model/measured_pk.dat", home_directory);
	
	z_current = a_to_z(0.7);
	z_target = 0.75;
	
	//current simulation data input
	Catalogue *current_catalogue = input_catalogue_file(sim_current, 1, mocks_format);	
	
	//min and max masses, here for current catalogue, but actually must be for target catalogue
	mass_min_max(current_catalogue);
	printf("mmin: %le,  mmax: %le\n", Mmin,  Mmax);
	R1_primed = m_to_R(Mmin, z_current);
	R2_primed = m_to_R(Mmax, z_current);
	printf("R1': %le, R2': %le \n", R1_primed, R2_primed);		
	
	overdensity_allocateMemory();
	haloes_measure_Pk(current_catalogue, outFile, 1.0, CIC);
	overdensity_freeMemory();
	
	Pk_current = input_spline_file(outFile, k_Pkmono_format, NORMAL);
	
	Pk_target = input_spline_file(Pk_target_path, Pk_model_format, NORMAL);	
	
	printf("current kmin: %lf, kmax: %lf \n", Pk_current->xmin, Pk_current->xmax);
	
	//find optimal s, z		
	vary_z_current = true;
	
	k_bins = 300;
	z_bins = 200;
	Dplus_bins = 300;
	R_bins = 300;
	k_min = 0.02;
	k_max = 0.2;
	z_min = 0.0;
	z_max = 4.0;	
	current_gravity = GR;
	target_gravity = GR;	
	
	
	
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	rescaling_R_bin_info = prep_bins(0.6*R1_primed, 1.5*R2_primed, R_bins, LOG_BIN); 
	//HOD_mass_bin_info = prep_bins(Mmin, Mmax, 100, NORMAL_BIN);
		
	
	//integral workspaces
	z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	OLV_workspace = gsl_integration_workspace_alloc(OLV_workspace_size);
	
	//redshift-time relation reversed:	
	prep_redshift_spline(rescaling_z_bin_info);
	
	//for current simulation:		
	Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	Pk_splines_zBins_extended = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins_extended));
	variance_splines_zBins = malloc((size_t)z_bins * sizeof(*variance_splines_zBins));
	Pk_current_extended = extend_spline(Pk_current);	
	
	//for target simulation:	
	Pk_target_extended = extend_spline_model(Pk_target);	
	variance_spline_target = prep_variances(rescaling_R_bin_info, Pk_target_extended);

	//finding a redshift to start minimising from
	double z_cur, z_guess, sq_diff, min_sq_diff;
	min_sq_diff = DBL_MAX;
	for (int i = 0; i < z_bins; i++) {	
		z_cur = bin_to_x(rescaling_z_bin_info, i);
		Pk_splines_zBins[i] = prep_Pk_constz(current_gravity, z_current, z_cur, Pk_current, rescaling_k_bin_info);	
		
		//check more points
		sq_diff = pow(splint_generic(Pk_splines_zBins[i], 0.1) - splint_generic(Pk_target, 0.1), 2.0);
		if (sq_diff < min_sq_diff) {
			min_sq_diff = sq_diff;
			z_guess = z_cur;
		}			
	}	
	
	printf("z guess: %lf \n", z_guess);			
	
	//prepare OLV bins at different redshifts
	for (int i = 0; i < z_bins; i++) {
		z_cur = bin_to_x(rescaling_z_bin_info, i);
		printf("zbin %d, z: %lf \n", i, z_cur);		

		Pk_splines_zBins_extended[i] = extend_spline(Pk_splines_zBins[i]);		
		variance_splines_zBins[i] = prep_variances(rescaling_R_bin_info, Pk_splines_zBins_extended[i]);
	}
	
	//current 
	char current_outPath[100];
	sprintf(current_outPath, "%s/data/rescaling/cat_to_model/delsq_current.dat", home_directory);
	print_delsq(Pk_current_extended, current_outPath, rescaling_k_bin_info, false);
	
	//target
	char target_outPath[100];
	sprintf(target_outPath, "%s/data/rescaling/cat_to_model/delsq_target.dat", home_directory);
	print_delsq(Pk_target_extended, target_outPath, rescaling_k_bin_info, false);	
	
	//minimization
	dsq_multimin(vary_z_current, z_guess, variance_spline_target, variance_splines_zBins);
	int z_ind_final = x_to_bin(rescaling_z_bin_info, z_rescaled);
	BinInfo* rescaled_k_bins = prep_bins(k_min/s, k_max/s, k_bins, LOG_BIN);
	
	//current s scaled
	char current_outPath_s[100];
	sprintf(current_outPath_s, "%s/data/rescaling/cat_to_model/delsq_current_s.dat", home_directory);
	print_delsq(Pk_current_extended, current_outPath_s, rescaling_k_bin_info, true);
	
	//scale box
	volume_limits[0] *= s;
	volume_limits[1] *= s;
	volume_limits[2] *= s;
	volume = volume_limits[0]*volume_limits[1]*volume_limits[2];
	
	//ensure PBC
	for (int i = 0; i < current_catalogue->particle_no; i++) {
		PBC(&(current_catalogue->particles[i].x), &(current_catalogue->particles[i].y), &(current_catalogue->particles[i].z));
	}	
	
	//scale masses, OMEGA M DEPENDS ON Z!!!!!!!!!!!
	for (int i = 0; i < particle_no; i++) {		
		current_catalogue->particles[i].mass *= pow(s, 3.0) * (omega_m_primed/omega_m_z(z_rescaled));
	}
	
	mass_min_max(current_catalogue);
	printf("AFTER RESCALING. M min: %le, M max: %le\n", Mmin, Mmax);
	
	overdensity_allocateMemory();
	haloes_measure_Pk(current_catalogue, outFile_scaled, 1.0, CIC);			
	overdensity_freeMemory();	
	
	SplineInfo* Pk_current_rescaled = input_spline_file(outFile_scaled, k_Pkmono_format, NORMAL);
	//print_spline(Pk_current_rescaled, rescaling_k_bin_info);
	//exit(0);
	
	SplineInfo_Extended* Pk_current_rescaled_extended = extend_spline(Pk_current_rescaled);
	
	
	//current s scaled, measured
	char current_outPath_s_measured[100];
	sprintf(current_outPath_s_measured, "%s/data/rescaling/cat_to_model/delsq_current_s_measured.dat", home_directory);
	print_delsq(Pk_current_rescaled_extended, current_outPath_s_measured, rescaling_k_bin_info, false);
	
	
	
	
	exit(0);		
	
	
	//scale velocities. Current z is z after rescaling? Parameters are now those of the target cosmology ? (H specifically)
	scale_velocities(s, current_catalogue, &cur_params, &targ_params, z_current, z_rescaled);
	toRedshift(current_catalogue, &targ_params, z_rescaled);
	
	double** individual_displacements;
	calloc2D_double(&individual_displacements, particle_no, 3);
	//dont forget to change to fractional displacements for ZA
	//ZA_displacements(current_catalogue, &individual_displacements, Pk_target_extended, NGP);
	mass_bias_displacements(current_catalogue, &individual_displacements);
	apply_ZA_displacements(current_catalogue, &individual_displacements, REAL_SPACE);
	ZA_velocities(current_catalogue, &individual_displacements, &targ_params, z_rescaled);
	

	//calc_b_eff(z_rescaled);
	//bias_plot();
	
	

	
	//scale masses, velocities with s
	//does this come before or after HOD??
	//NEED SOME PARAMETERS HERE!!
	//scale_masses();
	//scale_velocities();
	
	
	//Pk_out(out_final);

	
	free(Pk_splines_zBins);
	free(Pk_splines_zBins_extended);
	free(variance_splines_zBins);
	free(rescaling_k_bin_info);
	free(rescaling_R_bin_info);
	free(rescaling_z_bin_info);
	gsl_integration_workspace_free(z_to_t_workspace);	
	gsl_integration_workspace_free(OLV_workspace);	
	return 0;		
}
*/



