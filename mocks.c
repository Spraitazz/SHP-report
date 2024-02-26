Catalogue* populateHaloes_oneMass(int halo_no, double halo_mass) {	
	
	Catalogue* toReturn = malloc(sizeof(*toReturn));	
	Particle **particles = malloc((size_t)halo_no * sizeof(*particles));
	printf("Populating catalogue with %d haloes of mass %le in a box of length %lf\n", halo_no, halo_mass, volume_limits[0]);
	double x, y, z;
	for (int i = 0; i < halo_no; i++) {	
		x = randomDouble(0.0, volume_limits[0]);	
		y = randomDouble(0.0, volume_limits[1]);
		z = randomDouble(0.0, volume_limits[2]);
		char *label = malloc(20 * sizeof(*label));
		sprintf(label, "halo %d", i);				
		particles[i] = new_particle(x, y, z, 0.0, 0.0, 0.0, halo_mass, label, HALO, i);		
	}	
	
	toReturn->particle_no = halo_no;	
	toReturn->particles = particles;			
	return toReturn;
}

int ZA_mocks(char outPath[], int mock_no, int particle_no, Pk_Spline *pk_spline, Spline *variance_HOD, Parameters *params, double mass_lims[]) {
	
	printf("\npreparing ZA MOCKS \n");	
	
	double Mmin = mass_lims[0];
	double Mmax = mass_lims[1];
	double R1 = m_to_R(Mmin, params);
	double R2 = m_to_R(Mmax, params);
	BinInfo *R_binInfo = prep_bins(R1, R2, 300, LOG_BIN);
	Spline *variance = prep_variances(R_binInfo, pk_spline);
	printf("linear variance min: %lf, linear variance max: %lf \n", splint_generic(variance, variance->xmin), splint_generic(variance, variance->xmax));	
	
	/*
	double tmp_m = 1.4e+14;
	printf("i want cdm: %lf, b(M): %lf\n", HOD_params(tmp_m, variance_HOD, params)->cdm, b_m(tmp_m, params, variance));
	exit(0);
	*/
	
	//prepare cumsums from fvdv/m integral
	double m = 0.0;
	double v = 0.0;
	int massBins = 200;
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
	//normalise last element to 1.0
	for (int i = 0; i < massBins; i++) cumsums[i] /= cumsums[massBins-1];	
	//reversed spline needed
	Spline *fv_cumsums_reverse = new_spline(massBins, cumsums, ms);
	free(ms);
	free(cumsums);
	
	//effective bias, should match the one given
	double b_eff = calc_b_eff(params, variance, Mmin, Mmax);	
	printf("b eff from integral (in ZA_mocks() ): %lf \n\n", b_eff);
	//displacements from linear spectrum
	double R_nl = 0.0;		

	//displacement storage, catalogue
	double **disps_init = NULL;
	calloc2D_double(&disps_init, particle_no, 3);
	Catalogue *current_catalogue = NULL;
	Catalogue *satellites = NULL;
	Catalogue *centrals = NULL;
	char cur_outPath[200];
	//char dm_outPath[200];	

	for (int i = 0; i < mock_no; i++) {	
	
		printf("\nMaking mock %d\n", i); 			

		current_catalogue = populateHaloes_oneMass(particle_no, 0.0);		
		
		//set ST masses
		catalogue_setMasses(current_catalogue, fv_cumsums_reverse);		
				
		//generate ZA disp field
		ZA_displacements(current_catalogue, &disps_init, pk_spline, NULL, false, R_nl);
		
		//large scale velocity field not biased tcdm pg. 3
		//apply_ZA_velocities(current_catalogue, &disps_init, params);
		
		//biased displacements
		apply_biased_ZA_displacements(current_catalogue, &disps_init, params, variance);	
		
		//centrals = HOD_centrals(current_catalogue);			
		//satellites = HOD_satellites(current_catalogue, variance_HOD, params);
		
		//sprintf(dm_outPath, "%s/data/ZA_mocks_report/before/ZA_mock_catalogue_%d.dat", home_directory, i);	
		//catalogue_to_file(current_catalogue, dm_outPath);		
		//free_catalogue(current_catalogue);	
			
		//combine makes hard copy
		//current_catalogue = combine_catalogues(centrals, satellites);		
		//free_catalogue(centrals);
		//free_catalogue(satellites);
		
		sprintf(cur_outPath, outPath, i);
		//printf("before randomising \n");
		//randomise_catalogue(current_catalogue);	
		catalogue_to_file(current_catalogue, cur_outPath);
		free_catalogue(current_catalogue);
	}
	
	free2D_double(&disps_init, particle_no);	
	return 0;
}


/*
int HOD_mocks_test2() {

	char Pk_current_path[100];
	sprintf(Pk_current_path, "%s/data/models/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	Spline *Pk_current_model = input_spline_file(Pk_current_path, Pk_model_format, MODEL);	
	//this will go from redshift 0.75 to redshift 3.5 later	
		
	Parameters *cur_params = malloc(sizeof(*cur_params));		
	cur_params->H0 = 67.8;
	cur_params->omega_m_0 = 0.24;
	cur_params->omega_v_0 = 0.76;
	cur_params->omega_r_0 = 0.0;
	cur_params->omega = 1.0;	
	cur_params->gamma = 0.545;
	cur_params->z = 1.5; //actually 0.75 to start off with, check filename	
	
	char mocks_in[200];
	char cov_out[200];	
	sprintf(mocks_in, "%s/data/ZA_mocks_report/before/HOD/ZA_mock_catalogue_HOD_%s.dat", home_directory, "%d");
	sprintf(cov_out, "%s/data/cov/mocks_before/300_mocks_cov_v3.dat", home_directory);	
	covariance_catalogue(mocks_in, 300, 0, cur_params, cov_out);
	//generate_covariance_oneFile(mocks_in, cov_out, 305, 1, "%*e \t %le \t %le \t %*d");
	exit(0);
	
	
	char mocks_in[200];
	char Pk_out[200];	
	sprintf(mocks_in, "%s/data/ZA_mocks_HOD_z=1.5/DM/ZA_mock_catalogue_DM_only_%s.dat", home_directory, "%d");
	sprintf(Pk_out, "%s/data/ZA_mocks_HOD_z=1.5/Pk_DM/ZA_mock_catalogue_DM_only_Pk_%s.dat", home_directory, "%d");
	catalogues_measure_Pk(mocks_in, Pk_out, mocks_format_withIds, 200, 0, 1.0, CIC, 1.0, cur_params);
	exit(0);
	
	int particle_no = 300000;
	double Mmin_cur = 5e11;
	double Mmax_cur = 5e14;
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
	BinInfo *k_binInfo_curModel = prep_bins(Pk_current_model->xmin, Pk_current_model->splineInfo->xmax, Pk_current_model->splineInfo->lines, LOG_BIN);
	Spline *Pk_current_model_z0 = prep_Pk_constz(GR, Dplus_z_cur, 0.75, 0.0, Pk_current_model, k_binInfo_curModel, z_binInfo, cur_params);
	Pk_current_model = prep_Pk_constz(GR, Dplus_z_cur, 0.75, cur_params->z, Pk_current_model, k_binInfo_curModel, z_binInfo, cur_params);	
	free(k_binInfo_curModel);
	
	Spline *variance_HOD = prep_variances(R_binInfo, Pk_current_model_z0);
	
	
	
	double mass_lims_mocks[] = {Mmin_cur, Mmax_cur};
	char mocks_out[100];
	sprintf(mocks_out, "%s/data/ZA_mocks_HOD_z=1.5/DM/ZA_mock_catalogue_DM_only_%s.dat", home_directory, "%d");	
	ZA_mocks(mocks_out, 200, particle_no, Pk_current_model, variance_HOD, cur_params, mass_lims_mocks);
	
		
	//weighted SN - volume/sat no
	sat_no_tot /= 10;
	weighted_shotnoise_tot /= 10.0;

	printf("average no of sats: %d, average weighted shotnoise: %lf \n", sat_no_tot, weighted_shotnoise_tot);
	
	
	
	
	exit(0);



	return 0;
}
*/


//one mass
/*
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
	double Mmax_cur = 5e14;
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
	
	
	HOD_Parameters *HOD_parameters = HOD_params(5e13, variance_HOD, cur_params);
	double cdm = HOD_parameters->cdm;	
	double rvir = HOD_parameters->rvir;
	double rs = HOD_parameters->rs;	
	printf("\n FROM RESCALING PARAMS: cdm: %lf, rvir: %lf, rs: %lf \n", cdm, rvir, rs);

	
	char Pk_sats_model[200];
	sprintf(Pk_sats_model, "%s/data/rescaling/model_to_model/HOD3_simple_vels/Pk_sats_model.dat", home_directory);
	haloModel_RSD_onlySats(Pk_sats_model, k_binInfo, 194.035569, 1.0, Pk_current_model, 1.0, 0.0, HOD_parameters, 639362, cur_params);
	
	//exit(0);
	
	double mass_lims_mocks[] = {Mmin_cur, Mmax_cur};
	char mocks_out[100];
	sprintf(mocks_out, "%s/data/ZA_mocks_HOD/ZA_mock_catalogue_%s.dat", home_directory, "%d");	
	ZA_mocks_test(mocks_out, 10, particle_no, Pk_current_model, variance_HOD, cur_params, mass_lims_mocks);
	//exit(0);	

	
		
	//weighted SN - volume/sat no
	sat_no_tot /= 10;
	weighted_shotnoise_tot /= 10.0;

	printf("average no of sats: %d, average weighted shotnoise: %lf \n", sat_no_tot, weighted_shotnoise_tot);
	
	
	
	
	exit(0);



	return 0;
}
*/


/*

Catalogue* populateHaloes_twoMasses(double smallMass, double bigMass, int halo_no, int no_big) {
	Catalogue* toReturn = malloc(sizeof(*toReturn));
	toReturn->particle_no = halo_no;
	Particle* particles = initParticles(halo_no);	
	for (int i = 0; i < halo_no; i++) {	
		//positions		
		particles[i].x = randomDouble(0.0, volume_limits[0]);	
		particles[i].y = randomDouble(0.0, volume_limits[1]);
		particles[i].z = randomDouble(0.0, volume_limits[2]);
		
		// velocities
		particles[i].vx = 0.0;	
		particles[i].vy = 0.0;
		particles[i].vz = 0.0;	
		
		// mass
		if (i < no_big) {
			particles[i].mass = bigMass;
		} else {
			particles[i].mass = smallMass;
		}
	}		
	toReturn->particles = particles;			
	return toReturn;	
}

*/

/*
Catalogue* populateParticles_crystalLattice(int particle_no) {
	if (particle_no != cells[0]*cells[1]*cells[2]) {
		printf("for lattice particle number should equal cell number\n");
		exit(0);
	}	
	double cellSizes[3];
	cellSizes[0] = volume_limits[0] / (double) cells[0];
	cellSizes[1] = volume_limits[1] / (double) cells[1];
	cellSizes[2] = volume_limits[2] / (double) cells[2];
	
	Catalogue* toReturn = malloc(sizeof(*toReturn));
	toReturn->particle_no = particle_no;
	Particle* particles = initParticles(particle_no);	
	
	//placing particle in center of cell
	for (int i = 0; i < cells[0]; i++) {
		for (int j = 0; j < cells[1]; j++) {
			for (int k = 0; k < cells[2]; k++) {

				particles[arr_ind(i,j,k)].x = ((double)i + 0.5) * cellSizes[0];
				particles[arr_ind(i,j,k)].y = ((double)j + 0.5) * cellSizes[1];
				particles[arr_ind(i,j,k)].z = ((double)k + 0.5) * cellSizes[2];
								
				// velocities
				particles[arr_ind(i,j,k)].vx = 0.0;	
				particles[arr_ind(i,j,k)].vy = 0.0;
				particles[arr_ind(i,j,k)].vz = 0.0;	
		
				// mass
				particles[arr_ind(i,j,k)].mass = 0.0;
			}
		}
	}
	toReturn->particles = particles;
	return toReturn;
}
*/

/*

int populateHaloes_diffMass_lattice(double smallMass, double bigMass, int no_big) {	

	if (particle_no != cells[0]*cells[1]*cells[2]) {
		printf("for lattice particle number should equal cell number\n");
		exit(0);
	}	
	
	double cellSizes[3];
	cellSizes[0] = volume_limits[0] / (double) cells[0];
	cellSizes[1] = volume_limits[1] / (double) cells[1];
	cellSizes[2] = volume_limits[2] / (double) cells[2];
	int p_no = 0;
	
	//placing particle in center of cell
	for (int i = 0; i < cells[0]; i++) {
		for (int j = 0; j < cells[1]; j++) {
			for (int k = 0; k < cells[2]; k++) {
				particles[arr_ind(i,j,k)][0] = ((double)i + 0.5) * cellSizes[0];
				particles[arr_ind(i,j,k)][1] = ((double)j + 0.5) * cellSizes[1];
				particles[arr_ind(i,j,k)][2] = ((double)k + 0.5) * cellSizes[2];
								
				// velocities
				particles[arr_ind(i,j,k)][3] = 0.0;	
				particles[arr_ind(i,j,k)][4] = 0.0;
				particles[arr_ind(i,j,k)][5] = 0.0;			
				
				// mass
				if (p_no < no_big) {
					particles[arr_ind(i,j,k)][6] = bigMass;
				} else {
					particles[arr_ind(i,j,k)][6] = smallMass;
				}
				
				p_no += 1;

			}
		}
	}
	return 0;
}

int populateHaloes_randomMass(double M_min, double M_max) {	
	for (int i = 0; i < halo_no; i++) {	
		//positions		
		particles[i][0] = randomDouble(0.0, volume_limits[0]);	
		particles[i][1] = randomDouble(0.0, volume_limits[1]);
		particles[i][2] = randomDouble(0.0, volume_limits[2]);
		
		// velocities
		particles[i][3] = 0.0;	
		particles[i][4] = 0.0;
		particles[i][5] = 0.0;		
		
		double mass = pow(10.0, 12.0 + gsl_ran_gaussian(gsl_ran_r, 0.5));			
		
		particles[i][6] = mass;
		
	}				
	return 0;
}



int populateHaloes_catalogueMass(char inFile[]) {
	FILE *catalogue_input = fopen(inFile, "r");
	int num_lines = countLines(catalogue_input);
	
	particle_no = num_lines-1;
	initParticles();
	
	ph = fscanf(catalogue_input, "%*[^\n]\n");
	printf("adding %d particles \n", particle_no);
	
	for (int i = 0; i < particle_no; i++) {
		//positions		
		particles[i][0] = randomDouble(0.0, volume_limits[0]);	
		particles[i][1] = randomDouble(0.0, volume_limits[1]);
		particles[i][2] = randomDouble(0.0, volume_limits[2]);
		
		// velocities
		particles[i][3] = 0.0;	
		particles[i][4] = 0.0;
		particles[i][5] = 0.0;	
		
		// mass
		ph = fscanf(catalogue_input, "%*f \t %*f \t %*f \t %*f \t %*f \t %*f \t %le", &particles[i][6]);	
	}	

	mass_min_max();
}
*/

