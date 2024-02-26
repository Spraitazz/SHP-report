int test_random(bool verbose) {
	double min = 0.;
	double max = 10.;
	double r;
	for (int i = 0; i < 10000; i++) {
		r = randomDouble(min, max);
		if (r < min || r > max) {
			printf("Random number bounds error\n");
			return 1;
		}	
	}
	if (verbose) printf("random number generator OK\n");
	return 0;
}

int test_file_io(bool verbose) {
	char filename[100];
	sprintf(filename, "%s%s", home_directory, "testFile");
	FILE *writer = fopen(filename, "w");
	if (writer == NULL) {
		printf("error creating file\n");
		return 1;
	}
	for (int i = 0; i < 1000; i++) {
		fprintf(writer, "%d\n", i);
	}
	fclose(writer);
	FILE *reader = fopen(filename, "r");
	if (reader == NULL) {
		printf("error reading file\n");
		return 1;
	}
	int lines = countLines(reader);
	if (lines != 1000) {
		printf("wrong line count error\n");
		return 1;
	}
	rewind(reader);
	int zero;
	ph = fscanf(reader, "%d", &zero);
	if (zero != 0) {
		printf("IO mismatch\n");
		return 1;
	}
	fclose(reader);
	remove(filename);
	reader = fopen(filename, "r");
	if (reader != NULL) {
		printf("file remove error\n");
		return 1;
	}
	if (verbose) printf("File IO test OK\n");
	return 0;
}

int randoms_populate_test(bool verbose) {
	int save_no = particle_no;
	particle_no = 10000;
	initParticles();
	populateParticles();
	if (particles[9999][0] == particles[9999][1] && particles[9999][1] == particles[9999][2]) {
		printf("check random numbers in particle population with randoms\n");
		return 1;
	}
	particles[0][3] = 100.;
	double temp = particles[0][2];
	toRedshift();
	if (particles[0][2] == temp) {
		printf("toRedshift() error\n");
		return 1;
	}
	clearParticles();
	if (particles[0][0] != 0.) {
		printf("memory clear error at particles array\n");
		return 1;
	}
	particle_no = save_no;
	if (verbose) printf("randoms populate, toRedshift() OK\n");
	return 0;
}

int overdensity_pk_test(bool verbose) {
	int save_no = particle_no;
	particle_no = 100000;
	//unfolded test
	initGrid();
	fftw_allocateMemory();
	//double Pk_avg, k_maxmodes;
	for (int i = 0; i < 20; i++) {
		noise_randoms(particle_no, "", false, true);
	
	}
	//pick up data, check for consistency
	//run ~20 times
	//find bin with largest number of modes, make sure monopole average is 1/nbar and quadrupole is around 0
	
	//folded test
	noise_randoms(particle_no, "", true, true);	
	
	particle_no = save_no;
	if (verbose) printf("unfolded, folded randoms power spectra OK\n");
	return 0;

}

int memory_test(bool verbose) {
	int save_no = particle_no;
	particle_no = 1;
	initParticles();	
	populateParticles();
	
	initGrid();
	populateGrid();
	grid[100] = 1.5;
		
	fftw_allocateMemory();
	gridIntoOverdensity();
	overdensity_fourier[100][0] = 1.;
	
	overdensity_pk(false);
	
	clearMemory();
	
	if (particles[0][0] != 0. || particles[0][1] != 0. || particles[0][2] != 0.) {
		printf("particle clear failed\n");
		return 1;
	}
	if (overdensity_fourier[100][0] != 0.) {
		printf("overdensity clear failed\n");
		return 1;
	}
	if (grid[100] != 0.) {
		printf("grid clear failed\n");
		return 1;
	}
	
	particle_no = save_no;	
	if (verbose) printf("memory cleaning OK\n");
	return 0;
}

double delta_sq_int_func_test(double R, void *params) {
	params = params;
	return 1.0/R;
}

int delta_sq_rms_test() {
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	double error, result, R1_primed, R2_primed;
	R1_primed = 0.1;
	R2_primed = 10.0;
			
	//function under the integral
	gsl_function F;	
	F.function = &delta_sq_int_func_test;
	F.params = 0;        
	
    gsl_integration_qags(&F, R1_primed, R2_primed, 0, 1e-7, 1000, w, &result, &error);	
	result /= log(R2_primed / R1_primed);
	
	gsl_integration_workspace_free(w);
	printf("delta sq test result: %lf\n", result);
	return 0;	
}


int full_test(bool verbose) {
	//normally suppress printing, only print final result and errors, except "VERBOSE" mode
	init_rng();
	test_random(verbose);
	test_file_io(verbose);
	randoms_populate_test(verbose);	
	overdensity_pk_test(verbose);
	memory_test(verbose);
	return 0;
}

//full test run with randomly generated particles! CODE THIS UP
int test_run_randoms() {

	volume_limits[0] = 1500.;  // h^-1 Mpc 
	volume_limits[1] = 1500.;
	volume_limits[2] = 1500.;
	volume = volume_limits[0] * volume_limits[1] * volume_limits[2];
	
	//define grid size, EVEN NUMBERS ONLY TO AVOID PROBLEMS, 2^n best for FFT	
	cells[0] = 256;
	cells[1] = 256;
	cells[2] = 256;

	particle_no = 10000;	

	initGrid();	
	return 0;
}

int z_to_t_test() {
	//prep_redshift_spline();
	double z, Mpc, seconds_in_year, t, human_t;
	seconds_in_year = 3.154 * 1e7; //thx wolfram
	Mpc = 3.086e+22;
	for (int i = 0; i <= 20; i++) {
		z = 0.1*(double)i;
		t = z_to_t(z);
		human_t = t * (Mpc / 1000.0);
		human_t /= seconds_in_year;
		printf("z: %lf, t: %le, human t, years: %le\n", z, z_to_t(z), human_t);
	}
	return 0;
}

int Dplus_EdS_test() {
	omega_m = 1.0;
	omega_v = 0.0;
	z_min = 0.0;
	z_max = 10.0;
	Dplus_bins = 200;
	z_to_t_workspace = gsl_integration_workspace_alloc(1000);
	prep_redshift_spline(z_min, z_max, 1000);
	double* zs;
	double* Dpluses;
	double* zs = malloc((size_t)Dplus_bins * sizeof(*zs));
	double* Dpluses = malloc((size_t)Dplus_bins * sizeof(*Dpluses));
	Dplus_calc(GR, 0.0, "", &zs, &Dpluses);
}

int diffeq_test_function(double t, const double x[], double dxdt[], void *params) {
	params = params;
	double g = 9.81;
	double L = 1.0;	
	dxdt[0] = x[1];
	dxdt[1] = (-g/L)*sin(x[0]);
	return GSL_SUCCESS;
}

int diffeq_test() {
	double t_current = 0.0;
	double t_next;
	double dt = 0.01;	
	double x[2] = {1.0, 0.0}; 
	int status;
	
	gsl_odeiv2_system sys = {diffeq_test_function, NULL, 2, NULL}; //rk8pd
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, dt, 1e-7, 0.0);	
	printf("t: %lf, x: %le \n", t_current, x[0]);
		
	for (int i = 1; i < 200; i++) {	
		//the next iteration goes from t_curr to t_next in some number of steps chosen by the driver
		t_next = t_current + dt;
		printf("i: %d, t current: %lf, t next: %lf, x: %le \n", i, t_current, t_next, x[0]);
		status = gsl_odeiv2_driver_apply(d, &t_current, t_next, x);		
		//printf("i: %d, t current: %lf, t next: %lf, x: %le \n", i, t_current, t_next, x[0]);
				
		if (status != GSL_SUCCESS) {
			printf ("error, return value=%d\n", status);
			break;
		}	
	}	
	gsl_odeiv2_driver_free(d);
	return 0;
}

int lin_reg_test() {
	//model spline for Pk
	char Pk_in_model[] = "/home/jonas/Testing_GR/code/Jonas/linear_matter_pk_sig8_0.593_z_0.75.dat";
	//Variable to hold spline information for splinting
	Pk_spline_current = input_spline_file(Pk_in_model, Pk_model_format, false);	
	//also adds regression by powerlaw for small and high k
	prep_Pk_model_regression(&Pk_spline_current, &Pk_spline_current_model);
	printf("splined Pk\n");

	double kmin_tmp = 1e-5;
	double kmax_tmp = 15.0;
	int kbins_tmp = 1000;
	double logkmin_tmp = log(kmin_tmp);
	double logkmax_tmp = log(kmax_tmp);
	double logdk_tmp = (logkmax_tmp - logkmin_tmp) / (double) (kbins_tmp-1);
	double logk_tmp, k_tmp, Pk_tmp;
	char out_tmp[] = "/home/jonas/Testing_GR/code/Jonas/output/tmp/Pk_regression_test.dat";
	FILE *f = fopen(out_tmp, "w");
	for (int i = 0; i < kbins_tmp; i++) {
		logk_tmp = logkmin_tmp + (double)i * logdk_tmp;
		k_tmp = exp(logk_tmp);
		Pk_tmp = splint_Pk_model(&Pk_spline_current_model, k_tmp);
		fprintf(f, "%le \t %le \n", k_tmp, volume*Pk_tmp);
	}
	fclose(f);
	return 0;
}

//works
int linalg_test() {
	double margin, determinant;	
	int signum;
	margin = 1e-6;
	
	gsl_matrix *M = gsl_matrix_alloc(4,4);
	gsl_vector *f = gsl_vector_alloc(4);
	gsl_vector *multipoles = gsl_vector_alloc(4);

	gsl_matrix_set(M, 0, 0, 0.18);
	gsl_matrix_set(M, 0, 1, 0.60);
	gsl_matrix_set(M, 0, 2, 0.57);
	gsl_matrix_set(M, 0, 3, 0.96);
	
	gsl_matrix_set(M, 1, 0, 0.41);
	gsl_matrix_set(M, 1, 1, 0.24);
	gsl_matrix_set(M, 1, 2, 0.99);
	gsl_matrix_set(M, 1, 3, 0.58);
	
	gsl_matrix_set(M, 2, 0, 0.14);
	gsl_matrix_set(M, 2, 1, 0.30);
	gsl_matrix_set(M, 2, 2, 0.97);
	gsl_matrix_set(M, 2, 3, 0.66);
	
	gsl_matrix_set(M, 3, 0, 0.51);
	gsl_matrix_set(M, 3, 1, 0.13);
	gsl_matrix_set(M, 3, 2, 0.19);
	gsl_matrix_set(M, 3, 3, 0.85);
	
	gsl_vector_set(f, 0, 1.0);
	gsl_vector_set(f, 1, 2.0);
	gsl_vector_set(f, 2, 3.0);
	gsl_vector_set(f, 3, 4.0);
	
	//https://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html#Linear-Algebra-Examples		
	gsl_permutation * p = gsl_permutation_alloc(4);

	gsl_linalg_LU_decomp(M, p, &signum);
	
	determinant = gsl_linalg_LU_det(M, signum);
	
	//only has a meaningful value for a nonzero determinant
	if (determinant > margin || determinant < -margin) {
		gsl_linalg_LU_solve(M, p, f, multipoles);
		printf("ans: \n %lf \n %lf \n %lf \n %lf \n", gsl_vector_get(multipoles, 0), gsl_vector_get(multipoles, 1), gsl_vector_get(multipoles, 2), gsl_vector_get(multipoles, 3));
	} else {
		printf("det: %lf \n", determinant);
	}
	
	gsl_permutation_free(p);
		
	
	gsl_matrix_free(M);
	gsl_vector_free(f);		
	gsl_vector_free(multipoles);


}
