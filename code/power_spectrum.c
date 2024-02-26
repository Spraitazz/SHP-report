int toRedshift(Catalogue *cat, Parameters *params) {
	double s_prefac = z_to_a(params->z)*(params->H0);
	printf("Going to redshift, z + vz/%lf\n", s_prefac);
	for (int i = 0; i < cat->particle_no; i++) {		
		cat->particles[i]->z += cat->particles[i]->vz / s_prefac;	//z-axis line of sight		
		PBC(cat->particles[i]);
	}	
	return 0;
}

int toRealspace(Catalogue *cat, Parameters *params) {
	for (int i = 0; i < cat->particle_no; i++) {		
		cat->particles[i]->z -= cat->particles[i]->vz / (z_to_a(params->z)*(params->H0));	//z-axis line of sight		
		PBC(cat->particles[i]);
	}	
	return 0;
}

int populateGrid(Catalogue *catalogue, int grid_func) {
	if (grid_func == CIC) {
		populateGridCIC(catalogue);
	} else if (grid_func == NGP) {
		populateGridNGP(catalogue);	
	} else {
		printf("grid func neither NGP nor CIC in populateGrid() \n");
		exit(0);
	}
	return 0;
}

int populateGridNGP(Catalogue *catalogue) {
	int grid_x, grid_y, grid_z, index;		
	for (int i = 0; i < catalogue->particle_no; i++) {

		grid_x = (int) floor((catalogue->particles[i]->x / volume_limits[0]) * (double) cells[0]);
		grid_y = (int) floor((catalogue->particles[i]->y / volume_limits[1]) * (double) cells[1]);
		grid_z = (int) floor((catalogue->particles[i]->z / volume_limits[2]) * (double) cells[2]);		
		index = arr_ind(grid_x, grid_y, grid_z);
		grid[index] += 1.0;	
		//printf("gridding, index: %d, grid[index]: %lf ", index, grid[index]);
		//print_particle(catalogue->particles[i]);	
	}
	WindowFunc_pow = 1.0;
	return 0;
}
/*
int populateGrid_massRange(double M_min, double M_max) {
	int grid_x, grid_y, grid_z, index, particles_gridded;
	particles_gridded = 0;			
	for (int i = 0; i < particle_no; i++) {
		if (M_max*(1.000001) > particles[i][6] && particles[i][6] > M_min*(0.999999)) {
			particles_gridded += 1;
			grid_x = (int) floor((particles[i][0] / volume_limits[0]) * (double) cells[0]);
			grid_y = (int) floor((particles[i][1] / volume_limits[1]) * (double) cells[1]);
			grid_z = (int) floor((particles[i][2] / volume_limits[2]) * (double) cells[2]);		
			index = arr_ind(grid_x, grid_y, grid_z);
			grid[index] += 1.0;	
		}	
	}
	WindowFunc_pow = 1.0;
	return particles_gridded;
}
*/
//http://background.uchicago.edu/~whu/Courses/Ast321_11/pm.pdf p14
int populateGridCIC(Catalogue *catalogue) {

	int grid_x, grid_y, grid_z;
	double xc, yc, zc, dx, dy, dz, tx, ty, tz;
	double cellSizes[3];
	cellSizes[0] = volume_limits[0] / (double) cells[0];
	cellSizes[1] = volume_limits[1] / (double) cells[1];
	cellSizes[2] = volume_limits[2] / (double) cells[2];	
	
	for (int i = 0; i < catalogue->particle_no; i++) {
		grid_x = (int) floor((catalogue->particles[i]->x / volume_limits[0]) * (double) cells[0]);
		grid_y = (int) floor((catalogue->particles[i]->y / volume_limits[1]) * (double) cells[1]);
		grid_z = (int) floor((catalogue->particles[i]->z / volume_limits[2]) * (double) cells[2]);		
		
		xc = cellSizes[0] * ((double) grid_x + 0.5);
		yc = cellSizes[1] * ((double) grid_y + 0.5);
		zc = cellSizes[2] * ((double) grid_z + 0.5);
		
		//must be from 0 to 1
		dx = fabs(catalogue->particles[i]->x - xc) / cellSizes[0];
		dy = fabs(catalogue->particles[i]->y - yc) / cellSizes[1];
		dz = fabs(catalogue->particles[i]->z - zc) / cellSizes[2];	
		
		tx = 1.0 - dx;
		ty = 1.0 - dy;
		tz = 1.0 - dz;
		
		//ensure PBC
		int xplus, yplus, zplus;
		xplus = (grid_x + 1) % cells[0];
		yplus = (grid_y + 1) % cells[1];
		zplus = (grid_z + 1) % cells[2];

		grid[arr_ind(grid_x, grid_y, grid_z)] += tx*ty*tz;	
		grid[arr_ind(xplus, grid_y, grid_z)] += dx*ty*tz;
		grid[arr_ind(xplus, yplus, grid_z)] += dx*dy*tz;	
		grid[arr_ind(xplus, grid_y, zplus)] += dx*ty*dz;
		grid[arr_ind(xplus, yplus, zplus)] += dx*dy*dz;			
		grid[arr_ind(grid_x, yplus, grid_z)] += tx*dy*tz;		
		grid[arr_ind(grid_x, grid_y, zplus)] += tx*ty*dz;
		grid[arr_ind(grid_x, yplus, zplus)] += tx*dy*dz;		
	}
	WindowFunc_pow = 2.0;
	return 0;
}
/*
int populateGridCIC_massRange(Catalogue* cat, double M_min, double M_max) {
	int grid_x, grid_y, grid_z;
	double xc, yc, zc, dx, dy, dz, tx, ty, tz;
	double cellSizes[3];
	cellSizes[0] = volume_limits[0] / (double) cells[0];
	cellSizes[1] = volume_limits[1] / (double) cells[1];
	cellSizes[2] = volume_limits[2] / (double) cells[2];
	
	int particles_gridded = 0;			
	
	for (int i = 0; i < cat->particle_no; i++) {
	
		if (M_max*(1.000001) > cat->particles[i].mass && cat->particles[i].mass > M_min*(0.999999)) {
			
			particles_gridded += 1;
			
			grid_x = (int) floor((cat->particles[i].x / volume_limits[0]) * (double) cells[0]);
			grid_y = (int) floor((cat->particles[i].y / volume_limits[1]) * (double) cells[1]);
			grid_z = (int) floor((cat->particles[i].z / volume_limits[2]) * (double) cells[2]);			
		
			xc = cellSizes[0] * ((double) grid_x + 0.5);
			yc = cellSizes[1] * ((double) grid_y + 0.5);
			zc = cellSizes[2] * ((double) grid_z + 0.5);
		
			//must be from 0 to 1
			dx = fabs(cat->particles[i].x - xc) / cellSizes[0];
			dy = fabs(cat->particles[i].y - yc) / cellSizes[1];
			dz = fabs(cat->particles[i].z - zc) / cellSizes[2];	
		
			tx = 1.0 - dx;
			ty = 1.0 - dy;
			tz = 1.0 - dz;
		
			//ensure PBC
			int xplus, yplus, zplus;
			xplus = (grid_x + 1) % cells[0];
			yplus = (grid_y + 1) % cells[1];
			zplus = (grid_z + 1) % cells[2];		

			grid[arr_ind(grid_x, grid_y, grid_z)] += tx*ty*tz;	
			grid[arr_ind(xplus, grid_y, grid_z)] += dx*ty*tz;
			grid[arr_ind(xplus, yplus, grid_z)] += dx*dy*tz;	
			grid[arr_ind(xplus, grid_y, zplus)] += dx*ty*dz;
			grid[arr_ind(xplus, yplus, zplus)] += dx*dy*dz;			
			grid[arr_ind(grid_x, yplus, grid_z)] += tx*dy*tz;		
			grid[arr_ind(grid_x, grid_y, zplus)] += tx*ty*dz;
			grid[arr_ind(grid_x, yplus, zplus)] += tx*dy*dz;
		
		}
			
	}
	WindowFunc_pow = 2.0;
	return particles_gridded;
}
*/

int gridIntoOverdensity(int particle_no) {
	//average density (particles per grid box) simply particle_no / size of box
	double N_avg = (double) particle_no / ((double) cells[0]*cells[1]*cells[2]);	
	int index;	
	for (int i = 0; i < cells[0]; i++) {
		for (int j = 0; j < cells[1]; j++) {
			for (int k = 0; k < cells[2]; k++) {	
				//density of particles of a single box = number of particles in that box
				//overdensity = number of particles in box / rho_avg - 1	
				index = arr_ind(i, j, k);				
				overdensity[index][0] = (grid[index] / N_avg) - 1.0;			
			}
		}
	}
	return 0;
}

//gaussian here just means use overdensity fourier coeffs from a gaussian random field (NOT IMPLEMENTED)
int overdensity_pk(double multip_fact) {
    double kx, ky, kz, k, logdk, logkmin, logkmax, WindowFunc, mu, kmin, kmax;
    double fundmodes[3], nyquist[3], cellSizes[3];
    double phase, amplitude;
    bool flag;
    int index;  
    
    cellSizes[0] = volume_limits[0] / (double) cells[0];
	cellSizes[1] = volume_limits[1] / (double) cells[1];
	cellSizes[2] = volume_limits[2] / (double) cells[2];	

    fundmodes[0] = 2.0*pi/volume_limits[0];
    fundmodes[1] = 2.0*pi/volume_limits[1];
    fundmodes[2] = 2.0*pi/volume_limits[2];
    
    nyquist[0] = pi / cellSizes[0];
    nyquist[1] = pi / cellSizes[1];
    nyquist[2] = pi / cellSizes[2];
       
    //should be min(fundmodes), different law for freqs 
    kmin = fundmodes[0] * sqrt(3.0);
    kmax = nyquist[0] * sqrt(3.0) / 4.0;
    
    //get fourier overdensities
	fftw_execute(overdensity_plan);		
    
    BinInfo* Pk_binInfo = prep_bins(kmin, kmax, spectrum_size, LOG_BIN);
    
    //init_spectrum_storage(spectrum_size);
    
    printf("Pk bins: %d, kmin: %lf, kmax: %lf, cell size: %lf\n", spectrum_size, kmin, kmax, cellSizes[0]);

    for (int aa = 0; aa < cells[0]; aa++) {
    	for (int bb = 0; bb < cells[1]; bb++) {
    		for (int cc = 0; cc < cells[2]; cc++) {
    			kx = (double)aa * fundmodes[0];
    			ky = (double)bb * fundmodes[1];
				kz = (double)cc * fundmodes[2];  			
				
				//convention - up to N/2 positive frequencies, N/2 is fc, then -fc+ffund, ...
				if (aa > cells[0]/2) kx -= ((double) cells[0]) * fundmodes[0]; //or -= 2*kx_nyquist
				if (bb > cells[1]/2) ky -= ((double) cells[1]) * fundmodes[1];
				if (cc > cells[2]/2) kz -= ((double) cells[2]) * fundmodes[2];
				
				index = arr_ind(aa, bb, cc);					
				
				overdensity_fourier[index][0] /= ((double)cells[0]*cells[1]*cells[2]);
				overdensity_fourier[index][1] /= ((double)cells[0]*cells[1]*cells[2]);
				
				//at different redshift, for rescaling
				overdensity_fourier[index][0] *= multip_fact;
				overdensity_fourier[index][1] *= multip_fact;
								
				//Window function (convolution in real -> multip in fourier) to reduce aliased power
				WindowFunc = 1.0;
				if (kx != 0.0) WindowFunc *= sin(pi*kx*0.5/nyquist[0])/(pi*kx*0.5/nyquist[0]);
				if (ky != 0.0) WindowFunc *= sin(pi*ky*0.5/nyquist[1])/(pi*ky*0.5/nyquist[1]);
				if (kz != 0.0) WindowFunc *= sin(pi*kz*0.5/nyquist[2])/(pi*kz*0.5/nyquist[2]);
									
				// Correct for mass assignment of randoms.
				// Cloud in cell pow 2.0. NGP corresponds to WindowFunc rather than WindowFunc^2. 
				overdensity_fourier[index][0] /= pow(WindowFunc, WindowFunc_pow);
				overdensity_fourier[index][1] /= pow(WindowFunc, WindowFunc_pow);
				
				k = sqrt(kx*kx + ky*ky + kz*kz);				
				//unit vector along k dot n, chosen along polar axis kz
				mu = kz / k;			
				
				//one half independent. 
				flag = false;				
				if (kz > 0.0) {	flag = true; }
				else if ((kz == 0.0) && (ky > 0.0)) { flag = true; } 
				else if ((kz == 0.0) && (ky == 0.0) && (kx > 0.0)) { flag = true; } 
				
				if (flag && k >= kmin && k < kmax) pk_logbin(index, k, mu, Pk_binInfo);				
											  		
    		}    	
    	}    
	}
								        
	//find the average k for each bin - where the P(k) values are plotted, check for empty bins
	for (int i = 0; i < spectrum_size; i++) {
		if (k_number[i] != 0) {
			bin_k_average[i] /= (double) k_number[i];
		} else {
			toRemove_bins[i] = true;
		}
	}		
	//finally, calculate multipoles (mono, quadro)
	multipole_calc();	
	return 0;
}

int pk_logbin(int overdensity_index, double k, double mu, BinInfo* Pk_binInfo) {
	double fk_squared, coef_real, coef_im, L2, L4;
	int index_log_pk;	

	coef_real = overdensity_fourier[overdensity_index][0];
	coef_im = overdensity_fourier[overdensity_index][1];		
	fk_squared = coef_real * coef_real + coef_im * coef_im;	
	
	//plot is VP(k) versus k
	fk_squared *= volume;
	
	//SHOT NOISE, USE CAREFULLY, only monopole has shotnoise
	//fk_squared -= volume/particle_no;	
		
	L2 = 0.5 * (3.0 * mu * mu - 1.0);
	L4 = (1.0/8.0) * (35.0 * pow(mu, 4.0) - 30.0 * pow(mu, 2.0) + 3.0);	
							
	index_log_pk = x_to_bin(Pk_binInfo, k);	
	
	fk_abs_squares[index_log_pk] += fk_squared;	
	k_number[index_log_pk] += 1;	
	bin_k_average[index_log_pk] += k;
	
	if (k > bin_k_max[index_log_pk]) bin_k_max[index_log_pk] = k;
	if (k < bin_k_min[index_log_pk]) bin_k_min[index_log_pk] = k;	
	
	binsums_L2[index_log_pk][0] += L2;
	binsums_L2[index_log_pk][1] += L2 * L2;	
	binsums_L2[index_log_pk][2] += L2 * fk_squared;	
	
	binsums_L4[index_log_pk][0] += L4;
	binsums_L4[index_log_pk][1] += L4 * L4;
	binsums_L4[index_log_pk][2] += L2 * L4;
	binsums_L4[index_log_pk][3] += L4 * fk_squared;

	return 0;
}

//monopole + qudrupole + hexadecapole (pages 87, 88 of thesis)
//M p(vec) = f(vec)
int multipole_calc() {	
	double margin, determinant;	
	int signum;
	margin = 1e-6;
	
	gsl_matrix *M = gsl_matrix_alloc(3,3);
	gsl_vector *f = gsl_vector_alloc(3);
	gsl_vector *multipoles = gsl_vector_alloc(3);
	
	for (int i = 0; i < spectrum_size; i++) {			
	
		gsl_matrix_set(M, 0, 0, (double) k_number[i]);
		gsl_matrix_set(M, 0, 1, binsums_L2[i][0]);
		gsl_matrix_set(M, 0, 2, binsums_L4[i][0]);
		
		gsl_matrix_set(M, 1, 0, binsums_L2[i][0]);
		gsl_matrix_set(M, 1, 1, binsums_L2[i][1]);
		gsl_matrix_set(M, 1, 2, binsums_L4[i][2]);
		
		gsl_matrix_set(M, 2, 0, binsums_L4[i][0]);
		gsl_matrix_set(M, 2, 1, binsums_L4[i][2]);
		gsl_matrix_set(M, 2, 2, binsums_L4[i][1]);
		
		gsl_vector_set(f, 0, fk_abs_squares[i]);
		gsl_vector_set(f, 1, binsums_L2[i][2]);
		gsl_vector_set(f, 2, binsums_L4[i][3]);
		
		//https://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html#Linear-Algebra-Examples		
		gsl_permutation * p = gsl_permutation_alloc(3);

		gsl_linalg_LU_decomp(M, p, &signum);
		
		determinant = gsl_linalg_LU_det(M, signum);
		
		//only has a meaningful value for a nonzero determinant
		if (determinant > margin || determinant < -margin) {
			gsl_linalg_LU_solve(M, p, f, multipoles);
			monopole[i] = gsl_vector_get(multipoles, 0);
			quadrupole[i] = gsl_vector_get(multipoles, 1);
			hexadecapole[i] = gsl_vector_get(multipoles, 2);			
		} else {
			toRemove_bins[i] = true;
		}
		
		gsl_permutation_free(p);
		
	}
	gsl_matrix_free(M);
	gsl_vector_free(f);		
	gsl_vector_free(multipoles);
	return 0;
}


//prints k, monopole (simple calculation), measurements in bin
int just_print_monopole(char filename[]) {
	FILE *pk_out_file = fopen(filename, "w");	
	for (int i = 0; i < spectrum_size; i++) {
		if (toRemove_bins[i] == false) {
			fprintf(pk_out_file, "%lf \t %le \t %d \t %lf \t %lf\n", bin_k_average[i], fk_abs_squares[i]/(double)k_number[i], k_number[i], bin_k_min[i], bin_k_max[i]);
		}
	}
	fclose(pk_out_file);
	return 0;
}


//prints k, monopole, quadrupole, measurements in bin
int pk_to_file(const char fileName[], double shotnoise) {
	FILE *f = fopen(fileName, "w");
	if (f == NULL) {
		printf("no such file %s. pk_to_file() \n", fileName);
		exit(0);
	}
	printf("Printing P(k) to %s\n", fileName);
	for (int i = 0; i < spectrum_size; i++) {
		if (toRemove_bins[i] == false && monopole[i] > 1e-6) {			
			fprintf(f, pk_format, bin_k_average[i], monopole[i] - shotnoise, quadrupole[i], hexadecapole[i], k_number[i], bin_k_min[i], bin_k_max[i]);
			fprintf(f, "\n");
		} 
	}	
	fclose(f);
	return 0;
}

int delsq_to_file(const char filename[]) {
	FILE *f = fopen(filename, "w");	
	double delsq;
	for (int i = 0; i < spectrum_size; i++) {				
		if (toRemove_bins[i] == false) {
			delsq = pow(bin_k_average[i], 3.0)*monopole[i]/(2.0*pi*pi);
			fprintf(f, pk_format, bin_k_average[i], delsq, quadrupole[i], hexadecapole[i], k_number[i], bin_k_min[i], bin_k_max[i]);
			fprintf(f, "\n");
		}	
	}
	fclose(f);
	return 0;
}


