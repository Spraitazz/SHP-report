int ZA_displ_field_1D(int dim, Pk_Spline* Pk_current, Pk_Spline* Pk_target, bool fractional, double **phases, double R_nl) {

	if (fractional && Pk_target == NULL) {
		printf("target Pk spline must be supplied for fractional displacements\n");
		exit(0);
	}
	
	//printf("\nchange ZA from Pk current to Pk target reminder \n");
	
	double kx, ky, kz, k, k_sq, coef_real, coef_im;
	double phase, amplitude, Pk, WindowFunc, prefactor, R_nl_sq;
    int index, Index, negkIndex;
 	double fundmodes[3], nyquist[3], cellSizes[3];  
 	
 	cellSizes[0] = volume_limits[0] / (double) cells_displ[0];
	cellSizes[1] = volume_limits[1] / (double) cells_displ[1];
	cellSizes[2] = volume_limits[2] / (double) cells_displ[2];

    fundmodes[0] = 2.0*pi/volume_limits[0];
    fundmodes[1] = 2.0*pi/volume_limits[1];
    fundmodes[2] = 2.0*pi/volume_limits[2];
    
    nyquist[0] = pi / cellSizes[0];
    nyquist[1] = pi / cellSizes[1];
    nyquist[2] = pi / cellSizes[2];
 	
 	R_nl_sq = R_nl * R_nl;
 	
	for (int aa = 0; aa < cells_displ[0]; aa++) {
		for (int bb = 0; bb < cells_displ[1]; bb++) {
			for (int cc = 0; cc < cells_displ[2]; cc++) {
				kx = (double)aa * fundmodes[0];
				ky = (double)bb * fundmodes[1];
				kz = (double)cc * fundmodes[2];  			
			
				//convention - up to N/2 positive frequencies, N/2 is fc, then -fc+ffund, ...
				if (aa > cells_displ[0]/2) kx -= ((double) cells_displ[0]) * fundmodes[0]; //or -= 2*kx_nyquist
				if (bb > cells_displ[1]/2) ky -= ((double) cells_displ[1]) * fundmodes[1];
				if (cc > cells_displ[2]/2) kz -= ((double) cells_displ[2]) * fundmodes[2];
			
				index = arr_ind_displ(aa, bb, cc);					
				k_sq = kx*kx + ky*ky + kz*kz;												
				
				if (k_sq > 0.0) {
				
					k = sqrt(k_sq);								
					Pk = splint_Pk(Pk_current, k);
					
					/*					
					if (Pk*volume > 1e10 || Pk*volume < 1e-10) {
						printf("\n Are you sure the pk spline is working? ZA. k: %le, Pk: %le, vol: %le, VPk: %le \n", k, Pk, volume, Pk*volume);
						exit(0);
					}
					*/
					
					//looping over dimensions, so phases prepared in advance
					phase = (*phases)[index];
								
					//amplitude from expectation, Gaussian filter in fourier space
					amplitude = sqrt(Pk) * exp(-k_sq*R_nl_sq/2.0);
					
					coef_real = amplitude * cos(phase);
					coef_im = amplitude * sin(phase);	
					
					//ints defined in main.c for preprocessor							
					switch (dim) {					
						case X:
							displacements_fourier[index][0] = kx * coef_im / k_sq;							
							displacements_fourier[index][1] = -kx * coef_real / k_sq;
							break;
						case Y:						
							displacements_fourier[index][0] = ky * coef_im / k_sq;
							displacements_fourier[index][1] = -ky * coef_real / k_sq;										
							break;
						case Z:						
							displacements_fourier[index][0] = kz * coef_im / k_sq;
							displacements_fourier[index][1] = -kz * coef_real / k_sq;							
							break;
					}
					
					if (fractional) {
						//FRACTIONAL for rescaling
						if (Pk_current != NULL) {					
							//df(k) and dv(k) here, Pk current MUST BE WITH THE RESCALED K VALUES
							prefactor = sqrt(splint_Pk(Pk_target, k)/splint_Pk(Pk_current, k)) - 1.0;					
							displacements_fourier[index][0] *= prefactor;	
							displacements_fourier[index][1] *= prefactor;	
						} else {
							printf("Pk current spline must be provided for rescaling\n");
							exit(0);							
						}
					}
								
				}						
			}
		}		
	}
	
	// be safe: intialise k=0 H_k to zero. 
	displacements_fourier[0][0] = 0.0;
	displacements_fourier[0][1] = 0.0;
	
	// Hermitian condition. One hemi-sphere is independent, e.g. k_z >= 0.
    for (int u = cells_displ[2]-1; u >= cells_displ[2]/2; u--) {
        for (int j = 0; j < cells_displ[1]; j++) {
            for (int i = 0; i < cells_displ[0]; i++) {		        
		        negkIndex = u*cells_displ[1]*cells_displ[0] + j*cells_displ[0] + i;
		        Index = 0;		           
		        // zero maps to zero on reflection through the origin.
		        if(i!=0) Index += (cells_displ[0] - i);
		        if(j!=0) Index += (cells_displ[1] - j)*cells_displ[0];
		        Index += (cells_displ[2] - u)*cells_displ[1]*cells_displ[0];			
				displacements_fourier[negkIndex][0] = displacements_fourier[Index][0];
				displacements_fourier[negkIndex][1] = -1.0*displacements_fourier[Index][1];
				if(negkIndex == Index) displacements_fourier[Index][1] = 0.0; // purely real 		    		        
			}
		}
	}
   
    // in the k_z=0 plane one semi-circle is independent, k_y>0.        
    for(int j = cells_displ[1]-1; j >= cells_displ[1]/2; j--){
        for(int i = 0; i < cells_displ[0]; i++){
            negkIndex = j*cells_displ[0] + i;                      
            Index = 0;               
            // zero maps to zero on reflection through the origin.
            if(i!=0) Index += (cells_displ[0] - i);
            Index += (cells_displ[1] - j)*cells_displ[0];
            displacements_fourier[negkIndex][0] = displacements_fourier[Index][0];
            displacements_fourier[negkIndex][1] = -1.0*displacements_fourier[Index][1];
            if(negkIndex == Index) displacements_fourier[Index][1] = 0.0;           
        }
    }
   
    // on the line k_z=k_y=0, one half is independent, k_x>=0.
    for(int i = cells_displ[0]-1; i >= cells_displ[0]/2; i--){
        negkIndex = i;                      
        Index = 0;               
        // zero maps to zero on reflection through the origin.
        Index += (cells_displ[0] - i);
        displacements_fourier[negkIndex][0] = displacements_fourier[Index][0];
        displacements_fourier[negkIndex][1] = -1.0*displacements_fourier[Index][1];
        if(negkIndex == Index) displacements_fourier[Index][1] = 0.0;       
    }  
    return 0;
}

int ZA_displacements(Catalogue *all_haloes, double*** individual_displacements, Pk_Spline *Pk_current, Pk_Spline *Pk_target, bool fractional, double R_nl) {		
	
	ZA_allocateMemory();
	int index, grid_x, grid_y, grid_z;

	//As doing each dimension separately, phases are allocated in advance. IMPORTANT!!
	double* phases = malloc((size_t)cells_displ[0]*cells_displ[1]*cells_displ[2] * sizeof(*phases));
	for (int i = 0; i < cells_displ[0]*cells_displ[1]*cells_displ[2]; i++) {
		phases[i] = randomDouble(0.0, 2.0*pi);
	}
	
	//store displacements for each particle and each dimension in a [][], apply later, check max for sanity
	double dispmax = 0.0;
	double dispmax_im = 0.0;
	for (int dimension = 0; dimension < 3; dimension++) {		
		// get the fourier displacements, inverse FFT
		ZA_displ_field_1D(dimension, Pk_current, Pk_target, fractional, &phases, R_nl);		
		fftw_execute(ZA_displacements_plan);
		
		//collect real displacements for each particle	
		for (int i = 0; i < all_haloes->particle_no; i++) {		
			grid_x = (int) floor((all_haloes->particles[i]->x / volume_limits[0]) * (double) cells_displ[0]);
			grid_y = (int) floor((all_haloes->particles[i]->y / volume_limits[1]) * (double) cells_displ[1]);
			grid_z = (int) floor((all_haloes->particles[i]->z / volume_limits[2]) * (double) cells_displ[2]);
			index = arr_ind_displ(grid_x, grid_y, grid_z);
						
			(*individual_displacements)[i][dimension] = displacements[index][0];   
			   
			//min, max displacement.     										
			if (displacements[index][0] > dispmax) dispmax = displacements[index][0];
			if (displacements[index][1] > dispmax_im) dispmax_im = displacements[index][1];
			
		}						
		ZA_clearMemory();
	}		
	
	if (dispmax_im < -1e-6 || 1e-6 < dispmax_im) printf("\n Something wrong in ZA. max disp im: %lf \n", dispmax_im);	
	free(phases);
	ZA_freeMemory();
	return 0;
}

/*
double displacement_variance(double*** individual_displacements, int particle_no) {
	//variance of |f|
	double exp_f_sq = 0.0; //E(|f|^2)
	double exp_sq_f = 0.0; //E(|f|)^2 
	double f_sq, var_f;
	
	for (int i = 0; i < particle_no; i++) {			
		f_sq = pow((*individual_displacements)[i][0], 2.0) + pow((*individual_displacements)[i][1], 2.0) + pow((*individual_displacements)[i][2], 2.0);
		exp_f_sq += f_sq;
		exp_sq_f += sqrt(f_sq);	
	}
	
	exp_f_sq /= particle_no;
	exp_sq_f /= particle_no;
	exp_sq_f = pow(exp_sq_f, 2.0);
	var_f = exp_f_sq - exp_sq_f;
	printf("PARTICLES sqrt var f: %lf \n", sqrt(var_f));

}
*/

int apply_ZA_displacements(Catalogue *catalogue, double*** individual_displacements) {
	
	double max_r_shell = 2.0*pow((volume/catalogue->particle_no)/(4.0*pi/3.0), 1.0/3.0);
	bool warned = false;
	double disp_mag = 0.0;
	double avg_disp_mag = 0.0;
	
	printf("Applying ZA displacements\n");

	for (int i = 0; i < catalogue->particle_no; i++) {			
	
		catalogue->particles[i]->x += (*individual_displacements)[i][0];
		catalogue->particles[i]->y += (*individual_displacements)[i][1];
		catalogue->particles[i]->z += (*individual_displacements)[i][2];	
		
		disp_mag = sqrt(pow((*individual_displacements)[i][0], 2.0) + pow((*individual_displacements)[i][1], 2.0) + pow((*individual_displacements)[i][2], 2.0));
		avg_disp_mag += disp_mag;	
		
		PBC(catalogue->particles[i]);	
	}
	
	printf("Particles displaced. average displacement: %lf, approximate shell cross length: %lf\n\n", avg_disp_mag/(double)catalogue->particle_no, max_r_shell);
	return 0;
}

double apply_biased_ZA_displacements(Catalogue *catalogue, double*** individual_displacements, Parameters *params, Spline *variance) {
	
	double max_r_shell = 2.0*pow((volume/catalogue->particle_no)/(4.0*pi/3.0), 1.0/3.0);
	bool warned = false;
	double disp_mag = 0.0;
	double avg_disp_mag = 0.0;
	double bias = 0.0;
	double beff_actual = 0.0;
	
	printf("Applying (BIASED) ZA displacements\n");

	for (int i = 0; i < catalogue->particle_no; i++) {
				
		bias = b_m(catalogue->particles[i]->mass, params, variance);	
		beff_actual += bias;
		
		catalogue->particles[i]->x += bias*(*individual_displacements)[i][0];
		catalogue->particles[i]->y += bias*(*individual_displacements)[i][1];
		catalogue->particles[i]->z += bias*(*individual_displacements)[i][2];	
		
		disp_mag = sqrt(pow((*individual_displacements)[i][0], 2.0) + pow((*individual_displacements)[i][1], 2.0) + pow((*individual_displacements)[i][2], 2.0));
		avg_disp_mag += disp_mag;	
		
		PBC(catalogue->particles[i]);	
	}
	
	beff_actual /= (double)catalogue->particle_no;
	
	printf("Particles displaced. average displacement: %lf, approximate shell cross length: %lf\n", avg_disp_mag/(double)catalogue->particle_no, max_r_shell);
	printf("Actual effective bias: %lf\n\n", beff_actual);
	return beff_actual;
}

int apply_ZA_velocities(Catalogue *cat, double*** individual_displacements, Parameters *params) {
	double a = z_to_a(params->z);
	double fg = pow(params->omega_m_0, params->gamma);
	double scale_factor = a * params->H0 * fg;
	printf("Applying linear velocity field, a = %lf, H0: %lf, fg = %lf\n", a, params->H0, fg);
	for (int i = 0; i < cat->particle_no; i++) {
		cat->particles[i]->vx += scale_factor * (*individual_displacements)[i][0];
		cat->particles[i]->vy += scale_factor * (*individual_displacements)[i][1];
		cat->particles[i]->vz += scale_factor * (*individual_displacements)[i][2];
	}

	return 0;
}


/*
int remove_ZA_displacements(Catalogue *catalogue, double*** individual_displacements, int redshift_space) {
	
	for (int i = 0; i < particle_no; i++) {

		catalogue->particles[i].x -= (*individual_displacements)[i][0];
		catalogue->particles[i].y -= (*individual_displacements)[i][1];
		catalogue->particles[i].z -= (*individual_displacements)[i][2];	
			
		if (redshift_space == REDSHIFT_SPACE) {		
			catalogue->particles[i].z -= fg * (*individual_displacements)[i][2];
		} else if (redshift_space != REAL_SPACE) {
			printf("expected either REAL_SPACE or REDSHIFT_SPACE for apply_ZA_displacements() \n");
			exit(0);
		}	
			
		PBC(&(catalogue->particles[i].x), &(catalogue->particles[i].y), &(catalogue->particles[i].z));		
	}

	return 0;
}

*/

/*
int ZA_RSD() {

	Foldings *foldings = NULL;

	//model spline for Pk
	char Pk_in_model[100];
	//sprintf(Pk_in_model, "%s/data/models/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	sprintf(Pk_in_model, "%s/data/models/matter_pk_om_m_055_z_40_sigma8_033.dat", home_directory);
	Spline* Pk_current = input_spline_file(Pk_in_model, Pk_model_format, NORMAL, MODEL);
	
	
	Parameters *cur_params = malloc(sizeof(*cur_params));		
	cur_params->H0 = 100.0*0.678;
	cur_params->omega_m_0 = 0.55;
	cur_params->omega_v_0 = 0.45;
	cur_params->omega_r_0 = 0.0;
	cur_params->omega = 1.0;	
	cur_params->gamma = 0.545;
	cur_params->z = 4.0;	
	
	
	Parameters *cur_params = malloc(sizeof(*cur_params));		
	cur_params->H0 = 100.0*0.678;
	cur_params->omega_m_0 = 0.24;
	cur_params->omega_v_0 = 0.76;
	cur_params->omega_r_0 = 0.0;
	cur_params->omega = 1.0;	
	cur_params->gamma = 0.545;
	cur_params->z = 3.5;	


	//binning
	int k_bins = 500;	
	Dplus_bins = 500;
	int z_bins = 500;
	double k_min = Pk_current->splineInfo->xmin;
	double k_max = Pk_current->splineInfo->xmax;	
	double z_min = 0.0;
	double z_max = 4.0;
	double R_min = 0.05;
	double R_max = 6.5;
	BinInfo *rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	BinInfo *rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);

	//redshift 0.75 -> 3
	//Pk_current = prep_Pk_constz(GR, 0.75, 3.5, Pk_current, rescaling_k_bin_info, rescaling_z_bin_info, cur_params);	
	
	BinInfo *rescaling_R_bin_info = prep_bins(R_min, R_max, 500, LOG_BIN);

	Spline *variance_spline_current = prep_variances(rescaling_R_bin_info, Pk_current);
	Spline *variance_spline_current_reversed = reverse_spline(variance_spline_current);
	
	R_nl = splint_generic(variance_spline_current_reversed, 1.0);	
	printf("Rnl is: %le \n", R_nl);
	
	//sigma_exp_smoothed = sqrt(expected_variance_smoothed(R_nl, Pk_current_extended));
	//printf("expected standard deviation of smoothed displacement: %lf \n", sigma_exp_smoothed);				
	
	char outFile_redshift[100], out_model_redshift[100];
	sprintf(out_model_redshift, "%s/data/rescaling/model_to_model/Pk_current_model.dat", home_directory);
		
	//k bins for models
	BinInfo* model_bin_info = prep_bins(0.02, 0.4, 200, LOG_BIN);
	
	int particle_no = 128*128*128;	

	haloModel_out(out_model_redshift, model_bin_info, 0.0, 1.0, REDSHIFT_SPACE, Pk_current, cur_params, 1.0);	
	
	Catalogue *all_haloes;
	
	double** individual_displacements;	
	
	for (int i = 0; i < 20; i++) {
	
		printf("\n RUN %d \n", i);	
		
		sprintf(outFile_redshift, "%s/data/rescaling/model_to_model/test/Pk_current_%d.dat", home_directory, i);		

		all_haloes = populateParticles_crystalLattice(particle_no);		
		
		calloc2D_double(&individual_displacements, particle_no, 3);
		ZA_allocateMemory();
		ZA_displacements(all_haloes, &individual_displacements, NULL, Pk_current, false);		
		ZA_freeMemory();		
	
		apply_ZA_displacements(all_haloes, &individual_displacements, REAL_SPACE);
		apply_ZA_velocities(all_haloes, &individual_displacements, cur_params);
		toRedshift(all_haloes, cur_params);
		
		overdensity_allocateMemory();
		haloes_measure_Pk(all_haloes, outFile_redshift, 1.0, CIC, 1.0);		
		overdensity_freeMemory();
		
		free2D_double(&individual_displacements, particle_no);	
		free(all_haloes->particles);
		free(all_haloes);		
	}

	
	return 0;
}
*/


/*

int ST_bias() {

	char Pk_current_path[100];
	sprintf(Pk_current_path, "%s/data/models/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	Spline *Pk_current = input_spline_file(Pk_current_path, Pk_model_format, NORMAL, MODEL);
	
	Parameters *cur_params = malloc(sizeof(*cur_params));		
	cur_params->H0 = 100.0*0.678;
	cur_params->omega_m_0 = 0.24;
	cur_params->omega_v_0 = 0.76;
	cur_params->omega_r_0 = 0.0;
	cur_params->omega = 1.0;	
	cur_params->gamma = 0.545;
	cur_params->z = 3.5;	

	//binning
	int k_bins = 500;	
	Dplus_bins = 500;
	int z_bins = 500;
	double k_min = Pk_current->splineInfo->xmin;
	double k_max = Pk_current->splineInfo->xmax;	
	double z_min = 0.0;
	double z_max = 4.0;
	
	double Mmin = 1e10;
	double Mmax = 1e12;
	double R1 = m_to_R(Mmin, cur_params);
	double R2 = m_to_R(Mmax, cur_params);
	BinInfo *rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	BinInfo *rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	BinInfo *rescaling_R_bin_info = prep_bins(R1, R2, 500, LOG_BIN);
	
	Spline *Dplus_z_cur = Dplus_spline(GR, 0.0, rescaling_z_bin_info, cur_params);
	Pk_current = prep_Pk_constz(GR, Dplus_z_cur, 0.75, cur_params->z, Pk_current, rescaling_k_bin_info, rescaling_z_bin_info, cur_params);

	Spline *variance_current = prep_variances(rescaling_R_bin_info, Pk_current);
	R_nl = R_nonlinear(variance_current);
	printf("Rnl is: %le \n", R_nl);
	
	//sigma_exp_smoothed = sqrt(expected_variance_smoothed(R_nl, Pk_current_extended));
	//printf("expected standard deviation of smoothed displacement: %lf \n", sigma_exp_smoothed);				
	
	char outFile[100], outFile2[100], outFile3[100];
	char out_model_beff[100], out_model_bvmin[100], out_model_bvmax[100];
	sprintf(out_model_beff, "%s/data/bias_disp/model_pk_beff.dat", home_directory);
	sprintf(out_model_bvmin, "%s/data/bias_disp/model_pk_bvmin.dat", home_directory);
	sprintf(out_model_bvmax, "%s/data/bias_disp/model_pk_bvmax.dat", home_directory);
	
	//k bins for models
	BinInfo* model_bin_info = prep_bins(0.02, 1.0, 500, LOG_BIN);
	
	int particle_no = 128*128*128;	
	
	double fv_min = f_v(m_to_v(Mmin, cur_params, variance_current));
	double fv_max = f_v(m_to_v(Mmax, cur_params, variance_current));
	int Nbig = (int) particle_no/(1.0 + (fv_min/fv_max));
	int Nsmall = particle_no - Nbig;	

	double b_eff = calc_b_eff(cur_params, variance_current, Mmin, Mmax);
	double bv_min = b_m(Mmin, cur_params, variance_current);
	double bv_max = b_m(Mmax, cur_params, variance_current);
	printf("frac: %lf, particles: %d, beff: %lf, Nbig: %d, bv_Mmax: %lf, Nsmall: %d, bv_Mmin: %lf \n", fv_min/fv_max, particle_no, b_eff, Nbig, bv_max, Nsmall, bv_min);	
	
	//haloModel_out(out_model_beff, model_bin_info, volume/particle_no, pow(b_eff, 2.0), REAL_SPACE, Pk_current_extended);	
	//haloModel_out(out_model_bvmin, model_bin_info, volume/Nsmall, pow(bv_min, 2.0), REAL_SPACE, Pk_current_extended);
	//haloModel_out(out_model_bvmax, model_bin_info, volume/Nbig, pow(bv_max, 2.0), REAL_SPACE, Pk_current_extended);
	
	double** individual_displacements;
	
	
	int particles_gridded;
	
	for (int i = 0; i < 1; i++) {
	
		printf("\n RUN %d \n", i);	
		sprintf(outFile, "%s/data/tmp/ZA_debiased_pk_%d.dat", home_directory, i);
		sprintf(outFile2, "%s/data/tmp/ZA_debiased_1e13_pk_%d.dat", home_directory, i);
		sprintf(outFile3, "%s/data/tmp/ZA_debiased_1e14_pk_%d.dat", home_directory, i);

		all_haloes = populateHaloes_twoMasses(Mmin, Mmax, particle_no, Nbig);
		printf("populated all haloes, all: %d \n", all_haloes->particle_no);
		
		calloc2D_double(&individual_displacements, particle_no, 3);

		ZA_displacements(all_haloes, &individual_displacements, NULL, Pk_current, false);

		mass_bias_displacements(all_haloes, &individual_displacements, cur_params, variance_current);		
		//printf("Applying ZA displacement field \n");
		//apply_ZA_displacements(all_haloes, &individual_displacements, REAL_SPACE);			
	
		//overdensity_allocateMemory();
		//printf("Measuring Pk \n");
		
		//all haloes
		haloes_measure_Pk(all_haloes, outFile, 1.0, CIC);		
		
		
		//1e13 mass Pk measurement
		particles_gridded = populateGridCIC_massRange(all_haloes, 0.99*Mmin, 1.01*Mmin);	
		gridIntoOverdensity_particleNo(particles_gridded);	
		overdensity_pk();
		pk_to_file(outFile2);
		clearGrid();
		overdensity_clearMemory();	
		free_spectrum_storage();
		
		//1e14 mass Pk measurement
		particles_gridded = populateGridCIC_massRange(all_haloes, 0.99*Mmax, 1.01*Mmax);	
		gridIntoOverdensity_particleNo(particles_gridded);	
		overdensity_pk();
		pk_to_file(outFile3);
		clearGrid();
		overdensity_clearMemory();	
		free_spectrum_storage();	
		

		overdensity_freeMemory();		
		free(all_haloes->particles);
		free(all_haloes);
		
		free2D_double(&individual_displacements, particle_no);				
	}
	

	return 0;	
}

*/

/*
int biased_displacements() {

	//model spline for Pk
	char Pk_in_model[100];
	sprintf(Pk_in_model, "%s/data/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	SplineInfo* Pk_current = input_spline_file(Pk_in_model, Pk_model_format, NORMAL);

	//binning
	k_bins = 500;	
	Dplus_bins = 500;
	z_bins = 500;
	k_min = Pk_current->xmin;
	k_max = Pk_current->xmax;	
	z_min = 0.0;
	z_max = 5.0;
	double R_min = 0.05;
	double R_max = 6.5;
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	
	//PARAMETERS
	z_current_global = 5.0;	

	//for Dplus things
	z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	prep_redshift_spline(rescaling_z_bin_info);

	//redshift 0.75 -> z_current_global
	Pk_current = prep_Pk_constz(GR, 0.75, z_current_global, Pk_current, rescaling_k_bin_info);	
	Pk_current_extended = extend_spline_model(Pk_current);
	
	rescaling_R_bin_info = prep_bins(R_min, R_max, 500, LOG_BIN);
	OLV_workspace = gsl_integration_workspace_alloc(OLV_workspace_size);
	variance_spline_current = prep_variances(rescaling_R_bin_info, Pk_current_extended);
	variance_spline_current_reversed = reverse_spline(variance_spline_current);
	
	R_nl = splint_generic(variance_spline_current_reversed, 1.0);	
	printf("Rnl is: %le \n", R_nl);
	
	sigma_exp_smoothed = sqrt(expected_variance_smoothed(R_nl, Pk_current_extended));
	printf("expected standard deviation of smoothed displacement: %lf \n", sigma_exp_smoothed);				
	
	char outFile[100], outFile2[100], out_model[100], out_model_debiased[100];
	sprintf(out_model, "%s/data/bias_disp/model_pk.dat", home_directory);
	sprintf(out_model_debiased, "%s/data/bias_disp/model_pk_debiased.dat", home_directory);
	
	particle_no = 500000;
	halo_no = particle_no;
	initParticles();
	
	Mmin = 1e12;
	Mmax = 1e13;
	
	double f_big = f_v(m_to_v(Mmax, z_current_global, variance_spline_current));
	double f_small = f_v(m_to_v(Mmin, z_current_global, variance_spline_current));
	
	double frac = f_small / f_big;
	
	int no_big = particle_no / (1.0 + frac);
	
	//double b_eff_actual = 0.5*(b_m(Mmin, z_current_global, variance_spline_current) + b_m(Mmax, z_current_global, variance_spline_current));
	
	//printf("expect b eff: %lf \n", b_eff_actual);	
	
	//calc_b_eff(z_current_global, variance_spline_current);
	//b_eff = b_eff_actual;
	
	//k bins for models
	BinInfo* model_bin_info = prep_bins(0.02, 1.0, 500, LOG_BIN);	
	
	double avg_bias = 0.0;
	
	double** individual_displacements;
	
	char inCat[100];
	sprintf(inCat, "%s/AllHalo_GR_1500Mpc_a0.7_1.txt", mock_directory);
	
	for (int i = 0; i < 50; i++) {
	
		printf("\n RUN %d \n", i);
	
		//sprintf(outFile, "%s/data/bias_disp/debiased/ZA_debiased_pk_%d.dat", home_directory, i);
		sprintf(outFile, "%s/data/bias_disp/debiased/ZA_debiased_pk_%d.dat", home_directory, i);
		sprintf(outFile2, "%s/data/bias_disp/debiased/ZA_debiased_pk_4fold_%d.dat", home_directory, i);
	
		//populateHaloes_catalogueMass(inCat);
		populateHaloes_diffMass(Mmin, Mmax, no_big);
		
		if (m_to_R(Mmin, z_current_global) < R_min || m_to_R(Mmax, z_current_global) > R_max) {
			printf("expand R bounds to Rmin: %lf, Rmax: %lf \n", m_to_R(Mmin, z_current_global), m_to_R(Mmax, z_current_global));
		}
		
		calc_b_eff(z_current_global, variance_spline_current);	
		
		calloc2D_double(&individual_displacements, particle_no, 3);
		
		ZA_allocateMemory();
		ZA_displacements(&individual_displacements, Pk_current_extended, NGP);		
		ZA_freeMemory();
		
		mass_debias_displacements(&individual_displacements);
		avg_bias += effective_bias_calculated;
		printf("effective bias here: %lf, avg bias: %lf \n", effective_bias_calculated, avg_bias);
		printf("Applying ZA displacement field \n");
		apply_ZA_displacements(&individual_displacements, REAL_SPACE);
		
		free2D_double(&individual_displacements, particle_no);
	
		overdensity_allocateMemory();
		printf("Measuring Pk \n");
		haloes_measure_Pk(outFile, 1.0, NGP);	
		haloes_measure_Pk(outFile2, 4.0, NGP);
		overdensity_freeMemory();
	
		clearParticles();
		//freeParticles();
	
	}
	
	avg_bias /= 50.0;
	
	//debiased
	haloModel_out(out_model_debiased, model_bin_info, volume/particle_no, pow(b_eff, -2.0), REAL_SPACE, Pk_current_extended);
	
	//normal
	haloModel_out(out_model, model_bin_info, volume/particle_no, pow(avg_bias, -2.0), REAL_SPACE, Pk_current_extended);		
	
	
}

*/

/*

int ZA_test() {

	//model spline for Pk
	char Pk_in_model[100];
	sprintf(Pk_in_model, "%s/data/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	SplineInfo* Pk_current = input_spline_file(Pk_in_model, Pk_model_format, NORMAL);		
	
	//PARAMETERS
	set_params();
	z_current_global = 3.5;	
	
	//binning
	k_bins = 300;	
	Dplus_bins = 1000;
	k_min = 0.99*Pk_current.xmin;
	k_max = 1.01*Pk_current.xmax;	
	z_min = 0.0;
	z_max = 5.0;
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, true); //logbin	
	
	//redshift 0.75 -> 3.5
	prep_Pk_constz(GR, 0.75, z_current_global, Pk_current, rescaling_k_bin_info);	
	
	//WHY USING REGRESSION EXTENSION WORKS, BUT WITHOUT IT DOESNT????
	Pk_current_extended = extend_spline_model(&Pk_current);
	//extend_spline_nonmodel(&Pk_current, &Pk_spline_current);
	
	rescaling_R_bin_info = prep_bins(0.001, 10.0, 300, true); //logbin
	variance_spline_current = prep_variances(rescaling_R_bin_info, Pk_current_extended);
	
	reverse_spline(variance_spline_current, variance_spline_current_reversed);
	//R_nl is defined for a given z such that sigma^2_lin(R_nl, z) = 1	
	R_nl = splint_generic(variance_spline_current_reversed, 1.0);	
	printf("Rnl is: %le \n", R_nl);	
	
	sigma_exp_smoothed = sqrt(expected_variance_smoothed(R_nl, Pk_current_extended));
	printf("expected standard deviation of smoothed displacement: %lf \n", sigma_exp_smoothed);			
	
	prep_oneHalo();				
	
	char outFile[100];
	char tmpFile[100];
	char tmpFile2[100];
	sprintf(tmpFile, "%s/output/HOD_rands/unfolded.dat", home_directory);
	sprintf(tmpFile2, "%s/output/HOD_rands/folded.dat", home_directory);
	
	double** individual_displacements;
	calloc2D_double(&individual_displacements, particle_no, 3);	

	
	for (int i = 0; i < 100; i++) {
		sprintf(outFile, "%s/output/HOD_rands/test/measured_pk_%d.dat", home_directory, i);
		printf("run %d \n", i);
		
		initParticles();
		//populateParticles_lattice();
		populateParticles();		

		//ZA_allocateMemory();
		ZA_displacements(&individual_displacements, &Pk_current, NGP);
		apply_ZA_displacements(&individual_displacements, REAL_SPACE);
		//apply_ZA(false);
		//ZA_freeMemory();
		
		halo_no = particle_no; //now particle number will change, stay different for measuring pk, then revert
	//	oneHalo();
		
		overdensity_allocateMemory();
		haloes_measure_Pk(tmpFile, 1.0, false);	
		haloes_measure_Pk(tmpFile2, 4.0, false);
		combine_folded(tmpFile, tmpFile2, pk_format_folding, outFile, 0.0, 1.1, 0.3, true);		
		overdensity_freeMemory();
		
		clearParticles();
		remove(tmpFile);
		remove(tmpFile2);
		
		particle_no = halo_no;
		free2D_double(&particles, particle_no);
		
	}
	
	free2D_double(&individual_displacements, particle_no);
	//gsl_integration_workspace_free(onehalo_workspace);
	

	return 0;
}
*/
