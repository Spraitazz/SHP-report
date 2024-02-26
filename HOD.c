
/*
Takes the halo mass, the splined linear variance in overdensity (make sure range is enough)
and the cosmology parameters.

Returns the cdm, rvir and rs
*/
HOD_Parameters *HOD_params(double halo_mass, Spline *variance, Parameters *params) {
	
	Spline *variance_reversed = reverse_spline(variance);
	double R_mstar = splint_generic(variance_reversed, pow(delta_c, 2.0));
	double Mstar = R_to_m(R_mstar, params);
	
	if (R_mstar < 1e-6 || isnan(R_mstar)) {
		printf("R mstar %lf, Mstar %le in HOD_Parameters(). Check variance spline range? Min variance: %lf, corresponding r: %lf, max variance: %lf, corresponding r: %lf \n", R_mstar, Mstar, variance_reversed->splineInfo->xmin, splint_generic(variance_reversed, variance_reversed->splineInfo->xmin), variance_reversed->splineInfo->xmax, splint_generic(variance_reversed, variance_reversed->splineInfo->xmax));
		exit(0);
	}
		
	double cdm, rvir, rs;
	//del_nl_z() ????
	cdm = (c0 / (1.0 + params->z)) * pow(halo_mass / Mstar, beta); // eqn. B6
	//if (cdm < 1.0) cdm = 1.0; //otherwise unphysical??
	rvir = pow(3.0*halo_mass/(4.0*M_PI*rho_bar_z(params)*del_nl), 1.0/3.0); // eqn. B5
	rs = rvir/cdm;	

	HOD_Parameters *toReturn = malloc(sizeof(*toReturn));
	toReturn->cdm = cdm;
	toReturn->rvir = rvir;
	toReturn->rs = rs;
	free(variance_reversed);
	return toReturn;
}

HOD_Parameters **HOD_params_catalogue(BinInfo *mass_binInfo, Spline *variance, Parameters *params) {
	double m = 0.0;
	HOD_Parameters **toReturn = malloc(mass_binInfo->bins * sizeof(*toReturn));
	for (int i = 0; i < mass_binInfo->bins; i++) {
		m = bin_to_x(mass_binInfo, i);
		toReturn[i] = HOD_params(m, variance, params);	
	}
	return toReturn;
}

//u(k_m) = fourier transform x*sin(kx)*(ln(1 + x/a) - (x/a)/(1 + x/a))  (+ the 1/( ln(1+c) - 1/(1+c))normalization )
//cumulative sum in rho NFW up to r, given a halo of mass m, r here in units of rs
double cumsum_nfw(double r, double cdm) {
	double cs = (log(1.0 + r) - r/(1.0 + r)) / (log(1.0 + cdm) - cdm/(1.0 + cdm));
	if (cs < 0.0 || cs > 1.000001 || isnan(cs)) {
		printf("cumsum wtf? cumsum: %le for r: %lf, cdm: %lf \n", cs, r, cdm);
		exit(0);
	}
	return cs;		
}

//outputs REVERSED [integral of rho_nfw * dV from 0 to r for multiple r going from 0 to r_virial]
Spline *prep_nfw_cumsums(double halo_mass, int bins, Spline *variance, Parameters *params) {

	double r = 0.0;
	double cumsum = 0.0;		
	double *Rs = malloc((size_t)bins * sizeof(*Rs));
	double *cumsums = malloc((size_t)bins * sizeof(*cumsums));		
	
	HOD_Parameters *HOD_parameters = HOD_params(halo_mass, variance, params);
	double cdm = HOD_parameters->cdm;		
	double rvir = HOD_parameters->rvir;
	double rs = HOD_parameters->rs;	
	
	BinInfo *r_binInfo = prep_bins(0.0, rvir, bins, NORMAL_BIN);
	for (int i = 0; i < bins; i++) {
		r = bin_to_x(r_binInfo, i);		
		//for cumsum r in units of rs
		cumsum = cumsum_nfw(r/rs, cdm);
		Rs[i] = r;
		cumsums[i] = cumsum;
		if (cumsum < 0.0 || cumsum > 1.000001 || isnan(cumsum)) {
			printf("cumsum: %lf in prep_nfw_cumsums(). check variance spline range \n", cumsum);
			exit(0);
		}
	}
		
	
	Spline *toReturn = input_spline_values(bins, cumsums, Rs, MEASURED);	
	free(Rs);
	free(cumsums);
	free(HOD_parameters);
	return toReturn;
}

//prepares nfw cumulative sums for a given catalogue, places into array of *Spline objects
Spline **nfw_cumsums_catalogue(BinInfo *mass_binInfo, int cumsum_bins, Spline *variance, Parameters *params) {
	double m = 0.0;
	Spline **NFW_cumsums = malloc(mass_binInfo->bins * sizeof(*NFW_cumsums));
	for (int i = 0; i < mass_binInfo->bins; i++) {
		m = bin_to_x(mass_binInfo, i);
		NFW_cumsums[i] = prep_nfw_cumsums(m, cumsum_bins, variance, params);	
	}
	return NFW_cumsums;
}


//adds centrals according to HOD paper
Catalogue *HOD_centrals(Catalogue *current_catalogue) {
	double mean_N_central = 0.0;
	double this_halo_mass = 0.0;
	int has_central = 0;
	int total_centrals = 0;

	//find out how many centrals and allocate them
	int *have_centrals = malloc((size_t)current_catalogue->particle_no * sizeof(*have_centrals));
	
	for (int i = 0; i < current_catalogue->particle_no; i++) {
		this_halo_mass = current_catalogue->particles[i]->mass;
		mean_N_central = 0.5 * (1.0 + gsl_sf_erf((log10(this_halo_mass) - log_Mmin)/sigma_logm));		
		has_central = (int) gsl_ran_bernoulli(gsl_ran_r, mean_N_central);
				
		// boolean list over haloes. 
		have_centrals[i] = has_central;
		if (has_central == 1) total_centrals += 1;		
	}	
	
	Particle **centrals = malloc((size_t)total_centrals * sizeof(*centrals));
	printf("Adding %d centrals\n\n", total_centrals);
	
	int counter = 0;	
	for (int i = 0; i < current_catalogue->particle_no; i++) {
		if (have_centrals[i] == 1) {
		
			char *label = malloc(20 * sizeof(*label));
			sprintf(label, "central %d", counter);
			
			Particle *this_central = new_particle(current_catalogue->particles[i]->x, current_catalogue->particles[i]->y, current_catalogue->particles[i]->z, current_catalogue->particles[i]->vx, current_catalogue->particles[i]->vy, current_catalogue->particles[i]->vz, current_catalogue->particles[i]->mass, label, CENTRAL, current_catalogue->particles[i]->parent_id);
					
			centrals[counter] = this_central;	
			counter += 1;			  
		}
	}
	
	Catalogue *toReturn = malloc(sizeof(*toReturn));
	toReturn->particle_no = total_centrals;
	toReturn->particles = centrals;
	
	free(have_centrals);
	return toReturn;
}


Catalogue *HOD_satellites(Catalogue *current_catalogue, Spline *variance, Parameters *params) {

	double mean_N_sats = 0.0;
	double sum_weights = 0.0;
	double sum_weights_sq = 0.0;
	double cos_theta_sat = 0.0;
 	double phi_sat = 0.0;
 	double theta_sat = 0.0;
 	double NFW_cumsum = 0.0;
	double r_sat = 0.0;
	double x_sat = 0.0;
	double y_sat = 0.0;
	double z_sat = 0.0;
	double vx_sat = 0.0;
	double vy_sat = 0.0;
	double vz_sat = 0.0;;
	double sigma_sq = 0.0;
	double sigma = 0.0;
	double c = 0.0;
	double velDisp = 0.0;
	double cdm_vel_par = 0.0;
	double mean_sigma_vel = 0.0;
	double mean_sigma_sq_vel = 0.0;
	
	double *mass_lims = mass_min_max(current_catalogue);
	double Mmin = mass_lims[0];
	double Mmax = mass_lims[1];
	int massBins = 500;
	BinInfo *mass_binInfo = prep_bins(Mmin, Mmax, massBins, LOG_BIN);	
	
	//rho_nfw cumulative sums (inverted)
	Spline **NFW_cumsums = nfw_cumsums_catalogue(mass_binInfo, 300, variance, params);
	//rvir, rs, cdm for each mass bin
	HOD_Parameters **HOD_parameters_catalogue = HOD_params_catalogue(mass_binInfo, variance, params);		
	
	//find how many sats per each halo, also find weighted shotnoise
	int *no_satellites = malloc((size_t)current_catalogue->particle_no * sizeof(*no_satellites));	
	int satellite_no = 0;		
	for (int i = 0; i < current_catalogue->particle_no; i++) {	
		mean_N_sats = pow((current_catalogue->particles[i]->mass - M0)/M1, alpha);
		no_satellites[i] = gsl_ran_poisson(gsl_ran_r, mean_N_sats);		
		sum_weights += (double) no_satellites[i];
		sum_weights_sq += (double) (no_satellites[i] * no_satellites[i]);
		satellite_no += no_satellites[i];
	}	
		
	//shotnoise for each catalogue is weighted by the stochastic number of satellites allocated to the haloes
	double SN = sum_weights_sq / pow(sum_weights, 2.0);
	double weighted_shotnoise_sats = SN*volume;
	
	//globally stored to find average for all mocks
	sat_no_tot += satellite_no;
	weighted_shotnoise_tot += weighted_shotnoise_sats;
	
	//satellite array
	Particle **satellites = malloc((size_t)satellite_no * sizeof(*satellites));	
	
	int satellite_index = 0;
	double mean_r_sat = 0.0;
	double mean_r_sq_sat = 0.0;
	int cur_massBin = 0;	
	for (int i = 0; i < current_catalogue->particle_no; i++) {
		//each central i will have no_satellites[i] as found before	
		cur_massBin = x_to_bin(mass_binInfo, current_catalogue->particles[i]->mass);
		//velocity dispersion, paper1 eqn. 35			
		c = HOD_parameters_catalogue[cur_massBin]->cdm;
		cdm_vel_par = (c*(1.0 - pow(1.0+c, -2.0) - 2.0*log(1.0+c)/(1.0+c))/(2.0*pow(log(1.0+c) - c/(1.0+c), 2.0)));
		sigma_sq = cdm_vel_par*(G*current_catalogue->particles[i]->mass)/(3.0*HOD_parameters_catalogue[cur_massBin]->rvir);		
		sigma = sqrt(sigma_sq);
		mean_sigma_vel += sigma;
		mean_sigma_sq_vel += sigma_sq;
		for (int j = 0; j < no_satellites[i]; j++) {
			//picking values uniformly randomly distributed on a sphere
			cos_theta_sat = randomDouble(-1.0, 1.0);
			theta_sat = acos(cos_theta_sat);
			phi_sat = randomDouble(0.0, 2.0*pi);		

			//inverse transform sampled r
			r_sat = splint_generic(NFW_cumsums[cur_massBin], randomDouble(0.0, 1.0));
				
			//track the mean radial displacement
			mean_r_sat += r_sat;
			mean_r_sq_sat += r_sat*r_sat;				
	
			//actual position is with respect to the i-th central
			x_sat = r_sat * sin(theta_sat) * cos(phi_sat) + current_catalogue->particles[i]->x;
			y_sat = r_sat * sin(theta_sat) * sin(phi_sat) + current_catalogue->particles[i]->y;
			z_sat = r_sat * cos_theta_sat + current_catalogue->particles[i]->z;				
			
			velDisp = gsl_ran_gaussian(gsl_ran_r, sigma);			
			//this is the relative velocity to the halo, so add together with the host halo velocity
			vx_sat = current_catalogue->particles[i]->vx + velDisp;
			vy_sat = current_catalogue->particles[i]->vy + velDisp;
			vz_sat = current_catalogue->particles[i]->vz + velDisp;	
		
			//make new particle
			char *label = malloc(20 * sizeof(*label));
			sprintf(label, "satellite %d", satellite_index);			
			Particle *satellite = new_particle(x_sat, y_sat, z_sat, vx_sat, vy_sat, vz_sat, 0.0, label, SATELLITE, current_catalogue->particles[i]->parent_id);
			
			//PBC, add to satellite catalogue
			PBC(satellite);			
			satellites[satellite_index] = satellite;			
			satellite_index += 1;
		}		
	}
	
	//sanity check
	mean_r_sat /= (double) satellite_no;
	mean_r_sq_sat /= (double) satellite_no;
	mean_sigma_vel /= (double) current_catalogue->particle_no;
	mean_sigma_sq_vel /= (double) current_catalogue->particle_no;
	printf("HOD %d satellites populated. mean r of sats: %lf, sigma r sat: %lf, mean sigma vel: %lf, sigma of sigma vel: %lf\n\n", satellite_no, mean_r_sat, sqrt(mean_r_sq_sat - pow(mean_r_sat, 2.0)), mean_sigma_vel, sqrt(mean_sigma_sq_vel - pow(mean_sigma_vel, 2.0))); 
	
	free(no_satellites);
	free(mass_binInfo);
	for (int i = 0; i < massBins; i++) {
		free_spline(NFW_cumsums[i]);
		free(HOD_parameters_catalogue[i]);
	}
	Catalogue *toReturn = malloc(sizeof(*toReturn));
	toReturn->particle_no = satellite_no;
	toReturn->particles = satellites;
	return toReturn;
}

Catalogue *full_catalogue_HOD(Catalogue *current_catalogue, Parameters *params) {

	double *mass_lims = mass_min_max(current_catalogue);
	double Mmin = mass_lims[0];
	double Mmax = mass_lims[1];
	double R1 = m_to_R(Mmin, params);
	double R2 = m_to_R(Mmax, params);
	int R_bins = 300;
	BinInfo *R_binInfo = prep_bins(R1, R2, R_bins, LOG_BIN);
	
	int z_bins = 300;
	double z_min = 0.0;
	double z_max = 4.0;	

	BinInfo *z_binInfo = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	
	//D+(z), GR only, k-independent
	Spline *Dplus_z = Dplus_spline(GR, 0.0, z_binInfo, params);	
	
	char Pk_out_temp[200];
	sprintf(Pk_out_temp, "%s/data/tmp/Pk_full_cat_HOD_tmp.dat", home_directory);
	haloes_measure_Pk(current_catalogue, Pk_out_temp, 0.0, 1.0, CIC, 1.0);
	Spline *Pk_current = input_spline_file(Pk_out_temp, k_Pkmono_format, MEASURED);
	
	int k_bins = Pk_current->splineInfo->lines;
	double k_min = Pk_current->splineInfo->xmin;
	double k_max = Pk_current->splineInfo->xmax;
	BinInfo *k_binInfo = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	
	Spline *Pk_current_z0 = prep_Pk_constz(GR, Dplus_z, params->z, 0.0, Pk_current, k_binInfo, z_binInfo, params);
	Spline *variance_z0 = prep_variances(R_binInfo, Pk_current_z0);
	
	Catalogue *centrals = HOD_centrals(current_catalogue);			
	Catalogue *satellites = HOD_satellites(current_catalogue, variance_z0, params);
	Catalogue *final_catalogue = combine_catalogues(centrals, satellites);
	
	free_catalogue(centrals);
	free_catalogue(satellites);	
	free_spline(Pk_current_z0);
	free_spline(variance_z0);
	free_spline(Dplus_z);	
	free(k_binInfo);
	free(z_binInfo);
	free(R_binInfo);
	return final_catalogue;
}

int full_catalogue_HOD_file(char DM_catalogue_in[], char HOD_catalogue_out[], char input_format[], Parameters *params) {
	Catalogue *current_catalogue = input_catalogue_file(DM_catalogue_in, 0, input_format);	
	Catalogue *final_catalogue = full_catalogue_HOD(current_catalogue, params);
	catalogue_to_file(final_catalogue, HOD_catalogue_out);	
	free_catalogue(final_catalogue);
	free_catalogue(current_catalogue);
	return 0;
}


