//size limits for box
double volume_limits[3];
double volume_limits_save[3];
double volume;

//jenkins folding
double Jenkins_foldfactor = 1.0;

//grid size, number of cells 
int cells[3];
int cells_displ[3];

double *grid;

//int sat_no_tot = 0;
//double weighted_shotnoise_tot = 0.0;

//data input formats
//x y z vx vy vz mass
char mocks_format[] = "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %le";
//x y z vx vy vz mass parent_halo_id
char mocks_format_withIds[] = "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %le \t %*d";


char galaxies_format[] = "%le \t %le \t %le \t %le \t %le \t %le";
char cube_format[] = "%lf \t %lf \t %lf \t %*f \t %*f \t %lf";
char spline_format_Pk[] = "%lf \t %lf \t %*f \t %*d \n";
char monopoles_rescaling_format[] = "%lf \t %lf \t %*f \n";
char Pk_model_format[] = "%le \t %le \n";
char variance_format[] = "%le \t %le \n";

char pk_format[] = "%lf \t %le \t %le \t %le \t %d \t %lf \t %lf";
char Pk_format_covMat[] = "%*f \t %le \t %le \t %*e \t %*d \t %*f \t %*f";


char k_Pkmono_format[] = "%lf \t %le \t %*e \t %*e \t %*d \t %*f \t %*f";
char pk_format_folding[] = "%lf \t %le \t %le \t %*d \t %*f \t %*f";

//fftw data
fftw_plan overdensity_plan, ZA_displacements_plan;
fftw_complex *overdensity, *overdensity_fourier;
fftw_complex *displacements_fourier, *displacements;

int ph;

//Pk calculation data
int spectrum_size;
double* fk_abs_squares;
double* monopole;
double* quadrupole;
double* hexadecapole;
int* k_number;
bool* toRemove_bins;
double* bin_k_average;
double* bin_k_min;
double* bin_k_max;
//1st element is sum of L2 for that bin, 2nd is sum for [L2 squared], 3rd is sum of [L2 * fourier coefficient abs squared]
double** binsums_L2;
double** binsums_L4;
double WindowFunc_pow;



//catalogue input files
char mock_directory[] = "/home/jonas/Testing_GR/Mocks/L1500_np1024_a0.7_halo";
char* home_directory;
char temp_directory[100];
char cube_file[] = "/home/jonas/Testing_GR/Mocks/cube_gal_-20.0_vel.dat";
//char cube_directory[];

//HOD population stuff
double alpha, log_M1, M0, M1, delta_c, c0;
double beta, del_nl;
double log_Mmin, sigma_logm, log_M0;

//cosmological parameters
double rho_cosm, G;

//rescaling
double lambda_F4, lambda_F5, lambda_F6, lambda_GR;
double fR0_F4, fR0_F5, fR0_F6, fR0_GR;
//double R1_primed, R2_primed;
double dsq_minimized_value, s, z_rescaled;
int Dplus_bins;

//for ZA
double sigma_exp_smoothed;

int prep_oneHalo(void);
int oneHalo(void);

//functions.c
double randomDouble(double min, double max);
int PBC(Particle *p);
int countLines(FILE *f);
int arr_ind(int x, int y, int z);
int arr_ind_displ(int x, int y, int z);
double* mass_min_max(Catalogue *catalogue);
double m_to_R(double mass, Parameters *params);
double R_to_m(double R, Parameters *params);


//memory.c
Particle **initParticles(int particle_no);
int init_rng(void);
int clearParticles(Particle** particles);
int clearGrid(void);
int overdensity_clearMemory(void);
int overdensity_allocateMemory(void);
int overdensity_freeMemory(void);
int ZA_clearMemory(void);
int ZA_freeMemory(void);
int ZA_allocateMemory(void);
int clearPk(void);
int clearMemory(void);
int memCheck(void);
int init_spectrum_storage(int bins);
void calloc2D_double(double*** arr, int dim1, int dim2);
void free2D_double(double*** arr, int dim2);

//power_spectrum.c
int populateParticles(Particle** particles, int particle_no);
int toRedshift(Catalogue *cat, Parameters *params);
int initGrid(void);
int populateGrid(Catalogue *catalogue, int grid_func);
int populateGridNGP(Catalogue *catalogue);
int populateGridCIC(Catalogue *catalogue);
int gridIntoOverdensity(int particle_no);
int fftw_overdensity_calculate(void);
int overdensity_pk();
int pk_logbin(int overdensity_index, double k, double mu, BinInfo* Pk_binInfo);
int multipole_calc(void);
int pk_to_file_logbin(const char filename[]);
int delsq_to_file_logbin(const char filename[]);


//halo_model.c
int scale_velocities(double s, Catalogue *cat, Parameters *cur_params, Parameters *targ_params);
int scale_masses(double s, Catalogue *current_catalogue, Parameters *cur_params, Parameters *targ_params);
Spline **nfw_cumsums_catalogue(BinInfo *mass_binInfo, int cumsum_bins, Spline *variance_reversed, Parameters *params);
double halo_model_Pk(double k, double shot_noise, double twohalo_prefactor, int redshift_space, int order, Spline *Pk_spline, double beta, double R_nl);
Spline* prep_nfw_cumsums(double halo_mass, int bins, Spline *variance_reversed, Parameters *params);
int haloModel_realSpace(char out[], BinInfo *k_bins, double shot_noise, double twohalo_prefactor, Spline *Pk_spline, double scaling_factor, double R_nl, HOD_Parameters *HOD_params);
int haloModel_out(char out[], BinInfo *k_bins, double shot_noise, double twohalo_prefactor, int redshift_space, Spline *Pk_spline, Parameters *par, double scaling_factor, double mass_bias, double R_nl, HOD_Parameters *HOD_params);
double cumsum_nfw(double r, double cdm);
int set_halo_params(double halo_mass, Spline *variance_spline, Parameters *params);
int set_ZA_params(void);
double ukm_nfw_profile(double k, double rs, double cdm);
double I_integral_HOD(double t, void *params);
HOD_Parameters *HOD_params(double halo_mass, Spline *variance_reversed, Parameters *params);

//catalogue_input.c


//HOD.c
int add_centrals(void);
int prep_HOD_Pk(bool sats);
int add_satellites(void);
double weighted_shotnoise(void);
int prep_bin_info(void);
int HOD_main(void);

//rescaling.c


int halo_mass_func(void);
int rescale(void);
int rescale_testRun(void);
int rescale_catalogue(void);
int rescale_model(void);

//rescaling_functions.c
double delta_sq_int_func(double R, void *params);
double OLV(double k, void *params);
double integrate_OLV(double R, Pk_Spline *pk_spline);
double delta_sq_rms(const gsl_vector *inputs, void *params);
int dsq_multimin(bool vary_z_cur, double z_init, Spline* variance_const, Spline** variances_varz, BinInfo *z_binInfo, double *R_lims);
int dsq_multimin_test(void);
//Spline* prep_variances_test(BinInfo *R_bin_info, Spline *Pk_spline, double s_test);
double delsq_lin(Pk_Spline *pk_spline, double k);
double expected_variance_smoothed(double R_nl_this, Pk_Spline *pk_spline);
double OLV_smoothed(double k, void *params);
int scale_box(double s);
int move_particles(Catalogue *current_catalogue, double s);

//Spline* prep_Pk_constz(int gravity, Spline* Dplus_spline, double z_from, double z_to, Spline *Pk_spline_in, BinInfo *k_bin_info, BinInfo *z_binInfo, Parameters *params);
//Spline *prep_Pk_constz_rescaled(int gravity, Spline *Dplus_inp, double z_from, double z_to, Spline *Pk_spline_in, BinInfo *k_bin_info, BinInfo *z_binInfo, Parameters *params, double s);

Spline* Dplus_spline(int gravity, double k, BinInfo *z_binInfo, Parameters *params);
//Spline* prep_variances(BinInfo *R_bin_info, Pk_Spline *pk_spline);
int generate_sigma_plots(void);


//Dplus.c
double z_to_t(double z, Parameters *params);
double z_to_a(double z);
double a_to_z(double a);
double a_to_t(double a, Parameters *params);
double H_a(double a, Parameters *params);
double z_to_t_int(double z, void * params);
Spline* prep_redshift_spline(BinInfo *z_binInfo, Parameters *params);
double t_to_z(double t, Spline *redshift_reverse);
double H_z(double z);
int growth_function(double t, const double y[], double f[], void *params);
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
int Dplus_calc(int gravity, double k, double* zs, double* Dpluses, BinInfo *z_binInfo, Parameters *params);

//bias
double m_to_v(double mass, Parameters *params, Spline* variance_spline);
double f_v(double v);
double ln_fv_deriv(double v);
double v_to_m(double v, Bias_Params *bias_params);
double b_v(double v);
double b_m(double m, Parameters *params, Spline* variance_spline);
double bf_over_m(double v, void * params);
double f_over_m(double v, void * params);
int bias_plot(double z);
double calc_b_eff(Parameters *params, Spline* variance_spline, double Mmin, double Mmax);

//fold.c
Catalogue *Jenkins_fold_volume(Catalogue *catalogue);
int Jenkins_restore(Catalogue *catalogue, Catalogue *catalogue_saved);
FoldingInformation* combine_folded(const char unfolded_path[], const char folded_path[], char pk_format[], char outFile[], double kmin, double kmax, double k_fold, int* foldPoints, bool print);

//covariance_matrices
int generate_covariance(char inPath[], char outPath[], int sim_number, int startIndex);

//runs.c
int pre_covariance_run(const char inputData[], const char input_format[], const char outFile[], bool save_files, bool fold, int skipLines);
int randoms_fullRun(int particle_number, const char outFile[], bool fold, bool test_mode);
int haloes_measure_Pk(Catalogue *catalogue, char outFile[], double shotnoise, double foldfactor, int grid_func, double multip_fact);
	
//tests.c
int test_run_randoms(void);
int full_test(bool verbose);
int memory_test(bool verbose);
int overdensity_pk_test(bool verbose);
int randoms_populate_test(bool verbose);
int test_file_io(bool verbose);
int test_random(bool verbose);


//#JustGSLThings
const gsl_rng_type* gsl_ran_T;
struct gls_rng* gsl_ran_r;
gsl_integration_workspace* z_to_t_workspace;
size_t z_to_t_workspace_size = 1000;
gsl_integration_workspace* dsq_workspace; 
gsl_integration_workspace* OLV_workspace;
size_t OLV_workspace_size = 2000;
gsl_integration_workspace* onehalo_workspace;
gsl_integration_cquad_workspace* dsq_workspace_cquad;






