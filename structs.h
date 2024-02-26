
typedef struct {
	double xmin;
	double xmax;
	int bins;
	double dx;
	int bin_type;
} BinInfo;

typedef struct {
	double *x_vals;
	double *y_vals;
	double *coeffs;
	int lines;
	double xmin;
	double xmax;
} Spline;

typedef struct {
	Spline *spline;
	double pk_loA;
	double pk_hiA;
	double pk_lon;
	double pk_hin;
	int model;
	double z;
} Pk_Spline;

typedef struct {
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double mass;
	char *label;
	int type;
	int parent_id;
} Particle;

typedef struct {
	Particle **particles;
	int particle_no;
} Catalogue;


typedef struct {	
	int uf_first_ind;
	int uf_last_ind;
	int f_first_ind;
	int f_last_ind;
	bool max_reached;
	bool success;	
} FoldingInformation;

typedef struct {
	int folds;
	int** foldPoints;
} Foldings;


typedef struct {
	double H0;
	double omega_m_0;
	double omega_v_0;// = 1.0 - omega_m;
	double omega_r_0;
	double omega;
	double gamma; //f = omega_m ^ gamma
	double z;	
} Parameters;


typedef struct {
	Parameters *parameters;
	Spline *variance;
	Spline *variance_reversed;
} Bias_Params;

typedef struct {
	Pk_Spline *pk_spline;
	double R;	
} OLV_parameters;

typedef struct {
	bool vary_z_current;
	double R1_primed;
	double R2_primed;
	Spline *variance_const;
	Spline **variances_varz;
	BinInfo *z_binInfo;
} Multimin_Params;

typedef struct {	
	double s;
	double z_var;
	bool vary_z_current;
	Spline *variance_const;
	Spline **variances_varz;
	BinInfo *z_binInfo;
} Dsq_Params;



typedef struct {
	double cdm;
	double rvir;
	double rs;	
} HOD_Parameters;
