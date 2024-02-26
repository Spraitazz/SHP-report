/*
Takes the input file path string, which should contain %d for iteration,
the output path string, which should also contain %d, the number of simulations
and the index (usually 0) of the first simulation.
INPUT FILE MUST HAVE THE FORMAT "%lf \t %lf \t %lf \t %d \n" indicating k, monopole,
quadrupole and bin number, respectively.
Generates covariance matrix and outputs it to the outPath, where the k-values
increase down and to the right in the file.
*/
/*
int generate_covariance(char inPath[], char outPath[], int sim_number, int startIndex) {

	//check indices
	if (sim_number < 1 || startIndex < 0) {
		printf("wrong indices\n");
		return 1;
	}	
	
	FILE *inFile, *outFile;
	char in_filename[100];
	char out_filename[100];
	int lineno, binCount;		
	double sum_mono, sum_quadro;	
	double temp[3] = {0.0, 0.0, 0.0};
	
	sprintf(in_filename, inPath, startIndex);
	inFile = fopen(in_filename, "r");
	if (inFile == NULL) {
			printf("invalid covariance matrix input file %s\n", in_filename);
			exit(0);		
	}
	binCount = countLines(inFile);
	fclose(inFile);
	
	//Format: multipoles[sim no][bin no][mono/quadro]
	double multipoles[sim_number][binCount][2];	
	
	//collect multipoles to array
	for (int i = startIndex; i < sim_number + startIndex; i++) {
		sprintf(in_filename, inPath, i);
		inFile = fopen(in_filename, "r");
		//file exists check
		if (inFile == NULL) {
			printf("invalid covariance matrix input file %s\n", in_filename);
			exit(0);		
		}
		
		lineno = countLines(inFile);

		//all files should have the same number of lines (can also check k-vals are equal)
		if (lineno != binCount) {
			printf("Input files have different lengths, file: %s\n", in_filename);
			exit(0);
		}
		
		for (int k = 0; k < lineno; k++) {
			ph = fscanf(inFile, Pk_format_covMat, &multipoles[i - startIndex][k][0], 
			&multipoles[i - startIndex][k][1]);	
			//multipoles[i - startIndex][k][1] = fabs(multipoles[i - startIndex][k][1]);	
		}		
		fclose(inFile);
	}	
	
	//variables that depended on binCount
	double averages_multipoles[binCount][2];
	//Matr(i,j) Format: [mono-mono, mono-quadro, quadro-quadro][i][j]
	double cov_matrix[3][binCount][binCount];
	double cov_matrix_final[3][binCount][binCount];
	
	//averages of multipoles and quadrupoles for all sims
	for (int bin = 0; bin < binCount; bin++) {
		sum_mono = 0.0;
		sum_quadro = 0.0;
		for (int sim = 0; sim < sim_number; sim++) {
			sum_mono += multipoles[sim][bin][0];
			sum_quadro += multipoles[sim][bin][1];
		}
		averages_multipoles[bin][0] = sum_mono / (double)sim_number;
		averages_multipoles[bin][1] = sum_quadro / (double)sim_number;
	}		

	//calculate covariance entries
	for (int i = 0; i < binCount; i++) {
		for (int j = 0; j < binCount; j++) {				
			temp[0] = 0.0;	
			temp[1] = 0.0;
			temp[2] = 0.0;
			for (int sim = 0; sim < sim_number; sim++) {
				//P0 P0				
				temp[0] += ((multipoles[sim][i][0] - averages_multipoles[i][0]) * (multipoles[sim][j][0] - averages_multipoles[j][0]));
				//P0 P2
			    temp[1] += ((multipoles[sim][i][0] - averages_multipoles[i][0]) * (multipoles[sim][j][1] - averages_multipoles[j][1]));
			    //P2 P2
			    temp[2]	+= ((multipoles[sim][i][1] - averages_multipoles[i][1]) * (multipoles[sim][j][1] - averages_multipoles[j][1]));				   	
			}		
			cov_matrix[0][i][j] = temp[0] / (double) sim_number;
			cov_matrix[1][i][j] = temp[1] / (double) sim_number;
			cov_matrix[2][i][j] = temp[2] / (double) sim_number;
		}	
	}		

	//normalize by STDEVS
	double sigma_j, sigma_k;
	for (int poleOrder = 0; poleOrder < 3; poleOrder++) {
		for (int j = 0; j < binCount; j++) {
		
			sigma_j = sqrt(fabs(cov_matrix[poleOrder][j][j]));
			
			for (int k = 0; k < binCount; k++) {	
						
				sigma_k = sqrt(fabs(cov_matrix[poleOrder][k][k]));	
							
				if (sigma_j == 0.0 || sigma_k == 0.0) {
					printf("zero variance error in cov matrices, [%d][%d][%d]\n", poleOrder, j, k);
					return 1;
				} else {
					cov_matrix_final[poleOrder][j][k] = cov_matrix[poleOrder][j][k] / (sigma_j * sigma_k);
				}
			}
		}	
	}	
	
	//print to files, 0 - mono/mono, 1 - mono/quadro, 2 - quadro/quadro
	for (int i = 0; i < 3; i++) {
	
		sprintf(out_filename, outPath, i);
		outFile = fopen(out_filename, "w");
				
		for (int j = 0; j < binCount; j++) {
			for (int k = 0; k < binCount; k++) {
				fprintf(outFile, "%lf \t", cov_matrix_final[i][j][k]);
			}
			fprintf(outFile, "\n");
		}	
		fclose(outFile);
	}
	return 0;
}
*/

int generate_covariance_oneFile(char inPath[], char outPath[], int sim_number, int startIndex, char Pk_format[]) {

	//check indices
	if (sim_number < 1 || startIndex < 0) {
		printf("wrong indices\n");
		return 1;
	}	
	
	FILE *inFile, *outFile;
	char in_filename[200];
	int lineno = 0;
	int binCount = 0;		
	double sum_mono = 0.0;
	double sum_quadro = 0.0;	
	double temp = 0.0;
	int sc = 0;
	
	sprintf(in_filename, inPath, startIndex);
	inFile = fopen(in_filename, "r");
	if (inFile == NULL) {
			printf("invalid covariance matrix input file %s\n", in_filename);
			exit(0);		
	}
	binCount = countLines(inFile);
	fclose(inFile);
	
	//Format: multipoles[sim no][bin no]
	double multipoles[sim_number][binCount*2];	
	
	//collect multipoles to array
	for (int i = startIndex; i < sim_number + startIndex; i++) {
		sprintf(in_filename, inPath, i);
		inFile = fopen(in_filename, "r");
		//file exists check
		if (inFile == NULL) {
			printf("invalid covariance matrix input file %s\n", in_filename);
			exit(0);		
		}
		
		lineno = countLines(inFile);

		//all files should have the same number of lines (can also check k-vals are equal)
		if (lineno != binCount) {
			printf("Input files have different lengths, file: %s\n", in_filename);
			exit(0);
		}
		
		for (int k = 0; k < lineno; k++) {
			sc = fscanf(inFile, Pk_format, &multipoles[i - startIndex][k], 
			&multipoles[i - startIndex][binCount + k]);	
			sc = fscanf(inFile, "\n");
		}	
			
		fclose(inFile);
	}
	
	
	//variables that depended on binCount
	double averages_multipoles[binCount*2];
	//Matr(i,j) Format: [mono-mono, mono-quadro, quadro-quadro][i][j]
	double cov_matrix[binCount*2][binCount*2];
	double cov_matrix_final[binCount*2][binCount*2];
	
	//averages of multipoles and quadrupoles for all sims
	for (int bin = 0; bin < binCount*2; bin++) {		
		for (int sim = 0; sim < sim_number; sim++) {
			averages_multipoles[bin] += multipoles[sim][bin] / (double)sim_number;
		}		
	}		

	//calculate covariance entries
	for (int i = 0; i < binCount*2; i++) {
		for (int j = 0; j < binCount*2; j++) {				
			temp = 0.0;				
			for (int sim = 0; sim < sim_number; sim++) {
				temp += ((multipoles[sim][i]-averages_multipoles[i])*(multipoles[sim][j]-averages_multipoles[j]));	
			}		
			cov_matrix[i][j] = temp / (double) sim_number;			
		}	
	}		

	//normalize by STDEVS
	double sigma_j, sigma_k;

	for (int j = 0; j < binCount*2; j++) {
	
		sigma_j = sqrt(fabs(cov_matrix[j][j]));
		
		for (int k = 0; k < binCount*2; k++) {	
					
			sigma_k = sqrt(fabs(cov_matrix[k][k]));	
						
			if (sigma_j < 1e-10 || sigma_k < 1e-10) {
				printf("zero variance error in cov matrices, [%d][%d]. sigma %d: %lf, sigma %d: %lf\n", j, k, j, sigma_j, k, sigma_k);
				return 1;
			} else {
				cov_matrix_final[j][k] = cov_matrix[j][k] / (sigma_j * sigma_k);
			}
		}
	}	
	
	
	//print to file
	outFile = fopen(outPath, "w");
			
	for (int j = 0; j < binCount*2; j++) {
		for (int k = 0; k < binCount*2; k++) {
			fprintf(outFile, "%lf \t", cov_matrix_final[j][k]);
		}
		fprintf(outFile, "\n");
	}	
	fclose(outFile);
	
	return 0;
}
/*
int covariance_sequence() {
	char inPath[100];
	sprintf(inPath, "/home/jonas/Testing_GR/Mocks/W1/mock_%s_zlim_0.6_0.9_Jf_1.dat", "%d");
	char outPath[100];
	char outDel1[100];
	char outDel2[100];
	
	//sprintf(outPath, "%s/data/cov/305_mocks_cov_%s.dat", home_directory, "%d");
	//generate_covariance(inPath, outPath, 305, 1);
	
	sprintf(outPath, "%s/data/cov/305_mocks_cov.dat", home_directory);
	generate_covariance_oneFile(inPath, outPath, 305, 1);
	
	//for (int i = 305; i < 306; i++) {
		//sprintf(outPath, "%s/data/cov/cov_%d_sim%s.dat", home_directory, i, "%d");
		//printf("generating covariance, in path: %s, out path: %s \n", inPath, outPath);
		//generate_covariance(inPath, outPath, i, 1);
		//sprintf(outDel1, "%s", outPath);
		//sprintf(outDel1, outDel1, 1);
		//sprintf(outDel2, "%s", outPath);
		//sprintf(outDel2, outDel2, 2);
		//remove(outDel1);
		//remove(outDel2);
	//}
	
	
	return 0;

}
*/

int covariance_catalogue(char cats_in[], int sim_number, int startIndex, Parameters *params, char outPath[]) {
	
	char mockPath[200];
	char Pk_temp_out[200];
	
	sprintf(Pk_temp_out, "%s/data/Pk_300_mocks/Pk_cov_tmp_%s.dat", home_directory, "%d");
	catalogues_measure_Pk(cats_in, Pk_temp_out, mocks_format_withIds, sim_number, startIndex, 1.0, CIC, 1.0, params);
	generate_covariance_oneFile(Pk_temp_out, outPath, sim_number, startIndex, Pk_format_covMat);
	
	for (int i = 0; i < sim_number; i++) {
		//sprintf(Pk_temp_out, "%s/data/tmp/Pk_cov_tmp_%d.dat", home_directory, i+startIndex);
		//remove(Pk_temp_out);
	}	

	return 0;
}
