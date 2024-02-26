bool equalsDouble(double x, double y) {
	if (fabs(x - y) < DBL_EPSILON) {
		return true;
	} else {
		return false;
	}
}

bool equalsDoubleArray(double *arr1, double *arr2, int len) {
	for (int i = 0; i < len; i++) {
		if (!equalsDouble(arr1[i], arr2[i])) {
			printf("mismatch at index %d. arr1[%d] = %le, arr2[%d] = %le\n", i, i, arr1[i], i, arr2[i]);
			return false;
		}
	}
	return true;
}

//void exit_err()
