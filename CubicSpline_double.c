void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]){
// Given arrays x[0..n-1] and y[0..n-1] containing a tabulated function, i.e y_i = f(x_i)
// with x_1 < x_2 < .. < x_n and given values yp1 and ypn for the first derivative
// of the interpolating function at points 0 and (n-1), respectively, this
// routine returns an array y2[0..n-1] that contains the second derivatives of the 
// interpolating function at the tabulated points x_i.  If yp1 and/or ypn are 
// equal to 10^30 or larger, the routine is signaled to set the corresponding
// boundary conditions for a natural spline, with zero second derivative on that
// boundary. 

    int i, k;
    
    double p, qn, sig, un;
    
    // u = vector(1,n);                  
    double u[n];

    if(yp1 >0.99e30)           // The lower boundary condition is set either
        y2[0] = u[0] =0.0;      // to be "natural" or else to have a specified first
                                // derivative.
    else{
        y2[0]=-0.5;
         u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0]) - yp1);
    }
    
    for(i=1; i<n-1; i++){                   // This is the decomposition loop of the tridiagonal    
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);  // algorithm. y2 and u are used for temporary storage 
        
        p=sig*y2[i-1]+2.0;                  // of the decomposed factors. 
        
        y2[i]=(sig-1.0)/p;
        
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) -(y[i]-y[i-1])/(x[i]-x[i-1]);
        
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }

    if(ypn >0.99e30)           // The upper boundary condition is set either to be "natural"
        qn = un = 0.0;              // or to have a specified first derivative. 
    else {
        qn=0.5;
        un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }

    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2] + 1.0);

    for(k=n-2; k>=0; k--){      // This is the backsubstitution loop of the tridiagonal algorithm.
        y2[k] = y2[k]*y2[k+1]+u[k];
    }
}


void splint(double xa[], double ya[], double y2a[], int n, double x, double *y){
// Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the x_a's in order)
// and given the array y2a[1..n], which is the output from spline above, and given a value of 
// x, this routine returns a cubic spline interpolated value y. 

  int klo, khi, k;
  double h, b, aa;

  klo=0;
  khi=n-1;

  while (khi-klo >1){              // We will find the right place in the table by means of 
    k=(khi+klo) >>1;               // bisection. This is optimal if sequential calls to this routine
    if (xa[k] > x) khi=k;          // are at random values of x. If sequential calls are in order, 
    else klo=k;                    // and closely spaced, one would do better to store previous values
  }                                // of klo and khi and test if they remain appropriate on the next call.
                                   // klo and khi now bracket the input value of x.

  h=xa[khi]-xa[klo];
  if (equalsDouble(h, 0.0)) printf("Bad xa input to routine splint");  // The xa's must be distinct. 
  aa = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;  // Cubic spline polynomial is now evaluated.
  *y= aa*ya[klo]+b*ya[khi]+((aa*aa*aa-aa)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

int powerlaw_regression(int N, double rmin, double rmax, double sign, double ri[], double Di[], double* norm, double* index){
    // W    = A r^(n + 3.).
    // ln W = A' + n*log(r).
    
    // A'   = ln(A)
    
    double detA           = 0.0;
    
    // parameters.
    double Adash          = 0.0;
    double A              = 0.0;
    double n              = 0.0;
    
    double sum_lnri       = 0.0;
    double sum_lnri2      = 0.0;
    double sum_lnDi       = 0.0;
    double sum_lnDi_lnri  = 0.0;  
    
    double b1             = 0.0;
    double b2             = 0.0;
    
    int goodN             =   0;

    for(int i=0; i<N; i++){                
        if((ri[i]>= rmin) && (ri[i] <= rmax) && (sign*Di[i] > 0.0)){
            goodN           += 1;
            
            // printf("\n%e \t %e", ri[i], sign*Di[i]);
            
            sum_lnri        += log(ri[i]);
            sum_lnri2       += log(ri[i])*log(ri[i]);
        
            sum_lnDi        += log(sign*Di[i]);
            sum_lnDi_lnri   += log(sign*Di[i])*log(ri[i]);
        }
    }
    
    // printf("\n\n%d good N", goodN);
    
    // For a matrix A, paramater vector (A', n)^T, and vector B
    // Required to invert AP  = B. 2x2 matrix inversion.
         	
    // A reads as (a b) for a = N, b = sum log(r_i), c = sum log(r_i), d = sum log(r_i)*log(r_i).
    //            (c d)        

    // and B    = (b1, b2)^T for b1 = sum[log(Di) - 3log(ri)], b2 = sum[ln(Di)*ln(ri) - 3ln(ri)^2]
     
    // det = ad - bc.
    detA   = goodN*sum_lnri2 - sum_lnri*sum_lnri;
      
    // (A', n)^T = (A^-1)B  = (1./detA)*( d*b1 - b*b2) 
    //                                  (-c*b1 + a*b2)
    
    b1    = sum_lnDi      - 3.*sum_lnri;
    b2    = sum_lnDi_lnri - 3.*sum_lnri2;
     
    Adash = pow(detA, -1.)*(sum_lnri2*b1 - sum_lnri*b2); 
    
    A     = exp(Adash);
    
    n     = pow(detA, -1.)*(-sum_lnri*b1 + goodN*b2); 
    
    // printf("\n\nPower law regression: %e \t %e", sign*A, n);
    
    *norm = sign*A;

    *index = n;

    return 0;
}
