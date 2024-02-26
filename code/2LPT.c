//copied from 2LPTic


/* ??????????????????????????????????????????????/
 hubble_a = Hubble * sqrt(Omega / pow(InitTime, 3) + (1 - Omega - OmegaLambda) / pow(InitTime, 2) + OmegaLambda);

  vel_prefac = InitTime * hubble_a * F_Omega(InitTime);
  vel_prefac2 = InitTime * hubble_a * F2_Omega(InitTime);

  vel_prefac /= sqrt(InitTime);	// converts to Gadget velocity 
  vel_prefac2 /= sqrt(InitTime);
  
  */

//fg used for ZA velocities
double F_Omega(double a)
{
  double omega_a;

  omega_a = Omega / (Omega + a * (1 - Omega - OmegaLambda) + a * a * a * OmegaLambda);

  return pow(omega_a, 0.6);
}

//fg used for 2LPT velocities
double F2_Omega(double a)
{
  double omega_a;

  omega_a = Omega / (Omega + a * (1 - Omega - OmegaLambda) + a * a * a * OmegaLambda);

  return 2 * pow(omega_a, 4./7.);
}
