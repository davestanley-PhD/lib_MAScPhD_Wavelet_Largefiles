#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  if((nrhs != 9))
    {
      printf("Error: Invalid input parameters\n");
      printf("Usage L = lmax(x,DIMM,tau,dt,bb,cc,mult,angmx1,fiduc1)\n");
      printf("DIMM,tau,mult and fiduc1 are integer parameters\n");
      return;
    }

  mxArray *data,*param;
  double *x,*out,*p;
  double L,bb,cc,angmx1, params[9];
  int DIMM,tau,fiduc1,mult;
  float dt;
  int n;
  int i,j;
  int nRow,nCol;

  data = prhs[0];
  x = mxGetPr(data);
  nCol = mxGetN(data);
  nRow = mxGetM(data);

  if(nRow != 1)
    {
      printf("Error: Provide a row vector for x\n");
      return;
    }
  n = nCol;


  for(i = 1; i < 9; i++)
    {
      param = prhs[i];
      nCol = mxGetN(param);
      nRow = mxGetM(param);
      if((nCol == 1) && (nRow == 1))
	{
	  p = mxGetPr(param);
	  params[i] = p[0];
	}
      else
	{
	  printf("Error: Invalid parameter number %d\n",i+1);
	  printf("Error: Use only scalar parameters\n");
	  return;
	}
    }

  DIMM = (int)params[1];
  tau = (int)params[2];
  dt = (float)params[3];
  bb = params[4];
  cc = params[5];
  mult = (int)params[6];
  angmx1 = params[7];
  fiduc1 = (int)params[8];


  lamd_(x,&n,&L,&DIMM,&tau,&dt,&bb,&cc,&mult,&angmx1,&fiduc1);

  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  out = mxGetPr(plhs[0]);
  out[0] = L;

  return;
}
