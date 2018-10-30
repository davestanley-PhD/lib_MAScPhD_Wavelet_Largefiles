#include "mex.h"

double *x,*out,*p;
double L,bb,cc,angmx1, params[9];
int DIMM,tau,fiduc1,mult;
float dt;
int n;
int i,j;
int nRow,nCol;

//void lamd_(double*,int*,double*,int*,int*,float*,double*,double*,int*,double*,int*);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //printf("STLMAX\n");
  //printf("nlhs=%d, nrhs=%d\n",nlhs,nrhs);
  if(/*(nlhs != 1) ||*/ (nrhs != 9))
    {
      printf("Error: Invalid input parameters\n");
      printf("Usage L = lmax(x,DIMM,tau,dt,bb,cc,mult,angmx1,fiduc1)\n");
      printf("DIMM,tau,mult and fiduc1 are integer parameters\n");
      return;
    }

  //data = prhs[0];
  x = mxGetPr(prhs[0]);
  nCol = mxGetN(prhs[0]);
  nRow = mxGetM(prhs[0]);
  //printf("nRow=%d, nCol=%d\n",nRow,nCol);
  if(nRow != 1)
    {
      printf("Error: Provide a row vector for x\n");
      return;
    }
  n = nCol;

  //Load remaining parameters
  for(i = 1; i < 9; i++)
    {
      //param = prhs[i];
      nCol = mxGetN(prhs[i]);
      nRow = mxGetM(prhs[i]);
      if((nCol == 1) && (nRow == 1))
	{
	  p = mxGetPr(prhs[i]);
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
  //printf("DIMM=%d,tau=%d,dt=%f,bb=%f,cc=%f\n",DIMM,tau,dt,bb,cc);
  //printf("mult=%d,angmx1=%f,fiduc1=%d\n",mult,angmx1,fiduc1);

  /*
    lamd_(double*,uunsigned int*,double*,int*,int*,float*,double*,double*,int*,double*,int*)
    SUBROUTINE LMAD(x,isize,zlyap,DIMM_,tau_,dt_,bb_,cc_,mult_,angmx1_,fiduc1_)
  */
  lamd_(x,&n,&L,&DIMM,&tau,&dt,&bb,&cc,&mult,&angmx1,&fiduc1);
  //lamd_(x,&n,&L);

  //return data to matlab
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  out = mxGetPr(plhs[0]);
  out[0] = L; //Setup STLMAX output value

  return;
}
