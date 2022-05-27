// ******************************************************************
// **                                                              **
// ** PUMA 3.0.1 (June 10th, 2009)                                 **
// **                                                              **
// ** Modified on June 3rd, 2011                                   **  
// **                                                              **
// ** Modification: Crystalline quartz substrate refractive index  **
// ** refractive (SUBSTRATE=30) was corrected and Amorphous Quartz **
// ** substrate (SUBSTRATE=60) was added.                          **
// **                                                              **
// ** The problem of estimating the thickness and the optical      **
// ** constants of thin films using transmission data only is very **
// ** challenging from the mathematical point of view, and has a   **
// ** technological and an economic importance. In many cases it   **
// ** represents a very ill-conditioned inverse problem with many  ** 
// ** local-nonglobal solutions. In a recent publication, we       **
// ** proposed nonlinear programming models for solving this       **
// ** problem. Well-known software for linearly constrained        **
// ** optimization was used with success for this purpose.         **
// **                                                              **
// ** The software PUMA implements an unconstrained formulation of **
// ** the nonlinear programming model, which solves the estimation **
// ** problem using a method based on repeated calls to a recently **
// ** introduced unconstrained minimization algorithm. Numerical   **
// ** experiments on computer-generated films show that the new    **
// ** procedure is reliable.                                       **
// **                                                              **
// ** References:                                                  **
// **                                                              **
// ** E.G. Birgin, I. Chambouleyron, J.M. Martinez,                **
// ** "Estimation of the optical constants and the thickness of    **
// ** thin films using unconstrained optimization", Journal of     **
// ** Computational Physics 151, pp. 862-880, 1999.                **
// **                                                              **
// ** M. Mulato, I. Chambouleyron, E.G. Birgin and J.M. Martinez,  **
// ** "Determination of thickness and optical constants of a-Si:H  **
// ** films from transmittance data", Applied Physics Letters 77,  **
// ** pp. 2133-2135, 2000.                                         **
// **                                                              **
// ** I. Chambouleyron, S. D. Ventura, E. G. Birgin, and           **
// ** J. M. Martinez, "Optical constants and thickness             **
// ** determination of very thin amorphous semiconductor films",   **
// ** Journal of Applied Physics 92, pp. 3093-3102, 2002.          **
// **                                                              **
// ** E. G. Birgin, I. Chambouleyron, and J. M. Martinez,          **
// ** "Optimization problems in the estimation or parameters of    **
// ** thin films and the elimination of the influence of the       **
// ** substrate", Journal of Computational and Applied             **
// ** Mathematics 152, pp. 35-50, 2003.                            **
// **                                                              **
// ** E. G. Birgin, I. Chambouleyron, J. M. Martinez, and          **
// ** S. D. Ventura, "Estimation of optical parameters of very     **
// ** thin films", to appear in Applied Numerical Mathematics 47,  **
// ** pp. 109-119, 2003.                                           **
// **                                                              **
// ** S. Ventura, E. G. Birgin, J. M. Martínez and                 **
// ** I. Chambouleyron, "Optimization techniques for the           **
// ** estimation of the thickness and the optical parameters of    **
// ** thin films using reflectance data", Journal of Applied       **
// ** Physics 97, 043512, 2005.                                    **
// **                                                              **
// ** R. Andrade, E. G. Birgin, I. Chambouleyron, J. M. Martínez   **
// ** and S. D. Ventura, "Estimation of the thickness and the      **
// ** optical parameters of several superimposed thin films using  **
// ** optimization", Applied Optics 47, pp. 5208-5220, 2008.       **
// **                                                              **
// ******************************************************************

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define ninobsmax                10000
#define nobsmax                    200
#define nlayersmax                   5
#define nfilmsmax     (nlayersmax - 3)

#define gamma                  1.0e-04
#define lmin                   1.0e-16
#define lmax                   1.0e+01

#define nmax (2 * nobsmax * nfilmsmax)
#define mmax                       100

#define TRUE                         1
#define FALSE                        0
#define PI      3.14159265358979311600                          

#define min(x,y) ((x)<(y)) ? (x) : (y)
#define max(x,y) ((x)>(y)) ? (x) : (y)

char   firstc   (char s[]);
int    Error    (char * buffer);
int    Iniv     (int n, double v[], double firste[]);
int    Incrv    (int n, double v[], double firste[], double laste[], double estep[]);
void   compn    (int nlayers, int slayer, int nobs, double z[][nobsmax], double n[][nobsmax]);
void   compdn   (int nlayers, int slayer, int nobs, double z[][nobsmax], 
		 double n[][nobsmax], double dndz[][nobsmax]);
void   compk    (int nlayers, int slayer, int nobs, double ip[], double lambda[], 
		 double v[][nobsmax], double k[][nobsmax]);
void   compdk   (int nlayers, int slayer, int nobs, double ip[], double lambda[], 
		 double v[][nobsmax], double k[][nobsmax], double dkdv[][nobsmax]);
void   compt    (int nobs, double n[], double k[], double d, double s[], 
                 double lambda[], double t[]);
void   compdt   (int nobs, double n[], double k[], double d, double s[], 
                 double lambda[], double t[], double dtdn[], double dtdk[]);
void   compr    (int nobs, double n[], double k[], double d, double s[], 
                 double lambda[], double r[]);
void   compdr   (int nobs, double n[], double k[], double d, double s[], 
                 double lambda[], double r[], double drdn[], double drdk[]);
void   comptr   (int nlayers, int slayer, int nobs, double n[][nobsmax], double k[][nobsmax], 
		 double d[], double lambda[], double t[], double r[]);
void   compdtr  (int nlayers, int slayer, int nobs, double n[][nobsmax], double k[][nobsmax], 
                 double d[], double lambda[], double t[], double r[],
		 double dtdn[][nobsmax],double dtdk[][nobsmax], double drdn[][nobsmax], 
                 double drdk[][nobsmax]);
double compqe   (int nobs, double tobs[], double t[]);
double compdqe  (int nobs, double tobs[], double t[], double dqedt[]);
double evalf    (int ydim, double y[]);
double evalg    (int ydim, double y[], double g[]);
double sg       (int n, double x[], int m, double eps, int maxit,
		 int maxfc);
void   IniFilm  (char id[],int nlayers,int slayer,int substrate,char datatype,
                 int nobs,double lini,double lfin, char datfile[]);
void   IniPoint (double n0[],double nf[],double k0[],double y[]);
void   GenReport(double bestqe, double bestd[], double bestip[], 
		 char solfile[], char inffile[]);

struct tfilm
{
  char   id     [50];

  int    nlayers;

  int    slayer;
  int    substrate;
  double s      [nobsmax];

  char   datatype;

  int    nobs;

  double lambda [nobsmax];
  double h;

  double tobs   [nobsmax];
  double robs   [nobsmax];

  double d      [nfilmsmax];
  double ip     [nfilmsmax];
} film;

//*******************************************************************
//*******************************************************************

int main(int argc, char **argv)
{
  int substrate,nlayers,slayer,nfilms,nobs,maxit,iniptype,i,lip,ld,
      ln0,lnf,lk0,n,m,maxfc,skip;
  double lini,lfin,qe,bestqe,eps;
  char datatype;
  FILE * fp;

  char id[50],solfile[50],inffile[50],datfile[50],message[255];
  double y[2*nfilmsmax*nobsmax],yini[2*nfilmsmax*nobsmax],bestd[nfilmsmax],
      bestip[nfilmsmax],dini[nfilmsmax],dfin[nfilmsmax],dstep[nfilmsmax],
      ipini[nfilmsmax],ipfin[nfilmsmax],ipstep[nfilmsmax],n0ini[nfilmsmax],
      n0fin[nfilmsmax],n0step[nfilmsmax],nfini[nfilmsmax],nffin[nfilmsmax],
      nfstep[nfilmsmax],k0ini[nfilmsmax],k0fin[nfilmsmax],k0step[nfilmsmax],
      n0[nfilmsmax],nf[nfilmsmax],k0[nfilmsmax];

  if ( argc < 8 )
    Error("Incorrect number of parameters (see documentation)");

  strcpy (id,argv[1]);

  nlayers   = atoi  ( argv[2] );
  slayer    = atoi  ( argv[3] );
  substrate = atoi  ( argv[4] );
  datatype  = firstc( argv[5] );
  nobs      = atoi  ( argv[6] );
  lini      = atof  ( argv[7] );
  lfin      = atof  ( argv[8] );

  nfilms = nlayers - 3;

  if ( argc != 12 + 15 * nfilms && argc != 12 + 6 * nfilms )
    Error("Incorrect number of parameters (see documentation)");

  strcpy(solfile,id);
  strcpy(inffile,id);
  strcpy(datfile,id);
  strcat(datfile,"-dat.txt");
  strcat(solfile,"-sol.txt");
  strcat(inffile,"-inf.txt");

  IniFilm(id,nlayers,slayer,substrate,datatype,nobs,lini,lfin,datfile);

  maxit    = atoi( argv[9] );
  bestqe   = atof( argv[10] );
  iniptype = atoi( argv[11] );

  if ( iniptype == 0 )
    {
    if ( argc != 12 + 15 * nfilms )
      Error("Incorrect number of parameters (see documentation)");

    for ( i = 0; i < nfilms; i ++ )
      {
 	dini[i]   = atof( argv[12+15*i] ); 
 	dfin[i]   = atof( argv[13+15*i] );
	dstep[i]  = atof( argv[14+15*i] );

	ipini[i]  = atof( argv[15+15*i] );
	ipfin[i]  = atof( argv[16+15*i] );
	ipstep[i] = atof( argv[17+15*i] );

	n0ini[i]  = atof( argv[18+15*i] );
	n0fin[i]  = atof( argv[19+15*i] );
	n0step[i] = atof( argv[20+15*i] );
    
	nfini[i]  = atof( argv[21+15*i] );
	nffin[i]  = atof( argv[22+15*i] );
	nfstep[i] = atof( argv[23+15*i] );

	k0ini[i]  = atof( argv[24+15*i] );
	k0fin[i]  = atof( argv[25+15*i] );
	k0step[i] = atof( argv[26+15*i] );
      }
    }
  else
    {
    if ( argc != 12 + 6 * nfilms )
      Error("Incorrect number of parameters (see documentation)");

    for ( i = 0; i < nfilms; i ++ )
      {
 	dini[i]   = atof( argv[12+6*i] ); 
 	dfin[i]   = atof( argv[13+6*i] );
	dstep[i]  = atof( argv[14+6*i] );

	ipini[i]  = atof( argv[15+6*i] );
	ipfin[i]  = atof( argv[16+6*i] );
	ipstep[i] = atof( argv[17+6*i] );
      }
    }

  // Set SG parameters

  n      = 2 * nobs * nfilms;
  maxfc  = 10 * maxit;
  m      = 60;
  eps    = 1.0e-10;
  
  // Solve the optimization problems from different initial trials

  switch(iniptype)
    {
    case 0:

      for ( lip = Iniv(nfilms,film.ip,ipini); !lip; 
            lip = Incrv(nfilms,film.ip,ipini,ipfin,ipstep) )
      for ( ld  = Iniv(nfilms,film.d,dini); !ld; 
            ld  = Incrv(nfilms,film.d,dini,dfin,dstep) )
      for ( ln0 = Iniv(nfilms,n0,n0ini); !ln0; 
            ln0 = Incrv(nfilms,n0,n0ini,n0fin,n0step) ) 
      for ( lnf = Iniv(nfilms,nf,nfini); !lnf; 
            lnf = Incrv(nfilms,nf,nfini,nffin,nfstep) )
      for ( lk0 = Iniv(nfilms,k0,k0ini); !lk0; 
            lk0 = Incrv(nfilms,k0,k0ini,k0fin,k0step) )
        {
	  skip = FALSE;
	  for ( i = 0; i < nfilms; i++ )
	    if ( n0[i] < nf[i] )
 	      skip = TRUE;

          if ( ! skip )
	    { 
              printf("\n");
 	      for ( i = 0; i < nfilms; i ++ )
	        printf("%.2f %.2f %.2f %.2f %.2f ",
                film.ip[i],film.d[i],n0[i],nf[i],k0[i]);

              IniPoint(n0,nf,k0,y);

              qe = sg(n,y,m,eps,maxit,maxfc);

	      printf("%e ",qe);

              if ( qe < bestqe )
	        {
	          bestqe = qe;

                  for ( i = 0; i < nfilms; i ++ )
		    {
		      bestd[i]  = film.d[i];
		      bestip[i] = film.ip[i];
		    }

	          printf("!!!");

                  fp = fopen(solfile,"wt");
	          for ( i = 0; i < n; i ++ )
		    fprintf(fp,"%.20e\n",y[i]);
	          fclose(fp);
	        }
	    }
        }

        break;

    case 9:

      fp = fopen(solfile,"rt");

      if ( fp == NULL )
	{
	  strcat(strcat(strcpy(message,"File "),solfile)," does not exists");
	  Error(message);
	}

      for ( i = 0; i < n; i ++ )
	fscanf(fp,"%lf", & yini[i]);

      fclose(fp);

      for ( lip = Iniv(nfilms,film.ip,ipini); !lip; 
            lip = Incrv(nfilms,film.ip,ipini,ipfin,ipstep) )
      for ( ld  = Iniv(nfilms,film.d,dini); !ld; 
            ld = Incrv(nfilms,film.d,dini,dfin,dstep) )
        {
          printf("\n");
          for ( i = 0; i < nfilms; i ++ )
	    printf("%.2f ",film.ip[i]);
	  for ( i = 0; i < nfilms; i ++ )
	    printf("%.2f ",film.d[i]);

	  for ( i = 0; i < n; i ++ )
	    y[i] = yini[i];

          qe = sg(n,y,m,eps,maxit,maxfc);

	  printf("%e ",qe);

          if ( qe < bestqe )
	    {
	      bestqe = qe;

              for ( i = 0; i < nfilms; i ++ )
		{
		  bestd[i]  = film.d[i];
		  bestip[i] = film.ip[i];
		}

	      printf("!!!");

              fp = fopen(solfile,"wt");
	      for ( i = 0; i < n; i ++ )
		fprintf(fp,"%.20e\n",y[i]);
	      fclose(fp);
	    }
        }

        break;

    default:
      Error("Invalid initial point type (choose 0 or 9)");
    }

  GenReport(bestqe,bestd,bestip,solfile,inffile);

  return 0;
}

//*******************************************************************
//*******************************************************************

double sg(int n,double x[],int m,double eps,int maxit,int maxfc)
{
  double alpha,atemp,f,fbest,fmax,fnew,ginfn,gtd,lambda,sts,sty;
  int i,iter,fcnt;

  double g[nmax],gnew[nmax],s[nmax],y[nmax],d[nmax],xbest[nmax],
    xnew[nmax],lastfv[mmax];

  iter = 0;
  fcnt = 0;

  for ( i = 0; i < m; i ++ )
    lastfv[i] = -1.0e+99;

  for ( i = 0; i < n; i ++ )
    xbest[i] = x[i];

  f = evalg(n,x,g);
  fcnt ++;

  lastfv[0] = f;
  fbest     = f;

  ginfn = 0;
  for ( i = 0; i < n; i ++ )
    ginfn = max( ginfn, fabs( g[i] ) );

  if ( ginfn != 0 )
    lambda =  min( lmax, max( lmin, 1 / ginfn ) );

  while ( ginfn > eps && iter < maxit && fcnt < maxfc )
    {
      iter ++;

      gtd = 0;
      for ( i = 0; i < n; i ++ )
	{
          d[i] = - lambda * g[i];
          gtd = gtd + g[i] * d[i];
	}

      fmax = lastfv[0];
      for ( i = 1; i < m; i ++ )
	fmax = max( fmax, lastfv[i] );

      alpha = 1;

      for ( i = 0; i < n; i ++ )
	xnew[i] = x[i] + d[i];

      fnew  = evalf(n,xnew);
      fcnt ++;

      while ( fnew > fmax + gamma * alpha * gtd && fcnt < maxfc )
	{
	  if ( alpha <= 0.1 ) 
	    alpha = alpha / 2;

          else
	    {
	      atemp = ( -gtd * alpha * alpha ) / 
		( 2.0 * ( fnew - f - alpha * gtd) );
              if ( atemp < 0.1 || atemp > 0.9 * alpha ) 
		atemp = alpha / 2;
	      alpha = atemp;
	    }

	  for ( i = 0; i < n; i ++ )
	    xnew[i] = x[i] + alpha * d[i];

	  fnew = evalf(n,xnew);
	  fcnt ++;
	}      

      if ( fcnt >= maxfc )
	break;
               
      f = fnew;
      lastfv[iter % m] = f;

      if ( f < fbest ) 
	{
          fbest = f;
          for ( i = 0; i < n; i ++ )
	    xbest[i] = xnew[i];
	}

      evalg(n,xnew,gnew);

      sts = 0;
      sty = 0;
      for ( i = 0; i < n; i ++ )
	{
          s[i] = xnew[i] - x[i];
          y[i] = gnew[i] - g[i];
          sts  = sts + s[i] * s[i];
          sty  = sty + s[i] * y[i];
          x[i] = xnew[i];
          g[i] = gnew[i];
	}

      ginfn = 0;
      for ( i = 0; i < n; i ++ )
	ginfn = max( ginfn, fabs( g[i] ) );

      if ( sty <= 0 ) 
	lambda = lmax;

      else
	lambda = min( lmax, max( lmin, sts / sty ) );

      //printf("\n%d %f %f",iter,f,ginfn);
    }

  f = fbest;

  for ( i = 0; i < n; i ++ )
    x[i] = xbest[i];

  return f;
}

//*******************************************************************
//*******************************************************************

void compn(int nlayers,int slayer,int nobs,
	   double z[][nobsmax],double n[][nobsmax])
{
  int i,j;
  double h;

  h = film.h;

  // Substrate and air refractive indices
  for ( j = 0; j < nobs; j ++ )
    {
      n[0][j] = 1;
      n[nlayers-1][j] = 1;
      n[slayer][j] = film.s[j];
    }

  // Films refractive indices
  for ( i = 1; i <= slayer-1; i ++ )
    {
      n[i][nobs-1] = 1 + z[i-1][nobs-1] * z[i-1][nobs-1];
      n[i][nobs-2] = n[i][nobs-1] + h * z[i-1][nobs-2] * z[i-1][nobs-2];

      for ( j = nobs-3; j >= 0; j -- )
	n[i][j] = 2 * n[i][j+1] - n[i][j+2] + h * h * z[i-1][j] * z[i-1][j];
    }

  for ( i = slayer; i <= nlayers-3; i ++ )
    {
      n[i+1][nobs-1] = 1 + z[i-1][nobs-1] * z[i-1][nobs-1];
      n[i+1][nobs-2] = n[i+1][nobs-1] + h * z[i-1][nobs-2] * z[i-1][nobs-2];

      for ( j = nobs-3; j >= 0; j -- )
	n[i+1][j] = 2 * n[i+1][j+1] - n[i+1][j+2] + h * h * z[i-1][j] * z[i-1][j];
    }
}

//*******************************************************************
//*******************************************************************

void compdn(int nlayers,int slayer,int nobs,
	    double z[][nobsmax],double n[][nobsmax], double dndz[][nobsmax])
{
  int i,j;
  double h;

  h = film.h;

  // Substrate and air refractive indices
  for ( j = 0; j < nobs; j ++ )
    {
      n[0][j] = 1;
      n[nlayers-1][j] = 1;
      n[slayer][j] = film.s[j];
    }

  // Films refractive indices
  for ( i = 1; i <= slayer-1; i ++ )
    {
      n[i][nobs-1] = 1 + z[i-1][nobs-1] * z[i-1][nobs-1];
      n[i][nobs-2] = n[i][nobs-1] + h * z[i-1][nobs-2] * z[i-1][nobs-2];

      dndz[i-1][nobs-1] = 2 * z[i-1][nobs-1];
      dndz[i-1][nobs-2] = 2 * h * z[i-1][nobs-2];

      for ( j = nobs-3; j >= 0; j -- )
	{
	  n[i][j] = 2 * n[i][j+1] - n[i][j+2] + h * h * z[i-1][j] * z[i-1][j];
	  dndz[i-1][j] = 2 * h * h * z[i-1][j];
	}
    }

  // Films refractive indices
  for ( i = slayer; i <= nlayers-3; i ++ )
    {
      n[i+1][nobs-1] = 1 + z[i-1][nobs-1] * z[i-1][nobs-1];
      n[i+1][nobs-2] = n[i+1][nobs-1] + h * z[i-1][nobs-2] * z[i-1][nobs-2];

      dndz[i-1][nobs-1] = 2 * z[i-1][nobs-1];
      dndz[i-1][nobs-2] = 2 * h * z[i-1][nobs-2];

      for ( j = nobs-3; j >= 0; j -- )
	{
	  n[i+1][j] = 2 * n[i+1][j+1] - n[i+1][j+2] + h * h * z[i-1][j] * z[i-1][j];
	  dndz[i-1][j] = 2 * h * h * z[i-1][j];
	}
    }

}

//*******************************************************************
//*******************************************************************

void compk(int nlayers,int slayer,int nobs,
	   double ip[],double lambda[],double v[][nobsmax],double k[][nobsmax])
{
  int i,ipind,j,jmax;
  double h;

  h = film.h;

  // Substrate and air absorption coefficient
  for ( j = 0; j < nobs; j ++ )
    {
      k[0][j] = 0;
      k[nlayers-1][j] = 0;
      k[slayer][j] = 0;
    }

  // Films absorption coefficient
  for ( i = 1; i <= slayer-1; i ++ )
    {
      ipind = 0;
      while ( ipind < nobs && lambda[ipind] < ip[i-1] )
	ipind ++;

      jmax = min( nobs-3, ipind-1 );

      k[i][nobs-1] = v[i-1][nobs-1] * v[i-1][nobs-1];
      k[i][nobs-2] = k[i][nobs-1] + h * v[i-1][nobs-2] * v[i-1][nobs-2];

      for ( j = nobs-3; j >= ipind; j -- )
	k[i][j] = 2 * k[i][j+1] - k[i][j+2] + h * h * v[i-1][j] * v[i-1][j];

      for ( j = jmax; j >= 0; j -- )
	k[i][j] = 2 * k[i][j+1] - k[i][j+2] - h * h * v[i-1][j] * v[i-1][j];
    }

  // Films absorption coefficient
  for ( i = slayer; i <= nlayers-3; i ++ )
    {
      ipind = 0;
      while ( ipind < nobs && lambda[ipind] < ip[i-1] )
	ipind ++;

      jmax = min( nobs-3, ipind-1 );

      k[i+1][nobs-1] = v[i-1][nobs-1] * v[i-1][nobs-1];
      k[i+1][nobs-2] = k[i+1][nobs-1] + h * v[i-1][nobs-2] * v[i-1][nobs-2];

      for ( j = nobs-3; j >= ipind; j -- )
	k[i+1][j] = 2 * k[i+1][j+1] - k[i+1][j+2] + h * h * v[i-1][j] * v[i-1][j];

      for ( j = jmax; j >= 0; j -- )
	k[i+1][j] = 2 * k[i+1][j+1] - k[i+1][j+2] - h * h * v[i-1][j] * v[i-1][j];
    }
}

//*******************************************************************
//*******************************************************************

void compdk(int nlayers,int slayer,int nobs,double ip[],double lambda[], 
double v[][nobsmax],double k[][nobsmax],double dkdv[][nobsmax])
{
  int i,ipind,j,jmax;
  double h;

  h = film.h;

  // Substrate and air absorption coefficient
  for ( j = 0; j < nobs; j ++ )
    {
      k[0][j] = 0;
      k[nlayers-1][j] = 0;
      k[slayer][j] = 0;
    }

  // Films absorption coefficient
  for ( i = 1; i <= slayer-1; i ++ )
    {
      ipind = 0;
      while ( ipind < nobs && lambda[ipind] < ip[i-1] )
	ipind ++;

      jmax = min( nobs-3, ipind-1 );

      k[i][nobs-1] = v[i-1][nobs-1] * v[i-1][nobs-1];
      k[i][nobs-2] = k[i][nobs-1] + h * v[i-1][nobs-2] * v[i-1][nobs-2];

      dkdv[i-1][nobs-1] = 2 * v[i-1][nobs-1];
      dkdv[i-1][nobs-2] = 2 * h * v[i-1][nobs-2];

      for ( j = nobs-3; j >= ipind; j -- )
	{
	  k[i][j] = 2 * k[i][j+1] - k[i][j+2] + h * h * v[i-1][j] * v[i-1][j];
	  dkdv[i-1][j] = 2 * h * h * v[i-1][j];
	}

      for ( j = jmax; j >= 0; j -- )
	{
	  k[i][j] = 2 * k[i][j+1] - k[i][j+2] - h * h * v[i-1][j] * v[i-1][j];
	  dkdv[i-1][j] = - 2 * h * h * v[i-1][j];
	}
    }

  // Films absorption coefficient
  for ( i = slayer; i <= nlayers-3; i ++ )
    {
      ipind = 0;
      while ( ipind < nobs && lambda[ipind] < ip[i-1] )
	ipind ++;

      jmax = min( nobs-3, ipind-1 );

      k[i+1][nobs-1] = v[i-1][nobs-1] * v[i-1][nobs-1];
      k[i+1][nobs-2] = k[i+1][nobs-1] + h * v[i-1][nobs-2] * v[i-1][nobs-2];

      dkdv[i-1][nobs-1] = 2 * v[i-1][nobs-1];
      dkdv[i-1][nobs-2] = 2 * h * v[i-1][nobs-2];

      for ( j = nobs-3; j >= ipind; j -- )
	{
	  k[i+1][j] = 2 * k[i+1][j+1] - k[i+1][j+2] + h * h * v[i-1][j] * v[i-1][j];
	  dkdv[i-1][j] = 2 * h * h * v[i-1][j];
	}

      for ( j = jmax; j >= 0; j -- )
	{
	  k[i+1][j] = 2 * k[i+1][j+1] - k[i+1][j+2] - h * h * v[i-1][j] * v[i-1][j];
	  dkdv[i-1][j] = - 2 * h * h * v[i-1][j];
	}
    }
}

//*******************************************************************
//*******************************************************************
void compt(int nobs, double n[], double k[], double d, double s[], double lambda[], 
double t[])
{
  int i;

  double n2,k2,s2,n2pk2,nps2,nms2,np1,nm1,s2p1,constl,phi,
    s16,s32,A,B1,B2,B,C1,C2,C3,C4,C5,C6,C7,C8,C,D1,D2,D,x,
    x2,q,xdq,Axdq;

  for ( i = 0; i < nobs; i ++ )
    {
      n2     = n[i] * n[i];
      k2     = k[i] * k[i];
      s2     = s[i] * s[i];
      n2pk2  = n2 + k2;
      nps2   = n[i] + s2;
      nms2   = n[i] - s2;
      np1    = n[i] + 1;
      nm1    = n[i] - 1;
      s2p1   = s2 + 1;
      constl = 4 * PI * d / lambda[i];
      phi    = constl * n[i];
      s16    = 16 * s[i];
      s32    = 2 * s16;
             
      A      = s16 * ( n2 + k2 );
             
      B1     = np1 * np1 + k2;
      B2     = np1 * ( n[i] + s2 ) + k2;
      B      = B1 * B2;
             
      C1     = n2pk2 - 1;
      C2     = n2pk2 - s2;
      C3     = 2 * k2;
      C4     = s2p1;
      C5     = 2 * cos( phi );
      C6     = 2 * sin( phi );
      C7     = C1 * C2 - C3 * C4;
      C8     = 2 * C2 + C4 * C1;
      C      = C5 * C7 - C6 * k[i] * C8;
             
      D1     = nm1 * nm1 + k2;
      D2     = nm1 * nms2 + k2;
      D      = D1 * D2;
             
      x      = exp( - constl * k[i] );
      x2     = x * x;
      q      = B - C * x + D * x2;
      xdq    = x / q;
      Axdq   = A * xdq;
  
      t[i]   = Axdq;
    }
}

//*******************************************************************
//*******************************************************************
void compdt(int nobs, double n[], double k[], double d, double s[], 
double lambda[], double t[], double dtdn[], double dtdk[])
{
  int i;

  double n2,k2,s2,n2pk2,nps2,nms2,np1,nm1,s2p1,constl,phi,s16,s32,
    A,B1,B2,B,C1,C2,C3,C4,C5,C6,C7,C8,C,D1,D2,D,x,x2,q,xdq,Axdq,
    dtdA,Axdq2,dtdB,dtdC,dtdD,dtdx,dAdn,dAdk,dB1dn,dB2dn,dBdn,kpk,
    dB1dk,dB2dk,dBdk,npn,dC1dn,dC2dn,dC5dn,dC6dn,dC7dn,dC8dn,dCdn,
    dC1dk,dC2dk,dC3dk,dC7dk,dC8dk,dCdk,dD1dn,dD2dn,dDdn,dD1dk,
    dD2dk,dDdk,dxdk;

 for ( i = 0; i < nobs; i ++ )
    {
      dtdn[i] = 0;
      dtdk[i] = 0;
    }

 for ( i = 0; i < nobs; i ++ )
    {
      n2     = n[i] * n[i];
      k2     = k[i] * k[i];
      s2     = s[i] * s[i];
      n2pk2  = n2 + k2;
      nps2   = n[i] + s2;
      nms2   = n[i] - s2;
      np1    = n[i] + 1;
      nm1    = n[i] - 1;
      s2p1   = s2 + 1;
      constl = 4 * PI * d / lambda[i];
      phi    = constl * n[i];
      s16    = 16 * s[i];
      s32    = 2 * s16;
             
      A      = s16 * ( n2 + k2 );
             
      B1     = np1 * np1 + k2;
      B2     = np1 * ( n[i] + s2 ) + k2;
      B      = B1 * B2;
             
      C1     = n2pk2 - 1;
      C2     = n2pk2 - s2;
      C3     = 2 * k2;
      C4     = s2p1;
      C5     = 2 * cos( phi );
      C6     = 2 * sin( phi );
      C7     = C1 * C2 - C3 * C4;
      C8     = 2 * C2 + C4 * C1;
      C      = C5 * C7 - C6 * k[i] * C8;
             
      D1     = nm1 * nm1 + k2;
      D2     = nm1 * nms2 + k2;
      D      = D1 * D2;
             
      x      = exp( - constl * k[i] );
      x2     = x * x;
      q      = B - C * x + D * x2;
      xdq    = x / q;
      Axdq   = A * xdq;
  
      t[i]   = Axdq;

      dtdA = xdq;
	
      Axdq2 = Axdq / q;
	
      dtdB = - Axdq2;
	
      dtdC = x * Axdq2;
	
      dtdD = - x2 * Axdq2;
	
      dtdx = A / q - Axdq2 * ( 2 * D * x - C );
	
      dAdn = s32 * n[i];
      dtdn[i] += dtdA * dAdn;

      dAdk = s32 * k[i];
      dtdk[i] += dtdA * dAdk;

      dB1dn = 2 * np1;
      dB2dn = nps2 + np1;
      dBdn  = dB1dn * B2 + B1 * dB2dn;
      dtdn[i] += dtdB * dBdn;

      kpk = 2 * k[i];
	
      dB1dk = kpk;
      dB2dk = kpk;
      dBdk  = dB1dk * B2 + B1 * dB2dk;
      dtdk[i] += dtdB * dBdk;

      npn = 2*n[i];
	
      dC1dn = npn;
      dC2dn = npn;
      dC5dn = - C6 * constl;
      dC6dn = C5 * constl;
      dC7dn = dC1dn * C2 + C1 * dC2dn;
      dC8dn = 2 * dC2dn + C4 * dC1dn;
      dCdn  = dC5dn * C7 + C5 * dC7dn - k[i] * ( dC6dn * C8 + C6 * dC8dn );
      dtdn[i] += dtdC * dCdn;

      dC1dk = kpk;
      dC2dk = kpk;
      dC3dk = 2 * kpk;
      dC7dk = dC1dk * C2 + C1 * dC2dk - C4 * dC3dk;
      dC8dk = 2 * dC2dk + C4 * dC1dk;
      dCdk  = C5 * dC7dk - C6 * ( C8 + k[i] * dC8dk );
      dtdk[i] += dtdC * dCdk;

      dD1dn = 2 * nm1;
      dD2dn = nms2 + nm1;
      dDdn  = dD1dn * D2 + D1 * dD2dn;
      dtdn[i] += dtdD * dDdn;

      dD1dk = kpk;
      dD2dk = kpk;
      dDdk  = dD1dk * D2 + D1 * dD2dk;
      dtdk[i] += dtdD * dDdk;

      dxdk = - constl * x;
      dtdk[i] += dtdx * dxdk;
    }
}

//*******************************************************************
//*******************************************************************
void compr(int nobs, double n[], double k[], double d, double s[], double lambda[], 
double r[])
{
  int i;

  double up,down,alpha,x,x2,phi,n2,k2,s2,sphi,cphi,
    c1,c2,c3,c4,n2m1,n2ms2,k2mn2,k2pn2,k2pn22,k2p2n2,
    k2p2n2p1,k2p2n2m1,nps2,nms2,t0,t2,ts,tc,r0,r2,rs,rc,ps,
    pc,p1,p,q0,q2,qs,qc,q;

  for ( i = 0; i < nobs; i ++ )
    {
      alpha= 4 * PI * k[i] / lambda[i];
      x= exp( -alpha * d );
      x2= x*x;
      phi= 4 * PI * n[i] * d / lambda[i];

      // COMMON
      
      n2= n[i]*n[i]; 
      k2= k[i]*k[i]; 
      s2= s[i]*s[i];
      sphi= sin( phi ); 
      cphi= cos( phi );
      n2m1=     -1     +   n2;
      n2ms2=    n2     -   s2;
      k2mn2=    k2     -   n2;
      k2pn2=    k2     +   n2;
      k2pn22=   k2pn2 * k2pn2;
      k2p2n2=   k2     + 2*n2;
      k2p2n2p1= k2p2n2 +    1;
      k2p2n2m1= k2p2n2 -    1;
      nps2= (n[i] + s[i]) * (n[i] + s[i]);
      nms2= (n[i] - s[i]) * (n[i] - s[i]);
      c1= k2 + ( 1 + n[i]) * ( 1 + n[i]);
      c2= k2 + (-1 + n[i]) * (-1 + n[i]);
      c3= k2 + ( 1 + n[i]) * (n[i] + s2); 
      c4= k2 + (-1 + n[i]) * (n[i] - s2);
      
      // DOWN
      
      t0= c1 * c3; 
      t2= c2 * c4;
      tc= n2m1 * n2ms2 + k2 * ( k2p2n2 - 3 * (1 + s2) );
      ts= k[i] * ( -1 + (3 + s2) * k2pn2 - 3*s2 );
      down= 1 / ( t0 - 2*x* ( tc*cphi - ts*sphi ) + t2*x2 );
      
      // UP
      
      r0= c2 * c3; 
      r2= c1 * c4;
      rc= n2m1 * n2ms2 + k2 * (k2p2n2p1 + s2);
      rs= k[i] * (1 + k2pn2) * (-1 + s2);
      up= r0 - 2*x* ( rc*cphi - rs*sphi ) + r2*x2;
      
      // p: real{ B_{11} B_{22} }
      
      ps= sphi * ( 2 * k[i] * n[i] * (k2pn22 - s2) );
      pc= cphi * ( ( 1 + k2mn2 )   
          * k2pn22 + ( n2 * n2m1     + k2 * k2p2n2p1 )*s2 );
      p1= ( 1 - k2mn2 + 4*s[i]) 
          * k2pn22 + ( n2 * (n2 + 1) + k2 * k2p2n2m1 )*s2;
      p= (pc - ps) + 2 * p1 * x + (pc + ps) * x2;
      
      // q: | B_{22} |^2
      
      q0= c1 * ( k2 + nps2 );
      q2= c2 * ( k2 + nms2 );
      qs= 2 * k[i] * (k2pn2 - s[i]) * (1 + s[i]);
      qc= n2m1 * n2ms2 + k2 * ( k2p2n2m1 - s[i]*(4 + s[i]) );
      q= q0 - 2*x* ( qc*cphi - qs*sphi ) + q2*x2;
      
      // R
      
      r[i]= down * (up + 8 * x * (s[i]-1)*(s[i]-1) * p / q);
    }
}

//*******************************************************************
//*******************************************************************
void compdr(int nobs, double n[], double k[], double d, double s[], 
double lambda[], double r[], double drdn[], double drdk[])
{
  int i;

  double np1,np1E2,nm1,nps,nms,nm1E2,sp1,s2p1,s2p3,s2m1,k2p1,
    n2p1,n2m1,n2ms2,k2mn2,k2pn2,k2pn2E2,k2p2n2,k2p2n2p1,
    k2p2n2m1,npsE2,nmsE2,sm1E2,z2n,z2k,z4n,z4s,z2n2,z2k2,c1,
    c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,sp,cp,n2,k2,s2,
    x2,u0,u2,uc,us,u1,u,v0,v2,vc,vs,v1,w,q0,q2,qc,qs,
    q1,q,ps,pc,p1,p2,p,y,y0,pq,yn,yk,beta,x,phi,qq,uu,
    u0n,u2n,ucn,usn,un,u0k,u2k,uck,usk,uk,v0n,v2n,
    vcn,vsn,vn,v0k,v2k,vck,vsk,vk,q0n,q2n,qcn,qsn,
    qn,q0k,q2k,qck,qsk,qk,psn,pcn,p1n,pn,psk,pck,
    p1k,pk;

  for ( i = 0; i < nobs; i ++ )
    {
      beta= 4 * PI * d / lambda[i];
      x= exp( -beta * k[i] );
      phi= beta * n[i];
      
      // COMMON
      
      n2= n[i]*n[i]; k2= k[i]*k[i]; s2= s[i]*s[i]; x2= x*x;
      sp= sin( phi ); cp= cos( phi );
      
      np1=      n[i]   +    1;
      np1E2=    np1    *  np1;
      nm1=      n[i]   -    1;
      nps=      n[i]   +    s[i];
      nms=      n[i]   -    s[i];
      nm1E2=    nm1    *  nm1;
      sp1=      s[i]   +    1;
      s2p1=     s2     +    1;
      s2p3=     s2     +    3;
      s2m1=     s2     -    1;
      k2p1=     k2     +    1;
      n2p1=      1     +   n2;
      n2m1=     -1     +   n2;
      n2ms2=    n2     -   s2;
      k2mn2=    k2     -   n2;
      k2pn2=    k2     +   n2;
      k2pn2E2=  k2pn2 * k2pn2;
      k2p2n2=   k2     + 2*n2;
      k2p2n2p1= k2p2n2 +    1;
      k2p2n2m1= k2p2n2 -    1;
      
      npsE2= nps * nps;
      nmsE2= nms * nms;
      sm1E2= ( s[i] - 1 ) * ( s[i] - 1 );
      
      z2n=  2 * n[i];
      z2k=  2 * k[i];
      z4n=  4 * n[i];
      z4s=  4 * s[i];
      z2n2= 2 * n2;
      z2k2= 2 * k2;
      
      c1=  k2 + np1E2;        
      c2=  k2 + nm1E2;
      c3=  k2 + np1 * (n[i] + s2);
      c4=  k2 + nm1 * (n[i] - s2);
      c5=  1 - 2*k2pn2 + s2;
      c6=  s2p1 + 2*s2;
      c7=  3*s2 + 1;
      c8=  z2k2 + s2p1;
      c9=  n2 + s2p1;
      c10= n2 - s2p1;
      c11= s2p1 - n2/2;
      c12= 2*k2 + k2pn2;
      
      // U
      
      u0= c1 * c3; 
      u2= c2 * c4;
      uc= n2m1 * n2ms2 + k2 * ( k2p2n2 - 3 * s2p1 );
      us= k[i] * ( -1 + s2p3 * k2pn2 - 3*s2 );
      u1= uc*cp - us*sp;
      u= u0 - 2 * x * u1 + u2 * x2;
      
      // UP
      
      v0= c2 * c3; 
      v2= c1 * c4;
      vc= n2m1 * n2ms2 + k2 * (k2p2n2p1 + s2);
      vs= k[i] * ( 1 + k2pn2 ) * s2m1;
      v1= vc*cp - vs*sp;
      w= v0 - 2 * x * v1 + v2 * x2;
      
      // q: | B_{22} |^2
      
      q0= c1 * ( k2 + npsE2 ); 
      q2= c2 * ( k2 + nmsE2 );
      qs= z2k * ( k2pn2 - s[i] ) * sp1;
      qc= n2m1 * n2ms2 + k2 * ( k2p2n2m1 - s2 - z4s );
      q1= qc*cp - qs*sp;
      q= q0 - 2 * x * q1 + q2 * x2;
      
      // p: real{ B_{11} B_{22} }
      
      ps=  z2k * n[i] * (k2pn2E2 - s2);
      pc= ( 1 + k2mn2 )      * k2pn2E2 + s2 * ( n2 * n2m1 + k2 * k2p2n2p1 );
      p1= ( 1 - k2mn2 + z4s) * k2pn2E2 + s2 * ( n2 * n2p1 + k2 * k2p2n2m1 );
      p2= pc*cp + ps*sp;
      p= ( pc*cp - ps*sp ) + 2 * p1 * x + p2 * x2;
      
      // R
      
      pq= p / q; y0= 8 * x * sm1E2; y= y0 * pq;
      r[i]= ( w + y ) / u;

      // derivatives of u
      
      u0n= np1E2 * ( z4n + c6 ) + k2 * ( z4n + ( 2 + s2p1 ) );
      u2n= nm1E2 * ( z4n - c6 ) + k2 * ( z4n - ( 2 + s2p1 ) );
      usk= -c7 + c12 * s2p3;
      u0k=  z2k * ( 2*k2 + np1 * ( z2n + s2p1 ) );
      u2k=  z2k * ( 2*k2 + nm1 * ( z2n - s2p1 ) );
      usn=  z2k * n[i] * s2p3;
      ucn= -z2n *   c5;
      uck= -z2k * ( c5 + 2*s2p1 );
      
      un= u0n - 2 * x * ( ( ucn - beta*us )*cp - ( usn + beta*uc )*sp ) + u2n*x2;
      uk= u0k - 2 * x * ( uck*cp - usk*sp - beta*u1 ) + x2 * ( u2k - 2*beta*u2 );
      
      // derivatives of v
      
      v0n= k2 * ( z4n + s2m1 ) + nm1 * ( s2m1 + n[i] * ( z4n + c7 ) );
      v2n= k2 * ( z4n - s2m1 ) + np1 * ( s2m1 + n[i] * ( z4n - c7 ) );
      vsk=  ( 1 + c12 ) * s2m1;
      vsn=  z2k * n[i] * s2m1;
      v0k=  z2k * ( c8 + n[i] * ( z2n + s2m1 ) );
      v2k=  z2k * ( c8 + n[i] * ( z2n - s2m1 ) );
      vcn= -z2n * ( s2p1 - 2*k2pn2 );
      vck=  z2k * ( s2p1 + 2*k2pn2 );
      
      vn= v0n - 2 * x * ( (vcn-beta*vs)*cp - (vsn+beta*vc)*sp ) + v2n*x2;
      vk= v0k - 2 * x * ( vck*cp - vsk*sp - beta*v1 ) + x2 * ( v2k - 2*beta*v2 );
      
      // derivatives of q
      
      q0n= 2 * ( 2*n[i] + sp1 ) * ( k2 + np1 * nps );
      q2n= 2 * ( 2*n[i] - sp1 ) * ( k2 + nm1 * nms );
      qsk= 2 * ( 3*k2 + n2 - s[i] ) * sp1;
      q0k=  z2k * ( z2k2 + np1E2 + npsE2 );
      q2k=  z2k * ( z2k2 + nm1E2 + nmsE2 );
      qsn=  z2k * 2 * n[i] * sp1;
      qcn= -z2n *   c5;
      qck= -z2k * ( c5 + z4s );
      
      qn= q0n - 2 * x * ( (qcn-beta*qs)*cp - (qsn+beta*qc)*sp ) + q2n*x2;
      qk= q0k - 2 * x * ( qck*cp - qsk*sp - beta*q1 ) + x2 * ( q2k - 2*beta*q2 );
      
      // derivatives of p
      
      psn=  z2k * ( k2 * n2p1 + 5 * n2 * k2pn2 - s2 );
      psk=  z2n * ( n2 * k2p1 + 5 * k2 * k2pn2 - s2 );
      pcn= -z2n * ( s2 + z2k2 * ( c10 -   k2/2 )       - z2n2 * ( s2p1 - 3*n2/2 ) );
      pck=  z2k * ( s2 + z2k2 * ( c9  + 3*k2/2 )       + z2n2 *   c11 );
      p1n=  z2n * ( s2 + z2k2 * ( c9  -   k2/2 + z4s ) + z2n2 * ( s2p1 + 3*n2/2 + z4s ) );
      p1k= -z2k * ( s2 + z2k2 * ( c10 + 3*k2/2 - z4s ) - z2n2 * ( c11 + z4s ) );
      
      pn= ( ( pcn - beta*ps )*cp - ( psn + beta*pc )*sp ) + \
	( ( pcn + beta*ps )*cp + ( psn - beta*pc )*sp ) * x2 + 2 * p1n * x;
      pk= ( pck*cp - psk*sp ) + 2 * ( p1k - beta*p1 ) * x + \
	( pck*cp + psk*sp - 2*beta*p2) * x2; 
      
      p= ( pc*cp - ps*sp ) + 2 * p1 * x + p2 * x2;
      
      // derivatives of R
      
      qq = q * q; uu = u * u;
      yn = y0 * ( pn*q - p*qn ) / qq;
      drdn[i] = ( (vn + yn) * u - (w + y) * un ) / uu;
      
      yk = y0 * ( -beta*pq + ( pk*q - p*qk ) / qq );
      drdk[i] = ( (vk + yk) * u - (w + y) * uk ) / uu;
    }
}

//*******************************************************************
//*******************************************************************

double compqe(int nobs,double datobs[],double dat[])
{
  int i;
  double diff,qe;

  qe = 0.0;

  for ( i = 0; i < nobs; i ++ )
    {
      diff = datobs[i] - dat[i];
      qe = qe + diff * diff;
    }

  return qe;
}

//*******************************************************************
//*******************************************************************

double compdqe(int nobs,double datobs[],double dat[],double dqeddat[])
{
  int i;
  double diff,qe;

  qe = 0.0;

  for ( i = 0; i < nobs; i ++ )
    {
      diff = datobs[i] - dat[i];
      qe = qe + diff * diff;
      dqeddat[i] = - 2 * diff;
    }

  return qe;
}

//*******************************************************************
//*******************************************************************

double evalf(int ydim,double y[])
{
  int i,j,nlayers,slayer,nfilms,nobs;
  double h,qe;
  char datatype;

  double ip[nfilmsmax],d[nfilmsmax],lambda[nobsmax],n[nlayersmax][nobsmax],
    k[nlayersmax][nobsmax],t[nobsmax],tobs[nobsmax],r[nobsmax],
    robs[nobsmax],z[nfilmsmax][nobsmax],v[nfilmsmax][nobsmax];

  nlayers = film.nlayers;

  slayer = film.slayer;

  nfilms = nlayers - 3;

  for ( i = 0; i < nfilms; i ++ )
    {
      ip[i] = film.ip[i];
      d[i]  = film.d[i];
    }

  nobs = film.nobs;

  datatype = film.datatype;

  if ( datatype == 'T' || datatype == 'B' )
    for ( j = 0; j < nobs; j ++ )
      tobs[j] = film.tobs[j];

  if ( datatype == 'R' || datatype == 'B' )
    for ( j = 0; j < nobs; j ++ )
      robs[j] = film.robs[j];

  h = film.h;

  for ( j = 0; j < nobs; j ++ )
    lambda[j] = film.lambda[j];

  for ( i = 1; i <= nfilms; i ++ )
    for ( j = 0; j < nobs; j ++ )
      {
	z[i-1][j] = y[2 * nobs * (i-1) + j];
	v[i-1][j] = y[2 * nobs * (i-1) + nobs + j];
      }

  compn(nlayers,slayer,nobs,z,n);
  compk(nlayers,slayer,nobs,ip,lambda,v,k);

  if ( nfilms == 1 && slayer == 2 )
    {
      if ( datatype == 'T' || datatype == 'B' )
        compt(nobs,n[1],k[1],d[0],n[2],lambda,t);
      if ( datatype == 'R' || datatype == 'B' )
        compr(nobs,n[1],k[1],d[0],n[2],lambda,r);
    }
  else
    comptr(nlayers,slayer,nobs,n,k,d,lambda,t,r);

  qe = 0.0;

  if ( datatype == 'T' || datatype == 'B' )
    qe += compqe(nobs,tobs,t);

  if ( datatype == 'R' || datatype == 'B' )
    qe += compqe(nobs,robs,r);

  return qe;
}

//*******************************************************************
//*******************************************************************

double evalg(int ydim,double y[],double g[])
{
  int i,j,nlayers,slayer,nfilms,nobs;
  double h,qe;
  char datatype;

  double ip[nfilmsmax],d[nfilmsmax],lambda[nobsmax],n[nlayersmax][nobsmax],
    k[nlayersmax][nobsmax],t[nobsmax],tobs[nobsmax],r[nobsmax],
    robs[nobsmax],z[nfilmsmax][nobsmax],v[nfilmsmax][nobsmax],
    dqedt[nobsmax],dtdn[nfilmsmax][nobsmax],dtdk[nfilmsmax][nobsmax],
    dqedr[nobsmax],drdn[nfilmsmax][nobsmax],drdk[nfilmsmax][nobsmax],
    dndz[nfilmsmax][nobsmax],dkdv[nfilmsmax][nobsmax],
    dqedn[nfilmsmax][nobsmax],dqedk[nfilmsmax][nobsmax];

  nlayers = film.nlayers;

  slayer = film.slayer;

  nfilms = nlayers - 3;

  for ( i = 0; i < nfilms; i ++ )
    {
      ip[i] = film.ip[i];
      d[i]  = film.d[i];
    }

  nobs = film.nobs;

  datatype = film.datatype;

  if ( datatype == 'T' || datatype == 'B' )
    for ( j = 0; j < nobs; j ++ )
      tobs[j] = film.tobs[j];

  if ( datatype == 'R' || datatype == 'B' )
    for ( j = 0; j < nobs; j ++ )
      robs[j] = film.robs[j];

  h = film.h;

  for ( j = 0; j < nobs; j ++ )
    lambda[j] = film.lambda[j];

  for ( i = 1; i <= nfilms; i ++ )
    for ( j = 0; j < nobs; j ++ )
      {
	z[i-1][j] = y[2 * nobs * (i-1) + j];
	v[i-1][j] = y[2 * nobs * (i-1) + nobs + j];
      }

  compdn(nlayers,slayer,nobs,z,n,dndz);
  compdk(nlayers,slayer,nobs,ip,lambda,v,k,dkdv);

  if ( nfilms == 1 && slayer == 2 )
    {
      if ( datatype == 'T' || datatype == 'B' )
        compdt(nobs,n[1],k[1],d[0],n[2],lambda,t,dtdn[0],dtdk[0]);
      if ( datatype == 'R' || datatype == 'B' )
        compdr(nobs,n[1],k[1],d[0],n[2],lambda,r,drdn[0],drdk[0]);
    }
  else
    compdtr(nlayers,slayer,nobs,n,k,d,lambda,t,r,dtdn,dtdk,drdn,drdk);

  qe = 0.0;

  if ( datatype == 'T' || datatype == 'B' )
    qe += compdqe(nobs,tobs,t,dqedt);

  if ( datatype == 'R' || datatype == 'B' )
    qe += compdqe(nobs,robs,r,dqedr);

  for ( i = 0; i < nfilms; i ++ )
    {
      for ( j = 0; j < nobs; j ++ )
        {
	  dqedn[i][j] = 0;
	  dqedk[i][j] = 0;
        }

      if ( datatype == 'T' || datatype == 'B' )
	for ( j = 0; j < nobs; j ++ )
        {
	  dqedn[i][j] += dqedt[j] * dtdn[i][j];
	  dqedk[i][j] += dqedt[j] * dtdk[i][j];
        }

      if ( datatype == 'R' || datatype == 'B' )
	for ( j = 0; j < nobs; j ++ )
        {
	  dqedn[i][j] += dqedr[j] * drdn[i][j];
	  dqedk[i][j] += dqedr[j] * drdk[i][j];
        }

      for ( j = 0; j < nobs-2; j ++ )
        {
	  dqedn[i][j+1] += dqedn[i][j] * 2;
	  dqedn[i][j+2] -= dqedn[i][j];

	  dqedk[i][j+1] += dqedk[i][j] * 2;
	  dqedk[i][j+2] -= dqedk[i][j];
        }
    
      dqedn[i][nobs-1] += dqedn[i][nobs-2];
      dqedk[i][nobs-1] += dqedk[i][nobs-2];

      for ( j = 0; j < nobs; j ++ )
        {
	  g[2 * nobs * i + j]        = dqedn[i][j] * dndz[i][j];
	  g[2 * nobs * i + nobs + j] = dqedk[i][j] * dkdv[i][j];
        }
    }

  return qe;
}

//*******************************************************************
//*******************************************************************

int Iniv(int n,double v[],double firste[])
{
  int i;

  for ( i = 0; i < n; i ++ )
    v[i] = firste[i];

  return 0;
}

//*******************************************************************
//*******************************************************************

int Incrv(int n,double v[],double firste[],double laste[],double estep[])
{
  int i;

  i = n - 1;

  while ( i >= 0 )
    {
      if ( v[i] < laste[i] )
	{
	  v[i] = v[i] + estep[i];
	  break;
	}
      else
	{
	  v[i] = firste[i];
	  i --;
	}
    }

  return i < 0;
}

//*******************************************************************
//*******************************************************************

void IniFilm(char id[],int nlayers,int slayer,int substrate,char datatype,
int nobs,double lini,double lfin,char datfile[])
  {
  int i,j,in_nobs;
  double wa,wb,a,b,w;
  FILE * fp;

  double in_lambda[ninobsmax],in_tobs[ninobsmax],in_robs[ninobsmax];
  char message[255];

  // Set film id, number of films, number of observations and substrate

  strcpy(film.id,id);

  if ( nlayers > nlayersmax )
    {
      sprintf(message,"The number of layers can not exceed nlayersmax = %d.",nlayersmax);
      Error(message);
    }

  film.nlayers = nlayers;

  if ( slayer < 0 || slayer > nlayers - 1 )
    {
      sprintf(message,"The layer of the substrate must be a number between %d and %d.",0,nlayers-1);
      Error(message);
    }

  film.slayer = slayer;

  film.substrate = substrate;

  if ( datatype != 'T' && datatype != 'R' && datatype != 'B' )
    Error("Invalid datatype (T=transmittance, R=reflectance, B=both)");

  film.datatype = datatype;

  if ( nobs > nobsmax )
    {
      sprintf(message,"The number of observations (used in the optimization process) can not exceed nobsmax = %d.",nobsmax);
      Error(message);
    }

  film.nobs = nobs;

  // Set wavelengths where the optical parameters will be estimated

  film.h = ( lfin - lini ) / ( nobs - 1 );

  for ( i = 0; i < nobs; i ++ )
    film.lambda[i] = lini + i * film.h;
  
  if ( film.lambda[nobs-1] != lfin )
    film.lambda[nobs-1] = lfin;
  
  // Set substrate refractive index

  // Glass substrate (transparent in the range [350nm,2000nm])
  if ( substrate == 10 )

      for ( i = 0; i < nobs; i ++ )
          film.s[i] =
          sqrt(1.0+1.0/(0.7568-7930.0/(film.lambda[i]*film.lambda[i])));

  // Crystalline silicon substrate (transparent in the range [1250nm,2600nm])
  else if ( substrate == 20 )

      for ( i = 0; i < nobs; i ++ )
          film.s[i]=
          3.71382 - ( 8.69123e-05 - ( 2.47125e-08 + 1.04677e-11 * 
          film.lambda[i] ) * film.lambda[i] ) * film.lambda[i];

  // Crystalline quartz substrate (transparent in the range [200nm,1500nm])
  else if ( substrate == 30 )

      for ( i = 0; i < nobs; i ++ )
	  film.s[i]=
	  sqrt( 1.0 + 1.34157 * film.lambda[i] * film.lambda[i] / 
                      ( film.lambda[i] * film.lambda[i] - 8385.79 ) );

  // Glass slides substrate (transparent in the range [360nm,800nm])
  else if ( substrate == 40 )

      for ( i = 0; i < nobs; i ++ )
          film.s[i]=
	  sqrt( 1.0 + 1.43739 * film.lambda[i] * film.lambda[i] / 
                      ( film.lambda[i] * film.lambda[i] - 44220.2 ) );

  // Borosilicate substrate (transparent in the range [300nm,2600nm])
  else if ( substrate == 50 )

      for ( i = 0; i < nobs; i ++ )
	 film.s[i]= 1.51507 + 0.019 * exp( - ( film.lambda[i] - 435.8 ) / 175.09383 );

  // Amorphous Quartz substrate [transparent in the range 200nm,1500nm]	 
  else if ( substrate == 60 )

      for ( i = 0; i < nobs; i ++ )
	 film.s[i]=
	  sqrt( 1.0 + 0.18394 * film.lambda[i] * film.lambda[i] / 
                      ( film.lambda[i] * film.lambda[i] - 18231.83 ) );

   else Error("Invalid substrate (10=Glass, 20=Crystalline silicon, 30=Crystalline quartz, 40=Glass slides, 50=Borosilicate, 60=Amorphous quartz)");

  // Read the observed data (transmittance and/or reflectance) and
  // interpolate if necessary.

  fp= fopen(datfile,"rt");

  if ( fp == NULL )
    {
      strcat(strcat(strcpy(message,"File "),datfile)," does not exists");
      Error(message);
    }

  fscanf(fp, "%d", &in_nobs);

  if ( in_nobs > ninobsmax )
    {
      sprintf(message,"The number of observations in file %s can not exceed ninobsmax = %d.",datfile,ninobsmax);
      Error(message);
    }

  switch(datatype)
    {
    case 'T':

      for ( i = 0; i < in_nobs; i ++ )
	fscanf(fp, "%lf %lf", & in_lambda[i], & in_tobs[i]);

      break;

    case 'R':

      for ( i = 0; i < in_nobs; i ++ )
	fscanf(fp, "%lf %lf", & in_lambda[i], & in_robs[i]);

      break;

    case 'B':

      for ( i = 0; i < in_nobs; i ++ )
	fscanf(fp, "%lf %lf %lf", & in_lambda[i], & in_tobs[i], & in_robs[i]);

      break;
    }

  fclose(fp);

  if ( film.lambda[0]           < in_lambda[0]         || 
       film.lambda[film.nobs-1] > in_lambda[in_nobs-1] )
      Error("Invalid data file (wavelenght out of range)");

  for ( i = 0; i < nobs; i ++ )
      {
      j = 0;
      while( j < in_nobs && in_lambda[j] < film.lambda[i] )
          j++;

      if ( in_lambda[j] == film.lambda[i] )
	{
          if ( datatype == 'T' || datatype == 'B') film.tobs[i] = in_tobs[j];
          if ( datatype == 'R' || datatype == 'B') film.robs[i] = in_robs[j];
	}

      else
          {
          wa = in_lambda[j-1];
          wb = in_lambda[j];
          w  = film.lambda[i];

          if ( datatype == 'T' || datatype == 'B')
	    { 
	      a  = in_tobs[j-1];
	      b  = in_tobs[j];
	      film.tobs[i]= ( ( wb - w ) * a +( w - wa ) * b ) / ( wb - wa );
	    }
	  
          if ( datatype == 'R' || datatype == 'B')
	    { 
	      a  = in_robs[j-1];
	      b  = in_robs[j];
	      film.robs[i]= ( ( wb - w ) * a +( w - wa ) * b ) / ( wb - wa );
	    }
          }
      }
  }

//*******************************************************************
//*******************************************************************

void IniPoint(double n0[],double nf[],double k0[],double y[])
  {
  int i,j,nobs,nlayers,nfilms;
  double n[nobsmax],k[nobsmax],z[nobsmax],v[nobsmax],nk0,nk1,nk2,k1,k2,h;

  h       = film.h;
  nobs    = film.nobs;
  nlayers = film.nlayers;
  nfilms  = nlayers - 3;

  nk0 = 0.00 * nobs;
  nk1 = 0.20 * nobs;
  nk2 = 1.00 * nobs;

  for ( i = 0; i < nfilms; i ++ )
    {
    k1 = 0.1   * k0[i];
    k2 = 1e-10 * k0[i];

    for ( j = 0; j < nobs; j++ )
      {
        n[j] = n0[i] - ( n0[i] - nf[i] ) * j / ( nobs - 1 );
        if ( j < nk1 ) 
          k[j] = k0[i] - ( k0[i] - k1 ) * ( j - nk0 ) / ( nk1 - nk0 );
        else 
          k[j] = k1    - ( k1    - k2 ) * ( j + 1 - nk1 ) / ( nk2 - nk1 );
      }

    z[nobs-1] = sqrt( n[nobs-1] - 1 );
    z[nobs-2] = sqrt( ( n[nobs-2] - n[nobs-1] ) / h );

    v[nobs-1] = sqrt( k[nobs-1] );
    v[nobs-2] = sqrt( ( k[nobs-2] - k[nobs-1] ) / h );

    for ( j = 0; j < nobs - 2; j ++ )
      {
      z[j] = sqrt( max( 0.0, ( n[j] - 2 * n[j+1] + n[j+2] ) / ( h * h ) ) );
      v[j] = sqrt( max( 0.0, ( k[j] - 2 * k[j+1] + k[j+2] ) / ( h * h ) ) );
      }

    for ( j = 0; j < nobs; j ++ )
      {
      y[2 * nobs * i + j]        = z[j];
      y[2 * nobs * i + nobs + j] = v[j];
      }
    }

  }

//*******************************************************************
//*******************************************************************
void GenReport(double bestqe,double bestd[],double bestip[], 
char solfile[],char inffile[])
  {
  int i,j,nlayers,slayer,nfilms,nobs;
  double h;
  char datatype;
  FILE * fp;

  double lambda[nobsmax],n[nlayersmax][nobsmax],k[nlayersmax][nobsmax],
    t[nobsmax],r[nobsmax],z[nfilmsmax][nobsmax],v[nfilmsmax][nobsmax],
    y[2*nfilmsmax*nobsmax];
  char message[255];

  h       = film.h;
  nlayers = film.nlayers;
  slayer  = film.slayer;
  nfilms  = nlayers - 3;
  nobs    = film.nobs;

  for ( i = 0; i < nobs; i ++ )
    lambda[i] = film.lambda[i];

  fp = fopen(solfile,"rt");

  if ( fp == NULL )
    {
      strcat(strcat(strcpy(message,"File "),solfile)," does not exists");
      Error(message);
    }

  for ( i = 0; i < 2 * nfilms * nobs; i ++ )
    fscanf(fp,"%lf",&y[i]);

  fclose(fp);

  fp= fopen(inffile, "wt");

  fprintf(fp, "ESTIMATED THICKNESSES\n");

  for ( i = 1; i <= nfilms; i ++ )
      fprintf(fp,"\nThickness of film %d = %.2f nm",i,bestd[i-1]);

  fprintf(fp, "\n\nESTIMATED INFLEXION POINTS\n");

  for ( i = 1; i <= nfilms; i ++ )
      fprintf(fp,"\nInflexion point of the absorption coefficient of film %d = %.2f nm",
	      i,bestip[i-1]);

  fprintf(fp, "\n\nESTIMATED REFRACTIVE INDICES AND ABSORPTION COEFFICIENTS\n");

  for ( i = 1; i <= nfilms; i ++ )
    for ( j = 0; j < nobs; j ++ )
      {
	z[i-1][j] = y[2 * nobs * (i-1) + j];
	v[i-1][j] = y[2 * nobs * (i-1) + nobs + j];
      }

  compn(nlayers,slayer,nobs,z,n);
  compk(nlayers,slayer,nobs,bestip,lambda,v,k);
  comptr(nlayers,slayer,nobs,n,k,bestd,lambda,t,r);

  // print front films constants
  for ( i = 1; i <= slayer-1; i ++ )
    {
      fprintf(fp,"\nFilm %d\n",i);
      fprintf(fp, "\nlambda                     n                          kappa\n");
      for ( j = 0; j < nobs; j ++)
	fprintf(fp,"%.20e %.20e %.20e\n",lambda[j],n[i][j],k[i][j]);
    }

  // print back films constanst
  for ( i = slayer+1; i <= nlayers-2; i ++ )
    {
      fprintf(fp,"\nFilm %d\n",i);
      fprintf(fp, "\nlambda                     n                          kappa\n");
      for ( j = 0; j < nobs; j ++)
	fprintf(fp,"%.20e %.20e %.20e\n",lambda[j],n[i][j],k[i][j]);
    }

  datatype = film.datatype;

  if ( datatype == 'T' )
    {
      fprintf(fp, "\n\nTRANSMITTANCE WITH THE ESTIMATED THICKNESSES AND OPTICAL PARAMETERS\n");
      fprintf(fp, "\nlambda                     t\n");
      for ( i = 0; i < nobs; i ++)
	fprintf(fp,"%.20e %.20e\n",lambda[i],t[i]);
    }
  else if ( datatype == 'R' )
    {
      fprintf(fp, "\n\nREFLECTANCE WITH THE ESTIMATED THICKNESSES AND OPTICAL PARAMETERS\n");
      fprintf(fp, "\nlambda                     r\n");
      for ( i = 0; i < nobs; i ++)
	fprintf(fp,"%.20e %.20e\n",lambda[i],r[i]);
    }
  else
    {
      fprintf(fp, "\n\nTRANSMITTANCE/REFLECTANCE WITH THE ESTIMATED THICKNESSES AND OPTICAL PARAMETERS\n");
      fprintf(fp, "\nlambda                     t                          r\n");
      for ( i = 0; i < nobs; i ++)
	fprintf(fp,"%.20e %.20e %.20e\n",lambda[i],t[i],r[i]);
    }

  fprintf(fp, "\n\nQUADRATIC ERROR = %e\n",bestqe);

  fclose(fp);
  }

//*******************************************************************
//*******************************************************************

int Error(char * buffer)
  {
  fprintf(stderr, "\n%s", buffer);
  fprintf(stderr, "\nHit any key to end...\n");

  getchar();
  exit(-1);

  return(TRUE);
  }

//*******************************************************************
//*******************************************************************

char firstc (char s[])
{
  return s[0];
}

//*******************************************************************
//*******************************************************************

void comptr(int nlayers,int slayer,int nobs,double n[][nobsmax], 
double k[][nobsmax],double d[],double lambda[],double t[],double r[])
{
  int i,j;
  double tmp0, tmp1, beta, phi, cph, sph, x, fat, fatB, fatBC,
    br00, br01, br10, br11, bi00, bi01, bi10, bi11, 
    cr00, cr01, cr10, cr11, ci00, ci01, ci10, ci11, 
    xr00, xr01, xr10, xr11, xi00, xi01, xi10, xi11,
    yr00, yr01, yr10, yr11, yi00, yi01, yi10, yi11,
    Id, Cd, Sd, Iu, Cu, Su, rho, p, f,
    mb11, mb01, mb10, mb00, mc11, mc10;

  for ( j = 0; j < nobs; j ++ )
    {
      // evaluate matrix B

      beta = 0.5  / ( n[1][j] * n[1][j] + k[1][j] * k[1][j] );
      tmp0 = beta * ( n[0][j] * n[1][j] + k[0][j] * k[1][j] );
      tmp1 = beta * ( n[0][j] * k[1][j] - n[1][j] * k[0][j] );

      br00 = 0.5  + tmp0;
      br01 = 0.5  - tmp0;
      br10 = br01;
      br11 = br00;

      bi00 =        tmp1;
      bi01 =      - tmp1;
      bi10 =      - tmp1;
      bi11 =        tmp1;

      fatB = 0;
      for ( i = 1; i <= slayer - 1; i ++ )
	{
	  beta  = 2 * PI * d[i-1] / lambda[j];
	  phi   = beta * n[i][j];
          cph   = cos ( phi );
          sph   = sin ( phi );

	  // avoid overflow

	  fat   = beta * k[i][j];
	  x     = exp( -2 * fat );
	  fatB += fat;

          xr00 = ( br00 * cph + bi00 * sph ) * x;
          xi00 = ( bi00 * cph - br00 * sph ) * x;
          xr01 = ( br01 * cph + bi01 * sph ) * x;
          xi01 = ( bi01 * cph - br01 * sph ) * x;

          xr10 =   br10 * cph - bi10 * sph;
          xi10 =   br10 * sph + bi10 * cph;
          xr11 =   br11 * cph - bi11 * sph;
          xi11 =   br11 * sph + bi11 * cph;

	  beta = 0.5  / ( n[i+1][j] * n[i+1][j] + k[i+1][j] * k[i+1][j] );
	  tmp0 = beta * ( n[i][j]   * n[i+1][j] + k[i][j]   * k[i+1][j] );
	  tmp1 = beta * ( n[i][j]   * k[i+1][j] - n[i+1][j] * k[i][j]   );

	  yr00 =  0.5 * ( xr00 + xr10 );
	  yi00 =  0.5 * ( xi00 + xi10 );
	  yr01 =  0.5 * ( xr01 + xr11 );
	  yi01 =  0.5 * ( xi01 + xi11 );

	  yr10 = tmp0 * ( xr00 - xr10 ) - tmp1 * ( xi00 - xi10 );
	  yi10 = tmp0 * ( xi00 - xi10 ) + tmp1 * ( xr00 - xr10 );
	  yr11 = tmp0 * ( xr01 - xr11 ) - tmp1 * ( xi01 - xi11 );
	  yi11 = tmp0 * ( xi01 - xi11 ) + tmp1 * ( xr01 - xr11 );

	  br00 = yr00 + yr10;
	  bi00 = yi00 + yi10;
	  br01 = yr01 + yr11;
	  bi01 = yi01 + yi11;
		      	                 	         	   	   
	  br10 = yr00 - yr10;
	  bi10 = yi00 - yi10;
	  br11 = yr01 - yr11;
	  bi11 = yi01 - yi11;
	}

      // evaluate matrix C

      beta = 0.5  / ( n[slayer+1][j] * n[slayer+1][j] + k[slayer+1][j] * k[slayer+1][j] );
      tmp0 = beta * ( n[slayer][j]   * n[slayer+1][j] + k[slayer][j]   * k[slayer+1][j] );
      tmp1 = beta * ( n[slayer][j]   * k[slayer+1][j] - n[slayer+1][j] * k[slayer][j] );

      cr00 = 0.5  + tmp0;
      cr01 = 0.5  - tmp0;
      cr10 = cr01;
      cr11 = cr00;
      
      ci00 =        tmp1;
      ci01 =      - tmp1;
      ci10 =      - tmp1;
      ci11 =        tmp1;

      fatBC= fatB;
      for ( i = slayer + 1; i <= nlayers - 2; i ++ )
	{
	  beta  = 2 * PI * d[i-2] / lambda[j];
	  phi   = beta * n[i][j];
          cph   = cos ( phi );
          sph   = sin ( phi );

	  // avoid overflow

	  fat    = beta * k[i][j];
	  x      = exp( -2 * fat );
	  fatBC += fat;

          xr00 = ( cr00 * cph + ci00 * sph ) * x;
          xi00 = ( ci00 * cph - cr00 * sph ) * x;
          xr01 = ( cr01 * cph + ci01 * sph ) * x;
          xi01 = ( ci01 * cph - cr01 * sph ) * x;

          xr10 =   cr10 * cph - ci10 * sph;
          xi10 =   cr10 * sph + ci10 * cph;
          xr11 =   cr11 * cph - ci11 * sph;
          xi11 =   cr11 * sph + ci11 * cph;

	  beta = 0.5  / ( n[i+1][j] * n[i+1][j] + k[i+1][j] * k[i+1][j] );
	  tmp0 = beta * ( n[i][j]   * n[i+1][j] + k[i][j]   * k[i+1][j] );
	  tmp1 = beta * ( n[i][j]   * k[i+1][j] - n[i+1][j] * k[i][j]   );

	  yr00 =  0.5 * ( xr00 + xr10 );
	  yi00 =  0.5 * ( xi00 + xi10 );
	  yr01 =  0.5 * ( xr01 + xr11 );
	  yi01 =  0.5 * ( xi01 + xi11 );

	  yr10 = tmp0 * ( xr00 - xr10 ) - tmp1 * ( xi00 - xi10 );
	  yi10 = tmp0 * ( xi00 - xi10 ) + tmp1 * ( xr00 - xr10 );
	  yr11 = tmp0 * ( xr01 - xr11 ) - tmp1 * ( xi01 - xi11 );
	  yi11 = tmp0 * ( xi01 - xi11 ) + tmp1 * ( xr01 - xr11 );

	  cr00 = yr00 + yr10;
	  ci00 = yi00 + yi10;
	  cr01 = yr01 + yr11;
	  ci01 = yi01 + yi11;
		      	                 	         	   	   
	  cr10 = yr00 - yr10;
	  ci10 = yi00 - yi10;
	  cr11 = yr01 - yr11;
	  ci11 = yi01 - yi11;
	}

      // assembling \rho

      mb00= br00 * br00 + bi00 * bi00;
      mb01= br01 * br01 + bi01 * bi01;
      mb10= br10 * br10 + bi10 * bi10;
      mb11= br11 * br11 + bi11 * bi11;
      p   = br00 * br11 - bi00 * bi11;

      mc11= cr11 * cr11 + ci11 * ci11;
      mc10= cr10 * cr10 + ci10 * ci10;

      rho = mc11 * mb11 - mc10 * mb01;

      // evaluate transmittance

      t[j] = exp(-2*fatBC) / rho;

      // evaluate reflectanec
      
      p    = mc11 * mb10 
	   - mc10 * ( mb00 - 2 * exp(-2*fatB) * p / ( n[slayer][j] * mb11 ) );
      r[j] = p / rho;
    }
}

//*******************************************************************
//*******************************************************************

void assemble(double Lr[][2], double Li[][2], 
	      double Rr[][2], double Ri[][2], 
	      double Xr[][2], double Xi[][2], 
	      double dBr[2][2], double dBi[2][2])
{
  double c1,c2,c3,c4,c5,c6,c7,c8;

  c5 = Xi[0][0] * Ri[0][0] + Xi[0][1] * Ri[1][0] 
     - Xr[0][0] * Rr[0][0] - Xr[0][1] * Rr[1][0];
  c6 = Xi[1][0] * Ri[0][0] + Xi[1][1] * Ri[1][0] 
     - Xr[1][0] * Rr[0][0] - Xr[1][1] * Rr[1][0];
  c7 = Xi[0][0] * Ri[0][1] + Xi[0][1] * Ri[1][1] 
     - Xr[0][0] * Rr[0][1] - Xr[0][1] * Rr[1][1];
  c8 = Xi[1][0] * Ri[0][1] + Xi[1][1] * Ri[1][1] 
     - Xr[1][0] * Rr[0][1] - Xr[1][1] * Rr[1][1];

  c1 = Li[0][0] * Xr[0][0] + Li[0][1] * Xr[1][0]; 
  c2 = Li[0][0] * Xr[0][1] + Li[0][1] * Xr[1][1];
  c3 = Li[0][0] * Xi[0][0] + Li[0][1] * Xi[1][0];
  c4 = Li[0][0] * Xi[0][1] + Li[0][1] * Xi[1][1];    
  
  dBr[0][0] = - Ri[0][0] * c1 - Ri[1][0] * c2 // dBr[0][0]
    - Rr[0][0] * c3 - Rr[1][0] * c4
    - Lr[0][0] * c5 - Lr[0][1] * c6;
  
  dBr[0][1] = - Ri[0][1] * c1 - Ri[1][1] * c2 // dBr[0][1]
    - Rr[0][1] * c3 - Rr[1][1] * c4
    - Lr[0][0] * c7 - Lr[0][1] * c8;
  
  c1 = Li[1][0] * Xr[0][0] + Li[1][1] * Xr[1][0];
  c2 = Li[1][0] * Xr[0][1] + Li[1][1] * Xr[1][1];
  c3 = Li[1][0] * Xi[0][0] + Li[1][1] * Xi[1][0];
  c4 = Li[1][0] * Xi[0][1] + Li[1][1] * Xi[1][1];
  
  dBr[1][0] = - Ri[0][0] * c1 - Ri[1][0] * c2  // dBr[1][0]
    - Rr[0][0] * c3 - Rr[1][0] * c4 
    - Lr[1][0] * c5 - Lr[1][1] * c6;
  
  dBr[1][1] = - Ri[0][1] * c1 - Ri[1][1] * c2 // dBr[1][1]
    - Rr[0][1] * c3 - Rr[1][1] * c4
    - Lr[1][0] * c7 - Lr[1][1] * c8;
  
  c1 = Lr[0][0] * Xr[0][0] + Lr[0][1] * Xr[1][0]; 
  c2 = Lr[0][0] * Xr[0][1] + Lr[0][1] * Xr[1][1];
  c3 = Lr[0][0] * Xi[0][0] + Lr[0][1] * Xi[1][0];
  c4 = Lr[0][0] * Xi[0][1] + Lr[0][1] * Xi[1][1];
  
  dBi[0][0] =   Ri[0][0] * c1 + Ri[1][0] * c2 // dBi[0][0]
    + Rr[0][0] * c3 + Rr[1][0] * c4
    - Li[0][0] * c5 - Li[0][1] * c6;
  
  dBi[0][1] =   Ri[0][1] * c1 + Ri[1][1] * c2 // dBi[0][1]
    + Rr[0][1] * c3 + Rr[1][1] * c4
    - Li[0][0] * c7 - Li[0][1] * c8;
  
  c1 = Lr[1][0] * Xr[0][0] + Lr[1][1] * Xr[1][0];
  c2 = Lr[1][0] * Xr[0][1] + Lr[1][1] * Xr[1][1];
  c3 = Lr[1][0] * Xi[0][0] + Lr[1][1] * Xi[1][0];
  c4 = Lr[1][0] * Xi[0][1] + Lr[1][1] * Xi[1][1];
  
  dBi[1][0] = + Ri[0][0] * c1 + Ri[1][0] * c2 // dBi[1][0]
    + Rr[0][0] * c3 + Rr[1][0] * c4
    - Li[1][0] * c5 - Li[1][1] * c6;
  
  dBi[1][1] =   Ri[0][1] * c1 + Ri[1][1] * c2 // dBi[1][1]
    + Rr[0][1] * c3 + Rr[1][1] * c4
    - Li[1][0] * c7 - Li[1][1] * c8;
}

//*******************************************************************
//*******************************************************************

void tance(int final, int nobs, double n[], double k[], 
	   double d[], double lambda, double* fatorM, 
	   double Mr[2][2], double Mi[2][2],
	   double Lr[nlayersmax-1][2][2], double Li[nlayersmax-1][2][2], 
	   double Rr[nlayersmax-1][2][2], double Ri[nlayersmax-1][2][2], 
	   double dXrdn[nlayersmax-1][2][2], double dXidn[nlayersmax-1][2][2],
	   double dXrdk[nlayersmax-1][2][2], double dXidk[nlayersmax-1][2][2] )
{
  int a,c,i,j,l,v;

  double Xr[2][2],Xi[2][2],dXr[2][2], dXi[2][2],
    rSp[2][2],rSq[2][2],rCp[2][2],rCq[2][2],
    iSp[2][2],iSq[2][2],iCp[2][2],iCq[2][2],
    drSp[2][2],drSq[2][2],drCp[2][2],drCq[2][2],
    diSp[2][2],diSq[2][2],diCp[2][2],diCq[2][2],
    tmp0, tmp1, beta, phi, cph, sph, x, fat, fatM, 
    Id, Cd, Sd, Iu, Cu, Su, rho, nplus, nminu, Kplus, Kminu, mdl, i1,
    br00, br01, br10, br11, bi00, bi01, bi10, bi11,
    nv, Kv, xq, Kv2, Kv3, nv2, nvp, nva, Kvp, Kva, KvaKvp, Kv2Kvp, 
    nv2PKv2, nvaPnv, nvaLnv, nvPnvp, nvLnvp, nvaP2nv, nvaL2nv, c2nvPnvp, 
    c2nvLnvp, c2Kvnv, nv2Pnvanvp, nv2Lnvanvp,
    c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, 
    c11, c12, c13, c14, c15, c16, c17, c18, c19, c20,
    c21, c22, c23, c24, c25, c26, c27, c28, c29, c30,
    c31, c32, c33, c34, c35, c36, c37, c38, c39, c40,
    c41, c42, c43;

  Rr[0][0][0] = 1; // R_{0} = Identity
  Rr[0][0][1] = 0;
  Rr[0][1][0] = 0;
  Rr[0][1][1] = 1;
  Ri[0][0][0] = 0;
  Ri[0][0][1] = 0;
  Ri[0][1][0] = 0;
  Ri[0][1][1] = 0;
  
  Lr[final-1][0][0] = 1; // R_{m-3} = Identity
  Lr[final-1][0][1] = 0;
  Lr[final-1][1][0] = 0;
  Lr[final-1][1][1] = 1;
  Li[final-1][0][0] = 0;
  Li[final-1][0][1] = 0;
  Li[final-1][1][0] = 0;
  Li[final-1][1][1] = 0;
  
  // layers looping
  
  fatM = 0;
  for ( v = 1; v <= final-1; v ++ )
    {
      // -------------------------------------------------------- matrices L ----------
      
      // assign: A_v
      
      a = final + 2 - v;
      beta = 2 * ( n[a] * n[a] + k[a] * k[a] );
      c01 = n[a] * n[a-1] + k[a] * k[a-1];
      c02 = k[a] * n[a-1] - n[a] * k[a-1];
      c03 = 0.5 + c01 / beta;
      c04 = 0.5 - c01 / beta;
      c05 =       c02 / beta;

      // multiply: tmp= L_{v-1} * A_v
      
      a= final - v;
      br00 = c03 * Lr[a][0][0] + c04 * Lr[a][0][1] - c05 * ( Li[a][0][0] - Li[a][0][1] );
      br01 = c03 * Lr[a][0][1] + c04 * Lr[a][0][0] + c05 * ( Li[a][0][0] - Li[a][0][1] );
      br10 = c03 * Lr[a][1][0] + c04 * Lr[a][1][1] - c05 * ( Li[a][1][0] - Li[a][1][1] );
      br11 = c03 * Lr[a][1][1] + c04 * Lr[a][1][0] + c05 * ( Li[a][1][0] - Li[a][1][1] );
      bi00 = c03 * Li[a][0][0] + c04 * Li[a][0][1] + c05 * ( Lr[a][0][0] - Lr[a][0][1] );
      bi01 = c03 * Li[a][0][1] + c04 * Li[a][0][0] - c05 * ( Lr[a][0][0] - Lr[a][0][1] );
      bi10 = c03 * Li[a][1][0] + c04 * Li[a][1][1] + c05 * ( Lr[a][1][0] - Lr[a][1][1] );
      bi11 = c03 * Li[a][1][1] + c04 * Li[a][1][0] - c05 * ( Lr[a][1][0] - Lr[a][1][1] );

      // multiply: L_v= tmp * D_{v-1} 
      
      a       = final + 1 - v; // ajust index
      beta    = 2 * PI * d[a-1] / lambda;
      phi     = beta * n[a];
      cph     = cos( phi );
      sph     = sin( phi );
      fat     = beta * k[a]; // to avoid overflow
      x       = exp( -2 * fat );
      
      a = final - 1 - v;
      Lr[a][0][1] =   br01 * cph - bi01 * sph;
      Lr[a][1][1] =   br11 * cph - bi11 * sph;
      Li[a][0][1] =   bi01 * cph + br01 * sph;
      Li[a][1][1] =   bi11 * cph + br11 * sph;
      Lr[a][0][0] = ( br00 * cph + bi00 * sph ) * x;
      Lr[a][1][0] = ( br10 * cph + bi10 * sph ) * x;
      Li[a][0][0] = ( bi00 * cph - br00 * sph ) * x;
      Li[a][1][0] = ( bi10 * cph - br10 * sph ) * x;

      // -------------------------------------------------------- matrices R ----------
      
      // assign: A_v
      
      a = v-1; // adjust index
      beta = 2 * ( n[v] * n[v] + k[v] * k[v] );
      c01 = n[v] * n[a] + k[v] * k[a];
      c02 = k[v] * n[a] - n[v] * k[a];
      c03 = 0.5 + c01 / beta;
      c04 = 0.5 - c01 / beta;
      c05 =       c02 / beta;
      
      // multiply: tmp= A_v * R_{v-1}
	  
      br00 = c03 * Rr[a][0][0] + c04 * Rr[a][1][0] - c05 * ( Ri[a][0][0] - Ri[a][1][0] ); 
      br01 = c03 * Rr[a][0][1] + c04 * Rr[a][1][1] - c05 * ( Ri[a][0][1] - Ri[a][1][1] );
      br10 = c03 * Rr[a][1][0] + c04 * Rr[a][0][0] - c05 * ( Ri[a][1][0] - Ri[a][0][0] );
      br11 = c03 * Rr[a][1][1] + c04 * Rr[a][0][1] - c05 * ( Ri[a][1][1] - Ri[a][0][1] );
      bi00 = c03 * Ri[a][0][0] + c04 * Ri[a][1][0] + c05 * ( Rr[a][0][0] - Rr[a][1][0] );
      bi01 = c03 * Ri[a][0][1] + c04 * Ri[a][1][1] + c05 * ( Rr[a][0][1] - Rr[a][1][1] );
      bi10 = c03 * Ri[a][1][0] + c04 * Ri[a][0][0] + c05 * ( Rr[a][1][0] - Rr[a][0][0] );
      bi11 = c03 * Ri[a][1][1] + c04 * Ri[a][0][1] + c05 * ( Rr[a][1][1] - Rr[a][0][1] );

      // multiply: R_v= D_v * tmp

      beta = 2 * PI * d[a] / lambda;
      phi= beta * n[v];
      cph= cos( phi );
      sph= sin( phi );
      fat   = beta * k[v]; // to avoid overflow
      x     = exp( -2 * fat );
      fatM += fat;
	  
      Rr[v][0][0] = ( br00 * cph + bi00 * sph ) * x;
      Rr[v][0][1] = ( br01 * cph + bi01 * sph ) * x;
      Ri[v][0][0] = ( bi00 * cph - br00 * sph ) * x;
      Ri[v][0][1] = ( bi01 * cph - br01 * sph ) * x;
      Ri[v][1][0] =   bi10 * cph + br10 * sph;
      Ri[v][1][1] =   bi11 * cph + br11 * sph;
      Rr[v][1][0] =   br10 * cph - bi10 * sph;
      Rr[v][1][1] =   br11 * cph - bi11 * sph;
    }

  // ---------------------------------------------------------- matrix B ----------

  for ( v = final; v <= final+1; v ++ )
    {
      // assign: A_v

      a = v-1;
      beta= 2 * ( n[v] * n[v] + k[v] * k[v] );
      c01 = n[v] * n[a] + k[v] * k[a];
      c02 = k[v] * n[a] - n[v] * k[a];
      c03 = 0.5 + c01 / beta;
      c04 = 0.5 - c01 / beta;
      c05 =       c02 / beta;

      // multiply: tmp= A_v * R_{v-1}

      br00 = c03 * Rr[a][0][0] + c04 * Rr[a][1][0] - c05 * ( Ri[a][0][0] - Ri[a][1][0] ); 
      br01 = c03 * Rr[a][0][1] + c04 * Rr[a][1][1] - c05 * ( Ri[a][0][1] - Ri[a][1][1] );
      br10 = c03 * Rr[a][1][0] + c04 * Rr[a][0][0] - c05 * ( Ri[a][1][0] - Ri[a][0][0] );
      br11 = c03 * Rr[a][1][1] + c04 * Rr[a][0][1] - c05 * ( Ri[a][1][1] - Ri[a][0][1] );
      bi00 = c03 * Ri[a][0][0] + c04 * Ri[a][1][0] + c05 * ( Rr[a][0][0] - Rr[a][1][0] );
      bi01 = c03 * Ri[a][0][1] + c04 * Ri[a][1][1] + c05 * ( Rr[a][0][1] - Rr[a][1][1] );
      bi10 = c03 * Ri[a][1][0] + c04 * Ri[a][0][0] + c05 * ( Rr[a][1][0] - Rr[a][0][0] );
      bi11 = c03 * Ri[a][1][1] + c04 * Ri[a][0][1] + c05 * ( Rr[a][1][1] - Rr[a][0][1] );

      // multiply: R_v= D_v * tmp

      beta = 2 * PI * d[a] / lambda;
      phi= beta * n[v];
      cph= cos( phi );
      sph= sin( phi );
      fat   = beta * k[v]; // to avoid overflow
      x     = exp( -2 * fat );
      fatM += fat;
	  
      Rr[v][0][0] = ( br00 * cph + bi00 * sph ) * x;
      Rr[v][0][1] = ( br01 * cph + bi01 * sph ) * x;
      Ri[v][0][0] = ( bi00 * cph - br00 * sph ) * x;
      Ri[v][0][1] = ( bi01 * cph - br01 * sph ) * x;
      Ri[v][1][0] =   bi10 * cph + br10 * sph;
      Ri[v][1][1] =   bi11 * cph + br11 * sph;
      Rr[v][1][0] =   br10 * cph - bi10 * sph;
      Rr[v][1][1] =   br11 * cph - bi11 * sph;
    }

  Mr[0][0] = br00;  
  Mr[0][1] = br01; 
  Mr[1][0] = br10; 
  Mr[1][1] = br11; 
  Mi[0][0] = bi00; 
  Mi[0][1] = bi01; 
  Mi[1][0] = bi10; 
  Mi[1][1] = bi11; 

  (*fatorM) = exp( -2 * fatM ); // to avoid overflow

  // ------------------------------------------------- adjust matrices L ----------

  for ( v = 0; v <= final-1; v ++ )
    {      
      a= v+2;
      beta= 0.25 / ( n[a] * n[a] + k[a] * k[a] );
      Lr[v][0][0] = ( n[a] * Lr[v][0][0] - k[a] * Li[v][0][0] ) * beta;
      Lr[v][0][1] = ( n[a] * Lr[v][0][1] - k[a] * Li[v][0][1] ) * beta;
      Lr[v][1][0] = ( n[a] * Lr[v][1][0] - k[a] * Li[v][1][0] ) * beta;
      Lr[v][1][1] = ( n[a] * Lr[v][1][1] - k[a] * Li[v][1][1] ) * beta;
      Li[v][0][0] = ( n[a] * Li[v][0][0] + k[a] * Lr[v][0][0] ) * beta;
      Li[v][0][1] = ( n[a] * Li[v][0][1] + k[a] * Lr[v][0][1] ) * beta;
      Li[v][1][0] = ( n[a] * Li[v][1][0] + k[a] * Lr[v][1][0] ) * beta;
      Li[v][1][1] = ( n[a] * Li[v][1][1] + k[a] * Lr[v][1][1] ) * beta;
    }

  // ----------------------------------------------- partial derivatives ----------

  for ( v = 1; v <= final; v ++ )
    {
      // simple common constant values

      nv      = n[v];
      Kv      = k[v];
      beta    = 2 * PI * d[v-1] / lambda;
      fat     = beta * Kv; 
      x       = exp( -2 * fat );
      phi     = beta * nv;
      sph     = sin( phi );
      cph     = cos( phi );
      nvp     = n[v+1];
      nva     = n[v-1];
      Kvp     = k[v+1];
      Kva     = k[v-1];

      // compounded common constant values

      Kv2           = Kv * Kv;
      Kv3           = Kv2 * Kv;
      nv2           = nv * nv;
      KvaKvp        = Kva * Kvp;
      Kv2Kvp        = Kv2 * Kvp;
      nv2PKv2       = Kv2 + nv2;
      nvaPnv        = nva + nv;
      nvaLnv        = nva - nv; 
      nvPnvp        = nv  + nvp; 
      nvLnvp        = nv  - nvp; 
      nvaP2nv       = nva + 2 * nv;
      nvaL2nv       = nva - 2 * nv;
      c2nvPnvp      = 2 * nv + nvp;
      c2nvLnvp      = 2 * nv - nvp;
      c2Kvnv        = 2 * Kv * nv;
      nv2Pnvanvp    = nv2 + nva * nvp;
      nv2Lnvanvp    = nv2 - nva * nvp;
      c37 = Kvp * nva + Kva * nvp;
      c38 = 3 * Kv2 + 2 * Kv * Kvp;
      c39 = 3 * Kv2 - 2 * Kv * Kvp;
      c40 = Kva * (2 * Kv + Kvp);
      c41 = Kva * (2 * Kv - Kvp);
      c42 = nv * (3 * nv + 2 * nvp);
      c43 = nv * (3 * nv - 2 * nvp);
      c09 = Kv * nv2Pnvanvp;
      c10 = Kv * nv2Lnvanvp;
      c11 = Kva * (Kv * (Kv + Kvp) + nv * nvPnvp);
      c12 = Kva * (Kv * (Kv - Kvp) + nv * nvLnvp);
      c13 = Kvp * nv * nvaPnv;
      c14 = Kvp * nv * nvaLnv;
      c15 = Kv3 + Kv2Kvp;
      c16 = Kv3 - Kv2Kvp;
      c17 = nv * (KvaKvp - nvaPnv * nvPnvp);
      c18 = nv * (KvaKvp + nvaLnv * nvLnvp);
      c19 = nv * (KvaKvp - nvaLnv * nvPnvp);
      c20 = nv * (KvaKvp + nvaPnv * nvLnvp);
      c21 = Kv * c37;
      c22 = 2 * Kv * (nvaPnv + nvp);
      c23 = 2 * Kv * (nvaPnv - nvp);
      c24 = 2 * Kv * (nvaLnv + nvp);
      c25 = 2 * Kv * (nvaLnv - nvp);
      c26 = Kv2 * (nvaPnv + nvp);
      c27 = Kv2 * (nvaPnv - nvp);
      c28 = Kv2 * (nvaLnv + nvp);
      c29 = Kv2 * (nvaLnv - nvp);
      c30 = Kvp * nvaP2nv;
      c31 = Kvp * nvaL2nv;
      c32 = Kv2 - KvaKvp;
      c33 = Kv2 + KvaKvp;
      c34 = nva * c2nvPnvp;
      c35 = nva * c2nvLnvp;
      c36 = Kva * c2nvPnvp;
      c36 = Kva * c2nvLnvp;

      // Real part of X

      rSp[0][0] = -c15 - c13 - c10 - c11; rSp[0][1] = -c15 + c14 - c09 + c11;
      rSp[1][0] =  c16 - c13 + c09 + c12; rSp[1][1] =  c16 + c14 + c10 - c12;
      rSq[0][0] = -c16 - c14 - c10 + c12; rSq[0][1] = -c16 + c13 - c09 - c12;
      rSq[1][0] =  c15 - c14 + c09 - c11; rSq[1][1] =  c15 + c13 + c10 + c11;
      rCp[0][0] =  c26 + c21 - c17; rCp[0][1] = -c29 - c21 + c19;
      rCp[1][0] = -c27 + c21 - c20; rCp[1][1] =  c28 - c21 + c18;
      rCq[0][0] = rCp[1][1]; rCq[0][1] = rCp[1][0];
      rCq[1][0] = rCp[0][1]; rCq[1][1] = rCp[0][0];

      // Imaginary part of X

      iSp[0][1] =  c29 + c21 - c19; iSp[1][0] =  c27 - c21 + c20;
      iSp[0][0] = -c26 - c21 + c17; iSp[1][1] = -c28 + c21 - c18;
      iSq[0][0] = rCp[1][1]; iSq[0][1] = rCp[1][0];
      iSq[1][1] = rCp[0][0]; iSq[1][0] = rCp[0][1];
      iCp[0][0] = rSp[0][0]; iCp[0][1] = rSp[0][1];
      iCp[1][0] = rSp[1][0]; iCp[1][1] = rSp[1][1];
      iCq[0][0] = rSp[1][1]; iCq[0][1] = rSp[1][0];
      iCq[1][0] = rSp[0][1]; iCq[1][1] = rSp[0][0];  

      // Assemble X (with its real and imaginary parts)

      for ( l = 0; l < 2; l ++ )
	for ( c = 0; c < 2; c ++ )
	  {
	    Xr[l][c] = x  * ( rSp[l][c] * sph + rCp[l][c] * cph ) + 
	      rSq[l][c] * sph + rCq[l][c] * cph;

	    Xi[l][c] = x  * ( iSp[l][c] * sph + iCp[l][c] * cph ) + 
	      iSq[l][c] * sph + iCq[l][c] * cph;
	  }

      // ------------------------------- derivatives of n (refractive index) ----------
   
      drSp[1][0] =  c2Kvnv - c30 + c36; drSq[1][1] =  c2Kvnv + c30 + c36;
      drSp[1][1] =  c2Kvnv + c31 - c36; drSq[1][0] =  c2Kvnv - c31 - c36;
      diCq[1][0] = -c2Kvnv + c31 + c36; drSp[0][0] = -c2Kvnv - c30 - c36;
      drSq[0][1] = -c2Kvnv + c30 - c36; drSq[0][0] = -c2Kvnv - c31 + c36;
      drCp[0][0] =  c32 + c34 + c42; drCp[0][1] =  c33 - c34 + c42;
      diSp[1][0] =  c33 + c35 + c43; diSp[1][1] =  c32 - c35 + c43;
      drCp[1][1] = -c32 + c35 - c43; drCp[1][0] = -c33 - c35 - c43;
      diSp[0][1] = -c33 + c34 - c42; diSp[0][0] = -c32 - c34 - c42;
      diCp[0][1] = diCq[1][0]; diCp[1][0] = drSp[1][0];
      diCq[0][1] = drSp[1][0]; diCp[1][1] = drSp[1][1];
      diCq[0][0] = drSp[1][1]; diCp[0][0] = drSp[0][0];
      diCq[1][1] = drSp[0][0]; drSp[0][1] = diCq[1][0];
      drCq[0][0] = drCp[1][1]; drCq[1][0] = drCp[0][1];
      drCq[1][1] = drCp[0][0]; diSq[1][0] = drCp[0][1];
      diSq[1][1] = drCp[0][0]; diSq[0][0] = drCp[1][1];
      drCq[0][1] = drCp[1][0]; diSq[0][1] = drCp[1][0];
	  
      for ( l = 0; l < 2; l ++ )
	for ( c = 0; c < 2; c ++ )
	  {
	    dXr[l][c] = beta *  ( cph * ( x *  rSp[l][c] +  rSq[l][c] ) - 
				  sph * ( x *  rCp[l][c] +  rCq[l][c] ) ) + 
	      sph * ( x * drSp[l][c] + drSq[l][c] ) + 
	      cph * ( x * drCp[l][c] + drCq[l][c] );
	  
	    dXi[l][c] = beta *  ( cph * ( x *  iSp[l][c] +  iSq[l][c] ) - 
				  sph * ( x *  iCp[l][c] +  iCq[l][c] ) ) + 
	      sph * ( x * diSp[l][c] + diSq[l][c] ) + 
	      cph * ( x * diCp[l][c] + diCq[l][c] );
   	  
	    dXr[l][c] = ( dXr[l][c]  *  nv2PKv2 - 2 * nv * Xr[l][c] ) / 
	      ( nv2PKv2 * nv2PKv2 );

	    dXi[l][c] = ( dXi[l][c]  *  nv2PKv2 - 2 * nv * Xi[l][c] ) / 
	      ( nv2PKv2 * nv2PKv2 );

	    dXrdn[v-1][l][c] = dXr[l][c];
	    dXidn[v-1][l][c] = dXi[l][c];
	  }

      // ------------------------- derivatives of k (absorption coefficient) ----------

      drSp[1][0] =  c39 + c41 + nv2Pnvanvp; drSp[1][1] =  c39 - c41 + nv2Lnvanvp;
      drSq[1][0] =  c38 - c40 + nv2Pnvanvp; drSq[1][1] =  c38 + c40 + nv2Lnvanvp;
      drSp[0][0] = -c38 - c40 - nv2Lnvanvp; drSp[0][1] = -c38 + c40 - nv2Pnvanvp;
      drSq[0][1] = -c39 - c41 - nv2Pnvanvp; drSq[0][0] = -c39 + c41 - nv2Lnvanvp;
      drCp[0][0] =  c37 + c22; diSp[1][1] =  c37 - c24;
      drCp[1][0] =  c37 - c23; diSp[0][1] =  c37 + c25;
      drCp[0][1] = -c37 - c25; drCp[1][1] = -c37 + c24;
      diSp[0][0] = -c37 - c22; diSp[1][0] = -c37 + c23;
      diSq[0][0] = drCp[1][1]; diSq[1][0] = drCp[0][1];
      drCq[1][1] = drCp[0][0]; diSq[1][1] = drCp[0][0];
      drCq[0][1] = drCp[1][0]; diSq[0][1] = drCp[1][0];
      drCq[0][0] = drCp[1][1]; drCq[1][0] = drCp[0][1];
      diCp[1][0] = drSp[1][0]; diCq[0][1] = drSp[1][0];
      diCp[1][1] = drSp[1][1]; diCq[0][0] = drSp[1][1];
      diCp[0][0] = drSp[0][0]; diCp[0][1] = drSp[0][1];
      diCq[1][0] = drSp[0][1]; diCq[1][1] = drSp[0][0];

      for ( l = 0; l < 2; l ++ )
	for ( c = 0; c < 2; c ++ )
	  {
	    dXr[l][c] = x * ( sph * ( drSp[l][c] - beta * rSp[l][c] ) + 
			      cph * ( drCp[l][c] - beta * rCp[l][c] ) ) + 
	      sph * ( drSq[l][c] + beta * rSq[l][c] ) + 
	      cph * ( drCq[l][c] + beta * rCq[l][c] );
		
	    dXi[l][c] = x * ( sph * ( diSp[l][c] - beta * iSp[l][c] ) + 
			      cph * ( diCp[l][c] - beta * iCp[l][c] ) ) + 
	      sph * ( diSq[l][c] + beta * iSq[l][c] ) + 
	      cph * ( diCq[l][c] + beta * iCq[l][c] );
		
	    dXr[l][c] = ( dXr[l][c]  *  nv2PKv2 - 2 * Kv * Xr[l][c] ) / 
	      ( nv2PKv2 * nv2PKv2 );
		
	    dXi[l][c] = ( dXi[l][c]  *  nv2PKv2 - 2 * Kv * Xi[l][c] ) / 
	      ( nv2PKv2 * nv2PKv2 );

	    dXrdk[v-1][l][c] = dXr[l][c];
	    dXidk[v-1][l][c] = dXi[l][c];
	  }

    } // films :: v ::

} // routine

//*******************************************************************
//*******************************************************************

void compdtr( int nlayers, int slayer, int nobs, 
	      double n[][nobsmax], double k[][nobsmax], 
	      double d[], double lambda[], double t[], double r[],
	      double dtdn[][nobsmax], double dtdk[][nobsmax], 
	      double drdn[][nobsmax], double drdk[][nobsmax] )
{
  int i, j, v;
  double nl[nlayersmax], kl[nlayersmax], 
    tn[nlayersmax], tk[nlayersmax], rn[nlayersmax], rk[nlayersmax],
    fator[2], h, h2, dh, tmp, r1, r2, r2e2, r3, r4, dr1, dr2, dr3, dr4,
    dmB01, dmB10, dmB11, dmB00, dmC01, dmC10, dmC11, dmC00, dB[nlayersmax], dC[nlayersmax],
    Br[2][2], Bi[2][2], Cr[2][2], Ci[2][2], mB[2][2], mC[2][2],
    dBr[2][2], dBi[2][2], dCr[2][2], dCi[2][2],
    dXrdn1[nlayersmax-1][2][2], dXidn1[nlayersmax-1][2][2], 
    dXrdk1[nlayersmax-1][2][2], dXidk1[nlayersmax-1][2][2], 
    dXrdn2[nlayersmax-1][2][2], dXidn2[nlayersmax-1][2][2], 
    dXrdk2[nlayersmax-1][2][2], dXidk2[nlayersmax-1][2][2], 
    R1r[nlayersmax-1][2][2], R1i[nlayersmax-1][2][2],
    L1r[nlayersmax-1][2][2], L1i[nlayersmax-1][2][2],
    R2r[nlayersmax-1][2][2], R2i[nlayersmax-1][2][2],
    L2r[nlayersmax-1][2][2], L2i[nlayersmax-1][2][2];

  //d[1]= d[0];

  for ( i= 0; i <= slayer-2; i ++ ) // set front films thicknesses
    dB[i] = d[i];

  for ( i= slayer+1; i <= nlayers-2; i ++ ) // set back films thicknesses
    dC[i-slayer-1] = d[i-2];

  for ( j = 0; j < nobs; j ++ )
    {
      /********** MATRIX B (FRONTAL FILMS) *******************************************/

      for ( i= 0; i <= slayer; i ++ ) // set 'n' and 'kappa'
	{
	  nl[i]  = n[i][j];
	  kl[i]  = k[i][j];
	}

      if (slayer > 1) // number of front films > 0
	tance( slayer-1, nobs, nl, kl, dB, lambda[j], &fator[1], 
	       Br, Bi, L1r, L1i, R1r, R1i, dXrdn1, dXidn1, dXrdk1, dXidk1 );
      else // no front films
	{
	  fator[1] = 1;
	  Br[0][0] = ( n[slayer][j] + 1 ) / ( 2 * n[slayer][j] );
	  Br[0][1] = ( n[slayer][j] - 1 ) / ( 2 * n[slayer][j] );
	  Br[1][0] = Br[0][1];
	  Br[1][1] = Br[0][0];
	  Bi[0][0] = 0; Bi[0][1] = 0; Bi[1][0] = 0; Bi[1][1] = 0; 
	}
      
      mB[0][0] = Br[0][0] * Br[0][0] + Bi[0][0] * Bi[0][0]; // || B_{11} ||^2
      mB[0][1] = Br[0][1] * Br[0][1] + Bi[0][1] * Bi[0][1]; // || B_{12} ||^2
      mB[1][0] = Br[1][0] * Br[1][0] + Bi[1][0] * Bi[1][0]; // || B_{21} ||^2
      mB[1][1] = Br[1][1] * Br[1][1] + Bi[1][1] * Bi[1][1]; // || B_{22} ||^2

      /********** MATRIX C (BACK FILMS) **********************************************/

      for ( i= slayer; i <= nlayers-1; i ++ ) // set 'n' and 'kappa'
	{
	  nl[i-slayer] = n[i][j];
	  kl[i-slayer] = k[i][j];
	}

      // filmes posteriores

      if (nlayers - slayer > 2) // number of back films > 0
	tance( nlayers-slayer-2, nobs, nl, kl, dC, lambda[j], &fator[2],
	       Cr, Ci, L2r, L2i, R2r, R2i, dXrdn2, dXidn2, dXrdk2, dXidk2 );
      else // no back films
	{
	  fator[2] = 1;
	  Cr[0][0] = ( 1 + n[slayer][j] ) / 2;
	  Cr[0][1] = ( 1 - n[slayer][j] ) / 2;
	  Cr[1][0] = Cr[0][1];
	  Cr[1][1] = Cr[0][0];
	  Ci[0][0] = 0; Ci[0][1] = 0; Ci[1][0] = 0; Ci[1][1] = 0; 
	}

      mC[0][0] = Cr[0][0] * Cr[0][0] + Ci[0][0] * Ci[0][0]; // || C_{11} ||^2
      mC[0][1] = Cr[0][1] * Cr[0][1] + Ci[0][1] * Ci[0][1]; // || C_{12} ||^2
      mC[1][0] = Cr[1][0] * Cr[1][0] + Ci[1][0] * Ci[1][0]; // || C_{21} ||^2
      mC[1][1] = Cr[1][1] * Cr[1][1] + Ci[1][1] * Ci[1][1]; // || C_{22} ||^2

      /********** TRANSMITTANCE & REFLECTANCE ****************************************/

      h  = mC[1][1] * mB[1][1] - mC[1][0] * mB[0][1];
      h2 = h * h; 
      
      // evaluate transmittance
      t[j] = fator[1] * fator[2] / h;

      // evaluate reflectance
      r1 = 2 * ( Br[0][0] * Br[1][1] - Bi[0][0] * Bi[1][1] ) * fator[1];
      r2 = n[slayer][j] * mB[1][1];
      r2e2= r2 * r2;
      r3 = r1 / r2;
      r4 = mC[1][1] * mB[1][0] + mC[1][0] * ( r3 - mB[0][0] );
      r[j] = r4 / h;

      /********** PARTIAL DERIVATIVES ************************************************/

      if (slayer > 1) // ------------------------------------ number of front films > 0
	for ( v = 1; v <= slayer-1; v ++ ) // derivatives looping
	  {
	    // ---------- partial derivatives of n -----------------------------------------
	    assemble( L1r[v-1],L1i[v-1],R1r[v-1],R1i[v-1],dXrdn1[v-1],dXidn1[v-1], dBr,dBi );
	    dmB00 = 2 * ( Br[0][0] * dBr[0][0] + Bi[0][0] * dBi[0][0] );
	    dmB01 = 2 * ( Br[0][1] * dBr[0][1] + Bi[0][1] * dBi[0][1] );
	    dmB10 = 2 * ( Br[1][0] * dBr[1][0] + Bi[1][0] * dBi[1][0] );
	    dmB11 = 2 * ( Br[1][1] * dBr[1][1] + Bi[1][1] * dBi[1][1] );
	    
	    // DTDN
	    dh = mC[1][1] * dmB11 - mC[1][0] * dmB01;
    	    dtdn[v-1][j] = -dh * fator[1] * fator[2] / h2; 

	    // DRDN
	    dr1 = 2 * fator[1] * ( dBr[0][0] * Br[1][1] + Br[0][0] * dBr[1][1] 
				 - dBi[0][0] * Bi[1][1] - Bi[0][0] * dBi[1][1] );
	    dr2 = n[slayer][j] * dmB11;
	    dr3 = dr1 / r2 - r1 * dr2 / r2e2;
	    dr4 = mC[1][1] * dmB10 + mC[1][0] * ( dr3 - dmB00 );
    	    drdn[v-1][j] = dr4 / h - r4 * dh / h2;
	    
	    // ---------- partial derivatives of kappa -------------------------------------
	    assemble( L1r[v-1],L1i[v-1],R1r[v-1],R1i[v-1],dXrdk1[v-1],dXidk1[v-1], dBr,dBi );
	    dmB00 = 2 * ( Br[0][0] * dBr[0][0] + Bi[0][0] * dBi[0][0] );
	    dmB01 = 2 * ( Br[0][1] * dBr[0][1] + Bi[0][1] * dBi[0][1] );
	    dmB10 = 2 * ( Br[1][0] * dBr[1][0] + Bi[1][0] * dBi[1][0] );
	    dmB11 = 2 * ( Br[1][1] * dBr[1][1] + Bi[1][1] * dBi[1][1] );

	    // DTDK
	    dh = mC[1][1] * dmB11 - mC[1][0] * dmB01;
    	    dtdk[v-1][j] = -dh * fator[1] * fator[2] / h2; 

	    // DRDK
	    dr1 = 2 * fator[1] * ( dBr[0][0] * Br[1][1] + Br[0][0] * dBr[1][1] 
				 - dBi[0][0] * Bi[1][1] - Bi[0][0] * dBi[1][1] );
	    dr2 = n[slayer][j] * dmB11;
	    dr3 = dr1 / r2 - r1 * dr2 / r2e2;
	    dr4 = mC[1][1] * dmB10 + mC[1][0] * ( dr3 - dmB00 );
    	    drdk[v-1][j] = dr4 / h - r4 * dh / h2;
	  }

      if (nlayers - slayer > 2) // ---------------------------- number of back films > 0
	for ( v = 1; v <= nlayers-slayer-2; v ++ ) // derivatives looping
	  {
	    // ---------- partial derivatives of n -----------------------------------------
	    assemble( L2r[v-1],L2i[v-1],R2r[v-1],R2i[v-1],dXrdn2[v-1],dXidn2[v-1], dCr,dCi );
	    dmC00 = 2 * ( Cr[0][0] * dCr[0][0] + Ci[0][0] * dCi[0][0] );
	    dmC01 = 2 * ( Cr[0][1] * dCr[0][1] + Ci[0][1] * dCi[0][1] );
	    dmC10 = 2 * ( Cr[1][0] * dCr[1][0] + Ci[1][0] * dCi[1][0] );
	    dmC11 = 2 * ( Cr[1][1] * dCr[1][1] + Ci[1][1] * dCi[1][1] );
	    
	    // DTDN
	    dh = mB[1][1] * dmC11 - mB[0][1] * dmC10;
    	    dtdn[slayer+v-2][j] = -dh * fator[1] * fator[2] / h2; 

	    // DRDN
	    dr4 = dmC11 * mB[1][0] + dmC10 * ( r3 - mB[0][0] );
    	    drdn[slayer+v-2][j] = dr4 / h - r4 * dh / h2;
	    
	    // ---------- partial derivatives of kappa -------------------------------------
	    assemble( L2r[v-1],L2i[v-1],R2r[v-1],R2i[v-1],dXrdk2[v-1],dXidk2[v-1], dCr,dCi );
	    dmC00 = 2 * ( Cr[0][0] * dCr[0][0] + Ci[0][0] * dCi[0][0] );
	    dmC01 = 2 * ( Cr[0][1] * dCr[0][1] + Ci[0][1] * dCi[0][1] );
	    dmC10 = 2 * ( Cr[1][0] * dCr[1][0] + Ci[1][0] * dCi[1][0] );
	    dmC11 = 2 * ( Cr[1][1] * dCr[1][1] + Ci[1][1] * dCi[1][1] );

	    // DTDK
	    dh = mB[1][1] * dmC11 - mB[0][1] * dmC10;
    	    dtdk[slayer+v-2][j] = -dh * fator[1] * fator[2] / h2; 

	    // DRDK
	    dr4 = dmC11 * mB[1][0] + dmC10 * ( r3 - mB[0][0] );
    	    drdk[slayer+v-2][j] = dr4 / h - r4 * dh / h2;
	  }
      
    } // lambda ::: j :::

}
