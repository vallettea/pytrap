#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>
#include "conformal_map.h"

static struct {
    int potentialonly;

    int poincare;
    size_t nbpmax;
    double t1;

    double dt;

    double r0;
    double E0;
    double L0;

    double V_ex;
    double omega_ex;

    double rmax;
} par;

static double V1, V2, V3, V4, V5;

#define NPOINTS	14
static double zi[NPOINTS] = {
    -0.1997, -0.1917, -0.187, -0.183, -0.1804, -0.1735, -0.1705, -0.165,
    -0.1619, -0.1554, -0.1525,
    -0.14655, -0.134, -0.1257
};
static double Vi[NPOINTS];
static struct trap *trap;
static struct trap *excit;


#define DIM 6
static int
func(double t, const double y[], double f[], void *params)
{
    struct trap *trap = (struct trap *) params;
    const double phi2z = trap_get_pot2(trap, y[2]) / 2.0;

    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    if (par.V_ex == 0.0) {
	f[3] = y[0] * phi2z;
	f[4] = y[1] * phi2z;
	f[5] = -trap_get_pot1(trap, y[2])
	    + (y[0]*y[0] + y[1]*y[1])/4 * trap_get_pot3(trap, y[2]);
    }
    else {
	const double M = sin(par.omega_ex * t);
	const double phi2z_tot = phi2z + M * trap_get_pot2(excit, y[2]) / 2.0;
	f[3] = y[0] * phi2z_tot;
	f[4] = y[1] * phi2z_tot;
	f[5] = -(trap_get_pot1(trap, y[2]) + M*trap_get_pot1(excit, y[2]))
	    + (y[0]*y[0] + y[1]*y[1])/4 * (trap_get_pot3(trap, y[2])
		    + M*trap_get_pot3(excit, y[2]));
    }
    return GSL_SUCCESS;
}
//  r = sqrt(x*x + y*y)
//  dotr = (x*px + y*py) / sqrt(x*x + y*y)
//  L = x*py - y*px

static inline double
energy(const double y[])
{
    return (y[3]*y[3] + y[4]*y[4] + y[5]*y[5]) / 2
	+ trap_get_pot(trap, y[2])
	- (y[0]*y[0] + y[1]*y[1])/4 * trap_get_pot2(trap, y[2]);
}
static inline double
kinmom(const double y[])
{
    return y[0]*y[4] - y[1]*y[3];
}

static void
usage(char **argv)
{
    fprintf(stderr, "Usage: %s [options] V1 V2 V3 V4 Vz E0 L0 r0\n"
	    "\n"
	    "Options:\n"
	    "  -P            Print the potential and its derivatives then exits\n"
	    "  -n nb_points  Select Poincare section mode\n"
	    "  -t time       Select trajectory mode: solve the ODE until time\n"
	    "  -d timestep   Specify the time increment for the integrator\n"
	    "  -V V_ex       Amplitude of excitation (default is 0.0)\n"
	    "  -W omega_ex   Pulsation of excitation (default is 0.0)\n"
	    "\n"
	    ,
	    argv[0]);
    exit(1);
}

int
main(int argc, char **argv)
{
    /* DEFAULT VALUES */
    par.potentialonly = 0;
    par.dt = 1e-4;
    par.poincare = 1;
    par.nbpmax = 500;
    par.t1 = 1e2;
    par.r0 = 0.0;
    par.L0 = 0.0;
    par.E0 = 0.0;

    par.V_ex = 0.0;		/* excitation */
    par.omega_ex = 0.0;

    par.rmax = 1.0;		/* cutoff */

    /* PARSE */
    {
	int c;
	extern char *optarg;
	extern int optind;
	while((c = getopt(argc, argv, "Pd:n:t:V:W:")) != -1) {
	    switch(c) {
		case 'P':
		    par.potentialonly = 1;
		    break;
		case 'd':
		    par.dt = atof(optarg);
		    break;
		case 'n':
		    par.nbpmax = atoi(optarg);
		    par.poincare = 1;
		    par.t1 = 1e99;
		    break;
		case 't':
		    par.t1 = atof(optarg);
		    par.poincare = 0;
		    break;
		case 'V':
		    par.V_ex = atof(optarg);
		    break;
		case 'W':
		    par.omega_ex = atof(optarg);
		    break;
		case ':':
		case '?':
		    usage(argv);
		default:
		    fprintf(stderr, "Unknown error\n");
	    }
	}
	if (optind + (par.potentialonly?5:8) > argc)
	    usage(argv);

	V1 = atof(argv[optind++]);
	V2 = atof(argv[optind++]);
	V3 = atof(argv[optind++]);
	V4 = atof(argv[optind++]);
	V5 = atof(argv[optind++]);
	if (!par.potentialonly) {
	    par.E0 = atof(argv[optind++]);
	    par.L0 = atof(argv[optind++]);
	    par.r0 = atof(argv[optind++]);
	}
    }

    /* INIT TRAP */
    {
	size_t k;
	for(k=0; k<NPOINTS; k++) {
	    zi[k] += 0.1525;
	}

	trap = trap_alloc(zi, NPOINTS);
	trap->n = 11;
	trap->R = 0.008;
	trap->Rz = 0.013;
	trap_init(trap);

	Vi[0] = 0.0;
	Vi[1] = Vi[2] = V1;
	Vi[3] = Vi[4] = V2;
	Vi[5] = Vi[6] = V3;
	Vi[7] = Vi[8] = V4;
	Vi[9] = Vi[10] = 0.0;
	Vi[11] = Vi[12] = V5;
	Vi[13] = 0.0;
	assert(NPOINTS == 14);
	trap_set_pot(trap, Vi);


	excit = trap_alloc(zi, NPOINTS);
	excit->n = trap->n;
	excit->R = trap->R;
	excit->Rz = trap->Rz;
	trap_init(excit);
	{
	    double Vex[NPOINTS];
	    for(k=0; k<NPOINTS; k++) {
		Vex[k] = 0.0;
	    }
	    Vex[9] = Vex[10] = par.V_ex;
	    trap_set_pot(excit, Vex);
	}
    }

    /* the -P switch */
    if (par.potentialonly) {
	double z;
	printf("# static potential\n");
	for(z=0.0; z<0.2111; z+=1e-4) {
	    printf("%e %e %e %e %e\n", z,
		    trap_get_pot(trap, z),
		    trap_get_pot1(trap, z),
		    trap_get_pot2(trap, z),
		    trap_get_pot3(trap, z)
		  );
	}
	printf("\n");

	printf("# excitation potential\n");
	for(z=0.0; z<0.2111; z+=1e-4) {
	    printf("%e %e %e %e %e\n", z,
		    trap_get_pot(excit, z),
		    trap_get_pot1(excit, z),
		    trap_get_pot2(excit, z),
		    trap_get_pot3(excit, z)
		  );
	}
	exit(0);
    }

    /* rkck is good for energy conservation */
#define GSLODEIVTYPE	gsl_odeiv_step_rkck
#define GSLODEIVEPSREL	1e-8

    gsl_odeiv_step *s = gsl_odeiv_step_alloc(GSLODEIVTYPE, DIM);
    gsl_odeiv_system sys = {&func, NULL, DIM, trap};
    gsl_odeiv_control *c;
    gsl_odeiv_evolve *e;
    double t = 0.0, tprev, dt = par.dt;
    double y[DIM], yprev[DIM];
    int status;
    size_t iter = 0, npp = 0;

    y[0] = par.r0;		// x_i
    y[1] = 0.0;			// y_i
    y[2] = 0.0;			// z_i

    y[3] = 0.0;						// dot x_i
    y[4] = (par.L0 == 0.0) ? 0.0 : par.L0/par.r0;	// dot y_i
    y[5] = sqrt(2*par.E0 - y[3]*y[3] - y[4]*y[4]);	// dot z_i

    printf("\n");
    printf("#V[i]: %g %g %g %g %g\n", V1, V2, V3, V4, V5);
    printf("#excitation: %g %g\n", par.V_ex, par.omega_ex);
    printf("#E0,L0,r0: %g %g %g\n", par.E0, par.L0, par.r0);
    printf("#y0: %e %e %e %e\n", y[0], y[1], y[2], y[3]);
    printf("# 1:t 2:x 3:y 4:z 5:px 6:py 7:pz\n");

    c = gsl_odeiv_control_y_new(GSLODEIVEPSREL, 0.0);
    e = gsl_odeiv_evolve_alloc(DIM);

    while (t < par.t1)
    {
	tprev = t;
	memcpy(yprev, y, sizeof(y));

	/* ONE STEP */
	status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, par.t1, &dt, y);
	if (status != GSL_SUCCESS) {
	    fprintf(stderr, "GSL failure. Exiting\n");
	    break;
	}
	if (fabs(y[1]) >= 0.013 || fabs(y[2]) >= 0.24 || gsl_isnan(y[2])) {
	    fprintf(stderr, "Trajectory goes out of bound. Exiting.\n");
	    break;
	}
	if (dt > par.dt) dt = par.dt;

	/* OUTPUT */
	if (par.poincare) {
	    if (yprev[2]*y[2] <= 0.0) {
		double h0 = -yprev[2] / (y[2] - yprev[2]);
		printf("%e %e %e %e %e %e %e\n",
			tprev + h0*dt,
			yprev[0] + h0*(y[0]-yprev[0]),
			yprev[1] + h0*(y[1]-yprev[1]),
			yprev[2] + h0*(y[2]-yprev[2]),
			yprev[3] + h0*(y[3]-yprev[3]),
			yprev[4] + h0*(y[4]-yprev[4]),
			yprev[5] + h0*(y[5]-yprev[5])
		      );
		fflush(stdout);
		npp++;
		if (npp % 10 == 0) {
		    fprintf(stderr, "npp:%zu t:%.3e iter:%zuk dt:%.2e\n",
			    npp, t, iter/1000, dt);
		}
		if (npp >= par.nbpmax) {
		    fprintf(stderr, "Enough points on section. Exiting.\n");
		    break;
		}
	    }
	}
	else {
	    printf("%e %e %e %e %e %e %e\n",
		    t, y[0], y[1], y[2], y[3], y[4], y[5]);
	}
	if (iter % 10000 == 0) {
	    if (par.V_ex == 0.0) {
		double dE = (energy(y) - par.E0) / par.E0;
		fprintf(stderr, "t:%e iter:%zuk dt:%.2e dL0:%+.2e dE0:%+.2e\n",
			t, iter/1000, dt, kinmom(y)-par.L0, dE);
		if (fabs(dE) > 1e-3) {
		    fprintf(stderr, "dE is now too high. Exiting.\n");
		    break;
		}
	    }
	    if (isnan(y[0]) || isnan(y[1]) || hypot(y[0],y[1]) > par.rmax) {
		fprintf(stderr, "Diverging (x:%g, y:%g). Exiting.\n", y[0], y[1]);
		break;
	    }
	}
	iter++;
    }
    fprintf(stderr, "END t:%e (%zu iters) r:%e z:%e\n",
	    t, iter, hypot(y[0],y[1]), y[2]);

    /* */
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    trap_free(trap);
    return 0;
}
