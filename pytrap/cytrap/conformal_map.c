#include <string.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_poly.h>
#include "conformal_map.h"

#define OMEGA	(1.3152/M_PI)

static const double ai[] = {
    -0.6646800942825573, 308.41773377707017, -3671.0848933224097,
    1.2975843811561745e6, -7.4587154945175685e6,
    -2.0357109763813026e10, 7.87621758339513e10, 1.0895422241161766e14
};
static const double z1 = -0.007;
static const double z2 = +0.007;

double complex
f2inva(const struct trap *trap, const double z)
{
    if (z < z1) 
	return -clog(trap->abp
		+ trap->abm * cexp(-trap->sqrtb * (M_PI*z/trap->Rz - trap->aa))
		);
    if (z > z2)
	return clog(trap->ac * (trap->ad
		    + cexp(M_PI*z/trap->Rz + trap->aa / trap->sqrtb)
		    ));
    return gsl_poly_eval(ai, sizeof(ai)/sizeof(ai[0]), z);
}

double complex
df2inva(const struct trap *trap, const double z)
{
    if (z < z1) {
	double Z = M_PI / trap->Rz;
	complex double K = cexp(-trap->sqrtb * (Z*z - trap->aa));
	return (+trap->abm * trap->sqrtb * Z * K)
	    / (trap->abp + trap->abm * K);
    }
    if (z > z2) {
	double Z = M_PI / trap->Rz;
	complex double K = cexp(Z*z + trap->aa / trap->sqrtb);
	return (Z*K) / (trap->ad + K);
    }

    double res[2];
    gsl_poly_eval_derivs(ai, sizeof(ai)/sizeof(ai[0]), z, res, 2);
    return res[1];
}


static double complex
map_W_to_Z(const struct trap *trap, const double complex w)
{
    const double b = gsl_pow_2(trap->Rz / trap->R);
    const double complex t = cexp(w);
    const double complex sqrt_tm1 = csqrt(t-1);
    const double complex sqrt_tmb = csqrt(t-b);
    const double complex sqrt_btmb = trap->Rz / trap->R * sqrt_tm1;

    return (
	  trap->Rz * clog((sqrt_tmb + sqrt_tm1) / (sqrt_tmb - sqrt_tm1))
	+ trap->R  * clog((sqrt_btmb - sqrt_tmb) / (sqrt_btmb + sqrt_tmb))
    ) / M_PI;
}


struct params {
    const struct trap *trap;
    complex double z0;
};

static int
func_z(const gsl_vector *w, void *ptr, gsl_vector *z)
{
    const struct params *params = (const struct params *) ptr;
    complex double zz = map_W_to_Z(params->trap,
	    gsl_vector_get(w,0) + I*gsl_vector_get(w,1));
    gsl_vector_set(z, 0, creal(zz) - creal(params->z0));
    gsl_vector_set(z, 1, cimag(zz) - cimag(params->z0));
    return GSL_SUCCESS;
}

static complex double 
map_Z_to_W(const struct trap *trap, const double complex z)
{
    gsl_multiroot_function f; 
    gsl_multiroot_fsolver *s;
    struct params params;
    gsl_vector *w; 
    complex double ww;
    size_t iter = 0;
    int status;

    if (z == I*trap->R)
	return 0.0;

    params.trap = trap;
    params.z0 = z;
    f.f = &func_z;
    f.n = 2;
    f.params = (void *) & params;
    s = gsl_multiroot_fsolver_alloc (FSOLVER, 2);

    w = gsl_vector_alloc(2);
    gsl_vector_set(w, 0, (creal(z)<0)?-10.0:4.0);
    gsl_vector_set(w, 1, 2*M_PI);
    gsl_multiroot_fsolver_set(s, &f, w);

    do {
	iter++;
	status = gsl_multiroot_fsolver_iterate(s);
	if (status)
	    break;
	status = gsl_multiroot_test_residual(s->f, FSOLVER_EPSILON);
    }
    while (status == GSL_CONTINUE && iter < 1000);

    ww = gsl_vector_get(s->x, 0) + I*gsl_vector_get(s->x, 1);
    /*
    fprintf(stderr, "%s: (%e %e) after %zu iterations (%i: %s)\n",
	    __FUNCTION__, creal(ww), cimag(ww), iter, status,
	    gsl_strerror(status));
	    */

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(w);
    return ww;
}

struct trap *
trap_alloc(const double *zi, const size_t ntot)
{
    struct trap *trap = malloc(sizeof(*trap));
    if (!trap)
	goto errmem;

    trap->ntot = ntot;
    trap->n = DEFAULT_TRAP_N;
    trap->R = DEFAULT_TRAP_R;
    trap->Rz = DEFAULT_TRAP_RZ;

    trap->z = malloc(ntot * sizeof(* trap->z));
    if (!trap->z) goto errmem;
    memcpy(trap->z, zi, ntot*sizeof(*zi));

    trap->V = malloc(ntot * sizeof(* trap->V));
    if (!trap->V) goto errmem;

    trap->xiw = malloc(ntot * sizeof(* trap->xiw));
    if (!trap->xiw) goto errmem;
    trap->dxiw = malloc(ntot * sizeof(* trap->xiw));
    if (!trap->dxiw) goto errmem;
    trap->Q = malloc(ntot * sizeof(* trap->V));
    if (!trap->Q) goto errmem;

    trap->npi = (INTERP_ZMAX - INTERP_ZMIN) / INTERP_DZ;

    trap->zki = malloc(trap->npi*sizeof(double));
    if (!trap->zki) goto errmem;
    trap->ffzki = malloc(trap->npi*sizeof(double));
    if (!trap->ffzki) goto errmem;
    trap->dffzki = malloc(trap->npi*sizeof(double));
    if (!trap->dffzki) goto errmem;

    trap->Vzki = malloc(trap->npi*sizeof(double));
    if (!trap->Vzki) goto errmem;
    trap->V1zki = malloc(trap->npi*sizeof(double));
    if (!trap->V1zki) goto errmem;

    trap->acc = gsl_interp_accel_alloc();
    trap->spline = gsl_spline_alloc(gsl_interp_cspline, trap->npi);
    trap->acc1 = gsl_interp_accel_alloc();
    trap->spline1 = gsl_spline_alloc(gsl_interp_cspline, trap->npi);

    return trap;

errmem:
    perror(__FUNCTION__);
    exit(1);
}

int
trap_init(struct trap *trap)
{
    size_t k;

    trap->sqrtb = trap->Rz / trap->R;
    trap->b = trap->sqrtb * trap->sqrtb;
    trap->aa = log((trap->Rz + trap->R) / (trap->Rz - trap->R));
    trap->abm = (1 - trap->b) / (4*trap->b);
    trap->abp = (1 + trap->b) / (4*trap->b);
    trap->ac = (1 - trap->b) / 4;
    trap->ad = (1 + trap->b) / (1 - trap->b);


    for(k=0; k<trap->n; k++) {
	complex double xiw = map_Z_to_W(trap, trap->z[k] + I*trap->R);
	trap->xiw[k] = creal(xiw);
    }
    for(k=trap->n; k<trap->ntot; k++) {
	// border
	complex double zk;
	zk = map_Z_to_W(trap, trap->z[k] + I*trap->Rz);
	zk = map_W_to_Z(trap, creal(zk) - I*M_PI);
	zk = creal(zk);

	complex double xiw = map_Z_to_W(trap, zk + I*trap->Rz);
	trap->xiw[k] = creal(xiw);
    }

    for(k=0; k<trap->ntot-1; k++) {
	trap->dxiw[k] = trap->xiw[k+1] - trap->xiw[k];
    }

    for(k=0; k<trap->npi; k++) {
	trap->zki[k] = INTERP_ZMIN + k*INTERP_DZ;
	//trap->ffzki[k] = creal(map_Z_to_W(trap, 0.1525 - trap->zki[k]));
	trap->ffzki[k] = creal(f2inva(trap, 0.1525 - trap->zki[k]));
	trap->dffzki[k] = -creal(df2inva(trap, 0.1525 - trap->zki[k]));
	// signe moins parce qu'on derive avec (0.1525 - z)
    }

if (1)
{
    for(k=0; k<trap->ntot; k++) trap->xiw[k] += 30.0;
    for(k=0; k<trap->npi; k++) trap->ffzki[k] += 30.0;
}

    return 0;
}

int
trap_set_pot(struct trap *trap, const double *Vi)
{
    size_t i, k;

    memcpy(trap->V, Vi, trap->ntot*sizeof(*trap->V));
    for(i=0; i<trap->ntot-1; i++) {
	trap->Q[i] = (Vi[i+1] - Vi[i]) / (trap->xiw[i+1] - trap->xiw[i]);
    }

    for(k=0; k<trap->npi; k++) {
	double w = trap->ffzki[k];
	//double w = creal(f2inva(trap, 0.1525 - trap->zki[k]));
	double V = 0.0;
	double dV = 0.0;
	for(i=0; i<trap->ntot-1; i++) {
	    double Ai = trap->Q[i] * trap->xiw[i] - trap->V[i];
	    double AAi = 0.5/OMEGA * trap->Q[i];
	    double Bi = cosh(OMEGA * trap->dxiw[i]);
	    double Ui = trap->xiw[i] + trap->xiw[i+1];

	    double Di = (cosh(OMEGA*(2*w+Ui)) + Bi);
	    double Ei = cosh(OMEGA * (w - trap->xiw[i]));
	    double Eip1 = cosh(OMEGA * (w - trap->xiw[i+1]));
	    double DDi = Eip1 / Ei;

	    V += Ai * Bi / Di - AAi * log(DDi);
	    dV += Ai*Bi*(-2*sinh(OMEGA*(2*w+Ui)))/(Di*Di)
		- AAi /DDi * (
			sinh(OMEGA*(w-trap->xiw[i+1])) * Ei
			- Eip1 *  sinh(OMEGA*(w-trap->xiw[i]))
			) / (Ei*Ei);
	}
	dV *= OMEGA * trap->dffzki[k];
	trap->Vzki[k] = V;
	trap->V1zki[k] = dV;
    }

    gsl_spline_init(trap->spline, trap->zki, trap->Vzki, trap->npi);
    gsl_spline_init(trap->spline1, trap->zki, trap->V1zki, trap->npi);

    return 0;
}

double
trap_get_pot(struct trap *trap, const double z)
{
    return gsl_spline_eval(trap->spline, fabs(z), trap->acc);
}

double
trap_get_pot1(struct trap *trap, const double z)
{
    return (z>=0) ? gsl_spline_eval(trap->spline1, z, trap->acc1)
	: -gsl_spline_eval(trap->spline1, -z, trap->acc1);
}

double
trap_get_pot2(struct trap *trap, const double z)
{
    return gsl_spline_eval_deriv(trap->spline1, fabs(z), trap->acc1);
}

double
trap_get_pot3(struct trap *trap, const double z)
{
    return (z>=0) ? gsl_spline_eval_deriv2(trap->spline1, z, trap->acc1)
	: -gsl_spline_eval_deriv2(trap->spline1, -z, trap->acc1);
}

// DEPRECATED:
double
trap1_get_pot1(struct trap *trap, const double z)
{
    return (z>=0) ? gsl_spline_eval_deriv(trap->spline, z, trap->acc)
	: -gsl_spline_eval_deriv(trap->spline, -z, trap->acc);
}
double
trap1_get_pot2(struct trap *trap, const double z)
{
    return gsl_spline_eval_deriv2(trap->spline, fabs(z), trap->acc);
}
double
trap1_get_pot3(struct trap *trap, const double z)
{
    double dz = (fabs(z) > 0.03) ? 1e-7 : fabs(z)/10.;
    if (dz < 1e-7)
	return 0.0;
    double ap = trap_get_pot2(trap, fabs(z)+dz);
    double am = trap_get_pot2(trap, fabs(z)-dz);
    return ((z>0) ? (ap - am) : (am - ap)) / (2*dz);
}
// END DEPRECATED

void
trap_free(struct trap *trap)
{
    gsl_spline_free(trap->spline);
    gsl_interp_accel_free(trap->acc);

    free(trap->zki);
    free(trap->ffzki);
    free(trap->Vzki);

    free(trap->xiw);

    free(trap->z);
    free(trap->V);

    free(trap);
}


