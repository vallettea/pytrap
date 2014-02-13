#ifndef _CONFORMAL_MAP_H_
#define _CONFORMAL_MAP_H_
#include <complex.h>
#include <gsl/gsl_spline.h>

#define DEFAULT_TRAP_NTOT	14
#define DEFAULT_TRAP_N		11
#define DEFAULT_TRAP_R		3.79146919431e-2	// 8mm / 211mm
#define DEFAULT_TRAP_RZ		6.16113744076e-2	// 13mm / 211mm

struct trap {
    size_t n;			// number of points for the radius R zone
    size_t ntot;
    double R;			// radius (adimensional)
    double Rz;			// radius (adimensional)
    double *z;			// coordinates of the electrodes
    double *V;			// potential at the electrodes

    //
    double z2;
    double b, sqrtb;
    double aa, abm, abp, ac, ad;
    //
    double *xiw;
    double *dxiw;
    double *Q;

    size_t npi;			// number of points for interpolation
    double *zki;		// coordinates of interpolated points
    double *ffzki;
    double *dffzki;
    double *Vzki;		// V(z)
    double *V1zki;		// V'(z)

    gsl_interp_accel *acc;
    gsl_interp_accel *acc1;
    gsl_spline *spline;
    gsl_spline *spline1;
};

extern struct trap *trap_alloc(const double *, const size_t);
extern void trap_free(struct trap *);

extern int trap_init(struct trap *);

extern int trap_set_pot(struct trap *, const double *);

extern double trap_get_pot(struct trap *, const double);
extern double trap_get_pot1(struct trap *, const double);
extern double trap_get_pot2(struct trap *, const double);
extern double trap_get_pot3(struct trap *, const double);

/* these three are deprecated. Use for checking only: */
extern double trap1_get_pot1(struct trap *, const double);
extern double trap1_get_pot2(struct trap *, const double);
extern double trap1_get_pot3(struct trap *, const double);

/* the conformal mapping and its derivative: */
extern double complex f2inva(const struct trap *, const double);
extern double complex df2inva(const struct trap *, const double);


#define FSOLVER			gsl_multiroot_fsolver_hybrids
#define FSOLVER_EPSILON		0.0
#define INTERP_ZMIN		0.0
#define INTERP_ZMAX		0.25
#define INTERP_DZ		1e-5

#endif
