#ifndef EDVARDSCURVE_H
#define EDVARDSCURVE_H

#include <gmp.h>

struct __edvards_curve_t{
    mpz_t a;
    mpz_t c;
    mpz_t d;
    mpz_t p;
    mpz_t hint;
    mpz_t r2;
};

typedef struct __edvards_curve_t edvards_curve_t[1];

struct __point_t{
    mpz_t x;
    mpz_t y;
};

typedef struct __point_t point_t[1];
typedef const struct __point_t srcpoint_t[1];

void edvards_init(edvards_curve_t curve, mpz_srcptr a, mpz_srcptr c, mpz_srcptr d, mpz_srcptr p);
void edvards_add(edvards_curve_t curve, point_t result, srcpoint_t P, srcpoint_t Q);
void edvards_neg(edvards_curve_t curve, point_t result, srcpoint_t P);
void edvards_mult(edvards_curve_t curve, point_t result, mpz_srcptr k, srcpoint_t x);
void edvards_clear(edvards_curve_t curve);

void point_init(point_t result);
void point_init_set(point_t result, mpz_srcptr x, mpz_srcptr y);
void point_set(point_t result, mpz_srcptr x, mpz_srcptr y);
void point_clear(point_t result);

#endif // EDVARDSCURVE_H
