#include <gmp.h>
#include <edvardscurve.h>
#include "montarith.h"

void edvards_init(edvards_curve_t curve, mpz_srcptr a, mpz_srcptr c, mpz_srcptr d, mpz_srcptr p)
{

    mpz_init_set_si(curve->r2, 1);
    mpz_mul_2exp(curve->r2, curve->r2, 2*mpz_sizeinbase(p, 2));
    mpz_mod(curve->r2, curve->r2, p);


    mpz_init(curve->hint);
    mont_p_inv_neg(curve->hint, p, mpz_sizeinbase(p, 2));

    mpz_init(curve->a);
    mont_redc(curve->a, a, curve->r2, p, mpz_sizeinbase(p, 2), curve->hint);

    mpz_init(curve->c);
    mont_redc(curve->c, c, curve->r2, p, mpz_sizeinbase(p, 2), curve->hint);

    mpz_init(curve->d);
    mont_redc(curve->d, d, curve->r2, p, mpz_sizeinbase(p, 2), curve->hint);


    mpz_init_set(curve->p, p);
}


void edvards_clear(edvards_curve_t curve)
{
    mpz_clear(curve->a);
    mpz_clear(curve->c);
    mpz_clear(curve->d);
    mpz_clear(curve->p);
    mpz_clear(curve->hint);
}


void edvards_add(edvards_curve_t curve, point_t result, srcpoint_t P, srcpoint_t Q)
{
    size_t r = mpz_sizeinbase(curve->p, 2);
    mpz_t p;
    mpz_init_set(p, curve->p);

    mpz_srcptr X1 = P->x;
    mpz_srcptr Y1 = P->y;

    mpz_srcptr X2 = Q->x;
    mpz_srcptr Y2 = Q->y;

    mpz_t X3;
    mpz_t Y3;
    mpz_t Z3;
    mpz_inits(X3, Y3, Z3, NULL);

    mpz_t one;
    mpz_init_set_si(one, 1);

    {
        mpz_t X1Y2, X2Y1, AX1X2, Y1Y2, DX1X2Y1Y2 , Z1Z2pow2, temp;
        mpz_inits(X1Y2, X2Y1, DX1X2Y1Y2, Z1Z2pow2, AX1X2, Y1Y2, temp,  NULL);
        mont_redc(X1Y2, X1, Y2, p, r, curve->hint);
        mont_redc(X2Y1, X2, Y1, p, r, curve->hint);
        mont_redc(Z1Z2pow2, one, one, p, r, curve->hint);
        mont_redc(Z1Z2pow2, Z1Z2pow2, Z1Z2pow2, p, r, curve->hint);


        mont_redc(DX1X2Y1Y2, X1Y2, X2Y1, p, r, curve->hint);
        mont_redc(DX1X2Y1Y2, curve->d, DX1X2Y1Y2, p, r, curve->hint);
        if(mpz_cmp_si(DX1X2Y1Y2, 0) != 0){
            mpz_sub(DX1X2Y1Y2, p, DX1X2Y1Y2);
        }

        mont_add(temp, Z1Z2pow2, DX1X2Y1Y2, p);
        mont_add(X3, X1Y2, X2Y1, p);
        mont_redc(X3, temp, X3, p, r, curve->hint);
        mont_redc(X3, X3, Z1Z2pow2, p, r, curve->hint);
        mont_redc(X3, X3, Z1Z2pow2, p, r, curve->hint);



        mont_redc(AX1X2, X1, X2, p, r, curve->hint);
        mont_redc(AX1X2, curve->a, AX1X2, p, r, curve->hint);
        mont_redc(Y1Y2, Y1, Y2, p, r, curve->hint);
        if(mpz_cmp_si(AX1X2, 0) != 0){
            mpz_sub(AX1X2, p, AX1X2);
        }
        mont_add(Y3, Y1Y2, AX1X2, p);


        if(mpz_cmp_si(DX1X2Y1Y2, 0) != 0){
            mpz_sub(DX1X2Y1Y2, p, DX1X2Y1Y2);
        }
        mont_add(temp, Z1Z2pow2, DX1X2Y1Y2, p);
        mont_redc(Y3, temp, Y3, p, r, curve->hint);
        mont_redc(Y3, Y3, Z1Z2pow2, p, r, curve->hint);
        mont_redc(Y3, Y3, Z1Z2pow2, p, r, curve->hint);



        mpz_set(Z3, curve->c);
        mont_redc(Z3, Z1Z2pow2, Z3, p, r, curve->hint);
        mont_redc(Z3, one, Z3, p, r, curve->hint);
        mont_redc(Z3, one, Z3, p, r, curve->hint);
        mont_redc(DX1X2Y1Y2, DX1X2Y1Y2, DX1X2Y1Y2, p, r, curve->hint);
        if(mpz_cmp_si(DX1X2Y1Y2, 0) != 0){
            mpz_sub(DX1X2Y1Y2, p, DX1X2Y1Y2);
        }
        mont_redc(temp, Z1Z2pow2, Z1Z2pow2, p, r, curve->hint);
        mont_add(temp, temp, DX1X2Y1Y2, p);
        mont_redc(Z3, Z3, temp, p, r, curve->hint);



        mpz_sub_ui(p, p, 2);
        mont_pow(Z3, Z3, p, curve->p, r, curve->hint);
        mont_redc(X3, X3, Z3, curve->p, r, curve->hint);
        mont_redc(Y3, Y3, Z3, curve->p, r, curve->hint);



        mont_redc(X3, X3, one, curve->p, r, curve->hint);
        mont_redc(Y3, Y3, one, curve->p, r, curve->hint);

        point_set(result, X3, Y3);
        mpz_clears(one, AX1X2, Y1Y2, DX1X2Y1Y2, p, Z3, X1Y2, X2Y1, Z1Z2pow2, temp, NULL);
    }
}


void point_init(point_t result)
{
    mpz_init(result->x);
    mpz_init(result->y);
}

void point_init_set(point_t result, mpz_srcptr x, mpz_srcptr y)
{
    mpz_init_set(result->x, x);
    mpz_init_set(result->y, y);
}


void point_clear(point_t result)
{
    mpz_clear(result->x);
    mpz_clear(result->y);
}


void point_set(point_t result, mpz_srcptr x, mpz_srcptr y)
{
    mpz_set(result->x, x);
    mpz_set(result->y, y);
}


void edvards_neg(edvards_curve_t curve, point_t result, srcpoint_t P)
{
    if(mpz_cmp_si(P->x, 0)!= 0){
        mpz_sub(result->x, curve->p, P->x);
    }
    mpz_set(result->y, P->y);
}


void edvards_mult(edvards_curve_t curve, point_t result, mpz_srcptr k, srcpoint_t x)
{
    point_t garbage, Q;
    struct __point_t *dest;
    mpz_t zero, c, one;
    mpz_init(zero);
    mpz_init_set_si(one, 1);
    mpz_init(c);
    mont_redc(c, curve->c, one, curve->p, mpz_sizeinbase(curve->p, 2), curve->hint);

    point_init_set(garbage, x->x, x->y);
    point_init_set(Q, x->x, x->y);
    point_set(result, zero, c);

    for(mp_bitcnt_t t=0; t != mpz_sizeinbase(k, 2); ++t){
        dest = (mpz_tstbit(k, t) == 1)?  result : garbage;
        if(t==0){
            edvards_add(curve, dest, dest, Q);
        }
        else{
            edvards_add(curve, Q, Q, Q);
            edvards_add(curve, dest, dest, Q);
        }
    }
    point_clear(garbage);
    point_clear(Q);
    mpz_clears(c, one, NULL);

}
