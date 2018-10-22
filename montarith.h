#ifndef MONTARITH_H
#define MONTARITH_H

#include<gmp.h>


void mont_add(mpz_t result, mpz_srcptr x, mpz_srcptr y, mpz_srcptr p);
void mont_redc(mpz_t result, mpz_srcptr x, mpz_srcptr y, mpz_srcptr p, unsigned long r, mpz_srcptr p_inv_hint);
void mont_addredc(mpz_t result, mpz_srcptr x, mpz_srcptr y, mpz_srcptr p, unsigned long r, mpz_srcptr p_inv_hint);
void mont_pow(mpz_t result, mpz_srcptr x, mpz_srcptr pow, mpz_srcptr p, unsigned long r, mpz_srcptr p_inv_hint);
void mont_p_inv_neg(mpz_t result, mpz_srcptr p, unsigned long r);

#endif // MONTARITH_H
