#include<gmp.h>
#include "montarith.h"

void mont_p_inv_neg(mpz_t result, mpz_srcptr p, unsigned long r_num){
    mpz_t r;
    mpz_init_set_si(r, 1);
    mpz_mul_2exp(r, r, r_num);


    mpz_t Q;
    mpz_init_set(Q, p);
    mpz_set_si(result, 1);

    for(unsigned long i=1; i+1<=r_num; i++){ //on each i-th iteration compute res = p^1*p^2*...*p^i.
                                   //After that res = p^1*p^2*...*p^{k-1}
        if(i==1){
            mpz_set(result, p);
            continue;
        }
        mpz_mul(Q, Q, Q);
        mpz_mod_2exp(Q, Q, r_num);
        mpz_mul(result, result, Q);
        mpz_mod_2exp(result, result, r_num);
    }
    if(mpz_cmp_si(result, 0)!= 0){
        mpz_sub(result, r, result);
    }
    mpz_clears(r, Q, NULL);
}

void mont_add(mpz_t result, mpz_srcptr x, mpz_srcptr y, mpz_srcptr p){
    mpz_add(result, x, y);
    if(mpz_cmp(result, p) >= 0){
        mpz_sub(result, result, p);
    }
}

void mont_redc(mpz_t result, mpz_srcptr x, mpz_srcptr y, mpz_srcptr p, unsigned long r, mpz_srcptr hint){
    mpz_t p_inv_neg;
    if(hint == NULL){
        mpz_init(p_inv_neg);
        mont_p_inv_neg(p_inv_neg, p, r);
    }
    else{
        mpz_init_set(p_inv_neg, hint);
    }
    mpz_mul(result, x, y);

    mpz_t s;
    mpz_init(s);
    mpz_mul(s, result, p_inv_neg);
    mpz_mod_2exp(s, s, r);
    mpz_mul(s, s, p);

    mpz_add(result, result, s);
    mpz_div_2exp(result, result, r);
    if(mpz_cmp(result, p) >= 0){
        mpz_sub(result, result, p);
    }
    mpz_clears(p_inv_neg, s, NULL);
}

void mont_pow(mpz_t result, mpz_srcptr x_gotten, mpz_srcptr pow_gotten, mpz_srcptr p, unsigned long r, mpz_srcptr hint){
    mpz_t x, pow;
    mpz_init_set(x, x_gotten);
    mpz_init_set(pow, pow_gotten);

    mpz_t p_inv_neg;
    if(hint == NULL){
        mpz_init(p_inv_neg);
        mont_p_inv_neg(p_inv_neg, p, r);
    }
    else{
        mpz_init_set(p_inv_neg, hint);
    }

    mpz_t garbage, Q, dest;

    mpz_init_set_si(garbage, 1);
    mpz_mul_2exp(garbage, garbage, r);
    mpz_mod(garbage, garbage, p);

    mpz_init_set_si(Q, 1);
    mpz_set(result, garbage);
    mpz_init(dest);


    for(mp_bitcnt_t t=0; t != mpz_sizeinbase(pow, 2); ++t){
        (mpz_tstbit(pow, t) == 1)? mpz_swap(dest, result):mpz_swap(dest, garbage);
        if (t==0){
            mpz_set(Q, x);
            mpz_set(dest, x);
        }
        else{
            mont_redc(Q, Q, Q, p, r, p_inv_neg);
            mont_redc(dest, dest, Q, p ,r, p_inv_neg);
        }
        (mpz_tstbit(pow, t) == 1)? mpz_swap(dest, result):mpz_swap(dest, garbage);
    }

    mpz_clears(x, pow, garbage, Q, dest, p_inv_neg, NULL);
}

void mont_addredc(mpz_t result, mpz_srcptr x, mpz_srcptr y, mpz_srcptr p, unsigned long r, mpz_srcptr hint)
{
    mpz_t temp;
    mpz_init(temp);
    mont_redc(temp, x, y, p, r, hint);
    mont_add(result, temp, result, p);
    mpz_clear(temp);
}
