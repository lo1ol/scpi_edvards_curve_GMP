#include<gmp.h>
#include "montarith.h"


static mpz_t s;

void mont_init(){
    mpz_init(s); // s = t*(-p^(-1)) в формуле умножения по монтгомери
                 // инициализиурся тут для повышения быстродейтсивя
                 // чтобы каждый раз память не выделять
}

void mont_clear(){
    mpz_clear(s); // высвобождаем память для s
}

void mont_p_inv_neg(mpz_t result, mpz_srcptr p, unsigned long r_num){
    mpz_t r;
    mpz_init_set_si(r, 1);
    mpz_mul_2exp(r, r, r_num); // r -> 2^r (т.е. число тут в явном виде)

    //используем обобщение малой теоремы Ферма (теорему Эйлера)
    //(p^(-1)) mod 2^r == p^(phi(2^r)-1) mod 2^r == p^(2^(r-1) - 1) mod 2^r


    mpz_t Q;
    mpz_init_set(Q, p);
    mpz_set_si(result, 1);

    for(unsigned long i=1; i+1<=r_num; i++){ //on each i-th iteration compute res = p^1*p^2*p^4*...*p^(2^i).
                                             //After that res = p^1*p^2*p^4*...*p^(2^(k-1)) (ну я же англичянин :) )
        if(i==1){
            mpz_set(result, p); // на первой итерации просто присваиваем result == p^1
            continue;
        }
        mpz_mul(Q, Q, Q); // возводим в степень Q = p^(2^i)
        mpz_mod_2exp(Q, Q, r_num); // Q == p^(2^i) mod 2^r
        mpz_mul(result, result, Q); // result == p^1*p^2*p^4*...*p^(2^i)
        mpz_mod_2exp(result, result, r_num); // result == p^1*p^2*p^4*...*p^(2^i) mod 2^r
    }
    if(mpz_cmp_si(result, 0)!= 0){ //result = -p^(-1)
        mpz_sub(result, r, result);
    }
    mpz_clears(r, Q, NULL);
}

void mont_add(mpz_t result, mpz_srcptr x, mpz_srcptr y, mpz_srcptr p){
    mpz_add(result, x, y);
    if(mpz_cmp(result, p) >= 0){ // в случае переполнения приводим по модулю p
        mpz_sub(result, result, p);
    }
}

void mont_redc(mpz_t result, mpz_srcptr x, mpz_srcptr y, mpz_srcptr p, unsigned long r, mpz_srcptr hint){
    mpz_ptr p_inv_neg;
    int clear_afeter = 0;
    if(hint == NULL){ // если -p^(-1) не переданно, вычиляем сами
        mpz_t p_inv_neg_val;
        mpz_init(p_inv_neg_val);
        mont_p_inv_neg(p_inv_neg_val, p, r);
        p_inv_neg = p_inv_neg_val;
        clear_afeter =1; // выставляем флаг, что в конце нужно высвободить память
    }
    else{
        p_inv_neg =  (mpz_ptr) hint; // иначе используем подсказку
    }
    mpz_mul(result, x, y); // result == t == x * y

    mpz_mul(s, result, p_inv_neg); // s == t * -p^(-1)
    mpz_mod_2exp(s, s, r); //s == t * -p^(-1) mod 2^r
    mpz_mul(s, s, p); //s == (t * -p^(-1) mod 2^r) * p

    mpz_add(result, result, s); // result == t + s
    mpz_div_2exp(result, result, r); // result == (t+s) / r
    if(mpz_cmp(result, p) >= 0){ // если произошло переполнение, приводим по модулю p
        mpz_sub(result, result, p);
    }
    if(clear_afeter){
        mpz_clear(p_inv_neg);
    }
}

void mont_pow(mpz_t result, mpz_srcptr x_gotten, mpz_srcptr pow_gotten, mpz_srcptr p, unsigned long r, mpz_srcptr hint){
    mpz_t x, pow;
    mpz_init_set(x, x_gotten); // защита от совпадения результата и параметров
    mpz_init_set(pow, pow_gotten);

    mpz_t p_inv_neg;
    if(hint == NULL){
        mpz_init(p_inv_neg);
        mont_p_inv_neg(p_inv_neg, p, r);
    }
    else{
        mpz_init_set(p_inv_neg, hint);
    }

    mpz_t garbage, Q;
    mpz_ptr dest;

    mpz_init_set_si(garbage, 1); // garbage == r mod p т.е. 1 единица в представлении монтгомери
    mpz_mul_2exp(garbage, garbage, r);
    mpz_mod(garbage, garbage, p);

    mpz_init_set(Q, garbage);
    mpz_set(result, garbage);


    for(mp_bitcnt_t t=0; t != mpz_sizeinbase(pow, 2); ++t){ // используем лесинку монтгомери адаптированную
                                                            // для возведения степени в форме монгомери

        dest = (mpz_tstbit(pow, t) == 1)? result : garbage; // будем изменять или мусор или результат
                                                            // в зависимости от бита степени
        if (t==0){
            mpz_set(Q, x);
            mpz_set(dest, x);
        }
        else{
            mont_redc(Q, Q, Q, p, r, p_inv_neg);
            mont_redc(dest, dest, Q, p ,r, p_inv_neg);
        }
    }

    mpz_clears(x, pow, garbage, Q, p_inv_neg, NULL);
}

void mont_addredc(mpz_t result, mpz_srcptr x, mpz_srcptr y, mpz_srcptr p, unsigned long r, mpz_srcptr hint)
{
    mpz_t temp; // используем временную переменную для хранения результата умножения
    mpz_init(temp);
    mont_redc(temp, x, y, p, r, hint);
    mont_add(result, temp, result, p);
    mpz_clear(temp);
}
