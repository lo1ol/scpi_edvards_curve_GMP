#include <gmp.h>
#include "edvardscurve.h"
#include "montarith.h"


void edvards_init(edvards_curve_t curve, mpz_srcptr a, mpz_srcptr c, mpz_srcptr d, mpz_srcptr p)
{

    mpz_t r2; // вычисляем r^2 для перевода всех констант заранее в представление монтгомери
    mpz_init_set_si(r2, 1); // это необходимо для сохранения эффективности работы прогаммы и опущения
    mpz_mul_2exp(r2, r2, 2*mpz_sizeinbase(p, 2)); //подобных вычислений в результате непосредственной работы программы
    mpz_mod(r2, r2, p);


    mpz_init(curve->hint); // -p^(1) для mod (2^r)
    mont_p_inv_neg(curve->hint, p, mpz_sizeinbase(p, 2));

    mpz_init(curve->a);
    mont_redc(curve->a, a, r2, p, mpz_sizeinbase(p, 2), curve->hint);

    mpz_init(curve->c);
    mont_redc(curve->c, c, r2, p, mpz_sizeinbase(p, 2), curve->hint);

    mpz_init(curve->d);
    mont_redc(curve->d, d, r2, p, mpz_sizeinbase(p, 2), curve->hint);

    mpz_init_set(curve->p, p);

    mpz_init_set_si(curve->one, 1); // единица необходима во многих вычислениях на элиптической кривой.
                                    // выжно заметить, что она в представлении монтгомери это r^(-1)

    mpz_init(curve->p_minus2); // для получения обратного элемента к заданному (используя малую теорему Ферма)
    mpz_sub_ui(curve->p_minus2, p, 2);

    mpz_inits(curve->Z3, curve->X1Y2, curve->X2Y1, curve->AX1X2, curve->Y1Y2,
              curve->DX1X2Y1Y2, curve->Z1Z2, curve->temp, NULL); //переменные для промежуточных вычислений

    mpz_clear(r2);
}


void edvards_clear(edvards_curve_t curve)
{
    //высвобождает занимаемую память под ппараметры и промежуточные вычисления
    mpz_clear(curve->a);
    mpz_clear(curve->c);
    mpz_clear(curve->d);
    mpz_clear(curve->p);
    mpz_clear(curve->hint);

    mpz_clears(curve->Z3, curve->X1Y2, curve->X2Y1, curve->AX1X2, curve->Y1Y2,
              curve->DX1X2Y1Y2, curve->Z1Z2, curve->temp, NULL);

    mpz_clears(curve->one, curve->p_minus2, NULL);
}


void edvards_add(edvards_curve_t curve, point_t result, srcpoint_t P, srcpoint_t Q)
{
    size_t r = mpz_sizeinbase(curve->p, 2); //вычисляем r, для заданного p. r = ceil(log_2(p))

    //задаем удобные имена для промежуточных переменных.
    //Чтобы не писать curve->something
    mpz_ptr p_minus_2 = curve->p_minus2;
    mpz_ptr p = curve->p;
    mpz_ptr hint = curve->hint;
    mpz_ptr one = curve->one;
    mpz_ptr X1Y2 = curve->X1Y2, X2Y1 = curve->X2Y1, AX1X2 = curve->AX1X2,
            Y1Y2 = curve->Y1Y2, DX1X2Y1Y2 = curve->DX1X2Y1Y2,
            Z1Z2 = curve->Z1Z2, temp = curve->temp;

    mpz_srcptr X1 = P->x;
    mpz_srcptr Y1 = P->y;

    mpz_srcptr X2 = Q->x;
    mpz_srcptr Y2 = Q->y;

    int replace_after = 0; // флаг показывающий, что нужно высвободить X3 и Y3 помле окончания вычислений
                           // выставляется когда один из операндов совпадает с результатом.
    mpz_ptr X3;
    mpz_ptr Y3;
    mpz_ptr Z3 = curve->Z3;

    if(result == P || result == Q){ // если результат совпал с одним из результатов, производим вычисления в
                                    // другом месте
        mpz_t xval, yval;
        mpz_inits(xval, yval, NULL);
        X3 = xval;
        Y3 = yval;
        replace_after = 1;
    }
    else{
        X3 = result->x;
        Y3 = result->y;
    }

    // формулы используемые в проективной форме имеют вид
    // X3 = Z1*Z2*(X1*Y2 + X2Y1)(Z1^2*Z2^2 + d*X1*X2*Y1*Y2)
    // Y3 = Z1*Z2*(Y1*Y2 - X1X2)(Z1^2*Z2^2 + d*X1*X2*Y1*Y2)
    // Z3 = (Z1^4*Z2^4 - (d*X1*X2*Y1*Y2)^2)

    // C этого момента считаем, что r равно 2^r

    // перед началом работы программы, мы неяно считаем, что перешли в передставление монтгомери
    //используя следующий переход
    // (x : y) -> (x : y : 1) = (phi(x') : phi(y') : phi(r^(-1))) = (x'*r : y'*2r : r^(-1)*r)

    mont_redc(X1Y2, X1, Y2, p, r, hint); // X1Y2 == X1'*Y2'*r
    mont_redc(X2Y1, X2, Y1, p, r, hint); // X2Y2 == X2'*Y1'*r
    mont_redc(Z1Z2, one, one, p, r, hint); // Z1Z2 == (r^(-1))^2 * r


    mont_redc(DX1X2Y1Y2, X1Y2, X2Y1, p, r, hint); // DX1X2Y1Y2 == d*X1'*X2'*Y1'*Y2'*r
    mont_redc(DX1X2Y1Y2, curve->d, DX1X2Y1Y2, p, r, hint);
    if(mpz_cmp_si(DX1X2Y1Y2, 0) != 0){     // меняем знак. DX1X2Y1Y2 == -d*X1'*X2'*Y1'*Y2'*r
        mpz_sub(DX1X2Y1Y2, p, DX1X2Y1Y2);
    }

    mont_redc(temp, Z1Z2, Z1Z2, p, r, hint); // temp == Z1'^2 * Z2'^2 r == (r^(-1))^4 * r
    mont_add(temp, temp, DX1X2Y1Y2, p); // temp == (Z1'^2 * Z2'^2 - d*X1'*X2'*Y1'*Y2') * r
    mont_add(X3, X1Y2, X2Y1, p); // X3 == (X1'*Y2' + X2'*Y1')*r
    mont_redc(X3, temp, X3, p, r, hint); // X3 == (X1'*Y2' + X2'*Y1')*(Z1'^2 * Z2'^2 - d*X1'*X2'*Y1'*Y2')*r
    mont_redc(X3, X3, Z1Z2, p, r, hint); // X3 == Z1'*Z2'*(X1'*Y2' + X2'*Y1')*(Z1'^2 * Z2'^2 - d*X1'*X2'*Y1'*Y2')*r



    mont_redc(Y1Y2, Y1, Y2, p, r, hint); // Y1Y2 == Y1'*Y2'*r
    mont_redc(AX1X2, X1, X2, p, r, hint); // X1X2 == -a*X1'*X2'*r
    mont_redc(AX1X2, curve->a, AX1X2, p, r, hint);
    if(mpz_cmp_si(AX1X2, 0) != 0){
        mpz_sub(AX1X2, p, AX1X2);
    }
    mont_add(Y3, Y1Y2, AX1X2, p); // Y3 == (Y1'*Y2' - a*X1'*X2')*r


    if(mpz_cmp_si(DX1X2Y1Y2, 0) != 0){ // DX1X2Y1Y2 == d*X1'*X2'*Y1'*Y2'*r
        mpz_sub(DX1X2Y1Y2, p, DX1X2Y1Y2);
    }
    mont_redc(temp, Z1Z2, Z1Z2, p, r, hint);  // temp == Z1'^2 * Z2'^2 r == (r^(-1))^4 * r
    mont_add(temp, temp, DX1X2Y1Y2, p); // temp == (Z1'^2 * Z2'^2 + d*X1'*X2'*Y1'*Y2') * r
    mont_redc(Y3, temp, Y3, p, r, hint); // Y3 == (Y1'*Y2' - a*X1'*X2')*(Z1'^2 * Z2'^2 + d*X1'*X2'*Y1'*Y2') * r
    mont_redc(Y3, Y3, Z1Z2, p, r, hint); // Y3 == Z1'*Z2'*(Y1'*Y2' - a*X1'*X2')*(Z1'^2 * Z2'^2 + d*X1'*X2'*Y1'*Y2') * r



    mpz_set(Z3, curve->c); // Z3 == c*r
    mont_redc(DX1X2Y1Y2, DX1X2Y1Y2, DX1X2Y1Y2, p, r, hint); // DX1X2Y1Y2 = -(d*X1'*X2'*Y1'*Y2')^2 * r
    if(mpz_cmp_si(DX1X2Y1Y2, 0) != 0){
        mpz_sub(DX1X2Y1Y2, p, DX1X2Y1Y2);
    }

    mont_redc(temp, Z1Z2, Z1Z2, p, r, hint); //temp == Z1'^2 * Z2'^2 * r
    mont_redc(temp, temp, temp, p, r, hint); //temp == Z1'^4 * Z2'^4 * r
    mont_add(temp, temp, DX1X2Y1Y2, p); //temp == (Z1'^2 * Z2'^2 - (d*X1'*X2'*Y1'*Y2')^2) * r
    mont_redc(Z3, Z3, temp, p, r, hint); // Z3 == c * (Z1'^2 * Z2'^2 - (d*X1'*X2'*Y1'*Y2')^2) * r



    mont_pow(Z3, Z3, p_minus_2, p, r, hint); // Z3 == (c * (Z1'^2 * Z2'^2 - (d*X1'*X2'*Y1'*Y2')^2))^(-1) * r
    mont_redc(X3, X3, Z3, p, r, hint); // переходим обрато в аффиную форму x3 == X3'/Z3' *r
    mont_redc(Y3, Y3, Z3, p, r, hint); // y3 == Y3'/Z3' * r



    mont_redc(X3, X3, one, p, r, hint); // Вопользуемся свойством перехода от пердставления в монтгомери к обычному
                                        // Используя операцию умноежения в монтгомери x3 = X3'/Z3' = X3/Z3
    mont_redc(Y3, Y3, one, p, r, hint); // y3 = Y3'/Z3' = Y3/Z3


    if (replace_after){ // Если изначально нельзя было производить вычисления в result
        point_set(result, X3, Y3); // устанавливаем результат и очищаем память
        mpz_clears(X3, Y3, NULL);
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


void edvards_mult(edvards_curve_t curve, point_t result, mpz_srcptr k, srcpoint_t P)
{
    point_t garbage, Q; // garbage -- переменная для вычисления ненужного хлама
                        // для сглаживания занятости процессора. Q -- точки вида 2^i * P

    struct __point_t *dest; //dest становится или result или garbage в зависимости от того, нужно ли класть
                            // очередну степень в результат

    mpz_t zero, c; // инициализируем нуль и переводим, c в исходное состояние (не в форме монтгомери)

    mpz_ptr one = curve->one;
    mpz_init(zero);
    mpz_init(c);
    mont_redc(c, curve->c, one, curve->p, mpz_sizeinbase(curve->p, 2), curve->hint);

    point_init_set(garbage, P->x, P->y); // garbage = P
    point_init_set(Q, P->x, P->y); //Q == P
    point_set(result, zero, c); // нулевая точка

    for(mp_bitcnt_t t=0; t != mpz_sizeinbase(k, 2); ++t){ //реализация алгоритма двочного возведения в степень (умножения)
                                                          //используя лестницу Монтгомери
        dest = (mpz_tstbit(k, t) == 1)?  result : garbage;
        if(t==0){
            edvards_add(curve, dest, dest, Q); // в первый раз просто прибавляем к dest точку P
        }
        else{
            edvards_add(curve, Q, Q, Q); // теперь возводим в степень и прибавляем к dest
            edvards_add(curve, dest, dest, Q);
        }
    }
    point_clear(garbage);
    point_clear(Q);
    mpz_clears(c, zero, NULL);

}
