#include <gmp.h>
#include "edvardscurve.h"
#include "montarith.h"


void edvards_init(edvards_curve_t curve, mpz_srcptr a, mpz_srcptr c, mpz_srcptr d, mpz_srcptr p)
{

    curve->r = mpz_sizeinbase(p, 2);
    mpz_t r2; // вычисляем r^2 для перевода всех констант заранее в представление монтгомери
    mpz_init_set_si(r2, 1); // это необходимо для сохранения эффективности работы прогаммы и опущения
    mpz_mul_2exp(r2, r2, 2*curve->r); //подобных вычислений в результате непосредственной работы программы
    mpz_mod(r2, r2, p);


    mpz_init(curve->hint); // -p^(1) для mod (2^r)
    mont_p_inv_neg(curve->hint, p, curve->r);

    mpz_init(curve->a);
    mont_redc(curve->a, a, r2, p, curve->r, curve->hint);

    mpz_init(curve->c);
    mont_redc(curve->c, c, r2, p, curve->r, curve->hint);

    mpz_init(curve->d);
    mont_redc(curve->d, d, r2, p, curve->r, curve->hint);

    mpz_init_set(curve->p, p);

    mpz_init_set_si(curve->_one, 1); // единица необходима во многих вычислениях на элиптической кривой.
                                    // выжно заметить, что она в представлении монтгомери это r^(-1)

    mpz_init(curve->_p_minus2); // для получения обратного элемента к заданному (используя малую теорему Ферма)
    mpz_sub_ui(curve->_p_minus2, p, 2);

    mpz_inits(curve->_Z3, curve->_X1Y2, curve->_X2Y1, curve->_AX1X2, curve->_Y1Y2,
              curve->_DX1X2Y1Y2, curve->_Z1Z2, curve->_temp, NULL); //переменные для промежуточных вычислений

    point_init(curve->temp_point);
    proj_point_init(curve->temp_proj_point);

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

    mpz_clears(curve->_Z3, curve->_X1Y2, curve->_X2Y1, curve->_AX1X2, curve->_Y1Y2,
              curve->_DX1X2Y1Y2, curve->_Z1Z2, curve->_temp, NULL);

    mpz_clears(curve->_one, curve->_p_minus2, NULL);
    point_clear(curve->temp_point);
    proj_point_clear(curve->temp_proj_point);
}

//сложение двух точек в афинном прдставлении на заданной кривой
void edvards_add(edvards_curve_t curve, point_t result, srcpoint_t P, srcpoint_t Q)
{
    proj_point_t proj_result, proj_P, proj_Q;
    proj_point_init(proj_result);
    proj_point_init(proj_P);
    proj_point_init(proj_Q);

    convert_to_proj_point(proj_P, P);
    convert_to_proj_point(proj_Q, Q);

    proj_point_add(curve, proj_result, proj_P, proj_Q);

    convert_to_point(curve, result, proj_result);

    proj_point_clear(proj_P);
    proj_point_clear(proj_Q);
    proj_point_clear(proj_result);
}

//инициализация точек в афинном представлении
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
    proj_point_t proj_result, proj_P;
    proj_point_init(proj_result);
    proj_point_init(proj_P);
    convert_to_proj_point(proj_P, P);


    proj_point_mult(curve, proj_result, k, proj_P);
    convert_to_point(curve, result, proj_result);

    proj_point_clear(proj_P);
    proj_point_clear(proj_result);
}


void proj_point_init(proj_point_t P)
{
    mpz_inits(P->X, P->Y, P->Z, NULL);
}


void proj_point_init_set(proj_point_t P, mpz_srcptr X, mpz_srcptr Y, mpz_srcptr Z)
{
    proj_point_init(P);
    proj_point_set(P, X, Y,Z);
}


void proj_point_clear(proj_point_t P)
{
    mpz_clears(P->X, P->Y, P->Z, NULL);
}


void proj_point_set(proj_point_t result, mpz_srcptr X, mpz_srcptr Y, mpz_srcptr Z)
{
    mpz_set(result->X, X);
    mpz_set(result->Y, Y);
    mpz_set(result->Z, Z);
}


void convert_to_proj_point(proj_point_t result, srcpoint_t x)
{
    mpz_set(result->X, x->x);
    mpz_set(result->Y, x->y);
    mpz_set_si(result->Z, 1);
}

void convert_to_point(edvards_curve_t curve, point_t result, srcproj_point_t P)
{
    mont_pow(curve->_temp, P->Z, curve->_p_minus2, curve->p, curve->r, curve->hint);
    mont_redc(result->x, P->X, curve->_temp, curve->p, curve->r, curve->hint);
    mont_redc(result->y, P->Y, curve->_temp, curve->p, curve->r, curve->hint);

    mont_redc(result->x, result->x, curve->_one, curve->p, curve->r, curve->hint);
    mont_redc(result->y, result->y, curve->_one, curve->p, curve->r, curve->hint);
}



void proj_point_add(edvards_curve_t curve, proj_point_t result, srcproj_point_t P, srcproj_point_t Q)
{
    size_t r = curve->r; //вычисляем r, для заданного p. r = ceil(log_2(p))

    //задаем удобные имена для промежуточных переменных.
    //Чтобы не писать curve->something
    mpz_ptr p = curve->p;
    mpz_ptr hint = curve->hint;
    mpz_ptr X1Y2 = curve->_X1Y2, X2Y1 = curve->_X2Y1, AX1X2 = curve->_AX1X2,
            Y1Y2 = curve->_Y1Y2, DX1X2Y1Y2 = curve->_DX1X2Y1Y2,
            Z1Z2 = curve->_Z1Z2, temp = curve->_temp;

    mpz_srcptr X1 = P->X;
    mpz_srcptr Y1 = P->Y;
    mpz_srcptr Z1 = P->Z;

    mpz_srcptr X2 = Q->X;
    mpz_srcptr Y2 = Q->Y;
    mpz_srcptr Z2 = Q->Z;

    int replace_after = 0; // флаг показывающий, что нужно высвободить X3 и Y3 помле окончания вычислений
                           // выставляется когда один из операндов совпадает с результатом.
    mpz_ptr X3;
    mpz_ptr Y3;
    mpz_ptr Z3;

    if(result == P || result == Q){ // если результат совпал с одним из результатов, производим вычисления в
                                    // другом месте
        X3 = curve->temp_proj_point->X;
        Y3 = curve->temp_proj_point->Y;
        Z3 = curve->temp_proj_point->Z;
        replace_after = 1;
    }
    else{
        X3 = result->X;
        Y3 = result->Y;
        Z3 = result->Z;
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
    mont_redc(Z1Z2, Z1, Z2, p, r, hint); // Z1Z2 == Z1'*Z2'


    mont_redc(DX1X2Y1Y2, X1Y2, X2Y1, p, r, hint); // DX1X2Y1Y2 == d*X1'*X2'*Y1'*Y2'*r
    mont_redc(DX1X2Y1Y2, curve->d, DX1X2Y1Y2, p, r, hint);
    if(mpz_cmp_si(DX1X2Y1Y2, 0) != 0){     // меняем знак. DX1X2Y1Y2 == -d*X1'*X2'*Y1'*Y2'*r
        mpz_sub(DX1X2Y1Y2, p, DX1X2Y1Y2);
    }

    mont_redc(temp, Z1Z2, Z1Z2, p, r, hint); // temp == Z1'^2 * Z2'^2 r
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
    mont_redc(temp, Z1Z2, Z1Z2, p, r, hint);  // temp == Z1'^2 * Z2'^2 r
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


    if (replace_after){ // Если изначально нельзя было производить вычисления в result
        proj_point_set(result, X3, Y3, Z3); // устанавливаем результат и очищаем память
    }
}

void proj_point_mult(edvards_curve_t curve, proj_point_t result, mpz_srcptr k, srcproj_point_t P)
{
    proj_point_t garbage, Q; // garbage -- переменная для вычисления ненужного хлама
                        // для сглаживания занятости процессора. Q -- точки вида 2^i * P

    struct __proj_point_t *dest; //dest становится или result или garbage в зависимости от того, нужно ли класть
                            // очередну степень в результат

    mpz_t zero, c; // инициализируем нуль и переводим, c в исходное состояние (не в форме монтгомери)

    mpz_ptr one = curve->_one;
    mpz_init(zero);
    mpz_init(c);
    mont_redc(c, curve->c, one, curve->p, curve->r, curve->hint);

    proj_point_init_set(garbage, P->X, P->Y, P->Z); // garbage = P
    proj_point_init_set(Q, P->X, P->Y, P->Z); //Q == P
    proj_point_set(result, zero, c, one); // нулевая точка

    for(mp_bitcnt_t t=0; t != mpz_sizeinbase(k, 2); ++t){ //реализация алгоритма двочного возведения в степень (умножения)
                                                          //используя лестницу Монтгомери
        dest = (mpz_tstbit(k, t) == 1)?  result : garbage;
        if(t==0){
            proj_point_add(curve, dest, dest, Q); // в первый раз просто прибавляем к dest точку P
        }
        else{
            proj_point_add(curve, Q, Q, Q); // теперь возводим в степень и прибавляем к dest
            proj_point_add(curve, dest, dest, Q);
        }
    }
    proj_point_clear(garbage);
    proj_point_clear(Q);
    mpz_clears(c, zero, NULL);
}


