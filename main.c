#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "montarith.h"
#include "edvardscurve.h"

static const char* gX = "12";

static const char* gY = "46 9A F7 9D 1F B1 F5 E1 6B 99 59 2B 77 A0 1E 2A"
                        "0F DF B0 D0 17 94 36 8D 9A 56 11 7F 7B 38 66 95"
                        "22 DD 4B 65 0C F7 89 EE BF 06 8C 5D 13 97 32 F0"
                        "90 56 22 C0 4B 2B AA E7 60 03 03 EE 73 00 1A 3D";

static const char* gA = "01";

static const char* gC = "01";

static const char* gD = "00 9E 4F 5D 8C 01 7D 8D 9F 13 A5 CF 3C DF 5B FE"
                        "4D AB 40 2D 54 19 8E 31 EB DE 28 A0 62 10 50 43"
                        "9C A6 B3 9E 0A 51 5C 06 B3 04 E2 CE 43 E7 9E 36"
                        "9E 91 A0 CF C2 BC 2A 22 B4 CA 30 2D BB 33 EE 75"
                        "50";

static const char* gP = "00 FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF"
                        "FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF"
                        "FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF"
                        "FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FD"
                        "C7";

static const char* gK = "3F FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF"
                        "FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF"
                        "C9 8C DB A4 65 06 AB 00 4C 33 A9 FF 51 47 50 2C"
                        "C8 ED A9 E7 A7 69 A1 26 94 62 3C EF 47 F0 23 ED";

/**
 * @brief initilize_curve Функия инициализирующая кривую в форме Эдвардса параметрами с 512 битными значениями взятыми из
 * <a href="http://ftp.lissi.ru/TK26/%CF%E0%F0%E0%EC%E5%F2%F0%FB%20%FD%EB%EB%E8%EF%F2%E8%F7%E5%F1%EA%E8%F5%20%EA%F0%E8%E2%FB%F5.pdf">Тут</a>
 *
 *
 * @param crv кривая, которую инициализируем
 */
void initilize_curve(edvards_curve_t crv);

/**
 * @brief work_with_user фукнция запрашивающая у пользователя число, кратная точка от которого вычисляется
 */
void work_with_user(void);

/**
 * @brief test_program функция, тестирующую корректоность работы вычислений на кривой Эдвардса (запускает все тесты)
 */
void test_program(void);

/**
 * @brief test_order фуункция, тестирующая вычисления на кривой Эдвардса на то, что заданная точка имеет указанный порядок
 * @return  логическое значение успешности пройденности теста
 */
int test_order(void);

/**
 * @brief test_associativity функия, тесттирующая вычисления на кривой Эдвардса на сохаенения ассоциаативности
 * @param k кол-во повторений
 * @param seed зерно для генератора случайных чисел
 * @return логическое значение успешности пройденности теста
 */
int test_associativity(int k, unsigned seed);

/**
 * @brief test_communicativeness функия, тесттирующая вычисления на кривой Эдвардса на сохаенения коммунтиативности
 * @param k кол-во повторений
 * @param seed зерно для генератора случайных чисел
 * @return логическое значение успешности пройденности теста
 */
int test_communicativeness(int k, unsigned seed);


int main(int argc, char* argv[])
{
    mont_init();
    if (argc == 2 && strcmp(argv[1], "-t") == 0) {
        test_program();
    } else {
        work_with_user();
    }

    mont_clear();
    return 0;
}

void initilize_curve(edvards_curve_t crv)
{
    mpz_t a, d, c, p;
    mpz_init_set_str(a, gA, 16);

    mpz_init_set_str(c, gC, 16);

    mpz_init_set_str(d, gD, 16);

    mpz_init_set_str(p, gP, 16);

    edvards_init(crv, a, c, d, p);
    mpz_clears(a, c, d, p, NULL);
}

void work_with_user()
{
    mpz_t x, y, k;
    edvards_curve_t crv;
    point_t p1, p2;

    initilize_curve(crv);
    mpz_init(k);
    mpz_init_set_str(x, gX, 16);
    mpz_init_set_str(y, gY, 16);
    point_init_set(p1, x, y);
    point_init(p2);

    gmp_printf("Enter k number: ");
    gmp_scanf("%Zd", k);

    if (mpz_cmp_si(k, 0) < 0) {
        printf("k must be non negative!\n");
        edvards_clear(crv);
        mpz_clears(x, y, k, NULL);
        point_clear(p1);
        point_clear(p2);
        exit(EXIT_FAILURE);
    }

    edvards_mult(crv, p2, k, p1);

    gmp_printf("\n%Zd*\n(%Zd, %Zd)=\n (%Zd, %Zd)\n mod %Zd\n",
               k, p1->x, p1->y, p2->x, p2->y, crv->p);

    edvards_clear(crv);
    mpz_clears(x, y, k, NULL);
    point_clear(p1);
    point_clear(p2);
}

void test_program()
{
    const char* k_success = "Success!";
    const char* k_bad = "Bad!";
    int n = 100;
    unsigned seed = 0;

    printf("Enter seed: ");
    scanf("%u", &seed);

    printf("Order Test: %s\n", test_order() == 1 ? k_success : k_bad);
    printf("Associativity Test for %d random elements: %s\n", n, test_associativity(n, seed) == 1 ? k_success : k_bad);
    printf("Communicativeness Test for %d random elements: %s\n", n, test_communicativeness(n, seed) == 1 ? k_success : k_bad);
}

int test_order()
{
    int result;
    mpz_t x, y, k;
    edvards_curve_t crv;
    point_t Q, P;

    initilize_curve(crv);
    mpz_init_set_str(x, gX, 16);
    mpz_init_set_str(y, gY, 16);
    mpz_init_set_str(k, gK, 16);
    point_init_set(P, x, y);
    point_init(Q);

    edvards_mult(crv, Q, k, P);

    if (mpz_cmp_si(Q->x, 0) == 0 && mpz_cmp_si(Q->y, 1) == 0) {
        result = 1;
    } else {
        result = 0;
    }

    mpz_clears(x, y, k, NULL);
    edvards_clear(crv);
    point_clear(P);
    point_clear(Q);

    return result;
}

int test_associativity(int k, unsigned seed)
{
    int result = 1;
    mp_bitcnt_t n;
    mpz_t x, y, r12, r1, r2;
    gmp_randstate_t rand_st;
    edvards_curve_t crv;
    point_t R, R1, R2, R12, R12_;

    initilize_curve(crv);
    n = crv->r;

    gmp_randinit_default(rand_st);
    gmp_randseed_ui(rand_st, seed);
    mpz_init_set_str(x, gX, 16);
    mpz_init_set_str(y, gY, 16);
    mpz_inits(r1, r2, r12, NULL);
    point_init_set(R, x, y);
    point_init(R1); point_init(R2); point_init(R12); point_init(R12_);

    for (int i=0; i != k; ++i) {
        mpz_rrandomb(r1, rand_st, n);
        mpz_rrandomb(r2, rand_st, n);
        mpz_add(r12, r1, r2);

        edvards_mult(crv, R1, r1, R);
        edvards_mult(crv, R2, r2, R);
        edvards_add(crv, R12_, R1, R2);
        edvards_mult(crv, R12, r12, R);

        if (mpz_cmp(R12_->x, R12->x) != 0 && mpz_cmp(R12_->y, R12->y) != 0) {
            result = 0;
        }
    }

    mpz_clears(x, y, r1, r2, r12, NULL);
    edvards_clear(crv);
    point_clear(R); point_clear(R1); point_clear(R2); point_clear(R12); point_clear(R12_);
    gmp_randclear(rand_st);

    return result;
}

int test_communicativeness(int k, unsigned seed)
{
    int result = 1;
    mp_bitcnt_t n;
    mpz_t x, y, r1, r2;
    gmp_randstate_t rand_st;
    edvards_curve_t crv;
    point_t R, R1, R2, R12, R12_;

    initilize_curve(crv);
    n = crv->r;

    gmp_randinit_default(rand_st);
    gmp_randseed_ui(rand_st, seed);
    mpz_init_set_str(x, gX, 16);
    mpz_init_set_str(y, gY, 16);
    mpz_inits(r1, r2, NULL);
    point_init_set(R, x, y);
    point_init(R1); point_init(R2); point_init(R12); point_init(R12_);

    for (int i=0; i != k; ++i) {
        mpz_rrandomb(r1, rand_st, n);
        mpz_rrandomb(r2, rand_st, n);
        edvards_mult(crv, R1, r1, R);
        edvards_mult(crv, R2, r2, R);
        edvards_add(crv, R12, R1, R2);
        edvards_add(crv, R12_, R2, R1);

        if (mpz_cmp(R12_->x, R12->x) != 0 && mpz_cmp(R12_->y, R12->y) != 0) {
            result = 0;
        }
    }

    mpz_clears(x, y, r1, r2, NULL);
    edvards_clear(crv);
    point_clear(R); point_clear(R1); point_clear(R2); point_clear(R12); point_clear(R12_);
    gmp_randclear(rand_st);

    return result;
}
