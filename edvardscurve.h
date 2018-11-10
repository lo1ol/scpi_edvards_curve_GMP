#ifndef EDVARDSCURVE_H
#define EDVARDSCURVE_H

#include <gmp.h>

struct __point_t{
    mpz_t x;
    mpz_t y;
};

typedef struct __point_t point_t[1];
typedef const struct __point_t srcpoint_t[1];


struct __proj_point_t{
    mpz_t X;
    mpz_t Y;
    mpz_t Z;
};

typedef struct __proj_point_t proj_point_t[1];
typedef const struct __proj_point_t srcproj_point_t[1];

struct __edvards_curve_t{
    /*
        Данная структура содержит параметры элиптической кривой и переменные для временных вычислений.
        Переменные для временных вычислений будут храниться именно тут, для того, чтобы уменьшить
        кол-во обращений к менеджеру кучи и для того, чтобы с каждй кривой были связананы свои
        временные вычисления

        Кривая эдвардса в данном случае имеет вид

        a * x^2 + y^2 == c*(1 + d * x^2 * y^2)

        Что обобщает представление простого эдвардса и искривленного
    */
    size_t r;
    mpz_t a;
    mpz_t c;
    mpz_t d;
    mpz_t p;
    mpz_t hint;
    mpz_t _p_minus2, _one;
    mpz_t _Z3, _X1Y2, _X2Y1, _AX1X2, _Y1Y2, _DX1X2Y1Y2 , _Z1Z2, _temp;
    point_t temp_point;
    proj_point_t temp_proj_point;
};

typedef struct __edvards_curve_t edvards_curve_t[1];


//иницизиализирует характеристики кривой. Должна вызываться при начале работе с кривой
void edvards_init(edvards_curve_t curve, mpz_srcptr a, mpz_srcptr c, mpz_srcptr d, mpz_srcptr p);

//функция для сложения двух точек на заданной элиптичесской кривой
void edvards_add(edvards_curve_t curve, point_t result, srcpoint_t P, srcpoint_t Q);

//функция для нахождения обратной точки к заданной точке на элиптической кривой
void edvards_neg(edvards_curve_t curve, point_t result, srcpoint_t P);

//функция вычисления кратной точки на заданной кривой
void edvards_mult(edvards_curve_t curve, point_t result, mpz_srcptr k, srcpoint_t x);

//функция для освобождения занимаемых рессурсов
void edvards_clear(edvards_curve_t curve);


//инициализрует точку (делает ее с координатами 0, 0, что не принадлежит любой элиптической кривой эдвардса!!!)
void point_init(point_t result);

//инициализирует точку на кривой эдвардса копируя ззаданыне значения координат
void point_init_set(point_t result, mpz_srcptr x, mpz_srcptr y);

//задает значения координат точек (копируя)
void point_set(point_t result, mpz_srcptr x, mpz_srcptr y);

//высвобождает память занимаемую кординатами переданной точки
void point_clear(point_t result);

//инициализирует точку в проективной форме
void proj_point_init(proj_point_t P);

//инициализиет точку и устанавливает значение
void proj_point_init_set(proj_point_t P, mpz_srcptr X, mpz_srcptr Y, mpz_srcptr Z);

//высообождает выделенную под проектиную точку
void proj_point_clear(proj_point_t P);

//устанавливает значение точки в проективной форме
void proj_point_set(proj_point_t P, mpz_srcptr X, mpz_srcptr Y, mpz_srcptr Z);

//складывает две точки в проективной форме на заданной кривой Эдвардса
void proj_point_add(edvards_curve_t curve, proj_point_t result, srcproj_point_t P, srcproj_point_t Q);

//вычисляет кратную точку в проективной форме
void proj_point_mult(edvards_curve_t curve, proj_point_t result, mpz_srcptr k, srcproj_point_t P);

//переводит точку из афинного представления в проетивную
void convert_to_proj_point(proj_point_t result, srcpoint_t x);

//перевоит точку из проективной в афинную форму
void convert_to_point(edvards_curve_t curve, point_t result, srcproj_point_t P);


// важно заметить, что входные параметры и результат в данных функциях
// могут быть одной и той же переменной. Функции будут работать с ними абсолютно нормально
// пример: edvards_add(crv, x, x, x);

#endif // EDVARDSCURVE_H
