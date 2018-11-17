#ifndef EDVARDSCURVE_H
#define EDVARDSCURVE_H

#include <gmp.h>

/**
 * @brief The __point_t struct Точка на элептической кривой в афинном представлении
 */
struct __point_t{
    mpz_t x;
    mpz_t y;
};

/**
 * @brief point_t Тип для реализации семантики ссылки для афинной точки
 */
typedef struct __point_t point_t[1];

/**
 * @brief srcpoint_t Тип для реализации семантики константной ссылки афинной точки
 */
typedef const struct __point_t srcpoint_t[1];

/**
 * @brief The __proj_point_t struct Тип точки на элептической кривой в финном представлении
 */
struct __proj_point_t{
    mpz_t X;
    mpz_t Y;
    mpz_t Z;
};

/**
 * @brief proj_point_t Тип для реализации семантики ссылки для проективной формы точки
 */
typedef struct __proj_point_t proj_point_t[1];

/**
 * @brief srcproj_point_t Тип для реализации семантики констаной ссылки для проективной формы точки
 */
typedef const struct __proj_point_t srcproj_point_t[1];

/**
 * @brief The __edvards_curve_t struct Тип представляющий кривую в форме Эдвардса
 *
 * Данная структура содержит параметры элиптической кривой в форме Монтгомери и переменные для временных вычислений.
 * Переменные для временных вычислений будут храниться именно тут, для того, чтобы уменьшить
 * кол-во обращений к менеджеру кучи и для того, чтобы с каждй кривой были связананы свои
 * временные вычисления
 *
 * Кривая эдвардса в данном случае имеет вид
 * \f$ a  x^2 + y^2 = c(1 + d  x^2 * y^2)\f$
 * Что обобщает представление простого эдвардса и искривленного Эдвардса
 *
 */
struct __edvards_curve_t{
    size_t r; ///< параметр \f$ r\f$, опредяющийся как наименьшее число такое, что \f$ 2^r \ge p\f$
    mpz_t a; ///< параметр \f$ a \f$ Элептической кривой
    mpz_t c; ///< параметр \f$ c \f$ Элептической кривой
    mpz_t d; ///< параметр \f$ d \f$ Элептической кривой
    mpz_t p; ///< простое число, по модулю которого производятся вычисления
    mpz_t hint; ///< \f$ -p^{-1} \f$ в кольце \f$ 2^r \f$ (используется в качестве подсказки при вычислениях в форме Монтгомери)
    mpz_t _p_minus2, _one;
    mpz_t _Z3, _X1Y2, _X2Y1, _AX1X2, _Y1Y2, _DX1X2Y1Y2 , _Z1Z2, _temp;
    point_t temp_point;
    proj_point_t temp_proj_point;
};

/**
 * @brief edvards_curve_t Тип для реализации семантики ссылки для кривой эдвардса
 */
typedef struct __edvards_curve_t edvards_curve_t[1];


/**
 * @brief edvards_init Иницизиализирует параметры кривой. Должна вызываться при начале работе с кривой
 * @param curve кривая которую нужно инициализировать
 * @param a параметр \f$ a \f$ Элептической кривой
 * @param c параметр \f$ c \f$ Элептической кривой
 * @param d параметр \f$ d \f$ Элептической кривой
 * @param p простое число, по модулю которого производятся вычисления
 */
void edvards_init(edvards_curve_t curve, mpz_srcptr a, mpz_srcptr c, mpz_srcptr d, mpz_srcptr p);

/**
 * @brief edvards_add функция для сложения двух точек на заданной элиптичесской кривой
 *
 * Вычисления происходят по следующей схеме.
 *  1. данная точка переходит в проективную форму, путем принятия координаты \f$ Z=1 \f$ в проективной форме
 *  2. потом данная точка переходит в представление Монтгомери (по сути ничего в этот момент не происходит, т.к. умножение проективной точки на число не меняет точку)
 *  3. после этого все вычисления производятся в представлении монтгомери с использованием следующих формул:
 *     \f$ X3 = Z1*Z2*(X1*Y2 + X2Y1)(Z1^2*Z2^2 + d*X1*X2*Y1*Y2) \f$
 *     \f$ Y3 = Z1*Z2*(Y1*Y2 - X1X2)(Z1^2*Z2^2 + d*X1*X2*Y1*Y2) \f$
 *     \f$ Z3 = (Z1^4*Z2^4 - (d*X1*X2*Y1*Y2)^2) \f$
 *     Данные формулы были полученны в результате подставления известных формул для афинных координат, замены \f$ x = \frac{X}{Y} \f$  и \f$ y = \frac{Y}{Z} \f$
 *
 *     Формулы для вычисления в афиином представлении следующие
 *     \f$ x_3 = \frac{x_1 y_2+y_1 x_2}{1+d x_1 x_2 y_1 y_2} \f$, \f$ y3 = \frac{y_1 y_2-a x_1 x_2}{1-d x_1 x_2 y_1 y_2} \f$
 *  4. после всех вычислений, мы переводим нашу точку обратно в афинное представление путем деления каждой координаты на Z
 *     Операция деления происходит на самом деле как операция умноежния на обратный элемент вычисляемый как \f$ Z^{-1} = Z^{p-2}\f$
 *
 *
 * @param curve кривая в группе которой происходят все вычисления
 * @param result точка, куда положится результат сложения двух точек \f$ P - Q \f$
 * @param P \f$ P \f$
 * @param Q \f$ Q \f$
 */
void edvards_add(edvards_curve_t curve, point_t result, srcpoint_t P, srcpoint_t Q);

/**
 * @brief edvards_neg функция для нахождения обратной точки к заданной точке на элиптической кривой
 * @param curve кривая в группе которой происходят все вычисления
 * @param result точка, куда положится результат  \f$ -P \f$
 * @param P \f$ P \f$
 */
void edvards_neg(edvards_curve_t curve, point_t result, srcpoint_t P);

/**
 * @brief edvards_mult Функция вычисления кратной точки на заданной кривой
 *
 * При вычислении кратной точки используется лесенка монтгомери для сглаживания производительности процессора при вычсилении
 *
 * @param curve кривая в которой происходят вычисления
 * @param result точка куда будет положен результат \f$ k*P \f$
 * @param k \f$ k \f$ числа, кратную точку от которого мы хотим получить
 * @param P \f$ P \f$ точка, из которой мы хотим получить кратную
 */
void edvards_mult(edvards_curve_t curve, point_t result, mpz_srcptr k, srcpoint_t P);

/**
 * @brief edvards_clear функция для освобождения занимаемых рессурсов
 * @param curve кривая, которую мы хотим очистить
 */
void edvards_clear(edvards_curve_t curve);


/**
 * @brief point_init инициализрует точку (делает ее с координатами 0, 0, что не принадлежит любой элиптической кривой эдвардса!!!)
 * @param P инициализируемая точка
 */
void point_init(point_t P);

/**
 * @brief point_init_set инициализирует точку на кривой эдвардса копируя ззаданыне значения координат
 * @param P результат, куда будут положены две точки
 * @param x координата \f$ x \f$ точки
 * @param y коориданата \f$ y \f$ точки
 */
void point_init_set(point_t P, mpz_srcptr x, mpz_srcptr y);

/**
 * @brief point_set задает значения координат точек (копируя)
 * @param P точка, координаты которой задаются
 * @param x координата \f$ x \f$ точки
 * @param y коориданата \f$ y \f$ точки
 */
void point_set(point_t P, mpz_srcptr x, mpz_srcptr y);

/**
 * @brief point_clear высвобождает память занимаемую кординатами переданной точки
 * @param P точка, которая будет очищена
 */
void point_clear(point_t P);

/**
 * @brief proj_point_init инициализирует точку в проективной форме
 * @param P точка, которая будет инициализирована
 */
void proj_point_init(proj_point_t P);

/**
 * @brief proj_point_init_set инициализиет точку и устанавливает значение
 * @param P точка, кооридинаты которая инициализиурется
 * @param X координата \f$ X \f$ точки
 * @param Y координата \f$ Y \f$ точки
 * @param Z координата \f$ Z \f$ точки
 */
void proj_point_init_set(proj_point_t P, mpz_srcptr X, mpz_srcptr Y, mpz_srcptr Z);

/**
 * @brief proj_point_clear высообождает выделенную под проектиную точку
 * @param P точка, которая будет очищена
 */
void proj_point_clear(proj_point_t P);

/**
 * @brief proj_point_set устанавливает значение точки в проективной форме
 * @param P точка, координаты которой устаналиваюстя
 * @param X координата \f$ X \f$ точки
 * @param Yкоордината \f$ Y \f$ точки
 * @param Z координата \f$ Z \f$ точки
 */
void proj_point_set(proj_point_t P, mpz_srcptr X, mpz_srcptr Y, mpz_srcptr Z);

/**
 * @brief proj_point_add складывает две точки в проективной форме на заданной кривой Эдвардса
 *
 * Все выкладки для вычисления кратной точки описаны в описании функции edvards_add()
 *
 * @param curve кривая, в которой происходят все вычисления
 * @param result точка, куда будет положен результат сложения точек \f$ P \f$ и \f$ Q \f$
 * @param P \f$ P \f$
 * @param Q \f$ Q \f$
 */
void proj_point_add(edvards_curve_t curve, proj_point_t result, srcproj_point_t P, srcproj_point_t Q);

/**
 * @brief proj_point_mult вычисляет кратную точку в проективной форме
 * @param curve кривая, в которой производятся все вычисления
 * @param result точка, куда будет положен результат вычисления кратной точки
 * @param k \f$ k \f$ числа, кратную точку от которого мы хотим получить
 * @param P \f$ P \f$ точка, из которой мы хотим получить кратную
 */
void proj_point_mult(edvards_curve_t curve, proj_point_t result, mpz_srcptr k, srcproj_point_t P);

/**
 * @brief convert_to_proj_point переводит точку из афинного представления в проетивную
 * @param result результат, куда будет положенно проективное представление
 * @param P точка, в афинной форме
 */
void convert_to_proj_point(proj_point_t result, srcpoint_t P);

/**
 * @brief convert_to_point перевоит точку из проективной в афинную форму
 * @param curve кривая, в которой происзводятся вычисления
 * @param result результат конвертирования
 * @param P точка, которая конвертируется в проективной форме
 */
void convert_to_point(edvards_curve_t curve, point_t result, srcproj_point_t P);


#endif // EDVARDSCURVE_H
