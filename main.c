#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include "montarith.h"
#include "edvardscurve.h"

int main()
{

    mpz_t x, y, p, a, d, one, zero, k;
    mpz_init_set_str(k, "3F FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF"
                        "FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF"
                        "C9 8C DB A4 65 06 AB 00 4C 33 A9 FF 51 47 50 2C"
                        "C8 ED A9 E7 A7 69 A1 26 94 62 3C EF 47 F0 23 ED", 16);
    mpz_init_set_si(one, 1);
    mpz_init_set_si(zero, 0);
    mpz_init_set_str(x, "12", 16);

    mpz_init_set_str(y, "46 9A F7 9D 1F B1 F5 E1 6B 99 59 2B 77 A0 1E 2A"
                        "0F DF B0 D0 17 94 36 8D 9A 56 11 7F 7B 38 66 95"
                        "22 DD 4B 65 0C F7 89 EE BF 06 8C 5D 13 97 32 F0"
                        "90 56 22 C0 4B 2B AA E7 60 03 03 EE 73 00 1A 3D", 16);


    mpz_init_set_str(a, "01", 16);

    mpz_init_set_str(d, "00 9E 4F 5D 8C 01 7D 8D 9F 13 A5 CF 3C DF 5B FE"
                        "4D AB 40 2D 54 19 8E 31 EB DE 28 A0 62 10 50 43"
                        "9C A6 B3 9E 0A 51 5C 06 B3 04 E2 CE 43 E7 9E 36"
                        "9E 91 A0 CF C2 BC 2A 22 B4 CA 30 2D BB 33 EE 75"
                        "50", 16);

    mpz_init_set_str(p, "00 FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF"
                        "FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF"
                        "FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF"
                        "FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FD"
                        "C7", 16);

    edvards_curve_t crv;
    edvards_init(crv, a, one, d, p);
    point_t p1, p2, p3;
    point_init_set(p1, x, y);
    point_init_set(p2, x, y);
    point_init(p3);
    int i = 0;
    while(++i!=100){
        edvards_mult(crv, p3, k, p1);
    }
    gmp_printf("%Zd*\n(%Zd, %Zd)=\n (%Zd, %Zd)\n mod %Zd\n",
               k, p1->x, p1->y, p3->x, p3->y, p);

    edvards_clear(crv);
    mpz_clears(x, y, p, NULL);
    point_clear(p1);
    point_clear(p2);
    point_clear(p3);
    return 0;
}
