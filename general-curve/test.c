#include <ELiPS/bls12.h>
mp_limb_t curve_b_tmp[FPLIMB];

int main(){
    bls12_init();
    mpz_t Xmod72;
    mpz_init(Xmod72);

    mpz_mod_ui(Xmod72,X_z,72);
    gmp_printf("a=%Zd",Xmod72);
    return 0;
}