#include<ELiPS/bls12.h>
#include<gmp.h>
#include<time.h>

mp_bitcnt_t fp_bit(fp_t *A){
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_mpn_size(tmp,A->x0,FPLIMB);
}