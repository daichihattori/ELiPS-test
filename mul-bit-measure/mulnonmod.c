#include<ELiPS/bls12.h>
#include<gmp.h>
#include<time.h>



int main(){
    bls12_init();
    mpz_t tmp_A,tmp_B;
    mp_bitcnt_t bit=477;
    fp_t A,B,C_mod;
    fpd_t C;
    int cnt=10000;
    float mul_time=0;
    float mod_time=0;
    struct timeval tv_A, tv_B;

    mpz_init(tmp_A);
    mpz_init(tmp_B);
    fp_init(&A);
    fp_init(&B);
    //fpd_init(&C);

    for(bit=461;bit<511;bit++){
        mul_time=0;
        mod_time=0;

        for(int i=0;i<cnt;i++){
            mpz_rrandomb(tmp_A,state,bit);
            mpz_rrandomb(tmp_B,state,bit);

            mpn_set_mpz_size(A.x0,tmp_A,FPLIMB);
            mpn_set_mpz_size(B.x0,tmp_B,FPLIMB);
            fp_to_montgomery(&A,&A);
            fp_to_montgomery(&B,&B);

            gettimeofday(&tv_A,NULL);
            fp_mul_nonmod(&C,&A,&B);
            gettimeofday(&tv_B,NULL);
            mul_time+=timedifference_msec(tv_A, tv_B);

            gettimeofday(&tv_A,NULL);
            mpn_mod_montgomery(C_mod.x0,FPLIMB,C.x0,FPLIMB2);
            gettimeofday(&tv_B,NULL);
            mod_time+=timedifference_msec(tv_A, tv_B);

            if(mpn_cmp(C_mod.x0,prime,FPLIMB)>0) printf("error?\n");
        }
        printf("fp mul %ld bit    %.7f[ms] %.7f[ms]\n",bit,mul_time/cnt,mod_time/cnt);

    }
    
    return 0;

    



}