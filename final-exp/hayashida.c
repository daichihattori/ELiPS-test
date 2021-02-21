#include<ELiPS/bls12.h>
// ELiPS
/* Cost: 4Px + Px/2 + S + 8M + 4I + 1F1 + 1F2 + 1F3*/
void g1g2_to_g3_hardpart_compress_montrick(fp12_t *ANS, fp12_t *A){
    fp12_t tmp,t0,t1,t2,t3,t4,t5, test,At;
    fp12_init(&tmp);
    fp12_init(&t0);
    fp12_init(&t1);
    fp12_init(&t2);
    fp12_init(&t3);
    fp12_init(&t5);
    fp12_init(&t4);
    fp12_init(&test);
    fp12_init(&At);
    //HARDPART
    fp12_sqr_GS_lazy_montgomery(&t0, A);
    bls12_fp12_pow_X_compress_montrick(&t1, &t0);

    bls12_fp12_pow_X2_compress_montrick(&t2,&t1);//t2:=t1^(u2);
    fp12_frobenius_map_p6_montgomery(&t3,A);//t3:=f^(-1);
    fp12_mul_lazy_montgomery(&t1,&t3,&t1);//t1:=t3*t1;
    fp12_frobenius_map_p6_montgomery(&t1,&t1);//t1:=t1^(-1);
    fp12_mul_lazy_montgomery(&t1,&t1,&t2);//t1:=t1*t2;

    bls12_fp12_pow_X_compress_montrick(&t2,&t1);//t2:=t1^(u);
    bls12_fp12_pow_X_compress_montrick(&t3,&t2);//t3:=t2^(u);
    fp12_frobenius_map_p6_montgomery(&t1,&t1);//t1:=t1^(-1);

    fp12_mul_lazy_montgomery(&t3,&t1,&t3);//t3:=t1*t3;
    fp12_frobenius_map_p6_montgomery(&t1,&t1);//t1:=t1^(-1);
    fp12_frobenius_map_p3_lazy_montgomery(&t1,&t1);//t1:=t1^(p^3);
    fp12_frobenius_map_p2_montgomery(&t2,&t2);//t2:=t2^(p^2);


    fp12_mul_lazy_montgomery(&t1,&t1,&t2);//t1:=t1*t2;
    bls12_fp12_pow_X_compress_montrick(&t2,&t3);//t2:=t3^(u);
    fp12_mul_lazy_montgomery(&t2,&t2,&t0);//t2:=t2*t0;
    fp12_mul_lazy_montgomery(&t2,&t2,A);//t2:=t2*f;
    fp12_mul_lazy_montgomery(&t1,&t1,&t2);//t1:=t1*t2;

    fp12_frobenius_map_p1_lazy_montgomery(&t2,&t3);//t2:=t3^p;
    fp12_mul_lazy_montgomery(ANS,&t1,&t2);//t1:=t1*t2;
    //fp12_mod_montgomery(ANS,ANS);
}

// Hayashida
/* Cost: 4Px + Px/2 + S + 7M + 2Ic + 1F1 + 1F2 */
void g1g2_to_g3_hardpart_compress_montrick(fp12_t *ANS, fp12_t *A){
    fp12_t t1,t2,t3,t4;
    fp12_init(&t1);
    fp12_init(&t2);
    fp12_init(&t3);
    fp12_init(&t4);
    //HARDPART
    fp12_sqr_GS_lazy_montgomery(&t1, A);
    bls12_fp12_pow_X_compress_montrick(&t2, &t1);
    bls12_fp12_pow_X2_compress_montrick(&t3,&t2);//t3:=t2^(u2);
    fp12_frobenius_map_p6_montgomery(&t2,&t2);//t4:=f^(-1);
    fp12_mul_lazy_montgomery(&t2,&t2,&t3);//t2:=t4*t2;
    fp12_mul_lazy_montgomery(&t2,&t2,A);//t2:=t4*t2;
    bls12_fp12_pow_X_compress_montrick(&t3, &t2);
    fp12_frobenius_map_p1_lazy_montgomery(&t2,&t2);
    fp12_mul_lazy_montgomery(&t2,&t2,&t3);//t2:=t4*t2;
    bls12_fp12_pow_X_compress_montrick(&t3, &t2);
    bls12_fp12_pow_X_compress_montrick(&t3, &t3);
    fp12_frobenius_map_p2_montgomery(&t4,&t2);
    fp12_frobenius_map_p6_montgomery(&t2,&t2);//t4:=f^(-1);
    fp12_mul_lazy_montgomery(&t2,&t2,&t3);//t2:=t4*t2;
    fp12_mul_lazy_montgomery(&t2,&t2,&t4);//t2:=t4*t2;
    fp12_mul_lazy_montgomery(&t1,&t1,A);//t2:=t4*t2;
    fp12_mul_lazy_montgomery(ANS,&t2,&t1);//t2:=t4*t2;
}
int main(){
    bls12_init();
    return 0;
}