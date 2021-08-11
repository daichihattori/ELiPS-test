#include <ELiPS/bls12.h>
#include <gmp.h>
void fp12_mul_lazy_montgomery_complete(fp12_t *ANS,fp12_t *A,fp12_t *B){
    static fp6_t tmp1_fp6,tmp2_fp6,tmp3_fp6;
    static fpd6_t tmp1_fpd6,tmp2_fpd6,tmp3_fpd6,buf_fpd6;
    static fp6_t out;
    //set
    fp6_mul_nonmod_montgomery(&tmp1_fpd6,&A->x1,&B->x1);//b*d
    fp6_add_nonmod_single(&tmp2_fp6,&A->x0,&A->x1);//a+b
    fp6_add_nonmod_single(&tmp3_fp6,&B->x0,&B->x1);//c+d
    fp6_mul_nonmod_montgomery(&tmp2_fpd6,&tmp2_fp6,&tmp3_fp6);//(a+b)(c+d)
        fp6_mod_montgomery_double(&out,&tmp2_fpd6);
        //fp6_println("out_A_0=",&out);printf("\n");
    fp6_mul_nonmod_montgomery(&tmp3_fpd6,&A->x0,&B->x0);//a*c
        fp6_mod_montgomery_double(&out,&tmp3_fpd6);
        //fp6_println("out_A_1=",&out);printf("\n");

    //x0
    fp6_mul_basis_nonmod_double(&buf_fpd6,&tmp1_fpd6);//b*d*v
        fp6_mod_montgomery_double(&out,&buf_fpd6);
        //fp6_println("out_A_2=",&out);printf("\n");
    fp6_add_nonmod_double(&buf_fpd6,&buf_fpd6,&tmp3_fpd6);//a*c+b*d*v
        fp6_mod_montgomery_double(&out,&buf_fpd6);
        //fp6_println("out_A_3=",&out);printf("\n");
    fp6_mod_montgomery_double(&ANS->x0,&buf_fpd6);
        //fp6_println("chk1=",&ANS->x0);printf("\n");

    //x1
    fp6_sub_nonmod_double(&buf_fpd6,&tmp2_fpd6,&tmp1_fpd6);
        fp6_mod_montgomery_double(&out,&buf_fpd6);
        //fp6_println("out_A_4=",&out);printf("\n");
    fp6_sub_nonmod_double(&buf_fpd6,&buf_fpd6,&tmp3_fpd6);
        fp6_mod_montgomery_double(&out,&buf_fpd6);
        //fp6_println("out_A_5=",&out);printf("\n");
    fp6_mod_montgomery_double(&ANS->x1,&buf_fpd6);
        //fp6_println("chk2=",&ANS->x1);printf("\n");

        //under to debug
    fp6_mul_lazy_montgomery(&tmp2_fp6,&A->x1,&B->x1);//b*d
    fp6_add(&tmp1_fp6,&A->x0,&A->x1);//a+b
    fp6_add(&ANS->x1,&B->x0,&B->x1);//c+d
    fp6_mul_lazy_montgomery(&ANS->x1,&tmp1_fp6,&ANS->x1);//(a+b)(c+d)
        //fp6_println("out_B_0=",&ANS->x1);printf("\n");
    fp6_mul_lazy_montgomery(&tmp1_fp6,&A->x0,&B->x0);//a*c
        //fp6_println("out_B_1=",&tmp1_fp6);printf("\n");

    //x0
    fp6_mul_basis(&ANS->x0,&tmp2_fp6);//b*d*v
        //fp6_println("out_B_2=",&ANS->x0);printf("\n");
    fp6_add(&ANS->x0,&ANS->x0,&tmp1_fp6);//a*c+b*d*v
        //fp6_println("out_B_3=",&ANS->x0);printf("\n");
        //fp6_println("chk1=",&ANS->x0);printf("\n");

    //x1
    fp6_sub(&ANS->x1,&ANS->x1,&tmp1_fp6);
        //fp6_println("out_B_4=",&ANS->x1);printf("\n");
    fp6_sub(&ANS->x1,&ANS->x1,&tmp2_fp6);
        //fp6_println("out_B_5=",&ANS->x1);printf("\n");
        //fp6_println("chk2=",&ANS->x1);printf("\n");
        //getchar();

}

void test_fp6_mul_test(int cnt){
    fp12_t A,B,C_1,C_2;
    for(int i=0;i<cnt;i++){
        fp12_set_random(&A,state);
        fp12_to_montgomery(&A,&A);
        fp12_set_random(&B,state);
        fp12_to_montgomery(&B,&B);

        fp12_mul_lazy_montgomery(&C_1,&A,&B);
        fp12_mul_lazy_montgomery_complete(&C_2,&A,&B);
        if(fp12_cmp(&C_1,&C_2)==0) printf("ok!¥n");
        else printf("ng!¥n");

    }
}
int main(){
    bls12_init();
    test_fp6_mul_test(100);
    return 0;
}