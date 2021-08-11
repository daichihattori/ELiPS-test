#include<ELiPS/bls12.h>
#include<gmp.h>
#include <assert.h>

// void fp_sub_correct(fp_t *ANS,fp_t *A,fp_t *B){
//     fp_t A_tmp,B_tmp;
//     if(mpn_cmp(A_tmp.x0,B_tmp.x0,FPLIMB)>0) mpn_sub(ANS->x0,A_tmp.x0,B_tmp->x0,FPLIMB);
//     else{

//     }
// }


//A<P,B<P
void fp_sub_carry(fp_t *ANS,fp_t *A,fp_t *B){
#ifdef DEBUG_ASSERT
  assert(mpn_cmp(A->x0, prime, FPLIMB) > 0);
  assert(mpn_cmp(B->x0, prime, FPLIMB) > 0);
#endif
    if(mpn_sub_n(ANS->x0,A->x0,B->x0,FPLIMB)){
        mpn_add_n(ANS->x0,ANS->x0,prime,FPLIMB);
    }
}
int fp_sub_carry_nonmod(fp_t *ANS,fp_t *A,fp_t *B){

    return mpn_sub_n(ANS->x0,A->x0,B->x0,FPLIMB);
}
void test_fp_sub_under_prime_carry(){
    //C=A-B(A<B,A,B<-461bit)
    int n=10000;
    for(int i=0;i<n;i++){
        fp_t A,B,C_1,C_2;
        while(1){
            fp_set_random(&A,state);
            fp_set_random(&B,state);
            if(mpn_cmp(A.x0,B.x0,FPLIMB)<0) break;
        }
        fp_sub_nonmod_single(&C_1,&A,&B);
        fp_sub_carry(&C_2,&A,&B);
        if(fp_cmp(&C_1,&C_2)==0) printf("C_1==C_2¥n");
        else printf("C_1!=C_2¥n");
    }
}
void test_fp_sub_over_prime_carry(){
    //C=A-B(A<B,A,B<-461bit)
    int n=100;
    for(int i=0;i<n;i++){
        fp_t A,B,C_1,C_2;
        while(1){
            fp_set_random(&A,state);
            fp_set_random(&B,state);
            if(mpn_cmp(A.x0,B.x0,FPLIMB)<0) break;
        }
        // mpn_add_n(A.x0,A.x0,prime,FPLIMB);
        // mpn_add_n(A.x0,A.x0,prime,FPLIMB);
        // mpn_add_n(A.x0,A.x0,prime,FPLIMB);
        // mpn_add_n(A.x0,A.x0,prime,FPLIMB);

        mpn_add_n(B.x0,B.x0,B.x0,FPLIMB);
        mpn_add_n(B.x0,B.x0,B.x0,FPLIMB);
        mpn_add_n(B.x0,B.x0,B.x0,FPLIMB);
        mpn_add_n(B.x0,B.x0,B.x0,FPLIMB);


    



        fp_sub_nonmod_single(&C_1,&A,&B);
        fp_sub_carry(&C_2,&A,&B);
        if(fp_cmp(&C_1,&C_2)==0) printf("");
        else printf("C_1!=C_2¥n");
    }
}
int main(){
    bls12_init();
    //test_fp_sub_under_prime_carry();
    test_fp_sub_over_prime_carry();

}