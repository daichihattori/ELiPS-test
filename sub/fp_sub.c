#include <ELiPS/bls12.h>
#include <gmp.h>
void fp_nonmod_single_test(fp_t *ans,fp_t *a,fp_t *b){
    // if(mpn_sub_n(ans->x0,a->x0,b->x0,FPLIMB)){
    //     mpn_add_n(ans->x0,ans->x0,prime,FPLIMB);
    // }

    //１：単純にBをmodとる
    mp_limb_t dumy[FPLIMB];
    mpn_tdiv_qr(dumy,b->x0,0,b->x0,FPLIMB,prime,FPLIMB);
    if(mpn_sub_n(ans->x0,a->x0,b->x0,FPLIMB)){
        mpn_add_n(ans->x0,ans->x0,prime,FPLIMB);
    }
}

void fp_mod(fp_t *ans,mp_limb_t *a,mp_size_t size_a){
    #ifdef DEBUG_COST_A
    cost_mod_nomal++;
    #endif
    mp_limb_t dumy[size_a];
    mpn_tdiv_qr(dumy,ans->x0,0,a,size_a,prime,FPLIMB);
}
int main(){
    bls12_init(); 
    int cnt=1000;
    fp_t A,B,ANS1,ANS2;
    //check1: A<P,B<P
    fp_set_random(&A,state);
    fp_set_random(&B,state);
    fp_sub_nonmod_single(&ANS1,&A,&B);
    fp_nonmod_single_test(&ANS2,&A,&B);
    int flag=1;
    for(int i=0;i<cnt;i++){
        if(fp_cmp(&ANS1,&ANS2)!=0){
            flag=0;
            break;
        }
    }
    if(flag) printf("ok\n");
    else printf("ng\n");

    //check2: A<P,B>P
    //mpn_add_n(B.x0,B.x0,prime);



    return 0;
}