#include<stdio.h>
#include<gmp.h>
#define FPLIMB 1
typedef struct{
    mp_limb_t x0[FPLIMB];
}fp_t;
mp_limb_t prime[FPLIMB];
int main(){

    //set prime
    prime[0]=3;

    //test fp
    fp_t A,B,C;
    A.x0[0]=0;
    B.x0[0]=3;
    int carry;
    carry=mpn_sub_n(C.x0,A.x0,B.x0,FPLIMB);
    gmp_printf("C=%Nd carry=%d",C.x0,FPLIMB,carry);

    
    return 0;
}