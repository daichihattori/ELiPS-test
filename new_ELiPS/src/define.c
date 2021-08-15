#include <ELiPS/define.h>

int bls12_X_binary[bls12_X_length + 1];
int bls12_X2_binary[bls12_X2_length + 1];

//int cost_add,cost_add_ui,cost_sub,cost_sub_ui,cost_mul,cost_mul_ui,cost_sqr,cost_inv,cost_mod;
int cost_add, cost_add_ui, cost_sub, cost_sub_ui, cost_mul, cost_set_neg, cost_sqr, cost_inv, cost_mod;
int cost_add_nonmod, cost_add_nonmod_double, cost_sub_nonmod, cost_sub_nonmod_double, cost_r1shift, cost_mod_nomal;

mp_limb_t buf[FPLIMB], tmp_mul[FPLIMB2], tmp1[FPLIMB], tmp2[FPLIMB];

/*============================================================================*/
/* Pairing functions                                                          */
/*============================================================================*/

gmp_randstate_t state;
mpz_t X_z, prime_z, order_z, trace_z, X_abs_z;
mp_limb_t X_abs[FXLIMB], X2[FXLIMB2], prime[FPLIMB], order[FRLIMB], trace[FPLIMB];
mp_limb_t prime_carry[FPLIMB];
mp_limb_t prime2[FPLIMB2];
mp_limb_t epsilon1[FPLIMB], epsilon2[FPLIMB];
fp_t epsilon1_montgomery, epsilon2_montgomery;
mpz_t efp_total, efp12_total;
mp_limb_t curve_b[FPLIMB];

mpz_t sqrt_power_z;
mpz_t g1_power;

//montgomery
mp_limb_t N[FPLIMB2], RmodP[FPLIMB], R3[FPLIMB];
mp_limb_t Ni_neg;  //Ni_neg=-N^(-1)
mpz_t X_mod_order_z;
/*============================================================================*/
/* Test functions                                                             */
/*============================================================================*/
struct timeval tv_start, tv_end;

float MILLER_OPT_AFFINE, FINALEXP_OPT_AFFINE;
float MILLER_OPT_PROJECTIVE, FINALEXP_OPT_PROJECTIVE;

cost MILLER_OPT_AFFINE_COST, FINALEXP_OPT_AFFINE_COST;
cost MILLER_OPT_PROJECTIVE_COST, FINALEXP_OPT_PROJECTIVE_COST;

mpz_t to_g1_expo, to_g2_expo;
fp_t curve_b_montgomery;
fp_t inv2_montgomery;
