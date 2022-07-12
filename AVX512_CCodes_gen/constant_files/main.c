/*********************************************

	add_mult_poly :

multiplication multiprécision vectorisée :

https://software.intel.com/sites/landingpage/IntrinsicsGuide/#!=undefined


*********************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <immintrin.h>
#include <time.h>

#include "ccount.h"
#include <gmp.h>

#include "add_mult_poly.h"


#define WORD 64

#define NSAMPLES 50

#define NTEST 1000

unsigned long long int START, STOP, START1,STOP1;





void afficheVect(long int *A, char *var, int size)
{
	int i;
	printf("%s := ",var);
	
	for(i=0;i<size;i++){
		printf("%16.16lX ",A[	i]);
	}
	printf("\n");
}


// Rho = 2**47
static inline void init_nA_nB(long int * nA, long int * nB)
{

	for(int j=0; j<NB_COEFF;j++){
		nA[j] = ( (((int64_t)((rand()&0xffff)|(rand()<<16))<<32)) )| (((int64_t)((rand()&0xffff)|(rand()<<16))&0xffffffffUL));
		nB[j] = ( (((int64_t)((rand()&0xffff)|(rand()<<16))<<32)) )| (((int64_t)((rand()&0xffff)|(rand()<<16))&0xffffffffUL));
		
		// pour |nA| < phi et |nB| < phi
		nA[j] = nA[j]&CONV_MASK;
		nB[j] = nB[j]&CONV_MASK;
	}
	
}

/********************************************************************************
*
* MAIN
*
*********************************************************************************/

#define CEIL_DIVIDE(a, b)  (((a)/(b)) + ((a) % (b) == 0 ? 0 : 1)) /*!< Divide a by b and ceil the result*/
#define NB_COEFF_AVX512 CEIL_DIVIDE(NB_COEFF,8)

int main(int argc, char* argv[]){


	int flag=0, counter=0;

	static int512 nA[NB_COEFF_AVX512], nB[NB_COEFF_AVX512], res[NB_COEFF_AVX512], resMul[NB_COEFF_AVX512], tmp[NB_COEFF_AVX512];
	uint64_t mini = (uint64_t)-1L, mini1 = (uint64_t)-1L;

	mpz_t modul_p;
	mpz_t gama_pow[POLY_DEG];

	
	srand(time(NULL));
	
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	
	int nb_limbs;
	mpz_t A, B, C, E, F, R2, R3, G, H;
	mpz_inits (A, B, C, E, F, R2, R3, G, H, NULL);
	
	mp_limb_t mip0;
	
	const mp_limb_t *p_limbs;
	mp_limb_t *r2_limbs, *scratch_limbs, *mip_limbs;
	mp_limb_t *a_limbs, *b_limbs, *am1_limbs, *am2_limbs, *bm1_limbs, *bm2_limbs;
	
	init_data(modul_p,gama_pow);
	
	mpz_urandomm(A, r, modul_p);
	mpz_urandomm(B, r, modul_p);
	mpz_set(E, A);

	
	nb_limbs = mpz_size (modul_p);
	
	p_limbs = mpz_limbs_read (modul_p);   
	
	binvert_limb (mip0, p_limbs[0]);
	mip0 = -mip0;
	
	mpz_setbit (R2, 2*nb_limbs*8*sizeof(mp_limb_t)); 
	mpz_mod(R2, R2, modul_p);
	r2_limbs = mpz_limbs_modify (R2, nb_limbs); 
	
	a_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	am1_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	am2_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	bm1_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	bm2_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	mip_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	scratch_limbs = (mp_limb_t*) calloc (2*nb_limbs, sizeof(mp_limb_t));

	b_limbs = mpz_limbs_modify (B, nb_limbs);
	copy_limbs(a_limbs, A, nb_limbs);
	copy_limbs(am1_limbs, A, nb_limbs);
	copy_limbs(am2_limbs, A, nb_limbs);
	copy_limbs(bm1_limbs, B, nb_limbs);
	copy_limbs(bm2_limbs, B, nb_limbs);
	
	mpn_binvert (mip_limbs, p_limbs, nb_limbs, scratch_limbs);
	

	/***********************************************/
	printf("\nGMP mult_mod :\n----------\n");
	gmp_printf("p       : %Zd\n\n", modul_p);
	gmp_printf("A       : %Zd\n", E);
	gmp_printf("B       : %Zd\n\n", B);

	//~ conversion to Montgomery domain (block mont)
	mpn_mont_mul_red_1(am1_limbs, am1_limbs, r2_limbs, p_limbs, mip0, nb_limbs);
	mpn_mont_mul_red_1(bm1_limbs, bm1_limbs, r2_limbs, p_limbs, mip0, nb_limbs);
	
	//~ Montgomery modular multiplication (block mont)
	mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
	clean_limbs(scratch_limbs, 2*nb_limbs);
	for(int i=0; i<nb_limbs; i++)
		scratch_limbs[i] = am1_limbs[i];
	mpn_redc_1 (am1_limbs, scratch_limbs, p_limbs, nb_limbs, mip0);
	from_limbs_to_mpz_t(G, am1_limbs, nb_limbs);
	gmp_printf("r_mbgmp : %Zd\n", G);
	
	//~ conversion to Montgomery domain (classic mont)
	mpn_mont_mul_red_n(am2_limbs, am2_limbs, r2_limbs, p_limbs, mip_limbs, nb_limbs);
	mpn_mont_mul_red_n(bm2_limbs, bm2_limbs, r2_limbs, p_limbs, mip_limbs, nb_limbs);
	
	//~ Montgomery modular multiplication (block mont)
	mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
	
	clean_limbs(scratch_limbs, 2*nb_limbs);
	for(int i=0; i<nb_limbs; i++)
		scratch_limbs[i] = am2_limbs[i];
	mpn_redc_n (am2_limbs, scratch_limbs, p_limbs, nb_limbs, mip_limbs);
	from_limbs_to_mpz_t(H, am2_limbs, nb_limbs);
	gmp_printf("r_mcgmp : %Zd\n", H);
	
	
	//~ Montgomery low level modular multiplication
	mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
	from_limbs_to_mpz_t(F, a_limbs, nb_limbs);
	gmp_printf("r_lgmp  : %Zd\n", F);
	
#ifdef TEST
	goto chrono;
#endif
	
	
	
		mpz_mul (E, A, B);
		mpz_mod (E, E, modul_p);
	gmp_printf("r_gmp   : %Zd\n", E);


	/***********************************************/
	printf("\nmult_mod_poly_AVX512 n=%d, phi=2^52 :\n------------------------------------\n",NB_COEFF);
	
	from_int_to_pmns(nA->i64, A);
	from_int_to_pmns(nB->i64, B);
	gmp_printf("A       : %Zd\n", A);
	gmp_printf("B       : %Zd\n\n", B);
	afficheVect(nA->i64,"nA",NB_COEFF);
	afficheVect(nB->i64,"nB",NB_COEFF);


	//goto fin;


	mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
		
	afficheVect(res->i64,"res (AVX512)",NB_COEFF);
	/*goto fin;//*/
		
	printf("\nComparison with mult_mod_poly n=%d, phi=2^52 :\n---------------------------------------------\n",NB_COEFF);
	
	mult_mod_poly(resMul->i64, nA->i64,nB->i64);

	afficheVect(resMul->i64,"resMul    ",NB_COEFF);
	
	for(int i=0; i<NB_COEFF;i++) tmp[i>>3].i64[i&0x7] = res[i>>3].i64[i&0x7]^resMul[i>>3].i64[i&0x7];
	
	printf("\n");
	afficheVect(tmp->i64,"cmp",NB_COEFF);
	
	printf("\n");
	
	for(int i=0; i<NB_COEFF;i++)
		if(res[i>>3].i64[i&0x7]^resMul[i>>3].i64[i&0x7]) flag++;

	printf("flag = %d ; ",flag);
	
	if(!flag) printf("Victory !!!!!!!!!!!!!!!!!!!\n\n");
	else printf("Too bad !!!!!!!!!!!!!!!!!!!\n\n");
	
	//~ should not modify the value which is represented
	//exact_coeffs_reduction(resMul->i64, resMul->i64);
	from_pmns_to_int(C, resMul->i64, modul_p, gama_pow);
	
	gmp_printf("r_pmns  : %Zd\n", C);
	
	flag=0;
	printf("\n\n");
	
	//goto fin;//*/

	printf("\nComparison with mult_mod_poly_n%d :\n----------------------------------\n",NB_COEFF);
	
	
	
	printf("\nAll results :\n-------------\n");
	
	gmp_printf("p       : %Zd\n\n", modul_p);
	gmp_printf("A       : %Zd\n", A);
	gmp_printf("B       : %Zd\n\n", B);
	gmp_printf("r_pmns (n56, AVX512) : %Zd\n", C);
	
	printf("\nGMP mult_mod :\n----------\n");
	gmp_printf("r_mbgmp             : %Zd\n", G);
	gmp_printf("r_mcgmp             : %Zd\n", H);
	gmp_printf("r_lgmp              : %Zd\n", F);
	gmp_printf("r_gmp               : %Zd\n\n\n", E);
	//goto fin;//*/

chrono:
	
	
	printf("\t  /****************************/\n");
	printf("\t /   Test on 1000 datasets    /\n");
	printf("\t/****************************/\n\n");

	
	for(int i=0;i<NTEST;i++)
	{

		init_nA_nB(nA->i64,nB->i64);
		

		mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);

		//mult_mod_poly_n35_pmns(res->i64,nA->i64,nB->i64);
		
		mult_mod_poly(resMul->i64, nA->i64,nB->i64);
		
		for(int i=0; i<NB_COEFF;i++)
			if(res[i>>3].i64[i&0x7]^resMul[i>>3].i64[i&0x7]) flag++;
		flag?counter++,flag=0:counter,flag=0;
	
	}
	if(counter) printf("%d errors !\nToo bad !!!!!!!!!!!!!!!!!!!\n\n",counter),counter=0;
	else printf("mult_mod_poly_AVX512 vs mult_mod_poly: Victory !!!!!!!!!!!!!!!!!!!\n\n");
	counter=0;


	
	/********************
	//Timings  !!!!!!!!!!!!!
	*********************/
	
	#define K_T 10
	
	
	
	printf("\t  /*********************/\n");
	printf("\t / Timings !!!!!!!!!!!!/\n");
	printf("\t/*********************/\n\n");

 
	unsigned long long int timer=0, timer1=0, timer4=0, timer5=0, timer6=0, timer7=0;

	
	// timer mult_mod_poly        :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini1 = (uint64_t)-1L;

		init_nA_nB(nA->i64,nB->i64);
		
		for(int i=0;i<NTEST;i++)
		{
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			

			STAMP(START1)
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);//*/
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			STAMP(STOP1)
			

			if(mini1>STOP1-START1) mini1 = STOP1-START1;
			
		}
		
		timer += mini1;
	}

	printf("timer mult_mod_poly                : %llu\n",timer/(K_T*NSAMPLES));


		
	// timer mult_mod_poly AVX512 :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;

		init_nA_nB(nA->i64,nB->i64);
		
		for(int i=0;i<NTEST;i++)
		{
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			STAMP(STOP)
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer1 += mini;
	}
	
	printf("timer mult_mod_poly AVX512         : %llu\n",timer1/(K_T*NSAMPLES));//*/	
	
	



	
	// timer gmp mpz                      :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			STAMP(STOP)
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer4 += mini;
	}
	
	printf("timer gmp mpz                      : %llu\n",timer4/(K_T*NSAMPLES));

	
	// timer gmp low level                :
		
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;

		for(int i=0;i<NTEST;i++)
		{
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			STAMP(STOP)
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer5 += mini;
	}
	
	printf("timer gmp low level                : %llu\n",timer5/(K_T*NSAMPLES));


	// timer gmp Montgomery               :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			STAMP(STOP)
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer6 += mini;
	}
	
	printf("timer gmp Montgomery               : %llu\n",timer6/(K_T*NSAMPLES));


	// timer gmp Montgomery bloc          :

	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			STAMP(STOP)
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer7 += mini;
	}
	printf("timer gmp Montgomery bloc          : %llu\n",timer7/(K_T*NSAMPLES));
	
/*********************************************************

	Instruction counting with rdpmc

*********************************************************/

	printf("\n\t  /**************************/\n");
	printf("\t / Instructions !!!!!!!!!!!!/\n");
	printf("\t/**************************/\n\n");
	
	
	timer=0, timer1=0, timer4=0, timer5=0, timer6=0, timer7=0;

	
	// #Instructions mult_mod_poly_n6_pmns        :
	
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini1 = (uint64_t)-1L;

		init_nA_nB(nA->i64,nB->i64);
		
		for(int i=0;i<NTEST;i++)
		{
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START1 = rdpmc_instructions();
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);//*/
			mult_mod_poly(resMul->i64, nA->i64,nB->i64);
			STOP1 = rdpmc_instructions();

			if(mini1>STOP1-START1) mini1 = STOP1-START1;
			
		}
		
		timer += mini1;
	}
		
	printf("#Instructions mult_mod_poly                : %llu\n",timer/(K_T*NSAMPLES));


	// #Instructions mult_mod_poly_N25_pmns AVX512 :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;

		init_nA_nB(nA->i64,nB->i64);
		
		for(int i=0;i<NTEST;i++)
		{
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			mult_mod_poly_AVX512(res->i64,nA->i64,nB->i64);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer1 += mini;
	}
	
	printf("#Instructions mult_mod_poly AVX512         : %llu\n",timer1/(K_T*NSAMPLES));//*/	
	
	


	
	// #Instructions gmp mpz                      :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer4 += mini;
	}
	
	printf("#Instructions gmp mpz                      : %llu\n",timer4/(K_T*NSAMPLES));

	
	// #Instructions gmp low level                :
		
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;

		for(int i=0;i<NTEST;i++)
		{
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer5 += mini;
	}
	
	printf("#Instructions gmp low level                : %llu\n",timer5/(K_T*NSAMPLES));


	// #Instructions gmp Montgomery               :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer6 += mini;
	}
	
	printf("#Instructions gmp Montgomery               : %llu\n",timer6/(K_T*NSAMPLES));


	// #Instructions gmp Montgomery bloc          :

	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer7 += mini;
	}
	
	printf("#Instructions gmp Montgomery bloc          : %llu\n",timer7/(K_T*NSAMPLES));


fin:
	free(a_limbs);
	free(am1_limbs);
	free(am2_limbs);
	free(bm1_limbs);
	free(bm2_limbs);
	free(mip_limbs);
	free_data(modul_p,gama_pow);

	printf("\n");
	printf("Au revoir et merci !\n\n");

}
