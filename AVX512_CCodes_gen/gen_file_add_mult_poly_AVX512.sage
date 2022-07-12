
'''def build_add_mult_poly_h_file(dir_path, small_int, big_int):
	with open(dir_path+"/add_mult_poly.h", "w") as f:
		
		f.write("#ifndef POLY_MULT_ADD\n")
		f.write("#define POLY_MULT_ADD\n\n\n")
		
		f.write("void sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void double_add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void double_sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void neg_poly(" + small_int + " *rop, " + small_int + " *op);\n")
		f.write("void scalar_mult_poly(" + small_int + " *rop, " + small_int + " *op, " + small_int + " scalar);\n")
		f.write("void double_poly_coeffs(" + small_int + " *rop, " + small_int + " *op);\n")
		f.write("void lshift_poly_coeffs(" + small_int + " *rop, " + small_int + " *op, int nb_pos);\n\n")
		
		f.write("void mult_mod_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n\n")
		f.write("void square_mod_poly(" + small_int + " *rop, " + small_int + " *pa);\n\n")	
		
		f.write("void internal_reduction(" + small_int + " *rop, " + big_int + " *op);\n\n")
		f.write("void exact_coeffs_reduction(" + small_int + " *rop, " + small_int + " *op);\n\n")
		
		f.write("#endif\n\n")


dir_path="c_codes"
n=11
upPow_distrib=[float(0)]#[Integer(0)]#'''
big_int="__m512i"

def build_add_mult_poly_cAVX512_file(dir_path, n, upPow_distrib):#(dir_path, n, mont_phi, upPow_distrib, small_int, unsigned_small_int, big_int, red_int_coeff, neg_inv_ri_rep_coeff, mask_for_redint):
	with open(dir_path+"/add_mult_poly_AVX512.c", "w") as f:
		
		#~ to check if it is an AMNS or not
		if (type(upPow_distrib[0]) == int) or (type(upPow_distrib[0]) == Integer) : 
			is_AMNS = True
			#upPow_distrib = upPow_distrib[0]
		else :
			is_AMNS = False


		#f.write("#define idx_a1 (__m512i){0x0UL,0x8UL,0x8UL,0x8UL,0x8UL,0x8UL,0x8UL,0x8UL}\n")

		f.write("#define zero512 (__m512i){0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL}\n")


		if is_AMNS :
			tmp="REDINT_MASK_N"+str(n)+"_AMNS"
			f.write("#define sh_one (__m128i){0x1UL,0UL}\n")
		else :
			f.write("#define idx_a0 (__m512i){0x1UL,0x2UL,0x3UL,0x4UL,0x5UL,0x6UL,0x7UL,0x8UL}\n")
			off=n&0x7
			if off!=0:
				f.write("#define idx_a1 (__m512i){")
				
				for i in range(off-1):
					f.write(hex(i+1)+'UL,')
				for i in range(off-1,7):
					f.write("0x8UL,")
				f.write("0x8UL}\n")
			else:
				f.write("#define idx_a1 (__m512i){0x1UL,0x2UL,0x3UL,0x4UL,0x5UL,0x6UL,0x7UL,0x8UL}\n")

			
			#tmp="REDINT_MASK_N"+str(n)+"_PMNS"
			tmp="REDINT_MASK"
		f.write("#define mask52 (__m512i){")
		for i in range(7):
			f.write(tmp+",")
		f.write(tmp+"}\n\n")
		f.write("#include \"avx_const_data.h\"\n\n\n")
		
		f.write("#define idx_0 (__m512i){0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL}\n")
		f.write("#define idx_1 (__m512i){0x1UL,0x1UL,0x1UL,0x1UL,0x1UL,0x1UL,0x1UL,0x1UL}\n")
		f.write("#define idx_2 (__m512i){0x2UL,0x2UL,0x2UL,0x2UL,0x2UL,0x2UL,0x2UL,0x2UL}\n")
		f.write("#define idx_3 (__m512i){0x3UL,0x3UL,0x3UL,0x3UL,0x3UL,0x3UL,0x3UL,0x3UL}\n")
		f.write("#define idx_4 (__m512i){0x4UL,0x4UL,0x4UL,0x4UL,0x4UL,0x4UL,0x4UL,0x4UL}\n")
		f.write("#define idx_5 (__m512i){0x5UL,0x5UL,0x5UL,0x5UL,0x5UL,0x5UL,0x5UL,0x5UL}\n")
		f.write("#define idx_6 (__m512i){0x6UL,0x6UL,0x6UL,0x6UL,0x6UL,0x6UL,0x6UL,0x6UL}\n")
		f.write("#define idx_7 (__m512i){0x7UL,0x7UL,0x7UL,0x7UL,0x7UL,0x7UL,0x7UL,0x7UL}\n\n")


		f.write("#define idx_opl (__m512i){0UL,2UL,4UL,6UL,0x8UL,0xAUL,0xCUL,0xEUL}\n")
		f.write("#define idx_oph (__m512i){0x1UL,0x3UL,0x5UL,0x7UL,0x9UL,0xBUL,0xDUL,0xFUL}\n")
		f.write("#define idx_l (__m512i){0x0UL,0x8UL,0x1UL,0x9UL,0x2UL,0xAUL,0x3UL,0xBUL}\n")
		f.write("#define idx_h (__m512i){0x4UL,0xCUL,0x5UL,0xDUL,0x6UL,0xEUL,0x7UL,0xFUL}\n\n")

		if is_AMNS :
			f.write(idx_b_code(n))
		else:
			f.write(idx_c_code(n))

		print('fin des macros')

		f.write("void mult_mod_poly_AVX512(int64_t *rop, int64_t *pa, int64_t *pb){\n\n")
		
		n512=ceil(n/8)
		if is_AMNS :
			prod_funct = build_prod_code_AMNS_AVX512
			#sqr_funct = build_square_code_v2
				
		else :
			prod_funct = build_prod_code_PMNS_AVX512
			#sqr_funct = build_square_code
				
		
		f.write(prod_funct(n))
		
		f.write('\n#define idx_opl (__m512i){0UL,2UL,4UL,6UL,0x8UL,0xAUL,0xCUL,0xEUL}\n')
		f.write('#define idx_oph (__m512i){0x1UL,0x3UL,0x5UL,0x7UL,0x9UL,0xBUL,0xDUL,0xFUL}\n\n')
	
		f.write("inline void internal_reduction_512(int64_t *rop, int128 *op){\n\n")
		
		f.write(int_red_funct_AVX512(n))



def int_red_funct_AVX512(n):

	n512=ceil(n/8)
	tmp = '\t__m512i carry;\n'
	tmp += '\n\t__m512i '
	
	for i in range(n512-1):
		tmp+='acc512l'+str(i)+' = zero512, '
		
	tmp+='acc512l'+str(n512-1)+' = zero512'
		
		
	tmp+=' ;\n\n'
	
	tmp += '\t__m512i '
	
	for i in range(n512-1):
		tmp+='acc512h'+str(i)+' = zero512, '
		
	tmp+='acc512h'+str(n512-1)+' = zero512'
		
		
	tmp+=' ;\n\n\n'
	tmp += '\t//Computation of Q\n'
	tmp+= '\t//~ computation of : op*neg_inv_ri_rep_coeff mod(E, mont_phi)\n'
	
	
	tmp+='\n\t__m512i opl512 = _mm512_set1_epi64((uint64_t)op[0]);\n'
	
	for i in range(n512):
		tmp+='\t__m512i q512_'+str(i)+' = _mm512_madd52lo_epu64(zero512,opl512,KQ_0_'+str(i)+');\n'
	
	for i in range(1,n):
		ii = i>>3
		tmp+='\n\topl512 = _mm512_set1_epi64((uint64_t)op['+str(i)+']);\n'
		for j in range(n512):
			tmp+='\tq512_'+str(j)+' = _mm512_madd52lo_epu64(q512_'+str(j)+',opl512,KQ_'+str(i)+'_'+str(j)+');\n'
	
	
	print('fin computation of Q, début computation of : Q*red_int_coeff mod(E)')
	
	tmp+='\n\t//~ computation of : Q*red_int_coeff mod(E)\n'
	
	tmp+='\n\t__m512i tmpq512 = _mm512_permutexvar_epi64(idx_0,q512_0);\n'
	for i in range(n512):
		tmp+='\t__m512i tmpZero512l'+str(i)+' = _mm512_madd52lo_epu64(acc512l'+str(i)+',KZ_0_'+str(i)+',tmpq512);\n'
	for i in range(n512):
		tmp+='\t__m512i tmpZero512h'+str(i)+' = _mm512_madd52hi_epu64(acc512h'+str(i)+',KZ_0_'+str(i)+',tmpq512);\n'

	for i in range(1,n):
		ii = i>>3
		tmp+='\n\ttmpq512 = _mm512_permutexvar_epi64(idx_'+str(i%8)+',q512_'+str(ii)+');\n'	
		for j in range(n512):
			tmp+='\ttmpZero512l'+str(j)+' = _mm512_madd52lo_epu64(tmpZero512l'+str(j)+',KZ_'+str(i)+'_'+str(j)+',tmpq512);\n'
		for j in range(n512):
			tmp+='\ttmpZero512h'+str(j)+' = _mm512_madd52hi_epu64(tmpZero512h'+str(j)+',KZ_'+str(i)+'_'+str(j)+',tmpq512);\n'

	
	print('fin computation of : Q*red_int_coeff mod(E), final reconstruction')
	
	tmp+='\n\t// Final reconstruction\n'

	tmp+='\n\t// Final reconstruction\n'
	
	for j in range(n512):	
		tmp+='\tcarry = _mm512_srli_epi64(tmpZero512l'+str(j)+',52);\n'
		tmp+='\ttmpZero512h'+str(j)+' = _mm512_add_epi64(carry,tmpZero512h'+str(j)+');\n'
	
	tmp+='\n\t_mm512_store_epi64(rop,tmpZero512h0);\n'
		
	for j in range(1,n512):	
		tmp+='\t_mm512_store_epi64(rop+'+str(j*8)+',tmpZero512h'+str(j)+');\n'
		
	tmp += '\n\n}\n'
	
	print('fin internal_reduction_512()')

	return tmp
















	

def build_prod_code_AMNS_AVX512(n):
	
	n512=ceil(n/8)
	tmp = '\n\t__m512i tmp512_a;\n'
		

	tmp+="\t__m512i pb512[NB_COEFF_N"+str(n)+"_AVX512];\n\n"
	tmp+="\tpb512[0] = _mm512_loadu_epi64(pa);\n"
	for i in range(1,n512):
		tmp+="\tpb512["+str(i)+"] = _mm512_loadu_epi64(pa+"+str(i*8)+");\n"
	#print(n)
	
	tmp += '\n\t__m512i '
	
	for i in range(n512-1):
		tmp+='acc512l'+str(i)+' = zero512 '
		
	tmp+='acc512l'+str(n512-1)+' = zero512'
		
		
	tmp+=' ;\n\n'
	
	tmp += '\t__m512i '
	
	for i in range(n512-1):
		tmp+='acc512h'+str(i)+' = zero512, '
		
	tmp+='acc512h'+str(n512-1)+' = zero512'
		
		
	tmp+=' ;\n'
	
	
	for i in range(n):
		#print(i)
		tmp+='\n\ttmp512_a = _mm512_set1_epi64(pa['+str(i)+']);\n'
		for j in range(n512-1):
			tmp+='\ttmp512_'+str(j)+' = _mm512_permutex2var_epi64(pa512['+str(j)+'],idx_a0,pa512['+str(j+1)+']);\n'
		tmp+= '\tpa512['+str(n512-1)+'] = _mm512_permutex2var_epi64(pa512['+str(n512-1)+'],idx_a0,pa512[0]);\n'
		for j in range(n512-1):
			tmp+= '\tpa512['+str(j)+'] = tmp512_'+str(j)+' ;\n'
		
		
		for j in range(n512):
			tmp+='\n\tacc512l'+str(j)+' = _mm512_madd52lo_epu64(acc512l'+str(j)+',pa512['+str(j)+'],tmp512_b);\n'
			tmp+='\tacc512h'+str(j)+' = _mm512_madd52hi_epu64(acc512h'+str(j)+',pa512['+str(j)+'],tmp512_b);\n'
		
		if i!=0:
			ii = (i>>3)
			if i%8==0:
				tmp+='\n\tcx512l'+str(ii)+' = _mm512_permutex2var_epi64(cx512l'+str(ii)+',idx_c'+str(i+n-1)+',acc512l'+str(ii-1)+');\n'
				tmp+='\tcx512h'+str(ii)+' = _mm512_permutex2var_epi64(cx512h'+str(ii)+',idx_c'+str(i+n-1)+',acc512h'+str(ii-1)+');\n'
			else :
				tmp+='\n\tcx512l'+str(ii)+' = _mm512_permutex2var_epi64(cx512l'+str(ii)+',idx_c'+str(i+n-1)+',acc512l'+str(ii)+');\n'
				tmp+='\tcx512h'+str(ii)+' = _mm512_permutex2var_epi64(cx512h'+str(ii)+',idx_c'+str(i+n-1)+',acc512h'+str(ii)+');\n'
			
	
	tmp+='\n\t// final carries\n'

	for j in range(n512):
			tmp+='\n\tacc512l'+str(j)+' = _mm512_add_epi64(acc512l'+str(j)+',cx512l'+str(j)+');\n'
			tmp+='\tacc512h'+str(j)+' = _mm512_add_epi64(acc512h'+str(j)+',cx512h'+str(j)+');\n'
		
	tmp+='\n\t// internal reduction\n'
	
	print('fin multiplication mod E(X), début réduction interne')
	
	
	tmp+='\n\t//Computation of Q\n'
	tmp+='\t//~ computation of : op*neg_inv_ri_rep_coeff mod(E, mont_phi)\n'
	
	tmp+='\n\t__m512i opl512 = _mm512_permutexvar_epi64(idx_0,acc512l0);\n'
	
	for i in range(n512):
		tmp+='\t__m512i q512_'+str(i)+' = _mm512_madd52lo_epu64(zero512,opl512,KQ_0_'+str(i)+');\n'
	
	for i in range(1,n):
		ii = i>>3
		tmp+='\n\topl512 = _mm512_permutexvar_epi64(idx_'+str(i%8)+',acc512l'+str(ii)+');\n'
		for j in range(n512):
			tmp+='\tq512_'+str(j)+' = _mm512_madd52lo_epu64(q512_'+str(j)+',opl512,KQ_'+str(i)+'_'+str(j)+');\n'
	
	
	print('fin computation of Q, début computation of : Q*red_int_coeff mod(E)')
	
	tmp+='\n\t//~ computation of : Q*red_int_coeff mod(E)\n'
	
	tmp+='\n\t__m512i tmpq512 = _mm512_permutexvar_epi64(idx_0,q512_0);\n'
	for i in range(n512):
		tmp+='\t__m512i tmpZero512l'+str(i)+' = _mm512_madd52lo_epu64(acc512l'+str(i)+',KZ_0_'+str(i)+',tmpq512);\n'
	for i in range(n512):
		tmp+='\t__m512i tmpZero512h'+str(i)+' = _mm512_madd52hi_epu64(acc512h'+str(i)+',KZ_0_'+str(i)+',tmpq512);\n'

	for i in range(1,n):
		ii = i>>3
		tmp+='\n\ttmpq512 = _mm512_permutexvar_epi64(idx_'+str(i%8)+',q512_'+str(ii)+');\n'	
		for j in range(n512):
			tmp+='\ttmpZero512l'+str(j)+' = _mm512_madd52lo_epu64(tmpZero512l'+str(j)+',KZ_'+str(i)+'_'+str(j)+',tmpq512);\n'
		for j in range(n512):
			tmp+='\ttmpZero512h'+str(j)+' = _mm512_madd52hi_epu64(tmpZero512h'+str(j)+',KZ_'+str(i)+'_'+str(j)+',tmpq512);\n'

	
	print('fin computation of : Q*red_int_coeff mod(E), final reconstruction')
	
	tmp+='\n\t// Final reconstruction\n'
	
	for j in range(n512):	
		tmp+='\tcarry = _mm512_srli_epi64(tmpZero512l'+str(j)+',52);\n'
		tmp+='\ttmpZero512h'+str(j)+' = _mm512_add_epi64(carry,tmpZero512h'+str(j)+');\n'
	
	tmp+='\n\t_mm512_store_epi64(rop,tmpZero512h0);\n'
		
	for j in range(1,n512):	
		tmp+='\t_mm512_store_epi64(rop+'+str(j*8)+',tmpZero512h'+str(j)+');\n'
		
	tmp += '\n\n}\n'
	
	print('fin mult_mod_poly_AVX512()')

	return tmp













def build_prod_code_PMNS_AVX512(n):
	
	n512=ceil(n/8)
	tmp = '\n\t__m512i tmp512_b, carry;\n'
		

	#tmp+="\t__m512i pa512[NB_COEFF_N"+str(n)+"_AVX512];\n\n"
	tmp+="\t__m512i pa512[NB_COEFF];\n\n"
	tmp+="\tpa512[0] = _mm512_loadu_epi64(pa);\n"
	for i in range(1,n512):
		tmp+="\tpa512["+str(i)+"] = _mm512_loadu_epi64(pa+"+str(i*8)+");\n"
	tmp+="\t__m512i tmp512_0"
	for i in range(1,n512-1):
		tmp+=", tmp512_"+str(i)

	tmp+=' ;\n\n'
	#print(n)
	
	tmp += '\n\t__m512i '
	
	for i in range(n512-1):
		tmp+='acc512l'+str(i)+' = zero512, cx512l'+str(i)+' = zero512, '
		
	tmp+='acc512l'+str(n512-1)+' = zero512, cx512l'+str(n512-1)+' = zero512'
		
		
	tmp+=' ;\n\n'
	
	tmp += '\t__m512i '
	
	for i in range(n512-1):
		tmp+='acc512h'+str(i)+' = zero512, cx512h'+str(i)+' = zero512, '
		
	tmp+='acc512h'+str(n512-1)+' = zero512, cx512h'+str(n512-1)+' = zero512'
		
		
	tmp+=' ;\n'
	
	
	for i in range(n-1,-1,-1):
		#print(i)
		tmp+='\n\ttmp512_b = _mm512_set1_epi64(pb['+str(i)+']);\n'
		for j in range(n512-1):
			tmp+='\ttmp512_'+str(j)+' = _mm512_permutex2var_epi64(pa512['+str(j)+'],idx_a0,pa512['+str(j+1)+']);\n'
		tmp+= '\tpa512['+str(n512-1)+'] = _mm512_permutex2var_epi64(pa512['+str(n512-1)+'],idx_a1,pa512[0]);\n'
		for j in range(n512-1):
			tmp+= '\tpa512['+str(j)+'] = tmp512_'+str(j)+' ;\n'
		
		
		for j in range(n512):
			tmp+='\n\tacc512l'+str(j)+' = _mm512_madd52lo_epu64(acc512l'+str(j)+',pa512['+str(j)+'],tmp512_b);\n'
			tmp+='\tacc512h'+str(j)+' = _mm512_madd52hi_epu64(acc512h'+str(j)+',pa512['+str(j)+'],tmp512_b);\n'
		
		if i!=0:
			ii = (i>>3)
			if i%8==0:
				tmp+='\n\tcx512l'+str(ii)+' = _mm512_permutex2var_epi64(cx512l'+str(ii)+',idx_c'+str(i+n-1)+',acc512l'+str(ii-1)+');\n'
				tmp+='\tcx512h'+str(ii)+' = _mm512_permutex2var_epi64(cx512h'+str(ii)+',idx_c'+str(i+n-1)+',acc512h'+str(ii-1)+');\n'
			else :
				tmp+='\n\tcx512l'+str(ii)+' = _mm512_permutex2var_epi64(cx512l'+str(ii)+',idx_c'+str(i+n-1)+',acc512l'+str(ii)+');\n'
				tmp+='\tcx512h'+str(ii)+' = _mm512_permutex2var_epi64(cx512h'+str(ii)+',idx_c'+str(i+n-1)+',acc512h'+str(ii)+');\n'
			
	
	tmp+='\n\t// final carries\n'

	for j in range(n512):
			tmp+='\n\tacc512l'+str(j)+' = _mm512_add_epi64(acc512l'+str(j)+',cx512l'+str(j)+');\n'
			tmp+='\tacc512h'+str(j)+' = _mm512_add_epi64(acc512h'+str(j)+',cx512h'+str(j)+');\n'
		
	tmp+='\n\t// internal reduction\n'
	
	print('fin multiplication mod E(X), début réduction interne')
	
	
	tmp+='\n\t//Computation of Q\n'
	tmp+='\t//~ computation of : op*neg_inv_ri_rep_coeff mod(E, mont_phi)\n'
	
	tmp+='\n\t__m512i opl512 = _mm512_permutexvar_epi64(idx_0,acc512l0);\n'
	
	for i in range(n512):
		tmp+='\t__m512i q512_'+str(i)+' = _mm512_madd52lo_epu64(zero512,opl512,KQ_0_'+str(i)+');\n'
	
	for i in range(1,n):
		ii = i>>3
		tmp+='\n\topl512 = _mm512_permutexvar_epi64(idx_'+str(i%8)+',acc512l'+str(ii)+');\n'
		for j in range(n512):
			tmp+='\tq512_'+str(j)+' = _mm512_madd52lo_epu64(q512_'+str(j)+',opl512,KQ_'+str(i)+'_'+str(j)+');\n'
	
	
	print('fin computation of Q, début computation of : Q*red_int_coeff mod(E)')
	
	tmp+='\n\t//~ computation of : Q*red_int_coeff mod(E)\n'
	
	tmp+='\n\t__m512i tmpq512 = _mm512_permutexvar_epi64(idx_0,q512_0);\n'
	for i in range(n512):
		tmp+='\t__m512i tmpZero512l'+str(i)+' = _mm512_madd52lo_epu64(acc512l'+str(i)+',KZ_0_'+str(i)+',tmpq512);\n'
	for i in range(n512):
		tmp+='\t__m512i tmpZero512h'+str(i)+' = _mm512_madd52hi_epu64(acc512h'+str(i)+',KZ_0_'+str(i)+',tmpq512);\n'

	for i in range(1,n):
		ii = i>>3
		tmp+='\n\ttmpq512 = _mm512_permutexvar_epi64(idx_'+str(i%8)+',q512_'+str(ii)+');\n'	
		for j in range(n512):
			tmp+='\ttmpZero512l'+str(j)+' = _mm512_madd52lo_epu64(tmpZero512l'+str(j)+',KZ_'+str(i)+'_'+str(j)+',tmpq512);\n'
		for j in range(n512):
			tmp+='\ttmpZero512h'+str(j)+' = _mm512_madd52hi_epu64(tmpZero512h'+str(j)+',KZ_'+str(i)+'_'+str(j)+',tmpq512);\n'

	
	print('fin computation of : Q*red_int_coeff mod(E), final reconstruction')
	
	tmp+='\n\t// Final reconstruction\n'
	
	for j in range(n512):	
		tmp+='\tcarry = _mm512_srli_epi64(tmpZero512l'+str(j)+',52);\n'
		tmp+='\ttmpZero512h'+str(j)+' = _mm512_add_epi64(carry,tmpZero512h'+str(j)+');\n'
	
	tmp+='\n\t_mm512_store_epi64(rop,tmpZero512h0);\n'
		
	for j in range(1,n512):	
		tmp+='\t_mm512_store_epi64(rop+'+str(j*8)+',tmpZero512h'+str(j)+');\n'
		
	tmp += '\n\n}\n'
	
	print('fin mult_mod_poly_AVX512()')


	return tmp


# ne marche que pour E = X^n-X-1
def idx_c_code(n):
	
	result=""
	n512=ceil(n/8)
	for i in range(n512):
		if n+i*8>2*n-2:
			break
		result += "#define idx_c" + str(n+i*8) + " (__m512i){0x0UL,0x8UL,0x2UL,0x3UL,0x4UL,0x5UL,0x6UL,0x7UL}\n"
		if n+i*8+1>2*n-2:
			break
		result += "#define idx_c" + str(n+i*8+1) + " (__m512i){0x0UL,0x0UL,0x9UL,0x3UL,0x4UL,0x5UL,0x6UL,0x7UL}\n"
		if n+i*8+2>2*n-2:
			break
		result += "#define idx_c" + str(n+i*8+2) + " (__m512i){0x0UL,0x0UL,0x0UL,0xAUL,0x4UL,0x5UL,0x6UL,0x7UL}\n"
		if n+i*8+3>2*n-2:
			break
		result += "#define idx_c" + str(n+i*8+3) + " (__m512i){0x0UL,0x0UL,0x0UL,0x0UL,0xBUL,0x5UL,0x6UL,0x7UL}\n"
		if n+i*8+4>2*n-2:
			break
		result += "#define idx_c" + str(n+i*8+4) + " (__m512i){0x0UL,0x0UL,0x0UL,0x0UL,0x0UL,0xCUL,0x6UL,0x7UL}\n"
		if n+i*8+5>2*n-2:
			break
		result += "#define idx_c" + str(n+i*8+5) + " (__m512i){0x0UL,0x0UL,0x0UL,0x0UL,0x0UL,0x0UL,0xDUL,0x7UL}\n"
		if n+i*8+6>2*n-2:
			break
		result += "#define idx_c" + str(n+i*8+6) + " (__m512i){0x0UL,0x0UL,0x0UL,0x0UL,0x0UL,0x0UL,0x0UL,0xEUL}\n"
		if n+i*8+7>2*n-2:
			break
		result += "#define idx_c" + str(n+i*8+7) + " (__m512i){0xFUL,0x1UL,0x2UL,0x3UL,0x4UL,0x5UL,0x6UL,0x7UL}\n"
			
			
	result +="\n"	
			
			
	
	return result
	
	
# ne marche que pour E = X^n-2
def idx_b_code(n):
	result="#define idx_b0 (__m512i){"
	
	off = n&0x7
	print("off =",off)
	
	result += hex(8+off-1)+'UL'
	for i in range(1,8):
		result +=","+hex(i-1)+'UL'
		
	result += '}\n'

	result="#define idx_bi (__m512i){"
	
	print("off =",off)
	
	result += hex(8+off-1)+'UL'
	for i in range(1,8):
		result +=","+hex(i-1)+'UL'
		
	result += '}\n'

	result +="#define idx_b1 (__m512i){0x7UL,0x8UL,0x9UL,0xAUL,0x0UL,0x0UL,0x0UL,0x0UL}\n"
				
	result +="\n"	
			
			
	
	return result
	

