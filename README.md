# C-code-generators

The PMNS used in this work have been geenrated using the geenrator provided at [PMNS generator](https://github.com/arithPMNS/low_memory_efficient_PMNS).

This repository gives two codes generators for PMNS. Given a PMNS parameter, obtened from the above generator:
 - 'standard_C_gen' directory allows to generate standard C codes, for all the necessary operations (including forward and backward conversion to PMNS, addition, multiplication);
 - 'AVX_C_gen' directory allows to generate C codes, which takes advantges of AVX512 instruction set extension.
 
 
 
 Note : To use the AMNS generator and the C code generator, you will need SageMath library which can be found here: http://www.sagemath.org/
