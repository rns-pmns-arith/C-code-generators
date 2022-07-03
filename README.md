# C-code-generators

The PMNS used in this work have been generated using the generator provided at [PMNS generator](https://github.com/arithPMNS/low_memory_efficient_PMNS/tree/main/pmns_generator). 

A standard C codes generator is provied at [C codes generator](https://github.com/arithPMNS/low_memory_efficient_PMNS/tree/main/C_codes_generator). Given a PMNS  (obtained from the preceding generator), it generates C codes for all the necessary operations (including forward and backward conversion to PMNS, addition, multiplication).

<br />
<br />
<br />

The directory 'AVX_C_gen' contains a C codes generator, which takes advantge of AVX512 instruction set extension.

 
 Note : All these generators have been implmented using SageMath library which can be found here: http://www.sagemath.org/
