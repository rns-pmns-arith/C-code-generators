OPTIONS= -O3 -Wall -g -lgmp -lcrypto -mavx512f -mavx512dq -mavx512vl -mavx512ifma

CC= gcc-10

OPTIONS+=$(MACROS)

main : main.o add_mult_poly.o
	$(CC) -o main main.o add_mult_poly.o $(OPTIONS)

add_mult_poly.o : add_mult_poly.c add_mult_poly.h 
	$(CC) -c add_mult_poly.c $(OPTIONS)


main.o : add_mult_poly.o
	$(CC) -c main.c $(OPTIONS)


	
all : main


clean :
	rm *.o -fv *~

cleanall : clean
	rm main -f

