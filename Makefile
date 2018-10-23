TARGET = edvards_curve
CC = gcc
CCFLAGS = -Wall -std=gnu99 -pedantic-errors
OFLAGS = -c 

all: main.o montarith.o edvardscurve.o 
	$(CC) $(CCFLAGS) main.o montarith.o edvardscurve.o -lgmp -o $(TARGET)

main.o: main.c
	$(CC) $(OFLAGS) $(CCFLAGS) main.c -o main.o

montarith.o: montarith.c
	$(CC) $(CC) $(OFLAGS) $(CCFLAGS) montarith.c -o monarith.o

edvards.o: edvardscurve.c
	$(CC) $(CC) $(OFLAGS) $(CCFLAGS) edvardscurve.c -o edvardscurve.o

clean:
	rm -rf *.o $(TARGET)



