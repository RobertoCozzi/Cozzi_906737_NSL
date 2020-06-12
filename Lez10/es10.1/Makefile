CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio10_1: esercizio10_1.o random.o individual.o population.o
	$(CC) random.o individual.o population.o esercizio10_1.o -o esercizio10_1 $(CFLAGS)
random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
individual.o: individual.cpp individual.h
	$(CC) -c individual.cpp -o individual.o $(CFLAGS)
population.0: population.cpp population.h
	$(CC) -c population.cpp -o population.o $(CFLAGS)
esercizio10_1.o: es10_1.cpp random.h individual.h population.h
	$(CC) -c es10_1.cpp -o esercizio10_1.o $(CFLAGS)