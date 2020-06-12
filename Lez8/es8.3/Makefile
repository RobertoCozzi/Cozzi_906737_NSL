LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
 
%.o : %.cpp
	g++ -Wall -c $< ${INCS}

esercizio8_3: es8_3.o
	g++ -O3 -Wall -o $@ $^ ${LIBS}

clean:
	rm *.o es8_3 potential.dat kinetic.dat probability.dat
