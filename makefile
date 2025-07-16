CC=g++
CPPFLAGS=-I.
LDFLAGS=-g -fopenmp -lsqlite3
DEPS = bigraph.h utility.h biclique.h
OBJ = bigraph.o main.o utility.o biclique.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y -c -O3 -o $@ $< $(CPPFLAGS) $(LDFLAGS)  

biclique: $(OBJ)
	$(CC) -std=c++1y -O3 -pthread -o $@ $^ $(CPPFLAGS) $(LDFLAGS) -lgomp 

clean:
	-rm -f biclique *.o
