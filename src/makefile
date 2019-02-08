SRC    = shallow_water.cpp
   
OBJ    = shallow_water.o
   
CFLAGS =  -std=c++11 

CC=g++ 

sw.x:	$(OBJ)
	$(CC) -o sw.x $(OBJ) $(LIBS)

shallow_water.o: shallow_water.cpp
	$(CC) $(CFLAGS) -c shallow_water.cpp

clean:
	rm -f shallow_water.o
	rm -f sw.x
