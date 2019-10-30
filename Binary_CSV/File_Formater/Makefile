CXX = g++
FLAGS = -std=c++11

DEPS = fileFormat.h
OBJ = fileFormat.o main.o

%.0: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(FLAGS)
        
fileForamt: $(OBJ)
	$(CXX) -o $@ $^ $(FLAGS)
        
run:

clean:
	@rm *.o
	@rm fileFormat


