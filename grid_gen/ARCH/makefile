include makefile.in

HDR = $(wildcard *.h)
SRC = $(wildcard *.cc)
OBJ = $(patsubst %.cc,%.o,$(SRC))


TARGET = grid_part

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) -o $(TARGET) $(OBJ) $(LDFLAGS)

$(OBJ): $(HDR)

.cc.o: $(HDR)
	$(CXX) $(CFLAGS) -c $*.cc -o $*.o

clean:
	rm -f $(OBJ) $(TARGET)


allclean:
	make clean
