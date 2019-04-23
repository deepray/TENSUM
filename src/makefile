include makefile.in

HDR = $(wildcard *.h)
SRC = $(wildcard *.cc)
SRC += ../extern/fparser/fparser.cc
OBJ = $(patsubst %.cc,%.o,$(SRC))
INC = -I../extern/fparser


TARGET = tensum

all: $(TARGET)
	

$(TARGET): $(OBJ) 
	$(MPICC) -o $(TARGET) $(OBJ) $(LDFLAGS)

$(OBJ): $(HDR)

.cc.o: $(HDR)
	$(MPICC) $(CFLAGS) $(INC) -c $*.cc -o $*.o

clean:
	rm -f $(OBJ) $(TARGET)


allclean:
	make clean

$(info )
