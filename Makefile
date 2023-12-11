CC     := g++
CFLAGS := -fopenmp -O3 -fPIC -fno-loop-optimize -fno-aggressive-loop-optimizations -lgfortran

DIR1   := src/c_gkmPWMlasso
SRC1   := $(wildcard src/c_gkmPWMlasso/*.c)
OBJ1   := $(patsubst %.c,%.o,$(SRC1))
LIB1   := $(wildcard src/c_gkmPWMlasso/*.a)

DIR2   := src/c_gkmPWM
SRC2   := $(wildcard src/c_gkmPWM/*.c)
OBJ2   := $(patsubst %.c,%.o,$(SRC2))
LIB2   := $(wildcard src/c_gkmPWM/*.a)

DIR3   := src/c_mapTF
SRC3   := $(wildcard src/c_mapTF/*.c)
OBJ3   := $(patsubst %.c,%.o,$(SRC3))
LIB3   := $(wildcard src/c_mapTF/*.a)


.PHONY: all clean gkmPWMlasso gkmPWM mapTF

all: gkmPWMlasso gkmPWM mapTF

$(DIR1)/%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@  
$(DIR2)/%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@
$(DIR3)/%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@

gkmPWMlasso: $(OBJ1) $(LIB1)
	$(CC) $(OBJ1) $(LIB1) $(CFLAGS) -o gkmPWMlasso
gkmPWM     : $(OBJ2) $(LIB2)
	$(CC) $(OBJ2) $(LIB2) $(CFLAGS) -o gkmPWM
mapTF      : $(OBJ3) $(LIB3)
	$(CC) $(OBJ3) $(LIB3) $(CFLAGS) -o mapTF
clean:
	$(RM) $(OBJ1) $(OBJ2) $(OBJ3) gkmPWMlasso gkmPWM mapTF
