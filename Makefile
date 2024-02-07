CC     := g++
CFLAGS := -fopenmp -O3 -fPIC -fno-loop-optimize -fno-aggressive-loop-optimizations -lgfortran

DIR1   := src/c_gkmPWMlasso
SRC1   := $(wildcard src/c_gkmPWMlasso/*.c)
OBJ1   := $(patsubst %.c,%.o,$(SRC1))
LIB1   := $(wildcard src/*.a)

DIR2   := src/c_gkmPWM
SRC2   := $(wildcard src/c_gkmPWM/*.c)
OBJ2   := $(patsubst %.c,%.o,$(SRC2))
LIB2   := $(wildcard src/*.a)

DIR3   := src/c_mapTF
SRC3   := $(wildcard src/c_mapTF/*.c)
OBJ3   := $(patsubst %.c,%.o,$(SRC3))
LIB3   := $(wildcard src/*.a)

DIR4   := src/c_getgkmweights
SRC4   := $(wildcard src/c_getgkmweights/*.c)
OBJ4   := $(patsubst %.c,%.o,$(SRC4))
LIB4   := $(wildcard src/*.a)

.PHONY: all clean gkmPWMlasso gkmPWM mapTF getgkmweights

all: gkmPWMlasso gkmPWM mapTF getgkmweights

$(DIR1)/%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@  
$(DIR2)/%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@
$(DIR3)/%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@
$(DIR4)/%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@
	
gkmPWMlasso: $(OBJ1) $(LIB1)
	$(CC) $(OBJ1) $(LIB1) $(CFLAGS) -o gkmPWMlasso
gkmPWM     : $(OBJ2) $(LIB2)
	$(CC) $(OBJ2) $(LIB2) $(CFLAGS) -o gkmPWM
mapTF      : $(OBJ3) $(LIB3)
	$(CC) $(OBJ3) $(LIB3) $(CFLAGS) -o mapTF
getgkmweights      : $(OBJ4) $(LIB4)
	$(CC) $(OBJ4) $(LIB4) $(CFLAGS) -o getgkmweights
clean:
	$(RM) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) gkmPWMlasso gkmPWM mapTF getgkmweights
