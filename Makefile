CC     := gcc
CFLAGS := -fopenmp -O3 -fPIC -fno-loop-optimize -fno-aggressive-loop-optimizations

DIR1   := c/gkmPWMlasso
SRC1   := $(wildcard c/gkmPWMlasso/*.c)
OBJ1   := $(patsubst %.c,%.o,$(SRC1))
LIB1   := $(wildcard c/*.a)

DIR2   := c/gkmPWM
SRC2   := $(wildcard c/gkmPWM/*.c)
OBJ2   := $(patsubst %.c,%.o,$(SRC2))
LIB2   := $(wildcard c/*.a)

DIR3   := c/mapTF
SRC3   := $(wildcard c/mapTF/*.c)
OBJ3   := $(patsubst %.c,%.o,$(SRC3))
LIB3   := $(wildcard c/*.a)

DIR4   := c/getgkmweights
SRC4   := $(wildcard c/getgkmweights/*.c)
OBJ4   := $(patsubst %.c,%.o,$(SRC4))
LIB4   := $(wildcard c/*.a)

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
	$(CC) $(OBJ1) $(LIB1) $(CFLAGS) -o gkmPWMlasso -lm
gkmPWM     : $(OBJ2) $(LIB2)
	$(CC) $(OBJ2) $(LIB2) $(CFLAGS) -o gkmPWM -lm
mapTF      : $(OBJ3) $(LIB3)
	$(CC) $(OBJ3) $(LIB3) $(CFLAGS) -o mapTF -lm
getgkmweights      : $(OBJ4) $(LIB4)
	$(CC) $(OBJ4) $(LIB4) $(CFLAGS) -o getgkmweights -lm
clean:
	$(RM) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) gkmPWMlasso gkmPWM mapTF getgkmweights
