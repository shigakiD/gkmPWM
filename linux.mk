CC      := gcc
CFLAGS  := -O3 -fPIC -fno-loop-optimize -fno-aggressive-loop-optimizations -Ic/include
LIB_DIR := c/lib
LDFLAGS := -L$(LIB_DIR) -l:libopenblas_haswellp-r0.3.20.a -lm

# Directories containing source files
DIRS    := c/gkmPWMlasso c/gkmPWM c/mapTF c/getgkmweights

# Object and executable files for each directory
OBJ1    := $(patsubst %.c,%.o,$(wildcard c/gkmPWMlasso/*.c))
OBJ2    := $(patsubst %.c,%.o,$(wildcard c/gkmPWM/*.c))
OBJ3    := $(patsubst %.c,%.o,$(wildcard c/mapTF/*.c))
OBJ4    := $(patsubst %.c,%.o,$(wildcard c/getgkmweights/*.c))
EXES    := $(notdir $(DIRS))

.PHONY: all clean $(EXES)

all: $(EXES)

# Define compilation pattern rule for all .c files
%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

# Specific targets for each executable
gkmPWMlasso: $(OBJ1)
	$(CC) $^ $(LDFLAGS) $(CFLAGS) -o $@

gkmPWM: $(OBJ2)
	$(CC) $^ $(LDFLAGS) $(CFLAGS) -o $@

mapTF: $(OBJ3)
	$(CC) $^ $(LDFLAGS) $(CFLAGS) -o $@

getgkmweights: $(OBJ4)
	$(CC) $^ $(LDFLAGS) $(CFLAGS) -o $@

clean:
	$(RM) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(EXES)
