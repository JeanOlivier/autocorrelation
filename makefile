# Toolchain, using mingw on windows
CC = $(OS:Windows_NT=x86_64-w64-mingw32-)gcc
PKG_CFG = $(OS:Windows_NT=x86_64-w64-mingw32-)pkg-config

# flags
CFLAGS = -O3 -march=native -Wall
OMPFLAGS = -fopenmp -fopenmp-simd
HWLOCFLAGS = $(shell $(PKG_CFG) --cflags hwloc)

# libraries
LDLIBS = -lmpfr
HWLOCLIBS = $(shell $(PKG_CFG) --libs hwloc)

SRC = autocorrelation.c
EXT = $(if $(filter $(OS),Windows_NT),.exe,.out)
TARGET = $(SRC:.c=_2018$(EXT))


$(TARGET): $(SRC)
	$(CC) $(SRC) -o $(TARGET) $(CFLAGS) $(OMPFLAGS) $(HWLOCFLAGS) $(LDLIBS) $(HWLOCLIBS)

all: $(TARGET)

force: 
	$(CC) $(SRC) -o $(TARGET) $(CFLAGS) $(OMPFLAGS) $(HWLOCFLAGS) $(LDLIBS) $(HWLOCLIBS)


