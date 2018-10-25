# Toolchain, using mingw on windows
CC = $(OS:Windows_NT=x86_64-w64-mingw32-)gcc
PKG_CFG = $(OS:Windows_NT=x86_64-w64-mingw32-)pkg-config
RM = rm

# flags
CFLAGS = -O3 -march=native -Wall
OMPFLAGS = -fopenmp -fopenmp-simd
HWLOCFLAGS = #$(shell $(PKG_CFG) --cflags hwloc)
SHRFLAGS = -fPIC -shared

# libraries
LDLIBS = -lmpfr
HWLOCLIBS = #$(shell $(PKG_CFG) --libs hwloc)

# filenames
SRC = autocorrelation.c
EXT = $(if $(filter $(OS),Windows_NT),.exe,.out)
TARGET = $(SRC:.c=_2018$(EXT))
SHREXT = $(if $(filter $(OS),Windows_NT),.dll,.so)
SHRTGT = $(SRC:.c=$(SHREXT))


all: $(TARGET) $(SHRTGT) 

$(TARGET): $(SRC)
	$(CC) $(SRC) -o $(TARGET) $(CFLAGS) $(OMPFLAGS) $(LDLIBS) $(HWLOCFLAGS) $(HWLOCLIBS)

$(SHRTGT): $(SRC)
	$(CC) $(SRC) -o $(SHRTGT) $(SHRFLAGS) $(CFLAGS) $(OMPFLAGS) $(LDLIBS) $(HWLOCFLAGS) $(HWLOCLIBS)

force: 
	$(CC) $(SRC) -o $(TARGET) $(CFLAGS) $(OMPFLAGS) $(LDLIBS) $(HWLOCFLAGS) $(HWLOCLIBS)
	$(CC) $(SRC) -o $(SHRTGT) $(SHRFLAGS) $(CFLAGS) $(OMPFLAGS) $(LDLIBS) $(HWLOCFLAGS) $(HWLOCLIBS)

clean:
	@[ -f $(TARGET) ] && $(RM) $(TARGET) || true
	@[ -f $(SHRTGT) ] && $(RM) $(SHRTGT) || true

.PHONY: all clean force clean_target clean_shrtgt
