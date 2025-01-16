
# Compilation target.
TARGET = massey

# Compiler.
CC = clang

# Compiler flags. The variable $PARI must be set to the PARI/GP installation directory.
CPPFLAGS   = -I. -I$(PARI)/GPDIR/include 
CFLAGS = -O3 -Wall -Wunused-function -fno-strict-aliasing -fomit-frame-pointer -pipe -flto=thin -march=native -pthread -g
LDFLAGS = -Wl,-O3 -pthread
LIBS = -L$(PARI)/GPDIR/lib -lpari -lm 
#STATIC = -static

# Compilation.
$(TARGET): $(TARGET).c
	$(CC) -o $@ $< $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)

clean:
	-$(RM) *.o $(ALL) massey

