CFLAGS = -g -Wall
LIBS = -lproj -lm

all: mktiles

MKTILES_OBJS = mktiles.o
mktiles: $(MKTILES_OBJS)
	$(CC) $(CFLAGS) -o mktiles $(MKTILES_OBJS) $(LIBS)
