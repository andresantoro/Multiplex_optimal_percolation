CFLAGS="-O2"
CC="gcc"
MFLAG=-lm

all: tune_overlap

tune_overlap: tune_overlap.c
	$(CC) $(CFLAGS) -o tune_overlap tune_overlap.c $(MFLAG) 

clean:
	rm tune_overlap