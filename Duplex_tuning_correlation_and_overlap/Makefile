CFLAGS="-O2"
CC="gcc"
MFLAG=-lm

all: tune_rho deg_seq

tune_rho: tune_rho.c rank_utils.c 
	$(CC) $(CFLAGS) -o tune_rho tune_rho.c rank_utils.c $(MFLAG) 


deg_seq: deg_seq.c  utils.c
	$(CC) $(CFLAGS) -o deg_seq deg_seq.c utils.c $(MFLAG) 

clean:
	rm tune_rho deg_seq
