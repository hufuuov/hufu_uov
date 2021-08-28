CC=gcc
LDFLAGS= -lcrypto

SOURCES= rng.c hufu_uov.c PQCgenKAT_sign.c circulant_toeplitz_matrix.c field256.c hash.c
HEADERS= rng.h hufu_uov.h hufu_uov_parameter.h  field256.h hash.h 

PQCgenKAT_sign: $(HEADERS) $(SOURCES)
	$(CC) -g -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: clean

clean:
	-rm PQCgenKAT_sign *.req *.rsp   PQCgenKAT_sign