CC=gcc
LDFLAGS= -lcrypto

SOURCES= rng.c hufu_uov.c PQCgenKAT_sign.c
HEADERS= rng.h hufu_uov.h

PQCgenKAT_sign: $(HEADERS) $(SOURCES)
	$(CC) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: clean

clean:
	-rm PQCgenKAT_sign *.req *.rsp 