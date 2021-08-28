#ifndef HUFU_UOV_H
#define HUFU_UOV_H

int hufu_uov_keypair(unsigned char *pk, unsigned char *sk);
int hufu_uov_sign(const unsigned char *sk, unsigned char *docHash, unsigned char *sign); 
int hufu_uov_verify(const unsigned char *pk, unsigned char *sign, unsigned char *docHash);

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);
int crypto_sign(unsigned char *sm, unsigned long long *smlen,
                           unsigned char *m, unsigned long long mlen,
                          const unsigned char *sk);
int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                  unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk);


#endif
