#ifndef HUFU_UOV_H
#define HUFU_UOV_H


#define HUFU_UOV_FIELD_SIZE 256
#define HUFU_UOV_O 48
#define HUFU_UOV_V 144
#define HUFU_UOV_PUBSEED_SIZE 32
#define HUFU_UOV_SECSEED_SIZE 32
#define HUFU_UOV_PK_SIZE (2*HUFU_UOV_O*HUFU_UOV_O-HUFU_UOV_O+HUFU_UOV_PUBSEED_SIZE)
#define HUFU_UOV_SK_SIZE (HUFU_UOV_PUBSEED_SIZE+HUFU_UOV_SECSEED_SIZE)

#define HUFU_UOV_HASH_SIZE 32
#define HUFU_UOV_BYTES_FROM_PUBSEED (3*HUFU_UOV_O*HUFU_UOV_V)
#define HUFU_UOV_BYTES_FROM_SECSEED (2*HUFU_UOV_V+HUFU_UOV_O)
#define HUFU_UOV_SALT_SIZE 16
#define HUFU_UOV_SIGNATURE_SIZE (HUFU_UOV_O+HUFU_UOV_V+HUFU_UOV_SALT_SIZE)

int hufu_uov_keypair(unsigned char *pk, unsigned char *sk);
int hufu_uov_sign(const unsigned char *sk, unsigned char *docHash, unsigned char *sign); 
int hufu_uov_verify(const unsigned char *pk, unsigned char *sign, unsigned char *docHash);


int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);
int crypto_sign(unsigned char *sm, unsigned long long *smlen,
                          const unsigned char *m, unsigned long long mlen,
                          const unsigned char *sk);
int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                  unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk);


#endif
