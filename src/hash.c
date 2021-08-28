/**
  * sha256 by using  openssl lib
  */

#include "openssl/sha.h"
#include "hufu_uov_parameter.h"
#include "hash.h"


int pqc_sha256 (const unsigned char * m , unsigned long long mlen, unsigned char * digest) {
#if 32 == HUFU_UOV_HASH_SIZE	
	  SHA256_CTX sha256;  
	  SHA256_Init( &sha256 ); 
	  SHA256_Update( &sha256 , m , mlen );
	  SHA256_Final( digest , &sha256 );
#elif 48 == HUFU_UOV_HASH_SIZE
	  SHA512_CTX sha384;
	  SHA384_Init( &sha384 );
	  SHA384_Update( &sha384 , m , mlen );
	  SHA384_Final( digest , &sha384 );
#elif 64 == HUFU_UOV_HASH_SIZE
	  SHA512_CTX sha512;
	  SHA512_Init( &sha512 );
	  SHA512_Update( &sha512 , m , mlen );
	  SHA512_Final( digest , &sha512 );
#else 
	  error: un-supported HUFU_UOV_HASH_SIZE
#endif	
	  return 0;
  }
 
