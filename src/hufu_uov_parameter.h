/**
  * HUFU-UOV parameters
  * Chengdong Tao
  */

#ifndef HUFU_UOV_PARAMETERS_H
#define HUFU_UOV_PARAMETERS_H


#define HUFU_UOV_FIELD_SIZE 256
#define HUFU_UOV_O 48
#define HUFU_UOV_V 144

#define HUFU_UOV_HASH_SIZE 32
#define HUFU_UOV_SALT_SIZE 16

#define HUFU_UOV_PK_SIZE (11*HUFU_UOV_O*HUFU_UOV_O-HUFU_UOV_O)
#define HUFU_UOV_SK_SIZE (10*HUFU_UOV_O*HUFU_UOV_O+6*HUFU_UOV_O-1)
#define HUFU_UOV_SIGNATURE_SIZE (HUFU_UOV_O+HUFU_UOV_V+HUFU_UOV_SALT_SIZE)



#endif

