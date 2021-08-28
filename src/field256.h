/**
  * field gf(2^8)
  * Chengdong Tao
  */

#ifndef FIELD_256_H
#define FIELD_256_H


unsigned char add(unsigned char a, unsigned char b);
unsigned char sub(unsigned char a, unsigned char b);
unsigned char mul(unsigned char a, unsigned char b);
unsigned char gf_div(unsigned char a, unsigned char b);
unsigned char inv(unsigned char in);

#endif

