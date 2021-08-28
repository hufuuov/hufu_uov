/**
  * circulant and toeplitz matrix in field gf(2^8)
  * Chengdong Tao
  */

#ifndef CIRCULANT_TOEPLITZ_MATRIX_H
#define CIRCULANT_TOEPLITZ_MATRIX_H

#include"hufu_uov_parameter.h"

typedef struct  {
	unsigned char first_row[HUFU_UOV_V];
} circ_matrix_st;

typedef struct  {
	unsigned char gen_elements[2*HUFU_UOV_O-1];
} toeplitz_matrix_st;

typedef struct  {
    circ_matrix_st q[3];
	toeplitz_matrix_st matToeplitz;
	unsigned char lambda;
} hufu_poly_st;

/**
  * general matrix in field gf(2^8) 
  */
void print_matrix(unsigned char *matrix, int row, int col);
void print_vector(unsigned char *v, int len);
void swap(unsigned char *a, unsigned char *b);
int compute_inverse_matrix(unsigned char *a, int n);
int clean_vector(unsigned char *vec, int veclen);
int vector_mul_vector(unsigned char *a,unsigned char *b,  int len,  unsigned char *c);
int matrix_mul_vector(unsigned char *matrix, int row, int col, unsigned char *v, unsigned char *res) ;
int vector_mul_matrix(unsigned char *matrix, int row, int col, unsigned char *v, unsigned char *res);
void matrix_mul(unsigned char *a, 
                             unsigned char *b, 
                             int m, 
                             int n,  
                             int d, 
                             unsigned char *c);
void matrix_add(unsigned char *a, 
                             unsigned char *b, 
                             int m, 
                             int n,  
                             unsigned char *c);

void matrix_transpose(unsigned char *a, int row, int col, unsigned char *b);



/**
  * general matrix in finite field gf(2^8) 
  */
int random_circulant_matrix(circ_matrix_st *circMat);
void zero_circulant_matrix(circ_matrix_st *circMat);
void print_circulant_matrix(circ_matrix_st *circMat);
/**
 * Description: shift right. Complexity is O(o).             
 * @param  unsigned char *array:   pointer to the array whose length is len.
 * @param  int len: length of the array
 * @param  int nShift: the length of shifting.
 **/
int vector_shift_right(unsigned char *array, int len, int nShift);
int circulant_matrix_add(circ_matrix_st *a, circ_matrix_st *b, circ_matrix_st *c);
int generate_classic_circulant_matrix(circ_matrix_st *a,  
	                                                 unsigned char *circulantMatrix);
void print_classic_circulant_matrix(circ_matrix_st *circM);
int circulant_matrix_transpose(circ_matrix_st *a, circ_matrix_st *at);
int scalar_mul_circulant_matrix(unsigned char *lambda, circ_matrix_st *cm, circ_matrix_st *res);
int circulant_matrix_mul_vector(circ_matrix_st *circMat, 
	                                         circ_matrix_st *vec, 
	                                         circ_matrix_st *res);
int vector_mul_circulant_matrix(circ_matrix_st *vec, 
	                                         circ_matrix_st *circMat,
	                                         circ_matrix_st *res);
int circulant_matrix_mul(circ_matrix_st *a, circ_matrix_st *b, circ_matrix_st *c);
int generate_linear_transformation(circ_matrix_st* circMat, unsigned char *linearTrans);



void zero_toeplitz_matrix(toeplitz_matrix_st *toepMat);
int random_toeplitz_matrix(toeplitz_matrix_st *toepMat);
void print_toeplitz_matrix(toeplitz_matrix_st *toepMat);
int bytes_to_toeplitz_matrix(unsigned char *buf, int row, int col, unsigned char *matToeplitz);
int toeplitz_matrix_to_bytes(unsigned char *matToeplitz, int row, int col, unsigned char *buf);
int circulant_to_toeplitz(circ_matrix_st *cm, toeplitz_matrix_st *bufToep);
int toeplitz_matrix_add(toeplitz_matrix_st *a, toeplitz_matrix_st *b, toeplitz_matrix_st *c) ;
int scalar_mul_toeplitz_matrix(unsigned char *lambda, 
	                                       toeplitz_matrix_st *cm, 
	                                       toeplitz_matrix_st *res);


#endif
