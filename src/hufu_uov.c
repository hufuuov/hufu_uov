#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#include "openssl/sha.h"
#include "rng.h"
#include "hufu_uov.h"

static unsigned char Logtable[HUFU_UOV_FIELD_SIZE] = 
{ 0, 0, 25, 1, 50, 2, 26, 198, 75, 199, 27, 104, 51, 238, 223, 3, 
100, 4, 224, 14, 52, 141, 129, 239, 76, 113, 8, 200, 248, 105, 28, 193, 
125, 194, 29, 181, 249, 185, 39, 106, 77, 228, 166, 114, 154, 201, 9, 120, 
101, 47, 138, 5, 33, 15, 225, 36, 18, 240, 130, 69, 53, 147, 218, 142, 
150, 143, 219, 189, 54, 208, 206, 148, 19, 92, 210, 241, 64, 70, 131, 56, 
102, 221, 253, 48, 191, 6, 139, 98, 179, 37, 226, 152, 34, 136, 145, 16, 
126, 110, 72, 195, 163, 182, 30, 66, 58, 107, 40, 84, 250, 133, 61, 186, 
43, 121, 10, 21, 155, 159, 94, 202, 78, 212, 172, 229, 243, 115, 167, 87, 
175, 88, 168, 80, 244, 234, 214, 116, 79, 174, 233, 213, 231, 230, 173, 232, 
44, 215, 117, 122, 235, 22, 11, 245, 89, 203, 95, 176, 156, 169, 81, 160, 
127, 12, 246, 111, 23, 196, 73, 236, 216, 67, 31, 45, 164, 118, 123, 183, 
204, 187, 62, 90, 251, 96, 177, 134, 59, 82, 161, 108, 170, 85, 41, 157, 
151, 178, 135, 144, 97, 190, 220, 252, 188, 149, 207, 205, 55, 63, 91, 209, 
83, 57, 132, 60, 65, 162, 109, 71, 20, 42, 158, 93, 86, 242, 211, 171, 
68, 17, 146, 217, 35, 32, 46, 137, 180, 124, 184, 38, 119, 153, 227, 165, 
103, 74, 237, 222, 197, 49, 254, 24, 13, 99, 140, 128, 192, 247, 112, 7};


static unsigned char Alogtable[HUFU_UOV_FIELD_SIZE] = 
{ 1, 3, 5, 15, 17, 51, 85, 255, 26, 46, 114, 150, 161, 248, 19, 53, 
95, 225, 56, 72, 216, 115, 149, 164, 247, 2, 6, 10, 30, 34, 102, 170, 
229, 52, 92, 228, 55, 89, 235, 38, 106, 190, 217, 112, 144, 171, 230, 49, 
83, 245, 4, 12, 20, 60, 68, 204, 79, 209, 104, 184, 211, 110, 178, 205, 
76, 212, 103, 169, 224, 59, 77, 215, 98, 166, 241, 8, 24, 40, 120, 136, 
131, 158, 185, 208, 107, 189, 220, 127, 129, 152, 179, 206, 73, 219, 118, 154, 
181, 196, 87, 249, 16, 48, 80, 240, 11, 29, 39, 105, 187, 214, 97, 163, 
254, 25, 43, 125, 135, 146, 173, 236, 47, 113, 147, 174, 233, 32, 96, 160, 
251, 22, 58, 78, 210, 109, 183, 194, 93, 231, 50, 86, 250, 21, 63, 65, 
195, 94, 226, 61, 71, 201, 64, 192, 91, 237, 44, 116, 156, 191, 218, 117, 
159, 186, 213, 100, 172, 239, 42, 126, 130, 157, 188, 223, 122, 142, 137, 128, 
155, 182, 193, 88, 232, 35, 101, 175, 234, 37, 111, 177, 200, 67, 197, 84, 
252, 31, 33, 99, 165, 244, 7, 9, 27, 45, 119, 153, 176, 203, 70, 202, 
69, 207, 74, 222, 121, 139, 134, 145, 168, 227, 62, 66, 198, 81, 243, 14, 
18, 54, 90, 238, 41, 123, 141, 140, 143, 138, 133, 148, 167, 242, 13, 23, 
57, 75, 221, 124, 132, 151, 162, 253, 28, 36, 108, 180, 199, 82, 246, 1};


static unsigned char add(unsigned char a, unsigned char b) {
	return a ^ b;
}

static unsigned char sub(unsigned char a, unsigned char b) {
	return a ^ b;
}

static unsigned char mul(unsigned char a, unsigned char b) {
/* multiply two elements of GF(2^m)*/
    if (a && b) 
        return Alogtable[(Logtable[a] + Logtable[b])%(HUFU_UOV_FIELD_SIZE-1)];
    else 
		return 0;
}

static signed char gf_div(unsigned char a, unsigned char b) {
	int j;
	if (b == 0) {
		printf( "Division by zero\n" );
		abort();
	}
	if (a == 0) return (0);

	if ((j = Logtable[a] - Logtable[b]) < 0) j += (HUFU_UOV_FIELD_SIZE-1);

	return (Alogtable[j]);
}

static unsigned char inv(unsigned char in) {
	/* 0 is self inverting */
	if(in == 0) 
		return 0;
	else
        return Alogtable[((HUFU_UOV_FIELD_SIZE-1) - Logtable[in])];
}


static void print_matrix(unsigned char *matrix, int row, int col) {
    int i, j;
	for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            printf("%d ",(int)matrix[i*col+j]);
		}
		printf("\n");
	}
} 

static void print_vector(unsigned char *v, int len) {
    int i;
	printf("\n");
	for (i = 0; i < len; i++) {
		//if (i % (3*HUFU_UOV_O) == 0)
		//	printf("\n");
        printf("%d ", (int)v[i]);
		
	}
	printf("\n");
}

static void swap(unsigned char *a, unsigned char *b) {
   unsigned char t;
 
   t  = *b;
   *b = *a;
   *a = t;
}


//compute  A^-1
static int compute_inverse_matrix(unsigned char *a, int n) { 
    int *is,*js,i,j,k,l,u,v;
    unsigned char d,p;
	
    is=(int*)malloc(n*sizeof(int));
    js=(int*)malloc(n*sizeof(int));
	
    for (k=0; k<=n-1; k++) { 
		d=0;
		
        for (i=k; i<=n-1; i++) {
            for (j=k; j<=n-1; j++) { 
			    l=i*n+j; 
			    p=(a[l]);
                if (p>d) { 
				    d=p; 
				    is[k]=i; 
				    js[k]=j;
			    }
            }
        }
		
        if (d==0) { 
			free(is); 
			free(js); 
			printf("err**not inv\n");
            return(0);
        }
		
        if (is[k]!=k) {
            for (j=0; j<=n-1; j++) { 
		      	u=k*n+j; 
			    v=is[k]*n+j;
                swap(&(a[u]), &(a[v]));
            }
        }
		
        if (js[k]!=k) {
            for (i=0; i<=n-1; i++) { 
		    	u=i*n+k; 
			    v=i*n+js[k];
                swap(&(a[u]), &(a[v]));
            }
        }
		
        l=k*n+k;
        a[l]=inv(a[l]);
        for (j=0; j<=n-1; j++) {
            if (j!=k) { 
		    	u=k*n+j; 
		    	a[u]=mul(a[u],a[l]);
		    }
        }
		
        for (i=0; i<=n-1; i++) {
            if (i!=k) {
                for (j=0; j<=n-1; j++) {
                    if (j!=k) { 
			  	        u=i*n+j;
                        a[u]=a[u]^(mul(a[i*n+k],a[k*n+j]));
                    }
                }
            }
        }
		
        for (i=0; i<=n-1; i++) {
            if (i!=k) { 
		  	    u=i*n+k; 
		        a[u]=mul(a[u],a[l]);
		    }
        }
    }
	
    for (k=n-1; k>=0; k--) { 
		if (js[k]!=k) {
            for (j=0; j<=n-1; j++) { 
		  	    u=k*n+j; 
			    v=js[k]*n+j;
                swap(&(a[u]), &(a[v]));
            }
		}
		
        if (is[k]!=k) {
            for (i=0; i<=n-1; i++) { 
		  	    u=i*n+k; 
			    v=i*n+is[k];
                swap(&(a[u]), &(a[v]));
            }
        }
    }
	
    free(is);
	free(js);
	
    return(1);
}

static int clean_vector(unsigned char *vec, int veclen) {
    int i, status;
	
	for (i = 0; i < veclen; i++) {
		  vec[i] = (unsigned char)0;
	}
	
	status = 1;
	
	return status;
  }


static int vector_mul_vector(unsigned char *a,unsigned char *b,  int len,  unsigned char *c) {
    int i;
	unsigned char sum = 0;
	for (i = 0; i < len; i++) {
        sum ^= mul(a[i], b[i]);
	}
	*c = sum;
	return 1;
}

static int matrix_mul_vector(unsigned char *matrix, int row, int col, unsigned char *v, unsigned char *res) {
    int i, j;
	unsigned char temp;
	for (i = 0; i < row; i++) {
		temp = 0;
        for (j = 0; j < col; j++) {
            temp ^= mul(matrix[i*col + j], v[j]);
		}
		res[i] = temp;
	}
	return 1;
}

static int vector_mul_matrix(unsigned char *matrix, int row, int col, unsigned char *v, unsigned char *res) {
    int i, j;
	unsigned char temp;
	for (i = 0; i < col; i++) {
		temp = 0;
        for (j = 0; j < row; j++) {
            temp ^= mul(matrix[j*col + i], v[j]);
		}
		res[i] = temp;
	}
	return 1;
}

static int pqc_sha256 (const unsigned char * m , unsigned long long mlen, unsigned char * digest) {
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

typedef struct  {
	unsigned char first_row[HUFU_UOV_V];
} circ_matrix_st;

typedef struct  {
    circ_matrix_st q[3];
	unsigned char matToeplitz[2*HUFU_UOV_O-1];
	unsigned char lambda;
} hufu_poly_st;

static int random_circulant_matrix(circ_matrix_st *circMat) {
    rand_bytes(circMat->first_row, HUFU_UOV_V);
    return 1;
}

static void zero_circulant_matrix(circ_matrix_st *circMat) {
    int i;
	
	for (i = 0; i < HUFU_UOV_V; i++) {
        circMat->first_row[i] = 0;
	}
}

static void print_circulant_matrix(circ_matrix_st *circMat) {
	int i;
	
	printf("\n");
	for (i = 0; i < HUFU_UOV_V; i++) {
        printf("%d ", circMat->first_row[i]);
	}
	printf("\n");
}

/**
 * Description: shift right. Complexity is O(o).             
 * @param  unsigned char *array:   pointer to the array whose length is len.
 * @param  int len: length of the array
 * @param  int nShift: the length of shifting.
 **/
static int vector_shift_right(unsigned char *array, int len, int nShift) {
    int i = 0;
    unsigned char temp;
	
    do {
        i = (i + nShift) % len;
        temp = array[i];
        array[i] = array[0];
        array[0] = temp;
    } while(i);
	
	return 1;
}


static int circulant_matrix_add(circ_matrix_st *a, circ_matrix_st *b, circ_matrix_st *c) {
	int  i;
	for (i = 0; i < HUFU_UOV_V; i++) {
		c->first_row[i] = add(a->first_row[i], b->first_row[i]);
	}
	return 1;
}

static int generate_classic_circulant_matrix(circ_matrix_st *a,  
	                                                         unsigned char *circulantMatrix) {
    int i, j, status;
	unsigned char *str = (unsigned char *)malloc(HUFU_UOV_V*sizeof(unsigned char));
    memcpy(str, a->first_row, HUFU_UOV_V);
	for (i = 0; i < HUFU_UOV_V; i++) {
        for (j = 0; j < HUFU_UOV_V; j++) {
            circulantMatrix[i*HUFU_UOV_V+j] = str[j];
		}
		status = vector_shift_right(str, HUFU_UOV_V, 1);
	}
	free(str);
	return status;
}

static void print_classic_circulant_matrix(circ_matrix_st *circM) {
    int i, j;
	unsigned char *str = (unsigned char *)malloc(HUFU_UOV_V*sizeof(unsigned char));
    memcpy(str, circM->first_row, HUFU_UOV_V);
	printf("\n");
	for (i = 0; i < HUFU_UOV_V; i++) {
        for (j = 0; j < HUFU_UOV_V; j++) {
            printf("%d ", str[j]);
		}
		int status = vector_shift_right(str, HUFU_UOV_V, 1);
		printf("\n");
	}
	printf("\n");
	free(str);

}

static int circulant_matrix_transpose(circ_matrix_st *a, circ_matrix_st *at) {
	int i;
	at->first_row[0] = a->first_row[0];
	for (i = 1; i < HUFU_UOV_V; i++) {
        at->first_row[i] = a->first_row[HUFU_UOV_V-i];
	}
	return 1;
}

static int scalar_mul_circulant_matrix(unsigned char *lambda, circ_matrix_st *cm, circ_matrix_st *res) {
	int i;
	for (i = 0; i < HUFU_UOV_V; i++) {
        res->first_row[i] = mul(*lambda, cm->first_row[i]); 
	}
    return 1;
}

#if 0
static int vector_mul_circulant_matrix(circ_matrix_st *a, circ_matrix_st *b, circ_matrix_st *c) {
    int i, j, status;
	unsigned char temp;
	
	unsigned char * matrix = (unsigned char *)malloc(HUFU_UOV_V*HUFU_UOV_V*sizeof(unsigned char));
	status = generate_classic_circulant_matrix(b, matrix);
	for (i = 0; i < HUFU_UOV_V; i++) {
        temp = 0;
		for (j = 0; j < HUFU_UOV_V; j++) {
            temp ^=mul(a->first_row[j], matrix[j*HUFU_UOV_V+i]);
		}
		c->first_row[i] = temp;
	}

	free(matrix);

	return status;
}
#endif





static int circulant_matrix_mul_vector(circ_matrix_st *a, circ_matrix_st *b, circ_matrix_st *c) {
    int i, status;
	unsigned char str[HUFU_UOV_V];
	
	memcpy(str, a->first_row, HUFU_UOV_V);
	
	for (i = 0; i < HUFU_UOV_V; i++) {
	    status = vector_mul_vector(str, b->first_row,HUFU_UOV_V, &(c->first_row[i]));
		vector_shift_right(str, HUFU_UOV_V, 1);
	}
    
	return status;    
}


static int vector_mul_circulant_matrix(circ_matrix_st *a, circ_matrix_st *b, circ_matrix_st *c) {
    int status;

	circ_matrix_st at;
	status = circulant_matrix_transpose(a, &at);
    status = circulant_matrix_mul_vector(&at, b, c);

	return status;
}

static int circulant_matrix_mul(circ_matrix_st *a, circ_matrix_st *b, circ_matrix_st *c) {
	int status;
    status = vector_mul_circulant_matrix(a, b, c);
	return status;
}



//A is m*n matrix, B is n*d matrix, return C=A*B is m *d matrix
static void matrix_mul(unsigned char *a, 
                             unsigned char *b, 
                             int m, 
                             int n,  
                             int d, 
                             unsigned char *c) { 
    int i, j, k, u;
	
    for (i = 0; i< m ; i++) {
        for (j = 0; j< d; j++) { 
		    u = i*d + j; 
		    c[u]=0;
            for (k = 0; k < n; k++) {
               c[u]=c[u]^(mul(a[i*n+k],b[k*d+j]));
            }
        }
    }
}

//compute transpose(A)
static void matrix_transpose(unsigned char *a, int row, int col, unsigned char *b) {
	int i,j;
	
	for (i = 0; i < row; i++){
		for(j = 0; j < col; j++) {
			b[j*row+i] = a[i*col+j];
		}
	}

}

static int bytes_to_toeplitz_matrix(unsigned char *buf, int row, int col, unsigned char *matToeplitz) {
    int i, j, index;
	for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
			if ((i - j) <  0) {
                index = (j-i) + row-1;
			} else {
                index = i - j ;
			}
            matToeplitz[j*row+i] = buf[index];
		}
	}

    matToeplitz[(row-1)*col] = buf[2*row-2];
	
	return 1;
}

static int toeplitz_matrix_to_bytes(unsigned char *matToeplitz, int row, int col, unsigned char *buf) {
    int i;
	for (i = 0; i < col; i++) {
        buf[i] = matToeplitz[i];
	}

	for (i = 0; i < row-1; i++) {
        buf[col+i] = matToeplitz[(i+1)*col];
	}

	return 1;
}

static int generate_p1(circ_matrix_st *p1) {
    int status;
	status = random_circulant_matrix(p1);
	return status;
}

static int generate_p2(circ_matrix_st *p2) {
    int status;
	status = random_circulant_matrix(p2);
	return status;
}

static int generate_p3(circ_matrix_st *p3) {
    int status;
	status = random_circulant_matrix(p3);
	return status;
}

static int generate_q1(circ_matrix_st *p1, circ_matrix_st *q1) {
	memcpy(q1->first_row, p1->first_row, HUFU_UOV_V);
	return 1;
}

static int generate_q2(circ_matrix_st *p2, circ_matrix_st *p1, circ_matrix_st *s, circ_matrix_st *q2) {
    int status;
	circ_matrix_st cm;
	
	status = circulant_matrix_mul_vector(p1, s, &cm);
	status = circulant_matrix_add(p2, &cm, q2);

	return status;
}

static int generate_q3(circ_matrix_st *p3, circ_matrix_st *p1, circ_matrix_st *s, circ_matrix_st *q3) {
    int status;
	circ_matrix_st cm;

	status = vector_mul_circulant_matrix(s, p1, &cm);
	status = circulant_matrix_add(p3, &cm, q3);

	return status;
}

static int generate_s(circ_matrix_st *s) {
    int status;
	
	status = random_circulant_matrix(s);

	return status;
}

static int generate_matrix_s(circ_matrix_st *s, unsigned char *S) {
    int i, j, status;
	unsigned char *str = (unsigned char *)malloc(HUFU_UOV_V*sizeof(unsigned char));
    memcpy(str, s->first_row, HUFU_UOV_V);
	for (i = 0; i < HUFU_UOV_O; i++) {
        for (j = 0; j < HUFU_UOV_V; j++) {
            S[i*HUFU_UOV_V+j] = str[j];
		}
		status = vector_shift_right(str, HUFU_UOV_V, 1);
	}
	free(str);
	return status;
}

#if 0
// for debug
static int generate_linear_transformation(circ_matrix_st* circMat, unsigned char *linearTrans) {
    int status, i, j;

    circ_matrix_st cm;
	memcpy(cm.first_row, circMat->first_row, HUFU_UOV_V);
    status = clean_vector(linearTrans,(HUFU_UOV_O+HUFU_UOV_V)*(HUFU_UOV_O+HUFU_UOV_V));
	for (i = 0; i < HUFU_UOV_O; i++ ) {
        for (j = 0; j < HUFU_UOV_V; j++) {
            linearTrans[j*(HUFU_UOV_O+HUFU_UOV_V)+(i+HUFU_UOV_V)] = cm.first_row[j];
		}
		status = vector_shift_right(cm.first_row, HUFU_UOV_V, 1);
	}

	for (i = 0; i < HUFU_UOV_O+HUFU_UOV_V; i++) {
        linearTrans[i*(HUFU_UOV_O+HUFU_UOV_V)+i] = (unsigned char)1;
	}

	return status;
}
#endif

static int generate_p4(circ_matrix_st *q1, circ_matrix_st *q2, 
	                          circ_matrix_st *q3, circ_matrix_st *q4,
	                          circ_matrix_st *s,  circ_matrix_st *p4) {
    int status;
	circ_matrix_st  am, bm, cm, dm, em, fm;
	status = vector_mul_circulant_matrix(s, q1, &am);
	status = circulant_matrix_mul_vector(s, &am,  &bm);
    status = circulant_matrix_mul_vector( s, q3, &cm);
	status = circulant_matrix_mul_vector( q2, s,   &dm);
	status = circulant_matrix_add(&bm, &cm, &em);
	status = circulant_matrix_add(&dm, &em, &fm);
	status = circulant_matrix_add(q4, &fm, p4);

	return status;
}

static int circulant_to_toeplitz(circ_matrix_st *cm, unsigned char *bufToep) {
	int i;
    memcpy(bufToep, cm->first_row, HUFU_UOV_O);
	for (i = 0; i < HUFU_UOV_O-1; i++) {
		bufToep[HUFU_UOV_O+i] = cm->first_row[HUFU_UOV_V-i-1];
   }
	return 1;
}

static int generate_public_bytes_from_seed(unsigned char *seedPub, unsigned char *buf) {
	int i, quotient, remainder;
    unsigned char in[HUFU_UOV_HASH_SIZE], out[HUFU_UOV_HASH_SIZE];
	memcpy(in, seedPub, HUFU_UOV_HASH_SIZE);

	quotient = (int)(HUFU_UOV_BYTES_FROM_PUBSEED/HUFU_UOV_HASH_SIZE);
	for (i = 0; i < quotient; i++) {
		pqc_sha256(in, HUFU_UOV_HASH_SIZE, out);
        memcpy(buf+i*HUFU_UOV_HASH_SIZE, out, HUFU_UOV_HASH_SIZE);
		memcpy(in, out, HUFU_UOV_HASH_SIZE);
	}
	pqc_sha256(in, HUFU_UOV_HASH_SIZE, out);
	remainder = HUFU_UOV_BYTES_FROM_PUBSEED % HUFU_UOV_HASH_SIZE;
	memcpy(buf+quotient*HUFU_UOV_HASH_SIZE, out, remainder);

	return 1;
}


static int generate_secret_bytes_from_seed(unsigned char *seedSec, unsigned char *buf) {
	int i,  quotient, remainder;
    unsigned char in[HUFU_UOV_HASH_SIZE], out[HUFU_UOV_HASH_SIZE];
	memcpy(in, seedSec, HUFU_UOV_HASH_SIZE);
	quotient = (int)(HUFU_UOV_BYTES_FROM_SECSEED/HUFU_UOV_HASH_SIZE);
	for (i = 0; i < quotient; i++) {
		pqc_sha256(in, HUFU_UOV_HASH_SIZE, out);
        memcpy(buf+i*HUFU_UOV_HASH_SIZE, out, HUFU_UOV_HASH_SIZE);
		memcpy(in, out, HUFU_UOV_HASH_SIZE);
	}
	pqc_sha256(in, HUFU_UOV_HASH_SIZE, out);
	remainder = HUFU_UOV_BYTES_FROM_SECSEED % HUFU_UOV_HASH_SIZE;
	memcpy(buf+quotient*HUFU_UOV_HASH_SIZE, out, remainder);

	return 1;
}

int hufu_uov_keypair(unsigned char *pk, unsigned char *sk) {
    int status,  k;
	circ_matrix_st s,  a, q4, p4;
	unsigned char scalar[HUFU_UOV_O];
	unsigned char pubseed[HUFU_UOV_HASH_SIZE];
	unsigned char secseed[HUFU_UOV_HASH_SIZE];
    unsigned char *pubBuf = (unsigned char *)malloc(HUFU_UOV_BYTES_FROM_PUBSEED*sizeof(unsigned char));
	unsigned char *secBuf = (unsigned char *)malloc(HUFU_UOV_BYTES_FROM_SECSEED*sizeof(unsigned char));
	
	rand_bytes(pubseed, HUFU_UOV_HASH_SIZE);	
	rand_bytes(secseed, HUFU_UOV_HASH_SIZE);

    memcpy(sk, pubseed, HUFU_UOV_HASH_SIZE);
	memcpy(sk+HUFU_UOV_HASH_SIZE, secseed, HUFU_UOV_HASH_SIZE);
	memcpy(pk, pubseed, HUFU_UOV_HASH_SIZE);

	status = generate_public_bytes_from_seed(pubseed, pubBuf);
	status = generate_secret_bytes_from_seed(secseed, secBuf);
	
	memcpy(s.first_row, secBuf, HUFU_UOV_V);
    memcpy(a.first_row, secBuf+HUFU_UOV_V, HUFU_UOV_V);
	memcpy(scalar, secBuf+2*HUFU_UOV_V, HUFU_UOV_O);
	//print_vector(pubBuf,HUFU_UOV_BYTES_FROM_PUBSEED);

    hufu_poly_st  p[HUFU_UOV_O];

	for (k = 0; k < HUFU_UOV_O; k++) {
		memcpy(p[k].q[0].first_row, pubBuf+k*3*HUFU_UOV_V, HUFU_UOV_V);
		memcpy(p[k].q[1].first_row, pubBuf+k*3*HUFU_UOV_V+HUFU_UOV_V, HUFU_UOV_V);
		memcpy(p[k].q[2].first_row, pubBuf+k*3*HUFU_UOV_V+2*HUFU_UOV_V, HUFU_UOV_V);
		//status = generate_q1(&(p[k].q[0]), &(f[k].q[0]));
        //status = generate_q2(&(p[k].q[1]), &(p[k].q[0]), &s, &(f[k].q[1]));
       // status = generate_q3(&(p[k].q[2]), &(p[k].q[0]), &s, &(f[k].q[2]));

        status = scalar_mul_circulant_matrix(&(scalar[k]), &a, &q4);
		status = generate_p4(&(p[k].q[0]),&(p[k].q[1]),&(p[k].q[2]), &q4, &s, &p4);

		//f[k].lambda = scalar[k];
		//p[k].lambda = 1;
       // status = circulant_to_toeplitz(&a, f[k].matToeplitz);
		status = circulant_to_toeplitz(&p4, p[k].matToeplitz);

		memcpy(pk+HUFU_UOV_HASH_SIZE+k*(2*HUFU_UOV_O-1), p[k].matToeplitz, 2*HUFU_UOV_O-1);
	}
    free(secBuf);
	free(pubBuf);

	return status;
}

static int compute_constant_vector(hufu_poly_st *f, 
	                                             const unsigned char *hash, 
	                                             circ_matrix_st *cv, 
	                                             unsigned char * res) {
    int k, status;
	circ_matrix_st cm;
	unsigned char temp = 0;
	
	for (k = 0; k < HUFU_UOV_O; k++) {
		status = circulant_matrix_mul_vector(&(f[k].q[0]), cv, &cm);
		status = vector_mul_vector(cm.first_row, cv->first_row, HUFU_UOV_V,  &temp);
		res[k] = add(temp, hash[k]);
	}

    return status;
}


static int compute_coefficients_matrix(hufu_poly_st *f, 
	                                                   circ_matrix_st *cv , 
	                                                   unsigned char *coefficientMatrix,
	                                                   unsigned char *y) {
	int i, k, status, counter = 0;
    circ_matrix_st  cm[2], sum;
	
    for (k = 0; k < HUFU_UOV_O; k++) {	
		status = circulant_matrix_mul_vector(&(f[k].q[1]), cv, &cm[0]);
        status = circulant_matrix_mul_vector(&(f[k].q[2]), cv, &cm[1]);
		status = circulant_matrix_add(&cm[0], &cm[1], &sum);
		coefficientMatrix[counter++] = f[k].lambda;
		for (i = 0; i < HUFU_UOV_O-1; i++) {
            coefficientMatrix[counter++] = sum.first_row[i];
		}	
		y[k] = sum.first_row[HUFU_UOV_O-1];
	}
	
    return status;
}


static int compute_univariate_poly_constant(unsigned char *alpha, 
	                                                            unsigned char *matA, 
	                                                            unsigned char *yFirst, 
	                                                            unsigned char *res) {
    int status;
	unsigned char vec[HUFU_UOV_O];
	unsigned char byteV;
	status = matrix_mul_vector(matA, HUFU_UOV_O, HUFU_UOV_O, alpha, vec);
	status = vector_mul_vector(vec, alpha, HUFU_UOV_O, &byteV);
	*res = *yFirst ^ byteV;

	return status;
}

static int compute_univariate_poly_linear(unsigned char *alpha, 
                                                          	unsigned char *beta,
	                                                            unsigned char *matA, 
	                                                            unsigned char *zLast, 
	                                                            unsigned char *res) {
    int status;
	unsigned char vec[HUFU_UOV_O];
	unsigned char byteV, tempB;
	status = matrix_mul_vector(matA, HUFU_UOV_O, HUFU_UOV_O, beta, vec);
	status = vector_mul_vector(vec, alpha, HUFU_UOV_O, &byteV);
	status = matrix_mul_vector(matA, HUFU_UOV_O, HUFU_UOV_O, alpha, vec);
	status = vector_mul_vector(vec, beta, HUFU_UOV_O, &tempB);
	*res = *zLast ^ byteV ^ tempB;

	return status;
}

static int compute_univariate_poly_quadratic(unsigned char *beta, unsigned char *matA, unsigned char *res){
    int status;
	unsigned char vec[HUFU_UOV_O];
	unsigned char byteV;
	status = matrix_mul_vector(matA, HUFU_UOV_O, HUFU_UOV_O, beta, vec);
	status = vector_mul_vector(vec, beta, HUFU_UOV_O, &byteV);
	*res = byteV;

	return status;
}

static int hufu_poly_to_matrix(hufu_poly_st *poly,
	                                      unsigned char *Q) {
    int status, i, j;
	circ_matrix_st cm1, cm2, cm3;

    unsigned char * matToeplitz = (unsigned char *)malloc(HUFU_UOV_O*HUFU_UOV_O*sizeof(unsigned char));

	//status = clean_vector(Q,(HUFU_UOV_O+HUFU_UOV_V)*(HUFU_UOV_O+HUFU_UOV_V));
	memcpy(cm1.first_row, poly->q[0].first_row, HUFU_UOV_V);
	memcpy(cm2.first_row, poly->q[1].first_row, HUFU_UOV_V);
	memcpy(cm3.first_row, poly->q[2].first_row, HUFU_UOV_V);
	
	for (i = 0; i < HUFU_UOV_V; i++) {
        for (j = 0; j < HUFU_UOV_V; j++) {
                Q[i*(HUFU_UOV_O+HUFU_UOV_V)+j] = cm1.first_row[j];
		}
		vector_shift_right(cm1.first_row, HUFU_UOV_V, 1);
	}

    for (i = 0; i < HUFU_UOV_O; i++ ) {
	    for (j = 0; j < HUFU_UOV_V; j++) {
		    Q[j*(HUFU_UOV_O+HUFU_UOV_V)+(i+HUFU_UOV_V)] = cm2.first_row[j];
	    }
	    status = vector_shift_right(cm2.first_row, HUFU_UOV_V, 1);
    }

    for (i = 0; i < HUFU_UOV_O; i++ ) {
	    for (j = 0; j < HUFU_UOV_V; j++) {
		        Q[(i+HUFU_UOV_V)*(HUFU_UOV_O+HUFU_UOV_V)+j] = cm3.first_row[j];
	    }
	    status = vector_shift_right(cm3.first_row, HUFU_UOV_V, 1);
    }

    status = bytes_to_toeplitz_matrix(poly->matToeplitz, HUFU_UOV_O, HUFU_UOV_O, matToeplitz);
	for (i = 0; i < HUFU_UOV_O; i++) {
        for (j = 0; j < HUFU_UOV_O; j++) {
            Q[(3*HUFU_UOV_O+i)*(HUFU_UOV_V+HUFU_UOV_O)+3*HUFU_UOV_O+j] = mul(poly->lambda, matToeplitz[i*HUFU_UOV_O+j]);
		}
	}

	free(matToeplitz);
	
	return status;
}

#if 0
//for debug
static void print_hufu_poly(hufu_poly_st *f) {
    print_circulant_matrix(&(f->q[0]));
    print_circulant_matrix(&(f->q[1]));
    print_circulant_matrix(&(f->q[2]));

	print_vector(f->matToeplitz, 2*HUFU_UOV_O-1);
	printf(" \n %d \n", f->lambda);
}
#endif

static int sk_to_central_map(const unsigned char *sk, hufu_poly_st *f, circ_matrix_st *s) {
	int status, k;
	circ_matrix_st a;
	unsigned char pubseed[HUFU_UOV_HASH_SIZE];
	unsigned char secseed[HUFU_UOV_HASH_SIZE];
	unsigned char *pubBuf = (unsigned char *)malloc(HUFU_UOV_BYTES_FROM_PUBSEED*sizeof(unsigned char));
	unsigned char *secBuf = (unsigned char *)malloc(HUFU_UOV_BYTES_FROM_SECSEED*sizeof(unsigned char));
			
	memcpy(pubseed, sk, HUFU_UOV_HASH_SIZE);
	memcpy(secseed, sk+HUFU_UOV_HASH_SIZE,  HUFU_UOV_HASH_SIZE);
		
	status = generate_public_bytes_from_seed(pubseed, pubBuf);
	status = generate_secret_bytes_from_seed(secseed, secBuf);

	memcpy(s->first_row, secBuf, HUFU_UOV_V);
	memcpy(a.first_row, secBuf+HUFU_UOV_V, HUFU_UOV_V);
		
	hufu_poly_st  p[HUFU_UOV_O];
		
	for (k = 0; k < HUFU_UOV_O; k++) {
		memcpy(p[k].q[0].first_row, pubBuf+k*3*HUFU_UOV_V, HUFU_UOV_V);
		memcpy(p[k].q[1].first_row, pubBuf+k*3*HUFU_UOV_V+HUFU_UOV_V, HUFU_UOV_V);
		memcpy(p[k].q[2].first_row, pubBuf+k*3*HUFU_UOV_V+2*HUFU_UOV_V, HUFU_UOV_V);
		status = generate_q1(&(p[k].q[0]), &(f[k].q[0]));
		status = generate_q2(&(p[k].q[1]), &(p[k].q[0]), s, &(f[k].q[1]));
		status = generate_q3(&(p[k].q[2]), &(p[k].q[0]), s, &(f[k].q[2]));
		f[k].lambda = secBuf[k+2*HUFU_UOV_V];
		status = circulant_to_toeplitz(&a, f[k].matToeplitz);
	}

    free(secBuf);
    free(pubBuf);

    return status;
}


static int hash_with_salt(unsigned char *msghash, unsigned char *salt, unsigned char *digest) {
	int status, i;
	unsigned int quotient, remainder;
	unsigned char gamma[HUFU_UOV_HASH_SIZE+HUFU_UOV_SALT_SIZE], in[HUFU_UOV_HASH_SIZE], out[HUFU_UOV_HASH_SIZE];
    memcpy(gamma, msghash, HUFU_UOV_HASH_SIZE);
	memcpy(gamma+HUFU_UOV_HASH_SIZE, salt, HUFU_UOV_SALT_SIZE);
    
	status = pqc_sha256(gamma,HUFU_UOV_HASH_SIZE+HUFU_UOV_SALT_SIZE,out);
	memcpy(in, out, HUFU_UOV_HASH_SIZE);
	quotient = (int)(HUFU_UOV_O/HUFU_UOV_HASH_SIZE);
	remainder = HUFU_UOV_O % HUFU_UOV_HASH_SIZE;
	for (i = 0; i < quotient; i++) {
		status = pqc_sha256(in, HUFU_UOV_HASH_SIZE, out);
        memcpy(digest+i*HUFU_UOV_HASH_SIZE, out, HUFU_UOV_HASH_SIZE);
		memcpy(in, out, HUFU_UOV_HASH_SIZE);
	}
	status = pqc_sha256(in, HUFU_UOV_HASH_SIZE, out);
	memcpy(digest+quotient*HUFU_UOV_HASH_SIZE, out, remainder);
	
	return status;
}

static int find_root_by_search(unsigned char *polyCoeffs, unsigned char *root) {
    int  i, success = 0;
	unsigned char x, x2, res;
	for (i = 0; i < HUFU_UOV_FIELD_SIZE; i++) {
		 x = (unsigned char)i;
		 x2 = mul(x, x);
         res = mul(polyCoeffs[2], x2) ^ mul(polyCoeffs[1], x) ^ polyCoeffs[0];
		 if (res == 0) { 
		 	*root = x;
	        if (x < 255) {
		        success = 1;
	        } else {
                success  = 0;
			}
			break;
		 }
	}
    return success;
}

static int inverse_central_map(hufu_poly_st *f, 
	                                      unsigned char *hashWithSalt, 
	                                      circ_matrix_st *vinegarVal, 
	                                      unsigned char *coeffMatFirstPart,
	                                      unsigned char *lastPart,
                                          unsigned char *invfValue) {
	int status, i, success = 0;
	unsigned char polyCoeffs[3], x = 0;
    unsigned char u[HUFU_UOV_O], y[HUFU_UOV_O], z[HUFU_UOV_O], alpha[HUFU_UOV_O], beta[HUFU_UOV_O];
	unsigned char *toeplitzMat = (unsigned char*)malloc(HUFU_UOV_O*HUFU_UOV_O*sizeof(unsigned char));

	status = compute_constant_vector(f, hashWithSalt, vinegarVal, u);
	status = matrix_mul_vector(coeffMatFirstPart, HUFU_UOV_O, HUFU_UOV_O, u, y);
	status = matrix_mul_vector(coeffMatFirstPart, HUFU_UOV_O, HUFU_UOV_O, lastPart, z);
    for (i = 0; i < HUFU_UOV_O - 1; i++) {
        alpha[i] = y[i+1];
		beta [i] = z[i+1];
	}
	alpha[HUFU_UOV_O - 1] = 0;
	beta[HUFU_UOV_O - 1] = 1;

	status = bytes_to_toeplitz_matrix(f[0].matToeplitz, HUFU_UOV_O, HUFU_UOV_O, toeplitzMat);
	status = compute_univariate_poly_constant(alpha, toeplitzMat, &(y[0]), &polyCoeffs[0]);
	status = compute_univariate_poly_linear(alpha, beta, toeplitzMat, &(z[0]),&polyCoeffs[1]);
	status = compute_univariate_poly_quadratic(beta, toeplitzMat, &polyCoeffs[2]);
	success = find_root_by_search(polyCoeffs, &x);
	
	for (i = 0; i < HUFU_UOV_V; i++) {
        invfValue[i] = vinegarVal->first_row[i];
	}
	for (i = 0; i < HUFU_UOV_O - 1; i++) {
        invfValue[HUFU_UOV_V+i] = y[i+1] ^ mul(z[i+1], x);
	}
	invfValue[HUFU_UOV_V+HUFU_UOV_O - 1] = x;
	
	free(toeplitzMat);

    return success;
}

static int inverse_linear_transformation(circ_matrix_st *s, unsigned char *y, unsigned char *sign) {
    int status, i;
	unsigned char  *vec = (unsigned char*)malloc(HUFU_UOV_V*sizeof(unsigned char));
	unsigned char *S = (unsigned char*)malloc(HUFU_UOV_O*HUFU_UOV_V*sizeof(unsigned char));
	
	status = generate_matrix_s(s,  S);
	status = vector_mul_matrix(S,HUFU_UOV_O,HUFU_UOV_V, y+HUFU_UOV_V, vec);
	for (i = 0; i < HUFU_UOV_V; i++) {
        sign[i] = y[i] ^ vec[i];
	}
	memcpy(sign+HUFU_UOV_V, y+HUFU_UOV_V, HUFU_UOV_O);

    free(vec);
	free(S);
	
	return status;
}


int hufu_uov_sign(const unsigned char *sk, unsigned char *docHash, unsigned char *sign) {
    int status;
	unsigned char  v[HUFU_UOV_O], salt[HUFU_UOV_SALT_SIZE],  hashWithSalt[HUFU_UOV_O];
	unsigned char invfValue[HUFU_UOV_SIGNATURE_SIZE];
	unsigned char *coeffMat = (unsigned char*)malloc(HUFU_UOV_O*HUFU_UOV_O*sizeof(unsigned char));
	circ_matrix_st cv, s;
    int success = 0;
	hufu_poly_st f[HUFU_UOV_O];

    status = sk_to_central_map(sk, f, &s);

	while (success == 0) {
	    status = random_circulant_matrix(&cv);
	    status = compute_coefficients_matrix(f, &cv, coeffMat, v);
        success = compute_inverse_matrix(coeffMat, HUFU_UOV_O);
	}
	
	success = 0;
	while (success == 0) {
        rand_bytes(salt, HUFU_UOV_SALT_SIZE);
        status = hash_with_salt(docHash, salt, hashWithSalt);
	    success = inverse_central_map(f, hashWithSalt, &cv, coeffMat, v, invfValue);
	}

    status = inverse_linear_transformation(&s, invfValue, sign);
    memcpy(sign+(HUFU_UOV_O+HUFU_UOV_V), salt, HUFU_UOV_SALT_SIZE);

	free(coeffMat);
   
	return status;
}


static int compute_hufu_poly_value(hufu_poly_st *f, unsigned char *x, unsigned char *v) {
    int status, k;
	unsigned char *Q = (unsigned char*)malloc((HUFU_UOV_O+HUFU_UOV_V)*(HUFU_UOV_O+HUFU_UOV_V)*sizeof(unsigned char));
	unsigned char *y = (unsigned char*)malloc((HUFU_UOV_O+HUFU_UOV_V) * sizeof(unsigned char));
	unsigned char c;
	for (k = 0; k < HUFU_UOV_O; k++) {
		status= hufu_poly_to_matrix(&f[k], Q);
        status = matrix_mul_vector(Q, (HUFU_UOV_O+HUFU_UOV_V), (HUFU_UOV_O+HUFU_UOV_V), x, y);
		status = vector_mul_vector(x, y, (HUFU_UOV_O+HUFU_UOV_V), &c);
		v[k] = c;
	}

	free(Q);
	free(y);

	return status;
}

static int pk_to_poly(const unsigned char *pk, hufu_poly_st *p) {
    int status, k;
	unsigned char pubseed[HUFU_UOV_HASH_SIZE];
    unsigned char *pubBuf = (unsigned char *)malloc(HUFU_UOV_BYTES_FROM_PUBSEED*sizeof(unsigned char));
	memcpy(pubseed, pk, HUFU_UOV_HASH_SIZE);
	status = generate_public_bytes_from_seed(pubseed, pubBuf);
	for (k = 0; k < HUFU_UOV_O; k++) {
		memcpy(p[k].q[0].first_row, pubBuf+k*3*HUFU_UOV_V, HUFU_UOV_V);
		memcpy(p[k].q[1].first_row, pubBuf+k*3*HUFU_UOV_V+HUFU_UOV_V, HUFU_UOV_V);
		memcpy(p[k].q[2].first_row, pubBuf+k*3*HUFU_UOV_V+2*HUFU_UOV_V, HUFU_UOV_V);

		p[k].lambda = 1;
		memcpy( p[k].matToeplitz, pk+HUFU_UOV_HASH_SIZE+k*(2*HUFU_UOV_O-1), 2*HUFU_UOV_O-1);
	}
    free(pubBuf);
    return status;
}


int hufu_uov_verify(const unsigned char *pk, unsigned char *sign, unsigned char *docHash) {
    int status, i;
	unsigned char hashWithSalt[HUFU_UOV_O], ver[HUFU_UOV_O];

    hufu_poly_st  p[HUFU_UOV_O];
	status= pk_to_poly(pk, p);
    status = hash_with_salt(docHash, sign+HUFU_UOV_V+HUFU_UOV_O,hashWithSalt);
	status = compute_hufu_poly_value(p, sign, ver);
    for (i = 0; i < HUFU_UOV_O; i++) {
        if (ver[i] != hashWithSalt[i]) {
            return 0;
		}
	}

	return 1;
}

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
	int status = hufu_uov_keypair(pk, sk);
	return 0;
}

int crypto_sign(unsigned char *sm, unsigned long long *smlen, const unsigned char *m, unsigned long long mlen, const unsigned char *sk) {
	unsigned char msghash[HUFU_UOV_HASH_SIZE];
	unsigned char sign[HUFU_UOV_SIGNATURE_SIZE];
	int status ;
	memcpy(sm , m , mlen);
	smlen[0] = mlen + HUFU_UOV_SIGNATURE_SIZE;
	status = pqc_sha256(m, mlen, msghash);
	status = hufu_uov_sign(sk, msghash, sign);
	memcpy(sm+mlen, sign, HUFU_UOV_SIGNATURE_SIZE);
	return 0;
}

int crypto_sign_open(unsigned char *m, unsigned long long *mlen, unsigned char *sm, unsigned long long smlen,const unsigned char *pk)
{
	if( HUFU_UOV_SIGNATURE_SIZE > smlen ) return -1;
	memcpy(m , sm , smlen-HUFU_UOV_SIGNATURE_SIZE );
	mlen[0] = smlen-HUFU_UOV_SIGNATURE_SIZE;
	
	unsigned char msghash[HUFU_UOV_HASH_SIZE];
	unsigned char sign[HUFU_UOV_SIGNATURE_SIZE];
	memcpy(sign, sm + mlen[0], HUFU_UOV_SIGNATURE_SIZE);
	int status ;	
	status = pqc_sha256(m, mlen[0], msghash);
	status = hufu_uov_verify( pk , sign , msghash );
	return status;
}


