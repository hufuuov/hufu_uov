/**
  *  general matrix, circulant matrix, and toeplitz matrix computation
  * Author: Chengdong Tao
  */
  
#include <stdio.h>
//#include <stdint.h>
#include <stdlib.h>
#include <string.h>
//#include <sys/time.h>
#include <math.h>

#include "circulant_toeplitz_matrix.h"
#include "field256.h"
#include "rng.h"



void print_matrix(unsigned char *matrix, int row, int col) {
    int i, j;
	for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            printf("%d ",(int)matrix[i*col+j]);
		}
		printf("\n");
	}
} 

void print_vector(unsigned char *v, int len) {
    int i;
	printf("\n");
	for (i = 0; i < len; i++) {
		if (i % (3*HUFU_UOV_O) == 0)
			printf("\n");
        printf("%d ", (int)v[i]);
	}
	printf("\n");
}

void swap(unsigned char *a, unsigned char *b) {
   unsigned char t;
 
   t  = *b;
   *b = *a;
   *a = t;
}


//compute  A^-1
int compute_inverse_matrix(unsigned char *a, int n) { 
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

int clean_vector(unsigned char *vec, int veclen) {
    int i, status;
	
	for (i = 0; i < veclen; i++) {
		  vec[i] = (unsigned char)0;
	}
	
	status = 1;
	
	return status;
  }


int vector_mul_vector(unsigned char *a,unsigned char *b,  int len,  unsigned char *c) {
    int i;
	unsigned char sum = 0;
	for (i = 0; i < len; i++) {
        sum ^= mul(a[i], b[i]);
	}
	*c = sum;
	return 1;
}

int matrix_mul_vector(unsigned char *matrix, int row, int col, unsigned char *v, unsigned char *res) {
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

int vector_mul_matrix(unsigned char *matrix, int row, int col, unsigned char *v, unsigned char *res) {
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

void matrix_add(unsigned char *a, 
                             unsigned char *b, 
                             int m, 
                             int n,  
                             unsigned char *c) {
   int i, j;
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         c[i*m+j] = a[i*m+j] + b[i*m+j];
	  }
   }

}



//A is m*n matrix, B is n*d matrix, return C=A*B is m *d matrix
void matrix_mul(unsigned char *a, 
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
void matrix_transpose(unsigned char *a, int row, int col, unsigned char *b) {
	int i,j;
	
	for (i = 0; i < row; i++){
		for(j = 0; j < col; j++) {
			b[j*row+i] = a[i*col+j];
		}
	}

}

/**********************************************
  *                circulant matrix computation                                *
  **********************************************/
  
int random_circulant_matrix(circ_matrix_st *circMat) {
    rand_bytes(circMat->first_row, HUFU_UOV_V);
    return 1;
}

void zero_circulant_matrix(circ_matrix_st *circMat) {
    int i;
	
	for (i = 0; i < HUFU_UOV_V; i++) {
        circMat->first_row[i] = 0;
	}
}

void print_circulant_matrix(circ_matrix_st *circMat) {
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
int vector_shift_right(unsigned char *array, int len, int nShift) {
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

int circulant_matrix_add(circ_matrix_st *a, circ_matrix_st *b, circ_matrix_st *c) {
	int  i;
	for (i = 0; i < HUFU_UOV_V; i++) {
		c->first_row[i] = add(a->first_row[i], b->first_row[i]);
	}
	return 1;
}

int generate_classic_circulant_matrix(circ_matrix_st *a,  
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

void print_classic_circulant_matrix(circ_matrix_st *circM) {
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

int circulant_matrix_transpose(circ_matrix_st *a, circ_matrix_st *at) {
	int i;
	at->first_row[0] = a->first_row[0];
	for (i = 1; i < HUFU_UOV_V; i++) {
        at->first_row[i] = a->first_row[HUFU_UOV_V-i];
	}
	return 1;
}

int scalar_mul_circulant_matrix(unsigned char *lambda, circ_matrix_st *cm, circ_matrix_st *res) {
	int i;
	for (i = 0; i < HUFU_UOV_V; i++) {
        res->first_row[i] = mul(*lambda, cm->first_row[i]); 
	}
    return 1;
}


int circulant_matrix_mul_vector(circ_matrix_st *circMat, 
	                                         circ_matrix_st *vec, 
	                                         circ_matrix_st *res) {
    int i, status;
	unsigned char str[HUFU_UOV_V];
	
	memcpy(str, circMat->first_row, HUFU_UOV_V);
	
	for (i = 0; i < HUFU_UOV_V; i++) {
	    status = vector_mul_vector(str, vec->first_row,HUFU_UOV_V, &(res->first_row[i]));
		vector_shift_right(str, HUFU_UOV_V, 1);
	}
    
	return status;    
}

#if 1
int vector_mul_circulant_matrix(circ_matrix_st *vec, 
	                                         circ_matrix_st *circMat,
	                                         circ_matrix_st *res) {
    int status;

	circ_matrix_st at;
	status = circulant_matrix_transpose(circMat, &at);
    status = circulant_matrix_mul_vector(&at, vec, res);

	return status;
}
#endif

#if 0
 int vector_mul_circulant_matrix(circ_matrix_st *a, circ_matrix_st *b, circ_matrix_st *c) {
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

int circulant_matrix_mul(circ_matrix_st *a, circ_matrix_st *b, circ_matrix_st *c) {
	int status;
    status = vector_mul_circulant_matrix(a, b, c);
	return status;
}

#if 1
// for debug
int generate_linear_transformation(circ_matrix_st* circMat, unsigned char *linearTrans) {
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


/******************************************
  *                    Toeplitz matrix computation                    *
  ******************************************/

void zero_toeplitz_matrix(toeplitz_matrix_st *toepMat) {
    int i;
	
	for (i = 0; i < 2*HUFU_UOV_O-1; i++) {
        toepMat->gen_elements[i] = 0;
	}
}

int random_toeplitz_matrix(toeplitz_matrix_st *toepMat) {
    rand_bytes(toepMat->gen_elements, 2*HUFU_UOV_O-1);
    return 1;
}

void print_toeplitz_matrix(toeplitz_matrix_st *toepMat) {
	int i;
	
	printf("\n");
	for (i = 0; i < 2*HUFU_UOV_O-1; i++) {
        printf("%d ", toepMat->gen_elements[i]);
	}
	printf("\n");
}


int toeplitz_matrix_add(toeplitz_matrix_st *a, toeplitz_matrix_st *b, toeplitz_matrix_st *c) {
	int  i;
	for (i = 0; i < 2*HUFU_UOV_O-1; i++) {
		c->gen_elements[i] = add(a->gen_elements[i], b->gen_elements[i]);
	}
	return 1;
}

int scalar_mul_toeplitz_matrix(unsigned char *lambda, 
	                                       toeplitz_matrix_st *cm, 
	                                       toeplitz_matrix_st *res) {
	int i;
	for (i = 0; i < 2*HUFU_UOV_O-1; i++) {
        res->gen_elements[i] = mul(*lambda, cm->gen_elements[i]); 
	}
    return 1;
}


int bytes_to_toeplitz_matrix(unsigned char *buf, int row, int col, unsigned char *matToeplitz) {
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

int toeplitz_matrix_to_bytes(unsigned char *matToeplitz, int row, int col, unsigned char *buf) {
    int i;
	for (i = 0; i < col; i++) {
        buf[i] = matToeplitz[i];
	}

	for (i = 0; i < row-1; i++) {
        buf[col+i] = matToeplitz[(i+1)*col];
	}

	return 1;
}

int circulant_to_toeplitz(circ_matrix_st *cm, toeplitz_matrix_st *bufToep) {
	int i;
    memcpy(bufToep->gen_elements, cm->first_row, HUFU_UOV_O);
	for (i = 0; i < HUFU_UOV_O-1; i++) {
		bufToep->gen_elements[HUFU_UOV_O+i] = cm->first_row[HUFU_UOV_V-i-1];
   }
	return 1;
}


