/**
  *    hufu_uov.c
  *    Chengdong Tao
  **/
#include <stdio.h>
//#include <stdint.h>
#include <stdlib.h>
#include <string.h>
//#include <sys/time.h>
#include <math.h>

#include "hufu_uov.h"
#include "hufu_uov_parameter.h"
#include "circulant_toeplitz_matrix.h"
#include "hash.h"
#include "field256.h"
#include "rng.h"


static int generate_G1(circ_matrix_st *Q1, circ_matrix_st *G1) {
   memcpy(G1->first_row, Q1->first_row, HUFU_UOV_V);
   return 1;
}

static int generate_G2(circ_matrix_st *Q2, 
	                        circ_matrix_st *Q1, 
	                        circ_matrix_st *S, 
	                        circ_matrix_st *G2) {
    int status;
	circ_matrix_st cm;
	
	status = circulant_matrix_mul_vector(Q1, S, &cm);
	status = circulant_matrix_add(Q2, &cm, G2);

	return status;
}

static int generate_G3(circ_matrix_st *Q3, 
	                        circ_matrix_st *Q1, 
	                        circ_matrix_st *S, 
	                        circ_matrix_st *G3) {
    int status;
	circ_matrix_st cm;

	status = vector_mul_circulant_matrix(S, Q1, &cm);
	status = circulant_matrix_add(Q3, &cm, G3);

	return status;
}

static int generate_G4(hufu_poly_st       *hufuPoly,
	                        circ_matrix_st     *S,  
	                        toeplitz_matrix_st *G4) {
    int status;
	circ_matrix_st  am, bm, cm, dm, em, fm, ST;
	toeplitz_matrix_st tmpToep, teop;

	status = vector_mul_circulant_matrix(S, &(hufuPoly->q[0]), &am);
	status = circulant_matrix_mul_vector( S, &am,  &bm);
    status = circulant_matrix_mul_vector( S, &(hufuPoly->q[2]),  &cm);
	status = circulant_matrix_mul_vector(  &(hufuPoly->q[1]), S,    &dm);
	status = circulant_matrix_add(&bm, &cm, &em);
	status = circulant_matrix_add(&dm, &em, &fm);
	status = scalar_mul_toeplitz_matrix(&(hufuPoly->lambda), &(hufuPoly->matToeplitz), &tmpToep);
	status = circulant_to_toeplitz(&fm, &teop);
	status = toeplitz_matrix_add(&teop, &tmpToep, G4);

	return status;
}

/**
  *	 generate linear transformation L1
  */
static int generate_L1(circ_matrix_st *S) {
    int status;
	status = random_circulant_matrix(S);
    return status;
}

/**
  *	 generate central map F
  */
static int generate_F(hufu_poly_st *F) {
    int status, k;
	circ_matrix_st q3;
	toeplitz_matrix_st A;
	unsigned char scalar[HUFU_UOV_O];
	rand_bytes(scalar, HUFU_UOV_O);	
    status = random_toeplitz_matrix(&A);

	for (k = 0; k < HUFU_UOV_O; k++) {
		status = random_circulant_matrix(&(F[k].q[0]));
		status = random_circulant_matrix(&(F[k].q[1]));
		status = random_circulant_matrix(&(F[k].q[2]));
		memcpy(F[k].matToeplitz.gen_elements, A.gen_elements, 2*HUFU_UOV_O-1);
        F[k].lambda = scalar[k];
	}

	return status;
}

/**
  *	 generate linear transformation L2
  */
static int generate_L2(unsigned char *L2) {
    int status;
	status = rand_bytes(L2, HUFU_UOV_O*HUFU_UOV_O);
    return status;
}

static int compute_F_L1(hufu_poly_st *F, circ_matrix_st *S, hufu_poly_st *G) {
    int status, k;    
	
    for (k = 0; k < HUFU_UOV_O; k++) {
        status = generate_G1(&(F[k].q[0]), &(G[k].q[0]));
        status = generate_G2(&(F[k].q[1]), &(F[k].q[0]), S, &(G[k].q[1]));
		status = generate_G3(&(F[k].q[2]), &(F[k].q[0]), S, &(G[k].q[2]));
		status = generate_G4(&(F[k]), S, &(G[k].matToeplitz));
		G[k].lambda = 1;
	}

    return status;
}

int compute_L2_F_L1(unsigned char *L2, hufu_poly_st *G, hufu_poly_st *P) {
    int status, i, j, k;    
	
    for (k = 0; k < HUFU_UOV_O; k++) {
		circ_matrix_st  W[3];
		toeplitz_matrix_st W3;
		for (i = 0; i < 3; i++) {
		    zero_circulant_matrix(&(P[k].q[i]));
		}
		zero_toeplitz_matrix(&(P[k].matToeplitz));
		
		for (j = 0; j < HUFU_UOV_O; j++) {
            status = scalar_mul_circulant_matrix(&(L2[k*HUFU_UOV_O+j]), &(G[j].q[0]), &W[0]);
			status = circulant_matrix_add(&(P[k].q[0]), &W[0], &(P[k].q[0]));
            status = scalar_mul_circulant_matrix(&(L2[k*HUFU_UOV_O+j]), &(G[j].q[1]), &W[1]);
			status = circulant_matrix_add(&(P[k].q[1]), &W[1], &(P[k].q[1]));
            status = scalar_mul_circulant_matrix(&(L2[k*HUFU_UOV_O+j]), &(G[j].q[2]), &W[2]);
			status = circulant_matrix_add(&(P[k].q[2]), &W[2], &(P[k].q[2]));
            status = scalar_mul_toeplitz_matrix(&(L2[k*HUFU_UOV_O+j]), &(G[j].matToeplitz), &W3);
            status = toeplitz_matrix_add(&(P[k].matToeplitz), &W3, &(P[k].matToeplitz));
		}
		P[k].lambda = 1;
        
	}

    return status;
}

int hufu_uov_keypair(unsigned char *pk, unsigned char *sk) {
   int status, k;
   circ_matrix_st S;
   hufu_poly_st F[HUFU_UOV_O], G[HUFU_UOV_O], P[HUFU_UOV_O];
   unsigned char L2[HUFU_UOV_O*HUFU_UOV_O];

   /**
            *   generate linear transformation L1
            */
   status = generate_L1(&S);
   
   /**
            *  generate central map F
            */
   status = generate_F(F);
   
   /**
            *   generate linear transformation L2
            */   
    int success = 0;
   while(success == 0){
      status = generate_L2(L2);
      success = compute_inverse_matrix(L2, HUFU_UOV_O);
   }
	
   memcpy(sk, L2, HUFU_UOV_O*HUFU_UOV_O);
   memcpy(sk+HUFU_UOV_O*HUFU_UOV_O, F[0].matToeplitz.gen_elements, 2*HUFU_UOV_O-1);   

   for (k = 0; k < HUFU_UOV_O; k++) {
       memcpy(sk+HUFU_UOV_O*HUFU_UOV_O+2*HUFU_UOV_O-1+k, &(F[k].lambda), 1);   
   }
   
   for (k = 0; k < HUFU_UOV_O; k++) {
       memcpy(sk+HUFU_UOV_O*HUFU_UOV_O+3*HUFU_UOV_O-1+k*3*HUFU_UOV_V, F[k].q[0].first_row, HUFU_UOV_V);
	   memcpy(sk+HUFU_UOV_O*HUFU_UOV_O+3*HUFU_UOV_O-1+k*3*HUFU_UOV_V+HUFU_UOV_V, F[k].q[1].first_row, HUFU_UOV_V);
	   memcpy(sk+HUFU_UOV_O*HUFU_UOV_O+3*HUFU_UOV_O-1+k*3*HUFU_UOV_V+2*HUFU_UOV_V, F[k].q[2].first_row,HUFU_UOV_V);
   }  
   
   memcpy(sk+10*HUFU_UOV_O*HUFU_UOV_O+3*HUFU_UOV_O-1, S.first_row, HUFU_UOV_V);

   status = compute_F_L1(F, &S, G);
   status = compute_L2_F_L1(L2, G, P);

   for (k = 0; k < HUFU_UOV_O; k++) {
      memcpy(pk + k*(11*HUFU_UOV_O-1), P[k].q[0].first_row, HUFU_UOV_V);   
	  memcpy(pk + k*(11*HUFU_UOV_O-1) + HUFU_UOV_V, P[k].q[1].first_row, HUFU_UOV_V); 
	  memcpy(pk + k*(11*HUFU_UOV_O-1) + 2*HUFU_UOV_V, P[k].q[2].first_row, HUFU_UOV_V); 
	  memcpy(pk + k*(11*HUFU_UOV_O-1) + 3*HUFU_UOV_V, P[k].matToeplitz.gen_elements, 2*HUFU_UOV_O-1);   
   }

   return 1;
}


static int compute_constant_vector(hufu_poly_st       *f, 
	                                         const unsigned char *hash, 
	                                         circ_matrix_st      *cv, 
	                                         unsigned char       * res) { 
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

static int compute_coefficients_matrix(hufu_poly_st     *f, 
	                                               circ_matrix_st    *cv , 
	                                               unsigned char     *coefficientMatrix,
	                                               unsigned char     *y) {
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

    status = bytes_to_toeplitz_matrix(poly->matToeplitz.gen_elements, HUFU_UOV_O, HUFU_UOV_O, matToeplitz);
	for (i = 0; i < HUFU_UOV_O; i++) {
        for (j = 0; j < HUFU_UOV_O; j++) {
            Q[(3*HUFU_UOV_O+i)*(HUFU_UOV_V+HUFU_UOV_O)+3*HUFU_UOV_O+j] = mul(poly->lambda, matToeplitz[i*HUFU_UOV_O+j]);
		}
	}

	free(matToeplitz);
	
	return status;
}

static int sk_to_central_map(const unsigned char *sk, 
                                     hufu_poly_st        *F, 
                                     circ_matrix_st      *S, 
                                     unsigned char       *L2) {
	int status, k;
	circ_matrix_st a;
		
    memcpy(L2, sk, HUFU_UOV_O*HUFU_UOV_O);

	for (k = 0; k < HUFU_UOV_O; k++) {
	    memcpy(F[k].matToeplitz.gen_elements, sk+HUFU_UOV_O*HUFU_UOV_O, 2*HUFU_UOV_O-1);   
	}
		
	for (k = 0; k < HUFU_UOV_O; k++) {
		memcpy(&(F[k].lambda), sk+HUFU_UOV_O*HUFU_UOV_O+2*HUFU_UOV_O-1+k, 1);	
	}
		
	for (k = 0; k < HUFU_UOV_O; k++) {
		memcpy(F[k].q[0].first_row, sk+HUFU_UOV_O*HUFU_UOV_O+3*HUFU_UOV_O-1+k*3*HUFU_UOV_V,  HUFU_UOV_V);
		memcpy(F[k].q[1].first_row, sk+HUFU_UOV_O*HUFU_UOV_O+3*HUFU_UOV_O-1+k*3*HUFU_UOV_V+HUFU_UOV_V,  HUFU_UOV_V);
		memcpy(F[k].q[2].first_row, sk+HUFU_UOV_O*HUFU_UOV_O+3*HUFU_UOV_O-1+k*3*HUFU_UOV_V+2*HUFU_UOV_V, HUFU_UOV_V);
	}  
		
	memcpy(S->first_row, sk+10*HUFU_UOV_O*HUFU_UOV_O+3*HUFU_UOV_O-1,  HUFU_UOV_V);

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

	status = bytes_to_toeplitz_matrix(f[0].matToeplitz.gen_elements, HUFU_UOV_O, HUFU_UOV_O, toeplitzMat);
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

int hufu_uov_sign(const unsigned char *sk, unsigned char *docHash, unsigned char *sign) {
    int status;
	unsigned char  v[HUFU_UOV_O], salt[HUFU_UOV_SALT_SIZE], hashWithSalt[HUFU_UOV_O], L2Y[HUFU_UOV_O];
	unsigned char invfValue[HUFU_UOV_SIGNATURE_SIZE];
	unsigned char *coeffMat = (unsigned char*)malloc(HUFU_UOV_O*HUFU_UOV_O*sizeof(unsigned char));
	circ_matrix_st cv, s;
    int success = 0;
	hufu_poly_st f[HUFU_UOV_O];
	unsigned char L2[HUFU_UOV_O*HUFU_UOV_O];

    status = sk_to_central_map(sk, f, &s, L2);
	status = compute_inverse_matrix(L2, HUFU_UOV_O);

	while (success == 0) {
		while (success == 0) {
	       status = random_circulant_matrix(&cv);
	       status = compute_coefficients_matrix(f, &cv, coeffMat, v);
           success = compute_inverse_matrix(coeffMat, HUFU_UOV_O);
	    }
	
	    success = 0;
        rand_bytes(salt, HUFU_UOV_SALT_SIZE);
        status = hash_with_salt(docHash, salt, hashWithSalt);
		status = matrix_mul_vector(L2, HUFU_UOV_O, HUFU_UOV_O, hashWithSalt,L2Y);
	    success = inverse_central_map(f, L2Y, &cv, coeffMat, v, invfValue);
	}
	
    status = inverse_linear_transformation(&s, invfValue, sign);
    memcpy(sign+(HUFU_UOV_O+HUFU_UOV_V), salt, HUFU_UOV_SALT_SIZE);

	free(coeffMat);
   
	return status;
}

static int pk_to_poly(const unsigned char *pk, hufu_poly_st *P) {
    int status, k;
	for (k = 0; k < HUFU_UOV_O; k++) {
	   memcpy(P[k].q[0].first_row, pk + k*(11*HUFU_UOV_O-1),  HUFU_UOV_V);	
	   memcpy(P[k].q[1].first_row, pk + k*(11*HUFU_UOV_O-1) + HUFU_UOV_V, HUFU_UOV_V); 
	   memcpy(P[k].q[2].first_row, pk + k*(11*HUFU_UOV_O-1) + 2*HUFU_UOV_V, HUFU_UOV_V); 
	   memcpy(P[k].matToeplitz.gen_elements, pk + k*(11*HUFU_UOV_O-1) + 3*HUFU_UOV_V, 2*HUFU_UOV_O-1);	
	   P[k].lambda = 1;
	}

    return status;
}

int hufu_uov_verify(const unsigned char *pk, unsigned char *sign, unsigned char *docHash) {
    int status, i;
	unsigned char hashWithSalt[HUFU_UOV_O], ver[HUFU_UOV_O], salt[HUFU_UOV_SALT_SIZE];

    hufu_poly_st  p[HUFU_UOV_O];
	status= pk_to_poly(pk, p);
    memcpy(salt, sign+HUFU_UOV_V+HUFU_UOV_O, HUFU_UOV_SALT_SIZE);
    status = hash_with_salt(docHash, salt, hashWithSalt);
	status = compute_hufu_poly_value(p, sign, ver);
	
    for (i = 0; i < HUFU_UOV_O; i++) {
        if (ver[i] != hashWithSalt[i]) {
			printf("refuse signature \n");
            return 0;
		}
	}
    printf("***************accept success************ \n");
	return 1;
}

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
	int status = hufu_uov_keypair(pk, sk);
	return 0;
}

int crypto_sign(unsigned char *sm, 
	                 unsigned long long *smlen, 
	                  unsigned char *m, 
	                 unsigned long long mlen, 
	                 const unsigned char *sk) {
	unsigned char msghash[HUFU_UOV_HASH_SIZE];
	unsigned char sign[HUFU_UOV_SIGNATURE_SIZE];
	int status ;
	memcpy(sm , m , mlen);
	smlen[0] = mlen + HUFU_UOV_SIGNATURE_SIZE;
	status = pqc_sha256(m, mlen, msghash);
	status = hufu_uov_sign(sk, m, sign);
	memcpy(sm+mlen, sign, HUFU_UOV_SIGNATURE_SIZE);
	return 0;
}

int crypto_sign_open(unsigned char *m, 
	                        unsigned long long *mlen, 
	                        unsigned char *sm, 
	                        unsigned long long smlen,
	                        const unsigned char *pk)
{
	if( HUFU_UOV_SIGNATURE_SIZE > smlen ) return -1;
	memcpy(m , sm , smlen-HUFU_UOV_SIGNATURE_SIZE );
	mlen[0] = smlen-HUFU_UOV_SIGNATURE_SIZE;
	
	unsigned char msghash[HUFU_UOV_HASH_SIZE];
	unsigned char sign[HUFU_UOV_SIGNATURE_SIZE];
	memcpy(sign, sm + mlen[0], HUFU_UOV_SIGNATURE_SIZE);
	int status ;	
	status = pqc_sha256(m, mlen[0], msghash);
	status = hufu_uov_verify( pk , sign , m );
	return status;
}


