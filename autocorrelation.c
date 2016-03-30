#ifdef UNIX
    #define __USE_LARGEFILE64
    #define _LARGEFILE_SOURCE
    #define _LARGEFILE64_SOURCE
#endif 

#include <stdio.h> 
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "gmp.h"
#include "mpfr.h"

#include <omp.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define handle_error(msg) \
  do { perror(msg); exit(EXIT_FAILURE); } while (0)

     
// Compile binary:
//   gcc -O3 autocorrelation.c -o autocorrelation.out -Wall -lmpfr -lgmp
//
// Compile shared library:
//   gcc -O3 -fPIC -shared autocorrelation.c -o autocorrelation.so -Wall -lmpfr -lgmp
//
// Linux -> add -DUNIX
// OpenMP paralellization -> add -fopenmp

double mean( uint8_t *buffer, uint64_t size)
{
    double m=0;
    for (uint64_t i=0; i<size; i++)
    {
        m += buffer[i];
    }
    m /= size;
    return m;
}

double var( uint8_t *buffer, uint64_t size)
{
    double m=0;
    double t;
    double mu = mean( buffer, size);
    for (uint64_t i=0; i<size; i++)
    {
        t = buffer[i];
        t -= mu;
        m += (t*t);
    }
    m /= size;
    return m;
}
    
void calc_corr(mpfr_t R , mpfr_t Rk, mpfr_t M, mpfr_t N, mpfr_t K, mpfr_t Bk, mpfr_rnd_t RND)
{  
    // Everything should have the same precision as R
    mpfr_prec_t p = mpfr_get_prec(R);

    mpfr_t R0; // Temp var
    mpfr_init2(R0, p);
    
    // r0 = m^2/size
    mpfr_pow_ui(R0, M, (uint64_t)2, RND);
    mpfr_div(R0, R0, N, RND);
    // r = rk - m^2/size
    mpfr_sub(R, Rk, R0, RND);
    // r += m*bk/size
    mpfr_div(R0, M, N, RND);
    mpfr_mul(R0, R0, Bk, RND);
    mpfr_add(R, R, R0, RND);
    // r -= m^2*k/size^2
    mpfr_pow_ui(R0, M, (uint64_t)2, RND);
    mpfr_mul(R0, R0, K, RND);
    mpfr_div(R0, R0, N, RND);
    mpfr_div(R0, R0, N, RND);
    mpfr_sub(R, R, R0, RND);
    // r /= (size-k)
    mpfr_sub(R0, N, K, RND);
    mpfr_div(R, R, R0, RND);
        
    mpfr_clear(R0);
}

void aCorrk_UnNormalized(mpfr_t R, uint8_t *buffer, uint64_t size, int k)
{
    // Accumulators
    uint64_t rk=0;
    uint64_t m=0;
    
    // Accumulating
    uint64_t i=0;
    for (; i<size-k; i++)
    {
        rk += (uint64_t)buffer[i]*(uint64_t)buffer[i+k];
        m += (uint64_t)buffer[i];
    }
    for (;i<size; i++)
    {
        m += (uint64_t)buffer[i]; // Don't forget the last k values
    }
    // Corrections
    uint64_t bk=0;
    for (uint64_t i=size-k ; i<size ; i++)
    {
        bk+=buffer[i];
    }
    for (uint64_t i=0 ; i<k ; i++)
    {
        bk+=buffer[i];
    }
    
    // Converting to "p"bits floats for precision on those big numbers 
    mpfr_prec_t p = mpfr_get_prec(R);
    mpfr_t Rk,M,N,K,Bk;//,R;
    mpfr_inits2(p, Rk, M, N, K, Bk, NULL);// R, NULL);   
   
    mpfr_set_uj(Rk, rk,   MPFR_RNDN);
    mpfr_set_uj(M,  m,    MPFR_RNDN);
    mpfr_set_uj(N,  size, MPFR_RNDN);
    mpfr_set_uj(K,  k,    MPFR_RNDN);
    mpfr_set_uj(Bk, bk,   MPFR_RNDN);
    
    calc_corr(R, Rk, M, N, K, Bk, MPFR_RNDN);
    
    mpfr_clears(Rk, M, N, K, Bk, NULL);
}

double aCorrk(uint8_t *buffer, uint64_t size, int k)
{
    mpfr_t R, V;
    mpfr_inits2(256, R, V, NULL);
    aCorrk_UnNormalized(R, buffer, size, k);
    aCorrk_UnNormalized(V, buffer, size, 0);
    mpfr_div(R, R, V, MPFR_RNDN);

    double r = mpfr_get_d(R, MPFR_RNDN);
    
    mpfr_clears(V, R, NULL);

    return r;
}

double *aCorrUpTo( uint8_t *buffer, uint64_t n, int k )
{
    // Accumulators
    uint64_t m = 0;
    uint64_t *rk = (uint64_t *) calloc(k, sizeof(uint64_t)); // Filled with 0s.
    
    double *r = (double *) malloc(k*sizeof(double)); // Return buffer
    uint64_t bk; // Temporary variable for corrections

    mpfr_t R,Rk,M,N,K,Bk,V; // High precision floats for final calculation
    mpfr_inits2(256, R, Rk, M, N, K, Bk, V, NULL); // 256bits floats, ~71 decimal digits
    
    // Accumulating...
    #pragma omp parallel // Mods made for omp can be suboptimal for single thread
    {
        // Local arrays avoid race conditions and segfault problems
        uint64_t *rk_local = (uint64_t *) calloc(k, sizeof(uint64_t));
        #pragma omp for reduction(+:m)
        for (uint64_t i=0; i<n-k; i++)
        {
            m += (uint64_t)buffer[i];
            for (uint64_t j=0; j<k; j++)
            {
                rk_local[j] += (uint64_t)buffer[i]*(uint64_t)buffer[i+j];
            }
        }
        #pragma omp for reduction(+:m)
        for (uint64_t i=n-k; i<n; i++)
        {
            m += (uint64_t)buffer[i];
            for ( uint64_t j=0; j<n-i; j++)
            {
                rk_local[j] += (uint64_t)buffer[i]*(uint64_t)buffer[i+j];
            }
        }
        #pragma omp critical // Manual reduction of local arrays to the main one
        for (uint64_t j=0; j<k; j++)
        {
            rk[j] += rk_local[j];
        }
    }

    // Computing correlations using MPFR 256bits precision floats ...
    mpfr_set_uj(M,  m,  MPFR_RNDN);
    mpfr_set_uj(N,  n,  MPFR_RNDN);
    for (uint64_t j=0; j<k; j++)
    {
        // Corrections
        bk = 0;
        for (uint64_t l=n-j; l<n; l++) 
        {
            bk += buffer[l];
        }
        for (uint64_t l=0; l<j; l++)
        {
            bk += buffer[l];
        }
        // Converting j-specific values to high-precision floats
        mpfr_set_uj(K,  j,     MPFR_RNDN);
        mpfr_set_uj(Rk, rk[j], MPFR_RNDN);
        mpfr_set_uj(Bk, bk,    MPFR_RNDN);
        // The real calculation, result stored in R.
        calc_corr(R, Rk, M, N, K, Bk, MPFR_RNDN);

        // Storing the variance in V for future usage
        if (j==0) 
        {
            mpfr_set(V,R, MPFR_RNDN);
        }
        // Dividing by the variance for normalisation
        mpfr_div(R, R, V, MPFR_RNDN);
        r[j] = mpfr_get_d(R, MPFR_RNDN);
    }
    
    mpfr_clears(R, Rk, M, N, K, Bk, V, NULL);
    free(rk);
    
    return r;
}

double *aCorrUpToBit( uint8_t *buffer, uint64_t n, int k , int NthB)
{
    #define BITMASK(A) (A>>NthB & 1)
    #define BITMASK_AND(A,B) (BITMASK(A) & BITMASK(B))

    // Accumulators
    uint64_t m = 0;
    uint64_t *rk = (uint64_t *) calloc(k, sizeof(uint64_t)); // Filled with 0s.
    
    double *r = (double *) malloc(k*sizeof(double)); // Return buffer
    uint64_t bk; // Temporary variable for corrections

    mpfr_t R,Rk,M,N,K,Bk,V; // High precision floats for final calculation
    mpfr_inits2(256, R, Rk, M, N, K, Bk, V, NULL); // 256bits floats, ~71 decimal digits
    
    // Accumulating...
    #pragma omp parallel // Mods made for omp can be suboptimal for single thread
    {
        // Local arrays avoid race conditions and segfault problems
        uint64_t *rk_local = (uint64_t *) calloc(k, sizeof(uint64_t));
        #pragma omp for reduction(+:m)
        for (uint64_t i=0; i<n-k; i++)
        {
            m += (uint64_t)BITMASK(buffer[i]);
            for (uint64_t j=0; j<k; j++)
            {
                rk_local[j] += (uint64_t)BITMASK_AND(buffer[i], buffer[i+j]);
            }
        }
        #pragma omp for reduction(+:m)
        for (uint64_t i=n-k; i<n; i++)
        {
            m += (uint64_t)BITMASK(buffer[i]);
            for ( uint64_t j=0; j<n-i; j++)
            {
                rk_local[j] += (uint64_t)BITMASK_AND(buffer[i], buffer[i+j]);
            }
        }
        #pragma omp critical // Manual reduction of local arrays to the main one
        for (uint64_t j=0; j<k; j++)
        {
            rk[j] += rk_local[j];
        }
    }

    // Computing correlations using MPFR 256bits precision floats ...
    mpfr_set_uj(M,  m,  MPFR_RNDN);
    mpfr_set_uj(N,  n,  MPFR_RNDN);
    for (uint64_t j=0; j<k; j++)
    {
        // Corrections
        bk = 0;
        for (uint64_t l=n-j; l<n; l++) 
        {
            bk += BITMASK(buffer[l]);
        }
        for (uint64_t l=0; l<j; l++)
        {
            bk += BITMASK(buffer[l]);
        }
        // Converting j-specific values to high-precision floats
        mpfr_set_uj(K,  j,     MPFR_RNDN);
        mpfr_set_uj(Rk, rk[j], MPFR_RNDN);
        mpfr_set_uj(Bk, bk,    MPFR_RNDN);
        // The real calculation, result stored in R.
        calc_corr(R, Rk, M, N, K, Bk, MPFR_RNDN);

        // Storing the variance in V for future usage
        if (j==0) 
        {
            mpfr_set(V,R, MPFR_RNDN);
        }
        // Dividing by the variance for normalisation
        mpfr_div(R, R, V, MPFR_RNDN);
        r[j] = mpfr_get_d(R, MPFR_RNDN);
    }
    
    mpfr_clears(R, Rk, M, N, K, Bk, V, NULL);
    free(rk);
    
    return r;
}

double aCorrk_double( uint8_t *buffer, uint64_t size, int k)
{
    double m=0;
    double t;
    double mu = mean( buffer, size );
    double s2 = var( buffer, size );
    uint64_t i=0;
    while (i<size-k)
    {
        t = buffer[i]-mu;
        t *= buffer[i+k]-mu;
        m += t;
        i++;
    }
    m /= ((size-k)*s2);
    return m;
}

double *aCorrUpTo_double( uint8_t *buffer, uint64_t size, int k )
{
    double m = mean( buffer, size );
    double *r = (double *) calloc(k, sizeof(double)); // Filled with 0s.
    double r0;

    uint64_t i=0;
    for (i=0; i<size-k; i++)
    {
        for (int j=0; j<k; j++)
        {
            r[j] += (buffer[i]-m)*(buffer[i+j]-m);
        }
    }
    for (; i<size; i++)
    {
        for ( uint64_t j=0; j<size-i; j++)
        {
            r[j] += (buffer[i]-m)*(buffer[i+j]-m);
        }
    }
    
    r0 = r[0]/size;
    for (uint64_t i=size; i>size-k; i--)
    {
        r[size-i] /= (i*r0);   // divide by (N-k)*variance
    }
    
    return r;
}

int readBig(uint8_t *buffer, uint64_t size, FILE *fp)
{
    uint64_t n;
    uint64_t r;
    uint64_t g = ((uint64_t) 1) << 30; // 2^30 = 1024^3 = 1 Giga
    n = size/g;
    r = size%g;
    for (uint64_t i=0; i<n; i++)
    {
        if ( g != fread(buffer+i*g, 1, g, fp) )
        {
            fprintf(stderr, "There was an error reading the input file (n).\n");
            return 0;
        }
    }
    if ( r != fread(buffer+n*g, 1, r, fp) )
    {
        fprintf(stderr, "There was an error reading the input file (r).\n");
        return 0;
    }
    return 1;
}

int main(int argc, char *argv[])
{
   FILE *fp;
   int fd;
   uint64_t size;

   fp = fopen(argv[1], "rb");
   #ifdef __CYGWIN__
       struct stat sb;
       fd = fileno(fp);
       fstat(fd, &sb);
   #elif UNIX
       struct stat64 sb;
       fd = fileno(fp);
       fstat64(fd, &sb);
   #else
       struct __stat64 sb;
       fd = _fileno(fp);
       _fstat64(fd, &sb);
   #endif
   size = (uint64_t)sb.st_size;

   uint8_t *buffer = (uint8_t *) malloc(sizeof(uint8_t)*size);
   if ( !readBig(buffer, size, fp) )
   {
       exit(1);
   }
   fclose(fp);

   #ifdef DEBUG
       for(uint64_t i = 0; i < 10; i++)
       {
         printf("[%"PRIu64"]=%X ", i, buffer[i]);
       }
       printf("[%"PRIu64"]=%X", size-1, buffer[size-1]);
       printf("\n");
   #endif
   
   int k = 1;
   if ( argc == 3)
   {
       k=atoi(argv[2]);
   }

   /*
   double g;
   for (int i=0; i<100; i++)
   {
       g = mean(buffer, size);
   }   
   printf("\nMean: %0.20f\n",g);
   */

   
   double *f;
   f = aCorrUpTo(buffer, size, k);
   // Printing an easy to paste into python format, for testing.
   printf("\nMPFR:\nf = array([");
   for (int i=0; i<k-1 ; i++)
   {
     printf("%0.15f, ", f[i]);
   }
   printf("%0.15f])\n",f[k-1]);
   
   /*
   f = aCorrUpToBit(buffer, size, k, 7);
   // Printing an easy to paste into python format, for testing.
   printf("\n7th bit:\nf = array([");
   for (int i=0; i<k-1 ; i++)
   {
     printf("%0.15f, ", f[i]);
   }
   printf("%0.15f])\n",f[k-1]);
   */

   
   
   /*f = aCorrUpTo_double(buffer, size, k);
   // Printing an easy to paste into python format, for testing.
   printf("\nDouble:\nf = array([");
   for (int i=0; i<k-1 ; i++)
   {
     printf("%0.15f, ", f[i]);
   }
   printf("%0.15f])\n",f[k-1]);
   
      
   free(f);*/

   return 0;
}
