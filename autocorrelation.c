#ifdef unix
    #define __USE_LARGEFILE64
    #define _LARGEFILE_SOURCE
    #define _LARGEFILE64_SOURCE
#endif 

//#ifdef WINNT
#ifdef __MINGW64__
    // see number from: sdkddkver.h
    // https://docs.microsoft.com/fr-fr/windows/desktop/WinProg/using-the-windows-headers
    #define _WIN32_WINNT 0x0602 // Windows 8
    #include <Processtopologyapi.h>
    #include <processthreadsapi.h>
#endif

#include <stdio.h> 
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

//#include "gmp.h"
#include "mpfr.h"
#include <omp.h>

// Timing for benchmarking
#include <time.h>


#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define handle_error(msg) \
  do { perror(msg); exit(EXIT_FAILURE); } while (0)

      
      
// Compile binary:
//   gcc -O3 -march=native autocorrelation.c -o autocorrelation.out -Wall -lmpfr -fopenmp -fopenmp-simp
//
// Compile shared library:
//   gcc -O3 -march=native -fPIC -shared autocorrelation.c -o autocorrelation.so -Wall -lmpfr -fopenmp -fopenmp-simp
//
// Linux -> *unix* should be set, otherwise use -Dunix

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

void aCorrUpTo( uint8_t *buffer, uint64_t n, double *r, int k )
{
    // Accumulators
    uint64_t m = 0;
    uint64_t *rk = (uint64_t *) calloc(k, sizeof(uint64_t)); // Filled with 0s.
    uint64_t bk; // Temporary variable for corrections

    mpfr_t R,Rk,M,N,K,Bk,V; // High precision floats for final calculation
    mpfr_inits2(256, R, Rk, M, N, K, Bk, V, NULL); // 256bits floats, ~71 decimal digits
    
    
    // Accumulating...
    #pragma omp parallel // Mods made for omp can be suboptimal for single thread
    {
        #ifdef __MINGW64__
        int tid = omp_get_thread_num(); // Internal omp thread number
        HANDLE thandle = GetCurrentThread();
        _Bool result;
        
        GROUP_AFFINITY group = {0x0000000FFFFFFFFF, 0};
        group.Group = (tid<36)?0:1;
        result = SetThreadGroupAffinity(thandle, &group, NULL);
        if(!result) fprintf(stderr, "Failed setting output for tid=%i\n", tid);
        #endif
        
        #pragma omp for reduction(+:m), reduction(+:rk[:k])
        for (uint64_t i=0; i<n-k; i++)
        {
            m += (uint64_t)buffer[i];
            for (uint64_t j=0; j<k; j++)
            {
                rk[j] += (uint64_t)buffer[i]*(uint64_t)buffer[i+j];
            }
        }
        #pragma omp for reduction(+:m), reduction(+:rk[:k])
        for (uint64_t i=n-k; i<n; i++)
        {
            m += (uint64_t)buffer[i];
            for ( uint64_t j=0; j<n-i; j++)
            {
                rk[j] += (uint64_t)buffer[i]*(uint64_t)buffer[i+j];
            }
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
}

void aCorrUpToBit( uint8_t *buffer, uint64_t n, double *r, int k , int NthB)
{
    #define BITMASK(A) (A>>NthB & 1)
    #define AND(A,B) (A & B)
    //#define BITMASK_AND(A,B) AND(BITMASK(A), BITMASK(B))
    #define BITMASK_AND(A,B) BITMASK(AND(A,B)) // Should be slightly faster

    // Accumulators
    uint64_t m = 0;
    uint64_t *rk = (uint64_t *) calloc(k, sizeof(uint64_t)); // Filled with 0s.
    uint64_t bk; // Temporary variable for corrections

    mpfr_t R,Rk,M,N,K,Bk,V; // High precision floats for final calculation
    mpfr_inits2(256, R, Rk, M, N, K, Bk, V, NULL); // 256bits floats, ~71 decimal digits
    
    // Accumulating...
    #pragma omp parallel // Mods made for omp can be suboptimal for single thread
    {
        #ifdef __MINGW64__
        int tid = omp_get_thread_num(); // Internal omp thread number
        HANDLE thandle = GetCurrentThread();
        _Bool result;
        
        GROUP_AFFINITY group = {0x0000000FFFFFFFFF, 0};
        group.Group = (tid<36)?0:1;
        result = SetThreadGroupAffinity(thandle, &group, NULL);
        if(!result) fprintf(stderr, "Failed setting output for tid=%i\n", tid);
        #endif

        #pragma omp for reduction(+:m), reduction(+:rk[:k])
        for (uint64_t i=0; i<n-k; i++)
        {
            m += (uint64_t)BITMASK(buffer[i]);
            for (uint64_t j=0; j<k; j++)
            {
                rk[j] += (uint64_t)BITMASK_AND(buffer[i], buffer[i+j]);
            }
        }
        #pragma omp for reduction(+:m), reduction(+:rk[:k])
        for (uint64_t i=n-k; i<n; i++)
        {
            m += (uint64_t)BITMASK(buffer[i]);
            for ( uint64_t j=0; j<n-i; j++)
            {
                rk[j] += (uint64_t)BITMASK_AND(buffer[i], buffer[i+j]);
            }
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

void aCorrUpTo_double( uint8_t *buffer, uint64_t size, double *r, int k )
{
    double m = mean( buffer, size );
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
    // 1 - Reading the whole file into memory
    FILE *fp;
    int fd;
    uint64_t size;
 
    fp = fopen(argv[1], "rb");
    #ifdef __CYGWIN__
        struct stat sb;
        fd = fileno(fp);
        fstat(fd, &sb);
    #elif unix
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
    
    // The number of autocorrelations to compute
    int k = 2;
    if ( argc >= 3)
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
 
    
    ////////////////////
    // aCorrUpTo TEST //
    //////////////////// 
 
    // Result buffer, filled with 0 for compatibility with *_double* versions.   
    double *f = (double *) calloc(k, sizeof(double)); 
    // Timing stuff
    struct timespec start, end;
    double cpu_time_used;
    int timeflag = 0;
    
    if ( strcmp(argv[argc-1], "time\0") == 0 )
    {
        timeflag = 1;
    }
    
    if ( timeflag )
    {
        clock_gettime(CLOCK_MONOTONIC, &start);
    }
    // The computation itself!
    aCorrUpTo(buffer, size, f, k);
    // End of timing stuff
    if ( timeflag )
    {
        clock_gettime(CLOCK_MONOTONIC, &end);
        cpu_time_used = 1e-9*(end.tv_sec*1e9+end.tv_nsec - (start.tv_sec*1e9+start.tv_nsec));
    }
    // Printing an easy to paste into python format, for testing.
    printf("\nByte:\nf = array([");
    for (int i=0; i<k-1 ; i++)
    {
        printf("%0.15f, ", f[i]);
    }
    printf("%0.15f])\n",f[k-1]);
    if ( timeflag )
    {
        printf("\nComputation took %0.9f seconds!\n", cpu_time_used);
    }
    
    free(f);
    
     
    ///////////////////////
    // aCorrUpToBit TEST //
    ///////////////////////
/* 
    // Result buffer, filled with 0 for compatibility with *_double* versions.   
    double *g = (double *) calloc(k, sizeof(double)); 
    aCorrUpToBit(buffer, size, g, k, 0);
    // Printing an easy to paste into python format, for testing.
    printf("\nbit0:\nf = array([");
    for (int i=0; i<k-1 ; i++)
    {
      printf("%0.15f, ", g[i]);
    }
    printf("%0.15f])\n",g[k-1]);
 
    free(g);
 
    // Result buffer, filled with 0 for compatibility with *_double* versions.   
    double *h = (double *) calloc(k, sizeof(double)); 
    aCorrUpToBit(buffer, size, h, k, 7);
    // Printing an easy to paste into python format, for testing.
    printf("\nbit7:\nf = array([");
    for (int i=0; i<k-1 ; i++)
    {
      printf("%0.15f, ", h[i]);
    }
    printf("%0.15f])\n",h[k-1]);
 
    free(h);
    */ 
   
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
 
    // To test using "exact" double method. 
    // The parallel method might be more precise thanks to mpfr.
    // Both should match for the first several decimals.
    /*if ( argc >= 1000)
    {
        double g=0;
        g = aCorrk_double(buffer, size, k-1);
        printf("\n\n%0.15f",g);
    }*/
    return 0;
}
