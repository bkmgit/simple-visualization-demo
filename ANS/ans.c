/* An implementation of the ranged Asymmetric Numerical System
 * for an alphabet of upto 256 letters using a base of 512.
 *
 * This implementation follows Mart Simiskers [implementation](https://github.com/Martsim/crypto_seminar_2017_fall/)
 *
 * There is a dependency on [LibTomMath](http://www.libtom.net/LibTomMath/), 
 * a free open source portable number theoretic multiple-precision integer 
 * (MPI) library written entirely in C.
 *
 * Compile the code using
 *
 * gcc ans.c libtomath.a
 *
 * References
 *
 * [0] https://en.wikipedia.org/wiki/Asymmetric_numeral_systems
 *
 * [1] J. Duda "Asymmetric numeral systems: entropy coding combining speed of 
 * Huffman coding with compression rate of arithmetic coding"
 * https://arxiv.org/abs/1311.2540
 *
 * [2] J. Duda "Asymmetric Numeral Systems"
 * https://arxiv.org/abs/0902.0271
 *
 * [3] M. Simisker "A Review of Asymmetric Numeral Systems" 
 * https://courses.cs.ut.ee/MTAT.07.022/2017_fall/uploads/Main/mart-report-f17.pdf
 *
 * [4] R. Cheplyaka "Understanding Asymmetric Numeral Systems"
 * https://ro-che.info/articles/2017-08-20-understanding-ans
 *
 * [5] J. Gibbons "Coding with Asymmetric Numeral Systems (long version)"
 * http://www.cs.ox.ac.uk/jeremy.gibbons/publications/asymm-long.pdf
 *
 * [6] https://gist.github.com/reaandrew/21c0bbcc75086c72edfd
 *
 * [7] https://en.wikipedia.org/wiki/Bitwise_operations_in_C
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <inttypes.h>
#include "tommath.h"

#define arsize 20
#define maxnumsymbols 256
#define m 512 // need m >= numsymbols, power of 2 makes for easy division
int main()
{
        unsigned char numsin[arsize] =  {0,1,2,4,5,6,6,6,9,2,0,1,2,4,5,6,6,6,9,2};
        unsigned char numsout[arsize] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        unsigned char numtosymbol[maxnumsymbols];
        unsigned char symboltonum[maxnumsymbols];
        unsigned char symbolsin[arsize];
        unsigned char symbolsout[arsize];
        int symbolpresent[maxnumsymbols];
        int ls[maxnumsymbols];  // length of symbols range
        int s[m];
        int bs[maxnumsymbols]; // beginning of symbols range
        int i;
        int j;
        int numsymbols=0;
        unsigned char symbolcount=0;
        unsigned char count=0;
        unsigned long long x;
        char xstr[10005];
        mp_int xmp;
        mp_int tmpmp;
        mp_int lsmp;
        mp_int bsmp;
        mp_int mmp;
        mp_err result;

        if ((result = mp_init_size(&xmp, 10000)) != MP_OKAY) {
                printf("Error initializing the number.  %s",
                                mp_error_to_string(result));
                return 1;
        }
        if ((result = mp_init_size(&tmpmp, 10000)) != MP_OKAY) {
                printf("Error initializing the number.  %s",
                                mp_error_to_string(result));
                return 1;
        }
        if ((result = mp_init_multi(&lsmp,&bsmp,&mmp, NULL)) != MP_OKAY) {
                printf("Error initializing the numbers.  %s",
                                mp_error_to_string(result));
                return 1;
        }

        mp_set_u32(&mmp,(unsigned int) m);

        // Print out input array
        printf("Array to be encoded\n");
        printf("[");
        for(i=0;i<arsize-1;i++) printf(" %d,",numsin[i]);
        printf(" %d]\n",numsin[arsize-1]);

        // Count number of symbols
        for(i=0;i<maxnumsymbols;i++) symbolpresent[i]=0;
        for(i=0;i<arsize;i++) symbolpresent[numsin[i]]=1;
        for(i=0;i<maxnumsymbols;i++) numsymbols+=symbolpresent[i];
        // Initialize symbol tables
        for(i=0;i<maxnumsymbols;i++)
        {
                if(symbolpresent[i]==1) 
                {
                        numtosymbol[i]=symbolcount++;
                        symboltonum[count++]=i;
                }else{
                        numtosymbol[i]=255;
                }
        }
        /*printf("printing symbol table\n");
        for(i=0;i<maxnumsymbols;i++) printf("i %d symbol %hhu\n",i,numtosymbol[i]);
        for(j=0;j<count;j++) printf("j %d number %hhu\n",j,symboltonum[j]);*/
        // Get symbols from numbers in
        for(i=0;i<arsize;i++) symbolsin[i]=numtosymbol[numsin[i]];
        printf("Symbols to be encoded\n");
        printf("[");
        for(i=0;i<arsize-1;i++) printf(" %d,",symbolsin[i]);
        printf(" %d]\n",symbolsin[arsize-1]);

        // initialize arrays
        for(i=0;i<numsymbols;i++)
        {
                ls[i]=0;
                bs[i]=0;
        }
        for(i=0;i<m;i++) s[i]=0;
        // Setup encoding tables
        for(i=0;i<arsize;i++) ls[symbolsin[i]]++;
        bs[0]=0;
        for(i=1;i<numsymbols;i++) bs[i]=bs[i-1]+ls[i-1];
        for(j=0;j<m;j++)
        {
                unsigned long long thissum = 0;
                for(i=0;i<=j;i++) 
                {
                        thissum += ls[i];
                        if(thissum>j)
                        {
                                s[j]=i;
                                break;
                        }
                }
        }
        printf("Setup encoding tables\n");

        // Encode
        for(i=0;i<arsize;i++)
        {
                int st = symbolsin[i];
                //x = m*(x/ls[st]) + bs[st] + x%ls[st];
                mp_set_u32(&lsmp,(unsigned int) ls[st]);
                mp_set_u32(&bsmp,(unsigned int) bs[st]);
                result = mp_div(&xmp,&lsmp,&xmp,&tmpmp);
                result = mp_mul(&xmp,&mmp,&xmp);
                result = mp_add(&xmp,&bsmp,&xmp);
                result = mp_add(&xmp,&tmpmp,&xmp);
                //result = mp_to_radix(&xmp, xstr, 10000, NULL, 10);
                //printf("i %d x %s st %d ls %d bs %d\n",i,xstr,st,ls[st],bs[st]);
        }
        printf("The encoded number is %s\n",xstr);
        int compsize = mp_ubin_size(&xmp);
        // Decode
        for(i=arsize-1;i>=0;i--)
        {
                result = mp_mod(&xmp,&mmp,&tmpmp);
                unsigned long xmodm = mp_get_ul(&tmpmp);
                int st = s[xmodm];
                symbolsout[i]=st;
                //x=ls[st]*(x/m) + x%m - bs[st];
                mp_set_u32(&lsmp,(unsigned int) ls[st]);
                mp_set_i32(&bsmp,(int) bs[st]);
                result = mp_div(&xmp,&mmp,&xmp,NULL);
                result = mp_mul(&xmp,&lsmp,&xmp);
                result = mp_add(&xmp,&tmpmp,&xmp);
                result = mp_sub(&xmp,&bsmp,&xmp);
                //result = mp_to_radix(&xmp, xstr, 10000, NULL, 10);
                //printf("i %d x %s st %d ls %d bs %d\n",i,xstr,st,ls[st],bs[st]);
        }

        // Get numbers from symbols in
        for(i=0;i<arsize;i++) numsout[i]=symboltonum[symbolsout[i]];

        // Print out the array
        printf("Decoded array\n");
        printf("[");
        for(i=0;i<arsize-1;i++) printf(" %d,",numsout[i]);
        printf(" %d]\n",numsout[arsize-1]);

        printf("The original array has a size of %d bytes\n",8*arsize);
        printf("The compressed array has a size of %d bytes\n",compsize);
        // free memory
        mp_clear(&xmp);
        mp_clear(&tmpmp);
        mp_clear_multi(&lsmp,&bsmp,&mmp,NULL);
        return 0;
}
