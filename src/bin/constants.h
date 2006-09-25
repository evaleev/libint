
#ifndef _libint_constants_h_
#define _libint_constants_h_

static const char cart_comp[] = "XYZ";
static const char am_letter[] = "0pdfghiklmnoqrtuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
static const char *number[] = {"zero","one","two","three","four","five","six","seven","eight","nine","ten","eleven",
			       "twelve","thirteen","fourteen","fifteen","sixteen","seventeen","eighteen","nineteen","twenty",
                               "twentyone", "twentytwo", "twentythree", "twentyfour", "twentyfive", "twentysix", "twentyseven",
                               "twentyeight", "twentynine", "thirty"};

static inline int io(int i) { return i*(i+1)/2; }

/*----------------------------------------------------------------------------------
  hash(a,b) returns the index of the (a[0] a[1]) type product within a doublet.
  a contains x y and z exponents of functions on centers A and B, and b contains
  their angular momenta
 ----------------------------------------------------------------------------------*/
static inline int hash(int a[2][3], int b[2])
{
  int c[2] = {0,0};
  int i;

  if(b[0]){
    i=b[0]-a[0][0];
    c[0]=i+io(i)-a[0][1];
    }
  if(b[1]){
    i=b[1]-a[1][0];
    c[1]=i+io(i)-a[1][1];
    }

  return c[0]*io(b[1]+1)+c[1];
}

#endif
