
#define min( A, B ) ( ((A) < (B)) ? (A) : (B) )

void bansol(kkk, b, r, neq, iband, ndim, mdim)

/*
      How to call this function from a C language program:

       bansol (kkk, b, r, neq, iband, ndim, mdin);
          where passed parameters are declared as follows,
          float b[ndim][mdim], r[ndim];
          int   kkk, neq, iband, ndim, mdim;

       and where the meaning of these arguments is,
          kkk   = 1, perform lu decomposition destructively
                  2, perform back substitution
                  3, perform both options 1 and 2
          b     = left-hand side matrix, declared [ndim][mdim] in
                  invoking routine
          r     = right-hand side vector, declared [ndim] in
                  invoking routine
          neq   = number of equations
          iband = diagonal (column 1) + half-bandwidth
          ndim  = number of rows
          mdim  = number of columns

      Note: this function is a C language, literal translation of 
            the Fortran source routine bansol.f.  For ease in the
            translation activity, Fortran behaviors were incorporated,
            e.g. base 1 indexing through arrays in "for" loops, etc.

*/

/*     Input and output passed parameters     */
       float b[], r[];

/*     Input passed parameters     */
       int kkk, neq, iband, ndim, mdim;



{
 /* local variables */
 int nrs, nr, n, m, mr, i, j, k, l;
 float pivot, cp;

 nrs = neq-1;
 nr = neq;

                 /* if want (a.) only destructive lu decomposition or */
                 /*         (b.) both destructive lu decomposition and */
                 /*              back substitution */
 if ( kkk != 2 )
 {
  /* triangularize matrix b using doolittle method */
  for ( n=1; n<=nrs; n=n+1 )
  {
   m = n-1;
   mr = min(iband, (nr-m));
   pivot = 1./b[ (n-1)*mdim ]; /* compute reciprocal; divides are slow */
   for ( l=2; l<=mr; l=l+1 )
   {
    cp = b[ (n-1)*mdim + (l-1) ] * pivot; /* multiply by pivot is faster */
    i = m+l;                              /* than divide by pivot in the */
    j = 0;                                /* Fortran version             */
    for ( k=l; k<=mr; k=k+1 )
    {
     j = j+1;
     b[ (i-1)*mdim + (j-1) ] = b[ (i-1)*mdim + (j-1) ] -
        cp*b[ (n-1)*mdim + (k-1)];
    }
    b[ (n-1)*mdim + (l-1) ] = cp;
   }
  }
 }  /* end if ( kkk != 2 ) */


                 /* if want (a.) only back substitution or */
                 /*         (b.) both destructive lu decomposition and */
                 /*              back substitution */
 if ( kkk != 1 )
 {
  /* modify load vector r */
  for ( n=1; n<=nrs; n=n+1 )
  {
   m = n-1;
   mr = min(iband, (nr-m));
   cp = r[n-1];
   r[n-1] = cp / b[ (n-1)*mdim ];
   for ( l=2; l<=mr; l=l+1 )
   {
    i = m+l;
    r[i-1] = r[i-1] - b[ (n-1)*mdim + (l-1) ] * cp;
   }
  }
  /* back solution */
  r[nr-1] = r[nr-1] / b[ (nr-1)*mdim ];
  for ( i=1; i<=nrs; i=i+1 )
  {
   n = nr-i;
   m = n-1;
   mr = min(iband, (nr-m));
   for ( k=2; k<=mr; k=k+1 )
   {
    l = m+k;                      /* store computed displacements in load */
    r[n-1] = r[n-1] - b[ (n-1)*mdim + (k-1) ] * r[l-1];       /* vector r */
   }
  }
 }  /* end if ( kkk != 1 ) */

 return;
}