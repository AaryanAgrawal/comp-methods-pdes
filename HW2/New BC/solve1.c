/*************************************************************************
*                                                                        *
*     function solve1                                                     *
*                                                                        *
*     Jay B. Perry                                                       *
*     Numerical Methods Lab                                              *
*     Thayer School of Engineering                                       *
*     Dartmouth College                                                  *
*     Hanover, NH 03755                                                  *
*     February, 1992 

*     2-22-99.  The best of the solve.c variants.  use this one!
*      (as modified 2/22/99, drl.  
*           changed the definition of mdim;
*           removed the unused argument ndim;
*           made the comments understandable.)
                                                   
*               
*     Solves b*x=r via LU decompositon
*                                                        
*     b = matrix.  It will be overwritten with its LU decomposition.                                   
*     r = vector.  It will be overwritten with x i.e. the answer.
*
*     mdim: 
*        declarations in calling program are 
*            float b[ndim][mdim], r[ndim] 
*        the parameter mdim is needed here.
*
*     neq        = number of equations                                  
*     ihalfb     = half bandwidth of b                                       
*     2*ihalfb+1 = bandwidth of b
*
*     r is stored in rows 1 thru neq;
*     b is stored in rows 1 thru neq,  
*                 columns 1 thru 2*ihalfb+1
*
*     kkk = 1 ==> Perform LU decomposition destructively                 *
*     kkk = 2 ==> Perform LU back substitution                           *
*     kkk = 3 ==> Perform both                                           *
*                                                                        *
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#define min( A, B ) ( ((A) < (B)) ? (A) : (B) )

void solve1( kkk, b, r, neq, ihalfb, mdim, text)
int kkk, neq, ihalfb, mdim;
float b[], r[];
  {
    int nrs, ihbp, k, kk, kc, ki, lim, j, jc;
    int nn, iband, i, ji, iback, jp, kr, mr;
    float c, sum, pivot;

    nrs = neq - 1;
    ihbp = ihalfb + 1;
    switch ( kkk )
      {
        case 1:  /***  Triangularize matrix b using the Doolittle    ***
                  ***  Method                                        ***/
          for ( k = 1; k <= nrs; k++ )
            {
              pivot = b[k*mdim+ihbp];
              kk = k + 1;
              kc = ihbp;
              for ( i = kk, kc--; (i <= neq) && (kc > 0); i++, kc-- )
                {
                  c = -b[i*mdim+kc] / pivot;
                  if ( c != 0.0 )
                    {
                      b[i*mdim+kc] = c;
                      ki = kc + 1;
                      lim = kc + ihalfb;
                      for ( j = ki; j <= lim; j++ )
                        {
                          jc = ihbp + j - kc;
                          b[i*mdim+j] = b[i*mdim+j] + c * b[k*mdim+jc];
                        }
                    }
                }
            }
          break;
        case 2:  /***  Modify Load Vector r ***/
          nn = neq + 1;
          iband = 2 * ihalfb + 1;
          for ( i = 2; i <= neq; i++ )
            {
              jc = ihbp - i + 1;
              ji = 1;
              if ( jc <= 0 )
                {
                  jc = 1;
                  ji = i - ihbp + 1;
                }
              sum = 0.0;
              for ( j = jc; j <= ihalfb; j++ )
                {
                  sum = sum + b[i*mdim+j] * r[ji];
                  ji++;
                }
              r[i] = r[i] + sum;
            }
        /***  Back solution  ***/
          r[neq] = r[neq] / b[neq*mdim+ihbp];
          for ( iback = 2; iback <= neq; iback++ )
            {
              i = nn - iback;
              jp = i;
              kr = ihbp + 1;
              mr = (int)min( iband, (ihalfb + iback) );
              sum = 0.0;
              for ( j = kr; j <= mr; j++ )
                {
                  jp++;
                  sum = sum + b[i*mdim+j] * r[jp];
                }
              r[i] = (r[i] - sum) / b[i*mdim+ihbp];
            }
          break;

        case 3:  /***  Triangularize matrix b using the Doolittle    ***
                  ***  Method                                        ***/
          for ( k = 1; k <= nrs; k++ )
            {
              pivot = b[k*mdim+ihbp];
              //printf("pivot = %f\t%d\n", pivot, k*mdim+ihbp);
              kk = k + 1;
              kc = ihbp;
              for ( i = kk, kc--; (i <= neq) && (kc > 0); i++, kc-- )
                {
                  c = -b[i*mdim+kc] / pivot;
                  if ( c != 0.0 )
                    {
                      b[i*mdim+kc] = c;
                      ki = kc + 1;
                      lim = kc + ihalfb;
                      for ( j = ki; j <= lim; j++ )
                        {
                          jc = ihbp + j - kc;
                          b[i*mdim+j] = b[i*mdim+j] + c * b[k*mdim+jc];
                        }
                    }
                }
            }
        /***  Modify Load Vector r ***/
          nn = neq + 1;
          iband = 2 * ihalfb + 1;
          for ( i = 2; i <= neq; i++ )
            {
              jc = ihbp - i + 1;
              ji = 1;
              if ( jc <= 0 )
                {
                  jc = 1;
                  ji = i - ihbp + 1;
                }
              sum = 0.0;
              for ( j = jc; j <= ihalfb; j++ )
                {
                  sum = sum + b[i*mdim+j] * r[ji];
                  ji++;
                }
              //printf("r[%d] = %e\n", i, r[i]);
              r[i] = r[i] + sum;
              //printf("r[%d] = %e\n", i, r[i]);
            }
        /***  Back solution  ***/
          r[neq] = r[neq] / b[(neq)*mdim+ihbp];
          //printf("r[%d] = %f\n", neq, r[neq]);
          for ( iback = 2; iback <= neq; iback++ )
            {
              i = nn - iback;
              jp = i;
              kr = ihbp + 1;
              mr = (int)min( iband, (ihalfb + iback) );
              sum = 0.0;
              for ( j = kr; j <= mr; j++ )
                {
                  jp++;
                  //printf("r[%d] = %f\n", jp, r[jp]);
                  sum = sum + b[i*mdim+j] * r[jp];
                  //printf("sum = %f\tb[%d][%d] = %f\tr[%d] = %f\n", sum, i, j, b[i*mdim+j], jp, r[jp]);
                }
              r[i] = (r[i] - sum) / b[i*mdim+ihbp];
              //printf("r[%d] = %f / %f\n", i, sum, b[i*mdim+ihbp]);
            }
          break;
      }
      //write r to file
      char buf[0x100];
      char stext[0x100];
      itoa(text, stext, 10);
      snprintf(buf, sizeof(buf), "%s.txt", stext);
      FILE *fp;
      printf("writing to %s\n", buf);
      fp = fopen(buf, "w");
      for(int ii=1; ii<=neq; ii++){
        fprintf(fp, "%f\n", r[ii]);
        //printf("%f\n", r[i]);
      }
      fclose(fp);
  }