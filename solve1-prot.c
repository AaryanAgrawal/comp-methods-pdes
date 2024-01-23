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

/************************************************************************
*                                                                       *
*     Update: 22-01-2024                                                *
*                                                                       *
*     Aaryan Agrawal
*     Class of '25 
*     Dartmouth College 
*     Hanover, NH 03755
* 
*     22-01-2024.
*           changed array indexing to work for:
*                 r is stored in rows 0 thru neq-1;
*                 b is stored in rows 0 thru neq-1,                    
*                       columns 0 thru 2*ihalfb;                       
*
*                                                                      *
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#define min( A, B ) ( ((A) < (B)) ? (A) : (B) )

void solve1( kkk, b, r, neq, ihalfb, mdim )
int kkk, neq, ihalfb, mdim;
float b[], r[];
{
    int nrs, ihbp, k, kk, kc, ki, lim, j, jc;
    int nn, iband, i, ji, iback, jp, kr, mr;
    float c, sum, pivot;

    nrs = neq - 1;
    ihbp = ihalfb;
    switch ( kkk )
    {
        case 1:  // Triangularize matrix b
            for (k = 0; k < nrs; k++)
            {
                pivot = b[k * mdim + ihbp];
                for (i = k + 1; (i < neq) && (ihbp > 0); i++)
                {
                    c = -b[i * mdim + ihbp] / pivot;
                    if (c != 0.0)
                    {
                        b[i * mdim + ihbp] = c;
                        for (j = 1; j <= ihalfb; j++)
                        {
                            if (ihbp + j < mdim) {
                                b[i * mdim + ihbp + j] += c * b[k * mdim + ihbp + j];
                            }
                        }
                    }
                }
            }
            break;
        case 2:  // Modify Load Vector r
            for (i = 1; i < neq; i++)
            {
                sum = 0.0;
                for (j = 0; j < ihalfb; j++)
                {
                    if (i - j > 0) {
                        sum += b[i * mdim + j] * r[i - j - 1];
                    }
                }
                r[i] += sum;
            }

            // Back solution
            r[neq - 1] /= b[(neq - 1) * mdim + ihbp];
            for (iback = 1; iback < neq; iback++)
            {
                i = neq - iback - 1;
                sum = 0.0;
                for (j = ihbp + 1; j < mdim && i + j - ihbp < neq; j++)
                {
                    sum += b[i * mdim + j] * r[i + j - ihbp];
                }
                r[i] = (r[i] - sum) / b[i * mdim + ihbp];
            }
            break;
        case 3:  // Triangularize matrix b using the Doolittle Method
            for ( k = 0; k < nrs; k++ )
            {
                pivot = b[k*mdim+ihbp];
                kk = k;
                kc = ihbp;
                for ( i = kk, kc--; (i < neq) && (kc > 0); i++, kc-- )
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
                            if (jc < mdim) // To avoid accessing beyond the last column
                            {
                                b[i*mdim+j] += (c * b[k*mdim+jc]);
                            }
                        }
                    }
                }
            }
            // Modify Load Vector r
            nn = neq + 1;
            iband = 2 * ihalfb + 1;
            for ( i = 1; i < neq; i++ )
              {
                jc = ihbp - i + 1;
                ji = 1;
                if ( jc <= 0 )
                  {
                    jc = 1;
                    ji = i - ihbp + 1;
                  }
                sum = 0.0;
                for ( j = jc; j < ihalfb; j++ )
                  {
                    sum = sum + b[i*mdim+j] * r[ji];
                    ji++;
                  }
                //printf("r[%d] = %e\n", i, r[i]);
                r[i] = r[i] + sum;
                //printf("r[%d] = %e\n", i, r[i]);
              }
            // Back solution
            r[neq-1] /= b[(neq-1)*mdim+ihbp];
            //printf("r[%d] = %f\n", neq, r[neq]);
          for ( iback = 1; iback < neq; iback++ )
            {
              i = nn - iback;
              jp = i-1;
              kr = ihbp;
              mr = (int)min( iband, (ihalfb + iback) );
              sum = 0.0;
              for ( j = kr; j < mr; j++ )
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

    // Write r to file
    char buf[0x100];
    snprintf(buf, sizeof(buf), "%d.txt", ihalfb);
    FILE *fp;
    printf("writing to %s\n", buf);
    fp = fopen(buf, "w");
    for (int ii = 0; ii < neq; ii++) {
        fprintf(fp, "%f\n", r[ii]);
    }
    fclose(fp);
}