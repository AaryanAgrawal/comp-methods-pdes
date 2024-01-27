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

//void jacobi( )
