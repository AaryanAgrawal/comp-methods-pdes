void basisf(phi,dpx,dpy,dj,xi,eta,x,y)

/*     How to call this function from a C language program:

       basisf (phi, dpx, dpy, &dj, xi, eta, x, y);
          where passed parameters are defined as follows,
          float phi[4], dpx[4], dpy[4], dj, xi, eta, x[4], y[4];
*/

/*     Input passed parameters     */
       float xi, eta, x[4], y[4];

/*     Output passed parameters     */
       float phi[4], dpx[4], dpy[4], *dj;



{
 /* local variables */
 float dpxi[4], dpe[4], dja[4], djrecip;

 /* basis functions */
 phi[0] = (1.-xi)*(1.-eta)*.25;
 phi[1] = (1.+xi)*(1.-eta)*.25;
 phi[2] = (1.+xi)*(1.+eta)*.25;
 phi[3] = (1.-xi)*(1.+eta)*.25;

 /* derivatives with respect to xi */
 dpxi[0] = -(1.-eta)*.25;
 dpxi[1] = -dpxi[0];
 dpxi[2] =  (1.+eta)*.25;
 dpxi[3] = -dpxi[2];

 /* derivatives with respect to eta */
 dpe[0] = -(1.-xi)*.25;
 dpe[3] = -dpe[0];
 dpe[1] = -(1.+xi)*.25;
 dpe[2] = -dpe[1];

 /* jacobian */
 dja[0] = 0.; dja[1] = 0.; dja[2]=0.; dja[3] = 0.;

 dja[0] = x[0]*dpxi[0] + x[1]*dpxi[1] + x[2]*dpxi[2] + x[3]*dpxi[3];
 dja[1] = y[0]*dpxi[0] + y[1]*dpxi[1] + y[2]*dpxi[2] + y[3]*dpxi[3];
 dja[2] = x[0]*dpe[0]  + x[1]*dpe[1]  + x[2]*dpe[2]  + x[3]*dpe[3];
 dja[3] = y[0]*dpe[0]  + y[1]*dpe[1]  + y[2]*dpe[2]  + y[3]*dpe[3];

 *dj = dja[0]*dja[3] - dja[1]*dja[2];
 djrecip = 1./(*dj);

 /* derivatives with respect to x,y */
 dpx[0] = (dja[3]*dpxi[0] - dja[1]*dpe[0])  * djrecip;
 dpx[1] = (dja[3]*dpxi[1] - dja[1]*dpe[1])  * djrecip;
 dpx[2] = (dja[3]*dpxi[2] - dja[1]*dpe[2])  * djrecip;
 dpx[3] = (dja[3]*dpxi[3] - dja[1]*dpe[3])  * djrecip;
 dpy[0] = (-dja[2]*dpxi[0] + dja[0]*dpe[0]) * djrecip;
 dpy[1] = (-dja[2]*dpxi[1] + dja[0]*dpe[1]) * djrecip;
 dpy[2] = (-dja[2]*dpxi[2] + dja[0]*dpe[2]) * djrecip;
 dpy[3] = (-dja[2]*dpxi[3] + dja[0]*dpe[3]) * djrecip;

 return;
}