#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <iostream>

 /*   Subroutine to generate a B-spline open knot vector with multiplicity
    equal to the order at the ends.
	
    c            = order of the basis function
    n            = the number of defining polygon vertices
    nplus2       = index of x() for the first occurence of the maximum knot vector value
    nplusc       = maximum value of the knot vector -- $n + c$
    x()          = array containing the knot vector
*/


//knot(n,c,x)

//int n,c;
//int x[];

void knot(int n,int c,int x[]) {
	int nplusc,nplus2,i;

	nplusc = n + c;
	nplus2 = n + 2;

	x[1] = 0;
		for (i = 2; i <= nplusc; i++){
		    if ( (i > c) && (i < nplus2) )
				x[i] = x[i-1] + 1;
	    else
				x[i] = x[i-1];
		}
}


 /*
	c        = order of the B-spline basis function
    d        = first term of the basis function recursion relation
    e        = second term of the basis function recursion relation
	h[]	     = array containing the homogeneous weights
    npts     = number of defining polygon vertices
    nplusc   = constant -- npts + c -- maximum number of knot values
    r[]      = array containing the rationalbasis functions
               r[1] contains the basis function associated with B1 etc.
    t        = parameter value
    temp[]   = temporary array
    x[]      = knot vector
*/	



//rbasis(c,t,npts,x,h,r)

//int c,npts;
//float t;
//int x[];
//float h[];
//float r[];

void rbasis(int c,float t,int npts,int x[],float h[],float r[]) {
	int nplusc;
	int i,j,k;
	float d,e;
	float sum;
	float temp[36];

	nplusc = npts + c;

/*		printf("knot vector is \n");
		for (i = 1; i <= nplusc; i++){
			printf(" %d %d \n", i,x[i]);
		}
		printf("t is %f \n", t);
*/

/* calculate the first order nonrational basis functions n[i]	*/

	for (i = 1; i<= nplusc-1; i++){
    	if (( t >= x[i]) && (t < x[i+1]))
			temp[i] = 1;
	    else
			temp[i] = 0;
	}

/* calculate the higher order nonrational basis functions */

	for (k = 2; k <= c; k++){
    	for (i = 1; i <= nplusc-k; i++){
        	if (temp[i] != 0)    /* if the lower order basis function is zero skip the calculation */
           		d = ((t-x[i])*temp[i])/(x[i+k-1]-x[i]);
	        else
				d = 0;

    	    if (temp[i+1] != 0)     /* if the lower order basis function is zero skip the calculation */
        		e = ((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1]);
	        else
    			e = 0;

    	    temp[i] = d + e;
		}
	}

	if (t == (float)x[nplusc]){		/*    pick up last point	*/
 		temp[npts] = 1;
	}
/*
	printf("Nonrational basis functions are \n");
	for (i=1; i<= npts; i++){
		printf("%f ", temp[i]);
	}
	printf("\n");
*/
/* calculate sum for denominator of rational basis functions */

	sum = 0.;
	for (i = 1; i <= npts; i++){
		    sum = sum + temp[i]*h[i];
	}

/* form rational basis functions and put in r vector */

	for (i = 1; i <= npts; i++){
    	if (sum != 0){
        	r[i] = (temp[i]*h[i])/(sum);}
		else
			r[i] = 0;
	}
}

/*   b[]         = array containing the defining polygon vertices
                  b[1] contains the x-component of the vertex
                  b[2] contains the y-component of the vertex
                  b[3] contains the z-component of the vertex
	h[]			= array containing the homogeneous weighting factors 
    k           = order of the B-spline basis function
    nbasis      = array containing the basis functions for a single value of t
    nplusc      = number of knot values
    npts        = number of defining polygon vertices
    p[,]        = array containing the curve points
                  p[1] contains the x-component of the point
                  p[2] contains the y-component of the point
                  p[3] contains the z-component of the point
    p1          = number of points to be calculated on the curve
    t           = parameter value 0 <= t <= npts - k + 1
    x[]         = array containing the knot vector
*/

//rbspline(npts,k,p1,b,h,p)

//int npts,k,p1;
//
//float b[];
//float h[];
//float p[];

void rbspline(int npts,int k,int p1,float b[],float h[],float p[]) {
	int i,j,icount,jcount;
	int i1;
	int x[30];		/* allows for 20 data points with basis function of order 5 */
	int nplusc;

	float step;
	float t;
	float nbasis[20];
	float temp;


	nplusc = npts + k;

/*  zero and redimension the knot vector and the basis array */

	for(i = 0; i <= npts; i++){
		 nbasis[i] = 0.;
	}

	for(i = 0; i <= nplusc; i++){
		 x[i] = 0.;
		}

/* generate the uniform open knot vector */

	knot(npts,k,x);

/*
	printf("The knot vector is ");
	for (i = 1; i <= nplusc; i++){
		printf(" %d ", x[i]);
	}
	printf("\n");
*/

	icount = 0;

/*    calculate the points on the rational B-spline curve */

	t = 0;
	step = ((float)x[nplusc])/((float)(p1-1));

	for (i1 = 1; i1<= p1; i1++){

		if ((float)x[nplusc] - t < 5e-6){
			t = (float)x[nplusc];
		}

	    rbasis(k,t,npts,x,h,nbasis);      /* generate the basis function for this value of t */
/*
		printf("t = %f \n",t);
		printf("nbasis = ");
		for (i = 1; i <= npts; i++){
			printf("%f  ",nbasis[i]);
		}
		printf("\n");
*/
		for (j = 1; j <= 3; j++){      /* generate a point on the curve */
			jcount = j;
			p[icount+j] = 0.;

			for (i = 1; i <= npts; i++){ /* Do local matrix multiplication */
				temp = nbasis[i]*b[jcount];
			    p[icount + j] = p[icount + j] + temp;
/*
				printf("jcount,nbasis,b,nbasis*b,p = %d %f %f %f %f\n",jcount,nbasis[i],b[jcount],temp,p[icount+j]);
*/
				jcount = jcount + 3;
			}
		}
/*
		printf("icount, p %d %f %f %f \n",icount,p[icount+1],p[icount+2],p[icount+3]);
*/
    	icount = icount + 3;
		t = t + step;
	}
}


int main(int argc , char** argv){

	int i;
	int npts,k,p1;

	float b[31];  /* allows for up to 10  control vertices */
	float h[11];  /* allows for up to 10  control vertices */
	float p[103]; /* allows for up to 100 points on curve */

	npts = 5;
	k = 3;     /* third order, change for other orders */
	p1 = 11;   /* eleven points on curve */

	for (i = 1; i <= 3*npts; i++){
		b[i] = 0.;
	}

/*	set all homogeneous weighting factros to 1.0 */

		for (i=1; i <= npts; i++){
		h[i] = 1.0;
	}

/*  vary the homogeneous weighting factor 0, 0.25, 1.0, 5.0 */

	h[3] = 1;

	for (i = 1; i <= 3*p1; i++){
		p[i] = 0.;
	}

/*
	Define the control polygon, Ex. 4.1, p. 140 in the z=1 plane because
    this is three dimensional routine. x=b[1], y=b[2], z=b[3], etc.
*/	
	b[1]=0;
	b[2]=0;
	b[3]=1;
	b[4]=1;
	b[5]=2;
	b[6]=1;
	b[7]=2.5;
	b[8]=0;
	b[9]=1;
	b[10]=4;
	b[11]=2;
	b[12]=1;
	b[13]=5;
	b[14]=0;
	b[15]=1;
	
	rbspline(npts,k,p1,b,h,p);

	printf("\nPolygon points\n");

	for (i = 1; i <= 3*npts; i=i+3){
		printf(" %f %f %f \n",b[i],b[i+1],b[i+2]);
	}

	printf("\nHomogeneous weighting vector is \n");
	for (i = 1; i <= npts; i++){
		printf(" %f ", h[i]);
	}
	printf("\n");
	

	printf("\nCurve points\n");

	for (i = 1; i <= 3*p1; i=i+3){
		printf(" %f %f %f \n",p[i],p[i+1],p[i+2]);
	}
	//getch();

	system("pause");
}