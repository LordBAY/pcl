/**
 * @dilw levmar_test4.cpp
 * @brief Curve fitting
 */
#include <levmar/levmar.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <SQ_utils.h>

void ex5_fx( double *p, double* x,
	     int m, int n, void *data );
void ex5_jac( double* p, double* jac,
	      int m, int n, void* data );


struct myData {
    double* ix;
    double* iy;
    double* iz;
    int num;
};


int main( int argc, char* argv[] ) {
    
 
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud( new pcl::PointCloud<pcl::PointXYZ>() );
    cloud = sampleSQ_uniform( 0.5, 0.9, 0.5, 0.8, 0.8, 50 );
    
    std::cout << "Num points: "<< cloud->points.size() << std::endl;
    
    


    int n = cloud->points.size(); // Size of measurements
    int m = 5; // Size of parameters

    double p[m]; double x[n];

    double opts[LM_OPTS_SZ];
    double info[LM_INFO_SZ];

    opts[0] = LM_INIT_MU;
    opts[1] = 1E-15;
    opts[2] = 1E-15;
    opts[3] = 1E-20;
    opts[4] = LM_DIFF_DELTA;


    struct myData adata;
    adata.ix = new double[n];
    adata.iy = new double[n];
    adata.iz = new double[n];
    adata.num = n;

    for( int i = 0; i < n; ++i ) {
	adata.ix[i] = cloud->points[i].x;
	adata.iy[i] = cloud->points[i].y;
	adata.iz[i] = cloud->points[i].z;
    }


    int i, ret;
    
    // We want the function to be closest to their measurements
    for( int i = 0; i < n; ++i ) {
	x[i] = 1.0;
    }

    // Check
    double xi, yi, zi, n1, n2;
    for( int i = 0; i < n; ++i ) {
	
	xi = adata.ix[i];
	yi = adata.iy[i];
	zi = adata.iz[i];

	n1 = pow(fabs(xi/0.5), 2.0/0.8 ) + pow(fabs(yi/0.9), 2.0/0.8 );
	n2 = pow(fabs(zi/0.5), 2.0/0.8 );
	
	double pin = pow(n1, (0.8/0.8)) + n2;
	if( fabs(pin -1) > 0.0001 ) { std::cout << "["<<i<<"]: Value: "<< pin << std::endl; }
    }



    // Initial values
    p[0] = 0.45; p[1] = 0.6; p[2] = 0.2;
    p[3] = 0.7; p[4] = 0.4;
 
    std::cout << "BEFORE: fit parameter: Axes: "<< p[0]<<", "<<p[1]<<", "<< p[2] << std::endl;
    std::cout << " E1, E2: "<< p[3]<<", "<< p[4] << std::endl;


    double dt; clock_t ts, tf;
    ts = clock();
    ret = dlevmar_der( ex5_fx, ex5_jac,
		       p, x,
		       m, n,
		       1000,
		       NULL, info,
		       NULL, NULL, (void*)&adata );
    tf = clock();
    
    dt = (double)(tf-ts)/(double)CLOCKS_PER_SEC;
    std::cout << "Calculation time: "<< dt << std::endl;
    std::cout << " Levenberg Marquardt returned in "<<info[5]<<" iterations "<<
	", reason: "<< info[6] << " sumsq: "<< info[1] <<"["<<info[0]<<"]"<<std::endl;

    std::cout << "AFTER: Best fit parameter: Axes: "<< p[0]<<", "<<p[1]<<", "<< p[2] << std::endl;
    std::cout << " E1, E2: "<< p[3]<<", "<< p[4] << std::endl;
    return 0;

}

/**
 * @function SQ
 */
void ex5_fx( double *p, double* x,
	     int m, int n, void *data ) {

    struct myData* dptr;
    dptr = (struct myData*) data;

    double n1, n2;
    double xi, yi, zi;

    for( int i = 0; i < n; ++i ) {
	
	xi = dptr->ix[i];
	yi = dptr->iy[i];
	zi = dptr->iz[i];
	
	n1 = pow(fabs(xi/p[0]), 2.0/p[4] ) + pow(fabs(yi/p[1]), 2.0/p[4] );
	n2 = pow(fabs(zi/p[2]), 2.0/p[3] );
	
	x[i] = pow(n1, (p[4]/p[3])) + n2;
    }
    
}

void ex5_jac( double* p, double* jac,
	      int m, int n, void* data ) {

    struct myData* dptr;
    dptr = (struct myData*) data;
    double xi, yi, zi;


    register double t2,t3,t4,t5,t6,t7,t8,t9,t10;
    register double t11, t12, t13, t14, t15, t16;
    register double t17, t18, t19, t20, t21, t22, t23;

    for( int i = 0; i < n; ++i ) {

	xi = dptr->ix[i];
	yi = dptr->iy[i];
	zi = dptr->iz[i];

	t2 = xi*xi;
	t3 = 1.0/p[4];
	t4 = 1.0/p[3];
	t5 = 1.0/(p[0]*p[0]);
	t6 = t2*t5;
	t7 = yi*yi;
	t8 = pow(t6,t3);
	t9 = 1.0/(p[1]*p[1]);
	t10 = t7*t9;
	t11 = pow(t10,t3);
	t12 = t8+t11;
	t13 = p[4]*t4;
	t14 = t13-1.0;
	t15 = pow(t12,t14);
	t16 = t3-1.0;
	t17 = zi*zi;
	t18 = 1.0/(p[2]*p[2]);
	t19 = t17*t18;
	t20 = 1.0/(p[3]*p[3]);
	t21 = log(t12);
	t22 = pow(t12,t13);
	t23 = 1.0/(p[4]*p[4]);
	
	jac[5*i] = 1.0/(p[0]*p[0]*p[0])*t2*t4*pow(t6,t16)*t15*-2.0;
	jac[5*i+1] = 1.0/(p[1]*p[1]*p[1])*t4*t7*pow(t10,t16)*t15*-2.0;
	jac[5*i+2] = 1.0/(p[2]*p[2]*p[2])*t4*t17*pow(t19,t4-1.0)*-2.0;
	jac[5*i+3] = -pow(t19,t4)*t20*log(t19)-p[4]*t20*t21*t22;
	jac[5*i+4] = t4*t21*t22-p[4]*t4*t15*(t8*t23*log(t6)+t11*t23*log(t10));
    }

}

