#ifndef SURVEYOR_KS_TEST_H
#define SURVEYOR_KS_TEST_H

// Copied from https://github.com/tumi8/topas/blob/master/detectionmodules/statmodules/wkp-module/ks-test.cpp

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>

void mMultiply(double *A,double *B,double *C,int m) { int i,j,k; double s;
    for(i=0;i<m;i++) for(j=0; j<m; j++)
        {s=0.; for(k=0;k<m;k++) s+=A[i*m+k]*B[k*m+j]; C[i*m+j]=s;}
}
void mPower(double *A,int eA,double *V,int *eV,int m,int n)
{ double *B;int eB,i;
    if(n==1) {for(i=0;i<m*m;i++) V[i]=A[i];*eV=eA; return;} mPower(A,eA,V,eV,m,n/2); B=(double*)malloc((m*m)*sizeof(double)); mMultiply(V,V,B,m); eB=2*(*eV); if(n%2==0){for(i=0;i<m*m;i++) V[i]=B[i]; *eV=eB;}
    else {mMultiply(A,B,V,m); *eV=eA+eB;}
    if(V[(m/2)*m+(m/2)]>1e140) {for(i=0;i<m*m;i++) V[i]=V[i]*1e-140;*eV+=140;} free(B);
}
double K(int n,double d)
{ int k,m,i,j,g,eH,eQ;
    double h,s,*H,*Q;
//OMIT NEXT LINE IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL
    s=d*d*n; if(s>7.24||(s>3.76&&n>99)) return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
    k=(int)(n*d)+1; m=2*k-1; h=k-n*d; H=(double*)malloc((m*m)*sizeof(double)); Q=(double*)malloc((m*m)*sizeof(double)); for(i=0;i<m;i++) for(j=0;j<m;j++)
            if(i-j+1<0) H[i*m+j]=0; else H[i*m+j]=1;
    for(i=0;i<m;i++) {H[i*m]-=pow(h,i+1); H[(m-1)*m+i]-=pow(h,(m-i));} H[(m-1)*m]+=(2*h-1>0?pow(2*h-1,m):0);
    for(i=0;i<m;i++) for(j=0;j<m;j++)
            if(i-j+1>0) for(g=1;g<=i-j+1;g++) H[i*m+j]/=g;
    eH=0; mPower(H,eH,Q,&eQ,m,n);
    s=Q[(k-1)*m+k-1];
    for(i=1;i<=n;i++) {s=s*i/n; if(s<1e-140) {s*=1e140; eQ-=140;}} s*=pow(10.,eQ); free(H); free(Q); return s;
}


double ks_test (std::vector<double> sample1, std::vector<double> sample2) {

    unsigned int n1, n2, n_approx;
    // sample sizes
    float d;
    // the value of Kolmogorov's statistic
    // for the particular values of sample1 and sample2
    int D, Dmin, Dmax, s;
    // used in computing this value d

    std::vector<double>::iterator it1, it2;

    // Determine sample sizes
    n1 = sample1.size();
    n2 = sample2.size();

    // Calculate a conservative n approximation
    n_approx = (unsigned) ceil(float(n1*n2)/(n1+n2));

    // Sort samples
    std::sort(sample1.begin(), sample1.end());
    std::sort(sample2.begin(), sample2.end());

    // We divide the range 0..1 into n1*n2 intervals of equal size 1/(n1*n2).
    //
    // Each item in sample1 makes the sample c.d.f of sample1
    // jump by a step of n2 intervals.
    // Each item in sample2 makes the sample c.d.f of sample2
    // jump by a step of n1 intervals.
    //
    // For each item we compute D, related to the distance between the two
    // sample c.d.f., s_cdf_1 - s_cdf_2, by:
    //
    //    D/(n1*n2) = s_cdf_1 - s_cdf_2
    //
    // We want to determine:
    //
    //    Dmin/(n1*n2) = min [s_cdf_1 - s_cdf_2] <= 0
    //    Dmax/(n1*n2) = max [s_cdf_1 - s_cdf_2] >= 0
    //
    // And then the value of Kolmorogov's statistic D_n1n2 is just:
    //
    //    D_n1n2 = sup |s_cdf_1 - s_cdf_2|
    //           = max [ |Dmin/(n1*n2)| ; |Dmax/(n1*n2)| ]

    D = 0; Dmin = 0; Dmax = 0;
    it1 = sample1.begin();
    it2 = sample2.begin();

    while ( (it1 != sample1.end()) && (it2 != sample2.end()) ) {
        if (*it1 == *it2) {
            // steps in both sample c.d.f., we need to perform all steps
            // in this point before comparing D to Dmin and Dmax

            s = *it1;
            // perform all steps in s_cdf_1 first
            do {
                D += n2;
                it1++;
            }
            while ( (it1 != sample1.end()) && (*it1 == s) );
            // perform all steps in s_cdf_2 now
            do {
                D -= n1;
                it2++;
            }
            while ( (it2 != sample2.end()) && (*it2 == s) );

            // now adapt Dmin, Dmax if necessary
            if (D > Dmax)
                Dmax = D;
            else if (D < Dmin)
                Dmin = D;

        }
        else if (*it1 < *it2) {
            // step in s_cdf_1, increase D by n2
            D += n2;
            it1++;

            if (D > Dmax)
                Dmax = D;
        }

        else {
            // step in F2, decrease D by n1
            D -= n1;
            it2++;

            if (D < Dmin)
                Dmin = D;
        }
    }

    // For two-sided test, take D = max (|Dmax|, |Dmin|) and compute
    // the value d of Kolmogorov's statistic (two-sided only)

    if (-Dmin > Dmax)
        D = -Dmin;
    else
        D = Dmax;

    // Hence the observed value of Kolmogorov's statistic:
    d = float(D)/(n1*n2);

    // Return p-value
    return 1 - K(n_approx,d);

}

#endif //SURVEYOR_KS_TEST_H
