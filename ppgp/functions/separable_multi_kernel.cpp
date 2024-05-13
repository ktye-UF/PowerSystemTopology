//#include <iostream>
#include <Eigen/Eigen>
//using Eigen::MatrixXd;
using namespace Eigen;

#include "mex.h"
#include "shared_functions.cpp"


#pragma comment(lib,"libmx.lib")
#pragma comment(lib,"libmex.lib")
#pragma comment(lib,"libmat.lib")

//  allow  to use some short  cuts
typedef   Eigen::VectorXi        iVec;
typedef   Eigen::Map<iVec>      MapiVec;
typedef   Eigen::MatrixXd         Mat;
typedef   Eigen::Map<Mat>        MapMat;
typedef   Eigen::VectorXd         Vec;
//typedef   Eigen::Map<Vec>        MapVec;


void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]){
     

    
    Eigen::VectorXd beta=Eigen::matlab2Eigen<double> (prhs[1]);
    int  p=beta.size();

    
    //test
    // auto R0_i_ker = Eigen::matlab2Eigen<double> (R0[0]);    

    //std::cout << R0_i_ker << std::endl;
    //std::cout << *R0[0] << std::endl;

    Eigen::VectorXd kernel_type=Eigen::matlab2Eigen<double> (prhs[2]);
    Eigen::VectorXd alpha=Eigen::matlab2Eigen<double> (prhs[3]);
    
     //   std::cout << beta << std::endl;
     //   std::cout << kernel_type << std::endl;
     //   std::cout << alpha << std::endl;

    //Eigen::MatrixXd R=separable_multi_kernel (R0, beta, kernel_type, alpha);

    //int n=R.cols();
    
    const mwSize* dim_num = mxGetDimensions (prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(dim_num[0], dim_num[1], mxREAL);

    auto R = Eigen::matlab2Eigen<double> (plhs[0]); //output
        
    //if(p>1){
                
        mxArray *R0[p];
        array_pointer_to_3d_array_pointer(R0, prhs[0],p);


        R=separable_multi_kernel (R0, beta, kernel_type, alpha);
    //}else{//p=1
        
    //    auto R0=Eigen::matlab2Eigen<double> (prhs[0]);

     //    R=separable_multi_kernel_1d (R0, beta, kernel_type, alpha);

    //}
}


