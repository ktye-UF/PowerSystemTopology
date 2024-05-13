/*
// construct_ppgasp.cpp
//  Mengyang Gu, July 2019
//
// Copyright (C)  2019 Mengyang Gu
//
// This file is a part of the RobustGaSP Package available at Matlab
//
// The R version of the RobustGaSP Package is available at CRAN
//
// RobustGaSP is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RobustGaSP is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// in file.path(R.home("share"), "licenses").  If not, see
// <http://www.gnu.org/licenses/>.
*/


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
     

    Eigen::VectorXd beta=Eigen::matlab2Eigen<double> (prhs[0]);
    int  p=beta.size();

    double nu = (double)*mxGetPr(prhs[1]);


    //const mwSize* dim_num = mxGetDimensions (prhs[2]);
        
    mxArray *R0[p];
    array_pointer_to_3d_array_pointer(R0, prhs[2],p);
    
    Eigen::MatrixXd X=Eigen::matlab2Eigen<double> (prhs[3]);
    bool zero_mean = (bool)*mxGetPr(prhs[4]);
    Eigen::MatrixXd output=Eigen::matlab2Eigen<double> (prhs[5]);
    Eigen::VectorXd kernel_type=Eigen::matlab2Eigen<double> (prhs[6]);
    Eigen::VectorXd alpha=Eigen::matlab2Eigen<double> (prhs[7]);

    
      //  std::cout << beta << std::endl;
      //  std::cout << kernel_type << std::endl;
      //  std::cout << alpha << std::endl;

    //Eigen::MatrixXd R=separable_multi_kernel (R0, beta, kernel_type, alpha);

    //int n=R.cols();
    int num_obs=output.rows();
    int k=output.cols();
    MatrixXd R= separable_multi_kernel(R0,beta,kernel_type,alpha);
    R=R+nu*MatrixXd::Identity(num_obs,num_obs);  // nu could be zero or nonzero

    Eigen::LLT<MatrixXd> lltOfR(R);             // compute the cholesky decomposition of R called lltofR
    //Eigen::MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
  
    plhs[0] = mxCreateDoubleMatrix(num_obs,num_obs, mxREAL);

    auto L = Eigen::matlab2Eigen<double> (plhs[0]);
    L=lltOfR.matrixL();
          //  std::cout << L << std::endl;

      if(zero_mean){
         plhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
		 plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
              
         double *ptr_1;
         double *ptr_2;
              
         ptr_1=mxGetPr(plhs[1]);
         ptr_2=mxGetPr(plhs[2]);
         *ptr_1=0;
         *ptr_2=0;
         MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
         
         plhs[3]=mxCreateDoubleMatrix(k,1,mxREAL);

         VectorXd S_2_all=Eigen::matlab2Eigen<double> (plhs[3]);
         S_2_all=VectorXd::Zero(k);
    
         for(int loc_i=0;loc_i<k;loc_i++){
           S_2_all[loc_i]=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0);
         }
         S_2_all=S_2_all/num_obs;
    
  }else{
    int q=X.cols();
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X)); //one forward and one backward to compute R.inv%*%X
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X; //Xt%*%R.inv%*%X
    
    Eigen::LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X); // cholesky decomposition of Xt_R_inv_X called lltOfXRinvX
    
    plhs[1]=mxCreateDoubleMatrix(q,q,mxREAL);

    auto LX = Eigen::matlab2Eigen<double> (plhs[1]);
    LX=lltOfXRinvX.matrixL();  //  retrieve factor LX  in the decomposition 
    
    //list_return[1]=LX; //second element to return
    
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    MatrixXd Xt_R_inv_y= X.transpose()*yt_R_inv.transpose();
    
    plhs[2]=mxCreateDoubleMatrix(q,k,mxREAL);

    auto theta_hat=Eigen::matlab2Eigen<double> (plhs[2]);
    theta_hat=LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(Xt_R_inv_y)); 
    //list_return[2]=theta_hat;
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));          //compute  R_inv_X_Xt_R_inv_X_inv_Xt_R_inv through one forward and one backward solver
    //MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
    plhs[3]=mxCreateDoubleMatrix(k,1,mxREAL);

    auto S_2_all=Eigen::matlab2Eigen<double> (plhs[3]);
    S_2_all=MatrixXd::Zero(k,1);

    //VectorXd S_2_all=VectorXd::Zero(k);
    
    for(int loc_i=0;loc_i<k;loc_i++){
      S_2_all(loc_i,0)=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0);
    }
    S_2_all=S_2_all/num_obs;

    //list_return[3]=S_2_all/(num_obs-q);
  }


}
