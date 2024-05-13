/*
// pred_ppgasp.cpp
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

    Eigen::MatrixXd input=Eigen::matlab2Eigen<double> (prhs[2]);

    //Eigen::MatrixXd X=Eigen::matlab2Eigen<double> (prhs[3]);



    //const mwSize* dim_num = mxGetDimensions (prhs[2]);
        
    //int  p=dim_num[2];
    //mxArray *R0[p];
    //array_pointer_to_3d_array_pointer(R0, prhs[2]);
    
    Eigen::MatrixXd X=Eigen::matlab2Eigen<double> (prhs[3]);
    bool zero_mean = (bool)*mxGetPr(prhs[4]);
    Eigen::MatrixXd output=Eigen::matlab2Eigen<double> (prhs[5]);
    Eigen::MatrixXd testing_input=Eigen::matlab2Eigen<double> (prhs[6]);
    Eigen::MatrixXd X_testing=Eigen::matlab2Eigen<double> (prhs[7]);
    Eigen::MatrixXd L=Eigen::matlab2Eigen<double> (prhs[8]);
    auto LX=Eigen::matlab2Eigen<double> (prhs[9]);
    auto theta_hat=Eigen::matlab2Eigen<double> (prhs[10]);
    VectorXd sigma2_hat=Eigen::matlab2Eigen<double> (prhs[11]);

    double qt_025 = (double)*mxGetPr(prhs[12]);
    double qt_975 = (double)*mxGetPr(prhs[13]);

    //const mwSize* dim_num = mxGetDimensions (prhs[14]);
        
    
    mxArray *r0[p];
    array_pointer_to_3d_array_pointer(r0, prhs[14],p);

    Eigen::VectorXd kernel_type=Eigen::matlab2Eigen<double> (prhs[15]);
    Eigen::VectorXd alpha=Eigen::matlab2Eigen<double> (prhs[16]);

    bool mean_only = (bool)*mxGetPr(prhs[17]);

    int num_testing=testing_input.rows();
      //int dim_inputs=input.cols();
    int num_obs=output.rows();
    int k=output.cols();
  
    MatrixXd r= separable_multi_kernel(r0,beta, kernel_type,alpha);

    
    
    MatrixXd rt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(r.transpose()))).transpose();
    VectorXd c_star_star(num_testing);
    MatrixXd rtR_inv_r;
  
     if(zero_mean){
         
        plhs[0] = mxCreateDoubleMatrix(num_testing,k, mxREAL);
        auto MU_testing = Eigen::matlab2Eigen<double> (plhs[0]);
        MU_testing=rt_R_inv*output;
        
        if(!mean_only){

           for(int i_loc=0; i_loc<num_testing;i_loc++){
              rtR_inv_r=(rt_R_inv.row(i_loc)*r.row(i_loc).transpose());
              c_star_star[i_loc]=1+nu-rtR_inv_r(0,0);
            }
            plhs[1] = mxCreateDoubleMatrix(num_testing,k, mxREAL);
            plhs[2] = mxCreateDoubleMatrix(num_testing,k, mxREAL);
            plhs[3] = mxCreateDoubleMatrix(num_testing,k, mxREAL);

            auto pred_LB = Eigen::matlab2Eigen<double> (plhs[1]);
            auto pred_UB = Eigen::matlab2Eigen<double> (plhs[2]);
            auto pred_sd = Eigen::matlab2Eigen<double> (plhs[3]);

            MatrixXd pred_sigma_2_star=MatrixXd::Zero(num_testing,k);


            for(int loc_i=0;loc_i<k;loc_i++){
              pred_sigma_2_star.col(loc_i)=  sigma2_hat[loc_i]*c_star_star.array().abs().matrix();
            }
            pred_LB=MU_testing+pred_sigma_2_star.array().sqrt().matrix()*qt_025;
            pred_UB=MU_testing+pred_sigma_2_star.array().sqrt().matrix()*qt_975;
            pred_sd=(pred_sigma_2_star*(num_obs)/(num_obs-2)).array().sqrt().matrix();

        }
         
     }else{
         
        plhs[0] = mxCreateDoubleMatrix(num_testing,k, mxREAL);
        auto MU_testing = Eigen::matlab2Eigen<double> (plhs[0]);
        MU_testing=X_testing*theta_hat+rt_R_inv*(output-X*theta_hat);

        if(!mean_only){
            int q=X.cols();
            MatrixXd diff2;
            MatrixXd  R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));  
            MatrixXd X_testing_X_R_inv_r_i;
            for(int i=0; i<num_testing;i++){
              X_testing_X_R_inv_r_i=X_testing.row(i)-r.row(i)*R_inv_X;
              diff2=X_testing_X_R_inv_r_i*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(X_testing_X_R_inv_r_i.transpose())));

              rtR_inv_r=(rt_R_inv.row(i)*r.row(i).transpose());
              c_star_star[i]=1+nu-rtR_inv_r(0,0)+diff2(0,0);
            }

            plhs[1] = mxCreateDoubleMatrix(num_testing,k, mxREAL);
            plhs[2] = mxCreateDoubleMatrix(num_testing,k, mxREAL);
            plhs[3] = mxCreateDoubleMatrix(num_testing,k, mxREAL);


            auto pred_LB = Eigen::matlab2Eigen<double> (plhs[1]);
            auto pred_UB = Eigen::matlab2Eigen<double> (plhs[2]);
            auto pred_sd = Eigen::matlab2Eigen<double> (plhs[3]);

            MatrixXd pred_sigma_2_star=MatrixXd::Zero(num_testing,k);
            for(int loc_i=0;loc_i<k;loc_i++){
              pred_sigma_2_star.col(loc_i)=  sigma2_hat[loc_i]*c_star_star.array().abs().matrix();
            }
             //   std::cout << qt_025 << std::endl;
             //   std::cout << qt_975 << std::endl;

              //   std::cout << sigma2_hat << std::endl;

            pred_LB=MU_testing+pred_sigma_2_star.array().sqrt().matrix()*qt_025;
            pred_UB=MU_testing+pred_sigma_2_star.array().sqrt().matrix()*qt_975;
            pred_sd=(pred_sigma_2_star*(num_obs)/(num_obs-2)).array().sqrt().matrix();
        }
     }
    
    
    
}
