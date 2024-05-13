/*
// ppgasp_post.cpp
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

double log_marginal_lik_ppgasp(const Eigen::VectorXd param,double nugget, const bool nugget_est, mxArray *R0[], const MapMat & X,const bool zero_mean,const MapMat & output, const MapMat & kernel_type,const MapMat &  alpha);
//double log_ref_marginal_post_ppgasp(const Eigen::VectorXd param,double nugget, bool nugget_est,  mxArray *R0[], const MapMat & X,const bool zero_mean,const MapMat & output, const MatrixXd kernel_type,const MapMat & alpha);
double log_approx_ref_prior(const Eigen::VectorXd param,double nugget, bool nugget_est, const Eigen::VectorXd CL,const double a,const double b );

Eigen::VectorXd log_marginal_lik_deriv_ppgasp(const Eigen::VectorXd param,double nugget,  bool nugget_est, mxArray *R0[], const MapMat & X,const bool zero_mean,const MapMat & output, const  Eigen::VectorXd kernel_type,const Eigen::VectorXd alpha);

Eigen::VectorXd log_approx_ref_prior_deriv(const Eigen::VectorXd param,double nugget, bool nugget_est, const Eigen::VectorXd CL,const double a,const double b );

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]){
  
  //if (nrhs != 1)
   // mexErrMsgTxt("Wrong number of input arguments.\n");
  
  //if (nlhs > 1)   
  //  mexErrMsgTxt("Too many output argumnents.\n");
  
//#define A_IN  prhs[0]
//  #define B_OUT plhs[0]
     
    auto param=Eigen::matlab2Eigen<double> (prhs[0]);

    double nugget = (double)*mxGetPr(prhs[1]);
    bool nugget_est = (bool)*mxGetPr(prhs[2]);

    int  p;
    if(!nugget_est){
        p= param.size();
    }else{
        p= param.size()-1;
    }

    
    //form  a  3d array pointer for R0
    //there must be  a better way  to do  it
    //const mwSize* dim_num = mxGetDimensions (prhs[3]);
        
    mxArray *R0[p];
    array_pointer_to_3d_array_pointer(R0, prhs[3],p);
    // done for the list item
    auto X=Eigen::matlab2Eigen<double> (prhs[4]);

    bool zero_mean = (bool)*mxGetPr(prhs[5]);
    
    auto output=Eigen::matlab2Eigen<double> (prhs[6]);
    auto kernel_type = Eigen::matlab2Eigen<double> (prhs[7]);
    auto alpha = Eigen::matlab2Eigen<double> (prhs[8]);

    auto CL=Eigen::matlab2Eigen<double> (prhs[9]);
    double a=(double)*mxGetPr(prhs[10]);
    double b=(double)*mxGetPr(prhs[11]);

    /*
    std::cout << param << std::endl;
    std::cout << nugget << std::endl;
    std::cout << nugget_est << std::endl;
    std::cout << R0 << std::endl;
    std::cout << X << std::endl;
    std::cout << zero_mean << std::endl;
    std::cout << output << std::endl;
    std::cout << kernel_type << std::endl;
    std::cout << alpha << std::endl;
*/
    
    double log_lik=log_marginal_lik_ppgasp(param, nugget, nugget_est, R0,  X, zero_mean,  output,  kernel_type, alpha );
    double log_prior=log_approx_ref_prior(param,nugget,nugget_est,CL,a,b);
    //double log_mar_post=log_lik+log_prior;
    
    Eigen::VectorXd log_lik_deriv=log_marginal_lik_deriv_ppgasp(param, nugget,   nugget_est, R0,  X, zero_mean, output, kernel_type, alpha);
    Eigen::VectorXd log_prior_deriv=log_approx_ref_prior_deriv(param,nugget,nugget_est,CL,a,b);
  
    //Eigen::VectorXd log_post_deriv=log_lik_deriv+log_prior_deriv;
    
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
		
    double *log_mar_post;
    log_mar_post=mxGetPr(plhs[0]);
    *log_mar_post=log_lik+log_prior;
    
            
    plhs[1] = mxCreateDoubleMatrix(param.size(), 1, mxREAL);
   
    auto log_post_deriv = Eigen::matlab2Eigen<double> (plhs[1]);

    log_post_deriv=log_lik_deriv+log_prior_deriv;

}




//marginal likelihood for  pp gasp
//I change zero_mean to be bool type
double log_marginal_lik_ppgasp(const Eigen::VectorXd param,double nugget, const bool nugget_est, mxArray *R0[], const MapMat & X,const bool zero_mean,const MapMat & output, const MapMat & kernel_type,const MapMat & alpha ){
  Eigen::VectorXd beta;
  double nu=nugget;
  int k=output.cols();
  int param_size=param.size();
  if(!nugget_est){
    beta= param.array().exp().matrix();
    // nu=0;
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  
  int num_obs=output.rows();
  Eigen::MatrixXd R= separable_multi_kernel(R0,beta, kernel_type,alpha);
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
  
  Eigen::LLT<MatrixXd> lltOfR(R);             // compute the cholesky decomposition of R called lltofR
  Eigen::MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
  
  if(zero_mean){
    
   Eigen::MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    
    
    
    double log_S_2=0;
    
    for(int loc_i=0;loc_i<k;loc_i++){
      log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0));
    }
    
    //double log_S_2=log(S_2);
    
    return (-k*(L.diagonal().array().log().matrix().sum())-(num_obs)/2.0*log_S_2);
    
  }else{
    
    int q=X.cols();
    
    Eigen::MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X)); //one forward and one backward to compute R.inv%*%X
    Eigen::MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X; //Xt%*%R.inv%*%X
    
    Eigen::LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X); // cholesky decomposition of Xt_R_inv_X called lltOfXRinvX
    Eigen::MatrixXd LX = lltOfXRinvX.matrixL();  //  retrieve factor LX  in the decomposition 
    Eigen::MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));          //compute  R_inv_X_Xt_R_inv_X_inv_Xt_R_inv through one forward and one backward solve
    Eigen::MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    
    
    double log_S_2=0;
    
    for(int loc_i=0;loc_i<k;loc_i++){
      log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0));
    }
    // double log_S_2=log(S_2);
    
    //MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
    //double log_S_2=log(S_2(0,0));
    return -k*(L.diagonal().array().log().matrix().sum())-k*(LX.diagonal().array().log().matrix().sum())-(num_obs-q)/2.0*log_S_2;
  }
  
  
}


Eigen::VectorXd log_marginal_lik_deriv_ppgasp(const Eigen::VectorXd param,double nugget,  bool nugget_est, mxArray *R0[], const MapMat & X,const bool zero_mean,const MapMat & output, const  Eigen::VectorXd kernel_type,const  Eigen::VectorXd alpha){
  
  Eigen::VectorXd beta;
  double nu=nugget;
  int k=output.cols();
  int param_size=param.size();
  if(nugget_est==false){//not sure about the logical stuff
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  int p=beta.size();
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta,kernel_type,alpha);
  MatrixXd R_ori=  R;  // this is the one without the nugget
  
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
  
  LLT<MatrixXd> lltOfR(R);
  MatrixXd L = lltOfR.matrixL();
  VectorXd ans=VectorXd::Ones(param_size);
  
  //String kernel_type_ti;
  
  if(zero_mean){
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    //MatrixXd S_2= (yt_R_inv*output);
    
    //double log_S_2=log(S_2(0,0));
    VectorXd S_2_vec=VectorXd::Zero(k);
    
    for(int loc_i=0;loc_i<k;loc_i++){
      S_2_vec[loc_i]=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0);
      
    }
    MatrixXd dev_R_i;
    MatrixXd Vb_ti;
    //allow different choices of kernels
    for(int ti=0;ti<p;ti++){
      //kernel_type_ti=kernel_type[ti];
      auto R0_ti = Eigen::matlab2Eigen<double> (R0[ti]);    
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0_ti,R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0_ti,R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0_ti,R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }
      Vb_ti=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      
      
      double ratio=0;
      
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Vb_ti*(yt_R_inv.transpose()).col(loc_i) )(0,0))/S_2_vec[loc_i];
      }
      ans[ti]=-0.5*k*Vb_ti.diagonal().sum()+num_obs/2.0*ratio;
      
      //ans[ti]=-0.5*Vb_ti.diagonal().sum()+(num_obs/2.0*output.transpose()*Vb_ti*yt_R_inv.transpose()/ S_2(0,0))(0,0) ;  
    }
    //the last one if the nugget exists
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Vb_ti=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      
      double ratio=0;
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Vb_ti*(yt_R_inv.transpose()).col(loc_i))(0,0))/S_2_vec[loc_i];
      }
      ans[p]=-0.5*k*Vb_ti.diagonal().sum()+num_obs/2.0*ratio;
      //ans[p]=-0.5*Vb_ti.diagonal().sum()+(num_obs/2.0*output.transpose()*Vb_ti*yt_R_inv.transpose()/ S_2(0,0))(0,0); 
      
    }
    
  }else{
    int q=X.cols();
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X;
    
    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X);
    MatrixXd LX = lltOfXRinvX.matrixL();
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose();
    MatrixXd dev_R_i;
    MatrixXd Wb_ti;
    //allow different choices of kernels
    
    
    VectorXd S_2_vec=VectorXd::Zero(k);
    
    for(int loc_i=0;loc_i<k;loc_i++){
      S_2_vec[loc_i]=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0);
    }
    
    
    // double log_S_2=0;
    
    //for(int loc_i=0;loc_i<k;loc_i++){
    //  log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0));
    //}
    
    
    for(int ti=0;ti<p;ti++){
      //kernel_type_ti=kernel_type[ti];
        auto R0_ti = Eigen::matlab2Eigen<double> (R0[ti]);    

      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0_ti,R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0_ti,R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0_ti,R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }
      Wb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      
      double ratio=0;
      
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Wb_ti.transpose()*(yt_R_inv.row(loc_i).transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i)))(0,0))/S_2_vec[loc_i];
      }
      
      
      
      //MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
      
      //MatrixXd Q_output= yt_R_inv.transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output;
      
      ans[ti]=-0.5*k*Wb_ti.diagonal().sum()+(num_obs-q)/2.0*ratio; 
      
      
      //ans[ti]=-0.5*Wb_ti.diagonal().sum()+(num_obs-q)/2.0*(output.transpose()*Wb_ti.transpose()*Q_output/S_2(0,0))(0,0); 
    }
    
    
    
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Wb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      
      //double S_2_dev=0;
      double ratio=0;
      
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Wb_ti.transpose()*(yt_R_inv.row(loc_i).transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i)))(0,0))/S_2_vec[loc_i];
      }
      ans[p]=-0.5*k*Wb_ti.diagonal().sum()  +(num_obs-q)/2.0*ratio; 
      
      //ans[p]=-0.5*Wb_ti.diagonal().sum()+(num_obs-q)/2.0*(output.transpose()*Wb_ti.transpose()*Q_output/S_2(0,0))(0,0); 
    }
    
    
  }
  return ans;
  
}


double log_approx_ref_prior(const Eigen::VectorXd param,double nugget, bool nugget_est, const Eigen::VectorXd CL,const double a,const double b ){

  Eigen::VectorXd beta;
  double nu=nugget;
  int param_size=param.size();
  if(!nugget_est){
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  double t=CL.cwiseProduct(beta).sum()+nu;
  double part_I=-b*t;
  double part_II= a*log(t);
  return part_I+part_II;
}

Eigen::VectorXd log_approx_ref_prior_deriv(const Eigen::VectorXd param,double nugget, bool nugget_est, const Eigen::VectorXd CL,const double a,const double b ){

  Eigen::VectorXd beta;
  Eigen::VectorXd return_vec;
  double nu=nugget;
  int param_size=param.size();
  if(!nugget_est){//not sure about the logical stuff. Previously (nugget_est==false)
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }

  //  double a=1/2.0;//let people specify
  // double b=(a+beta.size())/2.0;
  double t=CL.cwiseProduct(beta).sum()+nu;

  if(!nugget_est){
    return_vec=(a*CL/t- b*CL);
  }else{
    Eigen::VectorXd CL_1(param_size);
    CL_1.head(param_size-1)=CL;
    CL_1[param_size-1]=1;
    return_vec=(a*CL_1/t- b*CL_1);
  }
  return return_vec;

}

// it is  useful  in the future for  the reference  prior
/*

double log_ref_marginal_post_ppgasp(const Eigen::VectorXd param,double nugget, bool nugget_est,  mxArray *R0[], const MapMat & X,const bool zero_mean,const MapMat & output, const VectorXd kernel_type,const MapMat & alpha){
  
  Eigen::VectorXd beta;
  double nu=nugget;
  int k=output.cols();
  int param_size=param.size();
  if(nugget_est==false){//not sure about the logical stuff
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  int p=beta.size();
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta,kernel_type,alpha);
  MatrixXd R_ori=  R;  // this is the one without the nugget
  
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
  
  LLT<MatrixXd> lltOfR(R);
  MatrixXd L = lltOfR.matrixL();
  
  // String kernel_type_ti;
  
  if(zero_mean){
    
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    //MatrixXd S_2= (yt_R_inv*output);
    //double log_S_2=log(S_2(0,0));
    
    double log_S_2=0;
    
    for(int loc_i=0;loc_i<k;loc_i++){
      log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0));
    }
    
   // double log_S_2=log(S_2);
    
    VectorXd ans=VectorXd::Ones(param_size);
    MatrixXd dev_R_i;
    MatrixXd IR(param_size+1,param_size+1);
    IR(0,0)=num_obs;
    
    
    //I make  it a matrix to avoid the awkard definition for pointer array 
    MatrixXd Vb=MatrixXd::Zero(num_obs,num_obs*param_size);

    //List Vb(param_size);
    //allow different choices of kernels
    for(int ti=0;ti<p;ti++){
      //  kernel_type_ti=kernel_type[ti];
      auto R0_ti = Eigen::matlab2Eigen<double> (R0[ti]);    

      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0_ti,R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0_ti,R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0_ti,R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }
      //Vb[ti]=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      Vb.block(0,ti*num_obs,num_obs,num_obs )=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      //Vb.block(0,(ti-1)*num_obs,n,n )=
    }
    
    //the last one if the nugget exists
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      //Vb[param_size-1]=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
       Vb.block(0,p*num_obs,num_obs,num_obs )=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));

    }
    // int q=X.cols();
    
    
    for(int i=0;i<param_size;i++){
      MatrixXd Vb_i=Vb.block(0,i*num_obs,num_obs,num_obs );
            //  Vb[i];
      IR(0,i+1)=IR(i+1,0)= Vb_i.trace();
      for(int j=0;j<param_size;j++){
       MatrixXd Vb_j=Vb.block(0,j*num_obs,num_obs,num_obs );
        IR(i+1,j+1)=IR(j+1,i+1)=(Vb_i*Vb_j).trace();
        
      }
    }
    
    LLT<MatrixXd> lltOfIR(IR);
    MatrixXd LIR = lltOfIR.matrixL();
    
    return (-k*(L.diagonal().array().log().matrix().sum())-(num_obs)/2.0*log_S_2+ LIR.diagonal().array().log().matrix().sum());
  }else{
    int q=X.cols();
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X;
    
    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X);
    MatrixXd LX = lltOfXRinvX.matrixL();
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose();
    //MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
    
    double log_S_2=0;
    
    for(int loc_i=0;loc_i<k;loc_i++){
      log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0));
    }
    
    //double log_S_2=log(S_2);
    
   // MatrixXd Q_output= yt_R_inv.transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output;
    MatrixXd dev_R_i;

   // List Wb(param_size);
   
    MatrixXd Wb=MatrixXd::Zero(num_obs,num_obs*param_size);

    
    for(int ti=0;ti<p;ti++){
        auto R0_ti = Eigen::matlab2Eigen<double> (R0[ti]);    

      //  kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0_ti,R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0_ti,R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0_ti,R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }
        
        Wb.block(0,ti*num_obs,num_obs,num_obs )=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      //Wb[ti]=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
    }
    
    
    //the last one if the nugget exists
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
     // Wb[param_size-1]=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
         Wb.block(0,p*num_obs,num_obs,num_obs )=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;

    }
    MatrixXd IR(param_size+1,param_size+1);
    IR(0,0)=num_obs-q;
    for(int i=0;i<param_size;i++){
      //MatrixXd Wb_i=Wb[i];
      MatrixXd Wb_i= Wb.block(0,i*num_obs,num_obs,num_obs );

      IR(0,i+1)=IR(i+1,0)= Wb_i.trace();
      for(int j=0;j<param_size;j++){
        MatrixXd Wb_j= Wb.block(0,j*num_obs,num_obs,num_obs );

        //MatrixXd Wb_j=Wb[j];
        IR(i+1,j+1)=IR(j+1,i+1)=(Wb_i*Wb_j).trace();
        
      }
    }
    
    LLT<MatrixXd> lltOfIR(IR);
    MatrixXd LIR = lltOfIR.matrixL();
    
   // double log_S_2=log(S_2(0,0));
    
    return (-k*(L.diagonal().array().log().matrix().sum())-k*(LX.diagonal().array().log().matrix().sum())-(num_obs-q)/2.0*log_S_2+ LIR.diagonal().array().log().matrix().sum());
  }
  //  return (-(L.diagonal().array().log().matrix().sum())-(LX.diagonal().array().log().matrix().sum())-(num_obs-q)/2.0*log_S_2+1/2.0*log(IR.determinant()) );
}

*/
