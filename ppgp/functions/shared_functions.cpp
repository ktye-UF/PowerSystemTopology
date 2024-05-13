/*
// shared_functions.cpp
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


typedef   Eigen::VectorXi        iVec;
typedef   Eigen::Map<iVec>      MapiVec;
typedef   Eigen::MatrixXd         Mat;
typedef   Eigen::Map<Mat>        MapMat;
typedef   Eigen::VectorXd         Vec;
typedef   Eigen::Map<Vec>        MapVec;


//define a way to change mex matrix  to Eigen
namespace Eigen {
template<typename T> Map<Matrix<T, Dynamic, Dynamic, ColMajor>> matlab2Eigen (const mxArray * pMat, bool needTranspose = true) {
    Map< Matrix<T, Dynamic, Dynamic, ColMajor>> matrixMap ( (T*) mxGetPr (pMat), mxGetM (pMat), mxGetN (pMat) );
    return matrixMap;
}
}



Eigen::MatrixXd matern_5_2_funct (const MapMat &d, double beta_i){
  const double cnst = sqrt(5.0);
  Eigen::MatrixXd matOnes = Eigen::MatrixXd::Ones(d.rows(),d.cols());
  Eigen::MatrixXd result = cnst*beta_i*d;
  return ((matOnes + result +
	   result.array().pow(2.0).matrix()/3.0).cwiseProduct((-result).array().exp().matrix()));
  
}

 Eigen::MatrixXd matern_3_2_funct (const MapMat &d, double beta_i){
  const double cnst = sqrt(3.0);
  Eigen::MatrixXd matOnes = Eigen::MatrixXd::Ones(d.rows(),d.cols());
  Eigen::MatrixXd result = cnst*beta_i*d;
  return ((matOnes + result ).cwiseProduct((-result).array().exp().matrix()));
  
}

Eigen::MatrixXd pow_exp_funct (const MapMat &d, double beta_i,double alpha_i){
  
  return (-(beta_i*d).array().pow(alpha_i)).exp().matrix();

}

Eigen::MatrixXd  matern_5_2_deriv(const MapMat & R0_i,  const MatrixXd R, double beta_i){
   
  const double sqrt_5 = sqrt(5.0);

  Eigen::MatrixXd matOnes = Eigen::MatrixXd::Ones(R.rows(),R.cols());
  Eigen::MatrixXd R0_i_2=R0_i.array().pow(2.0).matrix();

  Eigen::MatrixXd part1= sqrt_5*R0_i+10.0/3*beta_i*R0_i_2;
  Eigen::MatrixXd part2=matOnes+sqrt_5*beta_i*R0_i+5.0*pow(beta_i,2.0)*R0_i_2/3.0 ;
  return ((part1.cwiseQuotient(part2)  -sqrt_5*R0_i).cwiseProduct(R));
}
  
Eigen::MatrixXd  matern_3_2_deriv(const MapMat & R0_i,  const MatrixXd R, double beta_i){
   
  const double sqrt_3 = sqrt(3.0);

  return(-sqrt(3)*R0_i.cwiseProduct(R)+sqrt_3*R0_i.cwiseProduct((-sqrt_3*beta_i*R0_i).array().exp().matrix()));
    
}

Eigen::MatrixXd pow_exp_deriv(const MapMat & R0_i,  const MatrixXd R, const double beta_i, const double alpha_i){
    
 return  -(R.array()*(R0_i.array().pow(alpha_i))).matrix()*alpha_i*pow(beta_i,alpha_i-1);

}

void array_pointer_to_3d_array_pointer(mxArray *R0[], const mxArray *array_pointer, const int p){
   

        double* array_3D = mxGetPr(array_pointer);
        const mwSize* Dim3Dmatrix = mxGetDimensions (array_pointer);

        int M= Dim3Dmatrix[0];
        int N =Dim3Dmatrix[1];
        //int p =Dim3Dmatrix[2];
        
        int dim_num [3];
        dim_num[0]=M;
        dim_num[1]=N;
        dim_num[2]=p;
        
        //mxArray *R0[p];
        
        int MN=M*N;
        int iN;
        int lMN;
        int iNj;
        double*  R0_pointer;
        //matrix a text matrix that  contain all information
        for(int l=0;l<p;l++){
           R0[l]= mxCreateDoubleMatrix(M,N,mxREAL);
           double*  R0_pointer=(double*)mxGetData(R0[l]);

           lMN=l*MN;
            
           for(int i=0;i<M;i++){
                iN=i*N;
               for(int j=0;j<N;j++){
                   iNj=iN+j;
                   R0_pointer[iNj]=array_3D[lMN+iNj];
               }
            }          
          }
        
        //plhs=R0;
        
    }

//note here for kernel type I  define it  to  be double for simplicity

Eigen::MatrixXd separable_multi_kernel (mxArray *R0[], Eigen::VectorXd beta,Eigen::VectorXd kernel_type, Eigen::VectorXd alpha ){
  
  // Eigen::MatrixXd R0element = R0[0];
   int Rnrow = mxGetM(R0[0]);
  int Rncol = mxGetN(R0[0]);

  //int Rnrow =dim_num[0];
  //int Rncol =dim_num[1];
 // int p=dim_num[2]; //this one  should be equal to beta.size


  Eigen::MatrixXd R = R.Ones(Rnrow,Rncol);
  //String kernel_type_i_ker;
  for (int i_ker = 0; i_ker < beta.size(); i_ker++){
   // kernel_type_i_ker=kernel_type[i_ker];
      //I transfer the  R0 to the  eigen data  type
      auto R0_i_ker = Eigen::matlab2Eigen<double> (R0[i_ker]);    

    if(kernel_type[i_ker]==3){
      R = (matern_5_2_funct(R0_i_ker,beta[i_ker])).cwiseProduct(R);
    }else if(kernel_type[i_ker]==2){
      R = (matern_3_2_funct(R0_i_ker,beta[i_ker])).cwiseProduct(R);
    }else if(kernel_type[i_ker]==1){
      R = (pow_exp_funct(R0_i_ker,beta[i_ker],alpha[i_ker])).cwiseProduct(R);
    }
  }
  return R;
}

// Eigen::MatrixXd separable_multi_kernel_1d ( const MapMat &R0,  const Eigen::VectorXd beta,const Eigen::VectorXd kernel_type, const Eigen::VectorXd alpha ){
//   
//      Eigen::MatrixXd  R;
//      if(kernel_type[0]==3){
//           R = (matern_5_2_funct(R0,beta[0]));
//      }else if(kernel_type[0]==2){
//           R = (matern_3_2_funct(R0,beta[0]));
//      }else if(kernel_type[0]==1){
//           R = (pow_exp_funct(R0,beta[0],alpha[0]));
//      }
// 
//      return R;
// }
// 
// 
