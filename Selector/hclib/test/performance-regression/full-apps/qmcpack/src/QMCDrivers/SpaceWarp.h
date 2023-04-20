#ifndef QMCPLUSPLUS_QMC_SPACEWARP_H
#define QMCPLUSPLUS_QMC_SPACEWARP_H
#include "Configuration.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"

//namespace qmcplusplus {

//  struct SinglePtclWarp : public QMCTraits {
//
//    int ncenter;
//    PosType nusum;
//
//    //quantities depending on nuclear position
//    vector<RealType> r,rinv,omega,g;
//    vector<PosType>  dr,d_omega,nu;
//    vector<TensorType> d2omega;
//
//    //quantities depending on Hamiltonian
//    vector<RealType> Jacobian;
//    vector<PosType>  Deltax;
//    vector<TensorType> Jacob_matrix,Jacob_cofactor;
//
//    //Unit Matrix
//    TensorType unitMatrix;
//
//    SinglePtclWarp(int nCenter, int nPsi){
//
//      ncenter=nCenter;
//
//      r.resize(ncenter);
//      rinv.resize(ncenter);
//      omega.resize(ncenter);
//      g.resize(ncenter);
//      dr.resize(ncenter);
//      d_omega.resize(ncenter);
//      nu.resize(ncenter);
//      d2omega.resize(ncenter);
//
//      Jacobian.resize(nPsi);
//      Deltax.resize(nPsi);
//      Jacob_matrix.resize(nPsi);
//      Jacob_cofactor.resize(nPsi);
//
//      unitMatrix=0.e0;
//      for (int i=0; i<3; i++) unitMatrix(i,i)=1.e0;
//
//    }
//
//    RealType warpFunction(RealType r, RealType rinv){ return pow(rinv,4);}
//    RealType warpLogGrad(RealType r, RealType rinv){ return -4*rinv;}
//    RealType warpLogGrad2(RealType r, RealType rinv){ return -rinv;}
//
//    void update_warp_disp(const vector<RealType>& R,const vector<RealType>& Rinv,const vector<PosType>& dR){
//      r=R; rinv=Rinv; dr=dR;
//      RealType d=0.e0;
//      for(int iat=0; iat<ncenter; iat++){
//        omega[iat]=warpFunction(r[iat],rinv[iat]);
//        g[iat]=warpLogGrad(r[iat],rinv[iat]);
//        d+=omega[iat];
//      }
//      d=1.e0/d;
//      for(int iat=0; iat<ncenter; iat++) omega[iat]*=d;
//    }
//
//
//    PosType get_displacement(int ipsi, const vector<PosType>& delta){
//      //compute the warp for all geometries
//      Deltax[ipsi]=0.e0;
//      for(int iat=0; iat<ncenter; iat++)
//        Deltax[ipsi]+=(omega[iat]*delta[iat]);
//      return Deltax[ipsi];
//    }
//
//
//    void update_warp_jacob(){
//      //compute d_omega_I/d_r for each nucleus
//      nusum=0.e0;
//      for(int iat=0; iat<ncenter; iat++){
//        nu[iat]=(omega[iat]*g[iat]*rinv[iat])*dr[iat];
//        nusum+=nu[iat];
//      }
//      for(int iat=0; iat<ncenter; iat++)
//        d_omega[iat]=nu[iat]-omega[iat]*nusum;
//    }
//
//
//    RealType get_Jacobian(int ipsi, const vector<PosType>& delta){
//      Jacob_matrix[ipsi]=unitMatrix;
//      for(int iat=0; iat<ncenter; iat++) Jacob_matrix[ipsi]+= outerProduct(delta[iat],d_omega[iat]);
//      //Compute the cofactor matrix
//      Jacob_cofactor[ipsi] = getCof(Jacob_matrix[ipsi]);
//      //cout << Jacob_cofactor[ipsi] << endl;
//      //Compute the jacobian
//      Jacobian[ipsi] = dot(Jacob_matrix[ipsi].getRow(0),Jacob_cofactor[ipsi].getRow(0));
//      return Jacobian[ipsi];
//    }
//
//
//    void update_warp_grad(){
//      //Compute second derivative of \omega
//      TensorType d_nusum(0.e0);
//      for(int iat=0; iat<ncenter; iat++){
//        RealType g2=warpLogGrad2(r[iat],rinv[iat]);
//        PosType sigma( (g2-rinv[iat])*rinv[iat]*dr[iat]+(1.0e0/omega[iat])*d_omega[iat] );
//        TensorType d_nu( outerProduct(sigma,nu[iat]) );
//        RealType aux=omega[iat]*g[iat]*rinv[iat];
//        for(int ibeta=0; ibeta<3; ibeta++)d_nu(ibeta,ibeta)+=aux;
//        d_nusum+=d_nu;
//        d2omega[iat]=d_nu-outerProduct(d_omega[iat],nusum);
//      }
//
//      for(int iat=0; iat<ncenter; iat++)d2omega[iat]-=omega[iat]*d_nusum;
//    }
//
//
//    PosType get_grad_ln_Jacob(int ipsi,const vector<PosType>& delta){
//      PosType GradlnJacob;
//      for(int igamma=0; igamma<3; igamma++){
//        TensorType d_Jacob_matrix=outerProduct(delta[0],d2omega[0].getRow(igamma));
//        for(int iat=1; iat<ncenter; iat++){
//          d_Jacob_matrix+=outerProduct(delta[iat],d2omega[iat].getRow(igamma));
//        }
//        RealType GradJacob=traceAtB(d_Jacob_matrix,Jacob_cofactor[ipsi]);
//        GradlnJacob[igamma]=GradJacob/Jacobian[ipsi];
//      }
//      return GradlnJacob;
//    }
//
//
//    inline TensorType getCof(TensorType a){
//      TensorType b( a(1,1)*a(2,2)-a(2,1)*a(1,2),
//                    a(1,2)*a(2,0)-a(2,2)*a(1,0),
//                    a(1,0)*a(2,1)-a(2,0)*a(1,1),
//                    a(2,1)*a(0,2)-a(0,1)*a(2,2),
//                    a(2,2)*a(0,0)-a(0,2)*a(2,0),
//                    a(2,0)*a(0,1)-a(0,0)*a(2,1),
//                    a(0,1)*a(1,2)-a(1,1)*a(0,2),
//                    a(0,2)*a(1,0)-a(1,2)*a(0,0),
//                    a(0,0)*a(1,1)-a(1,0)*a(0,1) );
//        return b;
//    }
//
//
//    /*inline TinyVector getCofRow(TensorType a, int i0){
//      TinyVector cof;
//      i1=(i0+1)%3;
//      i2=(i1+1)%3;
//      cof[0] = a(i1,1)*a(i2,2) - a(i2,1)*a(i1,2);
//      cof[1] = a(i1,2)*a(i2,0) - a(i2,2)*a(i1,0);
//      cof[2] = a(i1,0)*a(i2,1) - a(i2,0)*a(i1,1);
//      return cof;
//    }*/
//
//  };
//}

namespace qmcplusplus
{

struct SinglePtclWarp;

struct SpaceWarp : public QMCTraits
{

  DistanceTableData* dtPrimary;

  int nptcl,ncenter,npsi,SizeOfR;

  vector<RealType> r,rinv;

  vector<PosType> dr;

  vector<vector<PosType> > Delta;

  vector<ParticleSet*> PtclRefs;

  vector<SinglePtclWarp*> WarpVector;

  Matrix<RealType> one_ptcl_Jacob;

  inline void resize_one_ptcl_Jacob()
  {
    one_ptcl_Jacob.resize(npsi,nptcl);
  }

  void update_one_ptcl_Jacob(int iel);
  //inline void update_one_ptcl_Jacob(int iel){
  //  for(int ipsi=0; ipsi< npsi; ipsi++)
  //    one_ptcl_Jacob[ipsi][iel]=WarpVector[iel]->Jacobian[ipsi];
  //}

  void initialize(vector<ParticleSet*>& IonSets, DistanceTableData* dtprime);

  void warp_one(int iel, bool require_register );

  PosType get_displacement(int iptcl, int ipsi);

  RealType get_Jacobian(int iptcl, int ipsi);

  TensorType get_Jacob_matrix(int iptcl, int ipsi);

  TensorType get_Jacob_cofactor(int iptcl, int ipsi);

  PosType get_grad_ln_Jacob(int iptcl, int ipsi);

  PosType get_grad_ln_Jacob_num(int iptcl, int ipsi);

  void registerData(vector<ParticleSet*>& plist, PooledData<RealType>& buf);


  void updateBuffer(PooledData<RealType>& buf);
  void copyToBuffer(PooledData<RealType>& buf);
  void copyFromBuffer(PooledData<RealType>& buf);

};
}
#endif
