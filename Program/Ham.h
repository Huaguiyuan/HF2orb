/** \ingroup HF */

/*! \file Ham.h
 *
 *  Build Hamiltonian and diagonalize Hamiltonian
 *
 */
#ifndef HAM_H
#define HAM_H
#include "Utils.h"

namespace HF {
	template<typename EngineParamsType, typename ModelParamsType, typename MFParamsType, typename LatticeType, typename FieldType>
	class Ham{
		public:
	
		typedef std::complex<FieldType> ComplexType;
		typedef psimag::Matrix<ComplexType> CMatrixType;

		CMatrixType ham;
		std::vector<FieldType> eigs;

		Ham(EngineParamsType& engineParams, ModelParamsType& modelParams, MFParamsType& mfParams, LatticeType& lattice)
			: engineParams_(engineParams), modelParams_(modelParams), mfParams_(mfParams), lattice_(lattice),
				lx_(engineParams_.lx), ly_(engineParams_.ly), volume_(lx_*ly_), dof_(numOrb*numSpin*volume_),
				hoppings_(modelParams_.hoppings), uu(modelParams_.U), up(modelParams_.Up), jj(modelParams_.J), Pi(3.1415926)
		{	
			ham.resize(dof_,dof_);		
		}
		
		void BuildHam()
		{
			for (size_t irow=0;irow<ham.n_row();irow++) 
				for (size_t icol = 0; icol < ham.n_col(); icol++) 
					ham(irow,icol)=0;
			/*--------- hoppings ------------*/
			for (int ispin=0;ispin<numSpin;ispin++) 
				for (int iorb = 0; iorb < numOrb; iorb++) 
					for (int isite = 0; isite < volume_; isite++) {
						int ix=isite+(iorb+ispin*numOrb)*volume_;
						for(int idir=0; idir<8;idir++){
							int isite2=lattice_.getNeighbour(isite,idir);
							for (int iorb2 = 0; iorb2 < numOrb; iorb2++) {
								int idir2=idir/2; //--- idir2 : four kinds of directions (X, Y, X+Y, X-Y)
								int iy=isite2+(iorb2+ispin*numOrb)*volume_;
								ham(ix,iy)=ham(ix,iy)+hoppings_[idir2*numOrb*numOrb+iorb*numOrb+iorb2];
							}
						}
					}
			/*--------- interaction ------------*/
			for (int isite = 0; isite < volume_; isite++) {

				FieldType nxu=mfParams_.nxu[isite];
				FieldType nxd=mfParams_.nxd[isite];
				FieldType nyu=mfParams_.nyu[isite];
				FieldType nyd=mfParams_.nyd[isite];
				
				ComplexType tpu=mfParams_.tpu[isite];
				ComplexType tpd=mfParams_.tpd[isite];

				ComplexType spx=mfParams_.spx[isite];
				ComplexType spy=mfParams_.spy[isite];

				ComplexType opx=mfParams_.opx[isite];
				ComplexType opy=mfParams_.opy[isite];

				/*--------- diagonal ------------*/
				int spin1=0, spin2=0;
				int orb1=0, orb2=0;
				int ix=isite+(orb1+spin1*numOrb)*volume_;
				int iy=isite+(orb2+spin2*numOrb)*volume_;

				ham(ix,iy)=ham(ix,iy) + uu*nxd + (up - jj/2)*(nyu + nyd) - 0.5*jj*(nyu - nyd);

				spin1=0; spin2=0;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) + uu*nyd + (up - jj/2)*(nxu + nxd) - 0.5*jj*(nxu - nxd);

				spin1=1; spin2=1;
				orb1=0; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) + uu*nxu + (up - jj/2)*(nyu + nyd) + 0.5*jj*(nyu - nyd);

				spin1=1; spin2=1;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) + uu*nyu + (up - jj/2)*(nxu + nxd) + 0.5*jj*(nxu - nxd);
				/*--------- off-diagonal ------------*/
				spin1=0; spin2=0;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) - (up - jj/2)*tpu + jj*tpd + 0.5*jj*tpu + jj*conj(tpd);
				ham(iy,ix)=conj(ham(ix,iy));

				spin1=0; spin2=1;
				orb1=0; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) - uu*spx - jj*spy;
				ham(iy,ix)=conj(ham(ix,iy));

				spin1=0; spin2=1;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) - (up - jj/2)*opx - 0.5*jj*opx - jj*opy;
				ham(iy,ix)=conj(ham(ix,iy));

				spin1=0; spin2=1;
				orb1=1; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) - (up - jj/2)*opy - 0.5*jj*opy - jj*opx;
				ham(iy,ix)=conj(ham(ix,iy));

				spin1=0; spin2=1;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) - uu*spy - jj*spx;
				ham(iy,ix)=conj(ham(ix,iy));

				spin1=1; spin2=1;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) - (up - jj/2)*tpd + jj*tpu + 0.5*jj*tpd + jj*conj(tpu);
				ham(iy,ix)=conj(ham(ix,iy));

			}
		}

		enum myCrossType {NEG=-1, NOT=0, POS=1};		
		bool BuildHamTBC(CMatrixType& hamT, FieldType& phase_x, FieldType& phase_y)
		{
		  FieldType ph_x = phase_x;
		  FieldType ph_y = phase_y;	

			for (size_t irow=0;irow<hamT.n_row();irow++) 
				for (size_t icol = 0; icol < hamT.n_col(); icol++) 
					hamT(irow,icol)=0;
			//--------- hoppings ------------//
			for (int ispin=0;ispin<numSpin;ispin++) 
				for (int iorb = 0; iorb < numOrb; iorb++) 
					for (int isite = 0; isite < volume_; isite++) {
						int i=isite+(iorb+ispin*numOrb)*volume_;
						for(int idir=0; idir<8;idir++){
							int isite2=lattice_.getNeighbour(isite,idir);
							
							myCrossType crossX = Is_Cross_X_Boundary(isite,isite2);
							myCrossType crossY = Is_Cross_Y_Boundary(isite,isite2);

							for (int iorb2 = 0; iorb2 < numOrb; iorb2++) {
								int idir2=idir/2; //--- idir2 : four kinds of directions (X, Y, X+Y, X-Y)
								int j=isite2+(iorb2+ispin*numOrb)*volume_;
								hamT(i,j)=hamT(i,j)+hoppings_[idir2*numOrb*numOrb+iorb*numOrb+iorb2];

								if(crossX==POS) hamT(i,j) = hamT(i,j)*(ComplexType(cos(2.0*Pi*ph_x),sin(2.0*Pi*ph_x)));
								if(crossX==NEG) hamT(i,j) = hamT(i,j)*(ComplexType(cos(2.0*Pi*ph_x),-sin(2.0*Pi*ph_x)));
								if(crossY==POS) hamT(i,j) = hamT(i,j)*(ComplexType(cos(2.0*Pi*ph_y),sin(2.0*Pi*ph_y)));
								if(crossY==NEG) hamT(i,j) = hamT(i,j)*(ComplexType(cos(2.0*Pi*ph_y),-sin(2.0*Pi*ph_y)));

								//hamT(j,i)=conj(hamT(i,j));

							}
						}
					}
			//--------- interaction ------------//
			for (int isite = 0; isite < volume_; isite++) {

				FieldType nxu=mfParams_.nxu[isite];
				FieldType nxd=mfParams_.nxd[isite];
				FieldType nyu=mfParams_.nyu[isite];
				FieldType nyd=mfParams_.nyd[isite];
				
				ComplexType tpu=mfParams_.tpu[isite];
				ComplexType tpd=mfParams_.tpd[isite];

				ComplexType spx=mfParams_.spx[isite];
				ComplexType spy=mfParams_.spy[isite];

				ComplexType opx=mfParams_.opx[isite];
				ComplexType opy=mfParams_.opy[isite];

				//--------- diagonal ------------//
				int spin1=0, spin2=0;
				int orb1=0, orb2=0;
				int ix=isite+(orb1+spin1*numOrb)*volume_;
				int iy=isite+(orb2+spin2*numOrb)*volume_;

				hamT(ix,iy)=hamT(ix,iy) + uu*nxd + (up - jj/2)*(nyu + nyd) - 0.5*jj*(nyu - nyd);

				spin1=0; spin2=0;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				hamT(ix,iy)=hamT(ix,iy) + uu*nyd + (up - jj/2)*(nxu + nxd) - 0.5*jj*(nxu - nxd);

				spin1=1; spin2=1;
				orb1=0; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				hamT(ix,iy)=hamT(ix,iy) + uu*nxu + (up - jj/2)*(nyu + nyd) + 0.5*jj*(nyu - nyd);

				spin1=1; spin2=1;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				hamT(ix,iy)=hamT(ix,iy) + uu*nyu + (up - jj/2)*(nxu + nxd) + 0.5*jj*(nxu - nxd);
				//--------- off-diagonal ------------//
				spin1=0; spin2=0;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				hamT(ix,iy)=hamT(ix,iy) - (up - jj/2)*tpu + jj*tpd + 0.5*jj*tpu + jj*conj(tpd);
				hamT(iy,ix)=conj(hamT(ix,iy));

				spin1=0; spin2=1;
				orb1=0; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				hamT(ix,iy)=hamT(ix,iy) - uu*spx - jj*spy;
				hamT(iy,ix)=conj(hamT(ix,iy));

				spin1=0; spin2=1;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				hamT(ix,iy)=hamT(ix,iy) - (up - jj/2)*opx - 0.5*jj*opx - jj*opy;
				hamT(iy,ix)=conj(hamT(ix,iy));

				spin1=0; spin2=1;
				orb1=1; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				hamT(ix,iy)=hamT(ix,iy) - (up - jj/2)*opy - 0.5*jj*opy - jj*opx;
				hamT(iy,ix)=conj(hamT(ix,iy));

				spin1=0; spin2=1;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				hamT(ix,iy)=hamT(ix,iy) - uu*spy - jj*spx;
				hamT(iy,ix)=conj(hamT(ix,iy));

				spin1=1; spin2=1;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				hamT(ix,iy)=hamT(ix,iy) - (up - jj/2)*tpd + jj*tpu + 0.5*jj*tpd + jj*conj(tpu);
				hamT(iy,ix)=conj(hamT(ix,iy));

			}
			return true;
		}

		myCrossType Is_Cross_X_Boundary(int isite1, int isite2)
		{			
			std::vector<int> coor_site1(2);
			lattice_.index2Coor(coor_site1,isite1);
			std::vector<int> coor_site2(2);
			lattice_.index2Coor(coor_site2,isite2);
			//cout<<"("<<isite1<<","<<isite2<<")"<<" site1: "<<coor_site1[0]<<coor_site1[0]<<" site2: "<<coor_site2[0]<<coor_site2[0]<<endl;
			if((coor_site1[0]==0 || coor_site1[0]==1) && (coor_site2[0]==(lx_-1) || coor_site2[0]==(lx_-2))) return POS;
			else if((coor_site1[0]==(lx_-1) || coor_site1[0]==(lx_-2)) && (coor_site2[0]==0 || coor_site2[0]==1)) return NEG;
			else return NOT;
		}

		// ------- false: NOT on Boundary ------ true: on Y Boundary ------ //
		myCrossType Is_Cross_Y_Boundary(int isite1, int isite2)
		{			
			std::vector<int> coor_site1(2);
			lattice_.index2Coor(coor_site1,isite1);
			std::vector<int> coor_site2(2);
			lattice_.index2Coor(coor_site2,isite2);
			if((coor_site1[1]==0 || coor_site1[1]==1) && (coor_site2[1]==(ly_-1) || coor_site2[1]==(ly_-2))) return POS;
			else if((coor_site1[1]==(ly_-1) || coor_site1[1]==(ly_-2)) && (coor_site2[1]==0 || coor_site2[1]==1)) return NEG;
			else return NOT;
		}

		void Diagonalize()
		{
			char jobz='V';
			utils::diag(ham,eigs,jobz);
			//if (jobz!='V') sort(eigs.begin(), eigs.end(), std::less<FieldType>());
		}

		bool IsConjugate(){
			bool isConj=true;
			for(int ix;ix<dof_;ix++)
				for(int iy=ix;iy<dof_;iy++){
					if(ham(ix,iy)!=conj(ham(iy,ix))){
						isConj=false;		
					}
				}
			std::cout<<"Hamiltonian"<<ham<<std::endl;
			return isConj;

		}

		void TestDiag()
		{
			std::vector<FieldType> eigsTest;
			utils::diag(ham,eigsTest,'N');
			std::cout<<"EigenValues:\n"<<eigsTest<<std::endl;
		}

		int getLength()
		{
			return dof_;
		}

		int getSites()
		{
			return volume_;
		}

		int getOrbs()
		{
			return numOrb;
		}

		int getSpins()
		{
			return numSpin;
		}

		void printHam(std::ostream& fout)
		{
			fout<<"#Hamiltonian\n"<<ham<<std::endl;
		}
		void printEigenValues(std::ostream& fout)
		{
			fout<<"#EigenValues\n"<<eigs<<std::endl;
		}

		private:
		
		EngineParamsType& engineParams_;
		ModelParamsType& modelParams_;
		MFParamsType& mfParams_;
		LatticeType& lattice_;
		int lx_;
		int ly_;
		int volume_;
		static const int numSpin=2;
		static const int numOrb=2;
		int dof_;
		std::vector<FieldType>& hoppings_; 
		FieldType uu;
		FieldType up;
		FieldType jj;
		const FieldType Pi;

	}; //Ham
} //HF

#endif
