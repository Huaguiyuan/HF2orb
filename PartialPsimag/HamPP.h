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
				hoppings_(modelParams_.hoppings), uu(modelParams_.U), up(modelParams_.Up), jj(modelParams_.J)
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

				//FieldType nxu=mfParams_.nxu[isite];
				//FieldType nxd=mfParams_.nxd[isite];
				//FieldType nyu=mfParams_.nyu[isite];
				//FieldType nyd=mfParams_.nyd[isite];
				
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

/*				ham(ix,iy)=ham(ix,iy) + uu*nxd + (up - jj/2)*(nyu+ nyd) - 2*jj*(nyu - nyd);
				
				spin1=0; spin2=0;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) + uu*nyd + (up - jj/2)*(nxu + nxd) - 2*jj*(nxu - nxd);
				
				spin1=1; spin2=1;
				orb1=0; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) + uu*nxu + (up - jj/2)*(nyu + nyd) + 2*jj*(nyu - nyd);
				ham(ix,iy)=0;
				spin1=1; spin2=1;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) + uu*nyu + (up - jj/2)*(nxu + nxd) + 2*jj*(nxu - nxd);
*/			
				/*--------- off-diagonal ------------*/
				spin1=0; spin2=0;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) - (up - jj/2)*tpu + 4*jj*tpd + 2*jj*tpu + jj*conj(tpd);
				ham(iy,ix)=conj(ham(ix,iy));

				spin1=0; spin2=1;
				orb1=0; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)= ham(ix,iy) - uu*spx; // - 4*jj*spy; 
				ham(iy,ix)=conj(ham(ix,iy));

				spin1=0; spin2=1;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) - (up - jj/2)*opx - 2*jj*opx; //- jj*opy; //
				ham(iy,ix)=conj(ham(ix,iy));

				spin1=0; spin2=1;
				orb1=1; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) - (up - jj/2)*opy - 2*jj*opy; // - jj*opx; //
				ham(iy,ix)=conj(ham(ix,iy));

				spin1=0; spin2=1;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) - uu*spy; //-4*jj*spx; //
				ham(iy,ix)=conj(ham(ix,iy));

				spin1=1; spin2=1;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;
				ham(ix,iy)=ham(ix,iy) - (up - jj/2)*tpd + 4*jj*tpu + 2*jj*tpd + jj*conj(tpu);
				ham(iy,ix)=conj(ham(ix,iy));

			}
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
		void printout(std::ostream& fout)
		{
			//fout<<"#Hamiltonian\n"<<ham<<std::endl;
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

	}; //Ham
} //HF

#endif
