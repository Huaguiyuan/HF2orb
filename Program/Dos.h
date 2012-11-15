/** \ingroup HF */

/*! \file Dos.h
 *
 *  Calculate DOS
 *
 */
#ifndef DOS_H
#define DOS_H
#include "Utils.h"

namespace HF{
	template<typename EngineParamsType, typename HamType, typename LatticeType, typename FieldType>
	class Dos {
		typedef typename HamType::ComplexType ComplexType;
		typedef typename HamType::CMatrixType HamMatrixType;
		
	public:
		Dos(const EngineParamsType& engineParams, HamType& hamiltonian, LatticeType& lattice) :
			engineParams_(engineParams), hamiltonian_(hamiltonian), lattice_(lattice), lx_(engineParams_.lx), ly_(engineParams_.ly), sites_(lx_*ly_), hilbertSize_(hamiltonian_.getLength()), numOrb(hamiltonian_.getOrbs()), numSpin(hamiltonian_.getSpins()), Pi(3.1415926), nTBC(64)
		{
		}
				
		void calcDOS(std::ofstream& fout)
		{
			hamiltonian_.BuildHam();

			size_t totalSize = hilbertSize_ * nTBC * nTBC;
			std::vector<FieldType> eigTmp(hilbertSize_);
			std::vector<FieldType> eigAll(totalSize);

			long count = 0;

			fout << "#DOS" <<std::endl;
			//for(int ikx=lx_/2;ikx>-lx_/2;ikx--)
			for(int i_delta_kx=0;i_delta_kx<nTBC;i_delta_kx++) {

					//for(int iky=ly_/2;iky>-ly_/2;iky--)
				for(int i_delta_ky=0;i_delta_ky<nTBC;i_delta_ky++) {

							//FieldType kx=double(ikx)/lx_;
							//FieldType ky=double(iky)/ly_;
							FieldType phaseX = double(i_delta_kx)/(nTBC);
							FieldType phaseY = double(i_delta_ky)/(nTBC);

							getEigenValues(phaseX,phaseY,eigTmp);
							for(int ix=0;ix<hilbertSize_;ix++) 
								eigAll[count*hilbertSize_+ix]=eigTmp[ix];
							count++;			
					}
			}
			long size = totalSize;
			if(count*hilbertSize_ != size) {
				std::cout<<"sizeofeigAll: "<<eigAll.size()<<" count: "<<count<<" totalSize: "<<totalSize<<std::endl;
				throw std::runtime_error ("ERROR: calcDOS ");
			}
			//else fout << "EigenValuesAll:\n"<<eigAll<<std::endl;

			std::vector<int> histogram;
			std::vector<FieldType> leftPoints;
			utils::MakeHistogram(eigAll, 100*nTBC, histogram, leftPoints);
			long sum = 0;
			for(size_t ii = 0; ii < histogram.size()-1; ii++) {
				std::cout << ii << std::endl;
				sum += histogram[ii];
				fout << (leftPoints[ii]+leftPoints[ii+1])/2.0 - engineParams_.mu <<" "<<(1.0*histogram[ii]/totalSize)<<std::endl;
			}
			sum += histogram[histogram.size()-1];
			std::cout<<"sum: "<<sum<<" totalSize: "<<totalSize<<std::endl;
		}
		
	private:

		void getEigenValues(FieldType& phase_x, FieldType& phase_y, std::vector<FieldType>& eigsTBC)
		{
			HamMatrixType hamTBC; 			
			//FieldType mu = engineParams_.mu; 
			hamTBC.resize(hilbertSize_,hilbertSize_);
			bool test = creatHamTBC(hamTBC,phase_x,phase_y);

			if(test) diagonalize(hamTBC,eigsTBC);
			else throw std::runtime_error ("ERROR: creatHamTBC ");
		}

		enum myCrossType {NEG=-1, NOT=0, POS=1};

		bool creatHamTBC(HamMatrixType& hamTBC, FieldType& phase_x, FieldType& phase_y)
		{
		  FieldType ph_x = phase_x;
		  FieldType ph_y = phase_y;	
			for (int ispin1=0;ispin1<numSpin;ispin1++)
			for (int ispin2=0;ispin2<numSpin;ispin2++) {
				for (int iorb1=0;iorb1<numOrb;iorb1++)
				for (int iorb2=0;iorb2<numOrb;iorb2++) {
					for(int isite=0;isite<sites_;isite++)
					for(int jsite=0;jsite<sites_;jsite++) {
						int i = isite+(iorb1+ispin1*numOrb)*sites_;
						int j = jsite+(iorb2+ispin2*numOrb)*sites_;

						myCrossType crossX = Is_Cross_X_Boundary(isite,jsite);
						myCrossType crossY = Is_Cross_Y_Boundary(isite,jsite);

						if ( (crossX==POS) && (crossY==POS) ) hamTBC(i,j) = hamiltonian_.ham(i,j)*(ComplexType(cos(2.0*Pi*(ph_x+ph_y)),sin(2.0*Pi*(ph_x+ph_y))));
						else if	( (crossX==NEG) && (crossY==NEG) ) hamTBC(i,j) = hamiltonian_.ham(i,j)*(ComplexType(cos(2.0*Pi*(ph_x+ph_y)),-sin(2.0*Pi*(ph_x+ph_y))));
						else if ( (crossX==POS) && (crossY==NEG) ) hamTBC(i,j) = hamiltonian_.ham(i,j)*(ComplexType(cos(2.0*Pi*(ph_x-ph_y)),sin(2.0*Pi*(ph_x-ph_y))));
						else if	( (crossX==NEG) && (crossY==POS) ) hamTBC(i,j) = hamiltonian_.ham(i,j)*(ComplexType(cos(2.0*Pi*(ph_x-ph_y)),-sin(2.0*Pi*(ph_x-ph_y))));

						else if ( (crossX==POS) && (crossY==NOT) ) hamTBC(i,j) = hamiltonian_.ham(i,j)*(ComplexType(cos(2.0*Pi*ph_x),sin(2.0*Pi*ph_x)));
						else if ( (crossX==NEG) && (crossY==NOT) ) hamTBC(i,j) = hamiltonian_.ham(i,j)*(ComplexType(cos(2.0*Pi*ph_x),-sin(2.0*Pi*ph_x)));

						else if ( (crossX==NOT) && (crossY==POS) ) hamTBC(i,j) = hamiltonian_.ham(i,j)*(ComplexType(cos(2.0*Pi*ph_y),sin(2.0*Pi*ph_y)));
						else if ( (crossX==NOT) && (crossY==NEG) ) hamTBC(i,j) = hamiltonian_.ham(i,j)*(ComplexType(cos(2.0*Pi*ph_y),-sin(2.0*Pi*ph_y)));

						else if ( (crossX==NOT) && (crossY==NOT) ) hamTBC(i,j) = hamiltonian_.ham(i,j);
					} // end for (site)
				} // end for (orb)
			} // end for (spin)
			
			return true;
		}

		void diagonalize(HamMatrixType& ham, std::vector<FieldType>& eigs)
		{
			char jobz='V';
			utils::diag(ham,eigs,jobz);
			//if (jobz!='V') sort(eigs.begin(), eigs.end(), std::less<FieldType>());
		}
		
		inline FieldType deltaFun(FieldType var)
		{
			return engineParams_.eps[2]/(pow(engineParams_.eps[2],2)+pow(var,2));
		}

		inline void calcComponents(int i, int& ix, int& iy)
		{
			iy = int(i/lx_);
			ix = i % lx_;
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
	
		const EngineParamsType& engineParams_;		
		HamType& hamiltonian_;	
		LatticeType& lattice_;
		int lx_;
		int ly_;
		int sites_;
		int hilbertSize_;
		const int numOrb;
		const int numSpin;
		const FieldType Pi;
		const int nTBC;

	}; // AKW
	
} // namespace HF

#endif
