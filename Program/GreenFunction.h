
/** \ingroup HF */

/*! \file GreenFunction.h
 *
 *  Contains GreenFunction and the calculations of # of particles, spin-spin correlations and 
 *  charge-charge correlations
 *
 */
#ifndef GREENFUNCTION_H
#define GREENFUNCTION_H
#include "Utils.h"

namespace HF{
	template<typename EngineParamsType, typename LatticeType, typename HamType, typename FieldType>
	class GreenFunction {
		typedef typename HamType::ComplexType ComplexType;
	public:
		GreenFunction(const EngineParamsType& engineParams, LatticeType& lattice, HamType& hamiltonian) :
			engineParams_(engineParams), hamiltonian_(hamiltonian), lattice_(lattice), hilbertSize_(hamiltonian_.getLength()), sites_(hamiltonian_.getSites()), numOrb(hamiltonian_.getOrbs())
		{
		}
		
		ComplexType operator()(int lambda1, int lambda2)
		{
			return greenFunction(lambda1,lambda2);
		}
		
		FieldType calcNumber()
		{
			FieldType sum=0;
			for (int i=0;i<hilbertSize_;i++) {
				sum += utils::fermi((hamiltonian_.eigs[i]-engineParams_.mu)*engineParams_.beta);
			}
			return sum;
		}
		
		FieldType calcElectronicEnergy() const
		{
			FieldType sum=0;
			for (int i=0;i<hilbertSize_;i++) {
				sum += hamiltonian_.eigs[i] * utils::fermi((hamiltonian_.eigs[i]-engineParams_.mu)*engineParams_.beta);
			}
			return sum;
				
		}

		void calcSpinCorrelation(std::vector<FieldType> &sc)
		{
			int spin1, spin2, spin3, spin4;
			int i1, i2, i3, i4;			
			sc.resize(sites_);
			for (int d=0;d<sites_;d++) {
				ComplexType tmp=0;
				for (int isite=0;isite<sites_;isite++) {
					for (int orb1=0;orb1<numOrb;orb1++) {
						for (int orb2=0;orb2<numOrb;orb2++) {
							int isite2=lattice_.add(isite, d);
							spin1=spin4=0;
							spin2=spin3=1;
							i1=isite+(orb1+spin1*numOrb)*sites_;
							i2=isite+(orb1+spin2*numOrb)*sites_;
							i3=isite2+(orb2+spin3*numOrb)*sites_;
							i4=isite2+(orb2+spin4*numOrb)*sites_;
							tmp += 2.0 * quartic(i1,i2,i3,i4);

							spin1=spin4=1;
							spin2=spin3=0;
							i1=isite+(orb1+spin1*numOrb)*sites_;
							i2=isite+(orb1+spin2*numOrb)*sites_;
							i3=isite2+(orb2+spin3*numOrb)*sites_;
							i4=isite2+(orb2+spin4*numOrb)*sites_;
							tmp += 2.0 * quartic(i1,i2,i3,i4);

							spin1=spin2=spin3=spin4=0;
							i1=isite+(orb1+spin1*numOrb)*sites_;
							i2=isite+(orb1+spin2*numOrb)*sites_;
							i3=isite2+(orb2+spin3*numOrb)*sites_;
							i4=isite2+(orb2+spin4*numOrb)*sites_;
							tmp += quartic(i1,i2,i3,i4);

							spin1=spin2=1;
							spin3=spin4=0;
							i1=isite+(orb1+spin1*numOrb)*sites_;
							i2=isite+(orb1+spin2*numOrb)*sites_;
							i3=isite2+(orb2+spin3*numOrb)*sites_;
							i4=isite2+(orb2+spin4*numOrb)*sites_;
							tmp -= quartic(i1,i2,i3,i4);

							spin1=spin2=0;
							spin3=spin4=1;
							i1=isite+(orb1+spin1*numOrb)*sites_;
							i2=isite+(orb1+spin2*numOrb)*sites_;
							i3=isite2+(orb2+spin3*numOrb)*sites_;
							i4=isite2+(orb2+spin4*numOrb)*sites_;
							tmp -= quartic(i1,i2,i3,i4);

							spin1=spin2=spin3=spin4=1;
							i1=isite+(orb1+spin1*numOrb)*sites_;
							i2=isite+(orb1+spin2*numOrb)*sites_;
							i3=isite2+(orb2+spin3*numOrb)*sites_;
							i4=isite2+(orb2+spin4*numOrb)*sites_;
							tmp += quartic(i1,i2,i3,i4);

						}
					}
				}
				if(fabs(imag(tmp)) > engineParams_.eps[0]) {std::cout<<"sc:tmp\t"<<real(tmp)<<"\t"<<imag(tmp)<<std::endl; throw std::runtime_error("ERROR!!! sc");}
				sc[d]=real(tmp)/sites_;
			}			
		}
		
		void calcChargeCorrelation(std::vector<FieldType> &cc)
		{
			int spin1, spin2, spin3, spin4;
			int i1, i2, i3, i4;
			FieldType den=calcNumber()/sites_;			
			cc.resize(sites_);
			for (int d=0;d<sites_;d++) {
				ComplexType tmp=0;
				for (int isite=0;isite<sites_;isite++) {
					for (int orb1=0;orb1<numOrb;orb1++) {
						for (int orb2=0;orb2<numOrb;orb2++) {
							int isite2=lattice_.add(isite, d);
							spin1=spin2=0;
							spin3=spin4=0;
							i1=isite+(orb1+spin1*numOrb)*sites_;
							i2=isite+(orb1+spin2*numOrb)*sites_;
							i3=isite2+(orb2+spin3*numOrb)*sites_;
							i4=isite2+(orb2+spin4*numOrb)*sites_;
							tmp += quartic(i1,i2,i3,i4);

							spin1=spin2=0;
							spin3=spin4=1;
							i1=isite+(orb1+spin1*numOrb)*sites_;
							i2=isite+(orb1+spin2*numOrb)*sites_;
							i3=isite2+(orb2+spin3*numOrb)*sites_;
							i4=isite2+(orb2+spin4*numOrb)*sites_;
							tmp += quartic(i1,i2,i3,i4);

							spin1=spin2=1;
							spin3=spin4=0;
							i1=isite+(orb1+spin1*numOrb)*sites_;
							i2=isite+(orb1+spin2*numOrb)*sites_;
							i3=isite2+(orb2+spin3*numOrb)*sites_;
							i4=isite2+(orb2+spin4*numOrb)*sites_;
							tmp += quartic(i1,i2,i3,i4);

							spin1=spin2=1;
							spin3=spin4=1;
							i1=isite+(orb1+spin1*numOrb)*sites_;
							i2=isite+(orb1+spin2*numOrb)*sites_;
							i3=isite2+(orb2+spin3*numOrb)*sites_;
							i4=isite2+(orb2+spin4*numOrb)*sites_;
							tmp += quartic(i1,i2,i3,i4);
						}
					}
				}
				if(fabs(imag(tmp)) > engineParams_.eps[0]) {std::cout<<"cc:tmp\t"<<real(tmp)<<"\t"<<imag(tmp)<<std::endl; throw std::runtime_error("ERROR!!! cc");}
				else cc[d] = real(tmp)/sites_ - den*den;
			}
		}


		ComplexType matrix(int lambda1,int lambda2)
		{
			return hamiltonian_.ham(lambda1,lambda2);
		}
		
		
	private:
		
		ComplexType greenFunction(int lambda1,int lambda2)
		{
			ComplexType sum = 0;
			FieldType beta = engineParams_.beta;
			FieldType mu = engineParams_.mu;
			
			for (int lambda=0;lambda<hilbertSize_;lambda++) 
				sum += conj(hamiltonian_.ham(lambda1,lambda)) * hamiltonian_.ham(lambda2,lambda) *utils::fermi(-beta*(hamiltonian_.eigs[lambda]-mu));
			return sum;
		}

		inline ComplexType quartic(int i1, int i2, int i3, int i4)
		{
			ComplexType tmp=0;
			tmp += (deltaFun(i1,i2) - greenFunction(i2,i1)) * (deltaFun(i3,i4) - greenFunction(i4,i3));
			tmp += greenFunction(i2,i3) * (deltaFun(i1,i4) - greenFunction(i4,i1));
			return tmp;
		}

		inline FieldType deltaFun(int x, int y)
		{
			if(x==y) return 1.0;
			else return 0;
		}


		const EngineParamsType& engineParams_;		
		HamType& hamiltonian_;
		LatticeType& lattice_;
		int hilbertSize_;
		int sites_;
		const int numOrb;
		
	}; // GreenFunction
	
} // namespace HF

#endif
