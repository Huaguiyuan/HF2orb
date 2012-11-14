
/* \ingroup HF */

#ifndef ENGINE_H
#define ENGINE_H
#include "Utils.h"
#include "GreenFunction.h"
#include "Adjustments.h"
#include "Akw.h"
#include "Dos.h"

namespace HF{
	template<typename EngineParamsType, typename ModelParamsType, typename MFParamsType, typename LatticeType, typename HamType, typename FieldType>
	class Engine {
		typedef typename HamType::ComplexType ComplexType;
		typedef HF::GreenFunction<EngineParamsType, LatticeType, HamType, FieldType> GreenFunType;
		typedef Spf::Adjustments<EngineParamsType> AdjustType;
		typedef HF::Akw<EngineParamsType, HamType, LatticeType, FieldType> AkwType;
		typedef HF::Dos<EngineParamsType, HamType, LatticeType, FieldType> DosType;

	public:	
		
		Engine(EngineParamsType& engineParams, ModelParamsType& modelParams, MFParamsType& mfParams, LatticeType& lattice)		
			: engineParams_(engineParams), modelParams_(modelParams), lattice_(lattice), mfParamsOld(mfParams), mfParamsNew(mfParams), hamiltonian(engineParams_, modelParams_, mfParamsOld, lattice_), gf(engineParams_, lattice_, hamiltonian), akw(engineParams_, hamiltonian, lattice_), dos(engineParams_, hamiltonian, lattice_), numOrb(hamiltonian.getOrbs()), fout_(engineParams_.output.c_str()), adjustment(engineParams_,1000)
		{	
			writeHeader();
		}
		
		void run()
		{	
			initialize();
			iterate();
			finalize();			
		}

	private:
		
		void writeHeader()
		{
			fout_<<"#Hartree-Fock Approx. for Multiorbital Hubbard Model\n";
			timePrint("Start", fout_);
			fout_<<engineParams_;
			fout_<<modelParams_;
		}

		void writeFooter()
		{			
			timePrint("Stop", fout_);
			fout_<<"#EOF\n";
		}

		void initialize()
		{
			std::cout<<"Engine.Initialize"<<std::endl;
			mfParamsNew.reset();
		}

		void iterate()
		{
			std::cout<<"Engine.Self-consistent iteration"<<std::endl;
			bool converged = false;
			int counter;
			for(counter=0; (counter < engineParams_.iterMAX) && !converged; counter++){
				hamiltonian.BuildHam();
				hamiltonian.Diagonalize();
				
				if(engineParams_.density>0) adjustMu();
				//printProgress(counter);
				printProgress();
				calMFParams();
				converged = isConverged();
				//std::cout<<"converged "<<converged<<std::endl;
				linearMixing();
			}
			if(converged) {
				std::cout<<"\tconverged!!!"<<std::endl;
				fout_<<"Conv:\nConverged!!!\tStop at: "<<counter<<std::endl;
			}
			else {
				std::cout<<"\tNOT converged!!!"<<std::endl;
				fout_<<"Conv:\nNOT Converged!!!\tStop at: "<<counter<<std::endl;
			}
		}

		void finalize()
		{
			std::cout<<"Engine.Finalization"<<std::endl;
			hamiltonian.BuildHam();
			//hamiltonian.printHam(fout_);
			hamiltonian.Diagonalize();
			hamiltonian.printEigenValues(fout_);
			if(engineParams_.density>0) adjustMu();
			//printProgress(1); 
			calMFParams();
			printEnergy();
			calDensity(fout_);
			calMagnetic(fout_);
//			spinCorrelations(fout_);
//			chargeCorrelations(fout_);
			if(fabs(engineParams_.omegaStep) > engineParams_.eps[0]) { 
				std::cout<<"Calculating DOS"<<std::endl; 
				dos.calcDOS(fout_);
				//akw.calc_TBC_Print(fout_); 
			} 
			mfParamsNew.printout("New", fout_);
			writeFooter();		
		}	
		
		void printProgress(const int &iter)
		{
			//std::cout<<"Iteration: "<<iter<<std::endl;
			//std::cout<<"NumOfElectrons: "<<gf.calcNumber()<<std::endl;
			//std::cout<<"Mu: "<<engineParams_.mu<<std::endl;
			fout_<<"Energy: "<<iter<<"\t"<<calcElectronicEnergy()<<"\t"<<calcConstEnergy()<<"\t"<<calcEnergy()<<std::endl;
		}

		void printEnergy()
		{
			fout_<<"#Mu= "<<engineParams_.mu<<"\n#NumOfElectrons= "<<gf.calcNumber()<<std::endl;
			fout_<<"#ElecEnergy= "<<calcElectronicEnergy()<<"\n#ConstEnergy= "<<calcConstEnergy()<<"\n#TotalEnergy= "<<calcEnergy()<<std::endl;
      //fout_<<"\n#KineticEnergy= "<<calcKineticEnergy()<<"\n#InteractionEnergy= "<<calcInteractionEnergy()<<std::endl;
		}
		
		void printProgress(const char &mark='$')
		{
			std::cout<<mark;
			std::cout.flush();
		}	
		
		void adjustMu()
		{
			engineParams_.mu=adjustment.simpleAdjChemPot(hamiltonian.eigs);
			//std::cout<<"Mu: "<<engineParams_.mu<<std::endl;
		} 

		void calMFParams()
		{
			ComplexType tmp,tmp0;
			int spin1, spin2, orb1, orb2, ix, iy;
			const int volume_= hamiltonian.getSites();
			for(int isite=0; isite < volume_; isite++){
				spin1=0; spin2=0;
				orb1=0; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;		
				tmp= gf(iy,ix);
				if(fabs(imag(tmp))<engineParams_.eps[0]) mfParamsNew.nxu[isite] = real(deltaFun(ix, iy) - tmp);
				else throw std::runtime_error("ERROR!!! nxu should be real");
				 
				spin1=1; spin2=1;
				orb1=0; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;	
				tmp= gf(iy,ix);
				if(fabs(imag(tmp))<engineParams_.eps[0]) mfParamsNew.nxd[isite] = real(deltaFun(ix, iy) - tmp);
				else throw std::runtime_error("ERROR!!! nxd should be real");

				spin1=0; spin2=0;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;		
				tmp= gf(iy,ix);
				if(fabs(imag(tmp))<engineParams_.eps[0]) mfParamsNew.nyu[isite] = real(deltaFun(ix, iy) - tmp);
				else throw std::runtime_error("ERROR!!! nyu should be real");

				spin1=1; spin2=1;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;		
				tmp= gf(iy,ix);
				if(fabs(imag(tmp))<engineParams_.eps[0]) mfParamsNew.nyd[isite] = real(deltaFun(ix, iy) - tmp);
				else throw std::runtime_error("ERROR!!! nyd should be real");

				spin1=0; spin2=0;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;		
				tmp= gf(iy,ix);
				tmp0=gf(ix,iy);
				if(abs(conj(tmp) - tmp0) < engineParams_.eps[0]) mfParamsNew.tpu[isite] = deltaFun(ix, iy) - tmp;
				else throw std::runtime_error("ERROR!!! tpu should be Complex Conjugate of tqu");
				
				spin1=1; spin2=1;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;		
				tmp= gf(iy,ix);
				tmp0=gf(ix,iy);
				if(abs(conj(tmp) - tmp0) < engineParams_.eps[0]) mfParamsNew.tpd[isite] = deltaFun(ix, iy) - tmp;
				else throw std::runtime_error("ERROR!!! tpd should be Complex Conjugate of tqd");

				spin1=0; spin2=1;
				orb1=0; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;		
				tmp= gf(iy,ix);
				tmp0=gf(ix,iy);
				if(abs(conj(tmp) - tmp0) < engineParams_.eps[0]) mfParamsNew.spx[isite] = deltaFun(ix, iy) - tmp;
				else throw std::runtime_error("ERROR!!! spx should be Complex Conjugate of sqx");

				spin1=0; spin2=1;
				orb1=1; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;		
				tmp= gf(iy,ix);
				tmp0=gf(ix,iy);
				if(abs(conj(tmp) - tmp0) < engineParams_.eps[0]) mfParamsNew.spy[isite] = deltaFun(ix, iy) - tmp;
				else throw std::runtime_error("ERROR!!! spy should be Complex Conjugate of sqy");

				spin1=0; spin2=1;
				orb1=0; orb2=1;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;		
				tmp= gf(iy,ix);
				tmp0=gf(ix,iy);
				if(abs(conj(tmp) - tmp0) < engineParams_.eps[0]) mfParamsNew.opx[isite] = deltaFun(ix, iy) - tmp;
				else throw std::runtime_error("ERROR!!! opx should be Complex Conjugate of opx");

				spin1=0; spin2=1;
				orb1=1; orb2=0;
				ix=isite+(orb1+spin1*numOrb)*volume_;
				iy=isite+(orb2+spin2*numOrb)*volume_;		
				tmp= gf(iy,ix);
				tmp0=gf(ix,iy);
				if(abs(conj(tmp) - tmp0) < engineParams_.eps[0]) mfParamsNew.opy[isite] = deltaFun(ix, iy) - tmp;
				else throw std::runtime_error("ERROR!!! opy should be Complex Conjugate of opy");

			}
		}

		bool isConverged()
		{	
			int counter=0;
			const int volume_= hamiltonian.getSites();
			for(int isite=0; isite < volume_; isite++){
				if(fabs(mfParamsNew.nxu[isite] - mfParamsOld.nxu[isite]) > engineParams_.eps[1]) ++counter;
				if(fabs(mfParamsNew.nxd[isite] - mfParamsOld.nxd[isite]) > engineParams_.eps[1]) ++counter;
				if(fabs(mfParamsNew.nyu[isite] - mfParamsOld.nyu[isite]) > engineParams_.eps[1]) ++counter;
				if(fabs(mfParamsNew.nyd[isite] - mfParamsOld.nyd[isite]) > engineParams_.eps[1]) ++counter;
				
				if(abs(mfParamsNew.tpu[isite] - mfParamsOld.tpu[isite]) > engineParams_.eps[1]) ++counter;
				if(abs(mfParamsNew.tpd[isite] - mfParamsOld.tpd[isite]) > engineParams_.eps[1]) ++counter;
				if(abs(mfParamsNew.spx[isite] - mfParamsOld.spx[isite]) > engineParams_.eps[1]) ++counter;
				if(abs(mfParamsNew.spy[isite] - mfParamsOld.spy[isite]) > engineParams_.eps[1]) ++counter;
				if(abs(mfParamsNew.opx[isite] - mfParamsOld.opx[isite]) > engineParams_.eps[1]) ++counter;
				if(abs(mfParamsNew.opy[isite] - mfParamsOld.opy[isite]) > engineParams_.eps[1]) ++counter;
			}
			if(counter>0) return false;
			else return true;
		}

		void simpleUpdate()
		{	
			mfParamsOld = mfParamsNew;	
		}

		void linearMixing()
		{
			double alpha=engineParams_.alpha;
			double beta=1.0 - alpha;
			for(int isite=0; isite < hamiltonian.getSites(); isite++){
				mfParamsOld.nxu[isite] = beta * mfParamsOld.nxu[isite] + alpha * mfParamsNew.nxu[isite];
				mfParamsOld.nxd[isite] = beta * mfParamsOld.nxd[isite] + alpha * mfParamsNew.nxd[isite];
				mfParamsOld.nyu[isite] = beta * mfParamsOld.nyu[isite] + alpha * mfParamsNew.nyu[isite];
				mfParamsOld.nyd[isite] = beta * mfParamsOld.nyd[isite] + alpha * mfParamsNew.nyd[isite];

				mfParamsOld.tpu[isite] = beta * mfParamsOld.tpu[isite] + alpha * mfParamsNew.tpu[isite];
				mfParamsOld.tpd[isite] = beta * mfParamsOld.tpd[isite] + alpha * mfParamsNew.tpd[isite];
				mfParamsOld.spx[isite] = beta * mfParamsOld.spx[isite] + alpha * mfParamsNew.spx[isite];
				mfParamsOld.spy[isite] = beta * mfParamsOld.spy[isite] + alpha * mfParamsNew.spy[isite];
				mfParamsOld.opx[isite] = beta * mfParamsOld.opx[isite] + alpha * mfParamsNew.opx[isite];
				mfParamsOld.opy[isite] = beta * mfParamsOld.opy[isite] + alpha * mfParamsNew.opy[isite];
			}
		}
		
		inline FieldType deltaFun(int x, int y)
		{
			if(x==y) return 1.0;
			else return 0;
		}

		void calDensity(std::ostream& fout)
		{	
			const int volume= hamiltonian.getSites();
			std::vector<FieldType> nx, ny;
			nx.resize(volume);
			ny.resize(volume);
			for(int isite=0; isite < volume; isite++){
				nx[isite]=mfParamsNew.nxu[isite]+mfParamsNew.nxd[isite];
				ny[isite]=mfParamsNew.nyu[isite]+mfParamsNew.nyd[isite];
			}
			fout<<"#Mu: "<<engineParams_.mu<<"\n#NumOfElectrons: "<<gf.calcNumber()<<std::endl;
			fout<<"#LocalChargeDensity\n"<<volume<<std::endl;
			for(int isite=0; isite < volume; isite++){
				fout<<nx[isite]+ny[isite]<<std::endl;
			}
			fout<<"#nx\n"<<nx;
			fout<<"#ny\n"<<ny;
		}

		void calMagnetic(std::ostream& fout)
		{
			const int volume= hamiltonian.getSites();
			const FieldType Pi=3.1415926;
			std::vector<FieldType> mx, my, thetaX, thetaY, phiX, phiY;
			mx.resize(volume);
			my.resize(volume);
			thetaX.resize(volume);
			thetaY.resize(volume);
			phiX.resize(volume);
			phiY.resize(volume);
			FieldType tmp=0;
			for(int isite=0; isite < volume; isite++){
				tmp = mfParamsNew.nxu[isite] - mfParamsNew.nxd[isite];
				mx[isite] = sqrt(tmp*tmp + 4*norm(mfParamsNew.spx[isite]));
				thetaX[isite] = acos(tmp/mx[isite]);
				if(imag(mfParamsNew.spx[isite]) >= 0) phiX[isite] = acos(2*real(mfParamsNew.spx[isite])/(mx[isite]*sin(thetaX[isite])));	
				else phiX[isite] = 2*Pi - acos(2*real(mfParamsNew.spx[isite])/(mx[isite]*sin(thetaX[isite])));	
				
				tmp = mfParamsNew.nyu[isite] - mfParamsNew.nyd[isite];
				my[isite] = sqrt(tmp*tmp + 4*norm(mfParamsNew.spy[isite]));
				thetaY[isite] = acos(tmp/my[isite]);
				if(imag(mfParamsNew.spy[isite]) >= 0) phiY[isite] = acos(2*real(mfParamsNew.spy[isite])/(my[isite]*sin(thetaY[isite])));
				else phiY[isite] = 2*Pi - acos(2*real(mfParamsNew.spy[isite])/(my[isite]*sin(thetaY[isite])));
			}
			fout<<"#MagneticMoment\n"<<volume<<std::endl;
/*			for(int isite=0; isite < volume; isite++){
				fout<<isite<<"\t"<<mx[isite]<<"\t"<<my[isite]<<"\t"<<thetaX[isite]<<"\t"<<thetaY[isite]<<"\t"<<phiX[isite]<<"\t"<<phiY[isite]<<std::endl;
			}
*/
			fout<<"#mx\n"<<mx;
			fout<<"#my\n"<<my;
			fout<<"#thetaX\n"<<thetaX;
			fout<<"#thetaY\n"<<thetaY;
			fout<<"#phiX\n"<<phiX;
			fout<<"#phiY\n"<<phiY;
		}

		void spinCorrelations(std::ostream& fout)
		{
			std::vector<FieldType> sc;
			gf.calcSpinCorrelation(sc);
			fout<<"#SpinCorrelations"<<std::endl;
			fout<<sc<<std::endl;
		}

		void chargeCorrelations(std::ostream& fout)
		{
			std::vector<FieldType> cc;
			gf.calcChargeCorrelation(cc);
			fout<<"#ChargeCorrelations"<<std::endl;
			fout<<cc<<std::endl;
		}

		FieldType calcEnergy() const
		{
			return (gf.calcElectronicEnergy()+mfParamsOld.calcConst());
		}

		FieldType calcElectronicEnergy() const
		{
			return (gf.calcElectronicEnergy());
		}

		FieldType calcConstEnergy() const
		{
			return mfParamsOld.calcConst();
		}

		void timePrint(const std::string& st, std::ostream& fout)
		{			
			time_t t = time(0);
			fout<<"#"<<st<<"Time= "<<ctime(&t);
		}

		void test()
		{
			mfParamsNew.reset();
			std::cout<<"NewParametersOne!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
			mfParamsNew.printout();
			mfParamsNew = mfParamsOld;
			std::cout<<"NewParametersTwo!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
			mfParamsNew.printout();
		}

		EngineParamsType& engineParams_;
		ModelParamsType& modelParams_;
		LatticeType& lattice_;
		
		MFParamsType& mfParamsOld;
		MFParamsType mfParamsNew;
		HamType hamiltonian;
		GreenFunType gf;
		AkwType akw;
		DosType dos;

		const int numOrb;
		std::ofstream fout_;
		AdjustType adjustment;

	}; // Engine
	
} // namespace HF

#endif
