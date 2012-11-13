/** \ingroup HF */

/*! \file MFParams.h
 *
 *  Contains the Mean-field parameters for the Hubbard model
 *
 */
#ifndef MFPARAMS_H
#define MFPARAMS_H
#include "Utils.h"

namespace HF {
	template<typename EngineParamsType, typename ModelParamsType, typename RandomNumberGeneratorType, typename FieldType>
	class MFParams{ //! Mean-Field Parameters
		public:
	
		typedef std::complex<FieldType> ComplexType;
	
						/* x,y : orbital \alpha (xz) and \beta (yz), 
						    u,d : spin up and down	  	*/
		std::vector<FieldType> nxu;	//nxu=< c_{i,x,u}^+  c_{i,x,u} >
		std::vector<FieldType> nxd;	//nxd=< c_{i,x,d}^+  c_{i,x,d} >
		std::vector<FieldType> nyu;	//nyu=< c_{i,y,u}^+  c_{i,y,u} >
		std::vector<FieldType> nyd;	//nyd=< c_{i,y,d}^+  c_{i,y,d} >

		std::vector<ComplexType> tpu;	//tpu=< c_{i,x,u}^+  c_{i,y,u} >
		std::vector<ComplexType> tpd;	//tpd=< c_{i,x,d}^+  c_{i,y,d} >
		//std::vector<ComplexType> tqu;	//tqu=< c_{i,y,u}^+  c_{i,x,u} >=C.C.(tpu)
		//std::vector<ComplexType> tqd;	//tqd=< c_{i,y,d}^+  c_{i,x,d} >=C.C.(tpd)

		std::vector<ComplexType> spx;	//spx=< c_{i,x,u}^+  c_{i,x,d} >
		std::vector<ComplexType> spy;	//spy=< c_{i,y,u}^+  c_{i,y,d} >
		//std::vector<ComplexType> sqx;	//sqx=< c_{i,x,d}^+  c_{i,x,u} >=C.C.(spx)
		//std::vector<ComplexType> sqy;	//sqy=< c_{i,y,d}^+  c_{i,y,u} >=C.C.(spy)

		std::vector<ComplexType> opx;	//opx=< c_{i,x,u}^+  c_{i,y,d} >
		std::vector<ComplexType> opy;	//opy=< c_{i,y,u}^+  c_{i,x,d} >
		//std::vector<ComplexType> oqx;	//oqx=< c_{i,y,d}^+  c_{i,x,u} >=C.C.(opx)
		//std::vector<ComplexType> oqy;	//oqy=< c_{i,x,d}^+  c_{i,y,u} >=C.C.(opy)

		MFParams(EngineParamsType& engineParams, ModelParamsType& modelParams) 
			: engineParams_(engineParams), modelParams_(modelParams), lx_(engineParams_.lx), ly_(engineParams_.ly), volume_(lx_*ly_), ne(2.0*volume_), uu(modelParams_.U), up(modelParams_.Up), jj(modelParams_.J) 
		{
			initial();
		}
		
		void printout(const std::string& st, std::ostream& fout){
			fout<<"MFParameters: "<<st<<std::endl;
			fout<<"nxu\n"<<nxu<<std::endl;
			fout<<"nxd\n"<<nxd<<std::endl;
			fout<<"nyu\n"<<nyu<<std::endl;
			fout<<"nyd\n"<<nyd<<std::endl;

			fout<<"tpu\n"<<tpu<<std::endl;
			fout<<"tpd\n"<<tpd<<std::endl;
			fout<<"spx\n"<<spx<<std::endl;
			fout<<"spy\n"<<spy<<std::endl;
			fout<<"opx\n"<<opx<<std::endl;
			fout<<"opy\n"<<opy<<std::endl;
		}

		
		void initial()
		{
			int vol=volume_;
			RandomNumberGeneratorType rng_;
			rng_.seed(engineParams_.randomSeed);
			resizeParams(vol);
			for(int iter=0; iter<vol;iter++){
				nxu[iter]=rng_();
				nxd[iter]=rng_();
				nyu[iter]=rng_();
				nyd[iter]=rng_();
				FieldType tmp=nxu[iter]+nxd[iter]+nyu[iter]+nyd[iter];
				nxu[iter]=nxu[iter]/tmp*ne;
					//nxu[iter]=0;
				nxd[iter]=nxd[iter]/tmp*ne;
					//nxd[iter]=0;
				nyu[iter]=nyu[iter]/tmp*ne;
					//nyu[iter]=0;
				nyd[iter]=nyd[iter]/tmp*ne;
					//nyd[iter]=0;

				tpu[iter]=ComplexType(rng_(),rng_());
					//tpu[iter]=ComplexType(0,0);
				tpd[iter]=ComplexType(rng_(),rng_());
					//tpd[iter]=ComplexType(0,0);
				spx[iter]=ComplexType(rng_(),rng_());
					//spx[iter]=ComplexType(0,0);
				spy[iter]=ComplexType(rng_(),rng_());
					//spy[iter]=ComplexType(0,0);
				opx[iter]=ComplexType(rng_(),rng_());
					//opx[iter]=ComplexType(0,0);
				opy[iter]=ComplexType(rng_(),rng_());
					//opy[iter]=ComplexType(1,2);
			}
			
		}

		MFParams &operator = (MFParams& paramsOther)
		{
			nxu = paramsOther.nxu;
			nxd = paramsOther.nxd;
			nyu = paramsOther.nyu;
			nyd = paramsOther.nyd;

			tpu = paramsOther.tpu;
			tpd = paramsOther.tpd;

			spx = paramsOther.spx;
			spy = paramsOther.spy;

			opx = paramsOther.opx;
			opy = paramsOther.opy;

			return *this;
		}
		
		int getLength()
		{
			return volume_;
		}

		void reset()
		{
			for(int i=0; i<volume_; i++){
				nxu[i] = 0;
				nxd[i] = 0;
				nyu[i] = 0;
				nyd[i] = 0;

				tpu[i] = 0;
				tpd[i] = 0;
				spx[i] = 0;
				spy[i] = 0;
				opx[i] = 0;
				opy[i] = 0;
			}
		}

		FieldType calcConst() const
		{
			FieldType sum=0;
			for (int i=0;i<volume_;i++) {
				sum += -uu*(nxu[i]*nxd[i]+nyu[i]*nyd[i]) -(up - jj/2)*(nxu[i]+nxd[i])*(nyu[i]+nyd[i]) +2*jj*(2*real(spx[i]*conj(spy[i]) + spy[i]*conj(spx[i])) + nxu[i]*nyu[i]+nxd[i]*nyd[i]-nxu[i]*nyd[i]-nxd[i]*nyu[i]) + jj*(real(opx[i]*conj(opy[i]) + opy[i]*conj(opx[i])));// + uu*(norm(spx[i])+norm(spy[i])) + (up - jj/2)*(norm(tpu[i])+norm(tpd[i])+norm(opx[i])+norm(opy[i])) - 2*jj*( 2*real(tpu[i]*conj(tpd[i]) + tpd[i]*conj(tpu[i])) + norm(tpu[i])+norm(tpd[i])-norm(opx[i])-norm(opy[i])) -jj*(real(tpu[i]*tpd[i] + conj(tpd[i])*conj(tpu[i])));
			}
			return sum;
		}

		private:

		void resizeParams(int vol)
		{
			nxu.resize(vol);
			nxd.resize(vol);
			nyu.resize(vol);
			nyd.resize(vol);
			tpu.resize(vol);
			tpd.resize(vol);
			spx.resize(vol);
			spy.resize(vol);
			opx.resize(vol);
			opy.resize(vol);
		}

		//ParametersType& engineParams_;
		EngineParamsType& engineParams_;
		ModelParamsType& modelParams_;
		int lx_;
		int ly_;
		int volume_;
		RandomNumberGeneratorType rng_;
		FieldType ne;

		FieldType uu;
		FieldType up;
		FieldType jj;
	}; //MFParams

} //HF

#endif
