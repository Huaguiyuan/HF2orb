
/*! \file ParametersEngine.h
 *
 *  Contains the Engine EngineParam
 *
 */
#ifndef PARAMETERSENGINE_H
#define PARAMETERSENGINE_H

#include "Utils.h"
#include "SimpleReader.h"

namespace HF{
	
	//! Structure that contains the Engine EngineParam
	template<typename FieldType_>
	struct ParametersEngine {
		typedef FieldType_ FieldType;
		
		std::string output; // filename to save observables and continued fractions
		FieldType density;
		int lx; //lattice: x
		int ly; //lattice: y
	  std::string boundaryConditions; // boundary conditions	
		FieldType mu; // chemical potential
		FieldType beta; // inverse temperature
		int iterMAX; //maximum of allowed iterations
		std::vector<FieldType> eps; // eps: eps_MFA, eps_shell, eps_Akw, eps_DOS
		std::string mfparams; //how to initalize MF EngineParam
		long randomSeed; //random seed to initialize MF EngineParam 
		FieldType alpha; //parameter of linear mixing
		FieldType omegaStep; // interval of omega in DOS 
	};

	//! Read EngineParam from input file
	template<typename FieldType>
	ParametersEngine<FieldType>&
	operator <= (ParametersEngine<FieldType>& EngineParam,  Dmrg::SimpleReader& reader) 
	{
		reader.read(EngineParam.output);
		reader.read(EngineParam.density); 
		reader.read(EngineParam.lx);
		reader.read(EngineParam.ly);
		reader.read(EngineParam.boundaryConditions);
		reader.read(EngineParam.mu);
		reader.read(EngineParam.beta);
		reader.read(EngineParam.iterMAX);
		reader.read(EngineParam.eps);
		reader.read(EngineParam.mfparams);
		reader.read(EngineParam.randomSeed);
		reader.read(EngineParam.alpha);
		reader.read(EngineParam.omegaStep);
	
		return EngineParam;
	} 

	//! print EngineParam
	template<typename FieldType>
	std::ostream &operator<<(std::ostream &os,ParametersEngine<FieldType> const &EngineParam)
	{
		os<<"EngineParam.filename="<<EngineParam.output<<"\n";
		os<<"EngineParam.density="<<EngineParam.density<<"\n";
		os<<"EngineParam.mu="<<EngineParam.mu<<"\n";
		os<<"EngineParam.lx/y="<<EngineParam.lx<<"/"<<EngineParam.ly<<"\n";
		os<<"EngineParam.boundaryConditions="<<EngineParam.boundaryConditions<<"\n";
		os<<"EngineParam.beta="<<EngineParam.beta<<"\n";
		os<<"EngineParam.iterMAX="<<EngineParam.iterMAX<<"\n";
		os<<"EngineParam.eps\n"<<EngineParam.eps<<"\n";
		os<<"EngineParam.intialMFparams="<<EngineParam.mfparams<<"\n";
		os<<"EngineParam.seed="<<EngineParam.randomSeed<<"\n";
		os<<"EngineParam.alpha="<<EngineParam.alpha<<"\n";
		os<<"EngineParam.omegaStep="<<EngineParam.omegaStep<<"\n";
		return os;
	}
} // namespace HF

#endif
