/*	HF for pnictides by Q. Luo
 *	Data: Sep 3, 2010
 *	Platform: linux
 *	Model: Mean-Field Approximation of Multiorbital Hubbard Model 
 */

#include "SimpleReader.h"
#include "ParametersEngine.h"
#include "ParametersPnictidesTwoOrbitals.h"
#include "Lattice.h"
#include "MFParams.h"
#include "Ham.h"
#include "RandomNumberGenerator.h"
#include "Engine.h"
#include "GreenFunction.h"

typedef double FieldType;
typedef HF::ParametersEngine<FieldType> ParametersEngineType;
typedef HF::ParametersPnictidesTwoOrbitals<FieldType> ParametersModelType;
typedef HF::Lattice<ParametersEngineType> LatticeType;
typedef Spf::RandomNumberGenerator<FieldType> RandomNumberGeneratorType;
typedef HF::MFParams<ParametersEngineType, ParametersModelType, RandomNumberGeneratorType, FieldType> MFParamsType;
typedef HF::Ham<ParametersEngineType,ParametersModelType,MFParamsType,LatticeType,FieldType> HamiltonianType;
typedef HF::Engine<ParametersEngineType,ParametersModelType,MFParamsType,LatticeType,HamiltonianType, FieldType> EngineType;

int main(int argc,char *argv[])
{
	ParametersEngineType engineParams;
	ParametersModelType mp;
	Dmrg::SimpleReader reader(argv[1]);

	std::cout<<"Loading"<<std::endl;
	reader.load(engineParams);
	reader.load(mp);

	std::cout<<"Initialize Lattice"<<std::endl;
	LatticeType lattice(engineParams);

	std::cout<<"Initialize MF Parameters"<<std::endl;
	MFParamsType mfParams(engineParams,mp);
	
	EngineType engine(engineParams, mp, mfParams, lattice);
	engine.run();

}

/* ####### End of main.cpp ######## */

