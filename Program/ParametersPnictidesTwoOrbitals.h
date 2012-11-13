/** \ingroup HF */

/*! \file ParametersPnictidesTwoOrbitals.h
 *
 *  Contains the parameters for the Hamiltonian
 *
 */
#ifndef PARAMETERSPNICTIDESTWOORBITALS_H
#define PARAMETERSPNICTIDESTWOORBITALS_H
#include "Utils.h"
#include "SimpleReader.h"

namespace HF {
	//! Hubbard Model Parameters
	template<typename Field>
	struct ParametersPnictidesTwoOrbitals {

		std::vector<Field> hoppings; 
		Field U;
		Field Up;
		Field J;

	};

	//! Operator to read Model Parameters from inp file.
	template<typename FieldType>
	ParametersPnictidesTwoOrbitals<FieldType>&
	operator <= (ParametersPnictidesTwoOrbitals<FieldType>& parameters,  Dmrg::SimpleReader& reader) 
	{
		reader.read(parameters.hoppings);
		reader.read(parameters.U);
		//reader.read(parameters.Up);		
		reader.read(parameters.J);
		parameters.Up = parameters.U - 2.0 * parameters.J;
		return parameters;
	}
	
	//! Function that prints model parameters to stream os
	template<typename FieldType>
	std::ostream& operator<<(std::ostream &os,const ParametersPnictidesTwoOrbitals<FieldType>& parameters)
	{
		os<<"ModelParam.U="<<parameters.U<<"\n";
		os<<"ModelParam.Up="<<parameters.Up<<"\n";
		os<<"ModelParam.J="<<parameters.J<<"\n";
		os<<"ModelParam.hoppings\n";
		os<<parameters.hoppings;
		return os;
	}
} // namespace HF

#endif
