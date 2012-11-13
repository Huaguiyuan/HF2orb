
/*read model parameters 
 */

#ifndef SIMPLE_READER_H
#define SIMPLE_READER_H
#include "Utils.h"


namespace Dmrg {
	class SimpleReader {
	public:
		SimpleReader(const std::string& fileName): fileName_(fileName),fin_(fileName.c_str())
		{
		}

		~SimpleReader() { fin_.close(); }

		template<typename ParametersType>
		void load(ParametersType& parameters) 
		{
			parameters <= *this; // operator<=(parameters,*this)
		}

		void read(std::string& s)  { readALabel(); fin_>>s; }

		template<typename FieldType>
		void read(std::vector<FieldType>& v) 
		{
			// read label
			readALabel();
			// read vector length
			size_t size;
			fin_>>size;
			// read vector
			v.resize(size);
			for (size_t i=0;i<v.size();i++) {
				fin_>>v[i];
			}
		}

		template<typename FieldType>
		void read(psimag::Matrix<FieldType>& v) 
		{
			readALabel();	
			// read row and col
			size_t row,col;
			fin_>>row;
			fin_>>col;
			// read matrix
			v.resize(row,col);
			for (size_t i=0;i<v.n_row();i++) 
				for (size_t j=0;j<v.n_col();j++) 
					fin_>>v(i,j);
		}

		template<typename FieldType>
		void read(FieldType& x)  { readALabel(); fin_>>x; }

		bool findLabel(const std::string& label)
		{
			std::string tmp;
			while(!fin_.eof()) {
				fin_>>tmp;
				//std::cerr<<tmp<<std::endl;
				if (tmp == label) return true;
			}
			return false;
		}


	private:
		std::string    fileName_,label_;
		std::ifstream fin_;

		void readALabel()
		{
			while(true && !fin_.eof()) {
				fin_>>label_;
				if (label_[0]!='#') break;
			} 
		}

	};
} // namespace Dmrg

/*@}*/
#endif


