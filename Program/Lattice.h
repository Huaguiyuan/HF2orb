
/** \ingroup HF */

/*! \file Lattice.h
 *
 * A 2D lattice
 *
 */

#ifndef LATTICE_H
#define LATTICE_H
#include "Utils.h"

namespace HF {
	template<typename ParametersType>
	class Lattice {
		public:

		/*	nearest- and next-nearest-neighbors:

			XNYP		YP		XPYP
					 |
			XN ------------- O -------------- XP
					 |
			XNYN		YN		XPYN
		 */
		enum {XP=0,XN=1,YP=2,YN=3,XPYP=4,XNYN=5,XPYN=6,XNYP=7};

		//enum {PBC=0,OBC=1,POBC=2};
		//psimag::Matrix<int> neigh;
//		Bc bc_;
		//typedef std::pair<size_t,size_t> PairType;
		Lattice(ParametersType& engineParams) 
			: engineParams_(engineParams),lx_(engineParams_.lx), ly_(engineParams_.ly), volume_(lx_*ly_)
			//, boundaryConditions_(engineParams.boundaryConditions)
		{
			//bc_=BoundaryConditions(boundaryConditions_);
			//if(bc_==OBC||bc_==POBC) throw std::runtime_error("Sorry, cannot support OBC or POBC currently!!!");
			psimag::Matrix<int> neighbors(volume_, 8);
			buildNeighbors(neighbors);
			neigh.resize(lx_,ly_);
			neigh=neighbors;
		}
				
/*
		Bc BoundaryConditions(const std::string& boundaryConditions) 
		{
			if(boundaryConditions=="periodic" || boundaryConditions=="Periodic" || boundaryConditions=="PERIODIC") 
			return PBC;
			else if(boundaryConditions=="open" || boundaryConditions=="Open" || boundaryConditions=="OPEN") 
			return OBC;
			else if(boundaryConditions="semi" || boundaryConditions="Semi" || boundaryConditions="SEMI") 
			return POBC;
			else throw std::runtime_error("Cannot deal with this boundary condition!!!");
		}
*/		
		int dim() const { return 2; }

		int volume() const { return lx_*ly_; }
		
		int lengthX() const { return lx_; }

		int lengthY() const { return ly_; }

		int getNeighbour (int i, int dir)
		{
			return neigh(i,dir);
		}

		int getNeighbour (int x, int y, int dx, int dy) 
		{
			int i=x+y*lx_;
			int dir=getDir(dx, dy);
			return neigh(i,dir);
		}
		
		void printout(){
			int vol=volume();
			std::cout<<"lattice.volume="<<volume()<<std::endl;
			std::cout<<"lattice.neighbors\n";
			for(int i=0; i<vol;i++){
				for(int dir=0; dir<8;dir++)
					std::cout<<neigh(i,dir)<<"\t";
				std::cout<<"\n";
			}
		}
		
		int add(int ind,int ind2) const
		{
			std::vector<int> x(2),y(2);
			index2Coor(x,ind);
			index2Coor(y,ind2);
			for (size_t i=0;i<x.size();i++) {
				x[i] += y[i];
				//g_pbc(x[i],lx_);
			}
			g_pbc(x[0],lx_);
			g_pbc(x[1],ly_);
			return g_index(x);
		}

		void index2Coor(std::vector<int> &v,int i) const
		{
			int lx = lx_;
			v[0] = i%lx;
			v[1] = int(i/lx);
		}

		private:

		void buildNeighbors(psimag::Matrix<int>& neighbors)
		{
			int lx = lx_;
			int ly = ly_;
			int zz = 0;
			int  i=0;
			
			for (int y=0;y<ly;y++) {
				for (int x=0;x<lx;x++) {
					i = x + y*lx;
					int counter = 0;
					//---------------nearest-neighbors------------------
					int xx=x+1;
					int yy=y;
					neighbors(i,counter++) = g_index(xx,yy,zz);

					xx=x-1;					
					neighbors(i,counter++) = g_index(xx,yy,zz);

					xx=x; 
					yy=y+1;					
					neighbors(i,counter++) = g_index(xx,yy,zz);

					yy=y-1;					
					neighbors(i,counter++) = g_index(xx,yy,zz);

					//--------------next-nearest-neighbors-------------
					xx=x+1; 
					yy=y+1;
					neighbors(i,counter++) = g_index(xx,yy,zz);

					xx=x-1; yy=y-1;
					neighbors(i,counter++) = g_index(xx,yy,zz);

					xx=x+1; yy=y-1;
					neighbors(i,counter++) = g_index(xx,yy,zz);

					xx=x-1; yy=y+1;
					neighbors(i,counter++) = g_index(xx,yy,zz);

					if (counter != 8) throw std::runtime_error("ERROR: neighbors -> counter");
				}
			}

			if (i != (lx*ly-1)) throw std::runtime_error("ERROR: neighbors -> i");
			
		}
		
		bool g_pbc(int& x, int l) const
		{
			int L = l;
			bool r=false;
			if (x<0) r=true; 
			if (x>=L) r=true; 
			while(x<0) x+=L;
			while(x>=L) x-=L;
			return r;
		}
		
		int g_index(std::vector<int>& x) const
		{
			int zz=0;
			return g_index(x[0],x[1],zz);
		}
		
		int g_index(int& x,int& y,int& z) const
		{
			int lx = lx_;
			int ly = ly_;
			g_pbc(x,lx);
			g_pbc(y,ly);
			//g_pbc(z,lz);
			return x+y*lx; //+z*L*L;
		}
		
		int getDir(int dx, int dy) const
		{
			if (dx==1)
				if(dy==1) return XPYP;
				else if (dy==0) return XP;
				else if (dy==-1) return XPYN;
				else throw std::runtime_error("ERROR: neighbors -> dy"); 
			else if (dx==0)
				if(dy==1) return YP;
				else if (dy==-1) return XPYN;
				else throw std::runtime_error("ERROR: neighbors -> dy");
			else if (dx==-1)
				if(dy==1) return XNYP;
				else if (dy==0) return XN;
				else if (dy==-1) return XNYN;
				else throw std::runtime_error("ERROR: neighbors -> dy");
			else throw std::runtime_error("ERROR: neighbors -> dx");
		}

		ParametersType& engineParams_;
		int lx_;
		int ly_;
		int volume_;
		
		/*	the neighbors are stored in Matrix neighbors as:
				
			XP=0,XN=1,YP=2,YN=3,XPYP=4,XNYN=5,XPYN=6,XNYP=7	

		*/
		//const int numDir=8;	
		psimag::Matrix<int> neigh;
	}; //Lattice
	
} // namespace Hf

#endif
