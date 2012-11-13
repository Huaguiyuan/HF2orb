

!**********************************************************************************
        subroutine spin_structure(efi,etheta,nsize,ns,indx, indy , Nc,near, S_str)
!**********************************************************************************(G start)

  use diag_mod
  implicit none
  integer::nsize,i,j,ns,xr,yr
  integer::k,kx,ky,ss
  integer::Nc(nsize,nsize),near(ns,4)
  integer:: indx(ns),indy(ns)                     !tell you x,y from ns
  real(8),dimension(nsize,nsize)::efi,etheta,ssx,ssy
  real(8),dimension(nsize,nsize)::sx,sy,sz
  real(8),dimension(0:nsize-1,0:nsize-1):: A_str                    !spin correlation function in x y space
  real(8),dimension(0:nsize-1,0:nsize-1):: S_str                    !spin correlation function in k   space 

!  real(8),parameter::pi=3.14159265358979


 
!  write(*,*) 'JJJ',JJJ

 
 sx=0.0
 sy=0.0
 sz=0.0
 A_str=0.0
 s_str=0.0
!_________________________________________________________
         do i=1,nsize
         do j=1,nsize
          
           sx(i,j)=DCOS(efi(i,j))*DSIN(etheta(i,j))
           sy(i,j)=DSIN(efi(i,j))*DSIN(etheta(i,j))
           sz(i,j)=DCOS(etheta(i,j)) 

         enddo
         enddo
!________________________________________________________


 do xr=0,nsize-1
 do yr=0,nsize-1

         do i=1,nsize
         do j=1,nsize
    
   A_str(xr,yr)=A_str(xr,yr) + sx(i,j) * sx(mod(i+xr+2*nsize-1,nsize)+1,mod(j+yr+2*nsize-1,nsize)+1)
   A_str(xr,yr)=A_str(xr,yr) + sy(i,j) * sy(mod(i+xr+2*nsize-1,nsize)+1,mod(j+yr+2*nsize-1,nsize)+1)
   A_str(xr,yr)=A_str(xr,yr) + sz(i,j) * sz(mod(i+xr+2*nsize-1,nsize)+1,mod(j+yr+2*nsize-1,nsize)+1)
         enddo
         enddo   
   A_str(xr,yr)=A_str(xr,yr)/real(nsize**2,8)
 enddo
 enddo



   S_str=0.0
   do kx=0,nsize-1
   do ky=0,nsize-1
     do xr=0,nsize-1
     do yr=0,nsize-1
   
   S_str(kx,ky)=A_str(xr,yr)*Cos(2.0*pi*real(kx*xr+ky*yr,8)/real(nsize,8))+S_str(kx,ky)

     enddo
     enddo
   enddo
   enddo

 return
end  subroutine spin_structure
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(G end)

