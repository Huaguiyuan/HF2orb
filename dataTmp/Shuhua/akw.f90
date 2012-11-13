!**********************************************************************************
        subroutine Akw(ham,kdim,eigen,nsize,ndim,ns,indx,indy,Spec_a)
!**********************************************************************************(G start)
!  use table_site
  use diag_mod
  implicit none
  integer::nsize,ns,xr,yr,ndim,kdim
  integer::k,kx,ky,l
  integer::k_a, k_b, i_a, i_b,n
  integer:: indx(ns),indy(ns)                     !tell you x,y from ns
      real(8)::ei,ej,ai,aj       !four different anger
      complex(8)::Ham(ndim,ndim), E_ik,  complex_akw
      real(8)::eigen(ndim)
      integer::Nc(nsize,nsize)
      integer::i,j,ic,k1,k2
!      real(8)::u(ns,4)
!      real(8)::Q(ns,3)
      real(8)::phase_1,phase_2, q_a, q_b,a,b,c,d

!  real(8),dimension(0:nsize-1,0:nsize-1,1:ndim,0:3):: Spec                    !spin correlation function in k   space 
 real(8),dimension(0:nsize-1,0:nsize-1):: Spec_a
!  real(8),parameter::pi=3.14159265358979


 
!  write(*,*) 'JJJ',JJJ

 

!  Spec=0.0
  Spec_a=0.0
  d=0

 
!________________________________________________________
!do i=1,nsize
!do j=1,nsize
!write(*,*) i,j,sx(i,j),sy(i,j)
!enddo
!enddo

n=4*ndim/6
! do l=1,ndim
do l=n-1,n-1

   do kx=0,nsize-1
   do ky=0,nsize-1


      do k=0,2   !alfa
      complex_akw=cmplx(0.0,0.0)
      do i = 1, ns   !i
      do j = 1, ns   !j
           xr=indx(j)-indx(i)
           yr=indy(j)-indy(i)
           a=2.0*pi*real(kx*xr+ky*yr,8)/real(nsize,8)
           E_ik=cdexp(cmplx(0.0, a, 8) )
!write(*,*) E_ik
           complex_akw= complex_akw+E_ik* (ham(i+ns*k ,l) * conjg(ham(j+ns*k ,l)))    ! l: # of eigen state, k # orbital, i,j: site 
           complex_akw= complex_akw+E_ik* (ham(i+ns*k+ns*3 ,l) * conjg(ham(j+ns*k+ns*3 ,l)))    ! l: # of eigen state, k # orbital, i,j: site 
      enddo  
      enddo


!      do i = 1, ns   !i
!          complex_akw= complex_akw+(ham(i ,l) * conjg(ham(i ,l)))  
!           complex_akw= complex_akw+(ham(i+ns*k ,l) * conjg(ham(i+ns*k ,l)))    ! l: # of eigen state, k # orbital, i,j: site 
!           complex_akw= complex_akw+(ham(i+ns*k+ns*3 ,l) * conjg(ham(i+ns*k+ns*3 ,l)))    ! l: # of eigen state, k # orbital, i,j: site 
!      enddo  

!      Spec(kx,ky,l,k+1)=real(complex_akw,8)/real(ns,8)
!      Spec(kx,ky,l,0)=Spec(kx,ky,l,0)+Spec(kx,ky,l,k+1)
!      d=d+Spec(kx,ky,l,0)
      Spec_a(kx,ky)=Spec_a(kx,ky)+real(complex_akw,8)/real(ndim,8)
      d=d+real(complex_akw,8)/real(ns,8)

!_________________________________________________________
      enddo  ! do k=0,2
 enddo
 enddo

 enddo ! do l=1,ndim



   do kx=0,nsize-1
   do ky=0,nsize-1
! write(34,*) kx, ky, Spec(kx,ky,n,0)
 write(34,*) kx, ky, Spec_a(kx,ky)
 enddo
 write(34,*) " "
 enddo 
!________________________________________________________

write(*,*) "d",d

 return
end  subroutine akw
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(G end)














!**********************************************************************************
        subroutine order(ham,eigen,nsize,ndim,ns,indx,indy,n)
!**********************************************************************************(G start)
!  use table_site
  use diag_mod
  implicit none
  integer::nsize,ns,xr,yr,ndim
  integer::k,kx,ky,l
  integer::k_a, k_b, i_a, i_b,n
  integer:: indx(ns),indy(ns)                     !tell you x,y from ns
      real(8)::ei,ej,ai,aj       !four different anger
      complex(8)::Ham(ndim,ndim), E_ik,  complex_n
      real(8)::eigen(ndim)
      integer::Nc(nsize,nsize)
      integer::i,j,ic,k1,k2
!      real(8)::u(ns,4)
!      real(8)::Q(ns,3)
      real(8)::phase_1,phase_2, q_a, q_b,a,b,c,d,den(0:5),den_a(1:3,1:2)

!  real(8),dimension(0:nsize-1,0:nsize-1,1:ndim,0:3):: Spec                    !spin correlation function in k   space 
 real(8),dimension(0:nsize-1,0:nsize-1):: Spec_a
!  real(8),parameter::pi=3.14159265358979


 
!  write(*,*) 'JJJ',JJJ

 

!  Spec=0.0
  Spec_a=0.0
  d=0

 
!________________________________________________________
!do i=1,nsize
!do j=1,nsize
!write(*,*) i,j,sx(i,j),sy(i,j)
!enddo
!enddo


! do l=1,ndim

      do i = 1, ns   !i
      den=0.0
      den_a=0.0
      do k=0,2   !alfa
      complex_n=cmplx(0.0,0.0)
      do l=1,n

           complex_n= complex_n+ (ham(i+ns*k ,l) * conjg(ham(i+ns*k ,l)))    ! l: # of eigen state, k # orbital, i,j: site
           den(4)= den(4) +(ham(i+ns*k ,l) * conjg(ham(i+ns*k ,l)))  
           den_a(k+1,1)= den_a(k+1,1)+(ham(i+ns*k ,l) * conjg(ham(i+ns*k ,l))) 
           complex_n= complex_n+ (ham(i+ns*k+ns*3 ,l) * conjg(ham(i+ns*k+ns*3 ,l)))    ! l: # of eigen state, k # orbital, i,j: site 
            den(5)= den(5) +(ham(i+ns*k+ns*3 ,l) * conjg(ham(i+ns*k+ns*3 ,l))) 
           den_a(k+1,2)= den_a(k+1,1)+(ham(i+ns*k+ns*3 ,l) * conjg(ham(i+ns*k+ns*3 ,l))) 
      enddo ! do l=1,ndim
       den(k+1)=real(complex_n,8)
       den(0)= den(0)+den(k+1)
 
      enddo  ! do k=0,2
     write(33,*) indx(i), indy(i), den(1), den(2), den(3),den(0) !den(4),den(5)
 !       write(33,*) i, den_a(1,1), den_a(1,2), den_a(2,1), den_a(2,2), den_a(3,1), den_a(3,2)
 !   d= den_a(1,1)+den_a(1,2)+den_a(2,1)+den_a(2,2)+den_a(3,1)+den_a(3,2)
 !     write(33,*) i, den_a(1,1)+den_a(1,2), den_a(2,1)+den_a(2,2), den_a(3,1)+den_a(3,2),d
      enddo
!      Spec(kx,ky,l,k+1)=real(complex_akw,8)/real(ns,8)
!      Spec(kx,ky,l,0)=Spec(kx,ky,l,0)+Spec(kx,ky,l,k+1)
!      d=d+Spec(kx,ky,l,0)


!_________________________________________________________









 return
end  subroutine order
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(G end)





!**********************************************************************************
        subroutine Akw3(ham,kdim,eigen,nsize,ndim,ns,indx,indy,k_a,k_b,phase_1,phase_2,Spec_b)
!**********************************************************************************(G start)
!  use table_site
  use diag_mod
  implicit none
  integer::nsize,ns,xr,yr,ndim,kdim
  integer::k,l
  integer::k_a, k_b, i_a, i_b,n
  integer:: indx(ns),indy(ns)                     !tell you x,y from ns
      real(8)::ei,ej,ai,aj       !four different anger
      real(8)::xe
      real(8)::ye
      complex(8)::Ham(ndim,ndim), E_ik,  complex_akw
      real(8)::eigen(ndim),kx,ky
      integer::Nc(nsize,nsize)
      integer::i,j,ic,k1,k2
!      real(8)::u(ns,4)
!      real(8)::Q(ns,3)
      real(8)::phase_1,phase_2, q_a, q_b,a,b,c,d

!  real(8),dimension(0:nsize-1,0:nsize-1,1:ndim,0:3):: Spec                    !spin correlation function in k   space 
 real(8):: Spec_b(1:ndim,0:1)
!  real(8),parameter::pi=3.14159265358979


      kx=(2.0*pi*k_a-phase_1)/real(nsize,8)
      ky=(2.0*pi*k_b-phase_2)/real(nsize,8)
 
!  write(*,*) 'JJJ',JJJ

!        xe=1.0
!        ye=1.0
!        if(kx**2+ky**2/=0) then
!         xe=real(kx,8)/(real(kx**2+ky**2))**0.5
!         ye=real(ky,8)/(real(kx**2+ky**2))**0.5
!        endif

 

!  Spec=0.0
  Spec_b=0.0
  d=0

 
!________________________________________________________
!do i=1,nsize
!do j=1,nsize
!write(*,*) i,j,sx(i,j),sy(i,j)
!enddo
!enddo


 do l=1,ndim
!do l=n-1,n-1



      do k=0,1   !alfa
      complex_akw=cmplx(0.0,0.0)
      do i = 1, ns   !i
      do j = 1, ns   !j
           xr=indx(j)-indx(i)
           yr=indy(j)-indy(i)
           a=kx*real(xr,8)+ky*real(yr,8)
           E_ik=cdexp(cmplx(0.0, a, 8) )   

           complex_akw= complex_akw+E_ik* (ham(i+ns*k ,l) * conjg(ham(j+ns*k ,l)))    ! l: # of eigen state, k # orbital, i,j: site 
           complex_akw= complex_akw+E_ik* (ham(i+ns*k+ns*2 ,l) * conjg(ham(j+ns*k+ns*2 ,l)))    ! l: # of eigen state, k # orbital, i,j: site 
      enddo  
      enddo


      Spec_b(l,k)=Spec_b(l,k)+real(complex_akw,8)/real(ns,8)
      d=d+real(complex_akw,8)/real(ns,8)

!_________________________________________________________
      enddo  ! do k=0,1


 enddo ! do l=1,ndim




!________________________________________________________

!write(*,*) "d",d

 return
end  subroutine akw3
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(G end)















!**********************************************************************************
        subroutine Akw2(ham,kdim,eigen,nsize,ndim,ns,indx,indy,Spec)
!**********************************************************************************(G start)
!  use table_site
  use diag_mod
  implicit none
  integer::nsize,ns,xr,yr,ndim,kdim
  integer::k,kx,ky,l
  integer::k_a, k_b, i_a, i_b
  integer:: indx(ns),indy(ns)                     !tell you x,y from ns
      real(8)::ei,ej,ai,aj       !four different anger
      complex(8)::Ham(ndim,ndim), E_ik,  complex_akw
      real(8)::eigen(ndim)
      integer::Nc(nsize,nsize)
      integer::i,j,ic,k1,k2
!      real(8)::u(ns,4)
!      real(8)::Q(ns,3)
      real(8)::phase_1,phase_2, q_a, q_b,a,b,c,d

  real(8),dimension(0:nsize-1,0:nsize-1,1:ndim,0:3):: Spec                    !spin correlation function in k   space 

!  real(8),parameter::pi=3.14159265358979


 
!  write(*,*) 'JJJ',JJJ

 

  Spec=0.0


 
!________________________________________________________
!do i=1,nsize
!do j=1,nsize
!write(*,*) i,j,sx(i,j),sy(i,j)
!enddo
!enddo


 do l=1,ndim
      do k=0,2    !alfa
   do kx=0,nsize-1
   do ky=0,nsize-1
      complex_akw=cmplx(0.0,0.0)



      do i = 1, ns   !i
      do j = 1, ns   !j
           xr=indx(j)-indx(i)
           yr=indy(j)-indy(i)
           a=2.0*pi*real(kx*xr+ky*yr,8)/real(ns,8)
           E_ik=cdexp(cmplx(0.0, a, 8) )

           complex_akw= complex_akw+E_ik* (ham(i+ns*k ,l) * conjg(ham(j+ns*k ,l)))    ! l: # of eigen state, k # orbital, i,j: site 
           complex_akw= complex_akw+E_ik* (ham(i+ns*k+ns*3 ,l) * conjg(ham(j+ns*k+ns*3 ,l)))    ! l: # of eigen state, k # orbital, i,j: site 
      enddo  
      enddo
      Spec(kx,ky,l,k+1)=real(complex_akw,8)/real(ns,8)
      Spec(kx,ky,l,0)=Spec(kx,ky,l,0)+Spec(kx,ky,l,k+1)
!_________________________________________________________

 enddo
 enddo
      enddo  ! do k=0,2
 enddo ! do l=1,ndim
! do xr=0,nsize-1
! do yr=0,nsize-1
!   write(*,*) xr,yr,A_str(xr,yr)
! enddo
! enddo 
!________________________________________________________



 return
end  subroutine akw2
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(G end)


