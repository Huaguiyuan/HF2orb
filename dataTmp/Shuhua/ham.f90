




!**************************************************************************
      subroutine Hamiltonian2D(J_Hund,DeltaXY,ham,v_x,v_y,nsize,ndim,ns,near,nextnear,indx,indy,efi,etheta,phase_1,phase_2)
!**************************************************************************(E start)
      implicit none!   real*8 (a-h,o-z)

     integer::nsize,ndim,ns,phon_on_or_not
      real(8)::J_Hund
      real(8)::t_0,t1,t2,t3,t4,t5,t6,t7,t8,DeltaXY,pp
      complex(8)::t(1:ns,0:1,0:1,1:4),ttt(1:ns,0:1,0:1,1:ns)
      real(8)::ei,ej,ai,aj       !four different anger
      complex(8)::Ham(ndim,ndim),v_x(ndim,ndim),v_y(ndim,ndim)
      real(8)::etheta(nsize,nsize)
      real(8):: efi(nsize,nsize)
      integer::indx(ns),indy(ns)
      integer::near(ns,4),nextnear(ns,4)
      integer::i,j,k,ic,k1,k2,l
!      real(8)::u(ns,4)
!      real(8)::Q(ns,3)
      real(8)::phase_1,phase_2, l_i, aaa

!***: t(i,orbital_a, orbital_b, orbital_c, direction)  orbital a: xz   orbital b: yz  orbital c: xz
!***______________________________________

      t1=0.02               ! maria
      t2=0.06               ! maria
      t3=0.03               ! maria
      t4=-0.01               ! maria


      t1=1.3               ! maria
      t2=-1.0               ! maria
      t3=0.85               ! maria
      t4=0.85               ! maria

      t1=1.3               ! maria
      t2=-1.0               ! maria
      t3=-0.85               ! maria
      t4=0.85               ! maria

      t1=-0.13               ! maria
      t2=0.1               ! maria
      t3=0.085               ! maria
      t4=0.085               ! maria
   t=0.0
   ttt=0.0
   v_x=0.0
   v_y=0.0
   ham=0.0
!0x
!1y
!2+x+y
!3-x+y


do l=1,ns
   l_i=real((-1)**(indx(l)+indy(l)),8)

!   write(*,*) l_i,indx(l),indy(l),ns
   t(l,0,0,1)=-t2
   t(l,1,1,2)=-t2

   t(l,0,0,2)=-t1
   t(l,1,1,1)=-t1



   t(l,0,0,3)=-t3  
   t(l,0,0,4)=-t3
   t(l,1,1,3)=-t3
   t(l,1,1,4)=-t3

   t(l,0,1,3)= t4  
   t(l,1,0,4)=-t4
   t(l,0,1,4)=-t4
   t(l,1,0,3)= t4



enddo ! do l ns







!  t(l,:,:,:)=t_0
!***______________________________________

      do i = 1,ndim
      do j = 1,ndim
        Ham(i,j) = 0.0
      end do
      end do


!write(*,*) 'aaa' , ns,ndim

!___________________________________________________________(E1 start define T term)
      do i = 1, ns
      do ic = 1,4
     if (ic<=2)      j = near(i,ic)
     if (ic>=3)      j = nextnear(i,ic-2)



           do k1=0,1
           do k2=0,1

           ham(i+k1*ns,j+k2*ns) = -t(i,k1,k2,ic)

       if (indx(j)-indx(i)>2)  ham(i+k1*ns,j+k2*ns) =   ham(i+k1*ns,j+k2*ns)*cdexp(cmplx(0.0,  phase_1, 8) )
       if (indy(j)-indy(i)>2)  ham(i+k1*ns,j+k2*ns) =   ham(i+k1*ns,j+k2*ns)*cdexp(cmplx(0.0,  phase_2, 8) )

       if (indx(j)-indx(i)<-2)  ham(i+k1*ns,j+k2*ns) =   ham(i+k1*ns,j+k2*ns)*cdexp(cmplx(0.0,  -phase_1, 8) )
       if (indy(j)-indy(i)<-2)  ham(i+k1*ns,j+k2*ns) =   ham(i+k1*ns,j+k2*ns)*cdexp(cmplx(0.0,  -phase_2, 8) ) 
       

           ham(j+k2*ns,i+k1*ns) = conjg(ham(i+k1*ns,j+k2*ns))
!           ttt(i,k1,k2,j)=t(i,k1,k2,ic)
!           ttt(j,k2,k1,i)=t(i,k1,k2,ic)
           enddo
           enddo

      end do
      end do

!write(*,'(48f4.1)') real(ham(1,1:32),8)

      do i = 1, 2*ns
      do j = 1, 2*ns   
           ham(i+2*ns,j+2*ns) = ham(i,j)
      enddo
      enddo   


      do i = 1, ns   !i
           ei=etheta(indx(i),indy(i))
           ai=efi(indx(i),indy(i))
      do k=0,2    !alfa

           ham(i+k*ns,i+k*ns) =  ham(i+k*ns,i+k*ns)- J_Hund*( cos(ei))
           ham(i+k*ns+2*ns,i+k*ns+2*ns) =  ham(i+k*ns+2*ns,i+k*ns+2*ns)- J_Hund*(-cos(ei))

           ham(i+k*ns,i+k*ns+2*ns) =  ham(i+k*ns,i+k*ns+2*ns)- J_hund* sin(ei)*cmplx(cos(ai),-sin(ai))  !S-
           ham(i+k*ns+2*ns,i+k*ns) =  ham(i+k*ns+2*ns,i+k*ns)- J_hund* sin(ei)*cmplx(cos(ai),+sin(ai))  !S+
      enddo  
      enddo



  v_x=0.0
  v_y=0.0
aaa=0
v_x=0



do l=1,ns
   l_i=real((-1)**(indx(l)+indy(l)),8)
!   l_i=0.0

!   t5=0
!   write(*,*) l_i,indx(l),indy(l),ns
   t(l,0,0,1)=-t2
   t(l,1,1,2)=-t2

   t(l,0,0,2)=-t1
   t(l,1,1,1)=-t1

   t(l,0,0,3)=-t3  
   t(l,0,0,4)=-t3
   t(l,1,1,3)=-t3
   t(l,1,1,4)=-t3

   t(l,0,1,3)= t4  
   t(l,1,0,4)=-t4
   t(l,0,1,4)=-t4
   t(l,1,0,3)= t4



enddo ! do l ns

            do i = 1, ns
            do ic = 1,4
               if (ic<=2)      j = near(i,ic)
               if (ic>=3)      j = nextnear(i,ic-2)

            do k1=0,1
            do k2=0,1

           ttt(i,k1,k2,j)=t(i,k1,k2,ic)



       if (indx(j)-indx(i)>2)  ttt(i,k1,k2,j) =   ttt(i,k1,k2,j)*cdexp(cmplx(0.0,  phase_1, 8) )
       if (indy(j)-indy(i)>2)  ttt(i,k1,k2,j) =   ttt(i,k1,k2,j)*cdexp(cmplx(0.0,  phase_2, 8) )

       if (indx(j)-indx(i)<-2)  ttt(i,k1,k2,j) =  ttt(i,k1,k2,j)*cdexp(cmplx(0.0,  -phase_1, 8) )
       if (indy(j)-indy(i)<-2)  ttt(i,k1,k2,j) =  ttt(i,k1,k2,j)*cdexp(cmplx(0.0,  -phase_2, 8) ) 


           ttt(j,k2,k1,i)=conjg(ttt(i,k1,k2,j))
           enddo
           enddo

           enddo
           enddo


      do i = 1,ns
      
    !  if(indx(i)==1) then
           do k1=0,1
           do k2=0,1
           j=near(i,1)                       !x
           v_x(i+k1*ns,j+k2*ns) = -ttt(i,k1,k2,j)
           v_x(j+k2*ns,i+k1*ns) =  ttt(j,k2,k1,i)
           j=nextnear(i,1)                   !+x +y
           v_x(i+k1*ns,j+k2*ns) = -ttt(i,k1,k2,j)
           v_x(j+k2*ns,i+k1*ns) =  ttt(j,k2,k1,i)
           j=nextnear(i,2)                   !+x -y
           v_x(i+k1*ns,j+k2*ns) = -ttt(i,k1,k2,j)
           v_x(j+k2*ns,i+k1*ns) =  ttt(j,k2,k1,i)
           enddo
           enddo
   !    endif

 !   if(indy(i)==1) then
           do k1=0,1
           do k2=0,1
           j=near(i,2)                       !y
           v_y(i+k1*ns,j+k2*ns) = -ttt(i,k1,k2,j)
           v_y(j+k2*ns,i+k1*ns) =  ttt(j,k2,k1,i)
           j=nextnear(i,1)                   !+x +y
           v_y(i+k1*ns,j+k2*ns) = -ttt(i,k1,k2,j)
           v_y(j+k2*ns,i+k1*ns) =  ttt(j,k2,k1,i)
           j=nextnear(i,4)                    !-x -y
           v_y(i+k1*ns,j+k2*ns) = -ttt(i,k1,k2,j)
           v_y(j+k2*ns,i+k1*ns) =  ttt(j,k2,k1,i)
           enddo
           enddo
  !   endif


      end do
!write(*,*) aaa
!write(*,'(48f4.1)') real(ham(1,1:32),8)

      do i = 1, 2*ns
      do j = 1, 2*ns   
           v_x(i+2*ns,j+2*ns) = v_x(i,j)
           v_y(i+2*ns,j+2*ns) = v_y(i,j)
      enddo
      enddo 


!____________________________________________________________(E 1 end)



!    write(*,*) 'bbb'
 
     end subroutine Hamiltonian2D
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(E end)







!*************************************************************************
         subroutine HHH(efi,etheta,nsize,ns,J_x,J_y,J_next,Nc,near,nextnear,H)     
!H=J*S.S                ! wihtout phon,     phon_on_or_not=0
!*************************************************************************(C start)

 implicit none
 integer::nsize,i,j,k,ns
 integer::Nc(nsize,nsize),near(ns,4),nextnear(ns,4)
 real(8)::H
 real(8),intent(in)::J_x,J_y,J_next
 real(8),dimension(nsize,nsize)::efi,etheta
 real(8),dimension(ns)::sx,sy,sz


!  write(*,*) 'JJJ',JJJ

 
 H=0.0
 sx=0.0
 sy=0.0
 sz=0.0

 do i=1,nsize
 do j=1,nsize

    sx(Nc(i,j))=DCOS(efi(i,j))*DSIN(etheta(i,j))
    sy(Nc(i,j))=DSIN(efi(i,j))*DSIN(etheta(i,j))
    sz(Nc(i,j))=DCOS(etheta(i,j)) 
 
 enddo
 enddo


! do i=1,nsize
! do j=1,nsize
!    write(*,'(3f5.2)') sz(Nc(i,j)),sz(near(Nc(i,j),1)),sz(near(Nc(i,j),2))
! enddo
! enddo


 do i=1,ns
    H= H + J_x * (  sx(i) * sx(near(i,1)) + sy(i) * sy(near(i,1)) + sz(i) * sz(near(i,1))  )
    H= H + J_y * (  sx(i) * sx(near(i,2)) + sy(i) * sy(near(i,2)) + sz(i) * sz(near(i,2))  )
 !   write(*,*) sz(i) * sz(near(i,1)),JJJ
 enddo

 do i=1,ns
    H= H + J_next * (  sx(i) * sx(nextnear(i,1)) + sy(i) * sy(nextnear(i,1)) + sz(i) * sz(nextnear(i,1))  )
    H= H + J_next * (  sx(i) * sx(nextnear(i,2)) + sy(i) * sy(nextnear(i,2)) + sz(i) * sz(nextnear(i,2))  )
 !   write(*,*) sz(i) * sz(near(i,1)),JJJ
 enddo



!   write(*,*) 'H',H

 return
 end subroutine 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(C end)











