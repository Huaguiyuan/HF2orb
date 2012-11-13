!**********************************************************************
     subroutine Probability(P, muu, eigen, E, beta, ndim, E_total)
!**********************************************************************(A start)
!use table_site
    implicit none

      integer::i,k
      integer::ndim     
      real(8):: P
      real(8):: beta
      real(8):: eigen(ndim)
      real(8):: E, E_total   ! E: clasical spin,     E_total: total energy

      real(8):: fermi
      real(8):: muu
      real(8):: x

     x=0.0
     P=0.0
     do i=1, ndim
!    write(*,*) i,P
     x = Beta* ( muu - (eigen(i)) )

     if ( x > real(5,8)) then
         P= P + x
!         write(*,*) 'aaa', P
!         aaa=aaa+1
     else if ( x <real(0.001,8) .and. x  >real(-0.001,8)) then
        P= P + DLOG(2.d0 + x )
!         write(*,*) 'bbb',P
!          bbb=bbb+1
     else if ( x <real(-5,8)) then
         P= P + DEXP( x )
!         write(*,*) 'ccc',P
!          ccc=ccc+1
     else 
      
      P = P + DLOG(real(1.0,8) + DEXP( x ))
!         write(*,*) 'ddd' ,P,  'x', x
!          ddd=ddd+1
     endif

     enddo



!      write(*,*) 'P',P, 'beta*E', beta*E
   
      P = P - Beta*E  
       E_total=E 
     do k=1, ndim          ! No. j eigen state
     E_total=E_total+ eigen(k)* fermi( muu , eigen(k),beta) 
     enddo

       end subroutine Probability
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(A end) 



!**********************************************************************
     subroutine total_energy(E_total, eigen, E_classic, ndim,ns,N_optimize)  ! only for optimization
!**********************************************************************(A start)
!use table_site
    implicit none

      integer::i,a
      integer::ndim ,ns    
      real(8):: E_total
      real(8):: E_classic
      real(8):: N_optimize
      real(8):: eigen(ndim)


     E_total=E_classic
     a=int(N_optimize*real(ns,8)+0.00001)
     do i=1, a
     E_total=E_total+eigen(i)
     enddo

       end subroutine total_energy
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(A end) 

!**********************************************************************
     subroutine statistic(data_list, average, error, data_size )
!**********************************************************************(F start)
       implicit none
   
      integer::data_size
      integer::i
      real(8)::data_list(data_size)
      real(8)::average,error
      real(8)::data_sum, data_square
 
        data_sum=0.d0
        data_square=0.d0
        do i=1,data_size
            data_sum=data_sum+data_list(i)
        enddo
        average=data_sum/real(data_size,8)

        do i=1,data_size
            data_square=data_square+(data_list(i)-average)**2
        enddo
        error=(data_square/data_size)**0.5

       return
      end subroutine statistic     


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(F end)     


