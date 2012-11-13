!*********************************************************
!*********************************************************
!**** 2009/4/28 f90  2-eg t2g model                    ***                                                                 !
!**** eg use zheev dianize it, t29 use classical spin  ***
!**** total energy use MC simulation                   ***
!****      -------------------Shuhua Liang             ***
!****                                                  ***
!****                                                  ***
!*********************************************************
!*********************************************************




module table_site
  real(8),parameter::edown=0.25!0.38 !-2.405
  real(8),parameter::eup=0.3 !0.4 ! 2.0, 
!  real(8)::muu
   real(8)::muu!=0.272415!0.273895!0.392!155369570099635   !392
  integer::n_optimize!=256,N_optimize=256
  integer::kdim
  real(8)::J_Hund!=0.15 !3 !3 !0.25!0.28 !18!8
  real(8)::DeltaXY!=0.4 !0.4 !0.4
  integer::k_phase!=8
  integer::fix_chemical!=1  !$$$$$$$$$$$$$$$$$$$$$$ =0, don't fix;   =1, fix chemical potential!!!!
  integer::nsize!=8  ! x and y dimetion 
  real(8)::temperature!=0.01 
  real(8)::J_x!=0.015
  real(8)::J_y!=0.015
  real(8)::J_next!=0.01  !next nearest J
  integer::MCstable!=3!100
  integer::MCstep!=3!10000
  integer::MC_unit!=2!10000
  integer,parameter::want_2D_spin=1  !$$$   must check!  $$$ =0, 2D   ;  =1,  3D
  integer,parameter::fix_spins=1     !$$$$$$$$$$$$$$$$$$$$$$ =0, don't fix;   =1, fix spins!!!!
  real(8)::mu,mu1,mu2
  real(8):: input_mu,output_mu,twist_mu
  real(8)::N_optimize_changable,alfa


  real(8),parameter::spinJ=1.0

   character(len=1000):: filename,u, perl_length

  integer::measure       ! from 0, 1--ns
  integer::count_MCstable  !from 1---MCstable
  integer::step
  integer::stable
  integer::unstable
  integer::accept, reject
  integer::ns                    !total site
  integer::ndim                  !dimension of Hamiltonion (in this program ndim=ns)
  integer::Irand                 !random number seed
  integer::kx,ky,kz

  real(8)::MC_alpha  ,MC_beta            !acceptance 
  real(8)::beta                  !temperature
  real(8)::aveden(0:2)                !average density
  real(8)::efi1,etheta1
  real(8)::E,Etry
  real(8)::E_total,E_classic,E_total_try
  real(8)::SS
  real(8)::Pold,Ptry
  real(8)::tau_x,tau_z
  real(8)::ss_div,den_div
  real(8)::ss_ave,den_ave
  real(8)::den_i,den_d,M_z
  real(8)::haha1,haha2,haha3
  integer::a,b,c,i_1,i_2,MC_counting
  integer::aaa,bbb,ccc,ddd
  integer::i,j,k,l,lll,l1,l2,l3,l4,l5,k_a,k_b,k1,k2,i_3
  real(8)::x,y,z                   !random number  
  real(8)::lll_real
  real(8)::phase_1,phase_2,sudo_phase, x0,x1,x2,x3,x4,x5,x6,x7,x8,x9, a_i, a_j,conduct_x,conduct_y,ppp,qqq
 complex(8)::green1,green2,com_a  ,com_b ,com_c ,com_d ,com_e ,com_f,com_x ,com_y  

integer::IER                                  !used in zheev subroutine    *


integer,allocatable:: indx(:),indy(:)                     !tell you x,y from ns
integer,allocatable:: near(:,:),nextnear(:,:)             !neigbour position
integer,allocatable:: Nc(:,:)                             !tell you ns from x,y
real(8),allocatable,dimension(:)::ss_list,den_xy_list,list
real(8),allocatable:: efi(:,:)             !back ground big t2g efi
real(8),allocatable:: etheta(:,:)         !back ground big t2g etheta
real(8),allocatable:: den(:)
real(8),allocatable:: eigen(:), eigen1(:), twist_eigen(:),spec_b(:,:) 
real(8),allocatable:: S_str(:,:) ,ave_s_str(:,:) ,con_list(:,:),den_list(:,:)
complex(8),allocatable:: Ham(:,:),Ham1(:,:),eigenstates(:,:),TTx(:,:), TTy(:,:)


end module table_site
