!************************************************************************
      subroutine coordinate(Nc,indx,indy,near,nextnear,nsize,ns)
!************************************************************************(B start)
!use table_site
      implicit real*8 (a-h,o-z)
      integer::indx(ns),indy(ns)
      integer,dimension(ns,4):: near
      integer,dimension(ns,4):: nextnear
      integer:: Nc(nsize,nsize)
      integer:: nsize

      do is = 1, ns
        do j = 1, 4
          near(is,j) = 0
          nextnear(is,j)=0
        end do
      end do



      icount = 0
      do iy = 1, nsize
      do ix = 1, nsize
          icount = icount + 1
          Nc(ix,iy) = icount

          indx(icount) = ix
          indy(icount) = iy
      end do
      end do


      do is = 1, ns
        ix = indx(is)
        iy = indy(is)

        ix1 = ix + 1 + nsize
        iy1 = iy
        ix1 = mod(ix1 - 1,nsize) + 1

        ix2 = ix
        iy2 = iy + 1 + nsize
        iy2 = mod(iy2 - 1,nsize) + 1

        ix3 = ix - 1 + nsize
        iy3 = iy
        ix3 = mod(ix3 - 1,nsize) + 1

        ix4 = ix
        iy4 = iy - 1 + nsize
        iy4 = mod(iy4 - 1,nsize) + 1

        near(is,1) = Nc(ix1,iy1)
        near(is,2) = Nc(ix2,iy2)
        near(is,3) = Nc(ix3,iy3)
        near(is,4) = Nc(ix4,iy4)

!define next near neibour 4 site coodinate
        ix5 = ix +1 + nsize              !!(x=1,y=1)
        iy5 = iy +1 + nsize
        ix5 = mod(ix5 - 1,nsize) + 1
        iy5 = mod(iy5 - 1,nsize) + 1

        ix6 = ix - 1 + nsize             !!(x=-1,y=1)
        iy6 = iy + 1 + nsize
        ix6 = mod(ix6 - 1,nsize) + 1
        iy6 = mod(iy6 - 1,nsize) + 1

        ix7 = ix - 1 + nsize             !!(x=-1,y=-1)
        iy7 = iy - 1 + nsize
        ix7 = mod(ix7 - 1,nsize) + 1
        iy7 = mod(iy7 - 1,nsize) + 1

        ix8 = ix + 1 + nsize             !!(x=1,y=-1)
        iy8 = iy - 1 + nsize
        ix8 = mod(ix8 - 1,nsize) + 1
        iy8 = mod(iy8 - 1,nsize) + 1

        nextnear(is,1) = Nc(ix5,iy5)
!        nextnear(is,2) = Nc(ix6,iy6)
        nextnear(is,2) = Nc(ix8,iy8)
        nextnear(is,3) = Nc(ix7,iy7)
!        nextnear(is,4) = Nc(ix8,iy8)
       nextnear(is,4) = Nc(ix6,iy6)


      end do



       end subroutine coordinate
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(B end)

