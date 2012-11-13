MODULE diag_mod
  USE nrtype
  IMPLICIT NONE

  INTERFACE STEV95
     MODULE PROCEDURE DSTEV95 
  END INTERFACE

  INTERFACE
     SUBROUTINE DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
       IMPLICIT NONE
       INTEGER, PARAMETER       :: dp = KIND(1.0d0)
       CHARACTER, INTENT(in)    :: JOBZ
       INTEGER,   INTENT(in)    :: N, LDZ
       INTEGER,   INTENT(out)   :: INFO
       REAL(dp),  INTENT(inout) :: D(N), E(N), WORK(*),Z(LDZ,N)  
     END SUBROUTINE DSTEV
                    
     SUBROUTINE DSTEVX( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, &
     &                   M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )
       IMPLICIT NONE
       INTEGER, PARAMETER       :: dp = KIND(1.0d0)
       CHARACTER, INTENT(in)    :: JOBZ, range
       INTEGER,   INTENT(in)    :: N, LDZ
       INTEGER,   INTENT(out)   :: INFO
       REAL(DP),  INTENT(in)    :: VL, VU, abstol
       INTEGER,   INTENT(in)    :: IL, IU
       INTEGER,   INTENT(out)   :: M
       REAL(dp),  INTENT(inout) :: D(N), E(N), WORK(5*N)
       INTEGER, INTENT(inout)   :: iwork(5*N)
       INTEGER, INTENT(out  )   :: ifail(N)
       REAL(dp),  INTENT(out)   :: W(N) ,Z(LDZ,*) 
!!$          Note: the user must ensure that at least max(1,M) columns are
!!$          supplied in the array Z; if RANGE = 'V', the exact value of M
!!$          is not known in advance and an upper bound must be used.
     END SUBROUTINE DSTEVX
                    
     SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO ) 
       IMPLICIT NONE
       INTEGER, PARAMETER       :: dp = KIND(1.0d0)
       CHARACTER, INTENT(in)    :: JOBZ, UPLO
       INTEGER,   INTENT(in)    :: N, LDA, LWORK
       INTEGER,   INTENT(out)   :: INFO
       REAL(dp),  INTENT(inout) :: W(N), A(LDA,N)
       REAL(dp),  INTENT(out)   :: WORK(LWORK)
     END SUBROUTINE DSYEV

     SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, &
&      RWORK, LRWORK, IWORK, LIWORK, INFO ) 
       ! full Hermitian matrix
       IMPLICIT NONE
       INTEGER, PARAMETER       :: dp = KIND(1.0d0)
       CHARACTER, INTENT(in)    :: JOBZ, UPLO
       INTEGER,   INTENT(in)    :: N, LDA, LWORK, lrwork, LIWORK
       INTEGER,   INTENT(out)   :: INFO
       REAL(dp),  INTENT(inout) :: W(N)
       COMPLEX(dp),  INTENT(inout) :: A(LDA,N)
       REAL(dp),  INTENT(out)   :: rWORK(LrWORK)
       COMPLEX(dp),  INTENT(out)   :: WORK(LWORK)
       INTEGER,  INTENT(out)    :: iWORK(LiWORK) 
     END SUBROUTINE ZHEEVD

     SUBROUTINE ZHEEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
          & ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL, INFO )
       IMPLICIT NONE
       INTEGER, PARAMETER       :: dp = KIND(1.0d0)
       CHARACTER, INTENT(in)    :: JOBZ, UPLO, range
       INTEGER,   INTENT(in)    :: N, LDA, LWORK, IL, IU, LDZ, M
       REAL(dp),   INTENT(in)   :: ABSTOL, VL, VU
       INTEGER,   INTENT(out)   :: INFO
       REAL(dp),  INTENT(inout) :: W(N)
       COMPLEX(dp),  INTENT(inout) :: A(LDA,N)
       COMPLEX(dp),  INTENT(out) :: Z(LDZ, MAX(1,M))
       COMPLEX(dp),  INTENT(out)   :: WORK(LWORK)
       REAL(dp),  INTENT(out)   :: rWORK(7*N)
       INTEGER,  INTENT(out)    :: iWORK(5*N), ifail(N)
     END SUBROUTINE ZHEEVX



     SUBROUTINE ZHESV(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, &
          & LWORK, INFO)
       ! solves system of equations
       IMPLICIT NONE
       INTEGER, PARAMETER       :: dp = KIND(1.0d0)
       CHARACTER, INTENT(in)    :: UPLO
       INTEGER,   INTENT(in)    :: N, NRHS, LDA, LDB
       INTEGER,   INTENT(in)    :: LWORK
       INTEGER,   INTENT(out)   :: INFO, IPIV(N)
       COMPLEX(dp),  INTENT(inout) :: A(LDA,N), B(LDB,NRHS)
       COMPLEX(dp),  INTENT(out)   :: WORK(LWORK)
     END SUBROUTINE ZHESV

     SUBROUTINE ZGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       IMPLICIT NONE
       INTEGER, PARAMETER       :: dp = KIND(1.0d0)
       INTEGER,   INTENT(in)    :: N, NRHS, LDA, LDB
       INTEGER,   INTENT(out)   :: INFO, IPIV(N)
       COMPLEX(dp),  INTENT(inout) :: A(LDA,N), B(LDB,NRHS)
     END SUBROUTINE ZGESV


     SUBROUTINE ZHBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
       IMPLICIT NONE
       INTEGER, PARAMETER      :: dp = KIND(1.0d0)
       CHARACTER, INTENT(in)   ::   JOBZ, UPLO
       INTEGER, INTENT(in)     ::   KD, LDAB, LDZ, LIWORK, LRWORK, LWORK, N
       INTEGER, INTENT(out)    ::   INFO
       INTEGER, INTENT(out)    ::   IWORK(liwork)
       REAL(dp), INTENT(out)    :: RWORK(lrwork), W(N)
       COMPLEX(dp), INTENT(inout) ::    AB(LDAB, N)
       COMPLEX(dp), INTENT(out) ::    WORK(lwork), Z( LDZ, N)
     END SUBROUTINE ZHBEVD

  END INTERFACE

CONTAINS

  SUBROUTINE ZHEEVD95( A, W, jobz, uplo, INFO )
    USE nrtype
    ! A: Matrix, Eigenvectors
    ! W: eigenvalues
    IMPLICIT NONE
    INTEGER,   INTENT(out),OPTIONAL  :: INFO
    CHARACTER(len=1), INTENT(in)     :: JOBZ, uplo
    REAL(DP), DIMENSION(:),  INTENT(inout) :: W
    COMPLEX(DP), DIMENSION(:,:),  INTENT(inout) :: A

    INTEGER                  :: N, LDA, info_inner, lwork, lrwork, liwork
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rwork
    COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: work
    INTEGER, DIMENSION(:), ALLOCATABLE :: iwork

    N= SIZE(A,1)
    IF (SIZE(W)/=N) STOP 'ZHEEVD95: dimensions of matrix and array for EV do not match!'
    lda=N
    IF (.NOT. ((jobz == 'N') .OR. (JOBZ  == 'V'))) STOP 'ZHEEVD95: Invalid option, must be "N" or "V"!'

    lwork  =  2*N + N**2
    lrwork = 1 + 5*N + 2*N**2
    liwork = 3 + 5*N
    ALLOCATE(work(lwork),rwork(lrwork),iwork(liwork))
    CALL ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, -1, &
&      RWORK, -1, IWORK, -1, INFO_inner ) 
    lwork = work(1)
    liwork = iwork(1)
    lrwork = rwork(1)
    DEALLOCATE(work,rwork,iwork)
    ALLOCATE(work(lwork),rwork(lrwork),iwork(liwork))

    CALL ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, &
&      RWORK, LRWORK, IWORK, LIWORK, INFO_inner ) 
    DEALLOCATE(work, rwork, iwork)

    IF (PRESENT(info)) info = info_inner
  END SUBROUTINE ZHEEVD95

  SUBROUTINE ZHEEVX95( A, W, Z, jobz, uplo, il, iu, vl, vu, N_found, ifail, INFO )
    USE nrtype
    IMPLICIT NONE
    INTEGER,   INTENT(out),OPTIONAL  :: INFO
    INTEGER,   INTENT(in),OPTIONAL  :: il,iu
    INTEGER,   INTENT(out),OPTIONAL  :: N_found
    REAL(DP),   INTENT(in),OPTIONAL  :: vl,vu
    CHARACTER(len=1), INTENT(in), OPTIONAL :: JOBZ, uplo
    REAL(DP), DIMENSION(:),  INTENT(inout) :: W
    COMPLEX(DP), DIMENSION(:,:),  INTENT(inout) :: A
    COMPLEX(DP), DIMENSION(:,:),  INTENT(out), OPTIONAL :: Z
    INTEGER, DIMENSION(:),  INTENT(out), OPTIONAL :: ifail

    COMPLEX(DP), DIMENSION(:,:), ALLOCATABLE :: Z_inner
    INTEGER                  :: N, LDA, ldz, info_inner, lwork, M
    CHARACTER(len=1) :: JOBZ_inner, uplo_inner
    COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: work
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rwork
    REAL(DP) :: abstol = -1.d0
    INTEGER, DIMENSION(:), ALLOCATABLE :: iwork, ifail_inner

    IF (PRESENT(uplo)) THEN
       uplo_inner = uplo
    ELSE
       uplo_inner = 'U'
    END IF

    IF (PRESENT(jobz)) THEN
       jobz_inner = jobz
    ELSE
       jobz_inner = 'U'
    END IF

    N= SIZE(A,1)
    IF (SIZE(W)/=N) STOP 'ZHEEVX95: dimensions of matrix and array for EV do not match!'
    lda=N
    IF (PRESENT(Z)) THEN
       IF (SIZE(Z,1)/=N) STOP 'ZHEEVX95: dimensions of matrix and array for verctors do not match!'
       IF (.NOT. PRESENT(jobz)) THEN
          jobz_inner = 'V'
       END IF
       IF (ALLOCATED(z_inner)) DEALLOCATE(z_inner) !EA
       ALLOCATE(Z_inner(N,N))
       ldz = N
    ELSE
       jobz_inner = 'N'
       IF (ALLOCATED(z_inner)) DEALLOCATE(z_inner) !EA
       ALLOCATE(Z_inner(N,1))
       ldz = 1
    END IF

    lwork  =  MAX((N + 1)*N,2*N-1)
    IF (ALLOCATED(work)) DEALLOCATE(work) !EA
    IF (ALLOCATED(rwork)) DEALLOCATE(rwork) !EA
    IF (ALLOCATED(iwork)) DEALLOCATE(iwork) !EA
    IF (ALLOCATED(ifail_inner)) DEALLOCATE(ifail_inner) !EA
    ALLOCATE(work(lwork),rwork(7*N),iwork(5*N), ifail_inner(N))

    IF (PRESENT(il)) THEN
       IF (.NOT.PRESENT(iu)) STOP 'ZHEEVX95: if il is set, iu also has to be present!'
       IF (PRESENT(vl))  STOP 'ZHEEVX95: only il OR vl makes sense!'
       IF (PRESENT(vu))  STOP 'ZHEEVX95: only il OR vu makes sense!'
       IF (SIZE(Z_inner,2) < (iu-il+1)) STOP 'ZHEEVX95: number of columns in array for eigenvectors too small!'
       CALL ZHEEVX( JOBZ_INNER, 'I', UPLO_inner, N, A, LDA, 1.d0, 1.d0, IL, IU, &
            & ABSTOL, M, W, Z_INNER, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL_INNER, INFO_inner )
    ELSE IF (PRESENT(vl)) THEN
       IF (.NOT.PRESENT(vu)) STOP 'ZHEEVX95: if vl is set, vu also has to be present!'
       IF (PRESENT(il))  STOP 'ZHEEVX95: only il OR vl makes sense!'
       IF (PRESENT(iu))  STOP 'ZHEEVX95: only iu OR vl makes sense!'
       CALL ZHEEVX( JOBZ_INNER, 'V', UPLO_INNER, N, A, LDA, VL, VU, 1, 1, &
            & ABSTOL, M, W, Z_inner, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL_INNER, INFO_inner )
    ELSE
       IF (PRESENT(iu))  STOP 'ZHEEVX95: if iu is set, il also has to be present!'
       IF (PRESENT(vu)) STOP 'ZHEEVX95: if vu is set, vl also has to be present!'
       CALL ZHEEVX( JOBZ_INNER, 'A', UPLO_INNER, N, A, LDA, 1.d0, 1.d0, 1, 1, &
            & ABSTOL, M, W, Z_inner, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL_INNER, INFO_inner )
    END IF
    DEALLOCATE(work, rwork, iwork)

    IF (PRESENT(info)) info = info_inner
    IF (PRESENT(ifail)) ifail = ifail_inner
    IF (PRESENT(Z) .AND. (jobz_inner == 'V')) Z = Z_inner
    IF (PRESENT(N_found)) N_found = M
  END SUBROUTINE ZHEEVX95


  SUBROUTINE DSYEV95( A, W, jobz, uplo, INFO )
    USE nrtype
    IMPLICIT NONE
    INTEGER,   INTENT(out),OPTIONAL  :: INFO
    CHARACTER(len=1), INTENT(in)     :: JOBZ, uplo
    REAL(DP), DIMENSION(:),  INTENT(inout) :: W
    REAL(DP), DIMENSION(:,:),  INTENT(inout) :: A

    INTEGER                  :: N, LDA, info_inner, lwork
    REAL(DP), DIMENSION(:), ALLOCATABLE :: work


    N= SIZE(A,1)
    IF (SIZE(W)/=N) STOP 'DSYEV95: dimensions of matrix and array for EV do not match!'
    lda=N
    
    IF (.NOT. ((jobz == 'N') .OR. (JOBZ  == 'V'))) STOP 'DSYEV95: Invalid option, must be "N" or "V"!'
    
    lwork = MAX(1,3*N-1)
    ALLOCATE(work(lwork))
    CALL DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK,-1, INFO_inner )
    lwork = work(1)
    DEALLOCATE(work)
    ALLOCATE(work(lwork))
    CALL DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK,lwork, INFO_inner )
    DEALLOCATE(work)
    IF (PRESENT(info)) info = info_inner
  END SUBROUTINE DSYEV95

  SUBROUTINE DSTEV95( D, E, Z, INFO )
    USE nrtype
    IMPLICIT NONE
    INTEGER,   INTENT(out),OPTIONAL  :: INFO
    REAL(DP), DIMENSION(:),  INTENT(inout) :: D, E
    REAL(DP), DIMENSION(:,:),  INTENT(inout),OPTIONAL :: Z  

    CHARACTER(len=1)         :: JOBZ
    INTEGER                  :: N, info_inner, lwork
    REAL(DP),DIMENSION(SIZE(D),SIZE(D)) :: Z_inner
    REAL(DP), DIMENSION(:), ALLOCATABLE :: work


    N= SIZE(D)
    IF (SIZE(E) /= N) STOP 'DSTEV95: diagonal D and off-diagonal E must have same length! &
         & (Last element of E is not used.)'

    IF (PRESENT(z)) THEN
       IF (SIZE(Z,1)/= N) STOP 'First dimension of Z (array for eigenvectors) must be size(D)!' 
       IF (SIZE(Z,2)/= N) STOP 'Second dimension of Z (array for eigenvectors) must be size(D)!' 
       JOBZ='V'
    ELSE
       JOBZ='N'
    END IF

    lwork = MAX(1,2*N-2)
    ALLOCATE(work(lwork))
    CALL DSTEV( JOBZ, N, D, E, Z_inner, N, WORK, INFO_inner )
    DEALLOCATE(work)
    IF (PRESENT(z)) z=z_inner       
    IF (PRESENT(info)) info =info_inner
    
  END SUBROUTINE DSTEV95

  SUBROUTINE DSTEVX95( D, E, RANGE, M, VL, VU, IL, IU, Z, ifail, abstol, INFO )
    USE nrtype
    IMPLICIT NONE


    REAL(DP), DIMENSION(:),  INTENT(inout) :: D, E
    CHARACTER, INTENT(in)    :: range
    INTEGER,   INTENT(out)   :: M
    REAL(DP),  INTENT(in), OPTIONAL :: VL, VU, abstol
    INTEGER,   INTENT(in), OPTIONAL    :: IL, IU
    REAL(DP), DIMENSION(:,:),  INTENT(inout),OPTIONAL :: Z  
    INTEGER, DIMENSION(:), INTENT(out), OPTIONAL   :: ifail
    INTEGER,   INTENT(out),OPTIONAL  :: INFO

    CHARACTER(len=1)         :: JOBZ
    INTEGER                  :: N, info_inner, ldz
    REAL(DP),DIMENSION(SIZE(D),SIZE(D)) :: Z_inner
    REAL(DP), DIMENSION(SIZE(D)*5) :: work
    REAL(DP), DIMENSION(SIZE(D))   :: w
    INTEGER, DIMENSION(SIZE(D)*5)  :: iwork
    INTEGER, DIMENSION(SIZE(D))    :: ifail_inner
    REAL(DP) :: abstol_inner

    N= SIZE(D)
    ldz = SIZE(d)
    IF (SIZE(E) /= N) STOP 'DSTEVX95: diagonal D and off-diagonal E must have same length! &
         & (Last element of E is not used.)'

    abstol_inner = 1.d-7
    IF (PRESENT(abstol)) abstol_inner = abstol

    IF (PRESENT(z)) THEN
       IF (SIZE(Z,1)/= N) STOP 'First dimension of Z (array for eigenvectors) must be size(D)!' 
       IF (((range =='A').OR.(range =='a')) .AND.(SIZE(Z,1)/= N)) &
            & STOP 'Second dimension of Z (array for eigenvectors) must be size(D) if range = A!' 
       JOBZ='V'
    ELSE
       JOBZ='N'
    END IF

    IF ((RANGE == 'I') .OR. (RANGE == 'i')) THEN
       IF (.NOT.(PRESENT(IL) .AND. PRESENT(IU))) STOP 'DSTEVX95: if range ==I, then il and iu must be given!'
       IF (PRESENT(z)) THEN
          IF (SIZE(Z,2) < IU-IL+1) STOP &
               & 'Second dimension of Z (array for eigenvectors) must at least be number of desired eigenvectors!' 
       END IF
       CALL DSTEVX( JOBZ, 'I', N, D, E, 0.d0, 1.d0, IL, IU, ABSTOL_inner, &
            &                   M, W, Z_inner, LDZ, WORK, IWORK, IFAIL_inner, INFO_inner )
       D = W
       IF (PRESENT(z)) z = z_inner(1:SIZE(Z, 1), 1:SIZE(Z, 2)) 
       IF (PRESENT(info)) info =info_inner
       IF (PRESENT(ifail)) ifail =ifail_inner
    ELSE IF ((RANGE == 'V') .OR. (RANGE == 'v')) THEN
       IF (.NOT.(PRESENT(VL) .AND. PRESENT(VU))) STOP 'DSTEVX95: if range ==V, then vl and vu must be given!'
       CALL DSTEVX( JOBZ, 'V', N, D, E, VL, VU, 1, 2, ABSTOL_inner, &
            &                   M, W, Z_inner, LDZ, WORK, IWORK, IFAIL_inner, INFO_inner )
       D = W
       IF (PRESENT(z)) THEN
          IF (M > SIZE(Z, 2)) STOP 'DSTEVX95:More eigenvectors found than columns in Z!'
          z = z_inner(1:SIZE(Z, 1), 1:M)
       END IF
       IF (PRESENT(info)) info =info_inner
       IF (PRESENT(ifail)) ifail =ifail_inner
    ELSE IF ((RANGE == 'A') .OR. (RANGE == 'a')) THEN
       CALL DSTEVX( JOBZ, 'A', N, D, E, 0.d0, 1.d0, 1, 2, ABSTOL_inner, &
            &                   M, W, Z_inner, LDZ, WORK, IWORK, IFAIL_inner, INFO_inner )
       D = W
       IF (PRESENT(z)) z = z_inner
       IF (PRESENT(info)) info =info_inner
       IF (PRESENT(ifail)) ifail =ifail_inner
    ELSE
       STOP 'DSTEVX95: if range must be I, A or V!'
    END IF

    
  END SUBROUTINE DSTEVX95

  SUBROUTINE ZHESV95(Mat, vecs, uplo, info)
    USE nrtype
    IMPLICIT NONE
    INTEGER,   INTENT(out),OPTIONAL  :: INFO
    COMPLEX(DP), DIMENSION(:,:), INTENT(inout) :: Mat, vecs
    CHARACTER(len=1), INTENT(in)     :: uplo

    INTEGER                  :: N, LDA, LDB, NRHS, info_inner, lwork
    COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: work
    INTEGER,  DIMENSION(SIZE(Mat,2))   :: IPIV

    N= SIZE(Mat,2)
    lda=N
    
    NRHS = SIZE(vecs,2)
    ldb  = SIZE(vecs,1)

    IF (ldb/=N) STOP 'ZHESV95: second dimension of matrix and first dim. of vectors do not match!'
    lwork  =  N**2
    ALLOCATE(work(lwork))
    CALL ZHESV(UPLO, N, NRHS, Mat, LDA, IPIV, vecs, LDB, WORK, &
          & -1, INFO_inner)
    lwork  =  work(1)
    DEALLOCATE(work)
    ALLOCATE(work(lwork))
    CALL ZHESV(UPLO, N, NRHS, Mat, LDA, IPIV, vecs, LDB, WORK, &
          & LWORK, INFO_inner)

    IF (PRESENT(info)) info =info_inner

  END SUBROUTINE ZHESV95

  SUBROUTINE ZGESV95(Mat, vecs, info)
    USE nrtype
    IMPLICIT NONE
    INTEGER,   INTENT(out),OPTIONAL  :: INFO
    COMPLEX(DP), DIMENSION(:,:), INTENT(inout) :: Mat, vecs

    INTEGER                  :: N, LDA, LDB, NRHS, info_inner
    INTEGER,  DIMENSION(SIZE(Mat,2))   :: IPIV

    N= SIZE(Mat,2)
    lda=N
    
    NRHS = SIZE(vecs,2)
    ldb  = SIZE(vecs,1)

    IF (ldb/=N) STOP 'ZGESV95: second dimension of matrix and first dim. of vectors do not match!'
    CALL ZGESV(N, NRHS, Mat, LDA, IPIV, vecs, LDB, INFO_inner)
    IF (PRESENT(info)) info =info_inner

  END SUBROUTINE ZGESV95


  SUBROUTINE HBEVD95( AB,eig,uplo, psi,INFO )
    IMPLICIT NONE

    COMPLEX(DP),  DIMENSION(:,:),  INTENT(inout) :: AB
    REAL(DP),     DIMENSION(:),  INTENT(out)   :: eig
    CHARACTER(len=1),INTENT(in)  :: uplo
    COMPLEX(DP),  DIMENSION(:,:),  INTENT(inout), OPTIONAL :: psi
    INTEGER,   INTENT(out),OPTIONAL  :: INFO

    INTEGER  :: N, kd, info_inner, lwork,lrwork,liwork
    CHARACTER(len=1)         :: JOBZ
    COMPLEX(DP),DIMENSION(SIZE(ab,dim=2),SIZE(ab,dim=2)) :: Z_inner
    COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: work
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rwork
    INTEGER, DIMENSION(:), ALLOCATABLE :: iwork

    N=SIZE(ab,dim=2)
    IF (SIZE(eig) /= N) STOP 'Number of eigenvalues must equal dimesnion of matrix!' 
    IF (PRESENT(psi)) THEN
       IF (SIZE(psi,1)/= N) STOP 'First dimension of psi (array for eigenvectors) must be size(AB,2)!' 
       IF (SIZE(psi,2)/= N) STOP 'Second dimension of psi (array for eigenvectors) must be size(AB,2)!' 
       JOBZ='V'
    ELSE
       JOBZ='N'
    END IF
    kd=SIZE(ab,dim=1)-1

    ALLOCATE(work(MAX(1,2*N**2)),rwork(1 + 5*N + 2*N**2),iwork(3 + 5*N))
    CALL ZHBEVD( JOBZ, UPLO, N, KD, AB, kd+1, eig, z_inner, N, WORK, -1, &
         & RWORK, -1, IWORK,-1, INFO_inner )
    lwork = work(1)
    liwork = iwork(1)
    lrwork = rwork(1)
    DEALLOCATE(work,iwork,rwork)
    ALLOCATE(work(lwork),rwork(lrwork),iwork(liwork))

    CALL ZHBEVD( JOBZ, UPLO, N, KD, AB, kd+1, eig, z_inner, N, WORK, lwork, &
         & RWORK, lrwork, IWORK,liwork , INFO_inner )
    IF (PRESENT(psi)) psi=z_inner       
    IF (PRESENT(info)) info =info_inner
    DEALLOCATE(work,rwork,iwork)

  END SUBROUTINE HBEVD95

END MODULE diag_mod
