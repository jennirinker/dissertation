!**********************************************************************************************************************************
! Copyright (C) 2013-2014  National Renewable Energy Laboratory
!
! This code provides a wrapper for the LAPACK routines currently used at the NWTC (mainly codes in the FAST framework).
!
!**********************************************************************************************************************************
! File last committed: $Date: 2014-09-18 15:01:50 -0600 (Thu, 18 Sep 2014) $
! (File) Revision #: $Rev: 257 $
! URL: $HeadURL: https://windsvn.nrel.gov/NWTC_Library/branches/NetLib/NWTC_source/NWTC_LAPACK.f90 $
!**********************************************************************************************************************************
MODULE NWTC_LAPACK

   USE NWTC_Base        ! we only need the precision and error level constants
!   USE, INTRINSIC               :: ISO_C_Binding, only: C_FLOAT, C_DOUBLE          ! this is included in NWTC_Library

      ! Notes:

         ! Your project must include the following files:
         ! From the NWTC Subroutine Library:
         !     {Doub | Sing}Prec.f90 [from NWTC Library]
         !     Sys*.f90              [from NWTC Library]
         !     NWTC_Base.f90         [from NWTC Library]
         ! lapack library (preferable a binary, but available in source form from http://www.netlib.org/, too)
         ! This wrapper file:
         !     NWTC_LAPACK.f90

   !INTEGER, PARAMETER  :: Lib_ReKi = SiKi   !
   !INTEGER, PARAMETER  :: Lib_DbKi = R8Ki   ! DbKi
   !
   ! bjj: when using the built-in (or dynamic) lapack libraries, S=Real(SiKi); D=Real(R8Ki).
   !      if people are compiling the lapack source, S=real; D=double precision. (default real and doubles)
   !      we need to check this somehow to make sure the right routines are called.
   ! (or define a directive that uses

   IMPLICIT  NONE

   INTERFACE LAPACK_gbsv ! Computes the solution to system of linear equations A * X = B for GB matrices
      MODULE PROCEDURE LAPACK_dgbsv
      MODULE PROCEDURE LAPACK_sgbsv
   END INTERFACE

   INTERFACE LAPACK_gesv ! Computes the solution to system of linear equations A * X = B for GE matrices
      MODULE PROCEDURE LAPACK_dgesv
      MODULE PROCEDURE LAPACK_sgesv
   END INTERFACE

   INTERFACE LAPACK_getrf ! Factor matrix into A=PLU
      MODULE PROCEDURE LAPACK_dgetrf
      MODULE PROCEDURE LAPACK_sgetrf
   END INTERFACE

   INTERFACE LAPACK_getri ! Compute the inverse of a matrix using the LU factorization
      MODULE PROCEDURE LAPACK_dgetri
      MODULE PROCEDURE LAPACK_sgetri
   END INTERFACE

   INTERFACE LAPACK_getrs ! Solve system of linear equations Ax=PLUx=b
      MODULE PROCEDURE LAPACK_dgetrs
      MODULE PROCEDURE LAPACK_sgetrs
      MODULE PROCEDURE LAPACK_dgetrs1
      MODULE PROCEDURE LAPACK_sgetrs1
   END INTERFACE

   INTERFACE LAPACK_ggev ! Compute generalized eigenvalues and/or eigenvectors for a pair of N-by-N real nonsymmetric matrices (A,B)
      MODULE PROCEDURE LAPACK_dggev
      MODULE PROCEDURE LAPACK_sggev
   END INTERFACE

   INTERFACE LAPACK_posv ! Compute the solution to system of linear equations A * X = B for PO matrices
      MODULE PROCEDURE LAPACK_dposv
      MODULE PROCEDURE LAPACK_sposv
   END INTERFACE

   INTERFACE LAPACK_pptrf ! Compute the Cholesky factorization of a real symmetric positive definite matrix A stored in packed format.
      MODULE PROCEDURE LAPACK_dpptrf
      MODULE PROCEDURE LAPACK_spptrf
   END INTERFACE

      
CONTAINS

!=======================================================================
   SUBROUTINE LAPACK_DGBSV( N, KL, KU, NRHS, AB, IPIV, B, ErrStat, ErrMsg )


      ! passed parameters

      INTEGER,         intent(in   ) :: KL                ! The number of subdiagonals within the band of A.  KL >= 0.
      INTEGER,         intent(in   ) :: KU                ! The number of superdiagonals within the band of A.  KU >= 0.
      INTEGER,         intent(in   ) :: N                 ! The number of linear equations, i.e., the order of the matrix A.  N >= 0.
      INTEGER,         intent(in   ) :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.

      !     .. Array Arguments ..
      REAL(R8Ki)      ,intent(inout) :: AB( :, : )        ! On entry, the matrix A in band storage, in rows KL+1 to 2*KL+KU+1; rows 1 to KL of the array need not be set.
                                                          ! The j-th column of A is stored in the j-th column of the array AB as follows:
                                                          !    AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
                                                          ! On exit, details of the factorization: U is stored as an upper triangular band matrix with KL+KU superdiagonals in
                                                          ! rows 1 to KL+KU+1, and the multipliers used during the factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
      REAL(R8Ki)      ,intent(inout) :: B( :, : )         ! On entry, the N-by-NRHS right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.
      INTEGER,         intent(  out) :: IPIV( : )         ! The pivot indices that define the permutation matrix P; row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.
      INTEGER                        :: LDAB              ! The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
      INTEGER                        :: LDB               ! The leading dimension of the array B.   LDB  >= max(1,N).


      LDAB  = SIZE(AB,1)
      LDB   = SIZE(B, 1)



      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL dgbsv (N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGBSV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DGBSV: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGBSV
!=======================================================================
   SUBROUTINE LAPACK_SGBSV( N, KL, KU, NRHS, AB, IPIV, B, ErrStat, ErrMsg )


      ! passed parameters

      INTEGER,         intent(in   ) :: KL                ! The number of subdiagonals within the band of A.  KL >= 0.
      INTEGER,         intent(in   ) :: KU                ! The number of superdiagonals within the band of A.  KU >= 0.
      INTEGER,         intent(in   ) :: N                 ! The number of linear equations, i.e., the order of the matrix A.  N >= 0.
      INTEGER,         intent(in   ) :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.

      !     .. Array Arguments ..
      REAL(SiKi)      ,intent(inout) :: AB( :, : )        ! On entry, the matrix A in band storage, in rows KL+1 to 2*KL+KU+1; rows 1 to KL of the array need not be set.
                                                          ! The j-th column of A is stored in the j-th column of the array AB as follows:
                                                          !    AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
                                                          ! On exit, details of the factorization: U is stored as an upper triangular band matrix with KL+KU superdiagonals in
                                                          ! rows 1 to KL+KU+1, and the multipliers used during the factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
      REAL(SiKi)      ,intent(inout) :: B( :, : )         ! On entry, the N-by-NRHS right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.
      INTEGER,         intent(  out) :: IPIV( : )         ! The pivot indices that define the permutation matrix P; row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.

      INTEGER                        :: LDAB              ! The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
      INTEGER                        :: LDB               ! The leading dimension of the array B.   LDB  >= max(1,N).


      LDAB  = SIZE(AB,1)
      LDB   = SIZE(B, 1)


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SGBSV (N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGBSV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SGBSV: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGBSV
!=======================================================================
   SUBROUTINE LAPACK_DGESV ( N, A, IPIV, B, ErrStat, ErrMsg )


      ! passed parameters

      INTEGER,         intent(in   ) :: N                 ! The number of linear equations, i.e., the order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(R8Ki)      ,intent(inout) :: A( :, : )         ! On entry, the N-by-N coefficient matrix A.  On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored.
      REAL(R8Ki)      ,intent(inout) :: B( :, : )         ! On entry, the N-by-NRHS matrix of right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.
      INTEGER,         intent(  out) :: IPIV( : )         ! The pivot indices that define the permutation matrix P; row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).

      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = SIZE(B,2)

      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGESV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DGESV: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGESV
!=======================================================================
   SUBROUTINE LAPACK_SGESV ( N, A, IPIV, B, ErrStat, ErrMsg )


      ! passed parameters

      INTEGER,         intent(in   ) :: N                 ! The number of linear equations, i.e., the order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(SiKi)      ,intent(inout) :: A( :, : )       ! On entry, the N-by-N coefficient matrix A.  On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored.
      REAL(SiKi)      ,intent(inout) :: B( :, : )       ! On entry, the N-by-NRHS matrix of right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.
      INTEGER,         intent(  out) :: IPIV( : )         ! The pivot indices that define the permutation matrix P; row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = SIZE(B,2)


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGESV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SGESV: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGESV
!=======================================================================
   SUBROUTINE LAPACK_DGETRF( M, N, A, IPIV, ErrStat, ErrMsg )


      ! passed parameters

      INTEGER,         intent(in   ) :: M                 ! The number of rows of the matrix A.  M >= 0.
      INTEGER,         intent(in   ) :: N                 ! The number of columns of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(R8Ki)      ,intent(inout) :: A( :, : )         ! On entry, the M-by-N matrix to be factored. On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored.
      INTEGER,         intent(  out) :: IPIV( : )         ! The pivot indices; for 1 <= i <= min(M,N), row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero. The factor U is exactly singular.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,M).

      LDA  = SIZE(A,1)


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DGETRF( M, N, A, LDA, IPIV, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGETRF: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DGETRF: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGETRF
!=======================================================================
   SUBROUTINE LAPACK_SGETRF( M, N, A, IPIV, ErrStat, ErrMsg )


      ! passed parameters

      INTEGER,        intent(in   ) :: M                 ! The number of rows of the matrix A.  M >= 0.
      INTEGER,        intent(in   ) :: N                 ! The number of columns of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(SiKi)     ,intent(inout) :: A( :, : )         ! On entry, the M-by-N matrix to be factored. On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored.
      INTEGER,        intent(  out) :: IPIV( : )         ! The pivot indices; for 1 <= i <= min(M,N), row i of the matrix was interchanged with row IPIV(i).

      INTEGER(IntKi), intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),   intent(  out) :: ErrMsg            ! Message describing error

         ! local variables
      INTEGER                       :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero. The factor U is exactly singular.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,M).

      LDA  = SIZE(A,1)



      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SGETRF( M, N, A, LDA, IPIV, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGETRF: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SGETRF: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Factor U is exactly singular.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGETRF
!=======================================================================
   SUBROUTINE LAPACK_DGETRS( TRANS, N, A, IPIV, B, ErrStat, ErrMsg )

      ! passed parameters

      CHARACTER(1),    intent(in   ) :: TRANS             ! Specifies the form of the system of equations: = 'N':  A * X = B  (No transpose)
                                                          !                                                = 'T':  A**T* X = B  (Transpose)
                                                          !                                                = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
      INTEGER,         intent(in   ) :: N                 ! The order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      INTEGER,         intent(in   ) :: IPIV( : )         ! The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(R8Ki)      ,intent(in   ) :: A( :, : )         ! The factors L and U from the factorization A = P*L*U as computed by DGETRF.
      REAL(R8Ki)      ,intent(inout) :: B( :, : )         ! On entry, the right hand side matrix B. On exit, the solution matrix X.

      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value;
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = SIZE(B,2)


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) -INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGETRS: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DGETRS: unknown error -'//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGETRS
!=======================================================================
   SUBROUTINE LAPACK_DGETRS1( TRANS, N, A, IPIV, B, ErrStat, ErrMsg )

      ! passed parameters

      CHARACTER(1),    intent(in   ) :: TRANS             ! Specifies the form of the system of equations: = 'N':  A * X = B  (No transpose)
                                                          !                                                = 'T':  A**T* X = B  (Transpose)
                                                          !                                                = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
      INTEGER,         intent(in   ) :: N                 ! The order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      INTEGER,         intent(in   ) :: IPIV( : )         ! The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(R8Ki)      ,intent(in   ) :: A( :, : )         ! The factors L and U from the factorization A = P*L*U as computed by DGETRF.
      REAL(R8Ki)      ,intent(inout) :: B( :    )         ! On entry, the right hand side matrix B. On exit, the solution matrix X.

      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value;
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = 1


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) -INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGETRS1: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DGETRS1: unknown error -'//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGETRS1
!=======================================================================
   SUBROUTINE LAPACK_SGETRS( TRANS, N, A, IPIV, B, ErrStat, ErrMsg )

      ! passed parameters

      CHARACTER(1),   intent(in   )  :: TRANS             ! Specifies the form of the system of equations: = 'N':  A * X = B  (No transpose)
                                                          !                                                = 'T':  A**T* X = B  (Transpose)
                                                          !                                                = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
      INTEGER,        intent(in   )  :: N                 ! The order of the matrix A.  N >= 0.
                                     
      !     .. Array Arguments ..    
      INTEGER,        intent(in   )  :: IPIV( : )         ! The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(SiKi),     intent(in   )  :: A( :, : )         ! The factors L and U from the factorization A = P*L*U as computed by SGETRF.
      REAL(SiKi),     intent(inout)  :: B( :, : )         ! On entry, the right hand side matrix B. On exit, the solution matrix X.
                                     
      INTEGER(IntKi), intent(  out)  :: ErrStat           ! Error level
      CHARACTER(*),   intent(  out)  :: ErrMsg            ! Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value;
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = SIZE(B,2)


      ErrStat = ErrID_None
      ErrMsg  = ""

      !IF (ReKi == C_FLOAT) THEN
         CALL SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      !ELSEIF (ReKi == C_DOUBLE) THEN
      !   CALL DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      !ELSE
      !   ErrStat = ErrID_FATAL
      !   ErrMsg  = "LAPACK_SGETRS: Matrix A is an invalid type."
      !   RETURN
      !END IF

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) -INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGETRS: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SGETRS: unknown error -'//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGETRS
!=======================================================================
   SUBROUTINE LAPACK_SGETRS1( TRANS, N, A, IPIV, B, ErrStat, ErrMsg )

      ! passed parameters

      CHARACTER(1),   intent(in   )  :: TRANS             ! Specifies the form of the system of equations: = 'N':  A * X = B  (No transpose)
                                                          !                                                = 'T':  A**T* X = B  (Transpose)
                                                          !                                                = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
      INTEGER,        intent(in   )  :: N                 ! The order of the matrix A.  N >= 0.
                                     
      !     .. Array Arguments ..    
      INTEGER,        intent(in   )  :: IPIV( : )         ! The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(SiKi),     intent(in   )  :: A( :, : )         ! The factors L and U from the factorization A = P*L*U as computed by SGETRF.
      REAL(SiKi),     intent(inout)  :: B( :    )         ! On entry, the right hand side matrix B. On exit, the solution matrix X.
                                     
      INTEGER(IntKi), intent(  out)  :: ErrStat           ! Error level
      CHARACTER(*),   intent(  out)  :: ErrMsg            ! Message describing error
                                     
         ! local variables           
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value;
      INTEGER                        :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
      NRHS = 1


      ErrStat = ErrID_None
      ErrMsg  = ""

      
      !IF (ReKi == C_FLOAT) THEN
         CALL SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      !ELSEIF (ReKi == C_DOUBLE) THEN
      !   CALL DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      !ELSE
      !   ErrStat = ErrID_FATAL
      !   ErrMsg  = "LAPACK_SGETRS: Matrix A is an invalid type."
      !   RETURN
      !END IF

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) -INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGETRS1: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SGETRS1: unknown error -'//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGETRS1
!=======================================================================
   SUBROUTINE LAPACK_DGETRI( N, A, IPIV, WORK, LWORK, ErrStat, ErrMsg )

      ! passed parameters

      INTEGER,         intent(in   ) :: N                 ! The order of the matrix A.  N >= 0.
      INTEGER,         intent(in   ) :: LWORK             ! The dimension of the array WORK. LWORK >= max(1,N). For optimal performance LWORK >= N*NB, where NB is the optimal blocksize returned by ILAENV.
                                                          ! If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns this value as the first
                                                          ! entry of the WORK array, and no error message related to LWORK is issued by XERBLA.

      !     .. Array Arguments ..
      INTEGER,         intent(in   ) :: IPIV( : )         ! dimension (N). The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(R8Ki)      ,intent(inout) :: A( :, : )         ! On entry, the factors L and U from the factorization A = P*L*U as computed by DGETRF. On exit, if INFO = 0, the inverse of the original matrix A.
      REAL(R8Ki)      ,intent(inout) :: WORK( : )         ! On exit, if INFO=0, then WORK(1) returns the optimal LWORK.

      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero. The matrix is singular and its inverse could not be computed.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,M).

      LDA  = SIZE(A,1)

      ErrStat = ErrID_None
      ErrMsg  = ""

      !IF (DbKi == C_DOUBLE) THEN
         CALL DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      !ELSEIF (DbKi == C_FLOAT) THEN
      !   CALL DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      !ELSE
      !   ErrStat = ErrID_FATAL
      !   ErrMsg  = "LAPACK_DGETRI: Matrix A is an invalid type."
      !   RETURN
      !END IF

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGETRI: illegal value in argument "//TRIM(ErrMsg)//"."
        ELSE
            ErrMsg = 'LAPACK_DGETRI: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Matrix is singular and its inverse cannot be computed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGETRI
!=======================================================================
   SUBROUTINE LAPACK_SGETRI( N, A, IPIV, WORK, LWORK, ErrStat, ErrMsg )

      ! passed parameters

      INTEGER,        intent(in   )  :: N                 ! The order of the matrix A.  N >= 0.
      INTEGER,        intent(in   )  :: LWORK             ! The dimension of the array WORK. LWORK >= max(1,N). For optimal performance LWORK >= N*NB, where NB is the optimal blocksize returned by ILAENV.
                                                          ! If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns this value as the first
                                                          ! entry of the WORK array, and no error message related to LWORK is issued by XERBLA.
                                     
      !     .. Array Arguments ..    
      INTEGER,        intent(in   )  :: IPIV( : )         ! dimension (N). The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix was interchanged with row IPIV(i).
      REAL(SiKi),     intent(inout)  :: A( :, : )         ! On entry, the factors L and U from the factorization A = P*L*U as computed by SGETRF. On exit, if INFO = 0, the inverse of the original matrix A.
      REAL(SiKi),     intent(inout)  :: WORK( : )         ! On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
                                     
      INTEGER(IntKi), intent(  out)  :: ErrStat           ! Error level
      CHARACTER(*),   intent(  out)  :: ErrMsg            ! Message describing error
                                     
         ! local variables           
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero. The matrix is singular and its inverse could not be computed.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,M).

      LDA  = SIZE(A,1)

      ErrStat = ErrID_None
      ErrMsg  = ""


      CALL SGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGETRI: illegal value in argument "//TRIM(ErrMsg)//"."
        ELSE
            ErrMsg = 'LAPACK_SGETRI: U('//TRIM(ErrMsg)//','//TRIM(ErrMsg)//')=0. Matrix is singular and its inverse cannot be computed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGETRI
!=======================================================================
   SUBROUTINE LAPACK_DGGEV(JOBVL, JOBVR, N, A, B, ALPHAR, ALPHAI, BETA, VL, VR, WORK, LWORK, ErrStat, ErrMsg)

      ! passed variables/parameters:

      CHARACTER(1),    intent(in   ) :: JOBVL             ! = 'N':  do not compute the left generalized eigenvectors; = 'V':  compute the left generalized eigenvectors.
      CHARACTER(1),    intent(in   ) :: JOBVR             ! = 'N':  do not compute the right generalized eigenvectors; = 'V':  compute the right generalized eigenvectors.

      INTEGER,         intent(in   ) :: N                 ! The order of the matrices A, B, VL, and VR.  N >= 0.

      INTEGER,         intent(in   ) :: LWORK             ! The dimension of the array WORK.  LWORK >= max(1,8*N). For good performance, LWORK must generally be larger.
                                                          !   If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns
                                                          !   this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.


      REAL(R8Ki)      ,intent(inout) :: A( :, : )         ! dimension (LDA, N). On entry, the matrix A in the pair (A,B). On exit, A has been overwritten.
      REAL(R8Ki)      ,intent(inout) :: B( :, : )         ! dimension (LDB, N). On entry, the matrix B in the pair (A,B). On exit, B has been overwritten.

      REAL(R8Ki)      ,intent(  out) :: ALPHAR( : )       ! dimension (N). See comments for variable "Beta"
      REAL(R8Ki)      ,intent(  out) :: ALPHAI( : )       ! dimension (N). See comments for variable "Beta".
      REAL(R8Ki)      ,intent(  out) :: BETA( : )         ! On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will be the generalized eigenvalues.  If ALPHAI(j) is zero, then
                                                          !   the j-th eigenvalue is real; if positive, then the j-th and (j+1)-st eigenvalues are a complex conjugate pair, with
                                                          !   ALPHAI(j+1) negative.
                                                          !
                                                          !   Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j) may easily over- or underflow, and BETA(j) may even be zero.
                                                          !   Thus, the user should avoid naively computing the ratio alpha/beta.  However, ALPHAR and ALPHAI will be always less
                                                          !   than and usually comparable with norm(A) in magnitude, and BETA always less than and usually comparable with norm(B).


      REAL(R8Ki)      ,intent(  out) :: VL( :, : )        ! dimension (LDVL,N). If JOBVL = 'V', the left eigenvectors u(j) are stored one after another in the columns of VL, in the same
                                                          !   order as their eigenvalues. If the j-th eigenvalue is real, then u(j) = VL(:,j), the j-th column of VL. If the j-th and
                                                          !   (j+1)-th eigenvalues form a complex conjugate pair, then u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
                                                          !   Each eigenvector is scaled so the largest component has abs(real part)+abs(imag. part)=1. Not referenced if JOBVL = 'N'.
      REAL(R8Ki)      ,intent(  out) :: VR( :, : )        ! dimension (LDVR,N). If JOBVR = 'V', the right eigenvectors v(j) are stored one after another in the columns of VR, in the same
                                                          !   order as their eigenvalues. If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR. If the j-th and
                                                          !   (j+1)-th eigenvalues form a complex conjugate pair, then v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
                                                          !   Each eigenvector is scaled so the largest component has abs(real part)+abs(imag. part)=1. Not referenced if JOBVR = 'N'.

      REAL(R8Ki)      ,intent(inout) :: WORK( : )         ! dimension (MAX(1,LWORK)). On exit, if INFO = 0, WORK(1) returns the optimal LWORK.



      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error


         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit;
                                                          ! < 0:
                                                          !   = -i, the i-th argument had an illegal value;
                                                          ! > 0:
                                                          !   = 1,...,N: The QZ iteration failed.  No eigenvectors have been calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) should be correct for j=INFO+1,...,N.
                                                          !   = N+1: other than QZ iteration failed in DHGEQZ.
                                                          !   = N+2: error return from DTGEVC.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
      INTEGER                        :: LDVL              ! The leading dimension of the matrix VL. LDVL >= 1, and if JOBVL = 'V', LDVL >= N
      INTEGER                        :: LDVR              ! The leading dimension of the matrix VR. LDVR >= 1, and if JOBVR = 'V', LDVR >= N.
      CHARACTER(20)                  :: n_str
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)

      LDVL  = SIZE(VL,1)
      LDVR  = SIZE(VR,1)
      
      

      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DGGEV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSEIF (INFO <N) THEN
            !ErrStat = ErrID_Severe
            WRITE( ErrMsg, * ) INFO + 1
            WRITE( n_str, * ) n
            ErrMsg  = "LAPACK_DGGEV: The QZ iteration failed. No eigenvectors have been calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) should be correct for j="&
                       //TRIM(ErrMsg)//",...,"//TRIM(n_str)//"."
         ELSEIF (INFO == N ) THEN
            ErrMsg  = "LAPACK_DGGEV: The QZ iteration failed. No eigenvectors have been calculated."
         ELSEIF (INFO == N+1) THEN
            ErrMsg  = "LAPACK_DGGEV: other than QZ iteration failed in DHGEQZ."
         ELSEIF (INFO == N+2) THEN
            ErrMsg  = "LAPACK_DGGEV: error return from DTGEVC."
         ELSE
            ErrMsg = 'LAPACK_DGGEV: unknown error '//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DGGEV
!=======================================================================
   SUBROUTINE LAPACK_SGGEV(JOBVL, JOBVR, N, A, B, ALPHAR, ALPHAI, BETA, VL, VR, WORK, LWORK, ErrStat, ErrMsg)

      ! subroutine arguments

      CHARACTER(1),   intent(in   )  :: JOBVL             ! = 'N':  do not compute the left generalized eigenvectors; = 'V':  compute the left generalized eigenvectors.
      CHARACTER(1),   intent(in   )  :: JOBVR             ! = 'N':  do not compute the right generalized eigenvectors; = 'V':  compute the right generalized eigenvectors.
                                    
      INTEGER,        intent(in   )  :: N                 ! The order of the matrices A, B, VL, and VR.  N >= 0.
      INTEGER,        intent(in   )  :: LWORK             ! The dimension of the array WORK.  LWORK >= max(1,8*N). For good performance, LWORK must generally be larger.
                                                          !   If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns
                                                          !   this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.
                                    
                                    
      REAL(SiKi),     intent(inout)  :: A( :, : )         ! dimension (LDA, N). On entry, the matrix A in the pair (A,B). On exit, A has been overwritten.
      REAL(SiKi),     intent(inout)  :: B( :, : )         ! dimension (LDB, N). On entry, the matrix B in the pair (A,B). On exit, B has been overwritten.
                                    
                                    
      REAL(SiKi),     intent(  out)  :: ALPHAR( : )       ! dimension (N). See comments for variable "Beta"
      REAL(SiKi),     intent(  out)  :: ALPHAI( : )       ! dimension (N). See comments for variable "Beta".
      REAL(SiKi),     intent(  out)  :: BETA( : )         ! On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will be the generalized eigenvalues.  If ALPHAI(j) is zero, then
                                                          !   the j-th eigenvalue is real; if positive, then the j-th and (j+1)-st eigenvalues are a complex conjugate pair, with
                                                          !   ALPHAI(j+1) negative.
                                                          !
                                                          !   Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j) may easily over- or underflow, and BETA(j) may even be zero.
                                                          !   Thus, the user should avoid naively computing the ratio alpha/beta.  However, ALPHAR and ALPHAI will be always less
                                                          !   than and usually comparable with norm(A) in magnitude, and BETA always less than and usually comparable with norm(B).


      REAL(SiKi),     intent(  out)  :: VL( :, : )        ! dimension (LDVL,N). If JOBVL = 'V', the left eigenvectors u(j) are stored one after another in the columns of VL, in the same
                                                          !   order as their eigenvalues. If the j-th eigenvalue is real, then u(j) = VL(:,j), the j-th column of VL. If the j-th and
                                                          !   (j+1)-th eigenvalues form a complex conjugate pair, then u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
                                                          !   Each eigenvector is scaled so the largest component has abs(real part)+abs(imag. part)=1. Not referenced if JOBVL = 'N'.
      REAL(SiKi),     intent(  out)  :: VR( :, : )        ! dimension (LDVR,N). If JOBVR = 'V', the right eigenvectors v(j) are stored one after another in the columns of VR, in the same
                                                          !   order as their eigenvalues. If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR. If the j-th and
                                                          !   (j+1)-th eigenvalues form a complex conjugate pair, then v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
                                                          !   Each eigenvector is scaled so the largest component has abs(real part)+abs(imag. part)=1. Not referenced if JOBVR = 'N'.
                                     
      REAL(SiKi),     intent(inout)  :: WORK( : )         ! dimension (MAX(1,LWORK)). On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
                                     
                                     
                                     
      INTEGER(IntKi), intent(  out)  :: ErrStat           ! Error level
      CHARACTER(*),   intent(  out)  :: ErrMsg            ! Message describing error
                                     
         ! local variables           
      INTEGER                        :: INFO              ! = 0:  successful exit;
                                                          ! < 0:
                                                          !   = -i, the i-th argument had an illegal value;
                                                          ! > 0:
                                                          !   = 1,...,N: The QZ iteration failed.  No eigenvectors have been calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) should be correct for j=INFO+1,...,N.
                                                          !   = N+1: other than QZ iteration failed in SHGEQZ.
                                                          !   = N+2: error return from STGEVC.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).
      INTEGER                        :: LDVL              ! The leading dimension of the matrix VL. LDVL >= 1, and if JOBVL = 'V', LDVL >= N
      INTEGER                        :: LDVR              ! The leading dimension of the matrix VR. LDVR >= 1, and if JOBVR = 'V', LDVR >= N.
      CHARACTER(20)                  :: n_str
                        
      
      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)

      LDVL  = SIZE(VL,1)
      LDVR  = SIZE(VR,1)


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SGGEV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSEIF (INFO <=N) THEN
            !ErrStat = ErrID_Severe
            WRITE( ErrMsg, * ) INFO + 1
            WRITE( n_str, * ) n
            ErrMsg  = "LAPACK_SGGEV: The QZ iteration failed. No eigenvectors have been calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) should be correct for j="&
                       //TRIM(ErrMsg)//",...,"//TRIM(n_str)//"."
         ELSEIF (INFO == N ) THEN
            ErrMsg  = "LAPACK_SGGEV: The QZ iteration failed. No eigenvectors have been calculated."
         ELSEIF (INFO == N+1) THEN
            ErrMsg  = "LAPACK_SGGEV: other than QZ iteration failed in SHGEQZ."
         ELSEIF (INFO == N+2) THEN
            ErrMsg  = "LAPACK_SGGEV: error return from STGEVC."
         ELSE
            ErrMsg = 'LAPACK_SGGEV: unknown error '//TRIM(ErrMsg)//'.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SGGEV
!=======================================================================
   SUBROUTINE LAPACK_DPOSV (UPLO, N, NRHS, A, B, ErrStat, ErrMsg)


      ! passed parameters

      INTEGER,         intent(in   ) :: N                 ! The number of linear equations, i.e., the order of the matrix A.  N >= 0.
      INTEGER,         intent(in   ) :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.

      !     .. Array Arguments ..
      REAL(R8Ki)      ,intent(inout) :: A( :, : )       ! On entry, the symmetric matrix A.  If UPLO = 'U', the leading N-by-N upper triangular part of A contains the upper
                                                          ! triangular part of the matrix A, and the strictly lower triangular part of A is not referenced.  If UPLO = 'L', the
                                                          ! leading N-by-N lower triangular part of A contains the lower triangular part of the matrix A, and the strictly upper
                                                          ! triangular part of A is not referenced.
                                                          ! On exit, if INFO = 0, the factor U or L from the Cholesky factorization A = U**T*U or A = L*L**T.
      REAL(R8Ki)      ,intent(inout) :: B( :, : )       ! On entry, the N-by-NRHS matrix of right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.

      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error
      CHARACTER(1),    intent(in   ) :: UPLO              ! 'U':  Upper triangle of A is stored; 'L':  Lower triangle of A is stored.

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, the leading minor of order i of A is not positive definite, so the factorization could not be
                                                          ! completed, and the solution has not been computed.

      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).


      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)
                                                          
                                                          
      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DPOSV (UPLO, N, NRHS, A, LDA, B, LDB, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DPOSV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DPOSV: Leading minor order '//TRIM(ErrMsg)//' of A is not positive definite, so factorization could not be completed, and the solution has not been computed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DPOSV
!=======================================================================
   SUBROUTINE LAPACK_SPOSV (UPLO, N, NRHS, A, B, ErrStat, ErrMsg)


      ! passed parameters

      INTEGER,         intent(in   ) :: N                 ! The number of linear equations, i.e., the order of the matrix A.  N >= 0.
      INTEGER,         intent(in   ) :: NRHS              ! The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.

      !     .. Array Arguments ..
      REAL(SiKi)      ,intent(inout) :: A( :, : )         ! On entry, the symmetric matrix A.  If UPLO = 'U', the leading N-by-N upper triangular part of A contains the upper
                                                          ! triangular part of the matrix A, and the strictly lower triangular part of A is not referenced.  If UPLO = 'L', the
                                                          ! leading N-by-N lower triangular part of A contains the lower triangular part of the matrix A, and the strictly upper
                                                          ! triangular part of A is not referenced.
                                                          ! On exit, if INFO = 0, the factor U or L from the Cholesky factorization A = U**T*U or A = L*L**T.
      REAL(SiKi)      ,intent(inout) :: B( :, : )         ! On entry, the N-by-NRHS matrix of right hand side matrix B.  On exit, if INFO = 0, the N-by-NRHS solution matrix X.

      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error
      CHARACTER(1),    intent(in   ) :: UPLO              ! 'U':  Upper triangle of A is stored; 'L':  Lower triangle of A is stored.

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; > 0: if INFO = i, U(i,i) is exactly zero.
                                                          ! > 0: if INFO = i, the leading minor of order i of A is not positive definite, so the factorization could not be
                                                          ! completed, and the solution has not been computed.
      INTEGER                        :: LDA               ! The leading dimension of the array A.  LDA >= max(1,N).
      INTEGER                        :: LDB               ! The leading dimension of the array B.  LDB >= max(1,N).


      LDA  = SIZE(A,1)
      LDB  = SIZE(B,1)



      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SPOSV (UPLO, N, NRHS, A, LDA, B, LDB, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SPOSV: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SPOSV: Leading minor order '//TRIM(ErrMsg)//' of A is not positive definite, so factorization could not be completed, and the solution has not been computed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SPOSV
!=======================================================================
   SUBROUTINE LAPACK_DPPTRF (UPLO, N, AP, ErrStat, ErrMsg)

   ! DPPTRF computes the Cholesky factorization of a real symmetric
   ! positive definite matrix A stored in packed format.
   !
   ! The factorization has the form
   !      A = U**T * U,  if UPLO = 'U', or
   !      A = L  * L**T,  if UPLO = 'L',
   ! where U is an upper triangular matrix and L is lower triangular.
   

      ! passed parameters

      INTEGER,         intent(in   ) :: N                 ! The order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(R8Ki)      ,intent(inout) :: AP( : )           ! AP is REAL array, dimension (N*(N+1)/2)
                                                          ! On entry, the upper or lower triangle of the symmetric matrix A, packed columnwise in a linear array.  The j-th column of A
                                                          ! is stored in the array AP as follows:
                                                          !    if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
                                                          !    if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
                                                          ! See below for further details.
                                                          ! On exit, if INFO = 0, the triangular factor U or L from the Cholesky factorization A = U**T*U or A = L*L**T, in the same storage format as A.

                                                          ! Further details:      
                                                          !   The packed storage scheme is illustrated by the following example
                                                          !   when N = 4, UPLO = 'U':
                                                          ! 
                                                          !   Two-dimensional storage of the symmetric matrix A:
                                                          ! 
                                                          !      a11 a12 a13 a14
                                                          !          a22 a23 a24
                                                          !              a33 a34     (aij = aji)
                                                          !                  a44
                                                          ! 
                                                          !   Packed storage of the upper triangle of A:
                                                          ! 
                                                          !   AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
                                                                                                                    
                                                          
      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error
      CHARACTER(1),    intent(in   ) :: UPLO              ! 'U':  Upper triangle of A is stored; 'L':  Lower triangle of A is stored.
      
      

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; 
                                                          ! > 0:  if INFO = i, the leading minor of order i is not positive definite, and the factorization could not be completed.

                                                          
      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL DPPTRF (UPLO, N, AP, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_DPPTRF: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_DPPTRF: Leading minor order '//TRIM(ErrMsg)//' of A is not positive definite, so Cholesky factorization could not be completed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_DPPTRF   
!=======================================================================
   SUBROUTINE LAPACK_SPPTRF (UPLO, N, AP, ErrStat, ErrMsg)

   ! SPPTRF computes the Cholesky factorization of a real symmetric
   ! positive definite matrix A stored in packed format.
   !
   ! The factorization has the form
   !      A = U**T * U,  if UPLO = 'U', or
   !      A = L  * L**T,  if UPLO = 'L',
   ! where U is an upper triangular matrix and L is lower triangular.
   

      ! passed parameters

      INTEGER,         intent(in   ) :: N                 ! The order of the matrix A.  N >= 0.

      !     .. Array Arguments ..
      REAL(SiKi)      ,intent(inout) :: AP( : )           ! AP is REAL array, dimension (N*(N+1)/2)
                                                          ! On entry, the upper or lower triangle of the symmetric matrix A, packed columnwise in a linear array.  The j-th column of A
                                                          ! is stored in the array AP as follows:
                                                          !    if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
                                                          !    if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
                                                          ! See LAPACK_DPPTRF for further details.
                                                          ! On exit, if INFO = 0, the triangular factor U or L from the Cholesky factorization A = U**T*U or A = L*L**T, in the same storage format as A.

      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error
      CHARACTER(1),    intent(in   ) :: UPLO              ! 'U':  Upper triangle of A is stored; 'L':  Lower triangle of A is stored.

         ! local variables
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value; 
                                                          ! > 0:  if INFO = i, the leading minor of order i is not positive definite, and the factorization could not be completed


      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL SPPTRF (UPLO, N, AP, INFO)

      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "LAPACK_SPPTRF: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'LAPACK_SPPTRF: Leading minor order '//TRIM(ErrMsg)//' of A is not positive definite, so Cholesky factorization could not be completed.'
         END IF
      END IF


   RETURN
   END SUBROUTINE LAPACK_SPPTRF   
   !=======================================================================
END MODULE NWTC_LAPACK
