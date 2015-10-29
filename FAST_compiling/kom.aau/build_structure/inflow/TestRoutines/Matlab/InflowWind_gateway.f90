!****************************************************************************************************
!Gateway program to call WindInflow_mex
!
!  16 Nov 2009    B. Jonkman
!  National Renewable Energy Laboratory
!  National Wind Technology Center
!
!     To initialize the subroutine:
!  [Err] = InflowWind_mex( 0, FileName, FileType, ReferenceHeight, ReferenceWidth )
!
!     To get wind speeds:
!  [Err, V] = InflowWind_mex( 3, t, X )
!
!     To terminate:
!  [Err] = InflowWind_mex( 9 )
!
!     Err             = INTEGER(8)         error status
!     V               = REAL(8) ARRAY(3)   velocity
!     FileName        = CHARACTER(*)
!     FileType        = REAL(8)
!     ReferenceHeight = REAL(8)
!     ReferenceWidth  = REAL(8)
!     t               = REAL(8)            time
!     x               = REAL(8) ARRAY(3)   position
!
!====================================================================================================
SUBROUTINE mexFunction(nlhs, plhs, nrhs, prhs)
!====================================================================================================
!   Purpose:  Glue routine for making FORTRAN MEX-file systems and blocks
!====================================================================================================

   USE                           InflowWind
   USE                           SharedInflowDefns
       
   IMPLICIT                      NONE
   

      !----------------------------------------------------------------------------------------------
      ! define parameters
      !----------------------------------------------------------------------------------------------

   INTEGER(4), PARAMETER        :: mxREAL      = 0         ! MATLAB uses 0 for REAL numbers; 1 for COMPLEX
   INTEGER,    PARAMETER        :: mwSize      = 4         ! Size of pointer variables; replace 4 with 8 on DE! Alpha and the SGI 64-bit platforms
   INTEGER,    PARAMETER        :: mwPointer   = 4         ! Size of pointer variables; replace 4 with 8 on DE! Alpha and the SGI 64-bit platforms
   INTEGER,    PARAMETER        :: mxDB        = 8         ! Matlab requires double-precision reals

      !----------------------------------------------------------------------------------------------
      ! define right- and left-hand side arguments of the MATLAB function we're creating
      !----------------------------------------------------------------------------------------------

   INTEGER                      :: nlhs                    ! MATLAB's count of the number of left-hand (output) arguments
   INTEGER                      :: nrhs                    ! MATLAB's count of the number of right-hand (input) arguments
   INTEGER(mwPointer)           :: plhs(*)                 ! MATLAB's pointer(s) to left-hand (output) arguments
   INTEGER(mwPointer)           :: prhs(*)                 ! MATLAB's pointer(s) to right-hand (input) arguments

   INTEGER(mwPointer)           :: ptr_X                   ! pointer to input  (RHS) argument #2 (x)
   INTEGER(mwPointer)           :: ptr_Y                   ! pointer to output (LHS) argument #1 (outputs)

   INTEGER                      :: FLAG                     ! input  argument #1, FLAG

      !----------------------------------------------------------------------------------------------
      ! define internal variables
      !----------------------------------------------------------------------------------------------

   TYPE(InflInitInfo)           :: InitWindData
   TYPE(InflIntrpOut)           :: MyWindSpeed

   INTEGER                      :: ArgNum
   INTEGER                      :: ErrStat
   INTEGER                      :: ITime
   INTEGER                      :: IPos
   INTEGER(mwSize)              :: M                       ! Number of rows in array
   INTEGER(mwSize)              :: N                       ! Number of columns in array
   INTEGER(mwSize)              :: NT
   INTEGER(mwSize)              :: NX
   INTEGER(mwSize)              :: NV
   
   REAL(ReKi)                   :: InputPosition(3)
   REAL(ReKi)                   :: Time
   REAL(mxDB), ALLOCATABLE      :: Velocity(:,:)
   REAL(mxDB), ALLOCATABLE      :: TimeAry(:)
   REAL(mxDB), ALLOCATABLE      :: X (:,:)                  ! input  argument #3, INPUT POSITION ARRAY, for FLAG=3

      !----------------------------------------------------------------------------------------------
      ! define the EXTERNAL MATLAB procedures
      !----------------------------------------------------------------------------------------------

   EXTERNAL                     :: mxCopyPtrToReal8        ! MATLAB mex function to create REAL(8) array from pointer to an array
   EXTERNAL                     :: mxCopyReal8ToPtr        ! MATLAB mex function to create pointer to copy of a REAL(8) array
   INTEGER(mwPointer), EXTERNAL :: mxCreateDoubleMatrix    ! MATLAB mex pointer to new matrix [Replace integer by integer*8 on the DE! Alpha and the SGI 64-bit platforms]
   INTEGER(mwSize),    EXTERNAL :: mxGetM                  ! MATALB mex function get number of rows in array
   INTEGER(mwSize),    EXTERNAL :: mxGetN                  ! MATALB mex function to get number of columns in array
   INTEGER(mwPointer), EXTERNAL :: mxGetPr                 ! MATLAB mex function to get the address of the first real number [Replace integer by integer*8 on the DE! Alpha and the SGI 64-bit platforms]
   REAL(mxDB),         EXTERNAL :: mxGetScalar             ! MATLAB mex function to return a scalar ( like mxCopyPtrToReal8() for 1 element )
   INTEGER(mwPointer), EXTERNAL :: mxGetString             ! MATLAB mex function to get string from its pointer



   !=================================================================================================
   ! Check the function input arguments to determine if this is an initialization step, termination,
   ! or request to obtain the wind speed.
   !=================================================================================================
   
   IF (NRHS == 0) THEN
      CALL ProgAbort( ' Must have at least one input argument.' )      
   ELSE
   
      M = mxGetM(PRHS(1))
      N = mxGetN(PRHS(1))
      
      IF ((M .NE. 1) .OR. (N .NE. 1)) CALL ProgAbort('Input 1 must be a scalar variable.')
      
      FLAG = INT(mxGetScalar(PRHS(1)))
      
      SELECT CASE (FLAG)

         CASE ( 0 )
         !...........................................................................................
         ! Initialize:
         !          [Err] = WindInflow_mex( 0, FileName, FileType, ReferenceHeight, ReferenceWidth )
         !...........................................................................................
         
            IF (NRHS /= 5) CALL ProgAbort('Wrong number of input arguments for initialization.')
            IF (NLHS /= 1) CALL ProgAbort('Wrong number of output arguments.')

         
               !-------------------------------------------------------------------------------------
               !      *******    The FileName input parameter, argument #2    *******
               !-------------------------------------------------------------------------------------
            ArgNum = 2

            M = mxGetM(PRHS(ArgNum))
            N = mxGetN(PRHS(ArgNum))
            IF( M .NE. 1) CALL ProgAbort('FileName (input #2) must be a row vector.')
      
            ErrStat = mxGetString( PRHS(ArgNum), InitWindData%WindFileName, M*N)
            IF (ErrStat /= 0) CALL ProgAbort('Error getting FileName from input #2.')

               !-------------------------------------------------------------------------------------
               !      *******    The FileType input parameter, argument #3    *******
               !-------------------------------------------------------------------------------------
            ArgNum = 3
            
            M = mxGetM(PRHS(ArgNum))
            N = mxGetN(PRHS(ArgNum))
            
            IF ((M .NE. 1) .OR. (N .NE. 1)) THEN
               CALL ProgAbort('FileType (input #3) must be a scalar variable.')
            ENDIF
            InitWindData%WindFileType = INT( mxGetScalar(PRHS(ArgNum)) )
            
            
               !-------------------------------------------------------------------------------------
               !      *******    The ReferenceHeight input parameter, argument #4    *******
               !-------------------------------------------------------------------------------------
            ArgNum = 4
            
            M = mxGetM(PRHS(ArgNum))
            N = mxGetN(PRHS(ArgNum))
            
            IF ((M .NE. 1) .OR. (N .NE. 1)) THEN
               CALL ProgAbort('ReferenceHeight (input #4) must be a scalar variable.')
            ENDIF
            InitWindData%ReferenceHeight = REAL( mxGetScalar(PRHS(ArgNum)), ReKi)
            
            
               !-------------------------------------------------------------------------------------
               !      *******    The ReferenceWidth input parameter, argument #5    *******
               !-------------------------------------------------------------------------------------
            ArgNum = 5
            
            M = mxGetM(PRHS(ArgNum))
            N = mxGetN(PRHS(ArgNum))
            
            IF ((M .NE. 1) .OR. (N .NE. 1)) THEN
               CALL ProgAbort('ReferenceWidth (input #5) must be a scalar variable.')
            ENDIF
            InitWindData%Width = REAL( mxGetScalar(PRHS(ArgNum)), ReKi)
          
         
               !-------------------------------------------------------------------------------------
               ! Initialize the wind inflow
               !-------------------------------------------------------------------------------------
         
            CALL WindInf_Init( InitWindData, ErrStat )    
         
         
         
         CASE ( 3 )     
         !...........................................................................................
         ! Return wind speed:
         !                         [Err, V] = WindInflow_mex( 3, t, X )
         !...........................................................................................
         
            IF (NRHS /= 3) CALL ProgAbort('Wrong number of input arguments for returning velocities.')
            IF (NLHS /= 2) CALL ProgAbort('Wrong number of output arguments for returning velocities.')
         
         
               !-------------------------------------------------------------------------------------
               !      *******    The TIME (t) input parameter, argument #2    *******
               !-------------------------------------------------------------------------------------

            ArgNum = 2
            
            M = mxGetM(PRHS(ArgNum))
            N = mxGetN(PRHS(ArgNum))
            IF ((M .NE. 1) .AND. (N .NE. 1)) THEN
               CALL ProgAbort('TIME (input #2) must be a 1-D variable array.')
            ELSE
               NT = MAX(M, N)
               N = 1
            ENDIF

!bjj: this will not work with Matlab 2009b!!!!!
            IF ( .NOT. ALLOCATED(TimeAry) ) THEN
               ALLOCATE( TimeAry(NT), STAT=ErrStat )
               IF (ErrStat /= 0) CALL ProgAbort( ' Error allocating the TimeAry array.' )
            END IF

!            Time = REAL( mxGetScalar(PRHS(ArgNum)), ReKi)

            ptr_X = mxGetPr(PRHS(ArgNum))
            CALL mxCopyPtrToReal8( ptr_X,   TimeAry,   NT)
         
               !-------------------------------------------------------------------------------------
               !      *******    The X (position) input parameter, argument #3    *******
               !-------------------------------------------------------------------------------------
               
            ArgNum = 3
                           
            M = mxGetM(PRHS(ArgNum))   ! number of positions
            N = mxGetN(PRHS(ArgNum))   ! 3
            IF (  ( N .EQ. 3 ) ) THEN
                  
               NX = M               
                  
               IF ( NX /= NT .AND. NT /= 1 .AND. NX /= 1 ) THEN
                  CALL ProgAbort( ' Position and time arrays must be the same length (unless one has length 1)' )
               END IF
      
               IF ( .NOT. ALLOCATED(X) ) THEN
                  ALLOCATE( X(NX,N), STAT=ErrStat )
                  IF (ErrStat /= 0) CALL ProgAbort( ' Error allocating the Position array.' )
               END IF

                  
               ptr_X = mxGetPr(PRHS(ArgNum))
               CALL mxCopyPtrToReal8(ptr_X,   X,   NX*N)                                 

                    
               NV = MAX(NT,NX)
               IF ( .NOT. ALLOCATED(Velocity) ) THEN
                  ALLOCATE( Velocity(NV,N), STAT=ErrStat )
                  IF (ErrStat /= 0) CALL ProgAbort( ' Error allocating the Velocity array.' )
               END IF                        
                    
                                                
            ELSE
               CALL ProgAbort('Position (X, input #3) must have 3-columns array.')
            ENDIF
         
                  
               !-------------------------------------------------------------------------------------
               ! Get the wind speed at this time/location
               !-------------------------------------------------------------------------------------         
         
            
            DO ITime = 1,NV
                              
               Time          = TimeAry( MIN(ITime, NT) )   
                           
               InputPosition = X( MIN(ITime,NX), : ) 
               
               MyWindSpeed   = WindInf_GetVelocity( Time, InputPosition, ErrStat )
               
               Velocity(ITime,:) = MyWindSpeed%Velocity(:)               
            END DO
         
               !-------------------------------------------------------------------------------------
               ! Return the velocity
               !-------------------------------------------------------------------------------------         
         
            PLHS(2) = mxCreateDoubleMatrix(NV, N, mxREAL)
            ptr_Y   = mxGetPr(PLHS(2))                               ! get a pointer to the output MATLAB array
	         CALL mxCopyReal8ToPtr( Velocity, ptr_Y, N*NV)            ! Load the WindInflow output into the MATLAB array (PLHS).
         
            DEALLOCATE( X        )
            DEALLOCATE( TimeAry  )
            DEALLOCATE( Velocity )
         
         CASE ( 9 )  
         !...........................................................................................
         ! Terminate:
         !                             [Err] = WindInflow_mex( 9 )
         !...........................................................................................

            IF (NLHS /= 1) CALL ProgAbort('Wrong number of output arguments.')

            CALL WindInf_Terminate( ErrStat )
                                    
            
      END SELECT
      
         !-------------------------------------------------------------------------------------------
         ! Return the error status
         !-------------------------------------------------------------------------------------------        
   
      M       = 1
      N       = 1
      PLHS(1) = mxCreateDoubleMatrix(M, N, mxREAL)
      ptr_Y   = mxGetPr(PLHS(1))                                    ! get a pointer to the output MATLAB array
	   CALL mxCopyReal8ToPtr( REAL(ErrStat, mxDB), ptr_Y, M*N)       ! Load the error status into the MATLAB array (PLHS).
      
      
   END IF
      

END SUBROUTINE MexFunction
