!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014  National Renewable Energy Laboratory
!
!    This file is part of TurbSim.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
! 
!
!**********************************************************************************************************************************
MODULE TS_TemporalCoherence

   USE  TurbSim_Types
   
   IMPLICIT NONE

CONTAINS


!=======================================================================
!> This subroutine defines the temporal coherence parameters (location
!> and concentration) for the case with no temporal coherence characteristics.
SUBROUTINE TC_None ( p, Rho, Mu )


   IMPLICIT                NONE

         ! Passed variables

   TYPE(Grid_ParameterType),   INTENT(IN   )  :: p                         !<  TurbSim grid parameters
   REAL(ReKi),                 INTENT(OUT  )  :: Rho     (:)               !<  Temporal coherence parameters: concentration    
   REAL(ReKi),                 INTENT(OUT  )  :: Mu      (:)               !<  Temporal coherence parameters: location  

      ! Internal variables

   INTEGER(IntKi)                             :: J                         !< loop counter


   ! Loop over points

DO J = 1,p%NPoints

   Rho(J) = 0.0
   Mu(J)  = 0.0

ENDDO !J


END SUBROUTINE TC_None
!=======================================================================
!> This subroutine defines the temporal coherence parameters (location
!> and concentration) for the empirical temporal coherence model that was
!> fit to data recorded at the National Renewable Energy Laboratory near
!> Boulder, Colorado, USA.
SUBROUTINE TC_NREL ( p, Rho, Mu )


   IMPLICIT                NONE

         ! Passed variables

   TYPE(Grid_ParameterType),   INTENT(IN   )  :: p                         !<  TurbSim grid parameters
   REAL(ReKi),                 INTENT(OUT  )  :: Rho     (:)               !<  Temporal coherence parameters: concentration    
   REAL(ReKi),                 INTENT(OUT  )  :: Mu      (:)               !<  Temporal coherence parameters: location  

      ! Internal variables

   INTEGER(IntKi)                             :: J                         !< loop counter


   ! Loop over points

DO J = 1,p%NPoints

   Rho(J) = 0.2
   Mu(J)  = 0.0

ENDDO !J


END SUBROUTINE TC_NREL
!=======================================================================
END MODULE TS_TemporalCoherence
