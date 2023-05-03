!$Id: next_observation_pdaf.F90 886 2021-11-29 12:56:34Z lnerger $
!>  Initialize information on next observation
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all filters
!!
!! The subroutine is called before each forecast phase
!! by PDAF_get_state. It has to initialize the number 
!! of time steps until the next available observation 
!! (nsteps) and the current model time (time). In 
!! addition the exit flag (exit) has to be initialized.
!! It indicates if the data assimilation process is 
!! completed such that the ensemble loop in the model 
!! routine can be exited.
!!
!! The routine is called from PDAF_get_state by all processes
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

  USE mod_assimilation, &        ! Assimilation variables
       ONLY: delt_obs, ntimes, startn
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world
       
  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: stepnow  !< Number of the current time step
  INTEGER, INTENT(out) :: nsteps   !< Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   !< Whether to exit forecasting (1 for exit)
  REAL(8), INTENT(out)    :: time     !< Current model (physical) time
  
! *************************************************************
! *** Determine number of time steps until next observation ***
! *************************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'next_observation_pdaf.F90: Set number of time steps in forecast!'

  nsteps = 0 ! delt_obs

! *********************************
! *** Set current physical time ***
! *********************************

  time = 0.0

! *********************
! *** Set exit flag ***
! *********************
! Added by Changliang Shao on 2022.1228
  if (stepnow + delt_obs <= delt_obs * (ntimes + startn)) then
     ! *** During the assimilation process ***
     nsteps = delt_obs   ! This assumes a constant time step interval
	 if(stepnow == 0) nsteps = delt_obs * startn
     doexit = 0          ! Not used in this implementation
     IF (mype_world == 0) THEN
         write (*, *) &
          'WRF-PDAF', stepnow, 'Next observation at time step', stepnow + nsteps
         ! CALL PDAF_print_info(2)
         CALL PDAF_print_info(3)
     ENDIF
  else
     ! *** End of assimilation process ***
     nsteps = 0          ! No more steps
     doexit = 1          ! Not used in this implementation
     IF (mype_world == 0) THEN
         write (*, *) &
          'WRF-PDAF', stepnow, 'No more observations - end assimilation'
         ! CALL PDAF_print_info(2)
         CALL PDAF_print_info(3)
     ENDIF
  end if
  
END SUBROUTINE next_observation_pdaf
