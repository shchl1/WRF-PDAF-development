!$Id: get_obs_f_pdaf.F90 883 2021-11-27 14:16:40Z lnerger $
!BOP
!
! !ROUTINE: get_obs_f_pdaf --- get vector of synthetic observations from PDAF
!
! !INTERFACE:
SUBROUTINE get_obs_f_pdaf(step, dim_obs, observation)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! The routine is called when synthetic observations
! are generated with PDAF. With the call, the user
! is provided with a generated observation vetor. 
! This can then e.g. be written to a file.
!
! The routine is called by all filter processes.
!
! Version for the dummy model with domain 
! decomposition. Here, the state is fully observed.
!
! !REVISION HISTORY:
! 2019-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE obs_SOUND_T_pdafomi, &
       ONLY: file_syntobs_SOUND_T
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_obs              ! Dimension of obs. vector
  REAL(8), INTENT(out)   :: observation(dim_obs) ! Observation vector
  CHARACTER(len=2) :: stepstr          ! String for time step

! !CALLING SEQUENCE:
! Called by: PDAF_gen_obs   (as U_get_obs_f)
!EOP

  IF (step<10) THEN
     WRITE (stepstr, '(i1)') step
  ELSE
     WRITE (stepstr, '(i2)') step
  END IF  

! *************************
! *** store observation ***
! *************************

  IF (filtertype == 100) THEN
    CALL write_syn_obs(step, file_syntobs_SOUND_T, &
                       dim_obs, observation, 1)
  ENDIF

END SUBROUTINE get_obs_f_pdaf

