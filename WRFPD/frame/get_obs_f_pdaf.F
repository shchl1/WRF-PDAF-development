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
  USE obs_SOUND_U_pdafomi, &      ! Variables for observation SOUND_U
       ONLY: assim_SOUND_U, file_syntobs_SOUND_U
  USE obs_SOUND_V_pdafomi, &      ! Variables for observation SOUND_V
       ONLY: assim_SOUND_V, file_syntobs_SOUND_V
  USE obs_SOUND_T_pdafomi, &      ! Variables for observation SOUND_T
       ONLY: assim_SOUND_T, file_syntobs_SOUND_T
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_obs              ! Dimension of obs. vector
  REAL(8), INTENT(out)   :: observation(dim_obs) ! Observation vector
  CHARACTER(len=5) :: stepstr          ! String for time step

! !CALLING SEQUENCE:
! Called by: PDAF_gen_obs   (as U_get_obs_f)
!EOP

  IF (step>=0) THEN
      WRITE (stepstr, '(i5.2)') step
  ELSE
      WRITE (stepstr, '(i5.2)') -step
  END IF

! *************************
! *** store observation ***
! *************************
! Note: The observation generation should always be performed for a single 
! observation type at a time. Thus one generates separate observation files 
! for each observation type.

  IF (filtertype == 100) THEN
    IF (assim_SOUND_U) &
	  CALL write_syn_obs(step, file_syntobs_SOUND_U, dim_obs, observation, 1)
    IF (assim_SOUND_V) &
	  CALL write_syn_obs(step, file_syntobs_SOUND_V, dim_obs, observation, 1)
    IF (assim_SOUND_T) &
	  CALL write_syn_obs(step, file_syntobs_SOUND_T, dim_obs, observation, 1)
  ENDIF

END SUBROUTINE get_obs_f_pdaf

