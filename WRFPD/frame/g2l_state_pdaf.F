!$Id: g2l_state_pdaf.F90 871 2021-11-22 16:44:34Z lnerger $
!>  Restrict a model state to a local analysis domain
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during the loop over all
!! local analysis domains in PDAF_X_update
!! before the analysis on a single local analysis 
!! domain.  It has to initialize elements of the 
!! state vector for the local analysis domains from
!! the PE-local full state vector.
!!
!! Generic implementation using index vector 
!! ID_LSTATE_IN_PSTATE.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

  USE mod_assimilation, &
       ONLY: id_lstate_in_pstate, n_fields

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  REAL(8), INTENT(in)    :: state_p(dim_p) !< PE-local full state vector 
  REAL(8), INTENT(out)   :: state_l(dim_l) !< State vector on local analysis domain

! *** local variables ***
  INTEGER :: i, l, n                    ! Counter


! *************************************
! *** Initialize local state vector ***
! *************************************

  ! Allocate array
  ! IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
  ! ALLOCATE(id_lstate_in_pstate(dim_l+4)) ! SCHEME 1
  ! Generic initialization using ID_LSTATE_IN_PSTATE set in INIT_DIM_L_PDAF
  ! l = 1
  ! DO n = 1, n_fields
    ! SELECT CASE (n)
        ! CASE (1)
          ! state_l(n) = (state_p(id_lstate_in_pstate(l)) + &
		  ! state_p(id_lstate_in_pstate(l+1)))/2.
          ! l = l+2
        ! CASE (2)
          ! state_l(n) = (state_p(id_lstate_in_pstate(l)) + &
		  ! state_p(id_lstate_in_pstate(l+1)))/2.
          ! l = l+2
        ! CASE (3)
          ! state_l(n) = (state_p(id_lstate_in_pstate(l)) + &
		  ! state_p(id_lstate_in_pstate(l+1)))/2.
          ! l = l+2
        ! CASE (4)
          ! state_l(n) = (state_p(id_lstate_in_pstate(l)) + &
		  ! state_p(id_lstate_in_pstate(l+1)))/2.
          ! l = l+2
        ! CASE DEFAULT
          ! state_l(n) = state_p(id_lstate_in_pstate(l))
          ! l = l+1
    ! END SELECT
  ! ENDDO

  ! ALLOCATE(id_lstate_in_pstate(dim_l)) ! SCHEME 2、3
  state_l(:) = state_p(id_lstate_in_pstate(:))

END SUBROUTINE g2l_state_pdaf
