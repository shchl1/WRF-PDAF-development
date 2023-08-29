!$Id: callback_obs_pdafomi.F90 883 2021-11-27 14:16:40Z lnerger $
!> callback_obs_pdafomi
!!
!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specific routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!!   When adding an observation type, one has to add one module
!!   obs_TYPE_pdafomi (based on the template obs_TYPE_pdafomi_TEMPLATE.F90).
!!   In addition one has to add a call to the different routines include
!!   in this file. It is recommended to keep the order of the calls
!!   consistent over all files. 
!! 
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYPE.
!!
SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

  ! Include functions for different observations
  USE obs_SOUND_U_pdafomi, ONLY: assim_SOUND_U, init_dim_obs_SOUND_U
  USE obs_SOUND_V_pdafomi, ONLY: assim_SOUND_V, init_dim_obs_SOUND_V
  USE obs_SOUND_T_pdafomi, ONLY: assim_SOUND_T, init_dim_obs_SOUND_T
  USE obs_SOUND_PH_pdafomi, ONLY: assim_SOUND_PH, init_dim_obs_SOUND_PH
  USE obs_GNSS_PW_pdafomi, ONLY: assim_GNSS_PW, init_dim_obs_GNSS_PW
  USE obs_GNSS_ZTD_pdafomi, ONLY: assim_GNSS_ZTD, init_dim_obs_GNSS_ZTD
! added by Changliang Shao on 2022.1228
  USE module_domain, ONLY : head_grid, get_ijk_from_grid 

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_SOUND_U ! Observation dimensions
  INTEGER :: dim_obs_SOUND_V ! Observation dimensions
  INTEGER :: dim_obs_SOUND_T ! Observation dimensions
  INTEGER :: dim_obs_SOUND_PH ! Observation dimensions
  INTEGER :: dim_obs_GNSS_PW ! Observation dimensions
  INTEGER :: dim_obs_GNSS_ZTD ! Observation dimensions
  
! added by Changliang Shao on 2022.1228  
  ! INTEGER                   :: ids, ide, jds, jde, kds, kde,    &
                               ! ims, ime, jms, jme, kms, kme,    &
                               ! ips, ipe, jps, jpe, kps, kpe
  ! integer                   :: mm,dimx,dimy,dimz,dimp,     &
                               ! cnt1,cnt2,cnt3,cnt4,cnt5
  ! call get_ijk_from_grid (  head_grid ,                &
                            ! ids, ide, jds, jde, kds, kde,    &
                            ! ims, ime, jms, jme, kms, kme,    &
                            ! ips, ipe, jps, jpe, kps, kpe   ) 

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_SOUND_U = 0
  dim_obs_SOUND_V = 0
  dim_obs_SOUND_T = 0
  dim_obs_SOUND_PH = 0
  dim_obs_GNSS_PW = 0
  dim_obs_GNSS_ZTD = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_SOUND_U) CALL init_dim_obs_SOUND_U(step, dim_obs_SOUND_U)
  IF (assim_SOUND_V) CALL init_dim_obs_SOUND_V(step, dim_obs_SOUND_V)
  IF (assim_SOUND_T) CALL init_dim_obs_SOUND_T(step, dim_obs_SOUND_T)
  IF (assim_SOUND_PH) CALL init_dim_obs_SOUND_PH(step, dim_obs_SOUND_PH)
  IF (assim_GNSS_PW) CALL init_dim_obs_GNSS_PW(step, dim_obs_GNSS_PW)
  IF (assim_GNSS_ZTD) CALL init_dim_obs_GNSS_ZTD(step, dim_obs_GNSS_ZTD)

  dim_obs = dim_obs_SOUND_U + dim_obs_SOUND_V + dim_obs_SOUND_T + &
            dim_obs_SOUND_PH + dim_obs_GNSS_PW + dim_obs_GNSS_ZTD

END SUBROUTINE init_dim_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_OBSTYPE.
!!
SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  ! Include functions for different observations
  USE obs_SOUND_U_pdafomi, ONLY: assim_SOUND_U, obs_op_SOUND_U, init_dim_obs_SOUND_U
  USE obs_SOUND_V_pdafomi, ONLY: assim_SOUND_V, obs_op_SOUND_V, init_dim_obs_SOUND_V
  USE obs_SOUND_T_pdafomi, ONLY: assim_SOUND_T, obs_op_SOUND_T, init_dim_obs_SOUND_T
  USE obs_SOUND_PH_pdafomi, ONLY: assim_SOUND_PH, obs_op_SOUND_PH, init_dim_obs_SOUND_PH
  USE obs_GNSS_PW_pdafomi, ONLY: assim_GNSS_PW, obs_op_GNSS_PW, init_dim_obs_GNSS_PW
  USE obs_GNSS_ZTD_pdafomi, ONLY: assim_GNSS_ZTD, obs_op_GNSS_ZTD, init_dim_obs_GNSS_ZTD

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs              !< Dimension of full observed state
  REAL(8), INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL(8), INTENT(inout) :: ostate(dim_obs)      !< PE-local full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  ! The order of these calls is not relevant as the setup
  ! of the overall observation vector is defined by the
  ! order of the calls in init_dim_obs_pdafomi
  IF (assim_SOUND_U) CALL obs_op_SOUND_U(dim_p, dim_obs, state_p, ostate)
  IF (assim_SOUND_V) CALL obs_op_SOUND_V(dim_p, dim_obs, state_p, ostate)
  IF (assim_SOUND_T) CALL obs_op_SOUND_T(dim_p, dim_obs, state_p, ostate)
  IF (assim_SOUND_PH) CALL obs_op_SOUND_PH(dim_p, dim_obs, state_p, ostate)
  IF (assim_GNSS_PW) CALL obs_op_GNSS_PW(dim_p, dim_obs, state_p, ostate)
  IF (assim_GNSS_ZTD) CALL obs_op_GNSS_ZTD(dim_p, dim_obs, state_p, ostate)
	! print *, 'size(ostate), mean(ostate) = ', size(ostate), sum(ostate)/size(ostate)
  ! ENDIF

END SUBROUTINE obs_op_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

  USE PDAFomi, ONLY: PDAFomi_set_debug_flag
  USE mod_parallel_pdaf, ONLY: mype_filter
  ! Include functions for different observations
  USE obs_SOUND_U_pdafomi, ONLY: assim_SOUND_U, init_dim_obs_l_SOUND_U
  USE obs_SOUND_V_pdafomi, ONLY: assim_SOUND_V, init_dim_obs_l_SOUND_V
  USE obs_SOUND_T_pdafomi, ONLY: assim_SOUND_T, init_dim_obs_l_SOUND_T
  USE obs_SOUND_PH_pdafomi, ONLY: assim_SOUND_PH, init_dim_obs_l_SOUND_PH
  USE obs_GNSS_PW_pdafomi, ONLY: assim_GNSS_PW, init_dim_obs_l_GNSS_PW
  USE obs_GNSS_ZTD_pdafomi, ONLY: assim_GNSS_ZTD, init_dim_obs_l_GNSS_ZTD

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! IF (domain_p==10 .AND. mype_filter==0) THEN
    ! CALL PDAFomi_set_debug_flag(domain_p)
  ! ELSE
    ! CALL PDAFomi_set_debug_flag(0)
  ! ENDIF
  ! Call init_dim_obs_l specific for each observation
  IF (assim_SOUND_U) CALL init_dim_obs_l_SOUND_U(domain_p, step, dim_obs, dim_obs_l)
  IF (assim_SOUND_V) CALL init_dim_obs_l_SOUND_V(domain_p, step, dim_obs, dim_obs_l)
  IF (assim_SOUND_T) CALL init_dim_obs_l_SOUND_T(domain_p, step, dim_obs, dim_obs_l)
  IF (assim_SOUND_PH) CALL init_dim_obs_l_SOUND_PH(domain_p, step, dim_obs, dim_obs_l)
  IF (assim_GNSS_PW) CALL init_dim_obs_l_GNSS_PW(domain_p, step, dim_obs, dim_obs_l)
  IF (assim_GNSS_ZTD) CALL init_dim_obs_l_GNSS_ZTD(domain_p, step, dim_obs, dim_obs_l)

END SUBROUTINE init_dim_obs_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for localize_covar
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type to apply covariance
!! localization in the LEnKF.
!!
SUBROUTINE localize_covar_pdafomi(dim_p, dim_obs, HP_p, HPH)

  ! Include functions for different observations
  USE obs_SOUND_U_pdafomi, ONLY: assim_SOUND_U, localize_covar_SOUND_U
  USE obs_SOUND_V_pdafomi, ONLY: assim_SOUND_V, localize_covar_SOUND_V
  USE obs_SOUND_T_pdafomi, ONLY: assim_SOUND_T, localize_covar_SOUND_T
  USE obs_SOUND_PH_pdafomi, ONLY: assim_SOUND_PH, localize_covar_SOUND_PH
  USE obs_GNSS_PW_pdafomi, ONLY: assim_GNSS_PW, localize_covar_GNSS_PW
  USE obs_GNSS_ZTD_pdafomi, ONLY: assim_GNSS_ZTD, localize_covar_GNSS_ZTD
! added by Changliang Shao on 2022.1228
  USE module_domain, ONLY : head_grid, get_ijk_from_grid
  
  USE mod_parallel_pdaf, &      ! Include rank of filter process
       ONLY: mype_filter
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: deg2rad, n_fields, off_fields, dim_fields, twin_experiment, &
             startn, filtertype, local_range, &
             ids, ide, jds, jde, kds, kde,    &
             ims, ime, jms, jme, kms, kme,    &
             ips, ipe, jps, jpe, kps, kpe

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs               !< number of observations
  REAL(8), INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
  REAL(8), INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH

! *** local variables ***
  INTEGER :: i, j, cnt_p             ! Counters
  INTEGER :: off_nx                  ! Offset in x-direction for parallelization
  REAL(8), ALLOCATABLE :: coords_p(:,:) ! Coordinates of PE-local state vector entries
  integer                   :: k,l,mm,dimx,dimy,dimz,dimp,     &
                               cnt1,cnt2,cnt3,cnt4,cnt5
							   
! **********************
! *** INITIALIZATION ***
! **********************

  ! Initialize coordinate array

  ! One needs to provide the array COORDS_P holding the coordinates of each
  ! element of the process-local state vector. Each column of the array holds
  ! the information for one element. The array can be initialized here using
  ! information on the model grid.

  ALLOCATE(coords_p(3, dim_p))
   
   l = 1   
   DO j = jps, MIN(jpe,jde-1)
      DO k = kps, kpe-1
	     DO i = ips, MIN(ipe,ide-1)
     		 coords_p(1, l) = real(j)
     		 coords_p(2, l) = real(i)
     		 coords_p(3, l) = real(k)
			 l = l +1
		 ENDDO
	  ENDDO
   ENDDO
   
   ! DO j = jps, MIN(jpe,jde-1)
      ! DO k = kps, kpe-1
	     ! DO i = ips, ipe
     		 ! coords_p(1, l) = head_grid%XLONG(i,j) * deg2rad
     		 ! coords_p(2, l) = head_grid%XLAT(i,j) * deg2rad
     		 ! coords_p(3, l) = (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
                              ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81
			 ! l = l +1
		 ! ENDDO
	  ! ENDDO
   ! ENDDO
   
   ! DO j = jps, jpe
      ! DO k = kps, (kpe-1)
	     ! DO i = ips, MIN(ipe,ide-1)
     		 ! coords_p(1, l) = head_grid%XLONG(i,j) * deg2rad
     		 ! coords_p(2, l) = head_grid%XLAT(i,j) * deg2rad
     		 ! coords_p(3, l) = (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
                              ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81
			 ! l = l +1
		 ! ENDDO
	  ! ENDDO
   ! ENDDO 
			 
   ! DO j = jps, MIN(jpe,jde-1)
      ! DO k = kps, kpe
	     ! DO i = ips, MIN(ipe,ide-1)
     		 ! coords_p(1, l) = head_grid%XLONG(i,j) * deg2rad
     		 ! coords_p(2, l) = head_grid%XLAT(i,j) * deg2rad
     		 ! coords_p(3, l) = (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
                              ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81
			 ! l = l +1
		 ! ENDDO
	  ! ENDDO
   ! ENDDO
   
   ! DO j = jps, MIN(jpe,jde-1)
      ! DO k = kps, kpe
	     ! DO i = ips, MIN(ipe,ide-1)
     		 ! coords_p(1, l) = head_grid%XLONG(i,j) * deg2rad
     		 ! coords_p(2, l) = head_grid%XLAT(i,j) * deg2rad
     		 ! coords_p(3, l) = (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
                              ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81
			 ! l = l +1
		 ! ENDDO
	  ! ENDDO
   ! ENDDO
   
   ! DO j = jps, MIN(jpe,jde-1)
      ! DO k = kps, kpe-1
	     ! DO i = ips, MIN(ipe,ide-1)
     		 ! coords_p(1, l) = head_grid%XLONG(i,j) * deg2rad
     		 ! coords_p(2, l) = head_grid%XLAT(i,j) * deg2rad
     		 ! coords_p(3, l) = (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
                              ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81
			 ! l = l +1
		 ! ENDDO
	  ! ENDDO
   ! ENDDO
   
   ! DO mm = 1, 3   
   ! DO j = jps, MIN(jpe,jde-1)
      ! DO k = kps, kpe-1
	     ! DO i = ips, MIN(ipe,ide-1)
     		 ! coords_p(1, l) = head_grid%XLONG(i,j) * deg2rad
     		 ! coords_p(2, l) = head_grid%XLAT(i,j) * deg2rad
     		 ! coords_p(3, l) = (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
                              ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81
			 ! l = l +1
		 ! ENDDO
	  ! ENDDO
   ! ENDDO
   ! ENDDO
   
! *************************************
! *** Apply covariance localization ***
! *************************************

  ! Call localize_covar specific for each observation
  IF (assim_SOUND_U) CALL localize_covar_SOUND_U(dim_p, dim_obs, HP_p, HPH, coords_p)
  IF (assim_SOUND_V) CALL localize_covar_SOUND_V(dim_p, dim_obs, HP_p, HPH, coords_p)
  IF (assim_SOUND_T) CALL localize_covar_SOUND_T(dim_p, dim_obs, HP_p, HPH, coords_p)
  IF (assim_SOUND_PH) CALL localize_covar_SOUND_PH(dim_p, dim_obs, HP_p, HPH, coords_p)
  IF (assim_GNSS_PW) CALL localize_covar_GNSS_PW(dim_p, dim_obs, HP_p, HPH, coords_p)
  IF (assim_GNSS_ZTD) CALL localize_covar_GNSS_ZTD(dim_p, dim_obs, HP_p, HPH, coords_p)


! ****************
! *** Clean up ***
! ****************

  DEALLOCATE(coords_p)

END SUBROUTINE localize_covar_pdafomi
