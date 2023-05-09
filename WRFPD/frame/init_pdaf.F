!$Id: init_pdaf.F90 954 2021-12-16 15:15:09Z lnerger $
!>  Interface routine to call initialization of PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! to perform the internal initialization of PDAF.
!!
!! This variant is for the online mode of PDAF.
!!
!! This routine is generic. However, it assumes a constant observation
!! error (rms_obs_A, etc.). Further, with parallelization the local state
!! dimension dim_state_p is used.
!!
!! __Revision history:__
!! * 2008-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf()

  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAF_init, PDAF_get_state
  USE PDAFomi_obs_f, &
       ONLY: PDAFomi_set_domain_limits
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world, n_modeltasks, task_id, mype_filter, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: dim_state_p, dim_state, screen, filtertype, subtype, dim_ens, &
       incremental, covartype, type_forget, forget, filename, &
       rank_analysis_enkf, locweight, local_range, srange, &
       type_trans, type_sqrt, delt_obs, startn, ntimes, &
       n_fields, off_fields, dim_fields, nam_fields, deg2rad, &
	   ids, ide, jds, jde, kds, kde,    &
	   ims, ime, jms, jme, kms, kme,    &
	   ips, ipe, jps, jpe, kps, kpe
  USE obs_SOUND_U_pdafomi, &      ! Variables for observation SOUND_U
       ONLY: assim_SOUND_U, rms_obs_SOUND_U, file_syntobs_SOUND_U
  USE obs_SOUND_V_pdafomi, &      ! Variables for observation SOUND_V
       ONLY: assim_SOUND_V, rms_obs_SOUND_V, file_syntobs_SOUND_V
  USE obs_SOUND_T_pdafomi, &      ! Variables for observation SOUND_T
       ONLY: assim_SOUND_T, rms_obs_SOUND_T, file_syntobs_SOUND_T
  USE obs_GNSS_PW_pdafomi, &      ! Variables for observation GNSS_PW
       ONLY: assim_GNSS_PW, rms_obs_GNSS_PW, file_syntobs_GNSS_PW
  USE obs_GNSS_ZTD_pdafomi, &      ! Variables for observation GNSS_ZTD
       ONLY: assim_GNSS_ZTD, rms_obs_GNSS_ZTD, file_syntobs_GNSS_ZTD
!   USE mod_model, &                ! Model variables
!        ONLY: nx, ny
! added by Changliang Shao on 2022.1228
  USE module_domain, ONLY : head_grid, get_ijk_from_grid 
! if nested define as below. id is the nest number
     ! INTEGER, INTENT(IN) :: id
     ! TYPE(domain), POINTER :: grid
     ! CALL find_grid_by_id( id, head_grid, grid )

  IMPLICIT NONE
  
! *** Local variables ***
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL(8) :: filter_param_r(3) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: doexit, steps     ! Not used in this implementation
  REAL(8) :: timenow           ! Not used in this implementation
  CHARACTER(len=2) :: filestr
  
! added by Changliang Shao on 2022.1228
  REAL(8) :: lim_coords(2,2)  !< coordinate array (1: longitude, 2: latitude)
                              !< geographic range: longitude (-pi, pi), latitude (-pi/2, pi/2)
                              !< Cartesian range: (x) coordinate grows from left to right; (y) from bottom to top  
  INTEGER                   :: i,j,k,mm,dimx,dimy,dimz,dimp,    &
                               cnt1,cnt2,cnt3,cnt4,cnt5,ierr
  integer,dimension(6) :: m
! *** External subroutines ***
  EXTERNAL :: init_ens_pdaf            ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time, 
                                       ! and dimension of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       prepoststep_pdaf, &             ! User supplied pre/poststep routine
	   init_file_syn_obs
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF - ONLINE MODE'
     WRITE (*,*) 'init_pdaf.F90: Initialize state dimension here!'
  END IF

  write (*,*) 'head_grid%id= ', head_grid%id

  ! *** Define state dimension ***
  ! added by Changliang Shao on 2022.1228
  call get_ijk_from_grid (  head_grid ,                &
                            ids, ide, jds, jde, kds, kde,    &
                            ims, ime, jms, jme, kms, kme,    &
                            ips, ipe, jps, jpe, kps, kpe   ) 
  write (*,*) 'ids, ide, jds, jde, kds, kde= ', ids, ide, jds, jde, kds, kde
  write (*,*) 'ims, ime, jms, jme, kms, kme= ', ims, ime, jms, jme, kms, kme
  write (*,*) 'ips, ipe, jps, jpe, kps, kpe= ', ips, ipe, jps, jpe, kps, kpe
  ! ! scheme 1 grid c
  ! dim_state = (ide-ids+1)*(kde-kds)*(jde-jds) + & ! U
              ! (ide-ids)*(kde-kds)*(jde-jds+1) + & ! V
              ! (ide-ids)*(kde-kds+1)*(jde-jds) + & ! W
              ! (ide-ids)*(kde-kds+1)*(jde-jds) + & ! PH
              ! (ide-ids)*(kde-kds)*(jde-jds)*4     ! T,QV,QC,QR,QI,QS,QG
  ! cnt1 = (ipe-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) ! U
  ! cnt2 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (jpe-jps+1) ! V
  ! cnt3 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps+1) * (MIN(jpe,jde-1)-jps+1) !W
  ! cnt4 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps+1) * (MIN(jpe,jde-1)-jps+1) !PH
  ! cnt5 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) !T
  
  ! scheme 2 grid a
  dim_state = (ide-ids)*(kde-kds)*(jde-jds)*8     ! U,V,W,PH,T,QV,QC,QR,QI,QS,QG
  cnt5 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) !T
  cnt1 = cnt5 ! U
  cnt2 = cnt5 ! V
  cnt3 = cnt5 ! W
  cnt4 = cnt5 ! PH
  ! cnt5 * 6 'QV','QC','QR','QI','QS','QG'
  n_fields = 8
  allocate(off_fields(n_fields))
  allocate(dim_fields(n_fields))
  allocate(nam_fields(n_fields))
  dim_fields = (/cnt1,cnt2,cnt3,cnt4,cnt5,cnt5,cnt5,cnt5/)
  off_fields = (/0,cnt1,cnt1+cnt2,cnt1+cnt2+cnt3,cnt1+cnt2+cnt3+cnt4,    &
               cnt1+cnt2+cnt3+cnt4+cnt5*1,cnt1+cnt2+cnt3+cnt4+cnt5*2,    &
               cnt1+cnt2+cnt3+cnt4+cnt5*3/)
  dim_state_p = cnt1+cnt2+cnt3+cnt4+cnt5*4
  nam_fields = ['U','V','W','PH','T','QV','QC','QR']  
  data m /2,3,4,5,8,9/  !'QV','QC','QR','QI','QS','QG'
  

! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

! *** IO options ***
  screen      = 2  ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
  filtertype = 7    ! Type of filter
                    !   (1) SEIK
                    !   (2) EnKF
                    !   (3) LSEIK
                    !   (4) ETKF
                    !   (5) LETKF
                    !   (6) ESTKF
                    !   (7) LESTKF
                    !   (8) localized EnKF
                    !   (9) NETF
                    !  (10) LNETF
                    !  (12) PF
                    !  (100) GENOBS
  dim_ens = n_modeltasks   ! Size of ensemble for all ensemble filters
  subtype = 0       ! subtype of filter: 
                    !   SEIK:
                    !     (0) mean forecast; new formulation
                    !     (1) mean forecast; old formulation
                    !     (2) fixed error space basis
                    !     (3) fixed state covariance matrix
                    !     (4) SEIK with ensemble transformation
                    !   EnKF:
                    !     (0) analysis for large observation dimension
                    !     (1) analysis for small observation dimension
                    !   LSEIK:
                    !     (0) mean forecast;
                    !     (2) fixed error space basis
                    !     (3) fixed state covariance matrix
                    !     (4) LSEIK with ensemble transformation
                    !   ETKF:
                    !     (0) ETKF using T-matrix like SEIK
                    !     (1) ETKF following Hunt et al. (2007)
                    !       There are no fixed basis/covariance cases, as
                    !       these are equivalent to SEIK subtypes 2/3
                    !   LETKF:
                    !     (0) LETKF using T-matrix like SEIK
                    !     (1) LETKF following Hunt et al. (2007)
                    !       There are no fixed basis/covariance cases, as
                    !       these are equivalent to LSEIK subtypes 2/3
                    !   ESTKF:
                    !     (0) Standard form of ESTKF
                    !     (2) fixed ensemble perturbations
                    !     (3) fixed state covariance matrix
                    !   LESTKF:
                    !     (0) Standard form of LESTKF
                    !     (2) fixed ensemble perturbations
                    !     (3) fixed state covariance matrix
                    !   NETF:
                    !     (0) Standard form of NETF
                    !   LNETF:
                    !     (0) Standard form of LNETF
                    !   PF:
                    !     (0) Standard form of PF
  type_trans = 0    ! Type of ensemble transformation
                    !   SEIK/LSEIK and ESTKF/LESTKF:
                    !     (0) use deterministic omega
                    !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   ETKF/LETKF:
                    !     (0) use deterministic symmetric transformation
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   NETF/LNETF:
                    !     (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                    !     (1) use identity transformation
  forget  = 0.9     ! Forgetting factor
  type_forget = 0   ! Type of forgetting factor 
                    ! SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                    !   (0) fixed
                    !   (1) global adaptive
                    !   (2) local adaptive for LSEIK/LETKF/LESTKF
                    ! NETF/LNETF/PF
                    !   (0) apply inflation on forecast ensemble
                    !   (2) apply inflation on analysis ensemble
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
  covartype = 1     ! Definition of factor in covar. matrix used in SEIK
                    !   (0) for dim_ens^-1 (old SEIK)
                    !   (1) for (dim_ens-1)^-1 (real ensemble covariance matrix)
                    !   This parameter has also to be set internally in PDAF_init.
  rank_analysis_enkf = 0   ! rank to be considered for inversion of HPH
                    ! in analysis of EnKF; (0) for analysis w/o eigendecomposition


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Forecast length (time interval between analysis steps) ***
  delt_obs = 60      ! This should be set according to the data availability
  ! Added by Changliang Shao on 2022.1228
  ntimes = 0        ! times for assim. 0 for first step to pert. maxtime 144
  startn = 0        ! times for start assim. 0 for initialization to generate ensemble
                    ! perturbation and ensemble forecast at least 1 step (from 1).
					! no perturbation for 0. 
					! assimilation start from 1, no assimilation for 0

! *** Which observation type to assimilate
  assim_SOUND_T = .false.
  assim_SOUND_U = .false.
  assim_SOUND_V = .false.
  assim_GNSS_PW = .false.
  assim_GNSS_ZTD = .false.

! *** specifications for observations ***
  rms_obs_SOUND_T = 0.1    ! Observation error standard deviation
  rms_obs_SOUND_U = 0.1    ! Observation error standard deviation
  rms_obs_SOUND_V = 0.1    ! Observation error standard deviation
  rms_obs_GNSS_PW = 0.1    ! Observation error standard deviation
  rms_obs_GNSS_ZTD = 0.1    ! Observation error standard deviation

! *** Localization settings
  locweight = 0     ! Type of localizating weighting # locweight = 2 
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  local_range = 0.001    ! Range in grid points for observation domain in local filters
                    !   local_range = 3.0
  srange = local_range  ! Support range for 5th-order polynomial
                    ! or range for 1/e for exponential weighting

! *** File names
  filename = 'output.dat'  
  
  IF (mype_filter<10) THEN
     WRITE (filestr, '(i1)') mype_filter
  ELSE
     WRITE (filestr, '(i2)') mype_filter
  END IF  
  file_syntobs_SOUND_T = '../obsT'//trim(filestr)//'.nc'
  file_syntobs_SOUND_U = '../obsU'//trim(filestr)//'.nc'
  file_syntobs_SOUND_V = '../obsV'//trim(filestr)//'.nc'
  file_syntobs_GNSS_PW = '../obsPW'//trim(filestr)//'.nc'
  file_syntobs_GNSS_ZTD = '../obsZTD'//trim(filestr)//'.nc'

! ***********************************
! *** Some optional functionality ***
! ***********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

  call init_pdaf_parse()


! *** Initial Screen output ***
! *** This is optional      ***

  IF (mype_world == 0) call init_pdaf_info()


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** For all filters, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

  ! Here we only specify the required set of parameters
  ! For more possible settings see the documentation
  filter_param_i(1) = dim_state_p ! State dimension
  filter_param_i(2) = dim_ens     ! Size of ensemble
  filter_param_r(1) = forget      ! Forgetting factor
  
  CALL PDAF_init(filtertype, subtype, 0, &
       filter_param_i, 2,&
       filter_param_r, 1, &
       COMM_model, COMM_filter, COMM_couple, &
       task_id, n_modeltasks, filterpe, init_ens_pdaf, &
       screen, status_pdaf)

! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF


! **********************************
! *** Prepare ensemble forecasts ***
! **********************************

  CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
       distribute_state_pdaf, prepoststep_pdaf, status_pdaf)
	   
! Initialize file for synthetic observations
  IF (filtertype == 100) THEN
      IF (assim_SOUND_U) CALL init_file_syn_obs(dim_fields(1), file_syntobs_SOUND_U, 1)
      IF (assim_SOUND_V) CALL init_file_syn_obs(dim_fields(2), file_syntobs_SOUND_V, 1)
      IF (assim_SOUND_T) CALL init_file_syn_obs(dim_fields(5), file_syntobs_SOUND_T, 1)
      IF (assim_GNSS_PW) CALL init_file_syn_obs(dim_fields(5)/(kpe-kps), &
	                                            file_syntobs_GNSS_PW, 1)
      IF (assim_GNSS_ZTD) CALL init_file_syn_obs(dim_fields(5)/(kpe-kps), &
	                                             file_syntobs_GNSS_ZTD, 1)
  END IF

! ************************************************************************
! *** Set domain coordinate limits (for use with OMI's use_global_obs) ***
! ************************************************************************  
  ! added by Changliang Shao on 2023.0311  
             ! (2,1)
   ! (:,1)+---------+
        ! |         |           - (1,1) longitude of the western edge of the domain (or x-coordinate of the left edge)
        ! |         |           - (1,2) longitude of the eastern edge of the domain (or x-coordinate of the right edge)
  ! (1,1) |         | (1,2)     - (2,1) latitude of the northern edge of the domain (or y-coordinate of the upper edge)
        ! |         |           - (2,2) latitude of the southern edge of the domain (or y-coordinate of the lower edge)
        ! |         |
        ! +---------+(:,2)
           ! (2,2)
  
  lim_coords(1,1) = real(ips)*15.
  lim_coords(1,2) = real(MIN(ipe,ide-1))*15. !
  lim_coords(2,1) = real(MIN(jpe,jde-1))*15. !
  lim_coords(2,2) = real(jps)*15.
  
  CALL PDAFomi_set_domain_limits(lim_coords)
  
  write (*,*) 'init_pdaf, DONE!'

END SUBROUTINE init_pdaf
