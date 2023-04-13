!$Id: prepoststep_pdaf.F90 886 2021-11-29 12:56:34Z lnerger $
!>  Used-defined Pre/Poststep routine for PDAF
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all ensemble filters.
!! 
!! The routine is called for global filters (e.g. ESTKF)
!! before the analysis and after the ensemble transformation.
!! For local filters (e.g. LESTKF) the routine is called
!! before and after the loop over all local analysis
!! domains.
!!
!! The routine provides full access to the state 
!! estimate and the state ensemble to the user.
!! Thus, user-controlled pre- and poststep 
!! operations can be performed here. For example 
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analyzed, e.g. by 
!! computing the estimated variances. 
!!
!! If a user considers to perform adjustments to the 
!! estimates (e.g. for balances), this routine is 
!! the right place for it.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code based on offline_1D
!! * Later revisions - see repository log
!!
SUBROUTINE prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE mpi                      ! MPI
  USE netcdf
  USE mod_parallel_pdaf, &     ! Parallelization
       ONLY: mype_filter, npes_filter, COMM_filter, MPIerr, MPIstatus
  USE mod_assimilation, &      ! Assimilation variables
       ONLY: dim_state, n_fields, dim_fields, off_fields, nam_fields, &
	   startn, filtertype
  ! added by Changliang Shao on 2022.1228
  USE module_domain, ONLY : head_grid, get_ijk_from_grid 
  USE obs_SOUND_T_pdafomi, ONLY : assim_SOUND_T, file_syntobs_SOUND_T
  USE obs_SOUND_U_pdafomi, ONLY : assim_SOUND_U, file_syntobs_SOUND_U
  USE obs_SOUND_V_pdafomi, ONLY : assim_SOUND_V, file_syntobs_SOUND_V

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step (negative for call after forecast)
  INTEGER, INTENT(in) :: dim_p       !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     !< Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   !< PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   !< PE-local dimension of observation vector
  REAL(8), INTENT(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
  !< (The array 'state_p' is not generally not initialized in the case of SEIK.
  !< It can be used freely here.)
  REAL(8), INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  REAL(8), INTENT(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
  INTEGER, INTENT(in) :: flag        !< PDAF status flag


! *** local variables ***
  INTEGER :: i, j, member, domain    ! Counters
  LOGICAL, SAVE :: firstio = .TRUE.   ! File output is peformed for first time?
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL(8) :: invdim_ens                  ! Inverse ensemble size
  REAL(8) :: invdim_ensm1                ! Inverse of ensemble size minus 1
  REAL(8), ALLOCATABLE :: statmean_est(:)! estimated RMS error
  REAL(8), ALLOCATABLE :: spread_est(:)  ! estimated RMS error
  REAL(8), ALLOCATABLE :: rmserror_est(:)! estimated RMS error
  REAL(8) :: statmean_tmp, rmserror_tmp, spread_tmp  ! estimated RMS error
  REAL(8), ALLOCATABLE :: statmean_p(:)  ! estimated mean
  REAL(8), ALLOCATABLE :: rmserror_p(:)  ! estimated RMS error
  REAL(8), ALLOCATABLE :: spread_p(:)  ! estimated spread
  REAL(8), ALLOCATABLE :: variance_p(:)  ! model state variances
  REAL(8), ALLOCATABLE :: field(:,:,:)   ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  CHARACTER(len=2) :: stepstr         ! String for time step
  CHARACTER(len=3) :: anastr          ! String for call type (initial, forecast, analysis)
  ! Variables for parallelization - global fields
  INTEGER :: off_p,  offset   ! Row-offset according to domain decomposition
  REAL(8), ALLOCATABLE :: variance(:)    ! local variance
  REAL(8), ALLOCATABLE :: ens(:,:)       ! global ensemble
  REAL(8), ALLOCATABLE :: state(:)       ! global state vector
  REAL(8),ALLOCATABLE :: ens_p_tmp(:,:)  ! Temporary ensemble for some PE-domain
  REAL(8),ALLOCATABLE :: state_p_tmp(:)  ! Temporary state for some PE-domain  
  REAL(8),ALLOCATABLE :: obs_t(:),obs_u(:),obs_v(:),obs_tmp(:)    !< PE-local observation
  INTEGER,ALLOCATABLE :: local_dims(:)   ! PE(/=0)-local state dimension
  INTEGER             :: dimtmp, ierr, dim_fields_tmp
  integer                   :: k,mm,dimx,dimy,dimz,dimp,     &
                               cnt1,cnt2,cnt3,cnt4,cnt5
  integer,dimension(6) :: m
  character(len=100) :: out_pathf
                        !- 查询目录文件夹状态或使用目录文件夹，允许字符串中存在空格 
  character(len=100) :: out_patha
                        !- 查询目录文件夹状态或使用目录文件夹，允许字符串中存在空格 
  external :: read_syn_obs_prepost

! **********************
! *** INITIALIZATION ***
! **********************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'prepoststep_pdaf.F90: functionality of prepoststep '
!  write (*,*) 'prepoststep_pdaf, state_p(58650), ens_p(58650,2), dim_p= ', &
!      state_p(58650), ens_p(58650,2), dim_p  

  ! Allocate fields
  write (*,*) 'dim_state= ', dim_state

  ! Set string for time step
  IF (step>=0) THEN
      WRITE (stepstr, '(i2.2)') step
  ELSE
      WRITE (stepstr, '(i2.2)') -step
  END IF
  out_pathf='../output_fore'//trim(stepstr)//'.nc'  
  out_patha='../output_anal'//trim(stepstr)//'.nc' 
  
  IF (mype_filter == 0) THEN
     IF (firsttime) THEN
        WRITE (*, '(8x, a)') 'Analyze initial state ensemble'
        anastr = 'ini'
        ! call init_file_syn_obs(dim_state, out_pathf, 1) 
        ! call init_file_syn_obs(dim_state, out_patha, 1)  
     ELSE
        IF (step<0) THEN
           WRITE (*, '(8x, a)') 'Analyze and write forecasted state ensemble'
           anastr = 'for'
		   ! call write_syn_obs(step, out_pathf, dim_state/8, observation_f, 0)
        ELSE
           WRITE (*, '(8x, a)') 'Analyze and write assimilated state ensemble'
           anastr = 'ana'
		   ! call write_syn_obs(step, out_patha, dim_state/8, observation_f, 0)
        END IF
     END IF
  END IF
  
! !if(0==1) then  
  ALLOCATE(variance_p(dim_p))
  ALLOCATE(statmean_p(n_fields))  
  ALLOCATE(statmean_est(n_fields))  
  ALLOCATE(spread_p(n_fields))
  ALLOCATE(spread_est(n_fields)) 
  ALLOCATE(rmserror_p(n_fields))
  ALLOCATE(rmserror_est(n_fields))
  ALLOCATE(obs_t(dim_fields(5)))
  ALLOCATE(obs_u(dim_fields(1)))
  ALLOCATE(obs_v(dim_fields(2)))
  
  ! Initialize numbers
  variance_p = 0.0
  statmean_p  = 0.0
  rmserror_p  = 0.0
  spread_p = 0.0
  statmean_est  = 0.0
  rmserror_est  = 0.0
  spread_est = 0.0
  obs_t = 0.0
  obs_u = 0.0
  obs_v = 0.0
  invdim_ens    = 1.0 / REAL(dim_ens)  
  IF (dim_ens > 1) THEN
    invdim_ensm1  = 1.0 / REAL(dim_ens - 1)
  ELSE
    invdim_ensm1  = invdim_ens
  ENDIF
  
  IF (filtertype /= 100) THEN 
    CALL read_syn_obs_prepost(file_syntobs_SOUND_U, &
	                dim_fields(1), obs_u, abs(step), 1-mype_filter)
    CALL read_syn_obs_prepost(file_syntobs_SOUND_V, &
	                dim_fields(2), obs_v, abs(step), 1-mype_filter)
    CALL read_syn_obs_prepost(file_syntobs_SOUND_T, &
	                dim_fields(5), obs_t, abs(step), 1-mype_filter)
  ENDIF
  
  print*, 'sum(obs_t)/size(obs_t)= ', sum(obs_t)/size(obs_t)
  print*, 'sum(obs_u)/size(obs_u)= ', sum(obs_u)/size(obs_u)
  print*, 'sum(obs_v)/size(obs_v)= ', sum(obs_v)/size(obs_v)


! **************************************************************
! *** Perform prepoststep for SEIK with re-inititialization. ***
! *** The state and error information is completely in the   ***
! *** ensemble.                                              ***
! *** Also performed for SEIK without re-init at the initial ***
! *** time.                                                  ***
! **************************************************************

  ! *** Compute mean state
  IF (mype_filter == 0) WRITE (*, '(8x, a)') '--- compute ensemble mean'
!if (0==1) then
  state_p = 0.0
  DO member = 1, dim_ens
     state_p(:) = state_p(:) + ens_p(:, member)
  END DO
  state_p(:) = invdim_ens * state_p(:)
!endif
  ! *** Compute sampled variances ***
  DO member = 1, dim_ens
     variance_p(:) = variance_p(:) + (ens_p(:, member) - state_p(:)) &
                     * (ens_p(:, member) - state_p(:))
  END DO
  IF (dim_ens > 1) then
     variance_p(:) = invdim_ensm1 * variance_p(:)
  ELSE
     variance_p(:) = invdim_ens * variance_p(:)
  ENDIF

! ******************************************************
! *** Assemble global variance vector on filter PE 0 ***
! ******************************************************

! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! Total estimated RMS error
  ! This example is univariate - one should distinguish different fields
  DO i = 1, n_fields
      statmean_p(i) = SUM(state_p((1+off_fields(i)):(dim_fields(i)+off_fields(i))))
      spread_p(i) = SUM(variance_p((1+off_fields(i)):(dim_fields(i)+off_fields(i))))	  
	  statmean_tmp = 0.0
	  spread_tmp = 0.0
	  dim_fields_tmp = 0
      call MPI_ALLreduce(statmean_p(i), statmean_tmp, 1, &
	  MPI_REAL8, MPI_SUM, COMM_filter, MPIerr)
      call MPI_ALLreduce(spread_p(i), spread_tmp, 1, &
	  MPI_REAL8, MPI_SUM, COMM_filter, MPIerr)
      call MPI_ALLreduce(dim_fields(i), dim_fields_tmp, 1, &
	  MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
      statmean_est(i) = statmean_tmp/REAL(dim_fields_tmp,8)
      spread_est(i) = SQRT(spread_tmp/REAL(dim_fields_tmp,8))
	  if (i == 5) then
	    rmserror_tmp = 0.0
        rmserror_p(i) = SUM((state_p((1+off_fields(5)):(dim_fields(5)+off_fields(5))) &
		- obs_t(1:dim_fields(5))) ** 2)
        call MPI_ALLreduce(rmserror_p(i), rmserror_tmp, 1, &
		MPI_REAL8, MPI_SUM, COMM_filter, MPIerr)
		rmserror_est(i) = SQRT(rmserror_tmp/REAL(dim_fields_tmp,8))
	  endif
	  if (i == 1) then
	    rmserror_tmp = 0.0
        rmserror_p(i) = SUM((state_p((1+off_fields(1)):(dim_fields(1)+off_fields(1))) &
		- obs_u(1:dim_fields(1))) ** 2)
        call MPI_ALLreduce(rmserror_p(i), rmserror_tmp, 1, &
		MPI_REAL8, MPI_SUM, COMM_filter, MPIerr)
		rmserror_est(i) = SQRT(rmserror_tmp/REAL(dim_fields_tmp,8))
	  endif
	  if (i == 2) then
	    rmserror_tmp = 0.0
        rmserror_p(i) = SUM((state_p((1+off_fields(2)):(dim_fields(2)+off_fields(2))) &
		- obs_v(1:dim_fields(2))) ** 2)
        call MPI_ALLreduce(rmserror_p(i), rmserror_tmp, 1, &
		MPI_REAL8, MPI_SUM, COMM_filter, MPIerr)
		rmserror_est(i) = SQRT(rmserror_tmp/REAL(dim_fields_tmp,8))
	  endif
  ENDDO 
    
! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  IF (mype_filter == 0) THEN
     write (*, '(a,6x,a)') 'WRF-PDAF', 'MEAN according to sample variance'
     do i = 1, n_fields
	    write (*,'(a,4x,a10,4x,a2,2x,es12.4)') &
              'WRF-PDAF', 'MEAN-'//anastr, trim(nam_fields(i)), statmean_est(i)
     end do
     write (*, '(a,6x,a)') 'WRF-PDAF', 'RMSE according to sample variance'
     do i = 1, n_fields
	    write (*,'(a,4x,a10,4x,a2,2x,es12.4)') &
              'WRF-PDAF', 'RMSE-'//anastr, trim(nam_fields(i)), rmserror_est(i)
     end do
     write (*, '(a,6x,a)') 'WRF-PDAF', 'SPREAD according to sample variance'
     do i = 1, n_fields
	    write (*,'(a,4x,a10,4x,a2,2x,es12.4)') &
              'WRF-PDAF', 'SPREAD-'//anastr, trim(nam_fields(i)), spread_est(i)
     end do
  END IF
  
  ! DEALLOCATE(variance)
  DEALLOCATE(variance_p)
  DEALLOCATE(statmean_p)  
  DEALLOCATE(statmean_est)  
  DEALLOCATE(spread_p)
  DEALLOCATE(spread_est) 
  DEALLOCATE(rmserror_p)
  DEALLOCATE(rmserror_est)
  DEALLOCATE(obs_t)
  DEALLOCATE(obs_u)
  DEALLOCATE(obs_v)
  
! *******************
! *** File output ***
! *******************
! if (1==0) then
  ! ! Here, one could implement file output
    ! notfirst: IF (.not. firsttime) THEN

     ! ALLOCATE(ens(dim_state, dim_ens))
     ! ALLOCATE(state(dim_state))

     ! ! Gather full ensemble on process with rank 0 and write file
     ! mype00b: IF (mype_filter /= 0) THEN

        ! ! *** Send ensemble substates on filter-PEs with rank > 0 ***
        ! do i = 1, dim_ens
            ! CALL MPI_send(ens_p(1, i), dim_p, MPI_DOUBLE_PRECISION, 0, &
            ! mype_filter*100000+i, COMM_filter, MPIerr)
        ! enddo
     ! ELSE mype00b

        ! ! *** Initialize and receive sub-states on PE 0 ***

        ! ! Initialize sub-ensemble for PE 0
        ! ens(1:dim_p, :) = ens_p(1:dim_p, :)

        ! ! Define offset in state vectors
        ! off_p = 0

        ! DO domain = 2, npes_filter
           ! ! Initialize sub-ensemble for other PEs and send sub-arrays

           ! ! Increment offset
           ! off_p = off_p + local_dims(domain-1)

           ! ! Receive sub-arrays
           ! do i = 1, dim_ens
               ! CALL MPI_recv(ens((1+off_p), i), local_dims(domain), MPI_DOUBLE_PRECISION, &
               ! domain-1, (domain-1)*100000+i, COMM_filter, MPIstatus, MPIerr)
           ! enddo
        ! END DO


        ! ! *** Now write analysis ensemble ***

        ! WRITE (*, '(8x, a)') '--- write ensemble and state estimate'
! !if (1==0) then
        ! ! Set string for time step
        ! IF (step>=0) THEN
           ! WRITE (stepstr, '(i2.2)') step
        ! ELSE
           ! WRITE (stepstr, '(i2.2)') -step
        ! END IF

        ! !ALLOCATE(field(dimx, dimy, dimz))

        ! DO member = 1, dim_ens
        
           ! !field(:, :, :) = reshape(ens(:, member), [dimx, dimy, dimz])

           ! WRITE (ensstr, '(i2.2)') member

           ! OPEN(11, file = 'ens_'//TRIM(ensstr)//'_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', status = 'replace')
            
           ! WRITE (11, *) ens(:, member)

           ! CLOSE(11)
        ! END DO
! !endif
     ! END IF mype00b

     ! ! Gather full state vector on process with rank 0 and write to file
     ! mype00c: IF (mype_filter /= 0) THEN

        ! ! *** Send ensemble substates on filter-PEs with rank > 0 ***

        ! CALL MPI_send(state_p, dim_p, &
             ! MPI_DOUBLE_PRECISION, 0, mype_filter, COMM_filter, MPIerr)

     ! ELSE mype00c

        ! ! *** Initialize and receive sub-states on PE 0 ***

        ! ! Initialize sub-state for PE 0
        ! state(1: dim_p) = state_p(1: dim_p)

        ! ! Define offset in state vectors
        ! off_p = 0

        ! DO domain = 2, npes_filter
           ! ! Initialize sub-ensemble for other PEs and send sub-arrays

           ! ! Increment offset
           ! off_p = off_p + local_dims(domain-1)
           ! ! Receive sub-arrays
           ! CALL MPI_recv(state(1+off_p), local_dims(domain), MPI_DOUBLE_PRECISION, &
           ! domain-1, domain-1, COMM_filter, MPIstatus, MPIerr)
        ! END DO
     
        ! ! *** Now write analysis state estimate ***
! !if (1==0) then        
        ! !field = reshape(state, [dimx, dimy, dimz])

        ! OPEN(11, file = 'state_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', status = 'replace')
 
        ! WRITE (11, *) state

        ! CLOSE(11)
! !endif
        ! !DEALLOCATE(field)
     ! END IF mype00c

     ! DEALLOCATE(ens, state)

     ! IF(mype_filter==0) DEALLOCATE ( local_dims )

    ! END IF notfirst
! endif
! ********************
! *** finishing up ***
! ********************
!endif
  firsttime = .FALSE.

END SUBROUTINE prepoststep_pdaf
