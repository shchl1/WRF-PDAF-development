!$Id: init_ens_pdaf.F90 883 2021-11-27 14:16:40Z lnerger $
!> Initialize ensemble
!!
!! User-supplied routine for PDAF.
!!
!! Used for all filters except SEEK, EnKF, LEnKF
!!
!! The routine is called when the filter is
!! initialized in PDAF_filter_init.  
!!
!! This template shows how the ensemble of dim_ens 
!! states can be initializated by exact 2nd order 
!! sampling.  State vectors of the form
!!   $x_i = x + sqrt(FAC) eofV (\Omega C^{-1})^T$
!! fulfill the condition
!!   $P = 1/(FAC)  \sum_{i=1}^{dim\_ens} (x_i - x)(x_i - x)^T$
!! The matrix is initialized in the form of
!! singular values and singular vectors.
!!
!! More information on this initialization variant can be
!! found on the PDAF web site on ensemble initialization.
!!
!! The routine is called by all filter processes and 
!! initializes the ensemble for the local domain.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! *Later revisions - see svn log
!!
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! added by Changliang Shao on 2022.1228
  USE module_domain, ONLY : head_grid, get_ijk_from_grid 
  USE mod_parallel_pdaf, &     ! Parallelization
       ONLY: mype_filter, npes_filter, COMM_filter, MPIerr, MPIstatus
  USE mod_assimilation, &      ! Assimilation variables
       ONLY: dim_state, n_fields, off_fields, dim_fields
  USE mpi
  USE netcdf

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: filtertype              !< Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 !< Size of ensemble
  REAL(8), INTENT(out)   :: state_p(dim_p)          !< PE-local model state
  REAL(8), INTENT(out)   :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for SEIK
  REAL(8), INTENT(out)   :: ens_p(dim_p, dim_ens)   !< PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 !< PDAF status flag


! *** local variables ***
  INTEGER :: ids, ide, jds, jde, kds, kde,    &
             ims, ime, jms, jme, kms, kme,    & 
			 ips, ipe, jps, jpe, kps, kpe
  INTEGER :: i, j, k, row, col, s, n  ! counters
  INTEGER, SAVE :: allocflag = 0      ! Flag for memory counting
  REAL(8), ALLOCATABLE :: ens(:,:)       ! global ensemble
  REAL(8), ALLOCATABLE :: state(:),meanstate(:) ! global state vector
  REAL(4), ALLOCATABLE :: U(:,:,:),V(:,:,:),W(:,:,:),QV(:,:,:),QC(:,:,:),QR(:,:,:)
  REAL(4), ALLOCATABLE :: T(:,:,:),PH(:,:,:)
  REAL(8), ALLOCATABLE :: EU(:,:,:,:),EV(:,:,:,:),EW(:,:,:,:),EPH(:,:,:,:)
  REAL(8), ALLOCATABLE :: ET(:,:,:,:),EQV(:,:,:,:),EQC(:,:,:,:),EQR(:,:,:,:)
  REAL(8), ALLOCATABLE :: eofV(:,:)      ! matrix of eigenvectors V 
  REAL(8), ALLOCATABLE :: svals(:)       ! singular values
  INTEGER :: dimu,dimv,dimw,dimph,dimt,dimqv,dimqc,dimqr
  INTEGER :: rank                     ! Rank of approximated covariance matrix
  REAL(8) :: fac                         ! Square-root of dim_ens-1 or dim_ens
  ! variables and arrays for domain decomposition
  INTEGER :: offset                   ! Row-offset according to domain decomposition
  INTEGER :: domain                   ! domain counter
  REAL(8),ALLOCATABLE :: ens_p_tmp(:,:)  ! Temporary ensemble for some PE-domain
  REAL(8),ALLOCATABLE :: state_p_tmp(:)  ! Temporary state vector for some PE-domain
  REAL(8),ALLOCATABLE :: sendtmp(:), recvtmp(:)  ! Temporary state vector
  INTEGER, ALLOCATABLE :: local_dims(:)   ! PE(/=0)-local state dimension
  INTEGER              :: IERR, dimtmp
  
  CHARACTER(len=120) :: inpath, outpath         ! File paths
  CHARACTER(len=120) :: infile, outfile         ! File names
  CHARACTER(len=150) :: ncfile_in, ncfile_out   ! File name including path
  INTEGER :: ncid_in,ipts,ipte,jpts,jpte        ! NC file IDs
  INTEGER :: id_sigma, id_svec, id_mstate       ! NC dimension and variable IDs
  INTEGER :: dimid_rank, dimid_state, dimid_one ! NC dimension and variable IDs
  INTEGER :: stat(100)                          ! Status for NC operations
  INTEGER :: countv(4), startv(4)               ! Vectors for NC operations
  integer, allocatable :: indx(:,:),cntx(:,:),indy(:)
  
! added by Changliang Shao on 2023.0210  
  integer                   :: mm,dimx,dimy,dimz,dimp,     &
                               cnt1,cnt2,cnt3,cnt4,cnt5
  integer,dimension(6) :: m
  EXTERNAL :: collect_state_pdaf      ! Collect a state vector from model fields

! **********************
! *** INITIALIZATION ***
! **********************

  ! added by Changliang Shao on 2022.1228
  data m /1,2,3,4,7,8/  !'QV','QC','QR','QI','QS','QG'  
  call get_ijk_from_grid (  head_grid ,                &
                            ids, ide, jds, jde, kds, kde,    &
                            ims, ime, jms, jme, kms, kme,    &
                            ips, ipe, jps, jpe, kps, kpe   ) 
                             
   ! cnt1 = (ipe-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) ! U GRID C
   ! cnt2 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (jpe-jps+1) ! V
   ! cnt3 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps+1) * (MIN(jpe,jde-1)-jps+1) !W
   ! cnt4 = cnt3 !PH
   cnt5 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) !T
   cnt1 = cnt5 ! GRID A 
   cnt2 = cnt5
   cnt3 = cnt5
   cnt4 = cnt5
   ! dimu = (ide-ids+1)*(kde-kds)*(jde-jds) ! U GRID C 
   ! dimv = (ide-ids)*(kde-kds)*(jde-jds+1) ! V
   ! dimw = (ide-ids)*(kde-kds+1)*(jde-jds) ! W
   ! dimph = dimw ! PH
   dimt = (ide-ids)*(kde-kds)*(jde-jds)  ! T
   dimu = dimt
   dimv = dimt
   dimw = dimt
   dimph = dimt
   dimqv = dimt
   dimqc = dimt
   dimqr = dimt     ! QV,QC,QR
   
   ! Path to and name of output file holding covariance matrix
   inpath = '/albedo/work/user/chshao001/PDAFTOOLS/'
   infile = 'covarnobubble.nc'  
   ncfile_in = TRIM(inpath)//TRIM(infile)
   WRITE (*,*) 'Read trajectory from file: ',TRIM(ncfile_in)

   WRITE (*, *) 'init_ens_pdaf: generate ensemble of model states'
   WRITE (*, '(9x, a)') '--- generate from model initial state_p'
    
   CALL collect_state_pdaf(dim_p, state_p)
 
  ! DO i=1, dim_ens
    ! ens_p(:,i) = state_p(:)
  ! END DO
  
   ! write (*,*) 'maxval(state_p), minval(state_p), sum(state_p), ubound(state_p), lbound(state_p)= ', &
                ! maxval(state_p), minval(state_p), sum(state_p), ubound(state_p), lbound(state_p)  

! if (0==1) then 
  ! *** Rank of matrix is ensemble size minus one
   rank = dim_ens - 1
   dimtmp = dim_p
   ALLOCATE ( local_dims(npes_filter), STAT=ierr ) 
   IF (ierr.NE.0) THEN
       write (*,*) 'Failed to allocate local_dims(npes_filter). '
   ENDIF
   local_dims = 0
   call MPI_Gather(dimtmp, 1, MPI_INTEGER, local_dims, 1, &
                   MPI_INTEGER, 0, COMM_filter, MPIerr)
   IF (mype_filter == 0) THEN
      if (sum( local_dims ) /= dim_state) write(*,*) 'error, &
      required dim_state= ',dim_state,'actually ',sum( local_dims )
   ENDIF
   
   ALLOCATE(indx(npes_filter, 4), STAT=ierr)
   IF (ierr.NE.0) write (*,*) 'Failed to allocate indx(npes_filter, 4). '
   ALLOCATE(cntx(npes_filter, 5), STAT=ierr)
   IF (ierr.NE.0) write (*,*) 'Failed to allocate cntx(npes_filter, 5). ' 
   indx = 0
   cntx = 0	  
   call MPI_Gather(ips, 1, MPI_INTEGER, indx(1,1), 1, &
                   MPI_INTEGER, 0, COMM_filter, MPIerr)
   call MPI_Gather(ipe, 1, MPI_INTEGER, indx(1,2), 1, &
                   MPI_INTEGER, 0, COMM_filter, MPIerr)
   call MPI_Gather(jps, 1, MPI_INTEGER, indx(1,3), 1, &
                   MPI_INTEGER, 0, COMM_filter, MPIerr)
   call MPI_Gather(jpe, 1, MPI_INTEGER, indx(1,4), 1, &
                   MPI_INTEGER, 0, COMM_filter, MPIerr)
   call MPI_Gather(cnt1, 1, MPI_INTEGER, cntx(1,1), 1, &
                   MPI_INTEGER, 0, COMM_filter, MPIerr)
   call MPI_Gather(cnt2, 1, MPI_INTEGER, cntx(1,2), 1, &
                   MPI_INTEGER, 0, COMM_filter, MPIerr)
   call MPI_Gather(cnt3, 1, MPI_INTEGER, cntx(1,3), 1, &
                   MPI_INTEGER, 0, COMM_filter, MPIerr)
   call MPI_Gather(cnt4, 1, MPI_INTEGER, cntx(1,4), 1, &
                   MPI_INTEGER, 0, COMM_filter, MPIerr)
   call MPI_Gather(cnt5, 1, MPI_INTEGER, cntx(1,5), 1, &
                   MPI_INTEGER, 0, COMM_filter, MPIerr)
				   
  ! *** Generate full ensemble on filter-PE 0 ***
!if (0==1) then
   mype0: IF (mype_filter == 0) THEN
   
      WRITE (*, '(/9x, a)') 'Generate state ensemble from covariance matrix'
      WRITE (*, '(9x, a)') &
           '--- use rank reduction and 2nd order exact sampling (SEIK type)'
      WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
      WRITE (*, '(9x, a, i5)') '--- number of EOFs: ', rank

     ! allocate memory for temporary fields
      ALLOCATE(eofV(dim_state, rank), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate eofV(dim_state, rank). '
      ALLOCATE(svals(rank), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate svals(rank). '

! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************

  ! We show an example here, in which the ensemble is
  ! generated by second-order exact sampling from a 
  ! state estimate together with an error estimate 
  ! based on singular values and vectors. The full ensemble
  ! is first generated on process 0 only. Subsequently the
  ! ensemble sub-states are distributed.

  ! For this initialization, the mean state, the singular 
  ! values, and the singluar vectors have to be initialized
   
     ! Allocate global ensemble and state
      ALLOCATE(ens(dim_state, dim_ens), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate ens(dim_state, dim_ens). '
      ALLOCATE(state(dim_state), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate state(dim_state). '
	  ALLOCATE(meanstate(dim_state), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate meanstate(dim_state). '
      ! ALLOCATE(U(ide-ids+1, kde-kds, jde-jds), STAT=ierr) ! grid c
      ! IF (ierr.NE.0) write (*,*) 'Failed to allocate U(ide-ids+1, kde-kds, jde-jds). '
      ! ALLOCATE(V(ide-ids, kde-kds, jde-jds+1), STAT=ierr)
      ! IF (ierr.NE.0) write (*,*) 'Failed to allocate V(ide-ids, kde-kds, jde-jds+1). '
      ! ALLOCATE(W(ide-ids, kde-kds+1, jde-jds), STAT=ierr)
      ! IF (ierr.NE.0) write (*,*) 'Failed to allocate W(ide-ids, kde-kds+1, jde-jds). '
      ! ALLOCATE(PH(ide-ids, kde-kds+1, jde-jds), STAT=ierr)
      ! IF (ierr.NE.0) write (*,*) 'Failed to allocate PH(ide-ids, kde-kds+1, jde-jds). '
      ALLOCATE(U(ide-ids, kde-kds, jde-jds), STAT=ierr) ! grid a
      IF (ierr.NE.0) write (*,*) 'Failed to allocate U(ide-ids, kde-kds, jde-jds). '
      ALLOCATE(V(ide-ids, kde-kds, jde-jds), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate V(ide-ids, kde-kds, jde-jds). '
      ALLOCATE(W(ide-ids, kde-kds, jde-jds), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate W(ide-ids, kde-kds, jde-jds). '
      ALLOCATE(PH(ide-ids, kde-kds, jde-jds), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate PH(ide-ids, kde-kds, jde-jds). '
      ALLOCATE(T(ide-ids, kde-kds, jde-jds), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate T(ide-ids, kde-kds, jde-jds). '
      ALLOCATE(QV(ide-ids, kde-kds, jde-jds), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate QV(ide-ids, kde-kds, jde-jds). '
      ALLOCATE(QC(ide-ids, kde-kds, jde-jds), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate QC(ide-ids, kde-kds, jde-jds). '
      ALLOCATE(QR(ide-ids, kde-kds, jde-jds), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate QR(ide-ids, kde-kds, jde-jds). '
      ! ALLOCATE(EU(ide-ids+1, kde-kds, jde-jds, dim_ens), STAT=ierr)
      ! IF (ierr.NE.0) write (*,*) 'Failed to allocate EU(ide-ids+1, kde-kds, jde-jds, dim_ens). '
      ! ALLOCATE(EV(ide-ids, kde-kds, jde-jds+1, dim_ens), STAT=ierr)
      ! IF (ierr.NE.0) write (*,*) 'Failed to allocate EV(ide-ids, kde-kds, jde-jds+1, dim_ens). '
      ! ALLOCATE(EW(ide-ids, kde-kds+1, jde-jds, dim_ens), STAT=ierr)
      ! IF (ierr.NE.0) write (*,*) 'Failed to allocate EW(ide-ids, kde-kds+1, jde-jds, dim_ens). '
      ! ALLOCATE(EPH(ide-ids, kde-kds+1, jde-jds, dim_ens), STAT=ierr)
      ! IF (ierr.NE.0) write (*,*) 'Failed to allocate EPH(ide-ids, kde-kds+1, jde-jds, dim_ens). '
	  ALLOCATE(EU(ide-ids, kde-kds, jde-jds, dim_ens), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate EU(ide-ids, kde-kds, jde-jds, dim_ens). '
      ALLOCATE(EV(ide-ids, kde-kds, jde-jds, dim_ens), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate EV(ide-ids, kde-kds, jde-jds, dim_ens). '
      ALLOCATE(EW(ide-ids, kde-kds, jde-jds, dim_ens), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate EW(ide-ids, kde-kds, jde-jds, dim_ens). '
      ALLOCATE(EPH(ide-ids, kde-kds, jde-jds, dim_ens), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate EPH(ide-ids, kde-kds, jde-jds, dim_ens). '
      ALLOCATE(ET(ide-ids, kde-kds, jde-jds, dim_ens), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate ET(ide-ids, kde-kds, jde-jds, dim_ens). '
      ALLOCATE(EQV(ide-ids, kde-kds, jde-jds, dim_ens), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate EQV(ide-ids, kde-kds, jde-jds, dim_ens). '
      ALLOCATE(EQC(ide-ids, kde-kds, jde-jds, dim_ens), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate EQC(ide-ids, kde-kds, jde-jds, dim_ens). '
      ALLOCATE(EQR(ide-ids, kde-kds, jde-jds, dim_ens), STAT=ierr)
      IF (ierr.NE.0) write (*,*) 'Failed to allocate EQR(ide-ids, kde-kds, jde-jds, dim_ens). '
      state = 0.0
      ens = 0.0
	  U = 0.0
	  V = 0.0
	  W = 0.0
	  PH = 0.0
	  T = 0.0
	  QV = 0.0
	  QC = 0.0
	  QR = 0.0
	  EU = 0.0
	  EV = 0.0
	  EW = 0.0
	  EPH = 0.0
	  ET = 0.0
	  EQV = 0.0
	  EQC = 0.0
	  EQR = 0.0 
   
   END IF mype0
   
  ! *** initialize state *** 
  
   PE0_200: IF (mype_filter /= 0) THEN
 
      ! send sub-fields from PEs /=0 grid c
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	      do i = ips,MIN(ipe,ide-1)
	      ! do i = ips,ipe
            CALL MPI_send(head_grid%u_2(i,k,j), 1, &
            MPI_REAL4, 0, 1111, COMM_filter, MPIerr)
          enddo
		enddo
      enddo
	  do j = jps,MIN(jpe,jde-1)
	  ! do j = jps,jpe
	    do k = kps,(kpe-1)
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_send(head_grid%v_2(i,k,j), 1, &
            MPI_REAL4, 0, 1112, COMM_filter, MPIerr)
          enddo
		enddo
      enddo
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	    ! do k = kps,kpe
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_send(head_grid%w_2(i,k,j), 1, &
            MPI_REAL4, 0, 1113, COMM_filter, MPIerr)
          enddo
		enddo
      enddo
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	    ! do k = kps,kpe
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_send(head_grid%ph_2(i,k,j), 1, &
            MPI_REAL4, 0, 1114, COMM_filter, MPIerr)
          enddo
		enddo
      enddo
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_send(head_grid%t_2(i,k,j), 1, &
            MPI_REAL4, 0, 1115, COMM_filter, MPIerr)
          enddo
		enddo
      enddo
	  do mm = 1,3
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_send(head_grid%moist(i,k,j,m(mm)), 1, &
            MPI_REAL4, 0, 1115+mm, COMM_filter, MPIerr)
          enddo
		enddo
      enddo
	  enddo
 
   ELSE PE0_200
      ! receive and assemble state field 	  
      DO n = 1, npes_filter
	  
      ips = indx(n,1)
      ipe = indx(n,2)
      jps = indx(n,3)
      jpe = indx(n,4)
	  
	  if (n == 1) then	  
	  
        ! ! On PE 0 init state directly grid c
        ! U(ips:ipe,kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
	    ! head_grid%u_2(ips:ipe,kps:(kpe-1),jps:MIN(jpe,jde-1)) 
	    ! V(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:jpe) = &
	    ! head_grid%v_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:jpe)
        ! W(ips:MIN(ipe,ide-1),kps:kpe,jps:MIN(jpe,jde-1)) = &
        ! head_grid%w_2(ips:MIN(ipe,ide-1),kps:kpe,jps:MIN(jpe,jde-1))   
        ! PH(ips:MIN(ipe,ide-1),kps:kpe,jps:MIN(jpe,jde-1)) = &
        ! head_grid%ph_2(ips:MIN(ipe,ide-1),kps:kpe,jps:MIN(jpe,jde-1))    
        ! T(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
        ! head_grid%t_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))	  
	    ! QV(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
	    ! head_grid%moist(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),1)  
	    ! QC(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
	    ! head_grid%moist(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),2)  
	    ! QR(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
	    ! head_grid%moist(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),3)
		
        ! On PE 0 init state directly grid a
		U(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
	    head_grid%u_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) 
	    V(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
	    head_grid%v_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))
        W(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
        head_grid%w_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))   
        PH(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
        head_grid%ph_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))    
        T(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
        head_grid%t_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))	  
	    QV(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
	    head_grid%moist(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),1)  
	    QC(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
	    head_grid%moist(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),2)  
	    QR(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) = &
	    head_grid%moist(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),3)
 
      ! Receive part of state field from PEs > 0 into 
      ! correct part of global state 

      else
	  
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	      do i = ips,MIN(ipe,ide-1)
	      ! do i = ips,ipe
            CALL MPI_recv(U(i,k,j), 1, &
            MPI_REAL4, n-1, 1111, COMM_filter, MPIstatus, MPIerr)
          enddo
		enddo
      enddo
	  ! do j = jps,jpe
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_recv(V(i,k,j), 1, &
            MPI_REAL4, n-1, 1112, COMM_filter, MPIstatus, MPIerr)
          enddo
		enddo
      enddo
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	    ! do k = kps,kpe
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_recv(W(i,k,j), 1, &
            MPI_REAL4, n-1, 1113, COMM_filter, MPIstatus, MPIerr)
          enddo
		enddo
      enddo
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	    ! do k = kps,kpe
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_recv(PH(i,k,j), 1, &
            MPI_REAL4, n-1, 1114, COMM_filter, MPIstatus, MPIerr)
          enddo
		enddo
      enddo
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_recv(T(i,k,j), 1, &
            MPI_REAL4, n-1, 1115, COMM_filter, MPIstatus, MPIerr)
          enddo
		enddo
      enddo
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_recv(QV(i,k,j), 1, &
            MPI_REAL4, n-1, 1115+1, COMM_filter, MPIstatus, MPIerr)
          enddo
		enddo
      enddo
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_recv(QC(i,k,j), 1, &
            MPI_REAL4, n-1, 1115+2, COMM_filter, MPIstatus, MPIerr)
          enddo
		enddo
      enddo
	  do j = jps,MIN(jpe,jde-1)
	    do k = kps,(kpe-1)
	      do i = ips,MIN(ipe,ide-1)
            CALL MPI_recv(QR(i,k,j), 1, &
            MPI_REAL4, n-1, 1115+3, COMM_filter, MPIstatus, MPIerr)
          enddo
		enddo
      enddo
	  
	  endif
	  
      END DO
       
	  ! *** initialize state ***
	  ! ! grid c
      ! state(1:dimu) = &
      ! reshape(U(ids:ide,kds:(kde-1),jds:(jde-1)), [dimu])  
      ! state((dimu+1):(dimu+dimv)) = &
      ! reshape(V(ids:(ide-1),kds:(kde-1),jds:jde), [dimv])
      ! state((dimu+dimv+1):(dimu+dimv+dimw)) = &
      ! reshape(W(ids:(ide-1),kds:kde,jds:(jde-1)), [dimw]) 
      ! state((dimu+dimv+dimw+1):(dimu+dimv+dimw+dimph)) = &
      ! reshape(PH(ids:(ide-1),kds:kde,jds:(jde-1)), [dimph])  
      ! state((dimu+dimv+dimw+dimph+1):(dimu+dimv+dimw+dimph+dimt)) = &
      ! reshape(T(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimt]) 
      ! state((dimu+dimv+dimw+dimph+dimt*1+1):(dimu+dimv+dimw+dimph+dimt+dimt*1)) = &
      ! reshape(QV(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimt])
      ! state((dimu+dimv+dimw+dimph+dimt*2+1):(dimu+dimv+dimw+dimph+dimt+dimt*2)) = &
      ! reshape(QC(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimt])
      ! state((dimu+dimv+dimw+dimph+dimt*3+1):(dimu+dimv+dimw+dimph+dimt+dimt*3)) = &
      ! reshape(QR(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimt])
	  
	  ! grid a
	  state(1:dimu) = &
      reshape(U(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimu])  
      state((dimu+1):(dimu+dimv)) = &
      reshape(V(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimv])
      state((dimu+dimv+1):(dimu+dimv+dimw)) = &
      reshape(W(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimw]) 
      state((dimu+dimv+dimw+1):(dimu+dimv+dimw+dimph)) = &
      reshape(PH(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimph])  
      state((dimu+dimv+dimw+dimph+1):(dimu+dimv+dimw+dimph+dimt)) = &
      reshape(T(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimt]) 
      state((dimu+dimv+dimw+dimph+dimt*1+1):(dimu+dimv+dimw+dimph+dimt+dimt*1)) = &
      reshape(QV(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimt])
      state((dimu+dimv+dimw+dimph+dimt*2+1):(dimu+dimv+dimw+dimph+dimt+dimt*2)) = &
      reshape(QC(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimt])
      state((dimu+dimv+dimw+dimph+dimt*3+1):(dimu+dimv+dimw+dimph+dimt+dimt*3)) = &
      reshape(QR(ids:(ide-1),kds:(kde-1),jds:(jde-1)), [dimt])
	  
      print *, 'mean(U)', sum(U)/size(U)
      print *, 'mean(V)', sum(V)/size(V)
      print *, 'mean(W)', sum(W)/size(W)
      print *, 'mean(PH)', sum(PH)/size(PH)
      print *, 'mean(T)', sum(T)/size(T)
      print *, 'mean(QV)', sum(QV)/size(QV)
      print *, 'mean(QC)', sum(QC)/size(QC)
      print *, 'mean(QR)', sum(QR)/size(QR)
	  
   END IF PE0_200

   mype0a: IF (mype_filter == 0) THEN 
   
     ! *** Open file
     s = 1
     stat(s) = NF90_OPEN(TRIM(ncfile_in), NF90_NOWRITE, ncid_in)	 
     ! *** read meanstate
     s = s + 1
     stat(s) = NF90_INQ_VARID(ncid_in, 'meanstate', id_mstate)
     s = s + 1
     stat(s) = NF90_GET_VAR(ncid_in, id_mstate, meanstate, &
               start = (/1/), count = (/dim_state/))
     DO k = 1,  s
         IF (stat(k) /= NF90_NOERR) &
             WRITE(*, *) 'NetCDF error in writing meanstate'
     END DO
     ! *** read singular values
     s = s + 1
     stat(s) = NF90_INQ_VARID(ncid_in, 'sigma', id_sigma)
     s = s + 1
     stat(s) = NF90_GET_VAR(ncid_in, id_sigma, svals, &
               start = (/1/), count = (/rank/))
     DO k = 1,  s
         IF (stat(k) /= NF90_NOERR) &
             WRITE(*, *) 'NetCDF error in writing singular values'
     END DO
     print *, svals   
     ! *** read singular vectors 
     s = 1
     stat(s) = NF90_INQ_VARID(ncid_in, 'u_svd', id_svec)
	 do i = 1, rank
	    s = s + 1
        stat(s) = NF90_GET_VAR(ncid_in, id_svec, eofV(:,i), &
	              start=(/1, i/), count=(/dim_state, 1/))
     enddo
     DO k = 1,  s
         IF (stat(k) /= NF90_NOERR) WRITE(*, *) 'NetCDF error in &
         writing singular vectors.'
     END DO
	 
     ! open (unit = 1111, file = trim(inpath)//trim('eofVt.log'))
     ! write (1111, *) eofV
     ! close(1111)

     ! Close file
     s = 1
     stat(s) = NF90_CLOSE(ncid_in)
     IF (stat(s) /= NF90_NOERR) WRITE(*, *) 'NetCDF error in closing file'

     WRITE (*,'(/1x,a/)') '------- END Reading svals and eofV-------------'
   
! *****************************************************
! *** DECOMPOSE COVARIANCE                          ***
! ***                                               ***
! *** P = eofV U eofV^T                             ***
! ***   = eofV C^(-1)^T Omega^T Omega C^(-1) eofV^T ***
! *** where U^(-1) = C C^T                          ***
! ***                                               ***
! *** Since the matrix is already initialized in    ***
! *** decomposed form we directly have the          ***
! *** inverses of C given by the singular values    ***
! *****************************************************


! *************************************************
! *** Generate ensemble of interpolating states ***
! *************************************************

     ! Very simple method here: We generate the full 
     ! ensemble on the filter PE with rank 0. Afterwards
     ! we distribute sub-states to other filter PEs
     if (dim_ens > 1) then 
       ! Generate ensemble using PDAF sampling routine
	   ! SCHEME 1: using state
       CALL PDAF_SampleEns(dim_state, dim_ens, eofV, svals, state, &
            ens, 1, flag)
       ! Generate ensemble using PDAF sampling routine
	   ! SCHEME 2: using meanstate
       ! CALL PDAF_SampleEns(dim_state, dim_ens, eofV, svals, meanstate, &
            ! ens, 1, flag)
     else
       ! added by Changliang Shao on 2022.1229 for using true state.
       do i = 1, dim_ens
           ens(:, i) = state(:)
       enddo
     endif
   END IF mype0a   

! ****************************
! *** Distribute substates ***
! ****************************
! if (0==1) then
			      
   mype0b: IF (mype_filter == 0) THEN
				   
     ! *** Initialize and send sub-state on PE 0 ***

      EU = reshape(ens(1:dimu,:), &
           ! [ide-ids+1,kde-kds,jde-jds,dim_ens])
           [ide-ids,kde-kds,jde-jds,dim_ens])
      EV = reshape(ens((dimu+1):(dimu+dimv),:), &
           ! [ide-ids,kde-kds,jde-jds+1,dim_ens])
           [ide-ids,kde-kds,jde-jds,dim_ens])
      EW = reshape(ens((dimu+dimv+1):(dimu+dimv+dimw),:), &
           ! [ide-ids,kde-kds+1,jde-jds,dim_ens])
           [ide-ids,kde-kds,jde-jds,dim_ens])
      EPH = reshape(ens((dimu+dimv+dimw+1):(dimu+dimv+dimw+dimph),:), &
           ! [ide-ids,kde-kds+1,jde-jds,dim_ens]) 
           [ide-ids,kde-kds,jde-jds,dim_ens]) 
      ET = reshape(ens((dimu+dimv+dimw+dimph+1):(dimu+dimv+dimw+dimph+dimt),:), &
           [ide-ids,kde-kds,jde-jds,dim_ens]) 
      EQV = reshape(ens((dimu+dimv+dimw+dimph+dimt*1+1):(dimu+dimv+dimw+dimph+dimt+dimt*1),:), &
           [ide-ids,kde-kds,jde-jds,dim_ens]) 
      EQC = reshape(ens((dimu+dimv+dimw+dimph+dimt*2+1):(dimu+dimv+dimw+dimph+dimt+dimt*2),:), &
           [ide-ids,kde-kds,jde-jds,dim_ens]) 
      EQR = reshape(ens((dimu+dimv+dimw+dimph+dimt*3+1):(dimu+dimv+dimw+dimph+dimt+dimt*3),:), &
           [ide-ids,kde-kds,jde-jds,dim_ens]) 
      do n = 1,dim_ens
      ! print *, 'mean(EU),n',n, sum(EU(:,:,:,n))/size(EU(:,:,:,n))
      ! print *, 'mean(EV),n',n, sum(EV(:,:,:,n))/size(EV(:,:,:,n))
      ! print *, 'mean(EW),n',n, sum(EW(:,:,:,n))/size(EW(:,:,:,n))
      ! print *, 'mean(EPH),n',n, sum(EPH(:,:,:,n))/size(EPH(:,:,:,n))
      ! print *, 'mean(ET),n',n, sum(ET(:,:,:,n))/size(ET(:,:,:,n))
      ! print *, 'mean(EQV),n',n, sum(EQV(:,:,:,n))/size(EQV(:,:,:,n))
      ! print *, 'mean(EQC),n',n, sum(EQC(:,:,:,n))/size(EQC(:,:,:,n))
      ! print *, 'mean(EQR),n',n, sum(EQR(:,:,:,n))/size(EQR(:,:,:,n))
      do j = jds,jde-1
	     do k = kds,kde-1
		    do i = ids,ide-1
			   if (EQV(i,k,j,n) < 0.0) EQV(i,k,j,n) = 0.0
			   if (EQC(i,k,j,n) < 0.0) EQC(i,k,j,n) = 0.0
			   if (EQR(i,k,j,n) < 0.0) EQR(i,k,j,n) = 0.0
            enddo
         enddo
      enddo
	  enddo

      !offset = 0
      DO domain = 1, npes_filter
	  
		 ips = indx(domain,1)
		 ipe = indx(domain,2)
		 jps = indx(domain,3)
		 jpe = indx(domain,4)
		 cnt1 = cntx(domain,1)
		 cnt2 = cntx(domain,2)
		 cnt3 = cntx(domain,3)
		 cnt4 = cntx(domain,4)
		 cnt5 = cntx(domain,5)
     
         whichdomain: IF (domain == 1) THEN
           ! Initialize sub-state and sub_ensemble for PE 0
           ! perform reordering of mode matrix for PE 0
			   
            do i = 1, dim_ens
            ens_p(1:cnt1,i) = &
            ! reshape(EU(ips:ipe,kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt1]) 
            reshape(EU(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt1]) 
            ens_p((cnt1+1):(cnt1+cnt2),i) = &
            ! reshape(EV(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:jpe,i), [cnt2])
            reshape(EV(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt2])
            ens_p((cnt1+cnt2+1):(cnt1+cnt2+cnt3),i) = &
            ! reshape(EW(ips:MIN(ipe,ide-1),kps:kpe,jps:MIN(jpe,jde-1),i), [cnt3])  
            reshape(EW(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt3])    
            ens_p((cnt1+cnt2+cnt3+1):(cnt1+cnt2+cnt3+cnt4),i) = &
            ! reshape(EPH(ips:MIN(ipe,ide-1),kps:kpe,jps:MIN(jpe,jde-1),i), [cnt4])  
            reshape(EPH(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt4])      
            ens_p((cnt1+cnt2+cnt3+cnt4+1):(cnt1+cnt2+cnt3+cnt4+cnt5),i) = &
            reshape(ET(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt5])   
            ens_p((cnt1+cnt2+cnt3+cnt4+cnt5*1+1):&
			(cnt1+cnt2+cnt3+cnt4+cnt5+cnt5*1),i) = &
            reshape(EQV(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt5])			
            ens_p((cnt1+cnt2+cnt3+cnt4+cnt5*2+1):&
			(cnt1+cnt2+cnt3+cnt4+cnt5+cnt5*2),i) = &
            reshape(EQC(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt5])			
            ens_p((cnt1+cnt2+cnt3+cnt4+cnt5*3+1):&
			(cnt1+cnt2+cnt3+cnt4+cnt5+cnt5*3),i) = &
            reshape(EQR(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i),[cnt5]) 
            enddo			
           
         ELSE whichdomain
           ! Initialize sub-state and sub_ensemble for other PEs
           ! and send sub-arrays

           ! allocate temporary sub-arrays
            ALLOCATE(ens_p_tmp(local_dims(domain), dim_ens))
            ! ALLOCATE(sendtmp(local_dims(domain)*dim_ens))

           ! perform reordering of mode matrix
           ! perform reordering of state
			
            do i = 1, dim_ens
			ens_p_tmp(1:cnt1,i) = &
            ! reshape(EU(ips:ipe,kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt1]) 
            reshape(EU(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt1]) 
            ens_p_tmp((cnt1+1):(cnt1+cnt2),i) = &
            ! reshape(EV(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:jpe,i), [cnt2])
            reshape(EV(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt2])
            ens_p_tmp((cnt1+cnt2+1):(cnt1+cnt2+cnt3),i) = &
            ! reshape(EW(ips:MIN(ipe,ide-1),kps:kpe,jps:MIN(jpe,jde-1),i), [cnt3])  
            reshape(EW(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt3])   
            ens_p_tmp((cnt1+cnt2+cnt3+1):(cnt1+cnt2+cnt3+cnt4),i) = &
            ! reshape(EPH(ips:MIN(ipe,ide-1),kps:kpe,jps:MIN(jpe,jde-1),i), [cnt4])  
            reshape(EPH(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt4])    
            ens_p_tmp((cnt1+cnt2+cnt3+cnt4+1):(cnt1+cnt2+cnt3+cnt4+cnt5),i) = &
            reshape(ET(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt5])   
            ens_p_tmp((cnt1+cnt2+cnt3+cnt4+cnt5*1+1):&
			(cnt1+cnt2+cnt3+cnt4+cnt5+cnt5*1),i) = &
            reshape(EQV(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt5])			
            ens_p_tmp((cnt1+cnt2+cnt3+cnt4+cnt5*2+1):&
			(cnt1+cnt2+cnt3+cnt4+cnt5+cnt5*2),i) = &
            reshape(EQC(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i), [cnt5])			
            ens_p_tmp((cnt1+cnt2+cnt3+cnt4+cnt5*3+1):&
			(cnt1+cnt2+cnt3+cnt4+cnt5+cnt5*3),i) = &
            reshape(EQR(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),i),[cnt5])
           ! Send sub-arrays
            CALL MPI_send(ens_p_tmp(1,i), local_dims(domain), &      
                 MPI_REAL8, domain-1, 2000+i, COMM_filter, MPIerr)
            enddo
            DEALLOCATE(ens_p_tmp)			
			
			! sendtmp = reshape(ens_p_tmp, [local_dims(domain)*dim_ens])
			
           ! ! Send sub-arrays
            ! CALL MPI_send(sendtmp, local_dims(domain)*dim_ens, &      
                 ! MPI_REAL8, domain-1, 2000, COMM_filter, MPIerr)
 
            ! DEALLOCATE(ens_p_tmp, sendtmp)

         END IF whichdomain

        ! Increment offset
        ! offset = offset + local_dims(domain)

      END DO
	  
   ELSE mype0b
     ! *** Receive substate on filter-PEs with rank > 0 ***
      do j = 1, dim_ens
	    CALL MPI_recv(ens_p(1,j), dim_p, &
             MPI_REAL8, 0, 2000+j, COMM_filter, MPIstatus, MPIerr)
      enddo
	  ! ALLOCATE(recvtmp(dim_p*dim_ens), STAT=ierr) 
      ! IF (ierr.NE.0) THEN
          ! write (*,*) 'Failed to allocate recvtmp(dim_p*dim_ens). '
      ! ENDIF
      ! CALL MPI_recv(recvtmp, dim_p*dim_ens, &
           ! MPI_REAL8, 0, 2000, COMM_filter, MPIstatus, MPIerr)	 
      ! ens_p = reshape(recvtmp, [dim_p, dim_ens])	  
      ! DEALLOCATE(recvtmp)
     
   END IF mype0b

  
! ****************************************************
! *** no pert added by Changliang Shao on 20230207 ***
! ****************************************************
   ! ALLOCATE(indy(8))
   ! indy = [4,6,7,8]
   ! DO i=1, dim_ens
       ! do j = 1, size(indy)
           ! ens_p((1+off_fields(indy(j))):(dim_fields(indy(j))+off_fields(indy(j))), i) = &
           ! state_p((1+off_fields(indy(j))):(dim_fields(indy(j))+off_fields(indy(j)))) 
       ! enddo
   ! END DO
   
! ****************
! *** clean up ***
! ****************
! endif
   DEALLOCATE(local_dims)
   DEALLOCATE(indx,cntx)
   IF (mype_filter == 0) THEN
      DEALLOCATE(svals, eofV)
      DEALLOCATE(ens, state)
      DEALLOCATE(U,V,W,PH,T,QV,QC,QR)
      DEALLOCATE(EU,EV,EW,EPH,ET,EQV,EQC,EQR)
   END IF
! endif
   ! write (*,*) 'maxval(state_p), minval(state_p), sum(state_p), ubound(state_p), lbound(state_p)= ', &
                ! maxval(state_p), minval(state_p), sum(state_p), ubound(state_p), lbound(state_p)  
   print *, 'complete init_ens_pdaf'

END SUBROUTINE init_ens_pdaf
