!$Id: init_dim_l_pdaf.F90 871 2021-11-22 16:44:34Z lnerger $
!>  Set dimension of local model state
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during analysis step
!! in PDAF_X_update in the loop over all local
!! analysis domains. It has to set the dimension
!! of the local model  state on the current analysis
!! domain.
!!
!! __Revision history:__
!! * 2005-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: coords_l, id_lstate_in_pstate, n_fields, off_fields, deg2rad
  USE module_domain, ONLY : head_grid, get_ijk_from_grid

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

! *** local variables ***
  INTEGER                   :: ids, ide, jds, jde, kds, kde,    &
                               ims, ime, jms, jme, kms, kme,    &
                               ips, ipe, jps, jpe, kps, kpe
  INTEGER                   :: i,j,k,l,m,n,x,cnt1,cnt2,cnt3,cnt4  !< Current time step


! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  ! Template reminder - delete when implementing functionality
  ! WRITE (*,*) 'init_dim_l_pdaf: Set local state dimension here!'

  ! dim_l = n_fields
  dim_l = 1

! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Global coordinates of local analysis domain
  call get_ijk_from_grid (  head_grid ,                &
                            ids, ide, jds, jde, kds, kde,    &
                            ims, ime, jms, jme, kms, kme,    &
                            ips, ipe, jps, jpe, kps, kpe   ) 
							
  ! don't search but calculate i,j,k
  ! n = (i-1)*(maxJ*maxK)+(j-1)*maxK + k
  ! A(S1:E1,S2:E2,S3:E3): A(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))
  ! D3=（E1-S1+1）*(E2-S2+1)* (E3-S3+1)=
  ! (MIN(ipe,ide-1)-ips+1)*(kpe-kps)*(MIN(jpe,jde-1)-jps+1)
  ! D2=（E1-S1+1）*(E2-S2+1) = 
  ! (MIN(ipe,ide-1)-ips+1)*(kpe-kps)
  ! D1=（E1-S1+1）= MIN(ipe,ide-1)-ips+1
  ! d(i,k,j): 
  ! 1+(i-S1)+(k-S2)*D1+(j-S3)*D2 = domain_p 
  ! cnt3 = (MIN(ipe,ide-1)-ips+1)*(kpe-kps)*(MIN(jpe,jde-1)-jps+1)
  ! cnt2 = (MIN(ipe,ide-1)-ips+1)*(kpe-kps)
  ! cnt1 = MIN(ipe,ide-1)-ips+1
  ! IF (1 <= domain_p .and. domain_p <= cnt1) THEN
      ! i = domain_p + ips -1
      ! k = 1
      ! j = 1
  ! ELSEIF (cnt1 < domain_p .and. domain_p <= cnt2) THEN 
      ! k = CEILING(REAL(domain_p)/REAL(cnt1)) + kps - 1
      ! i = domain_p - (k-kps)*cnt1 + ips - 1
      ! j = 1
  ! ELSEIF (cnt2 < domain_p .and. domain_p <= cnt3) THEN  
      ! j = CEILING(REAL(domain_p)/REAL(cnt2)) + jps - 1 
      ! k = CEILING(REAL(domain_p-(j-jps)*cnt2)/REAL(cnt1)) + kps - 1
      ! i = domain_p - (j-jps)*cnt2 - (k-kps)*cnt1 + ips - 1
  ! ELSE
      ! print *, 'domain_p loop error, domain_p= ', domain_p
  ! ENDIF
  ! ! ! for real case below
  ! coords_l(1) = real(i)
  ! coords_l(2) = real(j)
  ! coords_l(3) = real(k)
  
  ! print *, 'domain_p= ', domain_p, 'i,j,k= ', i, j, k
  ! for real case below
  ! coords_l(1) = head_grid%XLONG(i,j) * deg2rad
  ! coords_l(2) = head_grid%XLAT(i,j) * deg2rad
  ! coords_l(3) = (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
  ! ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81
  
  l = 1
  loop_j: DO j = jps, MIN(jpe,jde-1)
     DO k = kps, (kpe-1)
        DO i = ips, MIN(ipe,ide-1)
           IF ( l == domain_p ) THEN
              ! for ideal case on grid point
              coords_l(1) = real(i)
              coords_l(2) = real(j)
              coords_l(3) = real(k)
              EXIT loop_j
           ENDIF		   
           l = l + 1
        ENDDO
     ENDDO
  ENDDO loop_j
    
  ! cnt1 = (ipe-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1)
  ! cnt2 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (jpe-jps+1)
  ! cnt3 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps+1) * (MIN(jpe,jde-1)-jps+1)
  ! cnt4 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1)
    
  ! if (1 <= domain_p <= cnt1) then
  ! l = 1
  ! DO j = jps, MIN(jpe,jde-1)
     ! DO k = kps, (kpe-1)
        ! DO i = ips, jps
           ! IF ( l == domain_p ) THEN
              ! coords_l(1) = head_grid%XLONG_U(i,j) * deg2rad
              ! coords_l(2) = head_grid%XLAT_U(i,j) * deg2rad
              ! coords_l(3) = (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
              ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81
              ! write(*,*) '(head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
              ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81 = ', &
                          ! (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
              ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81
           ! ENDIF
           ! exit 		   
           ! l = l + 1
        ! ENDDO
     ! ENDDO
  ! ENDDO
  ! elseif ((cnt1+1) <= domain_p <= (cnt1+cnt2)) then
  ! l = 1
  ! DO j = jps, jpe
     ! DO k = kps, (kpe-1)
        ! DO i = ips, MIN(ipe,ide-1)
           ! IF ( l == (domain_p-cnt1) ) THEN
              ! coords_l(1) = head_grid%XLONG_V(i,j) * deg2rad
              ! coords_l(2) = head_grid%XLAT_V(i,j) * deg2rad
              ! coords_l(3) = (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
              ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81
           ! ENDIF
           ! exit 		   
           ! l = l + 1
        ! ENDDO
     ! ENDDO
  ! ENDDO
  ! elseif ((cnt1+cnt2+1) <= domain_p <= (cnt1+cnt2+cnt3)) then
  ! l = 1
  ! DO j = jps, MIN(jpe,jde-1)
     ! DO k = kps, kpe
        ! DO i = ips, MIN(ipe,ide-1)
           ! IF ( l == (domain_p-cnt1-cnt2) ) THEN
              ! coords_l(1) = head_grid%XLONG(i,j) * deg2rad
              ! coords_l(2) = head_grid%XLAT(i,j) * deg2rad
              ! coords_l(3) = (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j))/2./9.81
           ! ENDIF
           ! exit 		   
           ! l = l + 1
        ! ENDDO
     ! ENDDO
  ! ENDDO
  ! else
  ! l = 1
  ! DO j = jps, MIN(jpe,jde-1)
     ! DO k = kps, kpe-1
        ! DO i = ips, MIN(ipe,ide-1)
           ! IF ( l == (domain_p-cnt1-cnt2-cnt3) ) THEN
              ! coords_l(1) = head_grid%XLONG(i,j) * deg2rad
              ! coords_l(2) = head_grid%XLAT(i,j) * deg2rad
              ! coords_l(3) = (head_grid%ph_2(i,k,j)+head_grid%phb(i,k,j)+ &
              ! head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j))/2./9.81
           ! ENDIF
           ! exit 		   
           ! l = l + 1
        ! ENDDO
     ! ENDDO
  ! ENDDO
  ! endif

! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! Allocate array
  IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
  ! ALLOCATE(id_lstate_in_pstate(dim_l+4)) ! SCHEME 1

  ! l = 1  
  ! DO n = 1, n_fields
    ! SELECT CASE (n)
        ! CASE (1)
          ! id_lstate_in_pstate(l) = domain_p + off_fields(n)
          ! l = l+1
          ! id_lstate_in_pstate(l) = domain_p + off_fields(n) + 1
          ! l = l+1         
          ! ! if ((i /= MIN(ipe,ide-1)) .or. (ipe == ide)) then
        ! CASE (2)
          ! id_lstate_in_pstate(l) = domain_p + off_fields(n)
          ! l = l+1
          ! id_lstate_in_pstate(l) = domain_p + off_fields(n) &
              ! + (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * 1
          ! l = l+1
          ! ! if ((j /= MIN(jpe,jde-1)) .or. (jpe == jde)) then
        ! CASE (3)
          ! id_lstate_in_pstate(l) = domain_p + off_fields(n)
          ! l = l+1
          ! id_lstate_in_pstate(l) = domain_p + off_fields(n) + (MIN(ipe,ide-1)-ips+1) * 1
          ! l = l+1
        ! CASE (4)
          ! id_lstate_in_pstate(l) = domain_p + off_fields(n)
          ! l = l+1
          ! id_lstate_in_pstate(l) = domain_p + off_fields(n) + (MIN(ipe,ide-1)-ips+1) * 1
          ! l = l+1
        ! CASE DEFAULT
          ! id_lstate_in_pstate(l) = domain_p + off_fields(n)
          ! l = l+1
    ! END SELECT
  ! ENDDO

  ! ALLOCATE(id_lstate_in_pstate(dim_l)) ! SCHEME 2
  ! if (1 <= domain_p <= cnt1) then
  ! id_lstate_in_pstate = off_fields + domain_p
  ! elseif ((cnt1+1) <= domain_p <= (cnt1+cnt2)) then
  ! id_lstate_in_pstate = off_fields + domain_p - cnt1
  ! elseif ((cnt1+cnt2+1) <= domain_p <= (cnt1+cnt2+cnt3)) then
  ! id_lstate_in_pstate = off_fields + domain_p - cnt1 - cnt2
  ! else
  ! id_lstate_in_pstate = off_fields + domain_p - cnt1 - cnt2 - cnt3

  ALLOCATE(id_lstate_in_pstate(dim_l)) ! SCHEME 3
  ! id_lstate_in_pstate = off_fields + domain_p
  id_lstate_in_pstate = off_fields(5) + domain_p
  
END SUBROUTINE init_dim_l_pdaf
