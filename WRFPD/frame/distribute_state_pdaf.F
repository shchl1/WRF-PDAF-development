!$Id: distribute_state_pdaf.F90 871 2021-11-22 16:44:34Z lnerger $
!>  Initialize model fields from state vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all filters.
!!
!! During the forecast phase of the filter this
!! subroutine is called from PDAF_get_state
!! supplying a model state which has to be evolved. 
!! The routine has to initialize the fields of the 
!! model (typically available through a module) from 
!! the state vector of PDAF. With parallelization, 
!! MPI communication might be required to 
!! initialize all subdomains on the model PEs.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE distribute_state_pdaf(dim_p, state_p)


! added by Changliang Shao on 2022.1228
  USE module_domain, ONLY : head_grid, get_ijk_from_grid
  
  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL(8), INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

! *** local variables ***

! added by Changliang Shao on 2022.1228
  INTEGER                   :: ids, ide, jds, jde, kds, kde,    &
                               ims, ime, jms, jme, kms, kme,    &
                               ips, ipe, jps, jpe, kps, kpe
  integer                   :: i,j,k,l,mm,dimx,dimy,dimz,dimp,     &
                               cnt1,cnt2,cnt3,cnt4,cnt5
  integer,dimension(6) :: m

! *************************************************
! *** Initialize model fields from state vector ***
! *** for process-local model domain            ***
!**************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'distribute_state_pdaf: initialization of model fields here!'
  
  data m /1,2,3,4,7,8/  !'QV','QC','QR','QI','QS','QG'
  
  call get_ijk_from_grid (  head_grid ,                &
                            ids, ide, jds, jde, kds, kde,    &
                            ims, ime, jms, jme, kms, kme,    &
                            ips, ipe, jps, jpe, kps, kpe   ) 
                            
   cnt1 = (ipe-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) ! U
   cnt2 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (jpe-jps+1) ! V
   cnt3 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps+1) * (MIN(jpe,jde-1)-jps+1) ! W
   cnt4 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps+1) * (MIN(jpe,jde-1)-jps+1) ! PH
   cnt5 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) ! T
   ! qmoist cnt5 * 6
   dimp = cnt1+cnt2+cnt3+cnt4+cnt5*4   
   if (dimp /= dim_p) write(*,*) 'distribute dimmension of state_p error! ', &
               'required ', dim_p, 'actually ', dimp
               
   head_grid%u_2(ips:ipe, kps:(kpe-1), jps:MIN(jpe,jde-1)) =  &
   reshape(state_p(1:cnt1), &
   [ipe-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])

   head_grid%v_2(ips:MIN(ipe,ide-1), kps:(kpe-1), jps:jpe) =  &
   reshape(state_p((cnt1+1):(cnt1+cnt2)), &
   [MIN(ipe,ide-1)-ips+1,kpe-kps,jpe-jps+1])

   head_grid%w_2(ips:MIN(ipe,ide-1), kps:kpe, jps:MIN(jpe,jde-1)) =  &
   reshape(state_p((cnt1+cnt2+1):(cnt1+cnt2+cnt3)), &
   [MIN(ipe,ide-1)-ips+1,kpe-kps+1,MIN(jpe,jde-1)-jps+1])

   head_grid%ph_2(ips:MIN(ipe,ide-1), kps:kpe, jps:MIN(jpe,jde-1)) =  &
   reshape(state_p((cnt1+cnt2+cnt3+1):(cnt1+cnt2+cnt3+cnt4)), &
   [MIN(ipe,ide-1)-ips+1,kpe-kps+1,MIN(jpe,jde-1)-jps+1])

   ! l = 1
   ! DO j = jps, MIN(jpe,jde-1)
     ! DO k = kps, (kpe-1)
       ! DO i = ips, MIN(ipe,ide-1)
         ! head_grid%t_2(i,k,j) = state_p(l+cnt1+cnt2+cnt3+cnt4)
         ! l = l + 1
       ! ENDDO
     ! ENDDO
   ! ENDDO
   head_grid%t_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) =  &
   reshape(state_p((cnt1+cnt2+cnt3+cnt4+1):(cnt1+cnt2+cnt3+cnt4+cnt5)), &
   [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])

   DO mm = 1,3
      head_grid%moist(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),m(mm)) =  &
      reshape(state_p((cnt1+cnt2+cnt3+cnt4+cnt5*mm+1):(cnt1+cnt2+cnt3+cnt4+cnt5+cnt5*mm)), &
      [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   ENDDO
! if (1==0) then
   print *, maxval(head_grid%u_2),minval(head_grid%u_2),sum(head_grid%u_2),sum(head_grid%u_2)/size(head_grid%u_2)
   print *, maxval(head_grid%v_2),minval(head_grid%v_2),sum(head_grid%v_2),sum(head_grid%v_2)/size(head_grid%v_2)
   print *, maxval(head_grid%w_2),minval(head_grid%w_2),sum(head_grid%w_2),sum(head_grid%w_2)/size(head_grid%w_2)
   print *, maxval(head_grid%ph_2),minval(head_grid%ph_2),sum(head_grid%ph_2),sum(head_grid%ph_2)/size(head_grid%ph_2)
   print *, maxval(head_grid%t_2),minval(head_grid%t_2),sum(head_grid%t_2),sum(head_grid%t_2)/size(head_grid%t_2)
   print *, maxval(head_grid%moist(:,:,:,1)),minval(head_grid%moist(:,:,:,1)),sum(head_grid%moist(:,:,:,1))
   print *, maxval(head_grid%moist(:,:,:,2)),minval(head_grid%moist(:,:,:,2)),sum(head_grid%moist(:,:,:,2))
   print *, maxval(head_grid%moist(:,:,:,3)),minval(head_grid%moist(:,:,:,3)),sum(head_grid%moist(:,:,:,3))
! endif   
  WRITE (*,*) 'distribute_state_pdaf: DONE!'

END SUBROUTINE distribute_state_pdaf
