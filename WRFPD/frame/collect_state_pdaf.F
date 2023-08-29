!$Id: collect_state_pdaf.F90 871 2021-11-22 16:44:34Z lnerger $
!>  Initialize state vector from model fields
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all filters.
!!
!! This subroutine is called during the forecast 
!! phase from PDAF_put_state_X or PDAF_assimilate_X
!! after the propagation of each ensemble member. 
!! The supplied state vector has to be initialized
!! from the model fields (typically via a module). 
!! With parallelization, MPI communication might be 
!! required to initialize state vectors for all 
!! subdomains on the model PEs. 
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! added by Changliang Shao on 2022.1228
  USE module_domain, ONLY : head_grid
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: n_fields, off_fields, dim_fields,  &
	   ids, ide, jds, jde, kds, kde,    &
	   ims, ime, jms, jme, kms, kme,    &
	   ips, ipe, jps, jpe, kps, kpe

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL(8), INTENT(inout) :: state_p(dim_p)  !< local state vector
  
! *** local variables ***
! added by Changliang Shao on 2022.1228  
  integer                   :: i,j,k,l,mm,dimx,dimy,dimz,dimp,     &
                               cnt1,cnt2,cnt3,cnt4,cnt5,    &
	                           its, ite, jts, jte
  integer,dimension(6) :: m

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'collect_state_pdaf: initialization of state vector here!'
 
  ! added by Changliang Shao on 2022.1228
  data m /2,3,4,5,8,9/  !'QV','QC','QR','QI','QS','QG'  
   state_p = 0.
   
   its = ips
   ite = MIN(ipe,ide-1)
   jts = jps
   jte = MIN(jpe,jde-1)
   
   ! state_p(1:cnt1) = &
   ! reshape((head_grid%u_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) + &
   ! head_grid%u_2((ips+1):(MIN(ipe,ide-1)+1),kps:(kpe-1),jps:MIN(jpe,jde-1)))/2., [cnt1])
 
   ! state_p((cnt1+1):(cnt1+cnt2)) = &
   ! reshape((head_grid%v_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) + &
   ! head_grid%v_2(ips:MIN(ipe,ide-1),kps:(kpe-1),(jps+1):(MIN(jpe,jde-1)+1)))/2., [cnt2])   

   ! state_p((cnt1+cnt2+1):(cnt1+cnt2+cnt3)) = &
   ! reshape((head_grid%w_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) + &
   ! head_grid%w_2(ips:MIN(ipe,ide-1),(kps+1):kpe,jps:MIN(jpe,jde-1)))/2., [cnt3])  
    
   ! state_p((cnt1+cnt2+cnt3+1):(cnt1+cnt2+cnt3+cnt4)) = &
   ! reshape((head_grid%ph_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) + &   
   ! head_grid%ph_2(ips:MIN(ipe,ide-1),(kps+1):kpe,jps:MIN(jpe,jde-1)))/2., [cnt4])
      
   l = 1
   do j = jts, jte
     do k = kps, kpe-1
	   do i = its, ite
	     state_p(l+off_fields(1)) = (head_grid%u_2(i,k,j) + head_grid%u_2(i+1,k,j))/2.
	     state_p(l+off_fields(2)) = (head_grid%v_2(i,k,j) + head_grid%v_2(i,k,j+1))/2.
	     state_p(l+off_fields(3)) = (head_grid%w_2(i,k,j) + head_grid%w_2(i,k+1,j))/2.
	     state_p(l+off_fields(4)) = (head_grid%ph_2(i,k,j) + head_grid%ph_2(i,k+1,j))/2.
	     state_p(l+off_fields(5)) = head_grid%t_2(i,k,j)
	     state_p(l+off_fields(6)) = head_grid%moist(i,k,j,m(1))
	     state_p(l+off_fields(7)) = head_grid%moist(i,k,j,m(2))
	     state_p(l+off_fields(8)) = head_grid%moist(i,k,j,m(3))
	     state_p(l+off_fields(9)) = head_grid%moist(i,k,j,m(4))
	     state_p(l+off_fields(10)) = head_grid%moist(i,k,j,m(5))
	     state_p(l+off_fields(11)) = head_grid%moist(i,k,j,m(6))
	     state_p(l+off_fields(12)) = head_grid%p(i,k,j)
	     state_p(l+off_fields(13)) = head_grid%rho(i,k,j)
	     state_p(l+off_fields(14)) = (head_grid%ph_2(i,k+1,j)+head_grid%phb(i,k+1,j)- &
				       head_grid%ph_2(i,k,j)-head_grid%phb(i,k,j))/9.81  ! deta_h unit: m
		 l = l + 1
       enddo
	 enddo
   enddo
   
   l = 1
   do j = jts, jte
     do i = its, ite
	   state_p(l+off_fields(15)) = head_grid%psfc(i,j)
	   l = l + 1
	 enddo
   enddo
                  
   ! ! scheme 1
   ! cnt1 = (ipe-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) ! U
   ! cnt2 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (jpe-jps+1) ! V
   ! cnt3 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps+1) * (MIN(jpe,jde-1)-jps+1) !W
   ! cnt4 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps+1) * (MIN(jpe,jde-1)-jps+1) !PH
   ! cnt5 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) !T   
   
   ! state_p(1:cnt1) = &
   ! reshape(head_grid%u_2(ips:ipe,kps:(kpe-1), &
   ! jps:MIN(jpe,jde-1)), [cnt1])
 
   ! state_p((cnt1+1):(cnt1+cnt2)) = &
   ! reshape(head_grid%v_2(ips:MIN(ipe,ide-1), &
   ! kps:(kpe-1),jps:jpe), [cnt2])

   ! state_p((cnt1+cnt2+1):(cnt1+cnt2+cnt3)) = &
   ! reshape(head_grid%w_2(ips:MIN(ipe,ide-1), &
   ! kps:kpe,jps:MIN(jpe,jde-1)), [cnt3])   
    
   ! state_p((cnt1+cnt2+cnt3+1):(cnt1+cnt2+cnt3+cnt4)) = &
   ! reshape(head_grid%ph_2(ips:MIN(ipe,ide-1), &
   ! kps:kpe,jps:MIN(jpe,jde-1)), [cnt4])
         
   ! state_p((cnt1+cnt2+cnt3+cnt4+1):(cnt1+cnt2+cnt3+cnt4+cnt5)) = &
   ! reshape(head_grid%t_2(ips:MIN(ipe,ide-1), &
   ! kps:(kpe-1),jps:MIN(jpe,jde-1)), [cnt5])
   
   ! DO mm = 1,3
      ! state_p((cnt1+cnt2+cnt3+cnt4+cnt5*mm+1):(cnt1+cnt2+cnt3+cnt4+cnt5+cnt5*mm)) = &
      ! reshape(head_grid%moist(ips:MIN(ipe,ide-1), &
      ! kps:(kpe-1),jps:MIN(jpe,jde-1),m(mm)), [cnt5])
   ! ENDDO   
! if (1==0) then  
   print *, maxval(head_grid%u_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
     minval(head_grid%u_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
     sum(head_grid%u_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
     sum(head_grid%u_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)))/&
	 size(head_grid%u_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)))
   print *, maxval(head_grid%v_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
     minval(head_grid%v_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
	 sum(head_grid%v_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
	 sum(head_grid%v_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)))/&
	 size(head_grid%v_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)))
   print *, maxval(head_grid%w_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
     minval(head_grid%w_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
	 sum(head_grid%w_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
	 sum(head_grid%w_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)))/&
	 size(head_grid%w_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)))
   print *, maxval(head_grid%ph_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
     minval(head_grid%ph_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
	 sum(head_grid%ph_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1))),&
	 sum(head_grid%ph_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)))/&
	 size(head_grid%ph_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)))
   print *, maxval(head_grid%t_2),minval(head_grid%t_2),sum(head_grid%t_2),sum(head_grid%t_2)/size(head_grid%t_2)
   print *, maxval(head_grid%moist(:,:,:,2)),minval(head_grid%moist(:,:,:,2)),sum(head_grid%moist(:,:,:,2))
   print *, maxval(head_grid%moist(:,:,:,3)),minval(head_grid%moist(:,:,:,3)),sum(head_grid%moist(:,:,:,3))
   print *, maxval(head_grid%moist(:,:,:,4)),minval(head_grid%moist(:,:,:,4)),sum(head_grid%moist(:,:,:,4))
   ! print *, maxval(head_grid%moist(:,:,:,5)),minval(head_grid%moist(:,:,:,5)),sum(head_grid%moist(:,:,:,5))
   ! print *, maxval(head_grid%moist(:,:,:,8)),minval(head_grid%moist(:,:,:,8)),sum(head_grid%moist(:,:,:,8))
   ! print *, maxval(head_grid%moist(:,:,:,9)),minval(head_grid%moist(:,:,:,9)),sum(head_grid%moist(:,:,:,9))
! endif

  WRITE (*,*) 'collect_state_pdaf: DONE!'
  
END SUBROUTINE collect_state_pdaf
