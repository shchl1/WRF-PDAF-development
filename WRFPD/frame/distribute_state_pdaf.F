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
  USE module_domain, ONLY : head_grid
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: ids, ide, jds, jde, kds, kde,    &
	   ims, ime, jms, jme, kms, kme,    &
	   ips, ipe, jps, jpe, kps, kpe
  
  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL(8), INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

! *** local variables ***

! added by Changliang Shao on 2022.1228
  integer                   :: i,j,k,l,mm,dimx,dimy,dimz,dimp,     &
                               cnt1,cnt2,cnt3,cnt4,cnt5
  integer,dimension(6) :: m
  REAL(8), ALLOCATABLE :: state_tmp(:,:,:),state_tmp1(:,:,:), &
                          state_bak(:,:,:),state_bak1(:,:,:)  ! temp state
  REAL(8) :: state_p_bak(dim_p), state_p_amb(dim_p)

! *************************************************
! *** Initialize model fields from state vector ***
! *** for process-local model domain            ***
!**************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'distribute_state_pdaf: initialization of model fields here!'
  
  data m /1,2,3,4,7,8/  !'QV','QC','QR','QI','QS','QG'  
				
   ! scheme 2   			
   cnt5 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) ! T
   cnt1=cnt5
   cnt2=cnt5
   cnt3=cnt5
   cnt4=cnt5
   ! qmoist cnt5 * 4
   dimp = cnt1+cnt2+cnt3+cnt4+cnt5*4   
   if (dimp /= dim_p) write(*,*) 'distribute dimmension of state_p error! ', &
               'required ', dim_p, 'actually ', dimp
			      
   state_p_bak = real(0,8)
   state_p_amb = real(0,8)
   CALL collect_state_pdaf(dim_p, state_p_bak)
   
   state_p_amb = state_p - state_p_bak
   
   ALLOCATE(state_bak(ips-1:MIN(ipe,ide-1)+1, kps-1:kpe, jps-1:MIN(jpe,jde-1)+1))
   ALLOCATE(state_bak1(ips-1:MIN(ipe,ide-1)+1, kps-1:kpe, jps-1:MIN(jpe,jde-1)+1))
   ALLOCATE(state_tmp(ips-1:MIN(ipe,ide-1)+1, kps:kpe, jps-1:MIN(jpe,jde-1)+1))
   ALLOCATE(state_tmp1(ips-1:MIN(ipe,ide-1)+1, kps:kpe, jps-1:MIN(jpe,jde-1)+1))
   
   state_bak = real(0,8)
   state_tmp = real(0,8)
   state_bak(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     reshape(state_p_bak(1:cnt1), &
     [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     reshape(state_p_amb(1:cnt1), &
     [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   state_tmp(ips-1, kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     state_tmp(ips, kps:kpe-1, jps:MIN(jpe,jde-1))*real(2,8) - &
     state_tmp(ips+1, kps:kpe-1, jps:MIN(jpe,jde-1)) 
   state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) = &
   real(0.5,8)*(state_tmp(ips-1:MIN(ipe,ide-1)-1, kps:kpe-1, jps:MIN(jpe,jde-1)) + &
   state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)))
   do j = jps, MIN(jpe,jde-1)
     do k = kps, kpe-1
	   do i = MIN(ipe,ide-1), ips, -1
	     head_grid%u_2(i,k,j) = real(2,8)*state_bak(i,k,j) - head_grid%u_2(i+1,k,j) + &
		                        state_tmp(i,k,j)
       enddo
     enddo
   enddo
	    
   state_bak = real(0,8)
   state_tmp = real(0,8)
   state_bak(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     reshape(state_p_bak(cnt1+1:cnt1+cnt2), &
     [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     reshape(state_p_amb(cnt1+1:cnt1+cnt2), &
     [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1]) 
   state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps-1) =  &
     state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps)*real(2,8) - &
     state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps+1) 
   state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) = &
   real(0.5,8)*(state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps-1:MIN(jpe,jde-1)-1) + &
   state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)))
   do j = MIN(jpe,jde-1), jps, -1
     do k = kps, kpe-1
	   do i = ips, MIN(ipe,ide-1)
	     head_grid%v_2(i,k,j) = real(2,8)*state_bak(i,k,j) - head_grid%v_2(i,k,j+1) + &
		                        state_tmp(i,k,j)
       enddo
     enddo
   enddo
   	    
   ! state_bak = real(0,8)
   ! state_tmp = real(0,8)
   ! state_bak(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_bak(cnt1+cnt2+1:cnt1+cnt2+cnt3), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_amb(cnt1+cnt2+1:cnt1+cnt2+cnt3), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   ! state_tmp(ips:MIN(ipe,ide-1), kps-1, jps:MIN(jpe,jde-1)) =  &
     ! state_tmp(ips:MIN(ipe,ide-1), kps, jps:MIN(jpe,jde-1))*real(2,8) - &
     ! state_tmp(ips:MIN(ipe,ide-1), kps+1, jps:MIN(jpe,jde-1)) 
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) = &
   ! real(0.5,8)*(state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) + &
   ! state_tmp(ips:MIN(ipe,ide-1), kps-1:kpe-2, jps:MIN(jpe,jde-1)))
   
   ! state_bak1 = real(0,8)
   ! state_tmp1 = real(0,8)
   ! state_bak1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_bak(cnt1+cnt2+cnt3+1:cnt1+cnt2+cnt3+cnt4), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   ! state_tmp1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_amb(cnt1+cnt2+cnt3+1:cnt1+cnt2+cnt3+cnt4), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   ! state_tmp1(ips:MIN(ipe,ide-1), kps-1, jps:MIN(jpe,jde-1)) =  &
     ! state_tmp1(ips:MIN(ipe,ide-1), kps, jps:MIN(jpe,jde-1))*real(2,8) - &
     ! state_tmp1(ips:MIN(ipe,ide-1), kps+1, jps:MIN(jpe,jde-1)) 
   ! state_tmp1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) = &
   ! real(0.5,8)*(state_tmp1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) + &
   ! state_tmp1(ips:MIN(ipe,ide-1), kps-1:kpe-2, jps:MIN(jpe,jde-1)))
   
   ! do j = jps, MIN(jpe,jde-1)
     ! do k = kpe-1, kps+1, -1
	   ! do i = ips, MIN(ipe,ide-1)
	     ! head_grid%w_2(i,k,j) = real(2,8)*state_bak(i,k,j) - head_grid%w_2(i,k+1,j) + &
		                        ! state_tmp(i,k,j)
	     ! head_grid%ph_2(i,k,j) = real(2,8)*state_bak1(i,k,j) - head_grid%ph_2(i,k+1,j) + &
		                        ! state_tmp1(i,k,j)
       ! enddo
     ! enddo
   ! enddo
   
   state_bak = real(0,8)
   state_tmp = real(0,8)
   state_bak(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     reshape(state_p_bak(cnt1+cnt2+1:cnt1+cnt2+cnt3), &
     [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     reshape(state_p_amb(cnt1+cnt2+1:cnt1+cnt2+cnt3), &
     [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   state_tmp(ips:MIN(ipe,ide-1), kpe, jps:MIN(jpe,jde-1)) =  &
     state_tmp(ips:MIN(ipe,ide-1), kpe-1, jps:MIN(jpe,jde-1))*real(2,8) - &
     state_tmp(ips:MIN(ipe,ide-1), kpe-2, jps:MIN(jpe,jde-1)) 
   state_tmp(ips:MIN(ipe,ide-1), kps+1:kpe, jps:MIN(jpe,jde-1)) = &
   real(0.5,8)*(state_tmp(ips:MIN(ipe,ide-1), kps+1:kpe, jps:MIN(jpe,jde-1)) + &
   state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)))
   
   state_bak1 = real(0,8)
   state_tmp1 = real(0,8)
   state_bak1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     reshape(state_p_bak(cnt1+cnt2+cnt3+1:cnt1+cnt2+cnt3+cnt4), &
     [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   state_tmp1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     reshape(state_p_amb(cnt1+cnt2+cnt3+1:cnt1+cnt2+cnt3+cnt4), &
     [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   state_tmp1(ips:MIN(ipe,ide-1), kpe, jps:MIN(jpe,jde-1)) =  &
     state_tmp1(ips:MIN(ipe,ide-1), kpe-1, jps:MIN(jpe,jde-1))*real(2,8) - &
     state_tmp1(ips:MIN(ipe,ide-1), kpe-2, jps:MIN(jpe,jde-1)) 
   state_tmp1(ips:MIN(ipe,ide-1), kps+1:kpe, jps:MIN(jpe,jde-1)) = &
   real(0.5,8)*(state_tmp1(ips:MIN(ipe,ide-1), kps+1:kpe, jps:MIN(jpe,jde-1)) + &
   state_tmp1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)))
   
   do j = jps, MIN(jpe,jde-1)
     do k = kps+1, kpe
	   do i = ips, MIN(ipe,ide-1)
	     head_grid%w_2(i,k,j) = real(2,8)*state_bak(i,k-1,j) - head_grid%w_2(i,k-1,j) + &
		                        state_tmp(i,k,j)
	     head_grid%ph_2(i,k,j) = real(2,8)*state_bak1(i,k-1,j) - head_grid%ph_2(i,k-1,j) + &
		                        state_tmp1(i,k,j)
       enddo
     enddo
   enddo	  
	 
   ! head_grid%w_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) =  &
   ! reshape(state_p((cnt1+cnt2+1):(cnt1+cnt2+cnt3)), &
   ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   
   ! head_grid%ph_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) =  &
   ! reshape(state_p((cnt1+cnt2+cnt3+1):(cnt1+cnt2+cnt3+cnt4)), &
   ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   
   head_grid%t_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) =  &
   reshape(state_p((cnt1+cnt2+cnt3+cnt4+1):(cnt1+cnt2+cnt3+cnt4+cnt5)), &
   [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])

   DO mm = 1,3
      head_grid%moist(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),m(mm)) =  &
      reshape(state_p((cnt1+cnt2+cnt3+cnt4+cnt5*mm+1):(cnt1+cnt2+cnt3+cnt4+cnt5+cnt5*mm)), &
      [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   ENDDO
   
   DEALLOCATE(state_tmp, state_tmp1)
   DEALLOCATE(state_bak, state_bak1)
   
   ! ! scheme 1                            
   ! cnt1 = (ipe-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) ! U
   ! cnt2 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (jpe-jps+1) ! V
   ! cnt3 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps+1) * (MIN(jpe,jde-1)-jps+1) ! W
   ! cnt4 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps+1) * (MIN(jpe,jde-1)-jps+1) ! PH
   ! cnt5 = (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) ! T
   ! ! qmoist cnt5 * 6
   ! dimp = cnt1+cnt2+cnt3+cnt4+cnt5*4   
   ! if (dimp /= dim_p) write(*,*) 'distribute dimmension of state_p error! ', &
               ! 'required ', dim_p, 'actually ', dimp
               
   ! head_grid%u_2(ips:ipe, kps:(kpe-1), jps:MIN(jpe,jde-1)) =  &
   ! reshape(state_p(1:cnt1), &
   ! [ipe-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])

   ! head_grid%v_2(ips:MIN(ipe,ide-1), kps:(kpe-1), jps:jpe) =  &
   ! reshape(state_p((cnt1+1):(cnt1+cnt2)), &
   ! [MIN(ipe,ide-1)-ips+1,kpe-kps,jpe-jps+1])

   ! head_grid%w_2(ips:MIN(ipe,ide-1), kps:kpe, jps:MIN(jpe,jde-1)) =  &
   ! reshape(state_p((cnt1+cnt2+1):(cnt1+cnt2+cnt3)), &
   ! [MIN(ipe,ide-1)-ips+1,kpe-kps+1,MIN(jpe,jde-1)-jps+1])

   ! head_grid%ph_2(ips:MIN(ipe,ide-1), kps:kpe, jps:MIN(jpe,jde-1)) =  &
   ! reshape(state_p((cnt1+cnt2+cnt3+1):(cnt1+cnt2+cnt3+cnt4)), &
   ! [MIN(ipe,ide-1)-ips+1,kpe-kps+1,MIN(jpe,jde-1)-jps+1])

   ! ! l = 1
   ! ! DO j = jps, MIN(jpe,jde-1)
     ! ! DO k = kps, (kpe-1)
       ! ! DO i = ips, MIN(ipe,ide-1)
         ! ! head_grid%t_2(i,k,j) = state_p(l+cnt1+cnt2+cnt3+cnt4)
         ! ! l = l + 1
       ! ! ENDDO
     ! ! ENDDO
   ! ! ENDDO
   ! head_grid%t_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) =  &
   ! reshape(state_p((cnt1+cnt2+cnt3+cnt4+1):(cnt1+cnt2+cnt3+cnt4+cnt5)), &
   ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])

   ! DO mm = 1,3
      ! head_grid%moist(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),m(mm)) =  &
      ! reshape(state_p((cnt1+cnt2+cnt3+cnt4+cnt5*mm+1):(cnt1+cnt2+cnt3+cnt4+cnt5+cnt5*mm)), &
      ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
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
   print *, maxval(head_grid%moist(:,:,:,1)),minval(head_grid%moist(:,:,:,1)),sum(head_grid%moist(:,:,:,1))
   print *, maxval(head_grid%moist(:,:,:,2)),minval(head_grid%moist(:,:,:,2)),sum(head_grid%moist(:,:,:,2))
   print *, maxval(head_grid%moist(:,:,:,3)),minval(head_grid%moist(:,:,:,3)),sum(head_grid%moist(:,:,:,3))
! endif   

  WRITE (*,*) 'distribute_state_pdaf: DONE!'

END SUBROUTINE distribute_state_pdaf
