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
       ONLY: n_fields, off_fields, dim_fields,  &
	   ids, ide, jds, jde, kds, kde,    &
	   ims, ime, jms, jme, kms, kme,    &
	   ips, ipe, jps, jpe, kps, kpe
  
  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL(8), INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

! *** local variables ***

! added by Changliang Shao on 2022.1228
  integer                   :: i,j,k,l,mm,dimx,dimy,dimz,dimp,     &
                               cnt1,cnt2,cnt3,cnt4,cnt5,  &
							   its,ite,jts,jte
  integer,dimension(6) :: m
  REAL(8), ALLOCATABLE :: state_tmp(:,:,:),state_tmp1(:,:,:),state_tmp2(:,:,:), &
                          state_tmp3(:,:,:),state_tmp4(:,:,:), &
						  state_bak(:,:,:),state_bak1(:,:,:),state_bak2(:,:,:), &
                          state_bak3(:,:,:),state_bak4(:,:,:), &
                          state_tmpc1(:,:,:),state_tmpc2(:,:,:), &
						  state_tmpc3(:,:,:),state_tmpc4(:,:,:)  ! temp state
  REAL(8) :: state_p_bak(dim_p), state_p_amb(dim_p)

! *************************************************
! *** Initialize model fields from state vector ***
! *** for process-local model domain            ***
!**************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'distribute_state_pdaf: initialization of model fields here!'
  
  data m /2,3,4,5,8,9/  !'QV','QC','QR','QI','QS','QG', QI,QS,QG no value in ideal
				
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  	 
   ! scheme 2   			
			      
   state_p_bak = 0.
   state_p_amb = 0.
   CALL collect_state_pdaf(dim_p, state_p_bak)   
   state_p_amb = state_p - state_p_bak
      
   its = ips
   ite = MIN(ipe,ide-1)
   jts = jps
   jte = MIN(jpe,jde-1)
   
   ALLOCATE(state_tmp1(its:ite, kps:kpe, jts:jte)) ! increment in grid a
   ALLOCATE(state_tmp2(its:ite, kps:kpe, jts:jte))
   ALLOCATE(state_tmp3(its:ite, kps:kpe, jts:jte))
   ALLOCATE(state_tmp4(its:ite, kps:kpe, jts:jte))
   ALLOCATE(state_tmpc1(its:ite, kps:kpe, jts:jte)) ! increment in grid c
   ALLOCATE(state_tmpc2(its:ite, kps:kpe, jts:jte))
   ALLOCATE(state_tmpc3(its:ite, kps:kpe, jts:jte))
   ALLOCATE(state_tmpc4(its:ite, kps:kpe, jts:jte))
   
   state_tmp1 = 0.
   state_tmp2 = 0.
   state_tmp3 = 0.
   state_tmp4 = 0.
   state_tmpc1 = 0.
   state_tmpc2 = 0.
   state_tmpc3 = 0.
   state_tmpc4 = 0.
   
   l = 1 ! a = f + increment
   do j = jts, jte
     do k = kps, kpe-1
	   do i = its, ite
	     ! increment in grid a 
	     state_tmp1(i,k,j) = state_p_amb(off_fields(1)+l)
	     state_tmp2(i,k,j) = state_p_amb(off_fields(2)+l)
	     state_tmp3(i,k,j) = state_p_amb(off_fields(3)+l)
	     state_tmp4(i,k,j) = state_p_amb(off_fields(4)+l)
		 ! update T, QV, QC, QR
	     head_grid%t_2(i,k,j) = state_p(off_fields(5)+l)
	     head_grid%moist(i,k,j,m(1)) = state_p(off_fields(6)+l)
	     head_grid%moist(i,k,j,m(2)) = state_p(off_fields(7)+l)
	     head_grid%moist(i,k,j,m(3)) = state_p(off_fields(8)+l)
	     ! head_grid%moist(i,k,j,m(4)) = state_p(off_fields(9)+l)
	     ! head_grid%moist(i,k,j,m(5)) = state_p(off_fields(10)+l)
	     ! head_grid%moist(i,k,j,m(6)) = state_p(off_fields(11)+l)
		 if (head_grid%moist(i,k,j,m(1)) < 0) head_grid%moist(i,k,j,m(1)) = 0. 
		 if (head_grid%moist(i,k,j,m(2)) < 0) head_grid%moist(i,k,j,m(2)) = 0. 
		 if (head_grid%moist(i,k,j,m(3)) < 0) head_grid%moist(i,k,j,m(3)) = 0. 
		 ! if (head_grid%moist(i,k,j,m(4)) < 0) head_grid%moist(i,k,j,m(4)) = 0. 
		 ! if (head_grid%moist(i,k,j,m(5)) < 0) head_grid%moist(i,k,j,m(5)) = 0. 
		 ! if (head_grid%moist(i,k,j,m(6)) < 0) head_grid%moist(i,k,j,m(6)) = 0. 
		 l = l + 1
       enddo
     enddo
   enddo
   
   do j = jte, jts, -1
     do k = kps, kpe-1
	   do i = ite, its, -1
		 ! for increment U, V
		 if (i /= its) state_tmpc1(i,k,j) = (state_tmp1(i,k,j)+state_tmp1(i-1,k,j))/2.
		 if (j /= jts) state_tmpc2(i,k,j) = (state_tmp2(i,k,j)+state_tmp2(i,k,j-1))/2.
         ! for increment PH, W
		 if (k /= kps) then
		   state_tmpc3(i,k,j) = (state_tmp3(i,k,j)+state_tmp3(i,k-1,j))/2.
		   state_tmpc4(i,k,j) = (state_tmp4(i,k,j)+state_tmp4(i,k-1,j))/2.
		 endif
       enddo	   
     enddo
   enddo
   
   ! ! for edges increment
   state_tmpc1(its,kps:kpe-1,jts:jte) = &
	 2.*state_tmp1(its,kps:kpe-1,jts:jte) - state_tmpc1(its+1,kps:kpe-1,jts:jte)
   state_tmpc2(its:ite,kps:kpe-1,jts) = &
	 2.*state_tmp2(its:ite,kps:kpe-1,jts) - state_tmpc2(its:ite,kps:kpe-1,jts+1)
   ! state_tmpc3(its:ite,kps,jts:jte) = &
	 ! 2.*state_tmp3(its:ite,kps,jts:jte) - state_tmpc3(its:ite,kps+1,jts:jte)
   ! state_tmpc3(its:ite,kpe,jts:jte) = &
	 ! 2.*state_tmp3(its:ite,kpe,jts:jte) - state_tmpc3(its:ite,kpe-1,jts:jte)
   ! state_tmpc4(its:ite,kps,jts:jte) = &
	 ! 2.*state_tmp4(its:ite,kps,jts:jte) - state_tmpc4(its:ite,kps+1,jts:jte)
   ! state_tmpc4(its:ite,kpe,jts:jte) = &
	 ! 2.*state_tmp4(its:ite,kpe,jts:jte) - state_tmpc4(its:ite,kpe-1,jts:jte)	 
   
   ! update stagger U, V
   if (ipe == ide) then
     head_grid%u_2(ipe,kps:kpe-1,jts:jte) = &
	   head_grid%u_2(ipe,kps:kpe-1,jts:jte) + &
       2.*state_tmp1(ipe-1,kps:kpe-1,jts:jte) - &
	   state_tmpc1(ipe-1,kps:kpe-1,jts:jte)
   endif
   if (jpe == jde) then
     head_grid%v_2(its:ite,kps:kpe-1,jpe) = &
	   head_grid%v_2(its:ite,kps:kpe-1,jpe) + &
       2.*state_tmp2(its:ite,kps:kpe-1,jpe-1) - &
	   state_tmpc2(its:ite,kps:kpe-1,jpe-1)
   endif
   
   ! a = f + increment
   head_grid%u_2(its:ite,kps:kpe-1,jts:jte) = &
     state_tmpc1(its:ite,kps:kpe-1,jts:jte) + head_grid%u_2(its:ite,kps:kpe-1,jts:jte)
   head_grid%v_2(its:ite,kps:kpe-1,jts:jte) = &
     state_tmpc2(its:ite,kps:kpe-1,jts:jte) + head_grid%v_2(its:ite,kps:kpe-1,jts:jte)
   head_grid%w_2(its:ite,kps:kpe,jts:jte) = &
     state_tmpc3(its:ite,kps:kpe,jts:jte) + head_grid%w_2(its:ite,kps:kpe,jts:jte)
   head_grid%ph_2(its:ite,kps:kpe,jts:jte) = &
     state_tmpc4(its:ite,kps:kpe,jts:jte) + head_grid%ph_2(its:ite,kps:kpe,jts:jte)
	       	    
   DEALLOCATE(state_tmp1,state_tmp2,state_tmp3,state_tmp4, &              
              state_tmpc1,state_tmpc2,state_tmpc3,state_tmpc4)
			  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  	 
   ! scheme 2   			
			      
   ! state_p_bak = 0.
   ! state_p_amb = 0.
   ! CALL collect_state_pdaf(dim_p, state_p_bak)   
   ! state_p_amb = state_p - state_p_bak
      
   ! its = ips
   ! ite = MIN(ipe,ide-1)
   ! jts = jps
   ! jte = MIN(jpe,jde-1)
   
   ! ALLOCATE(state_tmp1(its:ite, kps:kpe, jts:jte)) ! increment in grid a
   ! ALLOCATE(state_tmp2(its:ite, kps:kpe, jts:jte))
   ! ALLOCATE(state_tmp3(its:ite, kps:kpe, jts:jte))
   ! ALLOCATE(state_tmp4(its:ite, kps:kpe, jts:jte))
   ! ALLOCATE(state_bak1(its:ite, kps:kpe, jts:jte)) ! f in grid a
   ! ALLOCATE(state_bak2(its:ite, kps:kpe, jts:jte))
   ! ALLOCATE(state_bak3(its:ite, kps:kpe, jts:jte))
   ! ALLOCATE(state_bak4(its:ite, kps:kpe, jts:jte))
   ! ALLOCATE(state_tmpc1(its:ite, kps:kpe, jts:jte)) ! increment in grid c
   ! ALLOCATE(state_tmpc2(its:ite, kps:kpe, jts:jte))
   ! ALLOCATE(state_tmpc3(its:ite, kps:kpe, jts:jte))
   ! ALLOCATE(state_tmpc4(its:ite, kps:kpe, jts:jte))
   
   ! state_tmp1 = 0.
   ! state_tmp2 = 0.
   ! state_tmp3 = 0.
   ! state_tmp4 = 0.
   ! state_bak1 = 0.
   ! state_bak2 = 0.
   ! state_bak3 = 0.
   ! state_bak4 = 0.
   ! state_tmpc1 = 0.
   ! state_tmpc2 = 0.
   ! state_tmpc3 = 0.
   ! state_tmpc4 = 0.
   
   ! l = 1 ! a = f + increment
   ! do j = jts, jte
     ! do k = kps, kpe-1
	   ! do i = its, ite
	     ! state_bak1(i,k,j) = state_p_bak(off_fields(1)+l)
	     ! state_bak2(i,k,j) = state_p_bak(off_fields(2)+l)
	     ! state_bak3(i,k,j) = state_p_bak(off_fields(3)+l)
	     ! state_bak4(i,k,j) = state_p_bak(off_fields(4)+l)
	     ! state_tmp1(i,k,j) = state_p_amb(off_fields(1)+l)
	     ! state_tmp2(i,k,j) = state_p_amb(off_fields(2)+l)
	     ! state_tmp3(i,k,j) = state_p_amb(off_fields(3)+l)
	     ! state_tmp4(i,k,j) = state_p_amb(off_fields(4)+l)
		 ! ! update T, QV, QC, QR
	     ! head_grid%t_2(i,k,j) = state_p(off_fields(5)+l)
	     ! head_grid%moist(i,k,j,m(1)) = state_p(off_fields(6)+l)
	     ! head_grid%moist(i,k,j,m(2)) = state_p(off_fields(7)+l)
	     ! head_grid%moist(i,k,j,m(3)) = state_p(off_fields(8)+l)
		 ! l = l + 1
       ! enddo
     ! enddo
   ! enddo
   
   ! ! for f
   ! do j = jte, jts, -1
     ! do k = kps, kpe-1
	   ! do i = ite, its, -1
	     ! head_grid%u_2(i,k,j) = 2.*state_bak1(i,k,j) - head_grid%u_2(i+1,k,j)
	     ! head_grid%v_2(i,k,j) = 2.*state_bak2(i,k,j) - head_grid%v_2(i,k,j+1)
	     ! head_grid%w_2(i,k+1,j) = 2.*state_bak3(i,k,j) - head_grid%w_2(i,k,j)
	     ! head_grid%ph_2(i,k+1,j) = 2.*state_bak4(i,k,j) - head_grid%ph_2(i,k,j)
		 ! ! for increment U, V
		 ! if (i /= its) state_tmpc1(i,k,j) = (state_tmp1(i,k,j)+state_tmp1(i-1,k,j))/2.
		 ! if (j /= jts) state_tmpc2(i,k,j) = (state_tmp2(i,k,j)+state_tmp2(i,k,j-1))/2.
         ! ! for increment PH, W
		 ! if (k /= kps) then
		   ! state_tmpc3(i,k,j) = (state_tmp3(i,k,j)+state_tmp3(i,k-1,j))/2.
		   ! state_tmpc4(i,k,j) = (state_tmp4(i,k,j)+state_tmp4(i,k-1,j))/2.
		 ! endif
       ! enddo	   
     ! enddo
   ! enddo
   
   ! ! for edges increment
   ! state_tmpc1(its,kps:kpe-1,jts:jte) = &
	 ! 2.*state_tmpc1(its+1,kps:kpe-1,jts:jte) - state_tmpc1(its+2,kps:kpe-1,jts:jte)
   ! state_tmpc2(its:ite,kps:kpe-1,jts) = &
	 ! 2.*state_tmpc2(its:ite,kps:kpe-1,jts+1) - state_tmpc2(its:ite,kps:kpe-1,jts+2)
   ! state_tmpc3(its:ite,kps,jts:jte) = &
	 ! 2.*state_tmpc3(its:ite,kps+1,jts:jte) - state_tmpc3(its:ite,kps+2,jts:jte)
   ! state_tmpc3(its:ite,kpe,jts:jte) = &
	 ! 2.*state_tmpc3(its:ite,kpe-1,jts:jte) - state_tmpc3(its:ite,kpe-2,jts:jte)
   ! state_tmpc4(its:ite,kps,jts:jte) = &
	 ! 2.*state_tmpc4(its:ite,kps+1,jts:jte) - state_tmpc4(its:ite,kps+2,jts:jte)
   ! state_tmpc4(its:ite,kpe,jts:jte) = &
	 ! 2.*state_tmpc4(its:ite,kpe-1,jts:jte) - state_tmpc4(its:ite,kpe-2,jts:jte)	 
   
   ! ! update stagger U, V
   ! if (ipe == ide) then
     ! head_grid%u_2(ipe,kps:kpe-1,jts:jte) = &
       ! 2.*state_bak1(ipe-1,kps:kpe-1,jts:jte) - &
	   ! head_grid%u_2(ipe-1,kps:kpe-1,jts:jte) + &
       ! 2.*state_tmpc1(ipe-1,kps:kpe-1,jts:jte) - &
	   ! state_tmpc1(ipe-2,kps:kpe-1,jts:jte)
   ! endif
   ! if (jpe == jde) then
     ! head_grid%v_2(its:ite,kps:kpe-1,jpe) = &
       ! 2.*state_bak2(its:ite,kps:kpe-1,jpe-1) - &
	   ! head_grid%v_2(its:ite,kps:kpe-1,jpe-1) + &
       ! 2.*state_tmpc2(its:ite,kps:kpe-1,jpe-1) - &
	   ! state_tmpc2(its:ite,kps:kpe-1,jpe-2)
   ! endif
   
   ! head_grid%u_2(its:ite,kps:kpe-1,jts:jte) = &
     ! state_tmpc1(its:ite,kps:kpe-1,jts:jte) + head_grid%u_2(its:ite,kps:kpe-1,jts:jte)
   ! head_grid%v_2(its:ite,kps:kpe-1,jts:jte) = &
     ! state_tmpc2(its:ite,kps:kpe-1,jts:jte) + head_grid%v_2(its:ite,kps:kpe-1,jts:jte)
   ! head_grid%w_2(its:ite,kps:kpe,jts:jte) = &
     ! state_tmpc3(its:ite,kps:kpe,jts:jte) + head_grid%w_2(its:ite,kps:kpe,jts:jte)
   ! head_grid%ph_2(its:ite,kps:kpe,jts:jte) = &
     ! state_tmpc4(its:ite,kps:kpe,jts:jte) + head_grid%ph_2(its:ite,kps:kpe,jts:jte)
	       	    
   ! DEALLOCATE(state_tmp1,state_tmp2,state_tmp3,state_tmp4, &
              ! state_bak1,state_bak2,state_bak3,state_bak4, &
              ! state_tmpc1,state_tmpc2,state_tmpc3,state_tmpc4)
			  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
   ! ALLOCATE(state_bak(ips-1:MIN(ipe,ide-1), kps:kpe, jps-1:MIN(jpe,jde-1)))
   ! ALLOCATE(state_tmp(ips-1:MIN(ipe,ide-1), kps:kpe, jps-1:MIN(jpe,jde-1)))
   ! ALLOCATE(state_tmp1(ips-1:MIN(ipe,ide-1), kps:kpe, jps-1:MIN(jpe,jde-1)))   
   
   ! CALL collect_state_pdaf(dim_p, state_p_bak)    
   ! state_p_amb = state_p - state_p_bak
   
   ! state_bak = 0.
   ! state_tmp = 0.
   ! state_tmp1 = 0.
   ! state_bak(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_bak(1+off_fields(1):dim_fields(1)+off_fields(1)), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_amb(1+off_fields(1):dim_fields(1)+off_fields(1)), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   ! state_tmp1(ips-1, kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! state_tmp(ips, kps:kpe-1, jps:MIN(jpe,jde-1))*2. - &
     ! state_tmp(ips+1, kps:kpe-1, jps:MIN(jpe,jde-1)) 
   ! state_tmp1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) = &
   ! (state_tmp(ips-1:MIN(ipe,ide-1)-1, kps:kpe-1, jps:MIN(jpe,jde-1)) + &
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)))/2.
   ! do j = jps, MIN(jpe,jde-1)
     ! do k = kps, kpe-1
	   ! do i = MIN(ipe,ide-1), ips, -1
	     ! head_grid%u_2(i,k,j) = 2.*state_bak(i,k,j) - head_grid%u_2(i+1,k,j) + &
		                        ! state_tmp1(i,k,j)
       ! enddo
     ! enddo
   ! enddo
	    
   ! state_bak = 0.
   ! state_tmp = 0.
   ! state_tmp1 = 0.
   ! state_bak(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_bak(1+off_fields(2):dim_fields(2)+off_fields(2)), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_amb(1+off_fields(2):dim_fields(2)+off_fields(2)), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1]) 
   ! state_tmp1(ips:MIN(ipe,ide-1), kps:kpe-1, jps-1) =  &
     ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps)*2. - &
     ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps+1) 
   ! state_tmp1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) = &
   ! (state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps-1:MIN(jpe,jde-1)-1) + &
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)))/2.
   ! do j = MIN(jpe,jde-1), jps, -1
     ! do k = kps, kpe-1
	   ! do i = ips, MIN(ipe,ide-1)
	     ! head_grid%v_2(i,k,j) = 2.*state_bak(i,k,j) - head_grid%v_2(i,k,j+1) + &
		                        ! state_tmp1(i,k,j)
       ! enddo
     ! enddo
   ! enddo
   
   ! state_bak = 0.
   ! state_tmp = 0.
   ! state_tmp1 = 0.
   ! state_bak(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_bak(1+off_fields(3):dim_fields(3)+off_fields(3)), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_amb(1+off_fields(3):dim_fields(3)+off_fields(3)), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   ! state_tmp1(ips:MIN(ipe,ide-1), kpe, jps:MIN(jpe,jde-1)) =  &
     ! state_tmp(ips:MIN(ipe,ide-1), kpe-1, jps:MIN(jpe,jde-1))*2. - &
     ! state_tmp(ips:MIN(ipe,ide-1), kpe-2, jps:MIN(jpe,jde-1)) 
   ! state_tmp1(ips:MIN(ipe,ide-1), kps+1:kpe, jps:MIN(jpe,jde-1)) = &
   ! (state_tmp(ips:MIN(ipe,ide-1), kps+1:kpe, jps:MIN(jpe,jde-1)) + &
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)))/2.   
   ! do j = jps, MIN(jpe,jde-1)
     ! do k = kps+1, kpe
	   ! do i = ips, MIN(ipe,ide-1)
	     ! head_grid%w_2(i,k,j) = 2.*state_bak(i,k-1,j) - head_grid%w_2(i,k-1,j) + &
		                        ! state_tmp1(i,k,j)
       ! enddo
     ! enddo
   ! enddo
   
   ! state_bak = 0.
   ! state_tmp = 0.
   ! state_tmp1 = 0.
   ! state_bak(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_bak(1+off_fields(4):dim_fields(4)+off_fields(4)), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_amb(1+off_fields(4):dim_fields(4)+off_fields(4)), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   ! state_tmp1(ips:MIN(ipe,ide-1), kpe, jps:MIN(jpe,jde-1)) =  &
     ! state_tmp(ips:MIN(ipe,ide-1), kpe-1, jps:MIN(jpe,jde-1))*2. - &
     ! state_tmp(ips:MIN(ipe,ide-1), kpe-2, jps:MIN(jpe,jde-1)) 
   ! state_tmp1(ips:MIN(ipe,ide-1), kps+1:kpe, jps:MIN(jpe,jde-1)) = &
   ! (state_tmp(ips:MIN(ipe,ide-1), kps+1:kpe, jps:MIN(jpe,jde-1)) + &
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)))/2.   
   ! do j = jps, MIN(jpe,jde-1)
     ! do k = kps+1, kpe
	   ! do i = ips, MIN(ipe,ide-1)
	     ! head_grid%ph_2(i,k,j) = 2.*state_bak(i,k-1,j) - head_grid%ph_2(i,k-1,j) + &
		                        ! state_tmp1(i,k,j)
       ! enddo
     ! enddo
   ! enddo	  
	    
   ! head_grid%t_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) =  &
   ! reshape(state_p(1+off_fields(5):dim_fields(5)+off_fields(5)), &
   ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])

   ! DO mm = 1,3
      ! head_grid%moist(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1),m(mm)) =  &
      ! reshape(state_p(1+off_fields(5+mm):dim_fields(5+mm)+off_fields(5+mm)), &
      ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   ! ENDDO
   
   ! DEALLOCATE(state_bak, state_tmp, state_tmp1)    
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
   
   ! state_bak = 0.
   ! state_tmp = 0.
   ! state_bak(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_bak(cnt1+cnt2+1:cnt1+cnt2+cnt3), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_amb(cnt1+cnt2+1:cnt1+cnt2+cnt3), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   ! state_tmp(ips:MIN(ipe,ide-1), kps-1, jps:MIN(jpe,jde-1)) =  &
     ! state_tmp(ips:MIN(ipe,ide-1), kps, jps:MIN(jpe,jde-1))*2. - &
     ! state_tmp(ips:MIN(ipe,ide-1), kps+1, jps:MIN(jpe,jde-1)) 
   ! state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) = &
   ! (state_tmp(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) + &
   ! state_tmp(ips:MIN(ipe,ide-1), kps-1:kpe-2, jps:MIN(jpe,jde-1)))/2.
   
   ! state_bak1 = 0.
   ! state_tmp1 = 0.
   ! state_bak1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_bak(cnt1+cnt2+cnt3+1:cnt1+cnt2+cnt3+cnt4), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])   
   ! state_tmp1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) =  &
     ! reshape(state_p_amb(cnt1+cnt2+cnt3+1:cnt1+cnt2+cnt3+cnt4), &
     ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   ! state_tmp1(ips:MIN(ipe,ide-1), kps-1, jps:MIN(jpe,jde-1)) =  &
     ! state_tmp1(ips:MIN(ipe,ide-1), kps, jps:MIN(jpe,jde-1))*2. - &
     ! state_tmp1(ips:MIN(ipe,ide-1), kps+1, jps:MIN(jpe,jde-1)) 
   ! state_tmp1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) = &
   ! (state_tmp1(ips:MIN(ipe,ide-1), kps:kpe-1, jps:MIN(jpe,jde-1)) + &
   ! state_tmp1(ips:MIN(ipe,ide-1), kps-1:kpe-2, jps:MIN(jpe,jde-1)))/2.
   
   ! do j = jps, MIN(jpe,jde-1)
     ! do k = kpe-1, kps+1, -1
	   ! do i = ips, MIN(ipe,ide-1)
	     ! head_grid%w_2(i,k,j) = 2.*state_bak(i,k,j) - head_grid%w_2(i,k+1,j) + &
		                        ! state_tmp(i,k,j)
	     ! head_grid%ph_2(i,k,j) = 2.*state_bak1(i,k,j) - head_grid%ph_2(i,k+1,j) + &
		                        ! state_tmp1(i,k,j)
       ! enddo
     ! enddo
   ! enddo
   
   ! head_grid%w_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) =  &
   ! reshape(state_p((cnt1+cnt2+1):(cnt1+cnt2+cnt3)), &
   ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   
   ! head_grid%ph_2(ips:MIN(ipe,ide-1),kps:(kpe-1),jps:MIN(jpe,jde-1)) =  &
   ! reshape(state_p((cnt1+cnt2+cnt3+1):(cnt1+cnt2+cnt3+cnt4)), &
   ! [MIN(ipe,ide-1)-ips+1,kpe-kps,MIN(jpe,jde-1)-jps+1])
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

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
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

! if (1==0) then
   print *, 'head_grid%u_2(ipe+1,kps,jpe)= ', head_grid%u_2(ipe+1,kps,jpe)
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

  WRITE (*,*) 'distribute_state_pdaf: DONE!'

END SUBROUTINE distribute_state_pdaf
