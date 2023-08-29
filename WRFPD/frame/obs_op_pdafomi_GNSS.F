!$Id: PDAFomi_obs_op.F90 552 2020-11-18 14:05:35Z lnerger $

!> PDAF-OMI TEMPLATE for observation operators
!!
!! This module contains a template routine for implementing
!! an observation operator using variables prepared in
!! init_dim_obs.
!!
!! Observation operators provided by PDAF can be found in
!! src/PDAFomi_obs_op.F90
!!
MODULE obs_op_pdafomi_GNSS

  USE PDAFomi_obs_f, ONLY: obs_f, PDAFomi_gather_obsstate, debug
  
  ! added by Changliang Shao on 2022.1228
  USE PDAFomi_obs_op, ONLY: PDAFomi_get_interp_coeff_lin  
  USE module_domain, ONLY : head_grid, get_ijk_from_grid
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: deg2rad, n_fields, off_fields, dim_fields, twin_experiment, &
             startn, filtertype, local_range, &
             ids, ide, jds, jde, kds, kde,    &
             ims, ime, jms, jme, kms, kme,    &
             ips, ipe, jps, jpe, kps, kpe

CONTAINS

!-------------------------------------------------------------------------------
!> Observation operator TEMPLATE
!!
!! The routine has to fill the part of the full observation 
!! vector OBS_F_ALL that represents the current observation
!! type.
!! 
!! The routine first applies the observation operator
!! for the current observation type and then calls
!! PDAFomi_gather_obsstate to gather the observation over
!! all processes and fills OBS_F_ALL.
!!
!! Thus, the observation operator has as least two parts
!! 1. process-local computation of observations of
!!    the current type (array OSTATE_P)
!! 2. gathering of full observation vector and filling
!!    of OBS_F_ALL by calling PDAFomi_gather_obsstate
!!
!! The operations in step 1 can include interpolation, 
!! averaging, or computing the observation values as
!! a function of different elements of the state vector STATE_P
!!
!! Please help to advance PDAF-OMI by providing new 
!! observation operator for inclusion into PDAF-OMI
!! under the LGPL license.
!!
!! __Revision history:__
!! * 2020-11 - Lars Nerger - Initial code for template
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obs_op_GNSS_ZTD(thisobs, nrows, state_p, obs_f_all)

    IMPLICIT NONE

! *** Arguments ***
    ! *** The interface can be adapted to the user's needs. 
    ! *** However, OBS_F_ALL is always required.
    ! *** THISOBS and STATE_P are usually required, but one could
    ! *** implement a multi-step observation operator that separates
    ! *** the functional computation from interpolations
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be averaged
    REAL(8), INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL(8), INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)

! *** Local variables ***
    INTEGER :: i, j, k, row, l, ite, jte, cnt1    ! Counter
    REAL(8), ALLOCATABLE :: ostate_p(:)       ! local observed part of state vector
    REAL(8) :: rrows                          ! REAL(8)-value for nrows
	REAL(8)    :: const, part, term1, term2, wzd, hzd, zf
	REAL(8)    :: a_ew, b_ew, coeff, dh, p, t, qv, psfc
	REAL(8)    :: wdk1, wdk2, zdk1, zdk2, zdk3, kappa, gas_constant, cp
	
	zdk1 = 2.2768e-5 
    zdk2 = 2.66e-3 
    zdk3 = 2.8e-7 
    !  Wet delay: 
    wdk1 = 2.21e-7 
    wdk2 = 3.73e-3  
    !  GPS Refractivity constant   
    a_ew = 0.622 
	b_ew = 0.378 
	coeff = (wdk2*1.e8) / 77.6
	gas_constant = 287.0
	cp = 7.0*gas_constant/2.0
	kappa = gas_constant / cp


! *********************************************
! *** Perform application of observation    ***
! *** operator on vector state_p            ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_TEMPLATE -- START'
          WRITE (*,*) '++ OMI-debug: ', debug, '  PDAFomi_obs_op_TEMPLATE -- DESCRIBE OBS OPERATOR'
          WRITE (*,*) '++ OMI-debug obs_op_TEMPLATE:', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_TEMPLATE:', debug, 'number of points to average', nrows
          WRITE (*,*) '++ OMI-debug obs_op_TEMPLATE:', debug, 'thisobs%id_obs_p', thisobs%id_obs_p
       END IF

       ! Consistency check
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_GNSS_ZTD - thisobs%id_obs_p is not allocated'
       END IF


       ! ************************************************************
       ! *** Process-local: Initialize observed part state vector ***
       ! ************************************************************

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
       END IF


       ! *** Initialize ostate_p e.g. by                      ***
       ! *** - selecting elements of the state vector state_p ***
       ! ***      according to thisobs%id_obs_p               ***
       ! *** - deriving an observed variable as a function    ***
       ! ***      of elements of the state vector             ***

       rrows = REAL(nrows)
	   ite = MIN(ipe,ide-1)
	   jte = MIN(jpe,jde-1)
	   cnt1 = ite - ips + 1
	   ! l = 1
       DO j = jps, jte
	       DO i = ips, ite
		       l = i-ips+1 + cnt1*(j-jps)
		       ostate_p(l) = 0.0
			   ! Wet delay:
			   wzd = 0.0
			   psfc = state_p(thisobs%id_obs_p(4*(kpe-kps)+1,l))
		       DO k = kps, kpe-1
			       dh = state_p(thisobs%id_obs_p(k,l))
				   p = state_p(thisobs%id_obs_p(k+(kpe-kps)*3,l))
				   p = p+head_grid%pb(i,k,j)
			       t = state_p(thisobs%id_obs_p(k+(kpe-kps)*1,l))
				   t = (t+300.)*(p/100000.)**kappa
				   qv = state_p(thisobs%id_obs_p(k+(kpe-kps)*2,l))
				   qv = qv/(1.0+qv)
			       const  = dh / a_ew
                   part   = p * qv / t
                   term1  = part * const * wdk1
                   term2  = part * const * wdk2 / t
                   wzd    = wzd + term1 + term2
			   ENDDO
			   ! Hydrostatic delay (Saastamoinen 1972):
			   ! IDEAL:
			   zf = (1.0 - zdk2*cos(2.0*(20.+(real(j)-101.)*15./111.)*deg2rad) - &
			       zdk3*head_grid%ht(i,j))
               hzd = zdk1 * psfc / zf
			   ! REAL:
			   ! zf = (1.0 - zdk2*cos(2.0*head_grid%XLAT(i,j)*deg2rad) - zdk3*head_grid%ht(i,j))
               ! hzd = zdk1 * head_grid%psfc(i,j) / zf
			   !-----To save the ZTD in cm to ztd:
               ostate_p(l) = (wzd + hzd) * 1.e2
			   ! for average or interpolation
               ! DO row = 1, nrows
                  ! ostate_p(l) = ostate_p(l) + state_p(thisobs%id_obs_p(row,l))
               ! END DO
               ! ostate_p(l) = ostate_p(l) / rrows
			   ! CALL PDAFomi_get_interp_coeff_lin(num_gp, n_dim, gpc, oc, icoeff)
			   
		       ! l = l + 1
		   ENDDO
	   ENDDO
	   
	   ! cnt1 = ite - ips + 1
       ! DO j = jps, jte
	     ! DO k = kps, kpe-1
	       ! DO i = ips, ite
		       ! l = i-ips+1 + cnt1*(j-jps)
		       ! ostate_p(l) = 0.0
			   ! ! Wet delay:
			   ! wzd = 0.0
			   ! psfc = state_p(thisobs%id_obs_p(4*(kpe-kps)+1,l))
			   ! dh = state_p(thisobs%id_obs_p(k,l))
			   ! t = state_p(thisobs%id_obs_p(k+(kpe-kps)*1,l))
			   ! qv = state_p(thisobs%id_obs_p(k+(kpe-kps)*2,l))
			   ! p = state_p(thisobs%id_obs_p(k+(kpe-kps)*3,l))
			   ! const  = dh/a_ew
               ! part   = (p+head_grid%pb(i,k,j))* qv/(1.0+qv)/ &
					   ! ((t+300.)*((p+head_grid%pb(i,k,j))/ &
					   ! 100000.)**kappa)
               ! term1  = part * const * wdk1
               ! term2  = part * const * wdk2 / (t+300.)
               ! wzd    = wzd + term1 + term2
			   
			   ! ! Hydrostatic delay (Saastamoinen 1972):
			   ! ! IDEAL:
			   ! zf = (1.0 - zdk2*cos(2.0*(20.-100.*15./111.+real(j)*15./111.)*deg2rad) - &
			       ! zdk3*head_grid%ht(i,j))
               ! hzd = zdk1 * psfc / zf
			   ! ! REAL:
			   ! ! zf = (1.0 - zdk2*cos(2.0*head_grid%XLAT(i,j)*deg2rad) - zdk3*head_grid%ht(i,j))
               ! ! hzd = zdk1 * head_grid%psfc(i,j) / zf
			   ! !-----To save the ZTD in cm to ztd:
               ! ostate_p(l) = (wzd + hzd) * 1.e2
			   ! ! for average or interpolation
               ! ! DO row = 1, nrows
                  ! ! ostate_p(l) = ostate_p(l) + state_p(thisobs%id_obs_p(row,l))
               ! ! END DO
               ! ! ostate_p(l) = ostate_p(l) / rrows
			   ! ! CALL PDAFomi_get_interp_coeff_lin(num_gp, n_dim, gpc, oc, icoeff)
		   ! ENDDO
		 ! ENDDO
	   ! ENDDO


       ! *************************************************
       ! *** Global: Gather full observed state vector ***
       ! ***            THIS IS MANDATORY!             ***
       ! *************************************************

       CALL PDAFomi_gather_obsstate(thisobs, ostate_p, obs_f_all)

       ! *** Clean up
       DEALLOCATE(ostate_p)

       ! Print debug information
       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_GNSS_ZTD -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_GNSS_ZTD
  
  
  SUBROUTINE PDAFomi_obs_op_GNSS_PW(thisobs, nrows, state_p, obs_f_all)

    IMPLICIT NONE

! *** Arguments ***
    ! *** The interface can be adapted to the user's needs. 
    ! *** However, OBS_F_ALL is always required.
    ! *** THISOBS and STATE_P are usually required, but one could
    ! *** implement a multi-step observation operator that separates
    ! *** the functional computation from interpolations
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be averaged
    REAL(8), INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL(8), INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)

! *** Local variables ***
    INTEGER :: i, j, k, row, l, ite, jte, cnt1    ! Counter
    REAL(8), ALLOCATABLE :: ostate_p(:)       ! local observed part of state vector
    REAL(8) :: rrows                          ! REAL(8)-value for nrows
	REAL(8)    :: pw, dh, qv, rho


! *********************************************
! *** Perform application of observation    ***
! *** operator on vector state_p            ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_TEMPLATE -- START'
          WRITE (*,*) '++ OMI-debug: ', debug, '  PDAFomi_obs_op_TEMPLATE -- DESCRIBE OBS OPERATOR'
          WRITE (*,*) '++ OMI-debug obs_op_TEMPLATE:', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_TEMPLATE:', debug, 'number of points to average', nrows
          WRITE (*,*) '++ OMI-debug obs_op_TEMPLATE:', debug, 'thisobs%id_obs_p', thisobs%id_obs_p
       END IF

       ! Consistency check
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_GNSS_PW - thisobs%id_obs_p is not allocated'
       END IF


       ! ************************************************************
       ! *** Process-local: Initialize observed part state vector ***
       ! ************************************************************

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
       END IF


       ! *** Initialize ostate_p e.g. by                      ***
       ! *** - selecting elements of the state vector state_p ***
       ! ***      according to thisobs%id_obs_p               ***
       ! *** - deriving an observed variable as a function    ***
       ! ***      of elements of the state vector             ***

       rrows = REAL(nrows)
	   ite = MIN(ipe,ide-1)
	   jte = MIN(jpe,jde-1)
	   ! l = 1
	   cnt1 = ite - ips + 1
       DO j = jps, jte
	       DO i = ips, ite
		       l = i-ips+1 + cnt1*(j-jps)
		       ostate_p(l) = 0.0
			   ! Wet delay:
			   pw = 0.0
		       DO k = kps, kpe-1
			       dh = state_p(thisobs%id_obs_p(k,l))
				   qv = state_p(thisobs%id_obs_p(k+(kpe-kps)*1,l))
				   qv = qv/(1.0+qv)
			       rho = state_p(thisobs%id_obs_p(k+(kpe-kps)*2,l))
			       pw = pw + dh * qv * rho
			   ENDDO
			   ! To convert the unit of PW to cm:
               ostate_p(l) = 0.1*pw
			   ! for average or interpolation
               ! DO row = 1, nrows
                  ! ostate_p(l) = ostate_p(l) + state_p(thisobs%id_obs_p(row,l))
               ! END DO
               ! ostate_p(l) = ostate_p(l) / rrows
			   ! CALL PDAFomi_get_interp_coeff_lin(num_gp, n_dim, gpc, oc, icoeff)
			   
		       ! l = l + 1
		   ENDDO
	   ENDDO
	   
	   ! cnt1 = ite - ips + 1
       ! DO j = jps, jte
	     ! DO k = kps, kpe-1
	       ! DO i = ips, ite
		       ! l = i-ips+1 + cnt1*(j-jps)
		       ! ostate_p(l) = 0.0
			   ! ! Wet delay:
			   ! pw = 0.0
			   ! dh = state_p(thisobs%id_obs_p(k,l))
			   ! qv = state_p(thisobs%id_obs_p(k+(kpe-kps)*1,l))
			   ! rho = state_p(thisobs%id_obs_p(k+(kpe-kps)*2,l))
			   ! pw = pw + dh * qv/(1.0+qv) * rho			   
			   ! ! To convert the unit of PW to cm:
               ! ostate_p(l) = 0.1*pw
		   ! ENDDO
		 ! ENDDO
	   ! ENDDO


       ! *************************************************
       ! *** Global: Gather full observed state vector ***
       ! ***            THIS IS MANDATORY!             ***
       ! *************************************************

       CALL PDAFomi_gather_obsstate(thisobs, ostate_p, obs_f_all)

       ! *** Clean up
       DEALLOCATE(ostate_p)

       ! Print debug information
       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_GNSS_PW -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_GNSS_PW

END MODULE obs_op_pdafomi_GNSS
