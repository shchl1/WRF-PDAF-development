!$Id: init_n_domains_pdaf.F90 883 2021-11-27 14:16:40Z lnerger $
!> Set number of local analysis domains
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called in PDAF_X_update 
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains. 
!! It has to set the number of local analysis 
!! domains for the PE-local domain.
!!
!! __Revision history:__
!! * 2005-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)

! added by Changliang Shao on 2022.1228
  USE module_domain, ONLY : head_grid, get_ijk_from_grid 
! added by Changliang Shao on 2023.0311
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: dim_fields
	   
  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step        !< Current time step
  INTEGER, INTENT(out) :: n_domains_p !< PE-local number of analysis domains

! ************************************
! *** Initialize number of domains ***
! ************************************
  ! added by Changliang Shao on 2023.0311
  ! SET n_domains_p. n_domains_p is the grid (not element) number.
  n_domains_p = dim_fields(5)
  ! n_domains_p = (ipe-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1) + &
  ! (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (jpe-jps+1) + &
  ! (MIN(ipe,ide-1)-ips+1) * (kpe-kps+1) * (MIN(jpe,jde-1)-jps+1) + &
  ! (MIN(ipe,ide-1)-ips+1) * (kpe-kps) * (MIN(jpe,jde-1)-jps+1)

END SUBROUTINE init_n_domains_pdaf
