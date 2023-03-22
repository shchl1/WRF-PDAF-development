program da_bdy

!----------------------------------------------------------------------
! Purpose: Generates boundary file by using wrfinput
!
! Input  : fg         -- first  time level wrfinput generated by real
!          fg02       -- second time level wrfinput generated by real
!          wrfbdy_ref -- reference boundary file generated by real
!
! Output : wrfbdy_out  -- the output boundary file
!
! Notes  : 1. variable name and attributes, dimension name, bdy_width
!             come from wrfbdy. 
!          2. domain size and time come from fg
!          3. boundary and tendency are calculated by using fg & fg02
!          4. the output boundary file only contain the 1st time level
!
! jliu@ucar.edu , 2011-12-15
!----------------------------------------------------------------------

  use netcdf

  implicit none

  integer :: i, n, offset, bdyfrq, domainsize, fg_jd, fg02_jd

  integer :: ncid, ncidfg, ncidfg02, ncidwrfbdy, ncidvarbdy, varid, varid_out, status
  integer :: nDims, nVars, nGlobalAtts, numsAtts
  integer :: dLen, attLen, xtype, unlimDimID
  integer :: bdy_width, varbdy_dimID, wrfbdy_dimID, fg_dimID, vTimes_ID, MSF_ID
  integer :: MU_fgID, MU_fg02ID, MUB_fgID, MUB_fg02ID, fg_varid, fg02_varid, tenid

  integer, dimension(4)  :: dsizes
  integer, dimension(4), target  :: start_u, start_v, start_mass
  integer, dimension(4)          :: cnt_4d, map_4d
  integer, dimension(3)          :: start_3d, cnt_3d, map_3d
  integer, dimension(3), target  :: start_msfu, start_msfv, cnt_msfu, cnt_msfv, map_msfu, map_msfv
  integer, dimension(:), pointer :: start_msf, cnt_msf, map_msf, start_4d

  integer :: south_north, south_north_stag
  integer :: west_east,   west_east_stag
  integer :: bottom_top,  bottom_top_stag

  integer, dimension(nf90_max_var_dims)    :: vDimIDs
  integer, dimension(:),       allocatable :: vdimsizes
  integer, dimension(:,:,:,:), allocatable :: iVar

  real,    dimension(:,:,:,:),  allocatable :: fVar_fg, fVar_fg02, Tend 
  real,    dimension(:,:,:),    allocatable , target :: MU_fg, MU_fg02, MUB_fg, MUB_fg02, MSF

  real,    dimension(:,:,:),    pointer :: MU_fgptr, MU_fg02ptr, MUB_fgptr, MUB_fg02ptr, MSF_ptr 
  
  character (len = 19), dimension(:),  allocatable :: times
  character (len = 19)                             :: fg_time, fg02_time
  character (len = 5)                              :: tenname
  character (len = NF90_MAX_NAME)                  :: vNam, dNam, attNam
  character (len = 9)                              :: MSF_NAME
  character (len = 255)                            :: err_msg=""
  character (len=8)                                :: i_char
  character (len=255)                              :: arg = ""
  character (len=255)                              :: appname =""
  character (len=255)                              :: fg      = "fg"
  character (len=255)                              :: fg02    = "fg02"
  character (len=255)                              :: wrfbdy  = "wrfbdy_ref"
  character (len=255)                              :: varbdy  = "wrfbdy_out"

  logical :: reverse, couple, stag

  integer iargc

  call getarg(0, appname)
  n=index(appname, '/', BACK=.true.)
  appname = trim(appname(n+1:))

  DO i = 1, iargc(), 2
    call getarg(i, arg)
    select case ( trim(arg) )
      case ("-fg")
        call getarg(i+1, arg)
        fg=trim(arg)
      case ("-fg02")
        call getarg(i+1, arg)
        fg02=trim(arg)
      case ("-bdy")
        call getarg(i+1, arg)
        wrfbdy=trim(arg)
      case ("-o")
        call getarg(i+1, arg)
        varbdy=trim(arg)
      case default
        Write(*,*) "Usage : "//trim(appname)//" [-fg filename] [-fg02 filename] [-bdy filename] [-o outputfile] [-h]"
        Write(*,*) "  -fg     Optional, 1st time levle first guess file,         default - fg"
        Write(*,*) "  -fg02   Optional, 2nd time levle first guess file,         default - fg02"
        Write(*,*) "  -bdy    Optional, reference boundary file comes from real, default - wrfbdy_ref"
        Write(*,*) "  -o      Optional, output boundary file,                    default - varbdy_out"
        Write(*,*) "  -h      Show this usage"
        call exit(0)
    end select
  END DO


  status = nf90_open(fg, NF90_NOWRITE, ncidfg)
  if ( status /= nf90_noerr ) then
    err_msg="Failed to open "//trim(fg)
    call nf90_handle_err(status, err_msg)
  endif

  status = nf90_open(fg02, NF90_NOWRITE, ncidfg02)
  if ( status /= nf90_noerr ) then
    err_msg="Failed to open "//trim(fg02)
    call nf90_handle_err(status, err_msg)
  endif

  status = nf90_inq_varid(ncidfg, "Times", vTimes_ID )
  if ( status /= nf90_noerr ) then
    err_msg="Please make sure fg has a vaild Times variable"
    call nf90_handle_err(status, err_msg)
  endif

  status = nf90_get_var(ncidfg, vTimes_ID, fg_time)
  if ( status /= nf90_noerr ) then
    err_msg="Please make sure fg has a vaild Time value"
    call nf90_handle_err(status, err_msg)
  endif

  status = nf90_inq_varid(ncidfg02, "Times", vTimes_ID )
  if ( status /= nf90_noerr ) then
    err_msg="Please make sure fg02 has a vaild Times variable"
    call nf90_handle_err(status, err_msg)
  endif

  status = nf90_get_var(ncidfg02, vTimes_ID, fg02_time)
  if ( status /= nf90_noerr ) then
    err_msg="Please make sure fg02 has a vaild Time value"
    call nf90_handle_err(status, err_msg)
  endif

  status = nf90_open(wrfbdy, NF90_NOWRITE, ncidwrfbdy)
  if ( status /= nf90_noerr ) then
    err_msg="Failed to open "//trim(wrfbdy)
    call nf90_handle_err(status, err_msg)
  endif

  status = nf90_create(varbdy, NF90_CLOBBER, ncidvarbdy)
  if ( status /= nf90_noerr ) then
    err_msg="Please make sure have write access"
    call nf90_handle_err(status, err_msg)
  endif

  bdyfrq = datediff(fg_time, fg02_time) 

  select case ( bdyfrq )
    case ( 0 )
      bdyfrq = 1
    case ( : -1 )
      Write (*,*) "***WARNNING : time levle of fg is LATER then fg02's.***"
  end select

  write(i_char, '(i8)') bdyfrq

  Write(*,*) " Input :"
  Write(*,*) "  fg            "//fg_time
  Write(*,*) "  fg02          "//fg02_time
  Write(*,*) "  Reference bdy "//trim(wrfbdy)
  Write(*,*) "Output : "
  Write(*,*) "  wrfbdy_out    "//fg_time
  Write(*,*) "  bdyfrq        ",adjustl(i_char)

  status = nf90_inquire(ncidfg, nAttributes=nGlobalAtts)
  do i=1, nGlobalAtts
    status = nf90_inq_attname(ncidfg, NF90_GLOBAL, i, attNam)
    status = nf90_copy_att(ncidfg, NF90_GLOBAL, attNam, ncidvarbdy, NF90_GLOBAL)
  end do

  status = nf90_inquire(ncidwrfbdy, nDims, nVars, nGlobalAtts, unlimDimID)
  if ( status /= nf90_noerr ) then
    err_msg="Please make sure have a valid wrf boundary file"
    call nf90_handle_err(status, err_msg)
  endif

  allocate (vdimsizes(nDims), stat=status)

  do i=1, nDims

    status = nf90_inquire_dimension(ncidwrfbdy, i, name=dNam, len = dLen)

    vdimsizes(i) = dLen
    select case (trim(dNam))
      case ("south_north")
        status = nf90_inq_dimid(ncidfg, dNam, fg_dimID)
        status = nf90_inquire_dimension(ncidfg, fg_dimID, len=dLen)
        vdimsizes(i) = dLen
        south_north = vdimsizes(i)
      case ("west_east")
        status = nf90_inq_dimid(ncidfg, dNam, fg_dimID)
        status = nf90_inquire_dimension(ncidfg, fg_dimID, len=dLen)
        vdimsizes(i) = dLen
        west_east = vdimsizes(i)
      case ("south_north_stag")
        status = nf90_inq_dimid(ncidfg, dNam, fg_dimID)
        status = nf90_inquire_dimension(ncidfg, fg_dimID, len=dLen)
        vdimsizes(i) = dLen
        south_north_stag = vdimsizes(i)
      case ("west_east_stag")
        status = nf90_inq_dimid(ncidfg, dNam, fg_dimID)
        status = nf90_inquire_dimension(ncidfg, fg_dimID, len=dLen)
        vdimsizes(i) = dLen
        west_east_stag = vdimsizes(i)
      case ("bottom_top")
        status = nf90_inq_dimid(ncidfg, dNam, fg_dimID)
        status = nf90_inquire_dimension(ncidfg, fg_dimID, len=dLen)
        vdimsizes(i) = dLen
        bottom_top = vdimsizes(i)
      case ("bottom_top_stag")
        status = nf90_inq_dimid(ncidfg, dNam, fg_dimID)
        status = nf90_inquire_dimension(ncidfg, fg_dimID, len=dLen)
        vdimsizes(i) = dLen
        bottom_top_stag = vdimsizes(i)
      case ("Time")
        vdimsizes(i) = 1
        allocate(times(vdimsizes(i)), stat=status)
      case ("bdy_width")
        bdy_width = dLen
    end select

    if  ( i == unlimDimID ) dLen = NF90_UNLIMITED

    status = nf90_def_dim(ncidvarbdy, dNam, dLen, varbdy_dimID)

  end do

  status = nf90_inq_varid(ncidfg  , "MU" , MU_fgID   )
  status = nf90_inq_varid(ncidfg  , "MUB", MUB_fgID  )
  status = nf90_inq_varid(ncidfg02, "MU" , MU_fg02ID )
  status = nf90_inq_varid(ncidfg02, "MUB", MUB_fg02ID)

  status = nf90_inq_varid(ncidfg, "Times", vTimes_ID )

  do varid=1, nVars

    status = nf90_inquire_variable(ncidwrfbdy,varid,name=vNam,xtype=xtype,ndims=nDims,dimids=vDimIDs,natts=numsAtts)
    status = nf90_def_var(ncidvarbdy, trim(vNam), xtype, vDimIDs(1:nDims), varid_out)
    if ( status /= nf90_noerr ) then
      err_msg="Failed to define variable : "//trim(vNam)
      call nf90_handle_err(status, err_msg)
    endif

    do i=1, numsAtts
      status = nf90_inq_attname(ncidwrfbdy, varid, i, attNam)
      status = nf90_copy_att(ncidwrfbdy, varid, trim(attNam), ncidvarbdy, varid_out)
      if ( status /= nf90_noerr ) then
        err_msg="Failed to copy att : "//trim(attNam)
        call nf90_handle_err(status, err_msg)
      endif
    end do

  end do

  status = nf90_enddef(ncidvarbdy)

  do varid=1, nVars

    status = nf90_inquire_variable(ncidwrfbdy,varid,name=vNam,xtype=xtype,ndims=nDims,dimids=vDimIDs)
    if ( status /= nf90_noerr ) then
      err_msg="Failed to inquire varialbe '"//trim(vNam)//"' for wrfbdy"
      call nf90_handle_err(status, err_msg)
    endif

    dsizes = 1
    do i = 1 , nDims
      dsizes(i) = vdimsizes(vDimIDs(i))
    end do

    offset = index(vNam, '_', BACK=.True.) 
    if ( offset <= 0 ) offset = Len(Trim(vNam))

    ! fg
    !     U         (west_east_stag,   south_north,      bottom_top,      time)
    !     V         (west_east,        south_north_stag, bottom_top,      time)
    !     T, QVAPOR (west_east,        south_north,      bottom_top,      time)
    !     PH        (west_east,        south_north,      bottom_top_stag, time)
    !     MU        (west_east,        south_north,                       time)
    !     MAPFAC_U  (west_east_stag,   south_north,                       time)
    !     MAPFAC_V  (west_east,        south_north_stag,                  time)
    ! bdy
    !   west & east   
    !     U         (south_north,      bottom_top,       bdy_width,  time)
    !     V         (south_north_stag, bottom_top,       bdy_width,  time)
    !     T, QVAPOR (south_north,      bottom_top,       bdy_width,  time)
    !     PH        (south_north,      bottom_top_stag,  bdy_width,  time)
    !     MU        (south_north,      bdy_width,                    time)
    !   north & south
    !     U         (west_east_stag,   bottom_top,       bdy_width,  time)
    !     V         (west_east,        bottom_top,       bdy_width,  time)
    !     T, QVAPOR (west_east,        bottom_top,       bdy_width,  time)
    !     PH        (west_east,        bottom_top_stag,  bdy_width,  time)
    !     MU        (west_east,        bdy_width,                    time)

    select case (Trim(vNam(offset:)))
      case ("_BXS") !  West  Boundary
        start_u    = (/1,1,1,1/)
        start_v    = (/1,1,1,1/)
        start_mass = (/1,1,1,1/)
        start_3d   = (/1,1,1/)
        start_msfu = (/1,1,1/)
        start_msfv = (/1,1,1/)

        cnt_4d     = (/dsizes(3),dsizes(1),dsizes(2),1/)
        cnt_3d     = (/bdy_width,south_north,1/)
        cnt_msfu   = (/bdy_width,south_north,1/)
        cnt_msfv   = (/bdy_width,south_north_stag,1/)

        map_4d     = (/dsizes(1)*dsizes(2), 1, dsizes(1), dsizes(1)*dsizes(2)*dsizes(3)/)
        map_3d     = (/south_north, 1, bdy_width*south_north/)
        map_msfu   = (/south_north, 1, bdy_width*south_north/)
        map_msfv   = (/south_north_stag, 1, bdy_width*south_north_stag/)

        reverse    = .False.
        tenname    = "_BTXS"
      case ("_BXE") !  East  Boundary
        start_u    = (/west_east_stag - bdy_width + 1, 1, 1, 1/)
        start_v    = (/west_east - bdy_width + 1, 1, 1, 1/)
        start_mass = (/west_east - bdy_width + 1, 1, 1, 1/)
        start_3d   = (/west_east - bdy_width + 1, 1, 1/)
        start_msfu = (/west_east_stag - bdy_width + 1, 1, 1/)
        start_msfv = (/west_east - bdy_width + 1, 1, 1/)

        cnt_4d     = (/dsizes(3),dsizes(1),dsizes(2),1/)
        cnt_3d     = (/bdy_width,south_north,1/)
        cnt_msfu   = (/bdy_width,south_north,1/)
        cnt_msfv   = (/bdy_width,south_north_stag,1/)

        map_4d     = (/dsizes(1)*dsizes(2), 1, dsizes(1), dsizes(1)*dsizes(2)*dsizes(3)/)
        map_3d     = (/south_north, 1, bdy_width*south_north/)
        map_msfu   = (/south_north, 1, bdy_width*south_north/)
        map_msfv   = (/south_north_stag, 1, bdy_width*south_north_stag/)

        reverse    = .True.
        tenname    = "_BTXE"
      case ("_BYE") !  North Boundary
        start_u    = (/1, south_north      - bdy_width + 1, 1, 1/)
        start_v    = (/1, south_north_stag - bdy_width + 1, 1, 1/)
        start_mass = (/1, south_north      - bdy_width + 1, 1, 1/)
        start_3d   = (/1, south_north      - bdy_width + 1, 1/)
        start_msfu = (/1, south_north      - bdy_width + 1, 1/)
        start_msfv = (/1, south_north_stag - bdy_width + 1, 1/)

        cnt_4d     = (/dsizes(1),dsizes(3),dsizes(2),1/)
        cnt_3d     = (/west_east, bdy_width,1/)
        cnt_msfu   = (/west_east_stag, bdy_width,1/)
        cnt_msfv   = (/west_east, bdy_width,1/)

        map_4d     = (/1, dsizes(1)*dsizes(2), dsizes(1), dsizes(3)*dsizes(1)*dsizes(2)/)
        map_3d     = (/1, west_east, west_east*bdy_width/)
        map_msfu   = (/1, west_east_stag, west_east_stag*bdy_width/)
        map_msfv   = (/1, west_east, west_east*bdy_width/)

        reverse    = .True.
        tenname    = "_BTYE"

      case ("_BYS") !  South Boundary
        start_u    = (/1, 1, 1, 1/)
        start_v    = (/1, 1, 1, 1/)
        start_mass = (/1, 1, 1, 1/)
        start_3d   = (/1, 1, 1/)
        start_msfu = (/1, 1, 1/)
        start_msfv = (/1, 1, 1/)

        cnt_4d     = (/dsizes(1),dsizes(3),dsizes(2),1/)
        cnt_3d     = (/west_east, bdy_width,1/)
        cnt_msfu   = (/west_east_stag, bdy_width,1/)
        cnt_msfv   = (/west_east, bdy_width,1/)

        map_4d     = (/1, dsizes(1)*dsizes(2), dsizes(1), dsizes(3)*dsizes(1)*dsizes(2)/)
        map_3d     = (/1, west_east, west_east*bdy_width/)
        map_msfu   = (/1, west_east_stag, west_east_stag*bdy_width/)
        map_msfv   = (/1, west_east, west_east*bdy_width/)

        reverse    = .False.
        tenname    = "_BTYS"

      case ("_BTXS", "_BTXE","_BTYS","_BTYE")
        cycle
    end select

    select case (nDims)
      case (2)
        if (vNam(1:offset) == "Times") then
          ncid = ncidfg
        else
          n = index(vNam, "bdytime") 
          if ( n <= 0 ) cycle
          select case (vNam(n-4:n-1))
            case ("this")
              ncid = ncidfg
            case ("next") 
              ncid = ncidfg02
            case default
              cycle
          end select
        end if
        status = nf90_get_var(ncid, vTimes_ID, times)
        status = nf90_put_var(ncidvarbdy, varid, times)
      case (3,4)

        Write(*,*) "Processing for "//trim(vNam)

        couple = .true.

        allocate(MU_fg   (dsizes(1),bdy_width,1), stat=status)
        allocate(MU_fg02 (dsizes(1),bdy_width,1), stat=status)
        allocate(MUB_fg  (dsizes(1),bdy_width,1), stat=status)
        allocate(MUB_fg02(dsizes(1),bdy_width,1), stat=status)
        allocate(MSF     (dsizes(1),bdy_width,1), stat=status)

        allocate(Tend(dsizes(1), dsizes(2), dsizes(3), dsizes(4)), stat=status)

        if ( dsizes(1) == west_east_stag .or. dsizes(1) == south_north_stag ) then
           MU_fgptr    => MU_fg   (2:,:,:)
           MU_fg02ptr  => MU_fg02 (2:,:,:)
           MUB_fgptr   => MUB_fg  (2:,:,:)
           MUB_fg02ptr => MUB_fg02(2:,:,:)
           stag = .True.
        else
           MU_fgptr    => MU_fg
           MU_fg02ptr  => MU_fg02
           MUB_fgptr   => MUB_fg
           MUB_fg02ptr => MUB_fg02
           stag = .False.
        end if

        err_msg="Failed to get variable : "//trim(vNam)
        status =  nf90_get_var(ncidfg, MU_fgID, MU_fgptr, start=start_3d, count=cnt_3d, map=map_3d)
        if ( status /= nf90_noerr ) call nf90_handle_err(status, err_msg)

        status =  nf90_get_var(ncidfg02, MU_fg02ID, MU_fg02ptr, start=start_3d,count=cnt_3d, map=map_3d)
        if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)

        status =  nf90_get_var(ncidfg, MUB_fgID, MUB_fgptr, start=start_3d, count=cnt_3d,map=map_3d)
        if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)

        status =  nf90_get_var(ncidfg02, MUB_fg02ID, MUB_fg02ptr, start=start_3d, count=cnt_3d, map=map_3d)
        if(status /= nf90_noerr) call nf90_handle_err(status, err_msg)

        err_msg="Failed to inquire tendency id for "//trim(vNam)//" for output file"
        status = nf90_inq_varid(ncidvarbdy, vNam(1:offset-1)//tenname, tenid)
        if(status /= nf90_noerr) call nf90_handle_err(status, err_msg)

        if ( reverse ) then
          MU_fg    = MU_fg   (:,bdy_width:1:-1,:)
          MU_fg02  = MU_fg02 (:,bdy_width:1:-1,:)
          MUB_fg   = MUB_fg  (:,bdy_width:1:-1,:)
          MUB_fg02 = MUB_fg02(:,bdy_width:1:-1,:)
        end if
         
        select case (vNam(1:offset))
          case ("U_", "V_")
            if ( stag ) then
              MU_fg   (1,:,:) = MU_fg   (2,:,:)
              MU_fg02 (1,:,:) = MU_fg02 (2,:,:)
              MUB_fg  (1,:,:) = MUB_fg  (2,:,:)
              MUB_fg02(1,:,:) = MUB_fg02(2,:,:)

              MU_fg   (2:dsizes(1)-1,:,:) = (MU_fg   (2:dsizes(1)-1,:,:) + MU_fg   (3:dsizes(1),:,:))*0.5
              MU_fg02 (2:dsizes(1)-1,:,:) = (MU_fg02 (2:dsizes(1)-1,:,:) + MU_fg02 (3:dsizes(1),:,:))*0.5
              MUB_fg  (2:dsizes(1)-1,:,:) = (MUB_fg  (2:dsizes(1)-1,:,:) + MUB_fg  (3:dsizes(1),:,:))*0.5
              MUB_fg02(2:dsizes(1)-1,:,:) = (MUB_fg02(2:dsizes(1)-1,:,:) + MUB_fg02(3:dsizes(1),:,:))*0.5
            else
              MU_fg   (:,2:bdy_width,:) = (MU_fg   (:,1:bdy_width-1,:) + MU_fg   (:,2:bdy_width,:))*0.5
              MU_fg02 (:,2:bdy_width,:) = (MU_fg02 (:,1:bdy_width-1,:) + MU_fg02 (:,2:bdy_width,:))*0.5
              MUB_fg  (:,2:bdy_width,:) = (MUB_fg  (:,1:bdy_width-1,:) + MUB_fg  (:,2:bdy_width,:))*0.5
              MUB_fg02(:,2:bdy_width,:) = (MUB_fg02(:,1:bdy_width-1,:) + MUB_fg02(:,2:bdy_width,:))*0.5
            end if

            if ( vNam(1:offset) == "U_"  ) then
              start_4d  => start_u
              start_msf => start_msfu
              cnt_msf   => cnt_msfu
              map_msf   => map_msfu
              MSF_NAME  = "MAPFAC_U"
            else
              start_4d  => start_v
              start_msf => start_msfv
              cnt_msf   => cnt_msfv
              map_msf   => map_msfv
              MSF_NAME  = "MAPFAC_V"
            end if

            status = nf90_inq_varid(ncidfg  , MSF_NAME , MSF_ID   )
            err_msg="Failed to get varialbe MSF"
            status =  nf90_get_var(ncidfg, MSF_ID, MSF, start=start_msf, count=cnt_msf, map=map_msf)
            if(status /= nf90_noerr) call nf90_handle_err(status, err_msg)

            if ( reverse ) MSF = MSF(:,bdy_width:1:-1,:)

          case ("T_","PH_","QVAPOR_")
            MSF = 1.0
            start_4d => start_mass
          case ("MU_")
            status = nf90_inq_varid(ncidvarbdy, "MU"//tenname, tenid)
            Tend(:,:,:,1) = ( MU_fg02 - MU_fg ) / bdyfrq
            status = nf90_put_var(ncidvarbdy, varid, MU_fg)
            !status = nf90_put_var(ncidvarbdy, varid, MU_fg02)
            err_msg="Failed to put variable "//trim(vNam)
            if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)
            status = nf90_put_var(ncidvarbdy, tenid, Tend(:,:,:,1))
            err_msg="Failed to put tendency for "//trim(vNam)
            if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)
            couple = .false.

          case default
            Tend = 0.0
            couple = .false.
            select case (xtype)
              case (nf90_float)
                allocate(fVar_fg( dsizes(1), dsizes(2), dsizes(3), dsizes(4) ), stat=status)
                fVar_fg = 0.0
                status = nf90_put_var(ncidvarbdy, varid, fVar_fg)
                err_msg="Failed to put variable "//trim(vNam)
                if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)
                status = nf90_put_var(ncidvarbdy, tenid, Tend)
                err_msg="Failed to put tendency for "//trim(vNam)
                if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)
                deallocate (fVar_fg)
              case (nf90_int)
                allocate(iVar( dsizes(1), dsizes(2), dsizes(3), dsizes(4) ), stat=status)
                iVar = 0
                status = nf90_put_var(ncidvarbdy, varid, iVar)
                err_msg="Failed to put variable "//trim(vNam)
                if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)
                status = nf90_put_var(ncidvarbdy, tenid, Tend)
                err_msg="Failed to put tendency for "//trim(vNam)
                if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)
                deallocate (iVar)
            end select ! end of xtype

        end select ! end of vNam

        if ( couple ) then

          allocate(  fVar_fg(dsizes(1), dsizes(2), dsizes(3), dsizes(4)), stat=status)
          allocate(fVar_fg02(dsizes(1), dsizes(2), dsizes(3), dsizes(4)), stat=status)

          err_msg="Failed to inquire variable id for "//vNam(1:offset-1)//" for fg"
          status = nf90_inq_varid(ncidfg, vNam(1:offset-1), fg_varid)
          if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)

          err_msg="Failed to inquire variable id for "//vNam(1:offset-1)//" for fg02"
          status = nf90_inq_varid(ncidfg02, vNam(1:offset-1), fg02_varid)
          if(status /= nf90_noerr) call nf90_handle_err(status, err_msg)

          err_msg="Failed to inquire tendency id for "//trim(vNam(1:offset-1))//" for output file"
          status = nf90_inq_varid(ncidvarbdy, vNam(1:offset-1)//tenname, tenid)
          if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)

          err_msg="Failed to get variable "//vNam(1:offset-1)//" from fg"
          status = nf90_get_var(ncidfg, fg_varid, fVar_fg, start=start_4d, count=cnt_4d, map=map_4d)
          if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)

          err_msg="Failed to get variable "//vNam(1:offset-1)//" from fg02"
          status = nf90_get_var(ncidfg02, fg02_varid, fVar_fg02, start=start_4d, count=cnt_4d, map=map_4d)
          if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)

          MU_fg   = MU_fg + MUB_fg
          MU_fg02 = MU_fg02 + MUB_fg
          !MU_fg02 = MU_fg02 + MUB_fg02

          if ( reverse ) then
            fVar_fg   = fVar_fg  (:,:,bdy_width:1:-1,:)
            fVar_fg02 = fVar_fg02(:,:,bdy_width:1:-1,:)
          end if

          do i = 1, dsizes(2)
            fVar_fg(:,i,:,:)   = (fVar_fg  (:,i,:,:) * MU_fg  ) / MSF
            fVar_fg02(:,i,:,:) = (fVar_fg02(:,i,:,:) * MU_fg02) / MSF
          end do

          Tend = ( fVar_fg02 - fVar_fg ) / bdyfrq

          err_msg="Failed to put variable "//trim(vNam)
          status = nf90_put_var(ncidvarbdy, varid, fVar_fg)
          !status = nf90_put_var(ncidvarbdy, varid, fVar_fg02)
          if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)

          err_msg="Failed to put tendency for "//trim(vNam)
          status = nf90_put_var(ncidvarbdy, tenid, Tend)
          if(status /= nf90_noerr) call nf90_handle_err(status,err_msg)

          deallocate (fVar_fg)
          deallocate (fVar_fg02)

        end if

        NULLIFY (MU_fgptr)
        NULLIFY (MU_fg02ptr)
        NULLIFY (MUB_fgptr)
        NULLIFY (MUB_fg02ptr)
        NULLIFY (MSF_ptr)

        deallocate (Tend)
        deallocate (MU_fg)
        deallocate (MU_fg02)
        deallocate (MUB_fg)
        deallocate (MUB_fg02)
        deallocate (MSF)
      case default
        cycle
    end select ! end of nDims

  end do

  deallocate (times)

  status = nf90_close(ncidfg)
  status = nf90_close(ncidfg02)
  status = nf90_close(ncidwrfbdy)
  status = nf90_close(ncidvarbdy)

  Write(*,*) "Boundary file generated successfully"

contains

  subroutine nf90_handle_err(status, err_msg)
    integer,         intent  (in) :: status
    character (len=*), intent(in) :: err_msg

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      print *, trim(err_msg)
      call exit(-1)
    end if
  end subroutine nf90_handle_err

  function jd(yyyy, mm, dd) result(ival)

    integer, intent(in)  :: yyyy
    integer, intent(in)  :: mm
    integer, intent(in)  :: dd
    integer              :: ival

    ! DATE ROUTINE JD(YYYY, MM, DD) CONVERTS CALENDER DATE TO
    ! JULIAN DATE.  SEE CACM 1968 11(10):657, LETTER TO THE
    ! EDITOR BY HENRY F. FLIEGEL AND THOMAS C. VAN FLANDERN.
    ! EXAMPLE JD(1970, 1, 1) = 2440588

    ival = dd - 32075 + 1461*(yyyy+4800+(mm-14)/12)/4 +  &
       367*(mm-2-((mm-14)/12)*12)/12 - 3*((yyyy+4900+(mm-14)/12)/100)/4

    return
  end function jd

  function datediff(date_1, date_2) result(ival)

    character(len=*), intent(in) :: date_1
    character(len=*), intent(in) :: date_2
    integer                      :: ival
    integer                      :: jd1, jd2
    integer                      :: yyyy,mm,dd
    integer                      :: hh1,nn1,ss1
    integer                      :: hh2,nn2,ss2
    

    ! date string : yyyy-mm-dd_hh:mm:ss
    ! calculate the difference between date_1 and date_2 in seconds

    read(date_1(1:19), '(i4,5(1x,i2))') &
         yyyy, mm, dd, hh1, nn1, ss1

    jd1=jd(yyyy,mm,dd)

    read(date_2(1:19), '(i4,5(1x,i2))') &
         yyyy, mm, dd, hh2, nn2, ss2

    jd2=jd(yyyy,mm,dd)

    ival=(jd2-jd1)*86400 + ( hh2-hh1)*3600 + (nn2-nn1)*60 + (ss2-ss1)

    return
  end function datediff

end program da_bdy
