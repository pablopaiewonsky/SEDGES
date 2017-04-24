      ! =======================
      ! main program
      ! =======================
      ! Drives SEDGES offline using external data input

      program main
      use globalvarmod
      implicit none

	  integer :: i = 1
	  integer :: j = 1

!     scale taus from years to seconds

      tau_veg  = tau_veg  * n_days_per_year * solar_day
      tau_soil = tau_soil * n_days_per_year * solar_day

      call readyearnl
	  
      call readdata2nl
	  
	  if (endofoldyearfileexist == 1) then
	    call readeoyfile
		! diagnostics
!		print *, 'call made to readeoyfile'
	  else ! starting from scratch
	    dwatc(:)  = 0.
		dwmax(:)  = 0.05
! bug fixed: used to be dwmax(:) = 0.1		
		where (dls(:) == 0)
		  dcveg(:)   = 0.
		  dcsoil(:)  = 0.
		  dforest(:) = 0.
		elsewhere
		  dcveg(:)  = riniveg
		  dcsoil(:) = rinisoil
		  dforest(:)= rinifor
		endwhere
	  endif
	  	  
	  open(11,file=swdir//'/'//TRIM(ADJUSTL(swprefix))//yrstring//'.sra',form='formatted')
!	  open(12,file=etdir//'/'//TRIM(ADJUSTL(etprefix))//yrstring//'.sra',form='formatted')
!	  open(13,file=soilwdir//'/'//TRIM(ADJUSTL(soilwprefix))//yrstring//'.sra',form='formatted')
!	  open(14,file=runoffdir//'/'//TRIM(ADJUSTL(runoffprefix))//yrstring//'.sra',form='formatted')
	  open(15,file=pressdir//'/'//TRIM(ADJUSTL(pressprefix))//yrstring//'.sra',form='formatted')
	  open(16,file=tsfcdir//'/'//TRIM(ADJUSTL(tsfcprefix))//yrstring//'.sra',form='formatted')
	  open(17,file=t25dir//'/'//TRIM(ADJUSTL(t25prefix))//yrstring//'.sra',form='formatted')
	  open(18,file=spechumdir//'/'//TRIM(ADJUSTL(spechumprefix))//yrstring//'.sra',form='formatted')
	  open(19,file=snowdir//'/'//TRIM(ADJUSTL(snowprefix))//yrstring//'.sra',form='formatted')
	  open(20,file=landseadir//'/'//TRIM(ADJUSTL(landseaprefix))//'.sra',form='formatted')
	  open(21,file=glacierdir//'/'//TRIM(ADJUSTL(glacierprefix))//'.sra',form='formatted')
      open(22,file=z0climodir//'/'//TRIM(ADJUSTL(z0climoprefix))//'.sra',form='formatted')
      open(23,file=precipdir//'/'//TRIM(ADJUSTL(precipprefix))//yrstring//'.sra',form='formatted')
      open(24,file=snowevapdir//'/'//TRIM(ADJUSTL(snowevapprefix))//yrstring//'.sra',form='formatted')
      open(25,file=snowfalldir//'/'//TRIM(ADJUSTL(snowfallprefix))//yrstring//'.sra',form='formatted')
      open(26,file=snowmeltdir//'/'//TRIM(ADJUSTL(snowmeltprefix))//yrstring//'.sra',form='formatted')
	  open(27,file=t60dir//'/'//TRIM(ADJUSTL(t60prefix))//yrstring//'.sra',form='formatted')
	  open(28,file=u60dir//'/'//TRIM(ADJUSTL(u60prefix))//yrstring//'.sra',form='formatted')
	  open(29,file=v60dir//'/'//TRIM(ADJUSTL(v60prefix))//yrstring//'.sra',form='formatted')

!	  open(monmeanunit,file=monmeansfile,form='formatted')
	  open(watcmonmeansunit,file=TRIM(ADJUSTL(watcmonmeansfile)),form='formatted')
	  open(romonmeansunit,file=TRIM(ADJUSTL(romonmeansfile)),form='formatted')
	  open(z0monmeansunit,file=TRIM(ADJUSTL(z0monmeansfile)),form='formatted')
	  open(albmonmeansunit,file=TRIM(ADJUSTL(albmonmeansfile)),form='formatted')
	  open(vegmonmeansunit,file=TRIM(ADJUSTL(vegmonmeansfile)),form='formatted')
	  open(laimonmeansunit,file=TRIM(ADJUSTL(laimonmeansfile)),form='formatted')
      open(forestmonmeansunit,file=TRIM(ADJUSTL(forestmonmeansfile)),form='formatted')
	  open(wmaxmonmeansunit,file=TRIM(ADJUSTL(wmaxmonmeansfile)),form='formatted')
	  open(nppmonmeansunit,file=TRIM(ADJUSTL(nppmonmeansfile)),form='formatted')
	  open(gpplmonmeansunit,file=TRIM(ADJUSTL(gpplmonmeansfile)),form='formatted')
	  open(gppwmonmeansunit,file=TRIM(ADJUSTL(gppwmonmeansfile)),form='formatted')
	  open(cvegmonmeansunit,file=TRIM(ADJUSTL(cvegmonmeansfile)),form='formatted')
	  open(csoilmonmeansunit,file=TRIM(ADJUSTL(csoilmonmeansfile)),form='formatted')
	  open(resmonmeansunit,file=TRIM(ADJUSTL(resmonmeansfile)),form='formatted')
	  open(watermonmeansunit,file=TRIM(ADJUSTL(watermonmeansfile)),form='formatted')
	  open(rhsmonmeansunit,file=TRIM(ADJUSTL(rhsmonmeansfile)),form='formatted')
      open(gamonmeansunit,file=TRIM(ADJUSTL(gamonmeansfile)),form='formatted')
      open(rcmonmeansunit,file=TRIM(ADJUSTL(rcmonmeansfile)),form='formatted')
      open(rcminmonmeansunit,file=TRIM(ADJUSTL(rcminmonmeansfile)),form='formatted')
      open(wsf_wuemonmeansunit,file=TRIM(ADJUSTL(wsf_wuemonmeansfile)),form='formatted')
      open(pgsvmonmeansunit,file=TRIM(ADJUSTL(pgsvmonmeansfile)),form='formatted')
      open(etoutmonmeansunit,file=TRIM(ADJUSTL(etoutmonmeansfile)),form='formatted')
      open(vftmonmeansunit,file=TRIM(ADJUSTL(vftmonmeansfile)),form='formatted')
      open(transpmonmeansunit,file=TRIM(ADJUSTL(transpmonmeansfile)),form='formatted')

	  call readmasks
	  close(20)
	  close(21)
	  close(22)

	  month = 1
	  day   = 1
      hour  = 0

      dz0climo(:)= 0.

	  do while (month<=12 .and. day<=31 .and. hour<=18)
	    oldmonth = month
	    call readdatetime
		if (month > oldmonth) then
		  call monthly_output(watcmonmeansunit,romonmeansunit,z0monmeansunit,albmonmeansunit       &
                             ,vegmonmeansunit,laimonmeansunit,forestmonmeansunit,wmaxmonmeansunit  &
							 ,nppmonmeansunit,gpplmonmeansunit,gppwmonmeansunit,cvegmonmeansunit   &
							 ,csoilmonmeansunit,resmonmeansunit,watermonmeansunit,rhsmonmeansunit  &
                             ,gamonmeansunit,rcmonmeansunit,rcminmonmeansunit,wsf_wuemonmeansunit  &
                             ,pgsvmonmeansunit,etoutmonmeansunit,vftmonmeansunit,transpmonmeansunit)
          call accumreset
		endif

        call readforcingdata

        call getpet

		call vegstep

        call getet

		call hydro

		call accumulatevariables
		
		! diagnostics
!		write(*,*) 'Here are the global max and min biomass values for this day'
!		write(*,*) maxval(dcveg)
!		write(*,*) minval(dcveg)

		if (month == 12 .and. day == 31 .and. hour==18) then
		  month=13
		  day  = 1
          call monthly_output(watcmonmeansunit,romonmeansunit,z0monmeansunit,albmonmeansunit &
                             ,vegmonmeansunit,laimonmeansunit,forestmonmeansunit,wmaxmonmeansunit&
                             ,nppmonmeansunit,gpplmonmeansunit,gppwmonmeansunit,cvegmonmeansunit &
                             ,csoilmonmeansunit,resmonmeansunit,watermonmeansunit,rhsmonmeansunit&
                             ,gamonmeansunit,rcmonmeansunit,rcminmonmeansunit,wsf_wuemonmeansunit&
                             ,pgsvmonmeansunit,etoutmonmeansunit,vftmonmeansunit,transpmonmeansunit)
		endif
      end do

	  close(11)
	  close(12)
!	  close(13)
!	  close(14)
	  close(15)
	  close(16)
	  close(17)
	  close(18)
	  close(19)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)

	  close(watcmonmeansunit)
	  close(romonmeansunit)
	  close(z0monmeansunit)
	  close(albmonmeansunit)
	  close(vegmonmeansunit)
	  close(laimonmeansunit)
	  close(forestmonmeansunit)
	  close(wmaxmonmeansunit)
	  close(nppmonmeansunit)
	  close(gpplmonmeansunit)
	  close(gppwmonmeansunit)
	  close(cvegmonmeansunit)
	  close(csoilmonmeansunit)
	  close(resmonmeansunit)
	  close(watermonmeansunit)
	  close(rhsmonmeansunit)
      close(gamonmeansunit)
      close(rcmonmeansunit)
      close(rcminmonmeansunit)
      close(wsf_wuemonmeansunit)
      close(pgsvmonmeansunit)
      close(etoutmonmeansunit)
      close(vftmonmeansunit)
      close(transpmonmeansunit)

	  call endofyear_output
	        
      end program main

!     =======================
!     SUBROUTINE READYEARNL
!     =======================

      subroutine readyearnl
	  use globalvarmod
      implicit none

      namelist /yearpar/ simyear,year,yrstring,ncveg,co2,endofoldyearfileexist,endofoldyearfile &
					   ,endofyearfile,watcmonmeansfile,romonmeansfile,z0monmeansfile            &
                       ,albmonmeansfile,vegmonmeansfile,laimonmeansfile,forestmonmeansfile      &
					   ,wmaxmonmeansfile,nppmonmeansfile,gpplmonmeansfile,gppwmonmeansfile      &
					   ,cvegmonmeansfile,csoilmonmeansfile,resmonmeansfile,watermonmeansfile    &
					   ,rhsmonmeansfile,gamonmeansfile,rcmonmeansfile,rcminmonmeansfile         &
                       ,wsf_wuemonmeansfile,pgsvmonmeansfile,etoutmonmeansfile,vftmonmeansfile  &
                       ,transpmonmeansfile

      open(11,file='year_namelist',form='formatted')
      read(11,yearpar)
      close(11)
      
      end subroutine readyearnl
	  

!     =====================================
!	  SUBROUTINE READDATA2NL
!     =====================================
      subroutine readdata2nl

      use globalvarmod
      implicit none

      namelist /data2par/ swdir,etdir                    &
                       ,pressdir,tsfcdir,t25dir          &
				       ,spechumdir,snowdir               &
					   ,landseadir,glacierdir,z0climodir &
                       ,precipdir,snowevapdir            &
                       ,snowfalldir,snowmeltdir,t60dir   &
                       ,u60dir,v60dir                    &
					   ,swprefix,etprefix                &
					   ,pressprefix                      &
					   ,tsfcprefix,t25prefix             &
					   ,spechumprefix,snowprefix         &
					   ,landseaprefix,glacierprefix      &
                       ,z0climoprefix,precipprefix       &
                       ,snowevapprefix                   &
                       ,snowfallprefix                   &
                       ,snowmeltprefix,t60prefix         &
                       ,u60prefix,v60prefix

      open(11,file='data2_namelist',form='formatted')
      read(11,data2par)
      close(11)

      end subroutine readdata2nl
	  
	  
!     ======================
!     SUBROUTINE READEOYFILE
!     ======================

	  
	  subroutine readeoyfile
      use globalvarmod
      implicit none

	  integer :: unitnum = 11
	!  integer :: head(8)
	  
	  open(unitnum,file=endofoldyearfile,form='formatted')

	  ! read soil water content (m) (dwatc)
	  ! header diagnostics
!	  read(unitnum,*) head(1), head(2), head(3),head(4),head(5), head(6), head(7), head(8)
!	  write(*,*) head(1), head(2), head(3),head(4),head(5), head(6), head(7), head(8)
      read(unitnum,*)
	  call readlanddata(dwatc,unitnum)
	  
	  ! read snow equivalent water content (m) (dsnow)
	  ! header diagnostics
	  !read(unitnum,*) head(1), head(2), head(3),head(4),head(5), head(6), head(7), head(8)
	  !write(*,*) head(1), head(2), head(3),head(4),head(5), head(6), head(7), head(8)
	  !call readlanddata(dsnow,unitnum)

      ! read vegetative cover fraction (dveg)
      read(unitnum,*)
      call readlanddata(dveg,unitnum)
	  
	  ! read field capacity (m) (dwmax)
	  ! header diagnostics
!	  read(unitnum,*) head(1), head(2), head(3),head(4),head(5), head(6), head(7), head(8)
!	  write(*,*) head(1), head(2), head(3),head(4),head(5), head(6), head(7), head(8)
      read(unitnum,*)
	  call readlanddata(dwmax,unitnum)
	  
	  ! read biomass (dcveg)
	  ! header diagnostics
!	  read(unitnum,*) head(1), head(2), head(3),head(4),head(5), head(6), head(7), head(8)
!	  write(*,*) head(1), head(2), head(3),head(4),head(5), head(6), head(7), head(8)
      read(unitnum,*)
	  call readlanddata(dcveg,unitnum)

	  ! read soil carbon content (dcsoil)
	  read(unitnum,*)
	  call readlanddata(dcsoil,unitnum)

      ! read surface aerodynamic conductance (dga)
      read(unitnum,*)
      call readlanddata(dga,unitnum)

      ! read canopy resistance (drc)
      read(unitnum,*)
      call readlanddata(drc,unitnum)

      ! read water stress factor for water use efficiency (dwsf_wue)
      read(unitnum,*)
      call readlanddata(dwsf_wue,unitnum)

      ! read surface short wave downward radiation (dfdoldsfc)
      read(unitnum,*)
      call readlanddata(dfdoldsfc,unitnum)

      ! read veg. temperature limitation function (dvft)
      read(unitnum,*)
      call readlanddata(dvft,unitnum)

      ! read surface roughness (dz0)
      read(unitnum,*)
      call readlanddata(dz0,unitnum)

	  close(unitnum)

	  end subroutine readeoyfile
	  
!     ===========================
!     SUBROUTINE READLANDDATA
!     ===========================

      subroutine readlanddata(var,unitnumber)
      use globalvarmod
      implicit none
	  
	  integer :: i
	  integer :: j
	  real, intent(inout) :: var(NHOR)
      integer, intent(in) :: unitnumber
	  
	  i = 1
	  	  
	  do j = 1,NUMROWS
	    read(unitnumber,*) var(i),var(i+1),var(i+2),var(i+3),var(i+4),var(i+5)
		i = i + 6
	  end do

      end subroutine readlanddata

!     ===========================
!     SUBROUTINE READFORCINGDATA
!     ===========================

      subroutine readforcingdata
      use globalvarmod
      implicit none

	  integer :: i
	  integer :: j

      i = 1
      j = 1
			  
      do j = 1,NUMROWS
        read(11,*) dfd(i,NLEP),dfd(i+1,NLEP),dfd(i+2,NLEP),dfd(i+3,NLEP)&
        ,dfd(i+4,NLEP),dfd(i+5,NLEP)
!        read(12,*) devap(i),devap(i+1),devap(i+2),devap(i+3),devap(i+4),devap(i+5)
!		  read(13,*) dsoilw(i),dsoilw(i+1),dsoilw(i+2),dsoilw(i+3),dsoilw(i+4),dsoilw(i+5)
!		  read(14,*) drunoff(i),drunoff(i+1),drunoff(i+2),drunoff(i+3),drunoff(i+4),drunoff(i+5)
        read(15,*) dp(i),dp(i+1),dp(i+2),dp(i+3),dp(i+4),dp(i+5)
        read(16,*) dt(i,NLEP),dt(i+1,NLEP),dt(i+2,NLEP),dt(i+3,NLEP),dt(i+4,NLEP),dt(i+5,NLEP)
        read(17,*) dtsoil(i),dtsoil(i+1),dtsoil(i+2),dtsoil(i+3),dtsoil(i+4),dtsoil(i+5)
        read(18,*) dq(i,NLEV),dq(i+1,NLEV),dq(i+2,NLEV),dq(i+3,NLEV),dq(i+4,NLEV),dq(i+5,NLEV)
        read(19,*) dsnow(i),dsnow(i+1),dsnow(i+2),dsnow(i+3),dsnow(i+4),dsnow(i+5)
!          read(22,*) dpevpr(i),dpevpr(i+1),dpevpr(i+2),dpevpr(i+3),dpevpr(i+4),dpevpr(i+5)
        read(23,*) dprecip(i),dprecip(i+1),dprecip(i+2),dprecip(i+3),dprecip(i+4),dprecip(i+5)
        read(24,*) dsnowevap(i),dsnowevap(i+1),dsnowevap(i+2),dsnowevap(i+3),dsnowevap(i+4),dsnowevap(i+5)
        read(25,*) dsnowfall(i),dsnowfall(i+1),dsnowfall(i+2),dsnowfall(i+3),dsnowfall(i+4),dsnowfall(i+5)
        read(26,*) dsnowmelt(i),dsnowmelt(i+1),dsnowmelt(i+2),dsnowmelt(i+3),dsnowmelt(i+4),dsnowmelt(i+5)
        read(27,*) dt(i,NLEV),dt(i+1,NLEV),dt(i+2,NLEV),dt(i+3,NLEV),dt(i+4,NLEV),dt(i+5,NLEV)
        read(28,*) du(i,NLEV),du(i+1,NLEV),du(i+2,NLEV),du(i+3,NLEV),du(i+4,NLEV),du(i+5,NLEV)
        read(29,*) dv(i,NLEV),dv(i+1,NLEV),dv(i+2,NLEV),dv(i+3,NLEV),dv(i+4,NLEV),dv(i+5,NLEV)

        i = i + 6

      end do

      end subroutine readforcingdata

!     ======================
!     SUBROUTINE READMASKS
!     ======================

	  
	  subroutine readmasks
      use globalvarmod
      implicit none
	  
	  integer :: i
	  integer :: j
	  
	  read(20,*)
	  read(21,*)
      read(22,*)
	  
	  i = 1
	  
      do j = 1,NUMROWS
		read(20,*) dls(i),dls(i+1),dls(i+2),dls(i+3),dls(i+4),dls(i+5)
		read(21,*) dglac(i),dglac(i+1),dglac(i+2),dglac(i+3),dglac(i+4),dglac(i+5)
		read(22,*) dz0climo(i),dz0climo(i+1),dz0climo(i+2),dz0climo(i+3),dz0climo(i+4),dz0climo(i+5)
        i = i + 6
      end do
	        
	  end subroutine readmasks

!     =======================
!     SUBROUTINE READDATETIME
!     =======================

	  
	  subroutine readdatetime
      use globalvarmod
      implicit none
	  
	  character (len=8)    :: str
      character (len=6)    :: str2
	  integer :: j

      read(11,*) (header(j), j=1,8)
 !     read(12,*) (header(j), j=1,8)
 !     read(13,*) (header(j), j=1,8)
 !     read(14,*) (header(j), j=1,8)
      read(15,*) (header(j), j=1,8)
      read(16,*) (header(j), j=1,8)
      read(17,*) (header(j), j=1,8)
      read(18,*) (header(j), j=1,8)
	  read(19,*) (header(j), j=1,8)
      read(23,*) (header(j), j=1,8)
      read(24,*) (header(j), j=1,8)
      read(25,*) (header(j), j=1,8)
      read(26,*) (header(j), j=1,8)
      read(27,*) (header(j), j=1,8)
      read(28,*) (header(j), j=1,8)
      read(29,*) (header(j), j=1,8)

      ! parse date into year, month, and day
      write(str,'(i8)') header(3)
      read(str,'(i4,2i2)') year,month,day

      ! extract hour from time
      write(str2,'(i6)') header(4)
      read(str2,'(i2)') hour
	  
	  ! diagnostics:
!	  write(*,*) 'the following comes from readdate subroutine:'
!	  write(*,*) year, month, day
	        
	  end subroutine readdatetime
	  
	  
!     =========================
!     SUBROUTINE MONTHLY_OUTPUT
!     =========================

	  subroutine monthly_output(watcmonmeansunum,romonmeansunum,z0monmeansunum,albmonmeansunum &
                             ,vegmonmeansunum,laimonmeansunum,forestmonmeansunum,wmaxmonmeansunum  &
							 ,nppmonmeansunum,gpplmonmeansunum,gppwmonmeansunum,cvegmonmeansunum   &
							 ,csoilmonmeansunum,resmonmeansunum,watermonmeansunum,rhsmonmeansunum  &
                             ,gamonmeansunum,rcmonmeansunum,rcminmonmeansunum,wsf_wuemonmeansunum  &
                             ,pgsvmonmeansunum,etoutmonmeansunum,vftmonmeansunum,transpmonmeansunum)
      use globalvarmod
      implicit none
	  
 	  integer :: date
      integer :: hhmmss
	  integer, intent(in) :: watcmonmeansunum
	  integer, intent(in) :: romonmeansunum
	  integer, intent(in) :: z0monmeansunum
	  integer, intent(in) :: albmonmeansunum
	  integer, intent(in) :: vegmonmeansunum
	  integer, intent(in) :: laimonmeansunum
	  integer, intent(in) :: forestmonmeansunum
	  integer, intent(in) :: wmaxmonmeansunum
	  integer, intent(in) :: nppmonmeansunum
	  integer, intent(in) :: gpplmonmeansunum
	  integer, intent(in) :: gppwmonmeansunum
	  integer, intent(in) :: cvegmonmeansunum
	  integer, intent(in) :: csoilmonmeansunum
	  integer, intent(in) :: resmonmeansunum
	  integer, intent(in) :: watermonmeansunum
	  integer, intent(in) :: rhsmonmeansunum
	  integer, intent(in) :: gamonmeansunum
	  integer, intent(in) :: rcmonmeansunum
	  integer, intent(in) :: rcminmonmeansunum
	  integer, intent(in) :: wsf_wuemonmeansunum
	  integer, intent(in) :: pgsvmonmeansunum
      integer, intent(in) :: etoutmonmeansunum
      integer, intent(in) :: vftmonmeansunum
      integer, intent(in) :: transpmonmeansunum

	  date   = year*10000+oldmonth*100+1
      hhmmss = hour*10000
	  
	  ! write soil water content (m) (awatc)
	  awatc(:)=awatc(:)/real(accumcount)
	  call writeheader(140,0,date,hhmmss,NLON,NLAT,simyear,watcmonmeansunum)
	  call writedata(awatc,watcmonmeansunum)
	  
      ! write soil runoff (aro)
	  aro(:)=aro(:)/real(accumcount)
	  call writeheader(160,0,date,hhmmss,NLON,NLAT,simyear,romonmeansunum)
	  call writedata(aro,romonmeansunum)
	  
	  ! write surface roughness (az0)
	  az0(:)=az0(:)/real(accumcount)
	  call writeheader(173,0,date,hhmmss,NLON,NLAT,simyear,z0monmeansunum)
	  call writedata(az0,z0monmeansunum)
	  
      ! write surface albedo (aalb)
	  aalb(:)=aalb(:)/real(accumcount)
	  call writeheader(175,0,date,hhmmss,NLON,NLAT,simyear,albmonmeansunum)
	  call writedata(aalb,albmonmeansunum)

      ! write simulated ET (aetout)
      aetout(:)=aetout(:)/real(accumcount)
      call writeheader(182,0,date,hhmmss,NLON,NLAT,simyear,etoutmonmeansunum)
      call writedata(aetout,etoutmonmeansunum)

	  ! write vegetative cover fraction (aveg)
	  aveg(:)=aveg(:)/real(accumcount)
	  call writeheader(199,0,date,hhmmss,NLON,NLAT,simyear,vegmonmeansunum)
	  call writedata(aveg,vegmonmeansunum)
	  
	  ! write leaf area index (alai)
	  alai(:)=alai(:)/real(accumcount)
	  call writeheader(200,0,date,hhmmss,NLON,NLAT,simyear,laimonmeansunum)
	  call writedata(alai,laimonmeansunum)
	  
	  ! write forest cover (aforest)
	  aforest(:)=aforest(:)/real(accumcount)
	  call writeheader(212,0,date,hhmmss,NLON,NLAT,simyear,forestmonmeansunum)
	  call writedata(aforest,forestmonmeansunum)
	  	  
	  ! write field capacity (m) (awmax)
	  awmax(:)=awmax(:)/real(accumcount)
	  call writeheader(229,0,date,hhmmss,NLON,NLAT,simyear,wmaxmonmeansunum)
	  call writedata(awmax,wmaxmonmeansunum)
	  
	  ! net primary production (anpp)
	  anpp(:)=anpp(:)/real(accumcount)
	  call writeheader(301,0,date,hhmmss,NLON,NLAT,simyear,nppmonmeansunum)
	  call writedata(anpp,nppmonmeansunum)
	  
	  ! write light-limited GPP (agppl)
	  agppl(:)=agppl(:)/real(accumcount)
	  call writeheader(302,0,date,hhmmss,NLON,NLAT,simyear,gpplmonmeansunum)
	  call writedata(agppl,gpplmonmeansunum)
	  
	  ! write water-limited GPP (agppw)
	  agppw(:)=agppw(:)/real(accumcount)
	  call writeheader(303,0,date,hhmmss,NLON,NLAT,simyear,gppwmonmeansunum)
	  call writedata(agppw,gppwmonmeansunum)
	  
	  ! write biomass (acveg)
	  acveg(:)=acveg(:)/real(accumcount)
	  call writeheader(304,0,date,hhmmss,NLON,NLAT,simyear,cvegmonmeansunum)
	  call writedata(acveg,cvegmonmeansunum)
	  
	  ! write soil carbon content (acsoil)
	  acsoil(:)=acsoil(:)/real(accumcount)
	  call writeheader(305,0,date,hhmmss,NLON,NLAT,simyear,csoilmonmeansunum)
	  call writedata(acsoil,csoilmonmeansunum)
	  	  	  
	  ! soil respiration (ares)
	  ares(:)=ares(:)/real(accumcount)
	  call writeheader(307,0,date,hhmmss,NLON,NLAT,simyear,resmonmeansunum)
	  call writedata(ares,resmonmeansunum)
	  
      ! write soil water time tendency (awater)
	  ! added new code, 310
	  awater(:)=awater(:)/real(accumcount)
	  call writeheader(310,0,date,hhmmss,NLON,NLAT,simyear,watermonmeansunum)
	  call writedata(awater,watermonmeansunum)
	  
      ! write evaporation restrictor (arhs)
	  ! added new code, 311
	  arhs(:)=arhs(:)/real(accumcount)
	  call writeheader(311,0,date,hhmmss,NLON,NLAT,simyear,rhsmonmeansunum)
	  call writedata(arhs,rhsmonmeansunum)

      ! write vegetation water limitation (apgsv)
      apgsv(:)=apgsv(:)/real(accumcount)
      call writeheader(313,0,date,hhmmss,NLON,NLAT,simyear,pgsvmonmeansunum)
      call writedata(apgsv,pgsvmonmeansunum)

      ! write surface aerodynamic conductance (aga)
      aga(:)=aga(:)/real(accumcount)
      call writeheader(316,0,date,hhmmss,NLON,NLAT,simyear,gamonmeansunum)
      call writedata(aga,gamonmeansunum)

      ! write canopy resistance (arc)
      arc(:)=arc(:)/real(accumcount)
      call writeheader(317,0,date,hhmmss,NLON,NLAT,simyear,rcmonmeansunum)
      call writedata(arc,rcmonmeansunum)

      ! write water stress factor for water use efficiency (awsf_wue)
      awsf_wue(:)=awsf_wue(:)/real(accumcount)
      call writeheader(318,0,date,hhmmss,NLON,NLAT,simyear,wsf_wuemonmeansunum)
      call writedata(awsf_wue,wsf_wuemonmeansunum)

      ! write minimum canopy resistance (arcmin)
      arcmin(:)=arcmin(:)/real(accumcount)
      call writeheader(319,0,date,hhmmss,NLON,NLAT,simyear,rcminmonmeansunum)
      call writedata(arcmin,rcminmonmeansunum)

      ! write transpiration (atransp)
      atransp(:)=atransp(:)/real(accumcount)
      call writeheader(389,0,date,hhmmss,NLON,NLAT,simyear,transpmonmeansunum)
      call writedata(atransp,transpmonmeansunum)

      ! write veg. temperature limitation function (avft)
      avft(:)=avft(:)/real(accumcount)
      call writeheader(391,0,date,hhmmss,NLON,NLAT,simyear,vftmonmeansunum)
      call writedata(avft,vftmonmeansunum)

	  end subroutine monthly_output

	  
!     ===============================
!     SUBROUTINE ACCUMULATEVARIABLES
!     ===============================

	  
	  subroutine accumulatevariables
      use globalvarmod
      implicit none
	  
	  acveg(:)   = acveg(:)+dcveg(:)
      acsoil(:)  = acsoil(:)+dcsoil(:)
	  aforest(:) = aforest(:)+dforest(:)
      agppl(:)   = agppl(:)+dgppl(:)
      agppw(:)   = agppw(:)+dgppw(:)
      alai(:)    = alai(:)+dlai(:)
      anpp(:)    = anpp(:)+dnpp(:)
	  ares(:)    = ares(:)+dres(:)
      aveg(:)    = aveg(:)+dveg(:)
      az0(:)     = az0(:)+dz0(:)
      awater(:)  = awater(:)+dwater(:)
      aro(:) = aro(:)+dro(:)
      awatc(:)   = awatc(:)+dwatc(:)
      awmax(:)   = awmax(:)+dwmax(:)
      arhs(:)    = arhs(:)+drhs(:)
      aalb(:)    = aalb(:)+dalb(:)
      apgsv(:)   = apgsv(:)+dpgsv(:)
      aga(:)     = aga(:)+dga(:)
      arc(:)     = arc(:)+drc(:)
      awsf_wue(:)=awsf_wue(:)+dwsf_wue(:)
      arcmin(:)  =arcmin(:)+drcmin(:)
      aetout(:)  =aetout(:)+devap(:)
      atransp(:) = atransp(:)+dtransp(:)
      avft(:)    =avft(:)+dvft(:)

	  accumcount=accumcount+1

	  end subroutine accumulatevariables
	  
!     ===============================
!     SUBROUTINE ACCUMRESET
!     ===============================

	  
	  subroutine accumreset
      use globalvarmod
      implicit none
	  
	  acveg(:)   = 0.
      acsoil(:)  = 0.
	  aforest(:) = 0.
      agppl(:)   = 0.
      agppw(:)   = 0.
      alai(:)    = 0.
      anpp(:)    = 0.
	  ares(:)    = 0.
      aveg(:)    = 0.
      az0(:)     = 0.
      awater(:)  = 0.
      aro(:)     = 0.
      awatc(:)   = 0.
      awmax(:)   = 0.
      arhs(:)    = 0.
      aalb(:)    = 0.
      apgsv(:)   = 0.
      aga(:)     = 0.
      arc(:)     = 0.
      awsf_wue(:)= 0.
      arcmin(:)  = 0.
      aetout(:)  = 0.
      atransp(:) = 0.
      avft(:)    = 0.

	  accumcount=0

	  end subroutine accumreset

		 
!     ===========================
!     SUBROUTINE ENDOFYEAR_OUTPUT
!     ===========================

	  subroutine endofyear_output
      use globalvarmod
      implicit none

	  integer :: date
      integer :: hhmmss
	  integer :: unitnum = 11
	
	  date = year*10000+12*100+31
      hhmmss = hour*10000
	  
	  open(unitnum,file=endofyearfile,form='formatted')

	  ! write soil water content [bucket model output] (m) (dwatc)
	  call writeheader(140,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
	  call writedata(dwatc,unitnum)
	  
	  ! write snow equivalent water content (m) (dsnow)
	  !call writeheader(141,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
	  !call writedata(dsnow,unitnum)
	  
      ! write vegetative cover fraction (dveg)
      call writeheader(199,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
      call writedata(dveg,unitnum)

	  ! write field capacity (m) (dwmax)
	  call writeheader(229,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
	  call writedata(dwmax,unitnum)
	  
	  ! write biomass (dcveg)
	  call writeheader(304,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
	  call writedata(dcveg,unitnum)
	  
	  ! write soil carbon content (dcsoil)
	  call writeheader(305,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
	  call writedata(dcsoil,unitnum)

      ! write surface aerodynamic conductance (dga)
      call writeheader(316,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
      call writedata(dga,unitnum)

      ! write canopy resistance (drc)
      call writeheader(317,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
      call writedata(drc,unitnum)

      ! write water stress factor for water use efficiency (dwsf_wue)
      call writeheader(318,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
      call writedata(dwsf_wue,unitnum)

      ! write surface short wave downward radiation (dfdoldsfc)
      call writeheader(076,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
      call writedata(dfdoldsfc,unitnum)

      ! write vegetation temperature limitation function (dvft)
      call writeheader(391,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
      call writedata(dvft,unitnum)

      ! write surface roughness (dz0)
      call writeheader(173,0,date,hhmmss,NLON,NLAT,simyear,unitnum)
      call writedata(dz0,unitnum)

	  close(unitnum)
	  
	  end subroutine endofyear_output
	  
!     ===========================
!     SUBROUTINE WRITEHEADER
!     ===========================

      subroutine writeheader(code,level,date,time,lonnum,latnum,sy,unitnumber)
      use globalvarmod
      implicit none
	  
	  integer :: i
	  integer, intent(in) :: code
	  integer, intent(in) :: level
	  integer, intent(in) :: date
	  integer, intent(in) :: time
	  integer, intent(in) :: lonnum
	  integer, intent(in) :: latnum
	  integer, intent(in) :: sy     ! simulation year
	  integer, intent(in) :: unitnumber
	  
	  headerout(1) = code
      headerout(2) = level
      headerout(3) = date
	  headerout(4) = time
      headerout(5) = lonnum
      headerout(6) = latnum
	  headerout(7) = sy
      headerout(8) = 0
	  
	  !diagnostics
!	  write(*,*) 'call to writeheader subroutine'
!      write(*,*) 'unitnumber = ', unitnumber
!	  write(*,*) 'month = ',month
!	  write(*,*) 'oldmonth = ',oldmonth
	  
	  write(unitnumber,*) (headerout(i), i=1,8)
	  
      end subroutine writeheader
	  
!     ===========================
!     SUBROUTINE WRITEDATA
!     ===========================

      subroutine writedata(var,unitnumber)
      use globalvarmod
      implicit none
	  
	  integer :: i
	  integer :: j
	  real, intent(in) :: var(NHOR)
      integer, intent(in) :: unitnumber
	  
	  i = 1
	  do j = 1,NUMROWS
	    write(unitnumber,*) var(i),var(i+1),var(i+2),var(i+3),var(i+4),var(i+5)
		i = i + 6
	  end do

      end subroutine writedata

!     ===========================
!     SUBROUTINE GETPET
!     ===========================

      subroutine getpet
      use globalvarmod
      implicit none

      integer :: jhor = 1
      real    :: zsvp = 1000.  ! saturation vapor pressure (Pa)
      real    :: zqs  = 0.0    ! surface saturation specific humidity (kg/kg)
      real    :: zqdiff
      real    :: zrho = 0.0         ! surface air density (kg/m^2)
      real    :: zexp
      real    :: zlnsig
      real    :: zumin      = 1.    ! minimum wind speed for PBL exhcange (m/s)
      real    :: ztransh
      real    :: vdiff_lamm = 160.  ! const. used in vdiff (see parameterization)
      real    :: vdiff_b    = 5.    !        "
      real    :: vdiff_c    = 5.    !        "
      real    :: vdiff_d    = 5.    !        "
      real    :: zabsu2(NHOR) = 0.       ! squared wind speed (m2/s2)
      real    :: znl(NHOR)    = 0.       ! z of lowermost level (m)
      real    :: zbz0(NHOR)   = 0.       ! z/z0
      real    :: zri(NHOR)    = 0.       ! bulk richardson number
      real    :: zrifh(NHOR)  = 0.       ! factor for heat flux transfer coeff.
      real    :: zkblnz2, zdenom

      zexp=-akap
      zlnsig=ALOG(sigma(NLEV))

      do jhor = 1,NHOR
        if (dls(jhor) > 0.0) then
          zsvp = ra1*exp(ra2*(dt(jhor,NLEP)-TMELT)/(dt(jhor,NLEP)-ra4))
          zqs  = rdbrv*zsvp/dp(jhor)
          zqs  = zqs/(1.-(1./rdbrv-1.)*zqs)
          zqdiff = zqs - dq(jhor,NLEV)
          zrho = dp(jhor)/dt(jhor,NLEP)/gascon
          dq(jhor,NLEP) = zqs

          ! virtual potential temperature (pot. = surface pressure) of lowest model level air
          dtsa(jhor)=dt(jhor,NLEV)*sigma(NLEV)**zexp                             &
          &         *(1.+(1./rdbrv-1.)*dq(jhor,NLEV))

          !*    calculate transfer coefficient
          !
          !     windspeed (squared):
          zabsu2(jhor)=AMAX1(zumin,du(jhor,NLEV)*du(jhor,NLEV)+dv(jhor,NLEV)*dv(jhor,NLEV))

          !     z off lowermost layer
          znl(jhor)=-gascon*0.5*(dt(jhor,NLEV)+dtsa(jhor))*zlnsig/ga

          !     z/z0
          zbz0(jhor)=znl(jhor)/dz0(jhor)

          !     bulk richardson number
          zri(jhor)=ga*znl(jhor)/(zabsu2(jhor)*dtsa(jhor))                             &
          &       *(dtsa(jhor)-dt(jhor,NLEP)*(1.+(1./rdbrv-1.)*dq(jhor,NLEP)))

          zkblnz2=(vonkarman/ALOG(zbz0(jhor)+1.))**2

          if(zri(jhor) <= 0.) then
            zdenom=1.+3.*vdiff_c*vdiff_b                                    &
         &        *sqrt(-1.*zri(jhor)*(zbz0(jhor)+1.))*zkblnz2
            zrifh(jhor)=1.-3.*vdiff_b*zri(jhor)/zdenom
          else
            zdenom=SQRT(1+vdiff_d*zri(jhor))
            zrifh(jhor)=1./(1.+3.*vdiff_b*zri(jhor)*zdenom)
          endif

          !
          !     transfer coeff.
          !

          ztransh=SQRT(zabsu2(jhor))*zkblnz2*zrifh(jhor)

          !     copy transfer coeff. to global array for SEDGES access
          !     and save old value
          dgaold(jhor)=dga(jhor)
          dga(jhor)=ztransh

          dpet(jhor)= -1.*zrho*dga(jhor)*zqdiff/1000.
          ! divide by 1000 converts from bulk formula units (kg/m^2/s) to Plasim units (m^3/m^2/s)
        endif
      end do

      end subroutine getpet

!     ===========================
!     SUBROUTINE GETET
!     ===========================

      subroutine getet
      use globalvarmod
      implicit none

      integer :: jhor = 1

      do jhor = 1,NHOR
          ! derive (explicit) Plasim ET
          if (nevap == 1) then
            devap(jhor)= dpet(jhor)*drhs(jhor)
          endif
      end do


      end subroutine getet
