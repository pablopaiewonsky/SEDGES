! =======================
! SUBROUTINE HYDRO    
! =======================
! MOSES86 version
! changes: criteria below are for dsnow wrt threshold, not 0.0 swe depth.

      subroutine hydro
	  use globalvarmod
      implicit none

	  integer :: jhor = 1

      ! preset
      dwater(:)=0.
      dro(:)=0.

      do jhor = 1,NHOR
        if (dls(jhor) > 0.0) then
          if (dsnow(jhor) <= snowthresh) then
            dwater(jhor)=dprecip(jhor)-dsnowfall(jhor)+devap(jhor)+dsnowmelt(jhor)
          elseif (dsnow(jhor) > snowthresh) then
          ! soil moisture can decrease only via snow-free bare soil
            dwater(jhor)=dprecip(jhor)-dsnowfall(jhor)+dsnowmelt(jhor)+desoil(jhor)
          endif
          if (dsnow(jhor) > snowthresh .and. dsnowmelt(jhor) <= 0.) then
            dwater(jhor)=min(0.,dwater(jhor)) ! soil moisture cannot increase
          endif
          dwatc(jhor)=dwatc(jhor)+deltsec*dwater(jhor)   ! in meters
          dro(jhor)=AMAX1(0.,dwatc(jhor)-dwmax(jhor))/deltsec
          dwatc(jhor)=AMAX1(AMIN1(dwmax(jhor),dwatc(jhor)),0.)
        endif
	  end do
      end subroutine hydro