! ========================
! VEGETATION MODULE SEDGES              ! by Pablo Paiewonsky
! ========================
! MOSES86

module vegmod
!use landmod
use globalvarmod
implicit none

integer :: jhor

! parameters for soil moisture stress
real, parameter :: vwfrac_crit_ss   = 1.00 ! critical soil wetness fraction for soil surface evaporation
real, parameter :: vwfrac_crit_lai  = 0.05 ! critical soil wetness fraction for leaf fall


! parameter for vegetation component
! all values are in m-k-s units

!real, parameter :: q10         = 4.29287155354      ! q10 value for low-temp soil resp

!     * parameters for land surface parameters

real, parameter :: zlaimin     =  0.05              ! max LAI for bare soil
real, parameter :: zlaimax     =  7.0              ! max possible LAI
real, parameter :: albveg      =  0.12              ! snow-free albedo for fully-vegetated surface
real, parameter :: albpeat     =  0.12              ! min albedo for dynamic snow-free bare soil
real, parameter :: albsand     =  0.32              ! max albedo for dynamic snow-free bare soil
real, parameter :: albdsnflatmin = 0.40              ! min albedo for deep, pure snow
real, parameter :: albdsnflatmax = 0.80              ! max albedo for deep, pure snow
real, parameter :: albsnformax   = 0.30              ! max snow albedo for forest
! real, parameter :: vwmax_max =  1.00              ! max SWHC for veg
real, parameter :: vwmax_min   =  0.05              ! SWHC for bare soil
real, parameter :: vz0_min     =  0.01              ! min roughness for bare soil
real, parameter :: cveg1      =  0.20              ! for forest cover to biomass conversion
real, parameter :: cveg2      =  1.0               ! biomass threshold for forest cover
real, parameter :: cveg3      =  1.570796326794897 ! Pi/2
real, parameter :: cveg4      =  1.5               ! affects snow albedo drop with biomass
real, parameter :: cveg5      =  1.5               ! biomass threshold for snow albedo drop
real, parameter :: cveg6      =  0.195             ! for biomass to maxLAI conversion
real, parameter :: cveg7      =  9.                ! soil carbon alb. sat. (Williamson et al. 2006, eqn. 30)
real, parameter :: cveg8      =  43.3336822        ! for normalizing soil resp. to old SimBA value at 10C
real, parameter :: cveg9      =  106.              ! for soil respiration (Clark et al., 2011)
real, parameter :: cveg10      = 227.13            ! for soil respiration (k32 in Williamson et al. 2006)
real, parameter :: cveg11      = 0.0012248818259   ! for soil resp. (k0 in Williamson et al. 2006, eqn. 38)
real, parameter :: cveg12      = 0.10              ! multiplier in root pipe model
real, parameter :: cveg13      = 0.63661977236758  ! 2/pi
real, parameter :: cveg14      = 0.7               ! canopy light extinct. coeffic. (1) * clumping index (0.7)
real, parameter :: cveg15      = 8.               ! for biomass to veg roughness conversion (kg/m^2)
real, parameter :: cveg16      = 0.5               ! for biomass to veg roughness conversion (m^2/kg^2)
real, parameter :: cveg17      =  2.5              ! approx. max roughness for veg
real, parameter :: tcrit       = 20.0              ! low temp at which productivity starts to drop (˚C)
real, parameter :: rcmin_min   = 0.0               ! minimum canopy resistance (s/m)
real, parameter :: rcmax       = 1.0E30            ! maximum canopy resistance (s/m)
real, parameter :: wsf_trmin   = 0.                ! minimum water stress factor for transpiration
real, parameter :: ci_ca       = 0.80              ! ci/ca (i.e. internal co2 / atmospheric co2)
real, parameter :: co2_ref     = 360.              ! co2 reference concentration in ppmv
real, parameter :: co2_comp    = 40.               ! light compensation point in ppmv
real, parameter :: zsnowfracf = 0.12               ! snow cover fraction of snow-covered forest
real, parameter :: rssmin     = 10.               ! minimum bare soil resistance (s/m)
real, parameter :: rssmax     = 1.0E30             ! max bare soil resistance (s/m)

!real,dimension(NHOR) :: dvsoil =   0.          ! soil cover

integer :: ibiomass
real    :: zlaim                 ! lai under moist soil conditions
real    :: zvegm                 ! veg cover under moist soil conditions
real    :: zvft_old               ! temperature limitation multiplier from previous time step
real    :: zsvp                  ! saturation vapor pressure
!real    :: zvpd                  ! vapor pressure deficit
!real    :: zbeta
real    :: zco2_mult
real    :: zco2_norm             ! normalizes zco2_mult to = 1 at the reference co2 level
real    :: zglacfree             ! glacier-free fraction

real :: zvwt                   =   0.0         ! vegetation weighting in dvrhs
real :: zveg_old               =   1.0         ! leaf cover fraction for a grid cell at top of loop
real :: zrc_old                =   rcmin_min   ! canopy resistance for a grid cell at top of loop
real :: zrcmin                 =   rcmin_min   ! min canopy resistance under water stress
real :: zwsf_ss                =   1.0         ! water stress factor for soil surface
real :: zwsf_tr                =   1.0         ! water stress factor for transpiration
real :: zfvegdry               =   1.0         ! water-limited leaf cover fraction
real :: zvalb                  =   1.          ! final albedo computed by vegetation module
real :: zalb0                  =   1.          ! snow-free albedo
real :: zalbsoil               =   1.          ! snow-free soil albedo
real :: zalbsnfor              =   1.          ! albedo of snow-covered forest
real :: zalbdsnflat            =   0.          ! deep snow albedo for flat soil/veg
real :: zalbsnflat             =   0.          ! snow albedo for flat soil/veg
real :: zforest                =   0.0         ! forest cover = dforest(:)
real :: zwmax                  =   0.0         ! bucket depth = dwmax(:)
real :: zvz0                   =   vz0_min     ! roughness length
real :: zra_old                =   rcmin_min   ! old aerodynamic resistance
real :: zra                    =   rcmin_min   ! current aerodynamic resistance
real :: zt                     =   0.           ! celsius sfc. temperature / tcrit
real :: zrss                   =   rssmax      ! bare soil resistance (s/m)
real :: zz0const

end module vegmod


!========!
! VEGINI !
!========!

subroutine vegini
use vegmod
implicit none

if (nrestart == 0) then
   where (dls(:) == 0) dforest(:) = 0.0
endif

end subroutine vegini


!=========!
! VEGSTOP !
!=========!

subroutine vegstop
use vegmod
implicit none

return
end subroutine vegstop


!=========!
! VEGSTEP !
!=========!

subroutine vegstep
use vegmod
implicit none

! Initialization and copy of puma fields needed

! zbeta  = max(0.,1.+co2_sens*log((co2-co2_comp)/(co2_ref-co2_comp)))
zco2_norm = (co2_ref + 2. * co2_comp)/(co2_ref - co2_comp)
zco2_mult = max(0.,zco2_norm * (co2 - co2_comp)/(co2 + 2. * co2_comp))
zz0const  = cveg17/(1.+exp(-1.*cveg16*(-cveg15))) - vz0_min

! Following plasim arrays are used but not modified
! -------------------------------------------------
! dwatc(:) = soil wetness [m]
! dswfl(:) = short wave radiation [W/m2]
! devap(:) = surface evaporation (negative)

! Initialize arrays (declared in plasimmod)

dgpp(:)     = 0.0  ! Gross primary production [kg C m-2 s-1]
dgppl(:)    = 0.0  ! GPP light limited        [kg C m-2 s-1]
dgppw(:)    = 0.0  ! GPP water limited        [kg C m-2 s-1]
dnpp(:)     = 0.0  ! Net primary production   [kg C m-2 s-1]
dwsf_wue(:) = 1. - ci_ca
if (nlaicarb == 0) then
  dveg(:) = 1.0 - exp(-cveg14 * dlai(:))
endif


do jhor = 1 , NHOR
  if (dls(jhor) > 0.0 .and. dglac(jhor) < 0.90) then ! land cell with < 90 % glacier

!    zsnowfree = 1.0 - max(dglac(jhor),tanh(100.*dsnow(jhor))) ! snow- and glacier- free fraction
    zglacfree = 1.0 - dglac(jhor) ! glacier free fraction

    ! Make local copies of some plasim arrays
    ! They are copied back if NBIOME is set to 1 (interactive)
    zforest = dforest(jhor) ! Forest cover (0.0 - 1.0)
    zwmax   = dwmax(jhor)   ! Bucket depth

    ! save old values
    zrc_old         = drc(jhor)            ! canopy resistance
	zveg_old        = dveg(jhor)           ! leaf cover fraction
    zvft_old        = dvft(jhor)           ! temperature limitation function

    ! get aerodynamic resistances
    zra_old         = 1./dgaold(jhor)
    zra             = 1./dga(jhor)

    ! diagnostic array: dpgsv
    if (zrc_old > 0.) then
      dpgsv(jhor) = drcmin(jhor)/zrc_old ! varies from small positive to 1: near 1 => water-limited
    else
      dpgsv(jhor) = 0.
    endif

    ! calculate gross primary productivity

    zvwt = zveg_old                            ! vegetation weighting in surface conductance

    ! light limited gpp [kg C / m2 /s]

    dgppl(jhor) = rlue * zco2_mult * zvft_old * dfdoldsfc(jhor) * zveg_old


! diagnostics
!    if (jhor == 8990) then
!      write(*,*) 'dgppl is ', dgppl(jhor)
!	  write(*,*) 'dveg is ', dveg(jhor)
!	endif

    ! water limited gpp [kg C / m2 /s]
    ! co2    : initialized with 360.0 [ppmv] - member of namelist $RADPAR
    ! co2conv:                               - member of namelist $LANDPAR

    dgppw(jhor) = co2conv/1000. * co2 * dwsf_wue(jhor) * zvwt * dp(jhor)/dt(jhor,NLEP)/gascon &
    &           / (1.6*zrc_old + zra_old)

    ! combine light & water limitations (eq. 6.1 PlaSim RM)

    dgpp(jhor) = min(dgppl(jhor),dgppw(jhor))

    ! net primary production [kg C / m2 /s]
    ! pgrow : initialized to 1.0 from variable forgrow from $LANDPAR

    dnpp(jhor) = pgrow(jhor) * 0.5 * dgpp(jhor) * zglacfree

    ! litterfall [kg C / m2 /s]
    ! tau_veg : initialized to 10 [years] - member of namelist $LANDPAR

    dlitter(jhor) = dcveg(jhor) / tau_veg

    ! soil respiration [kg C / m2 /s]
    ! tau_soil : initialized to 42 [years] - member of namelist $LANDPAR

!    dres(jhor) = q10**((dt(jhor,NLEP)-TMELT-10.)/10.) * dcsoil(jhor)/tau_soil
    if (dtsoil(jhor) > 256.1) then
	  dres(jhor) = dcsoil(jhor)/tau_soil * cveg8 / (1.+exp(cveg9/(dtsoil(jhor)-254.85)))
	else ! (dtsoil(jhor) <= 256.1)
	  dres(jhor) = 0.
	endif

!    if (dtsoil(jhor) >= TMELT) then
!      if(min(1.0,max(0.,dwatc(jhor)/zwmax)) > 0.5) then
!	    dres(jhor) = (1. - 0.8 * (min(1.0,max(0.,dwatc(jhor)/zwmax)) - 0.5)) * dres(jhor)
!      else ! soil wetness fraction <= 0.5
!        dres(jhor) = (0.2 + 0.8 * min(1.0,max(0.,dwatc(jhor)/zwmax)) / 0.5) * dres(jhor)
!      endif
!	else ! no unfrozen soil moisture
!	  dres(jhor) = 0.2 * dres(jhor)
!	endif


    ! update carbon stored in biomass (ncveg = time accelerator)

    dcveg(jhor) = dcveg(jhor) + (dnpp(jhor) - dlitter(jhor)) * deltsec * ncveg * nlaicarb !! edited by FL
	 
	dcveg(jhor)= max(0.,dcveg(jhor))

    ! update carbon stored in soil

    dcsoil(jhor) = dcsoil(jhor) + (dlitter(jhor) - dres(jhor)) * deltsec * ncveg * nlaicarb !! edited by FL
	 
	dcsoil(jhor)= max(0.,dcsoil(jhor))

    ! update forest cover

    zforest = 1.-exp(-1.*cveg1*(dcveg(jhor)-cveg2))
    zforest = max(0.0,zforest)
!    zforest = min(1.0,max(0.0,zforest))
		
	! update "bucket" depth
	
	!    zwmax(jhor) = vwmax_max * dvsoil(jhor) + vwmax_min * (1.0 - dvsoil(jhor))
    !    zwmax(jhor) = vwmax_min + (vwmax_max - vwmax_min) * atan(dcveg(jhor))/cveg3
    zwmax = max(vwmax_min, cveg12 * sqrt(dcveg(jhor)))

    ! structurally limited leaf area index (eq. 6.6 PlaSim RM)
    ! zlaimax : set to 12.0 - local constant
    ! zlaimin : set to 0.01 - minimal LAI for bare, wet soil

!    zlaim = zlaimax * zforest(jhor) + zlaimin * (1.0 - zforest(jhor))
    zlaim = zlaimin + cveg13 * atan(cveg6 * dcveg(jhor)) * (zlaimax-zlaimin)
     ! where cveg13 = 2/pi

    ! structurally limited vegetation cover (eq. 6.5 PlaSim RM)
    ! cveg14 = canopy light extinction coefficient (1) * clumping index (0.7) = 0.7

    zvegm = 1.0 - exp(-cveg14 * zlaim)

    ! calculate water stress factors for bare soil and transpiration
        ! water stress = soil wetness / (field capacity * critical soil wetness)
    ! calculate water-limited vegetation fraction
	
    if (zwmax > 0.0) then
!	  zwsf_ss  = min(1.0,max(0.0,dwatc(jhor)/(zwmax * vwfrac_crit_ss)))  ! range: 0 to 1
      zwsf_ss  = min(1.0,max(0.0,(dwatc(jhor)/(zwmax * vwfrac_crit_ss))**2))  ! range: 0 to 1
	  zwsf_tr   = min(1.0,max(wsf_trmin,dwatc(jhor)/zwmax))   ! range: wsf_trmin to 1
	  zfvegdry = min(1.0,max(0.0,dwatc(jhor)/(zwmax * vwfrac_crit_lai))) ! range: 0 to 1
    endif
	
	! update green leaf fraction with new water stress factor when LAI is dynamic
    if (nlaicarb > 0) then
	  dveg(jhor) = min(zvegm,zfvegdry)
    endif
    zvwt = dveg(jhor)

    ! update temperature limitation function
    dvft(jhor) = min(1.0,max(0.0,(dt(jhor,NLEP) - TMELT) / tcrit))
    if (dsnow(jhor) > snowthresh) then
      dvft(jhor) = 0.
    endif


	! update canopy resistance and minimum canopy resistance

    if (zwsf_tr > 0.) then
      zrcmin = min(rcmax,max(rcmin_min,(-1.*dpet(jhor)/zwsf_tr/tr_max - 1.)*zra))
    else
      zrcmin = rcmax
    endif

    if (dfd(jhor,NLEP) <= 0. .or. dwatc(jhor) <= 0. .or. dvft(jhor) <= 0. .or. zco2_mult <= 0.) then
      drc(jhor) = rcmax
    elseif (zveg_old <= 0.) then ! wetting of parched soil
      drc(jhor) = zrcmin
    elseif (dgppl(jhor) == 0.0) then
      drc(jhor) = (((1.6*zrc_old+zra_old)*dgppw(jhor)/rlue/zco2_mult/zveg_old                &
     &          /dvft(jhor)/dfd(jhor,NLEP)) - zra) / 1.6
      drc(jhor) = max(zrcmin, min(rcmax,drc(jhor)))
    else
      drc(jhor) = (((1.6*zrc_old+zra_old)*dgppw(jhor)*dfdoldsfc(jhor)*zvft_old               &
     &          /dgppl(jhor)/dvft(jhor)/dfd(jhor,NLEP)) - zra) / 1.6
      drc(jhor) = max(zrcmin, min(rcmax,drc(jhor)))
    endif

    ! update soil resistance

    if (zwsf_ss > 0.) then
      zrss=min(rssmax,rssmin/zwsf_ss)
    else
      zrss=rssmax
    endif

    ! derivation of remaining land surface parameters

	dvrhs(jhor) = zvwt/(1. + drc(jhor) * dga(jhor))                                          &
     &          + (1. - zvwt)/(1. + zrss * dga(jhor))
    zvz0        = cveg17/(1.+exp(-1.*cveg16*(dcveg(jhor)-cveg15))) - zz0const
!  zvz0        = min(vz0_max, vz0_min + (dcveg(jhor)**2)/cveg16)
!  zvz0        = vz0_min + (vz0_max - vz0_min)*zforest**2
    if (nlaicarb > 0) then
      dlai(jhor)  = -log(1.0 - dveg(jhor))/cveg14
    endif
    if (nsoilalb > 0) then
      zalbsoil = max(albpeat,(albpeat-albsand) / cveg7 * dcsoil(jhor) + albsand)
    else ! fixed soil albedo
      zalbsoil = 0.30
    endif
    zalb0       = albveg * dveg(jhor) + zalbsoil * (1.0 - dveg(jhor))

    ! soil cover diagnostic
!	dvsoil(jhor) = zwmax(jhor)/vwmax_max
!    dvsoil(jhor) = min(1.0,max(0.0,dvsoil(jhor)))

! transpiration diagnostic
    ! T = ET*T/ET
    dtransp(jhor)     = max(0.,-1.*dpet(jhor)*dveg(jhor)*zglacfree/(1.+drc(jhor)*dga(jhor)))

    ! discretization of vegetation state

    if (rnbiocats >= 2.0) then
      ibiomass      = zforest * rnbiocats
      zforest = min(1.0, real(ibiomass)/(rnbiocats-1.0))
!      ibiomass      = dvsoil(jhor) * rnbiocats
!      dvsoil(jhor)  = min(1.0, real(ibiomass)/(rnbiocats-1.0))
    endif


    ! modify albedo and surface conductance due to snow

    if (dsnow(jhor) > 0.0) then

      ! compute snow albedo (use ECHAM5 param. for snow fraction and temp. dependency)

      zalbdsnflat = (albdsnflatmax - albdsnflatmin) * (dt(jhor,NLEP)-268.16) / (TMELT-268.16)
!      zalbdsnflat = max(albsmin,min(albsmax,albsmax-zalbdsnflat))
      zalbdsnflat = max(albdsnflatmin,min(albdsnflatmax,albdsnflatmax-zalbdsnflat))
	  zalbsnflat  = zalb0+(zalbdsnflat-zalb0)*0.95*tanh(100.*dsnow(jhor))
      zalbsnfor   = min(zalbsnflat,albsnformax)
      zvalb       = (zalbsnflat-zalbsnfor) * exp(-cveg4 * max(0.,dcveg(jhor)-cveg5))              &
     &            +zalbsnfor
      if (dsnow(jhor) > snowthresh) then
        dvrhs(jhor) = (1.-zforest)*(1.-dveg(jhor))*(1.-tanh(100.*dsnow(jhor)))/(1.+zrss*dga(jhor))  &
     &              +(1.-zforest)*tanh(100.*dsnow(jhor)) + zforest * zsnowfracf
       ! soil evaporation
        desoil(jhor)  = (1.-zforest)*(1.-dveg(jhor))*(1.-tanh(100.*dsnow(jhor)))                    &
     &                /(1.+zrss*dga(jhor))*dpet(jhor)
      endif
    else
      zvalb       = zalb0
    endif

    ! interactive coupling

    if (nbiome == 1) then
      dz0(jhor)     = sqrt(zvz0*zvz0+dz0climo(jhor)*dz0climo(jhor))
      dwmax(jhor)   = zwmax
      drhs(jhor)    = dvrhs(jhor)
      dalb(jhor)    = zvalb
      dforest(jhor) = zforest
      drcmin(jhor)  = zrcmin	  
    endif

    ! update the old (i.e. previous time step's) surface short wave downward radiation
    dfdoldsfc(jhor) = dfd(jhor,NLEP)

  endif ! (dls(jhor) > 0.0 .and. dglac(jhor) < 0.9)
enddo ! jhor

! Send some arrays to GUI

!call guihor("DCVEG"   // char(0),dcveg  ,1,1000.0,0.0)
!call guihor("ZFOREST" // char(0),zforest,1,1000.0,0.0)

return
end subroutine vegstep
