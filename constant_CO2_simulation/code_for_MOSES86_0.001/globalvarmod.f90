! =======================
! MODULE GLOBALVARMOD    
! =======================

module globalvarmod

implicit none

integer, parameter :: NLAT  = 94                    ! number of latitudes
integer, parameter :: NLON  = 192                   ! number of longitudes
integer, parameter :: NHOR  = NLAT*NLON             ! number of gridpoints
integer, parameter :: NLEP  = 2                     ! number of vertical levels
integer, parameter :: NLEV  = NLEP-1                ! number of vertical levels above surface
integer, parameter :: NUMROWS = NHOR/6

real, parameter    :: TMELT  = 273.16               ! Melting point (H2O)
real, parameter    :: SOL_DAY_EARTH= 86400.         ! Solar day Earth    24h 00m 00s
real, parameter    :: RA1_LIQUID   = 610.78         ! Parameter for Magnus-Teten-Formula
real, parameter    :: RA2_LIQUID   =  17.2693882    ! for saturation vapor pressure
real, parameter    :: RA4_LIQUID   =  35.86         ! over liquid water (used on Earth)
real, parameter    :: akap         = 0.286          ! Kappa Earth
real, parameter    :: gascon       = 287.0          ! Gas constant for dry air on Earth
real, parameter    :: RV           = 461.51         ! Gas constant for water vapour
real, parameter    :: vonkarman    = 0.4            ! von karman const.
real, parameter    :: ga           = 9.81           ! gravity
real, parameter    :: snowthresh   = 0.001          ! snow depth threshold (m) below which
                                                    ! snow evaporation is neglected

character(len=44):: swdir
character(len=42):: etdir
character(len=44):: pressdir
character(len=43):: tsfcdir
character(len=46):: t25dir
character(len=44):: spechumdir
character(len=43):: snowdir
character(len=53):: landseadir
character(len=52):: glacierdir
character(len=47):: z0climodir
character(len=46):: precipdir
character(len=48):: snowevapdir
character(len=48):: snowfalldir
character(len=48):: snowmeltdir
character(len=43):: t60dir
character(len=43):: u60dir
character(len=43):: v60dir
	  
integer :: watcmonmeansunit    = 42
integer :: romonmeansunit      = 43
integer :: z0monmeansunit      = 44
integer :: albmonmeansunit     = 45
integer :: vegmonmeansunit     = 46
integer :: laimonmeansunit     = 47
integer :: forestmonmeansunit  = 48
integer :: wmaxmonmeansunit    = 49
integer :: nppmonmeansunit     = 50
integer :: gpplmonmeansunit    = 51
integer :: gppwmonmeansunit    = 52
integer :: cvegmonmeansunit    = 53
integer :: csoilmonmeansunit   = 54
integer :: resmonmeansunit     = 55
integer :: watermonmeansunit   = 56
integer :: rhsmonmeansunit     = 57
integer :: gamonmeansunit      = 58
integer :: rcmonmeansunit      = 59
integer :: rcminmonmeansunit   = 60
integer :: wsf_wuemonmeansunit = 61
integer :: pgsvmonmeansunit    = 62
integer :: etoutmonmeansunit   = 63
integer :: vftmonmeansunit     = 64
integer :: transpmonmeansunit  = 65

integer :: oldmonth = 1
integer :: month    = 1
integer :: day      = 1
integer :: hour     = -99
integer :: header(8)           ! for reading in date
integer :: headerout(8)        ! for writing date
integer :: accumcount = 0      ! monthly output counter

!     namelist parameters
integer :: year     = 1
integer :: simyear  = 1
integer :: endofoldyearfileexist = 0
character(len=12) :: endofoldyearfile
character(len=12) :: endofyearfile
character(len=4)  :: yrstring
character(len=511) :: swprefix
character(len=511) :: etprefix
character(len=511) :: pressprefix
character(len=511) :: tsfcprefix
character(len=511) :: t25prefix
character(len=511) :: spechumprefix
character(len=511) :: snowprefix
character(len=511) :: landseaprefix
character(len=511) :: glacierprefix
character(len=511) :: z0climoprefix
character(len=511) :: precipprefix
character(len=511) :: snowevapprefix
character(len=511) :: snowfallprefix
character(len=511) :: snowmeltprefix
character(len=511) :: t60prefix
character(len=511) :: u60prefix
character(len=511) :: v60prefix

character(len=511) :: watcmonmeansfile
character(len=511) :: romonmeansfile
character(len=511) :: z0monmeansfile
character(len=511) :: albmonmeansfile
character(len=511) :: vegmonmeansfile
character(len=511) :: laimonmeansfile
character(len=511) :: forestmonmeansfile
character(len=511) :: wmaxmonmeansfile
character(len=511) :: nppmonmeansfile
character(len=511) :: gpplmonmeansfile
character(len=511) :: gppwmonmeansfile
character(len=511) :: cvegmonmeansfile
character(len=511) :: csoilmonmeansfile
character(len=511) :: resmonmeansfile
character(len=511) :: watermonmeansfile
character(len=511) :: rhsmonmeansfile
character(len=511) :: gamonmeansfile
character(len=511) :: pgsvmonmeansfile
character(len=511) :: rcminmonmeansfile
character(len=511) :: rcmonmeansfile
character(len=511) :: wsf_wuemonmeansfile
character(len=511) :: etoutmonmeansfile
character(len=511) :: transpmonmeansfile
character(len=511) :: vftmonmeansfile

integer :: ncveg    = 1     ! compute new dcveg (0=keep initial state)
integer :: nlaicarb = 1     ! switch for LAI,dcveg,dcsoil (1/0: prog./clim)

integer :: nrestart =  0  ! 1 for true, 0 for false
integer :: nbiome   = 1     ! switch for vegetation model (1/0 : prog./clim)
integer :: nevap    = 1     ! switch for evaporation (1/0: prog./reanalysis)
integer :: nsoilalb = 1     ! switch for soil albedo (1/0: prog./fixed)
real    :: rlue     = 5.0E-10 ! light use efficiency coefficient
real    :: co2conv  = 4.15E-04 ! change this in future simulations to 4.15E-04!
real    :: tau_veg  = 10.0  ! [years] - gets scaled to seconds elsewhere
real    :: tau_soil = 42.0  ! [years] - gets scaled to seconds elsewhere
real    :: riniveg  =  0.0
real    :: rinifor  =  0.5
real    :: rinisoil =  0.0
real    :: rnbiocats=  0.0
real    :: forgrow  = 1.0
real    :: rlaigrow = 0.5
!real    :: rss      =  0.0     ! soil surface resistance to evaporation (s/m)
real    :: tr_max   =  2.78E-07 ! max possible transpiration rate (m/s) (Knorr 2000)

integer :: n_days_per_year =     365   ! set to 365 for real calendar
real :: solar_day    = SOL_DAY_EARTH   ! Length of solar day [sec]
real :: deltsec      = SOL_DAY_EARTH/4. ! length of each time step [sec]
real :: ra1          = RA1_LIQUID      !
real :: ra2          = RA2_LIQUID      !
real :: ra4          = RA4_LIQUID      !
real :: co2          = 360.0           ! atm. co2 concentration (ppmv)
real :: acpd         = gascon / akap   ! Specific heat for dry air
real :: rdbrv        = gascon / RV     ! rd / rv

real :: albsmin      = 0.4             ! min. albedo for snow
real :: albsmax      = 0.8             ! max. albedo for snow

real :: dfd(NHOR,NLEP) =0.       ! solar radiation downward
real :: dalb(NHOR)  = 0.11             ! albedo
real :: drhs(NHOR)  = 0.  ! surface wetness
real :: dls(NHOR)   = 1.  ! land(1)/sea(0) mask
real :: dz0(NHOR)   = 1.0E-04       ! surface roughness length
real :: dz0climo(NHOR) = 1.0E-04    ! orographic sfc roughness length
real :: devap(NHOR) = 0.            ! surface evaporation (m^3/m^2/s)
real :: dpet(NHOR) = 0.  ! surface potential evapotranspiration (m^3/m^2/s)
real :: dcsoil(NHOR)   =0.0  ! organic carbon stored in soil
real :: dcveg(NHOR)    =0.0  ! organic carbon stored in biomass
real :: dforest(NHOR) = 0.5  ! forest cover (fract.)
real :: dwmax(NHOR)   = 0.1  ! field capacity (m)
real :: dgpp(NHOR)    = 0.0  ! gross primary production (kg C /m2 /s)
real :: dgppl(NHOR)   = 0.0  ! light limited gpp (kg C /m2 /s)
real :: dgppw(NHOR)   = 0.0  ! water limited gpp (kg C /m2 /s)
real :: dlai(NHOR)    = 0.0  ! leaf area index (integer value)
real :: dlitter(NHOR) = 0.0  ! litterfall (kg C /m2 /s)
real :: dnogrow(NHOR) = 0.0  ! no growth allocation (kg C /m2 /s)
real :: dnpp(NHOR)    = 0.0  ! net primary production (kg C /m2 /s)
real :: dres(NHOR)    = 0.0  ! heterotrophic respiration (kg C /m2 /s)
real :: dveg(NHOR)    = 0.5  ! vegetation cover (fract.)
real :: dglac(NHOR)   = 0.   ! glacier mask (0.,1.)
real :: dp(NHOR)      = 100000. ! surface pressure
real :: plai(NHOR)    = 0.5     ! aboveground growth parameter
real :: pgrow(NHOR)   =1.0       ! biomass growth parameter
real :: pz0_max(NHOR) =2.0       ! maximum roughness
real :: pgs(NHOR)     = 1.0
real :: dprecip(NHOR) = 1.0E-08   ! surface precip (not in Plasim!)
!real :: dprs(NHOR)    = 0.  ! Snow Fall (in Plasim)          (m/s)
real :: dwatc(NHOR)   = 0.  ! soil wetness (m)
real :: dwater(NHOR)        ! soil water time tendency (m/s)
real :: dro(NHOR) = 0.      ! surface runoff (m/s) [SEDGES output]
!real :: drunoff(NHOR) = 0.  ! surface runoff (m) [data input]
!real :: dsoilw(NHOR)   = 0.  ! soil wetness (m)
!real :: dsoilwold(NHOR)   = 0.  ! soil wetness from last time step (m)
real :: dsnow(NHOR)   = 0.  ! snow depth (m)
real :: dsnowevap(NHOR)   = 0.  ! (liquid water equivalent) (m/s)
real :: dsnowfall(NHOR)   = 0.  ! (liquid water equivalent) (m/s)
real :: dsnowmelt(NHOR)   = 0.  ! (liquid water equivalent) (m/s)
!real :: dsnowold(NHOR)= 0.  ! snow depth from last time step (m)
!real :: dsmelt(NHOR)  = 0.  ! snow melt (in Plasim) (m/s water eq.)
!real :: dsndch(NHOR)  = 0.  ! snow depth change (m/s water eq.)
real :: dtsoil(NHOR)  = 273.  ! soil temperature at ~0.25m depth (K)
real :: dfdoldsfc(NHOR) = 0. ! sfc solar radiation downward from last time step
real :: dga(NHOR)   = 0.01  ! aerodynamic conductance for land
real :: dgaold(NHOR)= 0.01  ! aerodynamic conductance for land (old value)
real :: dpgsv(NHOR)   = 0.5  ! vegetation conductance (excl. water stress factor)
real :: dtransp(NHOR) = 0.   ! transpiration (m^3/m^2/s) [>= 0]
real :: dvrhs(NHOR)   = 0.5  ! surface conductance (same as drhs except not modified by landmod)
real :: drc(NHOR)     = 1.0 ! canopy resistance
real :: dwsf_wue(NHOR)= 0.2 ! water stress factor for water-use efficiency
real :: drcmin(NHOR)  = 0.0  ! minimum canopy resistance
real :: dvft(NHOR)    = 0.0  ! vegetation temperature limitation function
real :: dtsa(NHOR)    = 273. ! virtual potential temp. (pot. = surface pressure) of lowest model level
real :: desoil(NHOR)    = 0. ! soil evaporation

real :: dt(NHOR,NLEP)   = 273.    ! temperature (K)
real :: dq(NHOR,NLEP)   = 0.      ! spec. humidity
real :: du(NHOR,NLEV)   = 0.      ! u component of wind
real :: dv(NHOR,NLEV)   = 0.      ! v component of wind
real :: sigma(NLEV)     = 0.99881 ! ERA-Interim lowest atmospheric level (60)

! accumulated variables
real :: acveg(NHOR)   = 0.0  ! organic carbon stored in biomass
real :: acsoil(NHOR)  = 0.0  ! organic carbon stored in soil
real :: aforest(NHOR) = 0.0  ! forest cover (fract.)
real :: agppl(NHOR)   = 0.0  ! light limited gpp (kg C /m2 /s)
real :: agppw(NHOR)   = 0.0  ! water limited gpp (kg C /m2 /s)
real :: alai(NHOR)    = 0.0  ! leaf area index (integer value)
real :: anpp(NHOR)    = 0.0  ! net primary production (kg C /m2 /s)
real :: ares(NHOR)    = 0.0  ! heterotrophic respiration (kg C /m2 /s)
real :: aveg(NHOR)    = 0.0  ! vegetation cover (fract.)
real :: az0(NHOR)     = 0.     ! roughness length
real :: awater(NHOR)  = 0.      ! soil water time tendency (m/s)
real :: aro(NHOR) = 0.   ! surface runoff (m/s)
real :: awatc(NHOR)   = 0.   ! soil wetness (m)
real :: awmax(NHOR)   = 0.0  ! field capacity (m)
real :: arhs(NHOR)    = 0.     ! surface wetness
real :: aalb(NHOR)    = 0.       ! albedo

real :: apgsv(NHOR)   = 0.0  ! accum. vegetation conductance (excl. water stress factor)
real :: aga(NHOR)     = 0.0  ! accum. aerodynamic conductance for land
real :: arc(NHOR)     = 0.0  ! accum. canopy resistance
real :: awsf_wue(NHOR)= 0.0  ! accum. water stress factor for water-use efficiency
real :: arcmin(NHOR)  = 0.0  ! accum. minimum canopy resistance
real :: aetout(NHOR)  = 0.0  ! accum. simulated ET
real :: atransp(NHOR) = 0.0  ! accum. transpiration
real :: avft(NHOR)    = 0.0  ! accum. veg. temperature limitation function


end module globalvarmod
