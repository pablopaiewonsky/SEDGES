#!/bin/sh

# driver_MOSES86_0.001.sh
# 

WORKDIR=//nameofdirectorythatcontainscodeanddatasubdirectories
EXPNAME=MOSES86_0.001
MONMEANSBEGINYEAR=1471
MONMEANSFINALYEAR=1500
VEGACCEL=1
NCVEG=1
SETNUM=1
TOTALSETNUM=50
SIMYEAR=1
SY=0000
YRFROM=1981
YRTO=2010
ENDOFOLDYEARFILEEXIST=0
ENDOFOLDYEARFILE=eoy_${SY}.sra

F90EXECUTABLE=main_${EXPNAME}.out
DATADIR=$WORKDIR
SWDIR=${DATADIR}/dswr
ETDIR=${DATADIR}/ET
PRESSDIR=${DATADIR}/pres
TSFCDIR=${DATADIR}/skt
T25DIR=${DATADIR}/tsoil2
SPECHUMDIR=${DATADIR}/shum
SNOWDIR=${DATADIR}/swe
LANDSEADIR=${DATADIR}/land_sea_mask
GLACIERDIR=${DATADIR}/glacier_mask
Z0CLIMODIR=${DATADIR}/z0climo
PRECIPDIR=${DATADIR}/precip
SNOWEVAPDIR=${DATADIR}/snowevap
SNOWFALLDIR=${DATADIR}/snowfall
SNOWMELTDIR=${DATADIR}/snowmelt
T60DIR=${DATADIR}/t60
U60DIR=${DATADIR}/u60
V60DIR=${DATADIR}/v60

# note that prefixes are only for the yearly .sra data files
SWPREFIX=dswr_
ETPREFIX=ET_
PRESSPREFIX=pres_
TSFCPREFIX=skt_
T25PREFIX=tsoil2_
SPECHUMPREFIX=shum_
SNOWPREFIX=swe_
LANDSEAPREFIX=landsea_t62
GLACIERPREFIX=glacier_t62
Z0CLIMOPREFIX=z0climo_t62
PRECIPPREFIX=precip_
SNOWEVAPPREFIX=snowevap_
SNOWFALLPREFIX=snowfall_
SNOWMELTPREFIX=snowmelt_
T60PREFIX=t60_
U60PREFIX=u60_
V60PREFIX=v60_

cd $WORKDIR
mkdir -p ${WORKDIR}/${EXPNAME}/output
OUTPUTDIR=${WORKDIR}/${EXPNAME}/output

cd ${WORKDIR}/${EXPNAME}
mkdir -p ${WORKDIR}/${EXPNAME}/output

if [ -e $ENDOFOLDYEARFILE ]; then
  ENDOFOLDYEARFILEEXIST=1
fi

while [ $SETNUM -le $TOTALSETNUM ]; do

YEAR=$YRFROM

  while [ $YEAR -le $YRTO ]; do

  if [ $VEGACCEL -gt 0 ]; then
    if [ $SIMYEAR -le 261 ]; then NCVEG=1; fi
    if [ $SIMYEAR -le 260 ]; then NCVEG=2; fi
    if [ $SIMYEAR -le 190 ]; then NCVEG=5; fi
    if [ $SIMYEAR -le 120 ]; then NCVEG=10; fi
    if [ $SIMYEAR -le 60 ]; then NCVEG=20; fi
    if [ $SIMYEAR -le 3 ]; then NCVEG=1; fi
  fi

  SY=$SIMYEAR
  
  if [ $SIMYEAR -lt 10 ]; then SY="0"$SY; fi
  if [ $SIMYEAR -lt 100 ]; then SY="0"$SY; fi
  if [ $SIMYEAR -lt 1000 ]; then SY="0"$SY; fi
    
  ENDOFYEARFILE=eoy_${SY}.sra
  WATCMONMEANSFILE=monmeans_watc_${SY}.sra
  ROMONMEANSFILE=monmeans_ro_${SY}.sra
  Z0MONMEANSFILE=monmeans_z0_${SY}.sra
  ALBMONMEANSFILE=monmeans_alb_${SY}.sra
  VEGMONMEANSFILE=monmeans_veg_${SY}.sra
  LAIMONMEANSFILE=monmeans_lai_${SY}.sra
  FORESTMONMEANSFILE=monmeans_forest_${SY}.sra
  WMAXMONMEANSFILE=monmeans_wmax_${SY}.sra
  NPPMONMEANSFILE=monmeans_npp_${SY}.sra
  GPPLMONMEANSFILE=monmeans_gppl_${SY}.sra
  GPPWMONMEANSFILE=monmeans_gppw_${SY}.sra
  CVEGMONMEANSFILE=monmeans_cveg_${SY}.sra
  CSOILMONMEANSFILE=monmeans_csoil_${SY}.sra
  RESMONMEANSFILE=monmeans_res_${SY}.sra
  WATERMONMEANSFILE=monmeans_water_${SY}.sra
  RHSMONMEANSFILE=monmeans_rhs_${SY}.sra
  GAMONMEANSFILE=monmeans_ga_${SY}.sra
  RCMONMEANSFILE=monmeans_rc_${SY}.sra
  RCMINMONMEANSFILE=monmeans_rcmin_${SY}.sra
  WSF_WUEMONMEANSFILE=monmeans_wsf_wue_${SY}.sra
  PGSVMONMEANSFILE=monmeans_pgsv_${SY}.sra
  ETOUTMONMEANSFILE=monmeans_etout_${SY}.sra
  VFTMONMEANSFILE=monmeans_vft_${SY}.sra
  TRANSPMONMEANSFILE=monmeans_transp_${SY}.sra

  YRSTRING=$YEAR
  
    cat > year_namelist << EOF
&yearpar
 SIMYEAR=${SIMYEAR}
 YEAR=${YEAR}
 YRSTRING='${YRSTRING}'
 NCVEG=${NCVEG}
 ENDOFOLDYEARFILEEXIST=${ENDOFOLDYEARFILEEXIST}
 ENDOFOLDYEARFILE='${ENDOFOLDYEARFILE}'
 ENDOFYEARFILE='${ENDOFYEARFILE}' 
 WATCMONMEANSFILE='${WATCMONMEANSFILE}'
 ROMONMEANSFILE='${ROMONMEANSFILE}'
 Z0MONMEANSFILE='${Z0MONMEANSFILE}'
 ALBMONMEANSFILE='${ALBMONMEANSFILE}'
 VEGMONMEANSFILE='${VEGMONMEANSFILE}'
 LAIMONMEANSFILE='${LAIMONMEANSFILE}'
 FORESTMONMEANSFILE='${FORESTMONMEANSFILE}'
 WMAXMONMEANSFILE='${WMAXMONMEANSFILE}'
 NPPMONMEANSFILE='${NPPMONMEANSFILE}'
 GPPLMONMEANSFILE='${GPPLMONMEANSFILE}'
 GPPWMONMEANSFILE='${GPPWMONMEANSFILE}'
 CVEGMONMEANSFILE='${CVEGMONMEANSFILE}'
 CSOILMONMEANSFILE='${CSOILMONMEANSFILE}'
 RESMONMEANSFILE='${RESMONMEANSFILE}'
 WATERMONMEANSFILE='${WATERMONMEANSFILE}'
 RHSMONMEANSFILE='${RHSMONMEANSFILE}'
 GAMONMEANSFILE='${GAMONMEANSFILE}'
 RCMONMEANSFILE='${RCMONMEANSFILE}'
 RCMINMONMEANSFILE='${RCMINMONMEANSFILE}'
 WSF_WUEMONMEANSFILE='${WSF_WUEMONMEANSFILE}'
 PGSVMONMEANSFILE='${PGSVMONMEANSFILE}'
 ETOUTMONMEANSFILE='${ETOUTMONMEANSFILE}'
 VFTMONMEANSFILE='${VFTMONMEANSFILE}'
 TRANSPMONMEANSFILE='${TRANSPMONMEANSFILE}'
/
EOF

    cat > data2_namelist << EOF
&data2par
 SWDIR='${SWDIR}'
 ETDIR='${ETDIR}'
 PRESSDIR='${PRESSDIR}'
 TSFCDIR='${TSFCDIR}'
 T25DIR='${T25DIR}'
 SPECHUMDIR='${SPECHUMDIR}'
 SNOWDIR='${SNOWDIR}'
 LANDSEADIR='${LANDSEADIR}'
 GLACIERDIR='${GLACIERDIR}'
 Z0CLIMODIR='${Z0CLIMODIR}'
 PRECIPDIR='${PRECIPDIR}'
 SNOWEVAPDIR='${SNOWEVAPDIR}'
 SNOWFALLDIR='${SNOWFALLDIR}'
 SNOWMELTDIR='${SNOWMELTDIR}'
 T60DIR='${T60DIR}'
 U60DIR='${U60DIR}'
 V60DIR='${V60DIR}'
 SWPREFIX='${SWPREFIX}'
 ETPREFIX='${ETPREFIX}'
 PRESSPREFIX='${PRESSPREFIX}'
 TSFCPREFIX='${TSFCPREFIX}'
 T25PREFIX='${T25PREFIX}'
 SPECHUMPREFIX='${SPECHUMPREFIX}'
 SNOWPREFIX='${SNOWPREFIX}'
 LANDSEAPREFIX='${LANDSEAPREFIX}'
 GLACIERPREFIX='${GLACIERPREFIX}'
 Z0CLIMOPREFIX='${Z0CLIMOPREFIX}'
 PRECIPPREFIX='${PRECIPPREFIX}'
 SNOWEVAPPREFIX='${SNOWEVAPPREFIX}'
 SNOWFALLPREFIX='${SNOWFALLPREFIX}'
 SNOWMELTPREFIX='${SNOWMELTPREFIX}'
 T60PREFIX='${T60PREFIX}'
 U60PREFIX='${U60PREFIX}'
 V60PREFIX='${V60PREFIX}'
/
EOF

    ${WORKDIR}/$F90EXECUTABLE

	# convert to srv format, set year to simulation year, and transfer output data to output directory
	cdo -f srv setyear,$SIMYEAR -inputsrv eoy_${SY}.srv < $ENDOFYEARFILE
	cp $ENDOFYEARFILE $OUTPUTDIR
	mv eoy_${SY}.srv $OUTPUTDIR
	
	if [ -e $ENDOFOLDYEARFILE ]; then rm $ENDOFOLDYEARFILE; fi 
	ENDOFOLDYEARFILE=$ENDOFYEARFILE
	ENDOFOLDYEARFILEEXIST=1

#cdo -f srv setyear,1200 -inputsrv eoy_${SY}.srv < eoy_1200.sra

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_watc_${SY}.srv < $WATCMONMEANSFILE
	rm $WATCMONMEANSFILE
	mv monthly_watc_${SY}.srv $OUTPUTDIR
	
	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_ro_${SY}.srv < $ROMONMEANSFILE
	rm $ROMONMEANSFILE
	mv monthly_ro_${SY}.srv $OUTPUTDIR
	
	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_z0_${SY}.srv < $Z0MONMEANSFILE
	rm $Z0MONMEANSFILE
	mv monthly_z0_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_alb_${SY}.srv < $ALBMONMEANSFILE
	rm $ALBMONMEANSFILE
	mv monthly_alb_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_veg_${SY}.srv < $VEGMONMEANSFILE
	rm $VEGMONMEANSFILE
	mv monthly_veg_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_lai_${SY}.srv < $LAIMONMEANSFILE
	rm $LAIMONMEANSFILE
	mv monthly_lai_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_forest_${SY}.srv < $FORESTMONMEANSFILE
	rm $FORESTMONMEANSFILE
	mv monthly_forest_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_wmax_${SY}.srv < $WMAXMONMEANSFILE
	rm $WMAXMONMEANSFILE
	mv monthly_wmax_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_npp_${SY}.srv < $NPPMONMEANSFILE
	rm $NPPMONMEANSFILE
	mv monthly_npp_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_gppl_${SY}.srv < $GPPLMONMEANSFILE
	rm $GPPLMONMEANSFILE
	mv monthly_gppl_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_gppw_${SY}.srv < $GPPWMONMEANSFILE
	rm $GPPWMONMEANSFILE
	mv monthly_gppw_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_cveg_${SY}.srv < $CVEGMONMEANSFILE
	rm $CVEGMONMEANSFILE
	mv monthly_cveg_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_csoil_${SY}.srv < $CSOILMONMEANSFILE
	rm $CSOILMONMEANSFILE
	mv monthly_csoil_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_res_${SY}.srv < $RESMONMEANSFILE
	rm $RESMONMEANSFILE
	mv monthly_res_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_water_${SY}.srv < $WATERMONMEANSFILE
	rm $WATERMONMEANSFILE
	mv monthly_water_${SY}.srv $OUTPUTDIR

	cdo -f srv setyear,$SIMYEAR -inputsrv monthly_rhs_${SY}.srv < $RHSMONMEANSFILE
	rm $RHSMONMEANSFILE
	mv monthly_rhs_${SY}.srv $OUTPUTDIR

    cdo -f srv setyear,$SIMYEAR -inputsrv monthly_ga_${SY}.srv < $GAMONMEANSFILE
    rm $GAMONMEANSFILE
    mv monthly_ga_${SY}.srv $OUTPUTDIR

    cdo -f srv setyear,$SIMYEAR -inputsrv monthly_rc_${SY}.srv < $RCMONMEANSFILE
    rm $RCMONMEANSFILE
    mv monthly_rc_${SY}.srv $OUTPUTDIR

    cdo -f srv setyear,$SIMYEAR -inputsrv monthly_rcmin_${SY}.srv < $RCMINMONMEANSFILE
    rm $RCMINMONMEANSFILE
    mv monthly_rcmin_${SY}.srv $OUTPUTDIR

    cdo -f srv setyear,$SIMYEAR -inputsrv monthly_wsf_wue_${SY}.srv < $WSF_WUEMONMEANSFILE
    rm $WSF_WUEMONMEANSFILE
    mv monthly_wsf_wue_${SY}.srv $OUTPUTDIR

    cdo -f srv setyear,$SIMYEAR -inputsrv monthly_pgsv_${SY}.srv < $PGSVMONMEANSFILE
    rm $PGSVMONMEANSFILE
    mv monthly_pgsv_${SY}.srv $OUTPUTDIR

    cdo -f srv setyear,$SIMYEAR -inputsrv monthly_etout_${SY}.srv < $ETOUTMONMEANSFILE
    rm $ETOUTMONMEANSFILE
    mv monthly_etout_${SY}.srv $OUTPUTDIR

    cdo -f srv setyear,$SIMYEAR -inputsrv monthly_transp_${SY}.srv < $TRANSPMONMEANSFILE
    rm $TRANSPMONMEANSFILE
    mv monthly_transp_${SY}.srv $OUTPUTDIR

    cdo -f srv setyear,$SIMYEAR -inputsrv monthly_vft_${SY}.srv < $VFTMONMEANSFILE
    rm $VFTMONMEANSFILE
    mv monthly_vft_${SY}.srv $OUTPUTDIR

    YEAR=`expr $YEAR + 1`
    SIMYEAR=`expr $SIMYEAR + 1`

  done

   SETNUM=`expr $SETNUM + 1`

done

rm -f $ENDOFYEARFILE

# use cdo to create monthly mean files over all years
for VAR in watc ro z0 alb veg lai forest wmax npp gppl gppw cveg csoil res water rhs\
 ga rc rcmin pgsv etout transp vft ; do
SIMYEAR=$MONMEANSBEGINYEAR
touch ${OUTPUTDIR}/monmeans_${VAR}_${MONMEANSBEGINYEAR}_${MONMEANSFINALYEAR}.srv
MONMEANSFILE=${OUTPUTDIR}/monmeans_${VAR}_${MONMEANSBEGINYEAR}_${MONMEANSFINALYEAR}.srv
TEMPFILE1=${OUTPUTDIR}/monmeans_${VAR}_temp1.srv
TEMPFILE2=${OUTPUTDIR}/monmeans_${VAR}_temp2.srv
while [ $SIMYEAR -le $MONMEANSFINALYEAR ]; do
SY=$SIMYEAR
if [ $SIMYEAR -lt 10 ]; then SY="0"$SY; fi
if [ $SIMYEAR -lt 100 ]; then SY="0"$SY; fi
if [ $SIMYEAR -lt 1000 ]; then SY="0"$SY; fi
if [ -e $TEMPFILE1 ]; then
cdo mergetime ${OUTPUTDIR}/monthly_${VAR}_${SY}.srv $TEMPFILE1 $TEMPFILE2
else
cp ${OUTPUTDIR}/monthly_${VAR}_${SY}.srv $TEMPFILE2
fi
mv $TEMPFILE2 $TEMPFILE1
SIMYEAR=`expr $SIMYEAR + 1`
done
mv $TEMPFILE1 $MONMEANSFILE
done

exit
