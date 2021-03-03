#include <cstdlib>
#include <cmath>
#include "cvm/geodetic.hpp"
/*
!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!=====================================================================
!
!  UTM (Universal Transverse Mercator) projection from the USGS
!
!=====================================================================

! taken from the open-source package CAMx at http://www.camx.com/download
! converted by Dimitri Komatitsch to Fortran90 and slightly modified to add one parameter to the subroutine call
! and change the order of the arguments for compatibility with SPECFEM calls.
! Also converted to double precision for all calculations and all results.
! Also defined the UTM easting and northing in meters instead of kilometers for compatibility with SPECFEM calls.

! convert geodetic longitude and latitude to UTM, and back
! use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
! a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm

*/
namespace
{
void utm_geo(double *rlon4, double *rlat4,
             double *rx4,   double *ry4,
             const int UTM_PROJECTION_ZONE, const int iway,
             const bool SUPPRESS_UTM_PROJECTION)
{
/*
!
!---- CAMx v6.10 2014/04/02
!>!
!>    @brief
!>    UTM_GEO performs UTM to geodetic (long/lat) translation, and back.
!>
!>    This is a Fortran version of the BASIC program "Transverse Mercator
!>    Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
!>    Based on algorithm taken from "Map Projections Used by the USGS"
!>    by John P. Snyder, Geological Survey Bulletin 1532, USDI.
!>
!>    Portions Copyright 1996 - 2014
!>    ENVIRON International Corporation
!>
!>    Modifications:
!>     2012/12/02   Added logic for southern hemisphere UTM zones
!>                  Use zones +1 to +60 for the Northern hemisphere, -1 to -60 for the 
!>                  Southern hemisphere
!>                  Equator is defined as 0 km North for the Northern hemisphere, 10,000
!>                  km North for the Southern hemisphere
!>
!>     2016/03/24   Made it C compatible, moved hard to find constants inside,
!>                  introduced handling for doxygen, and ensured code does not 
!>                  go past column 90
!
!     Input/Output arguments:
!
!>    @param[inout] rlon4              Longitude (degrees, negative for West)
!>    @param[inout] rlat4              Latitude (degrees)
!>    @param[inout] rx4                UTM easting (meters)
!>    @param[inout] ry4                UTM northing (meters)
!>    @param[in] UTM_PROJECTION_ZONE   UTM zone
!>                                      The Northern hemisphere corresponds to zones +1 to +60
!>                                      The Southern hemisphere corresponds to zones -1 to -60
!>    @param[in] iway                 Conversion type
!>                                      ILONGLAT2UTM (0) = geodetic to UTM
!>                                      IUTM2LONGLAT (1) = UTM to geodetic
!>
!>    @reference Some general information about UTM:
!>               (for more details see e.g. http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system )
!
!>    @note There are 60 longitudinal projection zones numbered 1 to 60 starting at 180 degrees W.
!>          Each of these zones is 6 degrees wide, apart from a few exceptions around Norway and Svalbard.
!>          There are 20 latitudinal zones spanning the latitudes 80 degrees S to 84 degrees N and denoted
!>          by the letters C to X, ommitting the letter O.
!>          Each of these is 8 degrees south-north, apart from zone X which is 12 degrees south-north.
!>
!>          The UTM zone is described by the central meridian of that zone, i.e. the longitude at the
!>          midpoint of the zone, 3 degrees away from both zone boundary.
!>
!>    author Dimitri Komatitsch and Jeroen Tromp
!>
  USE ISO_C_BINDING
  !use constants, only: PI,ILONGLAT2UTM,IUTM2LONGLAT

  implicit none

! input/output parameters
  REAL(C_DOUBLE), intent(inout) :: rx4,ry4,rlon4,rlat4
  INTEGER(C_INT), intent(in) :: UTM_PROJECTION_ZONE,iway
  LOGICAL(C_BOOL), intent(in) :: SUPPRESS_UTM_PROJECTION

! local parameters
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793d0
  INTEGER, PARAMETER :: ILONGLAT2UTM = 0
  INTEGER, PARAMETER :: IUTM2LONGLAT = 1
*/
  const int ILONGLAT2UTM = 0;
  const int IUTM2LONGLAT = 1; 

/*
! From http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system :
! The Universal Transverse Mercator coordinate system was developed by the United States Army Corps of Engineers in the 1940s.
! The system was based on an ellipsoidal model of Earth. For areas within the contiguous United States
! the Clarke Ellipsoid of 1866 was used. For the remaining areas of Earth, including Hawaii, the International Ellipsoid was used.
! The WGS84 ellipsoid is now generally used to model the Earth in the UTM coordinate system,
! which means that current UTM northing at a given point can be 200+ meters different from the old one.
! For different geographic regions, other datum systems (e.g.: ED50, NAD83) can be used.

! Clarke 1866
! double precision, parameter :: SEMI_MAJOR_AXIS = 6378206.4d0, SEMI_MINOR_AXIS = 6356583.8d0

! WGS84 (World Geodetic System 1984)
  double precision, parameter :: SEMI_MAJOR_AXIS = 6378137.0d0, &
                                 SEMI_MINOR_AXIS = 6356752.314245d0
*/
  const double SEMI_MAJOR_AXIS = 6378137.0;
  const double SEMI_MINOR_AXIS = 6356752.314245;

/*
! Note that the UTM grids are actually Mercators which
! employ the standard UTM scale factor 0.9996 and set the Easting Origin to 500,000.
  double precision, parameter :: scfa=0.9996d0
  double precision, parameter :: north=0.d0, east=500000.d0

  double precision, parameter :: DEGREES_TO_RADIANS=0.017453292519943295d0
  double precision, parameter :: RADIANS_TO_DEGREES=57.29577951308232d0
*/
   const double DEGREES_TO_RADIANS=0.017453292519943295;
   const double RADIANS_TO_DEGREES=57.29577951308232;
   const double scfa=0.9996;
   const double north=0;
   const double east=500000.;

//! local variables
  int  zone;
  double rlon,rlat;
  double e2,e4,e6,ep2,xx,yy,dlat,dlon,cm,cmr,delam;
  double f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d;
  double rx_save,ry_save,rlon_save,rlat_save;
  bool lsouth;

  //! checks if conversion to UTM has to be done
  if (SUPPRESS_UTM_PROJECTION)
  {
    if (iway == ILONGLAT2UTM)
    {
      *rx4 = *rlon4;
      *ry4 = *rlat4;
    }
    else
    {
      *rlon4 = *rx4;
      *rlat4 = *ry4;
    } 
    return;
  }

  //! save original parameters
  rlon_save = *rlon4;
  rlat_save = *rlat4;
  rx_save = *rx4;
  ry_save = *ry4;

  e2=1.0-std::pow(SEMI_MINOR_AXIS/SEMI_MAJOR_AXIS, 2);
  e4=e2*e2;
  e6=e2*e4;
  ep2=e2/(1.0-e2);

//!
//!---- Set Zone parameters
//!

  lsouth = false;
  if (UTM_PROJECTION_ZONE < 0){lsouth = true;}
  zone = std::abs(UTM_PROJECTION_ZONE);
  cm = zone*6.0 - 183.0;
  cmr = cm*DEGREES_TO_RADIANS;

  dlat = 0.0; // ! removes an uninitialized warning - baker
/*
  ! this code block does nothing since i moved the variables down - baker
! if (iway == IUTM2LONGLAT) then
!   xx = rx4
!   yy = ry4
!   if (lsouth) yy = yy - 1.d7
! else
!   dlat = rlat4
!   dlon = rlon4
! endif

!
!---- Lat/Lon to UTM conversion
!
*/
  if (iway == ILONGLAT2UTM){

    dlat = *rlat4;// ! removes an uninitialized warning - baker
    dlon = *rlon4;// ! removes an uninitialized warning - baker

    rlon = DEGREES_TO_RADIANS*dlon;
    rlat = DEGREES_TO_RADIANS*dlat;

    delam = dlon - cm;
    if (delam < -180.0){delam = delam + 360.0;}
    if (delam > 180.0){delam = delam - 360.0;}
    delam = delam*DEGREES_TO_RADIANS;

    f1 = (1.0 - e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0)*rlat;
    f2 = 3.0*e2/8.0 + 3.0*e4/32.0 + 45.0*e6/1024.0;
    f2 = f2*std::sin(2.0*rlat);
    f3 = 15.0*e4/256.0*45.0*e6/1024.0;
    f3 = f3*std::sin(4.0*rlat);
    f4 = 35.0*e6/3072.0;
    f4 = f4*std::sin(6.0*rlat);
    rm = SEMI_MAJOR_AXIS*(f1 - f2 + f3 - f4);
    if (std::abs(dlat - 90.0) < 1.e-12 || std::abs(dlat + 90.0) < 1.e-12)
    {
      xx = 0.0;
      yy = scfa*rm;
    }else{
      rn = SEMI_MAJOR_AXIS/std::sqrt(1.0 - e2*std::pow(sin(rlat),2));
      t = std::pow(std::tan(rlat),2);//**2
      c = ep2*std::pow(std::cos(rlat),2);//**2
      a = std::cos(rlat)*delam;

      f1 = (1.0 - t + c)*std::pow(a,3)/6.0;
      f2 = 5.0 - 18.0*t + std::pow(t,2) + 72.0*c - 58.0*ep2;
      f2 = f2*std::pow(a,5)/120.0;
      xx = scfa*rn*(a + f1 + f2);
      f1 = std::pow(a,2)/2.0;
      f2 = 5.0 - t + 9.0*c + 4.0*std::pow(c,2);
      f2 = f2*std::pow(a,4)/24.0;
      f3 = 61.0 - 58.0*t + std::pow(t,2) + 600.0*c - 330.0*ep2;
      f3 = f3*std::pow(a,6)/720.0;
      yy = scfa*(rm + rn*std::tan(rlat)*(f1 + f2 + f3));
    } 
    xx = xx + east;
    yy = yy + north;

//!
//!---- UTM to Lat/Lon conversion
//!
  }else{

    xx = *rx4;//  ! removes uninitialized warning - baker
    yy = *ry4;//  ! removes uninititalted warning - baker
    if (lsouth){yy = yy - 1.e7;}// ! removes uninitialized warning - baker

    xx = xx - east;
    yy = yy - north;
    e1 = std::sqrt(1.0 - e2);
    e1 = (1.0 - e1)/(1.0 + e1);
    rm = yy/scfa;
    u = 1.0 - e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0;
    u = rm/(SEMI_MAJOR_AXIS*u);

    f1 = 3.0*e1/2.0 - 27.0*std::pow(e1,3.0)/32.0;
    f1 = f1*std::sin(2.0*u);
    f2 = 21.0*std::pow(e1,2)/16.0 - 55.0*std::pow(e1,4)/32.0;
    f2 = f2*std::sin(4.0*u);
    f3 = 151.0*pow(e1,3.0)/96.0;
    f3 = f3*std::sin(6.0*u);
    rlat1 = u + f1 + f2 + f3;
    dlat1 = rlat1*RADIANS_TO_DEGREES;
    if (dlat1 >= 90.0 || dlat1 <= -90.0){
      dlat1 = std::fmin(dlat1,90.0);
      dlat1 = std::fmax(dlat1,-90.0);
      dlon = cm;
    }else{
      c1 = ep2*std::pow(cos(rlat1),2);
      t1 = std::pow(std::tan(rlat1),2);
      f1 = 1.0 - e2*std::pow(std::sin(rlat1),2);
      rn1 = SEMI_MAJOR_AXIS/sqrt(f1);
      r1 = SEMI_MAJOR_AXIS*(1.0 - e2)/std::sqrt(std::pow(f1,3));
      d = xx/(rn1*scfa);

      f1 = rn1*std::tan(rlat1)/r1;
      f2 = std::pow(d,2)/2.0;
      f3 = 5.0*3.0*t1 + 10.0*c1 - 4.0*std::pow(c1,2) - 9.0*ep2;
      f3 = f3*std::pow(d,2)*std::pow(d,2)/24.0;
      f4 = 61.0 + 90.0*t1 + 298.0*c1 + 45.0*std::pow(t1,2) - 252.0*ep2 - 3.0*std::pow(c1,2);
      f4 = f4*std::pow((std::pow(d,2)),3.0)/720.0;
      rlat = rlat1 - f1*(f2 - f3 + f4);
      dlat = rlat*RADIANS_TO_DEGREES;

      f1 = 1.0 + 2.0*t1 + c1;
      f1 = f1*std::pow(d,2)*d/6.;
      f2 = 5.0 - 2.0*c1 + 28.0*t1 - 3.0*std::pow(c1,2) + 8.0*ep2 + 24.0*std::pow(t1,2);
      f2 = f2*std::pow((std::pow(d,2)),2)*d/120.0;
      rlon = cmr + (d - f1 + f2)/std::cos(rlat1);
      dlon = rlon*RADIANS_TO_DEGREES;
      if (dlon < -180.0){dlon = dlon + 360.0;}
      if (dlon > 180.0){dlon = dlon - 360.0;}
    } 
  } 

//!
//!----- output
//!
  if (iway == IUTM2LONGLAT){
    *rlon4 = dlon;
    *rlat4 = dlat;
    *rx4 = rx_save;
    *ry4 = ry_save;
  }else{
    *rx4 = xx;
    if (lsouth){yy = yy + 1.e7;}
    *ry4 = yy;
    *rlon4 = rlon_save;
    *rlat4 = rlat_save;
  } 
  return;
}
}
//  end subroutine utm_geo
/*
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Converts latitudes/longitudes to UTMs
!>
!>    @param[in] rlat4                 latitude (degrees)
!>    @param[in] rlon4                 longitude (degrees, negative for west)
!>    @param[in] UTM_PROJECTION_ZONE   UTM zone
!>                                      The Northern hemisphere corresponds to
!>                                      zones +1 to +60.
!>                                      The Southern hemisphere corresponds to
!>                                      zones -1 to -60
!>
!>    @param[out] rx4                  UTM easting (m)
!>    @param[out] ry4                  UTM northing (m)
!>
!>    @author Ben Baker, ISTI
!>
*/
std::pair<double, double> CVM::Geodetic::latitudeLongitudeToUTM(
      const std::pair<double, double> &latlon, // rlat4, const double rlon4,
      const int UTM_PROJECTION_ZONE)
{
    constexpr bool SUPPRESS_UTM_PROJECTION = false;
    const int iway = 0;// ! (lat, lon) -> UTM
    double rlon4_use, rlat4_use;
    double rx4, ry4;
    rlon4_use = latlon.second; //rlon4;
    rlat4_use = latlon.first; //rlat4;
    utm_geo(&rlon4_use, &rlat4_use, &rx4, &ry4, UTM_PROJECTION_ZONE, iway,
           SUPPRESS_UTM_PROJECTION);
    return std::pair(rx4, ry4);
}
/*
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Converts UTMs to latitudes/longitudes
!>
!>    @param[in] rx4                  UTM easting (m)
!>    @param[in] ry4                  UTM northing (m)
!>    @param[in] UTM_PROJECTION_ZONE  UTM zone
!>                                      The Northern hemisphere corresponds to
!>                                      zones +1 to +60.
!>                                      The Southern hemisphere corresponds to
!>                                      zones -1 to -60
!>
!>    @param[out] rlat4               corresponding latitude (degrees)
!>    @param[out] rlon4               corresponding longitude (degrees)
!>
!>    @author Ben Baker, ISTI
!>
*/
std::pair<double, double> CVM::Geodetic::utmToLatitudeLongitude(
    const std::pair<double, double> &xy, //double rx4, const double ry4,
    const  int UTM_PROJECTION_ZONE)
{
    const int iway = 1;// ! UTM -> (lat, lon)
    const bool SUPPRESS_UTM_PROJECTION = false;
    double rx4_use, ry4_use;
    double rlat4, rlon4;
    rx4_use = xy.first; //rx4;
    ry4_use = xy.second; //ry4;
    utm_geo(&rlon4, &rlat4, &rx4_use, &ry4_use, UTM_PROJECTION_ZONE, iway, 
            SUPPRESS_UTM_PROJECTION);
    return std::pair(rlat4, rlon4);
}
