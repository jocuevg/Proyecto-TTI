//$Header$
//------------------------------------------------------------------------------
// Sat_Const
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file Sat_Const.hpp
 * @brief Archivo cabecera de constantes globales.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#ifndef _SATCONST_
#define _SATCONST_
#include <cmath>

// Mathematical constants

//--------------
// const double pi 
//--------------
/** @brief Constante matemática π. */
const double pi        = M_PI;                // pi

//--------------
// const double pi2
//--------------
/** @brief Doble de la constante π (2π). */
const double pi2       = 2.0*pi;              // 2pi

//--------------
// const double Rad
//--------------
/** @brief Radianes por grado. */
const double Rad       = pi/180.0;            // Radians per degree

//--------------
// const double Deg
//--------------
/** @brief Grados por radian. */
const double Deg       = 180.0/pi;            // Degrees per radian

//--------------
// const double Arcs
//--------------
/** @brief Arcosegundos por radian. */
const double Arcs      = 3600.0*180.0/pi;     // Arcseconds per radian


// General

//--------------
// const double MJD_J2000
//--------------
/** @brief Fecha Juliana Modificada del J2000.0. */
const double MJD_J2000 = 51544.5;             // Modif. Julian Date of J2000.0

//--------------
// const double T_B1950
//--------------
/** @brief Diferencia de época respecto a B1950. */
const double T_B1950   = -0.500002108;        // Epoch B1950

//--------------
// const double c_light
//--------------
/** @brief Velocidad de la luz en el vacío [m/s]; DE200. */
const double c_light   = 299792458.0; // Speed of light  [m/s]; DE200

//--------------
// const double AU
//--------------
/** @brief Unidad astronómica [m]; DE200. */
const double AU        = 149597870700.0; // Astronomical unit [m]; DE200



// Physical parameters of the Earth, Sun and Moon

// Equatorial radius and flattening

//--------------
// const double R_Earth
//--------------
/** @brief Radio ecuatorial de la Tierra [m]; WGS-84. */
const double R_Earth   =   6378.1363e3;        // Radius Earth [m]; WGS-84

//--------------
// const double f_Earth
//--------------
/** @brief Aplanamiento de la Tierra; WGS-84. */
const double f_Earth   = 1.0/298.257223563;   // Flattening; WGS-84

//--------------
// const double R_Sun
//--------------
/** @brief Radio del Sol [m]; Seidelmann 1992. */
const double R_Sun     = 696000.0e3;          // Radius Sun [m]; Seidelmann 1992

//--------------
// const double R_Moon
//--------------
/** @brief Radio de la Luna [m]. */
const double R_Moon    =   1738.0e3;          // Radius Moon [m]


// Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)

//--------------
// const double omega_Earth
//--------------
/** @brief Velocidad angular de rotación de la Tierra [rad/s]; Aoki 1982, NIMA 1997. */
const double omega_Earth = 15.04106717866910/3600.0*Rad;   // [rad/s]; Aoki 1982, NIMA 1997


// Gravitational coefficients

//--------------
// const double GM_Earth
//--------------
/** @brief Constante gravitacional geocéntrica de la Tierra [m³/s²]; WGS-84. */
const double GM_Earth    = 398600.435436e9;                // [m^3/s^2]; WGS-84

//--------------
// const double GM_Sun
//--------------
/** @brief Constante gravitacional heliocéntrica del Sol [m³/s²]; DE200. */
const double GM_Sun      = 132712440041.939400e9;         // [m^3/s^2]; DE200

//--------------
// const double GM_Moon
//--------------
/** @brief Constante gravitacional de la Luna [m³/s²]; DE200. */
const double GM_Moon     = GM_Earth/81.30056907419062; // [m^3/s^2]; DE200

//--------------
// const double GM_Mercury
//--------------
/** @brief Constante gravitacional de Mercurio [m³/s²]; DE200. */
const double GM_Mercury  = 22031.780000e9;          // [m^3/s^2]; DE200

//--------------
// const double GM_Venus
//--------------
/** @brief Constante gravitacional de Venus [m³/s²]; DE200. */
const double GM_Venus    = 324858.592000e9;          // [m^3/s^2]; DE200

//--------------
// const double GM_Mars
//--------------
/** @brief Constante gravitacional de Marte [m³/s²]; DE200. */
const double GM_Mars     = 42828.375214e9;          // [m^3/s^2]; DE200

//--------------
// const double GM_Jupiter
//--------------
/** @brief Constante gravitacional de Júpiter [m³/s²]; DE200. */
const double GM_Jupiter  = 126712764.800000e9;          // [m^3/s^2]; DE200

//--------------
// const double GM_Saturn
//--------------
/** @brief Constante gravitacional de Saturno [m³/s²]; DE200. */
const double GM_Saturn   = 37940585.200000e9;          // [m^3/s^2]; DE200

//--------------
// const double GM_Uranus
//--------------
/** @brief Constante gravitacional de Urano [m³/s²]; DE200. */
const double GM_Uranus   = 5794548.600000e9;          // [m^3/s^2]; DE200

//--------------
// const double GM_Neptune
//--------------
/** @brief Constante gravitacional de Neptuno [m³/s²]; DE200. */
const double GM_Neptune  = 6836527.100580e9;          // [m^3/s^2]; DE200

//--------------
// const double GM_Pluto
//--------------
/** @brief Constante gravitacional de Plutón [m³/s²]; DE200. */
const double GM_Pluto    = 977.0000000000009e9;          // [m^3/s^2]; DE200
       

// Solar radiation pressure at 1 AU

//--------------
// const double P_Sol
//--------------
/** @brief Presión de radiación solar a 1 UA [N/m²]; (~1367 W/m²); IERS 96. */
const double P_Sol       = 1367/c_light;          // [N/m^2] (~1367 W/m^2); IERS 96


#endif
