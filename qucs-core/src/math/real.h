/*
 * real.h - some real valued function definitions
 *
 * Copyright (C) 2008 Stefan Jahn <stefan@lkcc.org>
 * Copyright (C) 2014 Guilheme Brondani Torri <guitorri@gmail.com>
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this package; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * $Id$
 *
 */

#ifndef __REAL_H__
#define __REAL_H__

#include <cmath>

/**It is prefered to add all used funcions into the qucs namespace.
 * Doing so one is forced do think about compatibility instead of using std directly.
 * Inline is optional at this moment
 * \todo test if inline indeed performace improves (optimization flags should inline them anyway)
 */

namespace qucs {

//
// trigonometric
//
/*! \brief Compute cosine of an angle
    \param[in] z angle in radians
    \return cosine of z
*/
inline nr_double_t    cos (const nr_double_t arg) {
  return std::cos(arg);
}
/*! \brief Compute sine of an angle
    \param[in] z angle in radians
    \return sine of z
*/
inline nr_double_t    sin (const nr_double_t arg) {
   return std::sin (arg);
};

/*! \brief Compute tangent of an angle
    \param[in] z angle in radians
    \return tangent of z
*/
inline nr_double_t  tan (const nr_double_t arg) {
  return std::tan (arg);
}

/*! \brief Compute arc cosine
    \param[in] z arc
    \return arc cosine of z
*/
inline nr_double_t  acos (const nr_double_t arg) {
  return std::acos (arg);
}

/*! \brief Compute arc sine
    \param[in] z arc
    \return arc sine of z
*/
inline nr_double_t  asin (const nr_double_t arg) {
  return std::asin (arg);
}

/*! \brief Compute arc tangent
    \param[in] z arc
    \return arc tangent of z
*/
inline nr_double_t  atan (const nr_double_t arg) {
  return std::atan (arg);
}

/*! \brief Compute arc tangent with two parameters (fortran like function)
    \param[in] x proportion of x-coordinate
    \param[in] y proportion of y-coordinate
    \return principal value of the arc tangent of y/x, expressed in radians.
*/
inline nr_double_t  atan2 (const nr_double_t x, const nr_double_t y) {
  return std::atan2 (x,y);
}

//
// hyperbolic
//
/*! \brief Compute hyperbolic cosine
    \param[in] z arc
    \return hyperbolic cosine of z
*/
inline nr_double_t  cosh (const nr_double_t arg) {
  return std::cosh (arg);
}

/*! \brief Compute hyperbolic sine
    \param[in] z arc
    \return hyperbolic sine of z
*/
inline nr_double_t  sinh (const nr_double_t arg) {
  return std::sinh (arg);
}

/*! \brief Compute hyperbolic tangent
    \param[in] z arc
    \return hyperbolic tangent of z
*/
inline nr_double_t  tanh (const nr_double_t arg) {
  return std::tanh (arg);
}

#ifdef HAVE_CXX_STD_NUMERICAL_ACOSH
inline nr_double_t  acosh (const nr_double_t arg) {
  return std::acosh (arg);
}
#else
nr_double_t  acosh (const nr_double_t);
#endif

#ifdef HAVE_CXX_STD_NUMERICAL_ASINH
inline nr_double_t  asinh (const nr_double_t arg) {
  return std::asinh (arg);
}
#else
 nr_double_t  asinh (const nr_double_t arg);
 #endif

#ifdef HAVE_CXX_STD_NUMERICAL_ATANH
inline nr_double_t  atanh (const nr_double_t arg) {
  return std::atanh (arg);
}
#else
 nr_double_t  atanh (const nr_double_t arg);
#endif


//
// exponential and logarithmic functions
//
inline nr_double_t exp (const nr_double_t arg) {
  return std::exp (arg);
}
inline nr_double_t log (const nr_double_t arg) {
   return std::log(arg);
}
inline nr_double_t log10 (const nr_double_t arg) {
  return std::log10(arg);
}


//
// power functions
//
inline nr_double_t pow (const nr_double_t a, const nr_double_t b) {
  return std::pow(a,b);
}
inline nr_double_t sqrt (const nr_double_t d) {
  return std::sqrt (d);
}

#ifdef HAVE_CXX_STD_NUMERICAL_HYPOT
inline nr_double_t xhypot (const nr_double_t a, const nr_double_t b) {
  return std::hypot(a,b);
}
#else
 nr_double_t xhypot (const nr_double_t, const nr_double_t );
#endif

//
// error functions
//
#ifdef HAVE_CXX_STD_NUMERICAL_ERF
 inline nr_double_t erf(const nr_double_t arg) {
   return std::erf (arg);
}
#else
nr_double_t erf(const nr_double_t arg);
#endif


//
// rounding and remainder functions
//
inline nr_double_t ceil(const nr_double_t arg) {
  return std::ceil(arg);
}
inline nr_double_t floor(const nr_double_t arg) {
  return std::floor(arg);
}
#ifdef HAVE_CXX_STD_NUMERICAL_FMOD
inline nr_double_t fmod(const nr_double_t arg1,const nr_double_t arg2) {
  return std::fmod(arg1,arg2);
}
#else
nr_double_t fmod(const nr_double_t arg1,const nr_double_t arg2);
#endif

#ifdef HAVE_CXX_STD_NUMERICAL_TRUNC
 inline nr_double_t trunc(const nr_double_t arg) {
   return std::trunc(arg);
 }
#else
 nr_double_t trunc(const nr_double_t arg);
#endif
#ifdef HAVE_CXX_STD_NUMERICAL_ROUND
 inline nr_double_t round(const nr_double_t arg) {
   return std::round(arg);
 }
#else
 nr_double_t round(const nr_double_t arg);
#endif
//
// Qucs extra trigonometric helper
//
nr_double_t coth (const nr_double_t );
nr_double_t sech (const nr_double_t );
nr_double_t cosech (const nr_double_t );


//
// Qucs extra math functions
//
nr_double_t  sqr (const nr_double_t );

unsigned int sqr (unsigned int);


nr_double_t  quadr (const nr_double_t );


//
// extra math functions
//
nr_double_t limexp (const nr_double_t);
nr_double_t signum (const nr_double_t);
nr_double_t   sign (const nr_double_t);
nr_double_t   sinc (const nr_double_t);
nr_double_t    fix (const nr_double_t);
nr_double_t   step (const nr_double_t);
unsigned int factorial (unsigned int);


//
// overload complex manipulations on reals
//
 inline nr_double_t   real (const nr_double_t d) { return d; };
 inline nr_double_t   imag (const nr_double_t d) { (void) d;return 0.0; };
nr_double_t   norm (const nr_double_t);
inline nr_double_t   conj (const nr_double_t d) { return d; };
 inline nr_double_t   abs (const nr_double_t d) { return std::abs(d); }

} // namespace qucs

#endif /* __REAL_H__ */
