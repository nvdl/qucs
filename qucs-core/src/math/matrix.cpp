/*
 * matrix.cpp - matrix class implementation
 *
 * Copyright (C) 2003-2009 Stefan Jahn <stefan@lkcc.org>
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
/*!\file matrix.cpp
   \brief Dense matrix class implementation

  References:

  [1] Power Waves and the Scattering Matrix
      Kurokawa, K.
      Microwave Theory and Techniques, IEEE Transactions on,
      Vol.13, Iss.2, Mar 1965
      Pages: 194- 202

  [2] A Rigorous Technique for Measuring the Scattering Matrix of
      a Multiport Device with a 2-Port Network Analyzer
      John C. TIPPET, Ross A. SPECIALE
      Microwave Theory and Techniques, IEEE Transactions on,
      Vol.82, Iss.5, May 1982
      Pages: 661- 666

  [3] Comments on "A Rigorous Techique for Measuring the Scattering
      Matrix of a Multiport Device with a Two-Port Network Analyzer"
      Dropkin, H.
      Microwave Theory and Techniques, IEEE Transactions on,
      Vol. 83, Iss.1, Jan 1983
      Pages: 79 - 81

  [4] Arbitrary Impedance
      "Accurate Measurements In Almost Any
      Impedance Environment"
      in Scropion Application note
      Anritsu
      online(2007/07/30) http://www.eu.anritsu.com/files/11410-00284B.pdf

  [5] Conversions between S, Z, Y, H, ABCD, and T parameters
      which are valid for complex source and load impedances
      Frickey, D.A.
      Microwave Theory and Techniques, IEEE Transactions on
      Vol. 42, Iss. 2, Feb 1994
      pages: 205 - 211
      doi: 10.1109/22.275248

  [6] Comments on "Conversions between S, Z, Y, h, ABCD,
      and T parameters which are valid for complex source and load impedances" [and reply]
      Marks, R.B.; Williams, D.F.; Frickey, D.A.
      Microwave Theory and Techniques, IEEE Transactions on,
      Vol.43, Iss.4, Apr 1995
      Pages: 914- 915
      doi: 10.1109/22.375247

  [7] Wave Techniques for Noise Modeling and Measurement
      S. W. Wedge and D. B. Rutledge,
      IEEE Transactions on Microwave Theory and Techniques,
      vol. 40, no. 11, Nov. 1992.
      pages 2004-2012,
      doi: 10.1109/22.168757
      Author copy online (2007/07/31)
      http://authors.library.caltech.edu/6226/01/WEDieeetmtt92.pdf

*/
#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <cmath>

#include "logging.h"
#include "object.h"
#include "complex.h"
#include "vector.h"
#include "matrix.h"

namespace qucs {


/*!\brief Computes magnitude in dB of each matrix element.
   \param[in] a matrix
*/
matrix dB (matrix a) {
  matrix res (a.rows (), a.cols ());
  for (int r = 0; r < a.rows (); r++)
    for (int c = 0; c < a.cols (); c++)
      res(r, c) = dB (a(r, c));
  return res;
}

/*!\brief Computes the argument of each matrix element.
   \param[in] a matrix
   \todo add arg in place
   \todo a is const
*/
matrix arg (matrix a) {
  matrix res (a.rows (), a.cols ());
  for (int r = 0; r < a.rows (); r++)
    for (int c = 0; c < a.cols (); c++)
      res(r, c) = arg (a(r, c));
  return res;
}

/*!\brief Multiply a matrix by itself
   \param[in] a matrix
*/
matrix sqr (matrix a) {
  return a * a;
}

/*!\brief Create a diagonal matrix from a vector
   \param[in] diag vector to write on the diagonal
   \todo diag is const
*/
matrix diagonal (qucs::vector diag) {
  int size = diag.getSize ();
  matrix res (size,size);
  for (int i = 0; i < size; i++) res (i, i) = diag (i);
  return res;
}

// Compute n-th power of the given matrix.
matrix pow (matrix a, int n) {
  matrix res;
  if (n == 0) {
    res = res.Identity(a.rows(),a.cols());
  }
  else {
    if(n <0) { 
      res = a.inverse();
      a = res;
    }
    else
      res = a;
    for (int i = 1; i < std::abs (n); i++)
      res = res * a;
  }
  return res;
}



/*!\brief Converts chain matrix to scattering parameters.

    Converts scattering parameters to chain matrix.
    Formulae are given by [5] tab 1. and are remembered here:

    \f{align*}
    A&=\frac{(Z_{01}^*+S_{11}Z_{01})(1-S_{22})
                 +S_{12}S_{21}Z_{01}}{\Delta} \\
    B&=\frac{(Z_{01}^*+S_{11}Z_{01})(Z_{02}^*+S_{22}Z_{02})
                 -S_{12}S_{21}Z_{01}Z_{02}}{\Delta} \\
    C&=\frac{(1-S_{11})(1-S_{22})
                 -S_{12}S_{21}}{\Delta} \\
    D&=\frac{(1-S_{11})(Z_{02}^*+S_{22}Z_{02})
                 +S_{12}S_{21}Z_{02}}{\Delta}
    \f}
    Where:
    \f[
    \Delta = 2 S_{21}\sqrt{\Re\text{e}\;Z_{01}\Re\text{e}\;Z_{02}}
    \f]
    \bug Do not need fabs
    \param[in] s Scattering matrix
    \param[in] z1 impedance at input 1
    \param[in] z2 impedance at input 2
    \return Chain matrix
    \note Assert 2 by 2 matrix
    \todo Why not s,z1,z2 const
*/
matrix stoa (matrix s, nr_complex_t z1, nr_complex_t z2) {
  nr_complex_t d = s (0, 0) * s (1, 1) - s (0, 1) * s (1, 0);
  nr_complex_t n = 2.0 * s (1, 0) * sqrt (fabs (real (z1) * real (z2)));
  matrix a (2,2);

  assert (s.rows () >= 2 && s.cols () >= 2);

  a(0, 0)= (conj (z1) + z1 * s (0, 0) -
		conj (z1) * s (1, 1) - z1 * d) / n;
  a(0, 1) = (conj (z1) * conj (z2) + z1 * conj (z2) * s (0, 0) +
		conj (z1) * z2 * s (1, 1) + z1 * z2 * d) / n;
  a(1, 0) = (1.0 - s (0, 0) - s (1, 1) + d) / n;
  a(1, 1) =  (conj (z2) - conj (z2) * s (0, 0) +
		z2 * s (1, 1) - z2 * d) / n;
  return a;
}


/*!\brief Converts chain matrix to scattering parameters.

    Converts chain matrix to scattering parameters
    Formulae are given by [5] and are remembered here:
    \f{align*}
    S_{11}&=\frac{AZ_{02}+B-CZ_{01}^*Z_{02}-DZ_{01}^*}{\Delta} \\
    S_{12}&=\frac{2(AD-BC)
                  (\Re\text{e}\;Z_{01}\Re\text{e}\;Z_{02})^{1/2}}
                {\Delta}\\
    S_{21}&=\frac{2(\Re\text{e}\;Z_{01}\Re\text{e}\;Z_{02})^{1/2}}{\Delta}\\
    S_{22}&=\frac{-AZ_{02}^*+B-CZ_{01}^*Z_{02}+DZ_{01}}{\Delta}
    \f}
    Where:
    \f[
    \Delta =AZ_{02}+B+CZ_{01}Z_{02}-DZ_{01}
    \f]
    \param[in] a Chain matrix
    \param[in] z1 impedance at input 1
    \param[in] z2 impedance at input 2
    \return Scattering matrix
    \bug Do not use fabs
    \todo a, z1, z2 const
*/
matrix atos (matrix a, nr_complex_t z1, nr_complex_t z2) {
  nr_complex_t d = 2.0 * sqrt (fabs (real (z1) * real (z2)));
  nr_complex_t n = a (0, 0) * z2 + a (0, 1) +
    a (1, 0) * z1 * z2 + a (1, 1) * z1;
  matrix s (2,2);

  assert (a.rows () >= 2 && a.cols () >= 2);

  s(0, 0)=  (a (0, 0) * z2 + a (0, 1)
                - a (1, 0) * conj (z1) * z2 - a (1, 1) * conj (z1)) / n;
  s(0, 1) = (a (0, 0) * a (1, 1) -
		a (0, 1) * a (1, 0)) * d / n;
  s(1, 0) = d / n;
  s(1, 1) = (a (1, 1) * z1 - a (0, 0) * conj (z2) +
		a (0, 1) - a (1, 0) * z1 * conj (z2)) / n;
  return s;
}

/*!\brief Converts scattering parameters to hybrid matrix.

    Converts chain matrix to scattering parameters
    Formulae are given by [5] and are remembered here:
    \f{align*}
    h_{11}&=\frac{(Z_{01}^*+S_{11}Z_{01})
                  (Z_{02}^*+S_{22}Z_{02})
		  -S_{12}S_{21}Z_{01}Z_{02}}{\Delta}\\
    h_{12}&=\frac{2S_{12}
                  (\Re\text{e}\;Z_{01} \Re\text{e}\;Z_{02})^\frac{1}{2}}
                 {\Delta} \\
    h_{21}&=\frac{-2S_{21}
                  (\Re\text{e}\;Z_{01} \Re\text{e}\;Z_{02})^\frac{1}{2}}
                 {\Delta} \\

    h_{22}&=\frac{(1-S_{11})(1-S_{22})-S_{12}S_{21}}{\Delta}
    \f}
    Where \f$\Delta\f$ is:
    \f[
    \Delta=(1-S_{11})(Z_{02}^*+S_{22}Z_{02})+S_{12}S_{21}Z_{02}
    \f]
    \bug{Programmed formulae are valid only for Z real}
    \param[in] s Scattering matrix
    \param[in] z1 impedance at input 1
    \param[in] z2 impedance at input 2
    \return hybrid matrix
    \note Assert 2 by 2 matrix
    \todo Why not s,z1,z2 const
 */
matrix stoh (matrix s, nr_complex_t z1, nr_complex_t z2) {
  nr_complex_t n = s (0, 1) * s (1, 0);
  nr_complex_t d = (1.0 - s (0, 0)) * (1.0 + s (1, 1)) + n;
  matrix h (2,2);

  assert (s.rows () >= 2 && s.cols () >= 2);

  h(0, 0) = ((1.0 + s (0, 0)) * (1.0 + s (1, 1)) - n) * z1 / d;
  h(0, 1) = +2.0 * s (0, 1) / d;
  h(1, 0) = -2.0 * s (1, 0) / d;
  h(1, 1) = ((1.0 - s (0, 0)) * (1.0 - s (1, 1)) - n) / z2 / d;
  return h;
}

/*!\brief Converts hybrid matrix to scattering parameters.

    Formulae are given by [5] and are remembered here:
    \f{align*}
    S_{11}&=\frac{(h_{11}-Z_{01}^*)(1+h_{22}Z_{02})-h_{12}h_{21}Z_{02}}
                 {\Delta}\\
    S_{12}&=\frac{-2h_{12}(\Re\text{e}\;Z_{01}\Re\text{e}\;Z_{02})^\frac{1}{2}}
                 {\Delta}\\
    S_{21}&=\frac{-2h_{21}(\Re\text{e}\;Z_{01}\Re\text{e}\;Z_{02})^\frac{1}{2}}
                 {\Delta}\\
    S_{22}&=\frac{(h_{11}+Z_{01})(1-h_{22}Z_{02}^*)-h_{12}h_{21}Z_{02}^*}
                 {\Delta}
   \f}
   Where \f$\Delta\f$ is:
   \f[
   \Delta=(Z_{01}+h_{11})(1+h_{22}Z_{02})-h_{12}h_{21}Z_{02}
   \f]
   \param[in] h hybrid matrix
   \param[in] z1 impedance at input 1
   \param[in] z2 impedance at input 2
   \return scattering matrix
   \note Assert 2 by 2 matrix
   \todo Why not h,z1,z2 const
*/
matrix htos (matrix h, nr_complex_t z1, nr_complex_t z2) {
  nr_complex_t n = h (0, 1) * h (1, 0);
  nr_complex_t d = (1.0 + h (0, 0) / z1) * (1.0 + z2 * h (1, 1)) - n;
  matrix s (2,2);

  assert (h.rows () >= 2 && h.cols () >= 2);

  s(0, 0) = ((h (0, 0) / z1 - 1.0) * (1.0 + z2 * h (1, 1)) - n) / d;
  s(0, 1) = +2.0 * h (0, 1) / d;
  s(1, 0) = -2.0 * h (1, 0) / d;
  s(1, 1) = ((1.0 + h (0, 0) / z1) * (1.0 - z2 * h (1, 1)) + n) / d;
  return s;
}

/*\brief Converts scattering parameters to second hybrid matrix.
  \bug{Programmed formulae are valid only for Z real}
  \bug{Not documented and references}
  \param[in] s Scattering matrix
  \param[in] z1 impedance at input 1
  \param[in] z2 impedance at input 2
  \return second hybrid matrix
  \note Assert 2 by 2 matrix
  \todo Why not s,z1,z2 const
*/
matrix stog (matrix s, nr_complex_t z1, nr_complex_t z2) {
  nr_complex_t n = s (0, 1) * s (1, 0);
  nr_complex_t d = (1.0 + s (0, 0)) * (1.0 - s (1, 1)) + n;
  matrix g (2,2);

  assert (s.rows () >= 2 && s.cols () >= 2);

  g(0, 0) = ((1.0 - s (0, 0)) * (1.0 - s (1, 1)) - n) / z1 / d;
  g(0, 1) = -2.0 * s (0, 1) / d;
  g(1, 0) = +2.0 * s (1, 0) / d;
  g(1, 1)= ((1.0 + s (0, 0)) * (1.0 + s (1, 1)) - n) * z2 / d;
  return g;
}

/*\brief Converts second hybrid matrix to scattering parameters.
  \bug{Programmed formulae are valid only for Z real}
  \bug{Not documented and references}
  \param[in] g second hybrid matrix
  \param[in] z1 impedance at input 1
  \param[in] z2 impedance at input 2
  \return scattering matrix
  \note Assert 2 by 2 matrix
  \todo Why not g,z1,z2 const
*/
matrix gtos (matrix g, nr_complex_t z1, nr_complex_t z2) {
  nr_complex_t n = g (0, 1) * g (1, 0);
  nr_complex_t d = (1.0 + g (0, 0) * z1) * (1.0 + g (1, 1) / z2) - n;
  matrix s (2,2);

  assert (g.rows () >= 2 && g.cols () >= 2);

  s(0, 0)= ((1.0 - g (0, 0) * z1) * (1.0 + g (1, 1) / z2) + n) / d;
  s(0, 1)= -2.0 * g (0, 1) / d;
  s(1, 0)= +2.0 * g (1, 0) / d;
  s(1, 1)= ((g (0, 0) * z1 + 1.0) * (g (1, 1) / z2 - 1.0) - n) / d;
  return s;
}


/*!\brief Generic conversion matrix

  This function converts 2x2 matrices from any of the matrix forms Y,
  Z, H, G and A to any other.  Also converts S<->(A, T, H, Y and Z)
  matrices.
  Convertion assumed:

  Y->Y, Y->Z, Y->H, Y->G, Y->A, Y->S,
  Z->Y, Z->Z, Z->H, Z->G, Z->A, Z->S,
  H->Y, H->Z, H->H, H->G, H->A, H->S,
  G->Y, G->Z, G->H, G->G, G->A, G->S,
  A->Y, A->Z, A->H, A->G, A->A, A->S,
  S->Y, S->Z, S->H, S->G, S->A, S->S,
  S->T,T->T,T->S
  \note assert 2x2 matrix
  \param[in] m base matrix
  \param[in] in matrix
  \param[in] out matrix
  \return matrix given by format out
  \todo m, in, out const
*/
matrix twoport (matrix m, char in, char out) {
  assert (m.rows () >= 2 && m.cols () >= 2);
  nr_complex_t d;
  matrix res (2,2);

  switch (in) {
  case 'Y':
    switch (out) {
    case 'Y': // Y to Y
      res = m;
      break;
    case 'Z': // Y to Z
      d = m (0, 0) * m (1, 1) - m (0, 1) * m (1, 0);
      res(0, 0)= m (1, 1) / d;
      res(0, 1)= -m (0, 1) / d;
      res(1, 0)= -m (1, 0) / d;
      res(1, 1)= m (0, 0) / d;
      break;
    case 'H': // Y to H
      d = m (0, 0);
      res(0, 0) = 1.0 / d;
      res(0, 1)= -m (0, 1) / d;
      res(1, 0)= m (1, 0) / d;
      res(1, 1)= m (1, 1) - m (0, 1) * m (1, 0) / d;
      break;
    case 'G': // Y to G
      d = m (1, 1);
      res(0, 0)= m (0, 0) - m (0, 1) * m (1, 0) / d;
      res(0, 1)= m (0, 1) / d;
      res(1, 0)= -m (1, 0) / d;
      res(1, 1)= 1.0 / d;
      break;
    case 'A': // Y to A
      d = m (1, 0);
      res(0, 0)= -m (1, 1) / d;
      res(0, 1)= -1.0 / d;
      res(1, 0)= m (0, 1) - m (1, 1) * m (0, 0) / d;
      res(1, 1)= -m (0, 0) / d;
      break;
    case 'S': // Y to S
      res = ytos (m);
      break;
    }
    break;
  case 'Z':
    switch (out) {
    case 'Y': // Z to Y
      d = m (0, 0) * m (1, 1) - m (0, 1) * m (1, 0);
      res(0, 0)= m (1, 1) / d;
      res(0, 1)= -m (0, 1) / d;
      res(1, 0)= -m (1, 0) / d;
      res(1, 1)= m (0, 0) / d;
      break;
    case 'Z': // Z to Z
      res = m;
      break;
    case 'H': // Z to H
      d = m (1, 1);
      res(0, 0)= m (0, 0) - m (0, 1) * m (1, 0) / d;
      res(0, 1)= m (0, 1) / d;
      res(1, 0)= -m (1, 0) / d;
      res(1, 1)= 1.0 / d;
      break;
    case 'G': // Z to G
      d = m (0, 0);
      res(0, 0)= 1.0 / d;
      res(0, 1)= -m (0, 1) / d;
      res(1, 0)= m (1, 0) / d;
      res(1, 1)= m (1, 1) - m (0, 1) * m (1, 0) / d;
      break;
    case 'A': // Z to A
      d = m (1, 0);
      res(0, 0)= m (0, 0) / d;
      res(0, 1)= m (0, 0) * m (1, 1) / d - m (0, 1);
      res(1, 0)= 1.0 / d;
      res(1, 1)= m (1, 1) / d;
      break;
    case 'S': // Z to S
      res = ztos (m);
      break;
    }
    break;
  case 'H':
    switch (out) {
    case 'Y': // H to Y
      d = m (0, 0);
      res(0, 0)=1.0 / d;
      res(0, 1)=-m (0, 1) / d;
      res(1, 0)= m (1, 0) / d;
      res(1, 1)= m (1, 1) - m (0, 1) * m(2, 1) / d;
      break;
    case 'Z': // H to Z
      d = m (1, 1);
      res(0, 0)= m (0, 0) - m (0, 1) * m (1, 0) / d;
      res(0, 1)= m (0, 1) / d;
      res(1, 0)= -m (1, 0) / d;
      res(1, 1)= 1.0 / d;
      break;
    case 'H': // H to H
      res = m;
      break;
    case 'G': // H to G
      d = m (0, 0) * m (1, 1) - m (0, 1) * m (1, 0);
      res(0, 0)= m (1, 1) / d;
      res(0, 1)= -m (0, 1) / d;
      res(1, 0)= -m (1, 0) / d;
      res(1, 1)= m (0, 0) / d;
      break;
    case 'A': // H to A
      d = m (1, 0);
      res(0, 0)= m (0, 1) - m (0, 0) * m (1, 1) / d;
      res(0, 1)= -m (0, 0) / d;
      res(1, 0)= -m (1, 1) / d;
      res(1, 1)= -1.0 / d;
      break;
    case 'S': // H to S
      res = htos (m);
      break;
    }
    break;
  case 'G':
    switch (out) {
    case 'Y': // G to Y
      d = m (1, 1);
      res(0, 0)= m (0, 0) - m (0, 1) * m (1, 0) / d;
      res(0, 1)= m (0, 1) / d;
      res(1, 0)= -m (1, 0) / d;
      res(1, 1)= 1.0 / d;
      break;
    case 'Z': // G to Z
      d = m (0, 0);
      res(0, 0)= 1.0 / d;
      res(0, 1)= -m (0, 1) / d;
      res(1, 0)= m (1, 0) / d;
      res(1, 1)= m (1, 1) - m (0, 1) * m (1, 0) / d;
      break;
    case 'H': // G to H
      d = m (0, 0) * m (1, 1) - m (0, 1) * m (1, 0);
      res(0, 0)= m (1, 1) / d;
      res(0, 1)= -m (0, 1) / d;
      res(1, 0)= -m (1, 0) / d;
      res(1, 1)= m (0, 0) / d;
      break;
    case 'G': // G to G
      res = m;
      break;
    case 'A': // G to A
      d = m (1, 0);
      res(0, 0)= 1.0 / d;
      res(0, 1)= m (1, 1) / d;
      res(1, 0)= m (0, 0) / d;
      res(1, 1)= m (0, 0) * m (1, 1) / d - m (0, 1);
      break;
    case 'S': // G to S
      res = gtos (m);
      break;
    }
    break;
  case 'A':
    switch (out) {
    case 'Y': // A to Y
      d = m (0, 1);
      res(0, 0)= m (1, 1) / d;
      res(0, 1)= m (1, 0) - m (0, 0) * m (1, 1) / d;
      res(1, 0)=-1.0 / d;
      res(1, 1)= m (0, 0) / d;
      break;
    case 'Z': // A to Z
      d = m (1, 0);
      res(0, 0)= m (0, 0) / d;
      res(0, 1)= m (0, 0) * m (1, 1) / d - m (0, 1);
      res(1, 0)= 1.0 / d;
      res(1, 1)= m (1, 1) / d;
      break;
    case 'H': // A to H
      d = m (1, 1);
      res(0, 0)= m (0, 1) / d;
      res(0, 1)= m (0, 0) - m (0, 1) * m (1, 0) / d;
      res(1, 0)= -1.0 / d;
      res(1, 1)= m (1, 0) / d;
      break;
    case 'G': // A to G
      d = m (0, 0);
      res(0, 0)=m (1, 0) / d;
      res(0, 1)=m (1, 0) * m (0, 1) / d - m (1, 1);
      res(1, 0)=1.0 / d;
      res(1, 1)=m(0, 1) / d;
      break;
    case 'A': // A to A
      res = m;
      break;
    case 'S': // A to S
      res = atos (m);
      break;
    }
    break;
  case 'S':
    switch (out) {
    case 'S': // S to S
      res = m;
      break;
    case 'T': // S to T
      d = m (1, 0);
      res(0, 0)= m (0, 1) - m (0, 0) * m (1, 1) / d;
      res(0, 1)= m (0, 0) / d;
      res(1, 0)= -m (1, 1) / d;
      res(0, 1)= 1.0 / d;
      break;
    case 'A': // S to A
      res = stoa (m);
      break;
    case 'H': // S to H
      res = stoh (m);
      break;
    case 'G': // S to G
      res = stog (m);
      break;
    case 'Y': // S to Y
      res = stoy (m);
      break;
    case 'Z': // S to Z
      res = stoz (m);
      break;
    }
    break;
  case 'T':
    switch (out) {
    case 'S': // T to S
      d = m (1, 1);
      res(0, 0)= m (0, 1) / d;
      res(0, 1)= m (0, 0) - m (0, 1) * m (1, 0) / d;
      res(1, 0)= 1.0 / d;
      res(0, 1)= -m (1, 0) / d;
      break;
    case 'T': // T to T
      res = m;
      break;
    }
    break;
  }
  return res;
}

/*!\brief Compute the Rollet stabilty factor

   The function returns the Rollet stability factor (\f$K\f) of the given
   S-parameter matrix:
   \[
   K=\frac{1-|S_{11}|^2-|S_{22}|^2+|\delta|^2}{2|S_{12}S_{21}|}
   \]
   Where:
   \[
   \Delta=S_{11}S_{22}-S_{12}S_{21}
   \]
   \param[in] m S parameter matrix
   \return Rollet factor
   \note Assert 2x2 matrix
   \todo m const?
   \todo Rewrite with abs and expand det. It is cleaner.
*/
nr_double_t rollet (matrix m) {
  assert (m.rows () >= 2 && m.cols () >= 2);
  nr_double_t res;
  res = (1 - norm (m (0, 0)) - norm (m (1, 1)) + norm (m.determinant())) /
    2 / abs (m (0, 1) * m (1, 0));
  return res;
}

/* Computes stability measure B1 of the given S-parameter matrix. */
nr_double_t b1 (matrix m) {
  assert (m.rows () >= 2 && m.cols () >= 2);
  nr_double_t res;
  res = 1 + norm (m (0, 0)) - norm (m (1, 1)) - norm (m.determinant());
  return res;
}


matrix rad2deg (matrix a) {
  matrix res (a.rows (), a.cols ());
  for (int r = 0; r < a.rows (); r++)
    for (int c = 0; c < a.cols (); c++)
      res(r, c)=rad2deg (a(r, c));
  return res;
}

matrix deg2rad (matrix a) {
  matrix res (a.rows (), a.cols ());
  for (int r = 0; r < a.rows (); r++)
    for (int c = 0; c < a.cols (); c++)
      res(r, c)= deg2rad (a(r, c));
  return res;
}

} // namespace qucs
