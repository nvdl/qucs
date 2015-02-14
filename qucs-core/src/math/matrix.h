/*
 * matrix.h - matrix class definitions
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

/*!\file matrix.h
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

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <cassert>

#include <Eigen/Dense>

#include "complex.h"
#include "vector.h"

namespace qucs {

class vector;
typedef Eigen::Matrix<nr_complex_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> matrix;
 
matrix dB (matrix);
matrix arg (matrix);
matrix sqr (matrix);
matrix diagonal (vector);
matrix pow (matrix, int);

/*!\brief S params to S params

  Convert scattering parameters with the reference impedance 'zref'
  to scattering parameters with the reference impedance 'z0'.

  Detail are given in [1], under equation (32)

  New scatering matrix \f$S'\f$ is:
  \f[
  S'=A^{-1}(S-\Gamma^+)(I-\Gamma S)^{-1}A^+
  \f]
  Where x^+ is the adjoint (or complex tranposate) of x,
  I the identity matrix and \f$A\f$ is diagonal the matrix such as:
  \f$   \Gamma_i= r_i \f$ and \f$A\f$ the diagonal matrix such
  as:
  \f[
  A_i =  \frac{(1-r_i^*)\sqrt{|1-r_ir_i^*|}}{|1-r_i|}
  \f]
  Where \f$x*\f$ is the complex conjugate of \f$x\f$
  and \f$r_i\f$ is wave reflexion coefficient of \f$Z_i'\f$ with respect
  to \f$Z_i^*\f$ (where \f$Z_i'\f$ is the new impedance and
  \f$Z_i\f$ is the old impedance), ie:
  \f[
  r_i = \frac{Z_i'-Z_i}{Z_i'-Z_i^*}
  \f]

  \param[in] s original S matrix
  \param[in] zref original reference impedance
  \param[in] z0 new reference impedance
  \bug This formula is valid only for real z!
  \todo Correct documentation about standing waves [1-4]
  \todo Implement Speciale implementation [2-3] if applicable
  \return Renormalized scattering matrix
*/
inline matrix stos (const matrix &s, const qucs::vector &zref, const qucs::vector &z0) {
  int d = s.rows ();
  matrix e, r, a;

  assert (d == s.cols () && d == z0.getSize () && d == zref.getSize ());

  e = e.Identity(d,d);
  r = diagonal ((z0 - zref) / (z0 + zref));
  a = diagonal (sqrt (z0 / zref) / (z0 + zref));
  return a.inverse() * (s - r) * (e - r * s).inverse() * a;
}

/*!\brief S renormalization with all part identic
   \param[in] s original S matrix
   \param[in] zref original reference impedance
   \param[in] z0 new reference impedance
   \return Renormalized scattering matrix
*/
inline matrix stos (const matrix &s, const nr_complex_t zref, const nr_complex_t z0 = 50.0) {
  int d = s.rows ();
  return stos (s, qucs::vector (d, zref), qucs::vector (d, z0));
}

/*!\brief S renormalization with all part identic and real
  \param[in] s original S matrix
  \param[in] zref original reference impedance
  \param[in] z0 new reference impedance
  \return Renormalized scattering matrix
*/
inline matrix stos (const matrix & s, const nr_double_t zref, const nr_double_t z0 = 50.0) {
  return stos (s, nr_complex_t (zref, 0), nr_complex_t (z0, 0));
}

/*!\brief S renormalization (variation)
   \param[in] s original S matrix
   \param[in] zref original reference impedance
   \param[in] z0 new reference impedance
*/
inline matrix stos (const matrix &s, const qucs::vector &zref, const nr_complex_t z0 = 50.0) {
  return stos (s, zref, qucs::vector (zref.getSize (), z0));
}

/*!\brief S renormalization (variation)
  \param[in] s original S matrix
  \param[in] zref original reference impedance
  \param[in] z0 new reference impedance
  \return Renormalized scattering matrix
*/
inline matrix stos (const matrix &s, const nr_complex_t zref, const qucs::vector &z0) {
  return stos (s, qucs::vector (z0.getSize (), zref), z0);
}

/*!\brief Scattering parameters to impedance matrix

  Convert scattering parameters to impedance matrix.
  According to [1] eq (19):
  \f[
  Z=F^{-1} (I- S)^{-1} (SG + G^+) F
  \f]
  Where \f$S\f$ is the scattering matrix, \f$x^+\f$
  is the adjoint of x, I the identity matrix. The matrix
  F and G are diagonal matrix defined by:
  \f{align*}
  F_i&=\frac{1}{2\sqrt{\Re\text{e}\; Z_i}} \\
  G_i&=Z_i
  \f}
  \param[in] s Scattering matrix
  \param[in] z0 Normalisation impedance
  \note We could safely drop the \f$1/2\f$ in \f$F\f$ because we compute
        \f$FXF^{-1}\f$ and therefore \f$1/2\f$ will simplify.
  \bug not correct if zref is complex
  \return Impedance matrix
*/
inline matrix stoz (const matrix &s, const qucs::vector &z0) {
  int d = s.rows ();
  matrix e, zref, gref;

  assert (d == s.cols () && d == z0.getSize ());

  e = e.Identity(d,d);
  zref = diagonal (z0);
  gref = diagonal (sqrt (real (1 / z0)));
  return gref.inverse() * (e - s).inverse() * (s * zref + zref) * gref;
}

/*!\brief Scattering parameters to impedance matrix identic case
   \param[in] s Scattering matrix
   \param[in] z0 Normalisation impedance
   \return Impedance matrix
*/
inline matrix stoz (const matrix &s, const nr_complex_t z0 = 50.0) {
  return stoz (s, qucs::vector (s.rows (), z0));
}
/*!\brief Convert impedance matrix scattering parameters.

   Convert scattering parameters to impedance matrix.
   According to [1] eq (18):
  \f[
  S=F(Z-G^+)(Z+G)^{-1} F^{-1}
  \f]
  Where \f$Z\f$ is the scattering matrix, \f$x^+\f$
  is the adjoint of x, I the identity matrix. The matrix
  F and G are diagonal matrix defined by:
  \f{align*}
  F_i&=\frac{1}{2\sqrt{\Re\text{e}\; Z_i}} \\
  G_i&=Z_i
  \f}
  \param[in] Z Impedance matrix
  \param[in] z0 Normalisation impedance
  \return Scattering matrix
  \note We could safely drop the \f$1/2\f$ in \f$F\f$ because we compute
        \f$FXF^{-1}\f$ and therefore \f$1/2\f$ will simplify.
  \bug not correct if zref is complex
*/
inline matrix ztos (const matrix &z, const qucs::vector &z0) {
  int d = z.rows ();
  matrix e, zref, gref;

  assert (d == z.cols () && d == z0.getSize ());

  e = e.Identity(d,d);
  zref = diagonal (z0);
  gref = diagonal (sqrt (real (1 / z0)));
  return gref * (z - zref) * (z + zref).inverse() * gref.inverse();
}

/*!\brief Convert impedance matrix to scattering parameters identic case
   \param[in] Z Impedance matrix
   \param[in] z0 Normalisation impedance
   \return Scattering matrix
 */
inline matrix ztos (const matrix &z, const nr_complex_t z0 = 50.0) {
  return ztos (z, qucs::vector (z.rows (), z0));
}

/*!\brief impedance matrix to admittance matrix.

   Convert impedance matrix to admittance matrix. By definition
   \f$Y=Z^{-1}\f$
   \param[in] z impedance matrix
   \return Admittance matrix
*/
inline matrix ztoy (const matrix& z) {
  assert (z.rows () == z.cols ());
  return z.inverse();
}
/*!\brief Scattering parameters to admittance matrix.

  Convert scattering parameters to admittance matrix.
  According to [1] eq (19):
  \f[
  Z=F^{-1} (I- S)^{-1} (SG + G^+) F
  \f]
  Where \f$S\f$ is the scattering matrix, \f$x^+\f$
  is the adjoint of x, I the identity matrix. The matrix
  F and G are diagonal matrix defined by:
  \f{align*}
  F_i&=\frac{1}{2\sqrt{\Re\text{e}\; Z_i}} \\
  G_i&=Z_i
  \f}
  Using the well know formula \f$(AB)^{-1}=B^{1}A^{1}\f$,
  we derivate:
  \f[
  Y=F^{-1} (SG+G^+)^{-1} (I-S) F
  \f]
  \param[in] s Scattering matrix
  \param[in] z0 Normalisation impedance
  \note We could safely drop the \f$1/2\f$ in \f$F\f$ because we compute
        \f$FXF^{-1}\f$ and therefore \f$1/2\f$ will simplify.
  \bug not correct if zref is complex
  \return Admittance matrix
*/
inline matrix stoy (const matrix &s, const qucs::vector &z0) {
  int d = s.rows ();
  matrix e, zref, gref;

  assert (d == s.cols () && d == z0.getSize ());

  e = e.Identity(d,d);
  zref = diagonal (z0);
  gref = diagonal (sqrt (real (1 / z0)));
  return gref.inverse() * (s * zref + zref).inverse() * (e - s) * gref;
}
 
/*!\brief Convert scattering pto adminttance parameters identic case
   \param[in] S Scattering matrix
   \param[in] z0 Normalisation impedance
   \return Admittance matrix
*/
inline matrix stoy (const matrix &s, const nr_complex_t z0 = 50.0) {
  return stoy (s, qucs::vector (s.rows (), z0));
}

/*!\brief Admittance matrix to scattering parameters

   Convert admittance matrix to scattering parameters.
   Using the same methodology as [1] eq (16-19), but writing
   (16) as \f$i=Yv\f$, ie
   \f[
   S=F(I-G^+Y)(I-GY)^{-1}F^{-1}
   \f]
   Where \f$S\f$ is the scattering matrix, \f$x^+\f$
   is the adjoint of x, I the identity matrix. The matrix
   F and G are diagonal matrix defined by:
   \f{align*}
   F_i&=\frac{1}{2\sqrt{\Re\text{e}\; Z_i}} \\
   G_i&=Z_i
   \f}
   Using the well know formula \f$(AB)^{-1}=B^{1}A^{1}\f$,
   we derivate:
   \f[
   Y=F^{-1} (SG+G^+)^{-1} (I-S) F
   \f]
   \param[in] y admittance matrix
   \param[in] z0 Normalisation impedance
   \note We could safely drop the \f$1/2\f$ in \f$F\f$ because we compute
         \f$FXF^{-1}\f$ and therefore \f$1/2\f$ will simplify.
   \bug not correct if zref is complex
   \return Scattering matrix
*/
inline matrix ytos (const matrix &y, const qucs::vector &z0) {
  int d = y.rows ();
  matrix e, zref, gref;

  assert (d == y.cols () && d == z0.getSize ());
  
  e = e.Identity(d,d);
  zref = diagonal (z0);
  gref = diagonal (sqrt (real (1 / z0)));
  return gref * (e - zref * y) * (e + zref * y).inverse() * gref.inverse();
}

/*!\brief Convert Admittance matrix to scattering parameters identic case
   \param[in] y Admittance matrix
   \param[in] z0 Normalisation impedance
   \return Scattering matrix
 */
inline matrix ytos (const matrix &y, const nr_complex_t z0 = 50.0) {
  return ytos (y, qucs::vector (y.rows (), z0));
}

/*!\brief Convert admittance matrix to impedance matrix.

  Convert \f$Y\f$ matrix to \f$Z\f$ matrix using well known relation
  \f$Z=Y^{-1}\f$
  \param[in] y admittance matrix
  \return Impedance matrix
  \note Check if y matrix is a square matrix
*/
inline matrix ytoz (const matrix &y) {
  assert (y.rows () == y.cols ());
  return y.inverse();
}

/*!\brief Admittance noise correlation matrix to S-parameter noise
   correlation matrix

   Converts admittance noise correlation matrix to S-parameter noise
   correlation matrix. According to [7] fig 2:
   \f[
   C_s=\frac{1}{4}(I+S)C_y(I+S)^+
   \f]
   Where \f$C_s\f$ is the scattering noise correlation matrix,
   \f$C_y\f$ the admittance noise correlation matrix, \f$I\f$
    the identity matrix and \f$S\f$ the scattering matrix
    of device. \f$x^+\f$ is the adjoint of \f$x\f$
   \warning cy matrix and s matrix are assumed to be normalized
   \param[in] cy Admittance noise correlation
   \param[in] s S parameter matrix of device
   \return S-parameter noise correlation matrix
   \note Assert compatiblity of matrix
*/
inline matrix cytocs (const matrix &cy, const matrix &s) {
  matrix e = e.Identity(s.rows (),s.rows());

  assert (cy.rows () == cy.cols () && s.rows () == s.cols () &&
	  cy.rows () == s.rows ());

  return 0.25 * (e + s) * cy * (e + s).adjoint();
}

/*!\brief Converts S-parameter noise correlation matrix to admittance noise
    correlation matrix.

    According to [7] fig 2:
    \f[
    C_y=(I+Y)C_s(I+Y)^+
    \f]
    Where \f$C_s\f$ is the scattering noise correlation matrix,
    \f$C_y\f$ the admittance noise correlation matrix, \f$I\f$
    the identity matrix and \f$S\f$ the scattering matrix
    of device. \f$x^+\f$ is the adjoint of \f$x\f$
    \warning cs matrix and y matrix are assumed to be normalized
    \param[in]  cs S parameter noise correlation
    \param[in] y Admittance matrix of device
    \return admittance noise correlation matrix
*/
inline matrix cstocy (const matrix &cs, const matrix &y) {
  matrix e = e.Identity(y.rows (),y.rows());

  assert (cs.rows () == cs.cols () && y.rows () == y.cols () &&
	  cs.rows () == y.rows ());

  return (e + y) * cs * (e + y).adjoint();
}

/*!\brief Converts impedance noise correlation matrix to S-parameter noise
   correlation matrix.

   According to [7] fig 2:
   \f[
   C_s=\frac{1}{4}(I-S)C_z(I-S)
   \f]
   Where \f$C_s\f$ is the scattering noise correlation matrix,
   \f$C_z\f$ the impedance noise correlation matrix, \f$I\f$
    the identity matrix and \f$S\f$ the scattering matrix
    of device. \f$x^+\f$ is the adjoint of \f$x\f$
   \warning Cz matrix and s matrix are assumed to be normalized
   \param[in] cz Impedance noise correlation
   \param[in] s S parameter matrix of device
   \return S-parameter noise correlation matrix
   \note Assert compatiblity of matrix
*/
inline matrix cztocs (const matrix &cz, const matrix &s) {
  matrix e = e.Identity (s.rows (),s.rows());

  assert (cz.rows () == cz.cols () && s.rows () == s.cols () &&
	  cz.rows () == s.rows ());

  return 0.25 * (e - s) * cz * (e - s).adjoint();
}

/*!\brief Converts S-parameter noise correlation matrix to impedance noise
    correlation matrix.

    According to [7] fig 2:
    \f[
    C_z=(I+Z)C_s(I+Z)^+
    \f]
    Where \f$C_s\f$ is the scattering noise correlation matrix,
    \f$C_z\f$ the impedance noise correlation matrix, \f$I\f$
    the identity matrix and \f$S\f$ the scattering matrix
    of device. \f$x^+\f$ is the adjoint of \f$x\f$
    \warning cs matrix and y matrix are assumed to be normalized
    \param[in]  cs S parameter noise correlation
    \param[in] z Impedance matrix of device
    \return Impedance noise correlation matrix
*/
inline matrix cstocz (const matrix &cs, const matrix &z) {
  assert (cs.rows () == cs.cols () && z.rows () == z.cols () &&
	  cs.rows () == z.rows ());
  matrix e = e.Identity (z.rows (),z.rows());
  return (e + z) * cs * (e + z).adjoint();
}

/*!\brief Converts impedance noise correlation matrix to admittance noise
    correlation matrix.

    According to [7] fig 2:
    \f[
    C_y=YC_zY^+
    \f]
    Where \f$C_z\f$ is the impedance correlation matrix,
    \f$I\f$ the identity matrix and \f$C_y\f$ the admittance noise
    correlation matrix.
    \f$x^+\f$ is the adjoint of \f$x\f$
    \warning cz matrix and y matrix are assumed to be normalized
    \param[in]  cz impedance noise correlation
    \param[in]  y Admittance matrix of device
    \return admittance noise correlation matrix
*/
inline matrix cztocy (const matrix &cz, const matrix &y) {
  assert (cz.rows () == cz.cols () && y.rows () == y.cols () &&
	  cz.rows () == y.rows ());

  return y * cz * y.adjoint();
}

/*!\brief Converts admittance noise correlation matrix to impedance noise
    correlation matrix.

    According to [7] fig 2:
    \f[
    C_z=ZC_yZ^+
    \f]
    Where \f$C_z\f$ is the impedance correlation matrix,
    \f$I\f$ the identity matrix and \f$C_y\f$ the admittance noise
    correlation matrix.
    \f$x^+\f$ is the adjoint of \f$x\f$
    \warning cy matrix and z matrix are assumed to be normalized
    \param[in]  cy Admittance noise correlation
    \param[in]  z Impedance matrix of device
    \return Impedance noise correlation matrix
    \todo cs, z const
*/
inline matrix cytocz (const matrix& cy, const matrix &z) {
  assert (cy.rows () == cy.cols () && z.rows () == z.cols () &&
	  cy.rows () == z.rows ());
  return z * cy * z.adjoint();
}

 
matrix stoa (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix atos (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix stoh (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix htos (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix stog (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix gtos (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix twoport (matrix, char, char);
nr_double_t rollet (matrix);
nr_double_t b1 (matrix);
matrix rad2deg     (matrix);
matrix deg2rad     (matrix);
} // namespace qucs

#endif /* __MATRIX_H__ */
