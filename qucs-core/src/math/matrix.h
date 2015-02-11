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
   \brief Dense matrix class header file
*/

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <Eigen/Dense>

namespace qucs {

class vector;
typedef Eigen::Matrix<nr_complex_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> matrix;
 
matrix eye (int);
matrix dB (matrix);
matrix arg (matrix);
matrix sqr (matrix);
matrix eye (int, int);
matrix diagonal (vector);
matrix pow (matrix, int);
nr_complex_t cofactor (matrix, int, int);
nr_complex_t detLaplace (matrix);
nr_complex_t detGauss (matrix);
nr_complex_t det (matrix);
matrix inverseLaplace (matrix);
matrix inverseGaussJordan (matrix);
matrix inverse (matrix);
matrix stos (matrix, nr_complex_t, nr_complex_t z0 = 50.0);
matrix stos (matrix, nr_double_t, nr_double_t z0 = 50.0);
matrix stos (matrix, vector, nr_complex_t z0 = 50.0);
matrix stos (matrix, nr_complex_t, vector);
matrix stos (matrix, vector, vector);
matrix stoz (matrix, nr_complex_t z0 = 50.0);
matrix stoz (matrix, vector);
matrix ztos (matrix, nr_complex_t z0 = 50.0);
matrix ztos (matrix, vector);
matrix ztoy (matrix);
matrix stoy (matrix, nr_complex_t z0 = 50.0);
matrix stoy (matrix, vector);
matrix ytos (matrix, nr_complex_t z0 = 50.0);
matrix ytos (matrix, vector);
matrix ytoz (matrix);
matrix stoa (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix atos (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix stoh (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix htos (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix stog (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix gtos (matrix, nr_complex_t z1 = 50.0, nr_complex_t z2 = 50.0);
matrix cytocs (matrix, matrix);
matrix cztocs (matrix, matrix);
matrix cztocy (matrix, matrix);
matrix cstocy (matrix, matrix);
matrix cytocz (matrix, matrix);
matrix cstocz (matrix, matrix);
matrix twoport (matrix, char, char);
nr_double_t rollet (matrix);
nr_double_t b1 (matrix);
matrix rad2deg     (matrix);
matrix deg2rad     (matrix);
} // namespace qucs

#endif /* __MATRIX_H__ */
