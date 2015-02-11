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
#if USE_EIGEN
typedef Eigen::Matrix<nr_complex_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> matrix;
#endif
 
#ifndef USE_EIGEN
/*!\class matrix
 * \brief Dense complex matrix class
 * This class defines a matrix object with its methods, operators and operations.
*/
class matrix
{
 private:
  Eigen::Matrix<nr_complex_t, Eigen::Dynamic, Eigen::Dynamic> m;
 public:
  matrix () = default;
  matrix (int r, int c) : m(Eigen::Matrix<nr_complex_t,Eigen::Dynamic,Eigen::Dynamic>::Zero(r,c)) {} ;
  matrix (const matrix &)= default;
  matrix (const Eigen::Matrix<nr_complex_t,Eigen::Dynamic,Eigen::Dynamic> &n):
    m(n) {};
  matrix (Eigen::Matrix<nr_complex_t,Eigen::Dynamic,Eigen::Dynamic> &&n):
    m(std::move(n)) {};
  ~matrix () = default;
  auto row(int i) -> decltype(m.row(i)) { return m.row(i); }
  nr_complex_t get (int r, int c) const = delete;
  void set(int r, int c, const nr_complex_t v) = delete;
  int cols (void) const { return m.cols(); }
  int rows (void) const { return m.rows(); }
  nr_complex_t * getData (void) = delete;
  auto data(void) -> decltype(m.data()) { return m.data(); }
  void print (void);
  void exchangeRows (int, int);
  void exchangeCols (int, int);
  matrix transpose() {
    decltype(this->m) temp=m.transpose();
    return temp;
  }
  matrix adjoint() {
    decltype(this->m) temp=m.adjoint();
    return temp;
  }
   matrix conjugate() {
    decltype(this->m) temp=m.conjugate();
    return temp;
  }

  auto array() -> decltype(this->m.array()) {
    return this->m.array();
  }
  
  // operator functions
  friend matrix operator + (const matrix&, const matrix&);
  friend matrix operator - (const matrix&, const matrix&);
  
  friend matrix operator * (const nr_complex_t, const matrix&);
  friend matrix operator * (const matrix&, const nr_complex_t);
  friend matrix operator * (const nr_double_t, const matrix&);
  friend matrix operator * (const matrix&, const nr_double_t);
  friend matrix operator * (const matrix&, const matrix&);

  // intrinsic operator functions
  matrix operator  - () { decltype(this->m) res = -this->m; return res;};

  // other operations
 
  /*! \brief Read access operator
      \param[in] r: row number (from 0 like usually in C)
      \param[in] c: column number (from 0 like usually in C)
      \return Cell in the row r and column c
      \todo: Why not inline
      \todo: Why not r and c not const
      \todo: Create a debug version checking out of bound (using directly assert)
  */
  nr_complex_t  operator () (int r, int c) const { return m(r,c); }
  /*! \brief Write access operator
      \param[in] r: row number (from 0 like usually in C)
      \param[in] c: column number (from 0 like usually in C)
      \return Reference to cell in the row r and column c
      \todo: Why not inline
      \todo: Why r and c not const
      \todo: Create a debug version checking out of bound (using directly assert)
  */
  nr_complex_t& operator () (int r, int c) { return m(r,c); }

};
#endif
inline matrix operator + (const matrix &a, const matrix &b) {
  decltype(a.m) res = a.m;
  res += b.m;
  return res;
}


inline matrix operator - (const matrix &a, const matrix &b) {
  decltype(a.m) res = a.m;
  res -= b.m;
  return res;
}

inline matrix operator * (const matrix&a, const nr_complex_t b) {
  decltype(a.m) res = b*a.m;
  return res;
}

inline matrix operator * (const matrix&a, const nr_double_t b) {
  decltype(a.m) res = b*a.m;
  return res;
}

inline matrix operator * (const nr_complex_t b, const matrix&a) {
  decltype(a.m) res = b*a.m;
  return res;
}

inline matrix operator * (const nr_double_t b,const matrix&a) {
  decltype(a.m) res = b*a.m;
  return res;
}
 
 
inline matrix operator * (const matrix&a, const matrix&b) {
  decltype(a.m) res = a.m;
  res *= b.m;
  return res;
}
 
matrix eye (int);
matrix abs (matrix);
matrix dB (matrix);
matrix arg (matrix);
matrix real (matrix);
matrix imag (matrix);
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
