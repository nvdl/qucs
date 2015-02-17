/*
 * matvec.cpp - matrix vector class implementation
 *
 * Copyright (C) 2004-2009 Stefan Jahn <stefan@lkcc.org>
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

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

#include "logging.h"
#include "object.h"
#include "complex.h"
#include "vector.h"
#include "matrix.h"
#include "matvec.h"

#if !HAVE_STRCHR
# define strchr  index
# define strrchr rindex
#endif

namespace qucs {

// Sets the name of the matvec object.
void matvec::setName (const std::string &n) {
  name = n;
}

// Returns the name of the matvec object.
std::string matvec::getName (void) const {
  return name;
}

/* This function saves the given vector to the matvec object with the
   appropriate matrix indices. */
void matvec::set (qucs::vector v, int r, int c) {
  assert (v.getSize () == data.size() &&
	  r >= 0 && r < rows && c >= 0 && c < cols);
  for (size_type i = 0; i < data.size(); i++)
    (data[i])(r, c)= v.get (i);
}

/* The function returns the vector specified by the given matrix
   indices.  If the matrix vector has a valid name 'A' the returned
   vector gets the name 'A[r,c]'. */
qucs::vector matvec::get (int r, int c) {
  assert (r >= 0 && r < rows && c >= 0 && c < cols);
  qucs::vector res;
  for (size_type i = 0; i < data.size(); i++) res.add ((data[i])(r, c));
  if (!name.empty()) {
    res.setName (createMatrixString (name.c_str(), r, c));
  }
  return res;
}

/* This function returns a static text representation with the
   'n[r,c]' scheme indicating a matrix (vector) entry. */
char * matvec::createMatrixString (const char * n, int r, int c) {
  static char str[256]; // hopefully enough
  sprintf (str, "%s[%d,%d]", n, r + 1, c + 1);
  return str;
}

/* This function also returns a static text representation with the
   'n[r,c]' scheme indicating a matrix (vector) entry but with
   different arguments. */
char * matvec::createMatrixString (char n, int r, int c) {
  static char str[256]; // hopefully enough
  sprintf (str, "%c[%d,%d]", n, r + 1, c + 1);
  return str;
}

/* The function investigates the given vectors name.  If this name
   matches the 'n[r,c]' pattern it returns the name 'n' and saves the
   row and column indices as well.  The caller is responsible to
   'free()' the returned string.  If the vectors name does not match
   the pattern the function returns NULL. */
char * matvec::isMatrixVector (const char * n, int& r, int& c) {
  const char * p; int len;
  char *pnew;
  if (n == NULL) return NULL;              // nothing todo here
  if ((p = strchr (n, '[')) != NULL) {     // find first '['
    r = atoi (p + 1) - 1;                  // get first index
    if ((p = strchr (p, ',')) != NULL) {   // find the ','
      c = atoi (p + 1) - 1;                // get second index
      if ((p = strchr (p, ']')) != NULL) { // find trailing ']'
	if (p[1] == '\0') {                // identifier must end in ']'
	  // parse actual identifier
	  if ((len = strchr (n, '[') - n) > 0) {
	    pnew = (char *) malloc (len + 1);
	    memcpy (pnew, n, len);
	    pnew[len] = '\0';
	    return pnew;
	  }
	}
      }
    }
  }
  return NULL;
}

/* This function looks through the vector list given in `data' to find
   matrix entries specified by `name' and returns the matrix vector
   dimensions. */
void matvec::getMatrixVectorSize (qucs::vector * data, char * name,
				  int& rs, int& cs, int& ss) {
  qucs::vector * v;
  char * n;
  const char *vn;
  int r, c, s;
  rs = cs = ss = -1;
  // go through vector list
  for (v = data; v != NULL; v = (qucs::vector *) v->getNext ()) {
    vn = v->getName ();
    // requested matrix name found?
    if (strstr (vn, name) == vn) {
      if ((n = matvec::isMatrixVector (vn, r, c)) != NULL) {
        if (rs < r) rs = r;
        if (cs < c) cs = c;
        s = v->getSize ();
	if (ss < s) ss = s;
        free (n);
      }
    }
  }
}

/* This function looks through the vector list given in `data' to find
   matrix entries specified by `name' and returns a matrix vector
   object.  If there are no such matrices the function returns
   NULL. */
matvec * matvec::getMatrixVector (qucs::vector * data, char * name) {

  // obtain matrix vector dimensions
  int rs, cs, ss;
  getMatrixVectorSize (data, name, rs, cs, ss);

  qucs::vector * v;
  const char * vn;
  char * n;
  int r, c;
  // valid matrix entries found
  if (rs >= 0 && cs >= 0 && ss > 0) {
    // create matrix vector
    matvec * mv = new matvec (ss, rs + 1, cs + 1);
    mv->setName (name);
    // go through vector list again and fill in matrix vectors
    for (v = data; v; v = (qucs::vector *) v->getNext ()) {
      vn = v->getName ();
      if (strstr (vn, name) == vn) {
        if ((n = matvec::isMatrixVector (vn, r, c)) != NULL) {
          mv->set (*v, r, c);
          free (n);
        }
      }
    }
    return mv;
  }
  return NULL;
}

/* This function saves the given matrix in the matrix vector at the
   specified position. */
void matvec::set (matrix m, int idx) {
  assert (m.rows () == rows && m.cols () == cols &&
	  idx >= 0 && idx < data.size());
  data[idx] = m;
}

/* The function returns the matrix stored within the matrix vector at
   the given position. */
matrix matvec::get (size_type idx) const {
  assert (idx < data.size());
  return data[idx];
}

// Matrix vector addition.
matvec operator + (const matvec &a, const matvec &b) {
  assert (a.getRows () == b.getRows () && a.getCols () == b.getCols () &&
	  a.size () == b.size ());
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++)
    res.set (a.get (i) + b.get (i), i);
  return res;
}

// Matrix vector addition with single matrix.
matvec operator + (const matvec &a, const matrix &b) {
  assert (a.getRows () == b.rows () && a.getCols () == b.cols ());
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (a.get (i) + b, i);
  return res;
}

// Matrix vector addition with vector.
matvec operator + (const matvec &a, const qucs::vector &b) {
  assert (a.size () == b.getSize ());
  matvec res (a);
  for (matvec::size_type i = 0; i < a.size (); i++) {
    res[i] = (res[i].array()+b.get(i)).matrix();
  }
  return res;
}

// Matrix vector addition with vector in different order.
matvec operator + (const qucs::vector &b, const matvec &a) {
  return a + b;
}

// Matrix vector addition with single matrix in different order.
matvec operator + (const matrix &a, const matvec &b) {
  return b + a;
}

// Matrix vector scalar addition.
matvec operator + (const matvec &a, nr_complex_t z) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) {
    matrix m((a.get(i).array()+z).matrix());
    res.set (m, i);
  }
  return res;
}

// Matrix vector scalar addition in different order.
matvec operator + (nr_complex_t z, const matvec &a) {
  return a+z;
}

// Matrix vector scalar addition.
matvec operator + (const matvec &a, nr_double_t d) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) {
    matrix m((a.get(i).array()+d).matrix());
    res.set (m, i);
  }
  return res;
}

// Matrix vector scalar addition in different order.
matvec operator + (nr_double_t d, const matvec &a) {
  return a+d;
}

// Matrix vector scalar subtraction.
matvec operator - (const matvec &a, nr_complex_t z) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) {
    matrix m((a.get(i).array()-z).matrix());
    res.set (m, i);
  }
  return res;
}

// Matrix vector scalar subtraction in different order.
matvec operator - (nr_complex_t z, const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) {
    matrix m((z-a.get(i).array()).matrix());
    res.set (m, i);
  }
  return res;
}

// Matrix vector scalar subtraction.
matvec operator - (const matvec &a, nr_double_t d) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) {
    matrix m((a.get(i).array()-d).matrix());
    res.set (m, i);
  }
  return res;
}

// Matrix vector scalar subtraction in different order.
matvec operator - (nr_double_t d, const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) {
    matrix m((d-a.get(i).array()).matrix());
    res.set (m, i);
  }
  return res;
}

// Intrinsic matrix vector addition.
matvec matvec::operator += (const matvec &a) {
  assert (a.getRows () == rows && a.getCols () == cols &&
	  a.size () == data.size());
  for (int i = 0; i < data.size(); i++) data[i] = data[i] + a.get (i);
  return *this;
}

// Matrix vector subtraction.
matvec operator - (const matvec &a, const matvec &b) {
  assert (a.getRows () == b.getRows () && a.getCols () == b.getCols () &&
	  a.size () == b.size ());
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (a.get (i) - b.get (i), i);
  return res;
}

// Matrix vector subtraction with single matrix.
matvec operator - (const matvec& a, const matrix& b) {
  assert (a.getRows () == b.rows () && a.getCols () == b.cols ());
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (a.get (i) - b, i);
  return res;
}

// Matrix vector subtraction with single matrix in different order.
matvec operator - (const matrix& a, const matvec &b) {
  matvec res(b);
  return -res + a;
}

// Matrix vector subtraction with vector.
matvec operator - (const matvec& a, const qucs::vector &b) {
  qucs::vector res(b);
  return -res + a;
}

// Matrix vector subtraction with vector in different order.
matvec operator - (const qucs::vector &b, const matvec &a) {
  matvec res(a);
  return -res + b;
}

// Unary minus.
matvec matvec::operator - () {
  matvec res (size (), getRows (), getCols ());
  for (int i = 0; i < size (); i++) res.set (-data[i], i);
  return res;
}

// Intrinsic matrix vector subtraction.
matvec matvec::operator -= (const matvec &a) {
  assert (a.getRows () == rows && a.getCols () == cols &&
	  a.size () == data.size());
  for (int i = 0; i < a.size (); i++) data[i] = data[i] - a.get (i);
  return *this;
}

// Matrix vector scaling.
matvec operator * (const matvec &a, nr_complex_t z) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (a.get (i) * z, i);
  return res;
}

// Matrix vector scaling in different order.
matvec operator * (nr_complex_t z, const matvec &a) {
  return a * z;
}

// Scalar matrix vector scaling.
matvec operator * (const matvec &a, nr_double_t d) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (a.get (i) * d, i);
  return res;
}

// Scalar matrix vector scaling in different order.
matvec operator * (nr_double_t d,const matvec &a) {
  return a * d;
}

// Matrix vector scaling by a second vector.
matvec operator * (const matvec &a, const qucs::vector &b) {
  assert (a.size () == b.getSize ());
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (a.get (i) * b.get (i), i);
  return res;
}

// Matrix vector scaling by a second vector in different order.
matvec operator * (const qucs::vector &a, const matvec &b) {
  return b * a;
}

// Matrix vector scaling.
matvec operator / (const matvec &a, const nr_complex_t z) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set ((1.0/z)*a.get (i) , i);
  return res;
}

// Scalar matrix vector scaling.
matvec operator / (const matvec &a, nr_double_t d) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set ((1.0/d)*a.get (i), i);
  return res;
}

// Matrix vector scaling by a second vector.
matvec operator / (const matvec& a, const qucs::vector &b) {
  assert (a.size () == b.getSize ());
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (a.get (i)* (1.0/b.get (i)), i);
  return res;
}

// Matrix vector multiplication.
matvec operator * (const matvec &a, const matvec &b) {
  assert (a.getCols () == b.getRows () && a.size () == b.size ());
  matvec res (a.size (), a.getRows (), b.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (a.get (i) * b.get (i), i);
  return res;
}

// Matrix vector multiplication with a single matrix.
matvec operator * (const matvec &a, const matrix &b) {
  assert (a.getCols () == b.rows ());
  matvec res (a.size (), a.getRows (), b.cols ());
  for (int i = 0; i < a.size (); i++) res.set (a.get (i) * b, i);
  return res;
}

// Matrix vector multiplication with a single matrix in different order.
matvec operator * (const matrix &a,const matvec &b) {
  return b * a;
}

// Compute determinants of the given matrix vector.
qucs::vector det (const matvec &a) {
  qucs::vector res (a.size ());
  for (int i = 0; i < a.size (); i++) res.set (a.get(i).determinant(), i);
  return res;
}

// Compute inverse matrices of the given matrix vector.
matvec inverse (const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) {
    res.set (a.get(i).inverse(), i);
  }
  return res;
}

// Compute inverse matrices of the given matrix vector.
matvec sqr (const matvec &a) {
  return a * a;
}

// Compute n-th power of the given matrix vector.
matvec pow (const matvec &a, int n) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (pow (a.get (i), n), i);
  return res;
}

// Compute n-th powers in the vector of the given matrix vector.
matvec pow (const matvec &a, const qucs::vector &v) {
  assert (a.size () == v.getSize ());
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++)
    res.set (pow (a.get (i), (int) real (v.get (i))), i);
  return res;
}

// Conjugate complex matrix vector.
matvec conj (const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set ((a.get (i)).conjugate(), i);
  return res;
}

// Computes magnitude of each matrix vector element.
matvec abs (const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) {
    matrix temp(a.get(i).cwiseAbs().cast<nr_complex_t>());
    res.set (temp, i);
  }
  return res;
}

// Computes magnitude in dB of each matrix vector element.
matvec dB (const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (dB (a.get (i)), i);
  return res;
}

// Computes the argument of each matrix vector element.
matvec arg (const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (arg (a.get (i)), i);
  return res;
}

// Real part matrix vector.
matvec real (const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++)
    {
      matrix m = a.get(i).real().cast<nr_complex_t>();
      res.set (m, i);
    }
  return res;
}

// Real part matrix vector.
matvec imag (const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++)
     {
       matrix m = a.get(i).imag().cast<nr_complex_t>();
       res.set (m, i);
    }
  return res;
}

/* The function returns the adjoint complex matrix vector.  This is
   also called the adjugate or transpose conjugate. */
matvec adjoint (const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set ((a.get (i)).adjoint(), i);
  return res;
}

// Transpose the matrix vector.
matvec transpose (const matvec &a) {
  matvec res (a.size (), a.getCols (), a.getRows ());
  for (int i = 0; i < a.size (); i++)
    res.set (a.get (i).transpose(), i);
  return res;
}

/* Convert scattering parameters with the reference impedance 'zref'
   to scattering parameters with the reference impedance 'z0'. */
matvec stos (const matvec &s, const qucs::vector &zref, const qucs::vector &z0) {
  assert (s.getCols () == s.getRows () &&
	  s.getCols () == zref.getSize () && s.getCols () == z0.getSize ());
  matvec res (s.size (), s.getCols (), s.getRows ());
  for (int i = 0; i < s.size (); i++)
    res.set (stos (s.get (i), zref, z0), i);
  return res;
}

matvec stos (const matvec &s, nr_complex_t zref, nr_complex_t z0) {
  int d = s.getRows ();
  return stos (s, qucs::vector (d, zref), qucs::vector (d, z0));
}

matvec stos (const matvec &s, nr_double_t zref, nr_double_t z0) {
  return stos (s, nr_complex_t (zref, 0), nr_complex_t (z0, 0));
}

matvec stos (const matvec &s, const qucs::vector &zref, nr_complex_t z0) {
  return stos (s, zref, qucs::vector (zref.getSize (), z0));
}

matvec stos (const matvec &s, nr_complex_t zref, const qucs::vector &z0) {
  return stos (s, qucs::vector (z0.getSize (), zref), z0);
}

// Convert scattering parameters to admittance matrix vector.
matvec stoy (const matvec &s, const qucs::vector &z0) {
  assert (s.getCols () == s.getRows () && s.getCols () == z0.getSize ());
  matvec res (s.size (), s.getCols (), s.getRows ());
  for (int i = 0; i < s.size (); i++) res.set (stoy (s.get (i), z0), i);
  return res;
}

matvec stoy (const matvec &s, nr_complex_t z0) {
  return stoy (s, qucs::vector (s.getCols (), z0));
}

// Convert admittance matrix to scattering parameter matrix vector.
matvec ytos (const matvec &y, const qucs::vector &z0) {
  assert (y.getCols () == y.getRows () && y.getCols () == z0.getSize ());
  matvec res (y.size (), y.getCols (), y.getRows ());
  for (int i = 0; i < y.size (); i++) res.set (ytos (y.get (i), z0), i);
  return res;
}

matvec ytos (const matvec &y, nr_complex_t z0) {
  return ytos (y, qucs::vector (y.getCols (), z0));
}

// Convert scattering parameters to impedance matrix vector.
matvec stoz (const matvec &s, const qucs::vector &z0) {
  assert (s.getCols () == s.getRows () && s.getCols () == z0.getSize ());
  matvec res (s.size (), s.getCols (), s.getRows ());
  for (int i = 0; i < s.size (); i++) res.set (stoz (s.get (i), z0), i);
  return res;
}

matvec stoz (const matvec &s, nr_complex_t z0) {
  return stoz (s, qucs::vector (s.getCols (), z0));
}

// Convert impedance matrix vector scattering parameter matrix vector.
matvec ztos (const matvec &z, const qucs::vector &z0) {
  assert (z.getCols () == z.getRows () && z.getCols () == z0.getSize ());
  matvec res (z.size (), z.getCols (), z.getRows ());
  for (int i = 0; i < z.size (); i++) res.set (ztos (z.get (i), z0), i);
  return res;
}

matvec ztos (const matvec &z, nr_complex_t z0) {
  return ztos (z, qucs::vector (z.getCols (), z0));
}

// Convert impedance matrix vector to admittance matrix vector.
matvec ztoy (const matvec &z) {
  assert (z.getCols () == z.getRows ());
  matvec res (z.size (), z.getCols (), z.getRows ());
  for (int i = 0; i < z.size (); i++) res.set (ztoy (z.get (i)), i);
  return res;
}

// Convert admittance matrix vector to impedance matrix vector.
matvec ytoz (const matvec &y) {
  assert (y.getCols () == y.getRows ());
  matvec res (y.size (), y.getCols (), y.getRows ());
  for (int i = 0; i < y.size (); i++) res.set (ytoz (y.get (i)), i);
  return res;
}

/* This function converts 2x2 matrix vectors from any of the matrix
   forms Y, Z, H, G and A to any other.  Also converts S<->(A, T, H, Y
   and Z) matrix vectors. */
matvec twoport (const matvec &m, char in, char out) {
  assert (m.getCols () >= 2 && m.getRows () >= 2);
  matvec res (m.size (), 2, 2);
  for (int i = 0; i < m.size (); i++)
    res.set (twoport (m.get (i), in, out), i);
  return res;
}

/* The function returns the Rollet stability factor vector of the
   given S-parameter matrix vector. */
qucs::vector rollet (const matvec &m) {
  assert (m.getCols () >= 2 && m.getRows () >= 2);
  qucs::vector res (m.size ());
  for (int i = 0; i < m.size (); i++) res.set (rollet (m.get (i)), i);
  return res;
}

/* The function returns the stability measure B1 vector of the given
   S-parameter matrix vector. */
qucs::vector b1 (const matvec &m) {
  assert (m.getCols () >= 2 && m.getRows () >= 2);
  qucs::vector res (m.size ());
  for (int i = 0; i < m.size (); i++) res.set (b1 (m.get (i)), i);
  return res;
}

matvec rad2deg (const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (rad2deg (a.get (i)), i);
  return res;
}

matvec deg2rad (const matvec &a) {
  matvec res (a.size (), a.getRows (), a.getCols ());
  for (int i = 0; i < a.size (); i++) res.set (deg2rad (a.get (i)), i);
  return res;
}

} // namespace qucs
