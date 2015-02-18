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


/* The function returns the vector specified by the given matrix
   indices.  If the matrix vector has a valid name 'A' the returned
   vector gets the name 'A[r,c]'. */
qucs::vector matvec::get (int r, int c) {
  assert (r >= 0 && r < this->rows() && c >= 0 && c < this->cols());
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
	  for(matvec::size_type ii=0;ii<(*mv).size();ii++)
	    ((*mv)[ii])(r,c)=(*v)(ii);
          free (n);
        }
      }
    }
    return mv;
  }
  return NULL;
}

// Matrix vector addition.
matvec operator + (const matvec &a, const matvec &b) {
  assert (a.rows () == b.rows () && a.cols () == b.cols () &&
	  a.size () == b.size ());
  matvec res (a);
  for (matvec::size_type i = 0; i < a.size (); i++)
    res[i]+=b[i];
  return res;
}

// Matrix vector addition with single matrix.
matvec operator + (const matvec &a, const matrix &b) {
  assert (a.rows () == b.rows () && a.cols () == b.cols ());
  matvec res (a);
  std::for_each(res.begin(),res.end(),
		[b](matrix &el) { el += b; });
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
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[z] (matrix& el) {
		  el = (el.array()+z).matrix();
		});
  return res;
}

// Matrix vector scalar addition in different order.
matvec operator + (nr_complex_t z, const matvec &a) {
  return a+z;
}

// Matrix vector scalar addition.
matvec operator + (const matvec &a, nr_double_t d) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[d] (matrix& el) {
		  el = (el.array()+d).matrix();
		});
  return res;
}

// Matrix vector scalar addition in different order.
matvec operator + (nr_double_t d, const matvec &a) {
  return a+d;
}

// Matrix vector scalar subtraction.
matvec operator - (const matvec &a, nr_complex_t z) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[z] (matrix& el) {
		  el = (el.array()-z).matrix();
		});
  return res;
}

// Matrix vector scalar subtraction in different order.
matvec operator - (nr_complex_t z, const matvec &a) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[z] (matrix& el) {
		  el = (z-el.array()).matrix();
		});
  return res;
}

// Matrix vector scalar subtraction.
matvec operator - (const matvec &a, nr_double_t d) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[d] (matrix& el) {
		  el = (el.array()-d).matrix();
		});
  return res;
}

// Matrix vector scalar subtraction in different order.
matvec operator - (nr_double_t d, const matvec &a) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[d] (matrix& el) {
		  el = (d-el.array()).matrix();
		});
  return res;
}

// Intrinsic matrix vector addition.
matvec matvec::operator += (const matvec &a) {
  assert (a.rows () == this->rows() && a.cols () == this->cols() &&
	  a.size () == data.size());
  for (matvec::size_type i = 0; i < data.size(); i++)
    data[i] += a[i];
  return *this;
}

// Matrix vector subtraction.
matvec operator - (const matvec &a, const matvec &b) {
  assert (a.rows () == b.rows () && a.cols () == b.cols () &&
	  a.size () == b.size ());
  matvec res (a);
  for (matvec::size_type i = 0; i < a.size (); i++)
    res[i] -= b[i];
  return res;
}

// Matrix vector subtraction with single matrix.
matvec operator - (const matvec& a, const matrix& b) {
  assert (a.rows () == b.rows () && a.cols () == b.cols ());
  matvec res (a);
  for (matvec::size_type i = 0; i < a.size (); i++)
    res[i] -= b;
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
  std::for_each(std::begin(*this),std::end(*this),
		[] (matrix& el) {
		  el = -el;
		});
  return *this;
}

// Intrinsic matrix vector subtraction.
matvec matvec::operator -= (const matvec &a) {
  assert (a.rows () == this->rows() && a.cols () == this->cols() &&
	  a.size () == data.size());
  for (matvec::size_type i = 0; i < a.size (); i++)
    data[i] -= a[i];
  return *this;
}

// Matrix vector scaling.
matvec operator * (const matvec &a, nr_complex_t z) {
  matvec res (a);
  std::for_each(std::begin(res),std::end(res),
		[z] (matrix& el) {
		  el = (z*el.array()).matrix();
		});
  return res;
}

// Matrix vector scaling in different order.
matvec operator * (nr_complex_t z, const matvec &a) {
  return a * z;
}

// Scalar matrix vector scaling.
matvec operator * (const matvec &a, nr_double_t d) {
  matvec res (a);
  std::for_each(std::begin(res),std::end(res),
		[d] (matrix& el) {
		  el = (d*el.array()).matrix();
		});
  return res;
}

// Scalar matrix vector scaling in different order.
matvec operator * (nr_double_t d,const matvec &a) {
  return a * d;
}

// Matrix vector scaling by a second vector.
matvec operator * (const matvec &a, const qucs::vector &b) {
  assert (a.size () == b.getSize ());
  matvec res (a);
  for (matvec::size_type i = 0; i < a.size (); i++)
    res[i]*=b.get(i);
  return res;
}

// Matrix vector scaling by a second vector in different order.
matvec operator * (const qucs::vector &a, const matvec &b) {
  return b * a;
}

// Matrix vector scaling.
matvec operator / (const matvec &a, const nr_complex_t z) {
  return (1.0/z)*a;
}

// Scalar matrix vector scaling.
matvec operator / (const matvec &a, nr_double_t d) {
  return (1.0/d)*a;
}

// Matrix vector scaling by a second vector.
matvec operator / (const matvec& a, const qucs::vector &b) {
  assert (a.size () == b.getSize ());
  matvec res (a);
  return (1.0/b)*res;
}

// Matrix vector multiplication.
matvec operator * (const matvec &a, const matvec &b) {
  assert (a.cols () == b.rows () && a.size () == b.size ());
  matvec res (a);
  for (matvec::size_type i = 0; i < a.size (); i++)
    res[i]*=b[i];
  return res;
}

// Matrix vector multiplication with a single matrix.
matvec operator * (const matvec &a, const matrix &b) {
  assert (a.cols () == b.rows ());
  matvec res (a);
  std::for_each(std::begin(res),std::end(res),
		[b] (matrix& el) {
		  el *= b;
		});
  return res;
}

// Matrix vector multiplication with a single matrix in different order.
matvec operator * (const matrix &a,const matvec &b) {
  assert (a.cols () == b.rows ());
  matvec res (b);
  std::for_each(std::begin(res),std::end(res),
		[a] (matrix& el) {
		  el = a*el;
		});
  return res;
}

// Compute determinants of the given matrix vector.
qucs::vector det (const matvec &a) {
  qucs::vector res(a.size());
  for(matvec::size_type i=0;i<a.size();i++)
    res(i)=a[i].determinant();
  return res;
}

// Compute inverse matrices of the given matrix vector.
matvec inverse (const matvec &a) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[] (matrix& el) {
		  el = el.inverse();
		});
  return res;
}

// Compute inverse matrices of the given matrix vector.
matvec sqr (const matvec &a) {
  return a * a;
}

// Compute n-th power of the given matrix vector.
matvec pow (const matvec &a, int n) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[n] (matrix& el) {
		  el = pow(el,n);
		});
  return res;
}

// Compute n-th powers in the vector of the given matrix vector.
matvec pow (const matvec &a, const qucs::vector &v) {
  assert (a.size () == v.getSize ());
  matvec res (a);
  for (matvec::size_type i = 0; i < a.size (); i++)
    res[i]=pow(a[i],static_cast<int>(v(i).real()));
  return res;
}

// Conjugate complex matrix vector.
matvec conj (const matvec &a) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[] (matrix& el) {
		  el = el.conjugate();
		});
  return res;
}

// Computes magnitude of each matrix vector element.
matvec abs (const matvec &a) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[] (matrix& el) {
		  el = el.cwiseAbs().cast<nr_complex_t>();
		});
  return res;
}

// Computes magnitude in dB of each matrix vector element.
matvec dB (const matvec &a) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[] (matrix& el) {
		  el = dB(el);
		});
  return res;
}

// Computes the argument of each matrix vector element.
matvec arg (const matvec &a) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[] (matrix& el) {
		  el = arg(el);
		});
  return res;
}

// Real part matrix vector.
matvec real (const matvec &a) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[] (matrix& el) {
		  el = el.real().cast<nr_complex_t>();
		});
  return res;
}

// Real part matrix vector.
matvec imag (const matvec &a) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[] (matrix& el) {
		  el = el.imag().cast<nr_complex_t>();
		});
  return res;
}

/* The function returns the adjoint complex matrix vector.  This is
   also called the adjugate or transpose conjugate. */
matvec adjoint (const matvec &a) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[] (matrix& el) {
		  el = el.adjoint();
		});
  return res;
}

// Transpose the matrix vector.
matvec transpose (const matvec &a) {
  matvec res(a);
  std::for_each(std::begin(res),std::end(res),
		[] (matrix& el) {
		  el = el.transpose();
		});
  return res;
}

/* Convert scattering parameters with the reference impedance 'zref'
   to scattering parameters with the reference impedance 'z0'. */
matvec stos (const matvec &s, const qucs::vector &zref, const qucs::vector &z0) {
  assert (s.cols () == s.rows () &&
	   s.cols () == zref.getSize () && s.cols () == z0.getSize ());
  matvec res (s);
  for (matvec::size_type i = 0; i < s.size (); i++)
    res[i]=stos (s[i], zref, z0);
  return res;
}

matvec stos (const matvec &s, nr_complex_t zref, nr_complex_t z0) {
  int d = s.rows ();
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
  assert (s.cols () == s.rows () && s.cols () == z0.getSize ());
  matvec res (s);
  for (matvec::size_type i = 0; i < s.size (); i++)
    res[i]=stoy (s[i], z0);
  return res;
}

matvec stoy (const matvec &s, nr_complex_t z0) {
  return stoy (s, qucs::vector (s.cols (), z0));
}

// Convert admittance matrix to scattering parameter matrix vector.
matvec ytos (const matvec &y, const qucs::vector &z0) {
  assert (y.cols () == y.rows () && y.cols () == z0.getSize ());
  matvec res (y);
  for (matvec::size_type i = 0; i < y.size (); i++)
    res[i] = ytos (y[i], z0);
  return res;
}

matvec ytos (const matvec &y, nr_complex_t z0) {
  return ytos (y, qucs::vector (y.cols (), z0));
}

// Convert scattering parameters to impedance matrix vector.
matvec stoz (const matvec &s, const qucs::vector &z0) {
  assert (s.cols () == s.rows () && s.cols () == z0.getSize ());
  matvec res (s);
  for (matvec::size_type i = 0; i < s.size (); i++)
    res[i] =stoz (s[i], z0);
  return res;
}

matvec stoz (const matvec &s, nr_complex_t z0) {
  return stoz (s, qucs::vector (s.cols (), z0));
}

// Convert impedance matrix vector scattering parameter matrix vector.
matvec ztos (const matvec &z, const qucs::vector &z0) {
  assert (z.cols () == z.rows () && z.cols () == z0.getSize ());
  matvec res (z);
  for (matvec::size_type i = 0; i < z.size (); i++)
    res[i] = ztos (z[i], z0);
  return res;
}

matvec ztos (const matvec &z, nr_complex_t z0) {
  return ztos (z, qucs::vector (z.cols (), z0));
}

// Convert impedance matrix vector to admittance matrix vector.
matvec ztoy (const matvec &z) {
  assert (z.cols () == z.rows ());
  matvec res (z);
  for (matvec::size_type i = 0; i < z.size (); i++)
    res[i] = ztoy (z[i]);
  return res;
}

// Convert admittance matrix vector to impedance matrix vector.
matvec ytoz (const matvec &y) {
  assert (y.cols () == y.rows ());
  matvec res (y);
  for (matvec::size_type i = 0; i < y.size (); i++)
    res[i] = ytoz (y[i]);
  return res;
}

/* This function converts 2x2 matrix vectors from any of the matrix
   forms Y, Z, H, G and A to any other.  Also converts S<->(A, T, H, Y
   and Z) matrix vectors. */
matvec twoport (const matvec &m, char in, char out) {
  assert (m.cols () >= 2 && m.rows () >= 2);
  matvec res (m);
  for (matvec::size_type i = 0; i < m.size (); i++)
    res[i] = twoport (m[i],in,out);
  return res;
}

/* The function returns the Rollet stability factor vector of the
   given S-parameter matrix vector. */
qucs::vector rollet (const matvec &m) {
  assert (m.cols () >= 2 && m.rows () >= 2);
  qucs::vector res (m.size ());
  for (matvec::size_type i = 0; i < m.size (); i++)
    res(i) = rollet(m[i]);
  return res;
}

/* The function returns the stability measure B1 vector of the given
   S-parameter matrix vector. */
qucs::vector b1 (const matvec &m) {
  assert (m.cols () >= 2 && m.rows () >= 2);
  qucs::vector res (m.size ());
  for (matvec::size_type i = 0; i < m.size (); i++)
    res(i) = b1 (m[i]);
  return res;
}

matvec rad2deg (const matvec &a) {
  matvec res (a);
  for (matvec::size_type i = 0; i < a.size (); i++)
    res[i] = rad2deg (a[i]);
  return res;
}

matvec deg2rad (const matvec &a) {
  matvec res (a);
  for (matvec::size_type i = 0; i < a.size (); i++)
    res[i] = deg2rad (a[i]);
  return res;
}

} // namespace qucs
