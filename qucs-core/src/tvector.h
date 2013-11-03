/*
 * tvector.h - simple vector template class definitions
 *
 * Copyright (C) 2004, 2005, 2006 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __TVECTOR_H__
#define __TVECTOR_H__

#include <vector>
#include <assert.h>
#include <limits>

#include <Eigen/Core>

#include "precision.h"


namespace qucs {



# if 0
template <typename nr_type_t>
using tvector = Eigen::Matrix<nr_type_t,Eigen::Dynamic,1>;

#else

template <class nr_type_t>
class tvector;

// Forward declarations of friend functions.
template <class nr_type_t>
nr_double_t maxnorm (tvector<nr_type_t>);
template <class nr_type_t>
nr_double_t norm (tvector<nr_type_t>);
template <class nr_type_t>
tvector<nr_type_t> operator + (tvector<nr_type_t>, tvector<nr_type_t>);
template <class nr_type_t>
tvector<nr_type_t> operator + (tvector<nr_type_t>, nr_type_t);
template <class nr_type_t>
tvector<nr_type_t> operator + (nr_type_t, tvector<nr_type_t>);
template <class nr_type_t>
tvector<nr_type_t> operator - (tvector<nr_type_t>, tvector<nr_type_t>);
template <class nr_type_t>
tvector<nr_type_t> operator * (const tvector<nr_type_t>&, const nr_type_t&);
template <class nr_type_t>
tvector<nr_type_t> operator * (const nr_type_t&, const tvector<nr_type_t>&);
template <class nr_type_t>
tvector<nr_type_t> operator * (tvector<nr_type_t>, tvector<nr_type_t>);
template <class nr_type_t>
tvector<nr_type_t> operator - (tvector<nr_type_t>);
template <class nr_type_t>
bool operator < (tvector<nr_type_t>, tvector<nr_type_t>);
template <class nr_type_t>
bool operator > (tvector<nr_type_t>, tvector<nr_type_t>);


template <class nr_type_t>
class tvector
{
 private:
  Eigen::Matrix<nr_type_t,Eigen::Dynamic,1> v;
 public:
  tvector () = default;
  tvector (const std::size_t n) = delete;
  tvector (const tvector &) = default;
  tvector (const Eigen::Matrix<nr_type_t,Eigen::Dynamic,1> &n):
    v(n) {};
  tvector (Eigen::Matrix<nr_type_t,Eigen::Dynamic,1> &&n):
    v(std::move(n)) {};
  ~tvector () = default;
  void set (int, nr_type_t) = delete;
  void setConstant (nr_type_t c) {
    this->v.setConstant(c);
  }
  void set (nr_type_t, int, int) = delete;
  void set (tvector, int, int) = delete;
  static tvector<nr_type_t> Zero(unsigned int row, unsigned int col) {
    return tvector(Eigen::Matrix<nr_type_t,Eigen::Dynamic,1>::Zero(row,1));
  }
  void setZero() { 
    if(this->size() > 0)
      for (unsigned int i = 0;i < this->size(); ++i)
	(*this)(i) = 0;
  }
  auto size (void) const -> decltype (v.size()) { return this->v.size(); }
  auto transpose() -> decltype(v.transpose()) { return this->v.transpose(); }
  auto adjoint() -> decltype(v.adjoint()) { return this->v.adjoint(); }
  tvector<nr_type_t> conjugate() const {
    auto res = tvector(this->v.conjugate());
    return res;
  }

  nr_type_t dot(const tvector<nr_type_t> &s) const {
    return this->v.dot(s.v);
  }
  
  nr_type_t sum() const {
    return this->v.sum();
  }

  nr_type_t * data (void) { return this->v.data(); }
  void exchangeRows (int, int);
  bool  isFinite (void) const {
    return this->v.isFinite();
  }
  void print ();
  void reorder (int *);
  int  contains (nr_type_t, nr_double_t eps = std::numeric_limits<nr_double_t>::epsilon());

  // some basic vector operations
#ifndef _MSC_VER
  friend tvector operator +<> (tvector, tvector);
  friend tvector operator -<> (tvector, tvector);
  friend tvector operator *<> (const tvector&, const nr_type_t&);
  friend tvector operator *<> (const nr_type_t&, const tvector&);
  friend tvector operator *<> (tvector, tvector);
  friend tvector operator -<> (tvector);
  friend tvector operator +<> (tvector, nr_type_t);
  friend tvector operator +<> (nr_type_t, tvector);
#endif

  // other operations
#ifndef _MSC_VER
  friend nr_double_t norm<> (tvector);
  friend nr_double_t maxnorm<> (tvector);
#endif

  // comparisons
#ifndef _MSC_VER
  friend bool operator < <> (tvector, tvector);
  friend bool operator > <> (tvector, tvector);
#endif

  // intrinsic operators
  tvector operator += (tvector);
  tvector operator -= (tvector);
  tvector operator *= (nr_double_t);

  // assignment operators
  tvector operator = (const nr_type_t) = delete;

  // easy accessor operators
  nr_type_t  operator () (int i) const {
    return this->v(i);
  }
  nr_type_t& operator () (int i) {
    return this->v(i); }
  nr_type_t  operator [] (int i) const {
    return this->v[i];
  }
  nr_type_t& operator [] (int i) {
    return this->v[i];
  }
};

  
} // namespace qucs
  
#include "tvector.cpp"
#endif

#endif /* __TVECTOR_H__ */
