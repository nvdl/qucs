/*
 * flags.h - flags helper
 *
 * Copyright (C) 2015 Bastien ROUCARIÈS
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
#include <type_traits>
#include <bitset>
#include <string>
#include <ostream>
#include <istream>
#include <iostream>


namespace qucs {
  /*! Fail if not enum
     
      undefined implementation on purpose : fail if T is not an enum
  */
  template<typename T, class Enable = void>
  class flags; // undefined
}

/*! forward declaration */
namespace std {
  template <class T> struct hash<qucs::flags<T>  >
  {
    size_t operator()(const qucs::flags<T> & x) const;
  };
}

namespace qucs {
  /*! forward declaration */
  template<class charT, class traits, class enumT>
  std::basic_istream<charT, traits>&
  operator>> (std::basic_istream<charT,traits>& is, flags<enumT>& rhs);
  
  /*! Class for flags

      define an enum and declare as flags<enum> A;
  */
  template<typename T> class flags<T, typename std::enable_if<std::is_enum<T>::value>::type > {
  private:
    constexpr static std::size_t endenum = static_cast<std::size_t>(T::_end);
    std::bitset<endenum> bit;
  public:
    /*! emulation of pointer */
    typedef typename std::bitset<endenum>::reference reference;
    /*! default constructor */
    constexpr flags() noexcept = default;
    /*! default destructor */
    ~flags() noexcept = default;
    /*! construct from std::string */
    template <class charT, class traits, class Alloc>
    explicit flags(const std::basic_string<charT,traits,Alloc>& str,
		   typename std::basic_string<charT,traits,Alloc>::size_type pos = 0,
		   typename std::basic_string<charT,traits,Alloc>::size_type n = std::basic_string<charT,traits,Alloc>::npos,
		   charT zero = charT('0'), charT one = charT('1')) {
      this->bits(str,pos,n,zero,one);
    };
    /*! construct from cstring */
    template <class charT>
    explicit flags (const charT* str,
		    typename std::basic_string<charT>::size_type n = std::basic_string<charT>::npos,
		    charT zero = charT('0'), charT one = charT('1')) {
      this->bits(str,n,zero,one);
    }
    /*! equality */
    bool operator== (const flags& rhs) const noexcept { return this->bit==rhs.bit; };
    /*! not equality */
    bool operator!= (const flags& rhs) const noexcept { return this->bit!=rhs.bit; };
    /*! access */
    constexpr bool operator[] (T pos) const {
      return this->bit[static_cast<std::size_t>(pos)];
    };
    /*! write access */
    reference operator[] (T pos) {
      return this->bit[static_cast<std::size_t>(pos)];
    };
    /*! count number of flag set */
    size_t count() const noexcept {
      return this->bit.count();
    };
    /*! size */
    constexpr size_t size() noexcept {
      return this->bit.size();
    };
    /*! test with exception if out of bound */
    bool test (T pos) const {
      return this->bit.test(static_cast<std::size_t>(pos));
    };
    /*! check if one flag is set */
    bool any() const noexcept {
      return this->bit.any();
    };
    /*! check if none flag is set */
    bool none() const noexcept {
      return this->bit.none();
    };
    /*! check if all flags are set */
    bool all() const noexcept {
      return this.bit.all();
    };
    /*! set all flags */
    flags<T>& set(void) noexcept {
      this->bit.set();
      return *this;
    };
    /*! set one flag (may raise an exception in out of bound */
    flags& set (T pos, bool val = true) {
      this->bit.set(static_cast<std::size_t>(pos),val);
      return *this;
    };
    /*! reset all flags */
    flags<T>& reset() noexcept {
      this->bit.reset();
      return *this;
    }
    /*! reset one flags */
    flags<T>& reset (T pos) {
      this->bit.reset(static_cast<std::size_t>(pos));
      return *this;
    }
    /*! flip all bit */
    flags<T>& flip() noexcept {
      this->bit.flip();
      return *this;
    }
    /*! flipt one bit */
    flags<T>& flip (T pos) {
      this->bit.flip(static_cast<std::size_t>(pos));
      return *this;
    };

    /*! convert to string */
    template <class charT = char,
	      class traits = std::char_traits<charT>,
	      class Alloc = std::allocator<charT> >
    std::basic_string<charT,traits,Alloc> to_string (charT zero = charT('0'),
						charT one  = charT('1')) const {
      return this->bit.to_string(zero,one);
    }
    /* friend function */
    template<class charT, class traits, class enumT>
    friend
    std::basic_istream<charT, traits>&
    operator>> (std::basic_istream<charT,traits>& is, flags<enumT>& rhs);    
    template<class charT, class traits, class enumT>
    friend
    std::basic_ostream<charT, traits>&
    operator<< (std::basic_ostream<charT,traits>& os, const flags<enumT>&  rhs);
    template<class enumT>
    friend std::size_t std::hash<flags<enumT> >::operator ()(const flags<enumT>&) const;
   
  };
  
  template<class charT, class traits, class enumT>
  std::basic_istream<charT, traits>&
  operator>> (std::basic_istream<charT,traits>& is, flags<enumT>& rhs) {
    return is >> rhs.bit;
  }
  template<class charT, class traits, class enumT>
  std::basic_ostream<charT, traits>&
  operator<< (std::basic_ostream<charT,traits>& os, const flags<enumT>&  rhs) {
    return os << rhs.bit;
  }
}

/*! flags specialisation */
namespace std {
  template <class T> size_t hash<qucs::flags<T>>::operator()(const qucs::flags<T> & x) const
  {
    return std::hash<std::bitset<static_cast<size_t>(T::_end)> >()(x.bit);
  }
}
