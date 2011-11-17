// -------------------------------------------------
// Simple Metadata Handling
// -------------------------------------------------

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#if !defined(LOOS_SIMPLEMETA_HPP)
#define LOOS_SIMPLEMETA_HPP


#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>


namespace loos {

  namespace DensityTools {

    //! Simple class for handling metadata
    /**
     * Metadata consists of multiple lines that begin with a hash-mark
     * ('#').  When reading into a SimpleMeta object, the hash-marks are
     * stripped and each line becomes a string in a vector.  When
     * writing out, the process is reversed.
     */
    class SimpleMeta {
    public:

      typedef std::string                     value_type;
      typedef std::vector<value_type>         container_type;
      typedef container_type::iterator        iterator;
      typedef container_type::const_iterator  const_iterator;

      SimpleMeta() { }
      SimpleMeta(const std::string& s) { data_.push_back(s); }
      SimpleMeta(const std::vector<std::string>& v) : data_(v) { }

      //! Direct access to stored container of data
      container_type& data() { return(data_); }
      const container_type& data() const { return(data_); }

      //! Allow STL-iteration
      iterator begin() { return(data_.begin()); }
      iterator end() { return(data_.end()); }
      const_iterator begin() const { return(data_.begin()); }
      const_iterator end() const { return(data_.end()); }
      bool empty() const { return(data_.empty()); }
      unsigned int size() const { return(data_.size()); }

      //! Clear all contained metadata
      void clear() { data_.clear(); }

      //! Set metadata to string (deletes existing metadata)
      void set(const std::string& s) { data_.clear(); data_.push_back(s); }
    
      //! Append metadata
      void add(const std::string& s) { data_.push_back(s); }

      friend std::ostream& operator<<(std::ostream& os, const SimpleMeta& m) {
        for (const_iterator i = m.data_.begin(); i != m.data_.end(); ++i)
          os << "# " << *i << std::endl;
        return(os);
      }

      friend std::istream& operator>>(std::istream& is, SimpleMeta& m) {
        std::string buf;

        m.data_.clear();
        while (true) {
          int c = is.peek();
          if (c != '#')
            break;
          std::getline(is, buf);
          if (is.fail() || is.eof())
            throw(std::runtime_error("Error while reading metadata"));
          m.data_.push_back(m.stripper(buf));
        }
        return(is);
      }


    private:

      // Strip leading space (skipping meta-marker)
      std::string stripper(const std::string& s) {
        unsigned int i;
        for (i=1; i < s.size() && s[i] == ' '; ++i) ;
        return(s.substr(i, s.size()));
      }

    private:
      std::vector<std::string> data_;
    };




  };


};




#endif
