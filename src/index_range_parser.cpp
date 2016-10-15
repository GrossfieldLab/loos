/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2016, Tod D. Romo, Alan Grossfield
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




#include <loos_defs.hpp>
#include <exceptions.hpp>

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

#include <boost/phoenix/object/new.hpp>

#include <boost/format.hpp>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>



using namespace std;

namespace loos {

  namespace internal {


    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    namespace phoenix = boost::phoenix;

    struct RangeItem {
      uint start;
      uint stop;
      int step;


      bool validate() const {
        if (start < stop && step > 0)
          return false;
        if (start < stop && step < 0)
          return false;

	return true;
      }


      vector<uint> generate() const {
        vector<uint> results;

        if (step == 0)
          results.push_back(start);
        else if (step > 0) {
          for (uint i=start; i<=stop; i += step)
            results.push_back(i);
        } else {
          uint i;
          for (i=start; i>= stop && i<=start; i += step)
            results.push_back(i);
        }

        return results;
      }
    
      RangeItem() : start(0), stop(0), step(0) { }
      RangeItem(const uint a) : start(a), stop(a), step(0) { }
      RangeItem(const uint a, const uint b) : start(a), stop(b), step(1) { }
      RangeItem(const uint a, const uint b, const int c) : start(a), stop(b), step(c) { }
    };



    typedef vector<RangeItem*>   RangeList;
    typedef RangeItem*           RangeItemPointer;
  
    template<typename Iterator>
    struct RangeParser : qi::grammar<Iterator, RangeList(), ascii::space_type>
    {
      RangeParser() : RangeParser::base_type(start), maxsize(0) {
        using qi::int_;
        using qi::uint_;
        using qi::_1;
        using qi::_2;
        using qi::_3;
        using qi::_val;
        using phoenix::ref;
        using qi::eps;


        endpoint = uint_ [ _val = _1 ]
          | eps [ _val = ref(maxsize) ];
      
        ranger = (uint_ >> ':' >> int_ >> ':' >> endpoint ) [ _val = phoenix::new_<RangeItem>(_1,_3, _2) ]
          | ( endpoint >> ':' >> int_ >> ':' >> uint_ ) [ _val = phoenix::new_<RangeItem>(_1,_3, _2) ]
          | ( uint_ >> ':' >> endpoint ) [ _val = phoenix::new_<RangeItem>(_1, _2) ]
          | ( endpoint >> ':' >> uint_ ) [ _val = phoenix::new_<RangeItem>(_1, _2, -1) ]
          | ( uint_ ) [ _val = phoenix::new_<RangeItem>(_1) ];

        start = ranger % ',';
      }

      uint maxsize;
    
      qi::rule<Iterator, uint, ascii::space_type> endpoint;
      qi::rule<Iterator, RangeItemPointer, ascii::space_type> ranger;
      qi::rule<Iterator, RangeList(), ascii::space_type> start;

    };


  } // internal
  

  vector<uint> parseIndexRange(const std::string& input, const uint maxsize = 0) {
    using boost::spirit::ascii::space;
    typedef std::string::const_iterator iterator_type;
    typedef internal::RangeParser<iterator_type> range_parser;

    range_parser parser;
    parser.maxsize = maxsize;

    vector<internal::RangeItemPointer> rangelist;
    string::const_iterator iter = input.begin();
    string::const_iterator end = input.end();

    vector<uint> result;
    bool ok = phrase_parse(iter, end, parser, space, rangelist);
    if (ok && iter == end) {
      for (uint i=0; i<rangelist.size(); ++i) {
        vector<uint> r = rangelist[i]->generate();
        copy(r.begin(), r.end(), back_inserter(result));
        delete rangelist[i];
      }

    } else {
      throw(ParseError("Could not parse range: " + input));
    }
    
    return result;
  }


}  // loos
