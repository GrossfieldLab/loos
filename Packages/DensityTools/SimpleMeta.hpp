/*
  smetad.hpp

  (c) 2009 Tod D. Romo, Grossfield Lab, URMC


  Simple meta-data handler...

*/




#if !defined(SMETAD_HPP)
#define SMETAD_HPP


#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>


namespace lab {

  class SMetaData {
  public:

    typedef std::string                     value_type;
    typedef std::vector<value_type>         container_type;
    typedef container_type::iterator        iterator;
    typedef container_type::const_iterator  const_iterator;

    SMetaData() { }
    SMetaData(const std::string& s) { data_.push_back(s); }
    SMetaData(const std::vector<std::string>& v) : data_(v) { }

    container_type& data() { return(data_); }
    const container_type& data() const { return(data_); }

    // For convenience...
    iterator begin() { return(data_.begin()); }
    iterator end() { return(data_.end()); }
    const_iterator begin() const { return(data_.begin()); }
    const_iterator end() const { return(data_.end()); }
    bool empty() const { return(data_.empty()); }
    unsigned int size() const { return(data_.size()); }

    void clear() { data_.clear(); }
    void set(const std::string& s) { data_.clear(); data_.push_back(s); }
    void add(const std::string& s) { data_.push_back(s); }

    friend std::ostream& operator<<(std::ostream& os, const SMetaData& m) {
      for (const_iterator i = m.data_.begin(); i != m.data_.end(); ++i)
        os << "# " << *i << std::endl;
      return(os);
    }

    friend std::istream& operator>>(std::istream& is, SMetaData& m) {
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
    std::string stripper(const std::string& s) {
      unsigned int i;
      for (i=1; i < s.size() && s[i] == ' '; ++i) ;
      return(s.substr(i, s.size()));
    }

  private:
    std::vector<std::string> data_;
  };







};




#endif
