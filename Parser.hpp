/*
  Parser.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Parser for atom selections using Bison/Flex

*/


#if !defined(PARSER_HPP)
#define PARSER_HPP

#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

#include <AtomicGroup.hpp>
#include <Kernel.hpp>
#include <ParserDriver.hpp>



class Parser {
  loos::Kernel krnl;
  ParserDriver driver;

public:
  Parser(const string& s) : driver(s, krnl) { }


  Kernel& kernel(void) { return(krnl); }
};





#endif

