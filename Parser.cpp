/*
  Parser.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Parser for atom selections...
*/

#include <sstream>

#include "Parser.hpp"

using namespace std;


#define GPTOKS(x)    (boost::get<0>((x)))
#define GLTOKS(x)    (boost::get<1>((x)))

namespace loos {

  namespace parser {

    void dump(const pares& p, const string s) {
      cout << s << endl;
      ptoks q = GPTOKS(p);
      Tokens r = GLTOKS(p);

      cout << "  o Parse tokens o\n";
      prptoks(q);
      cout << "  o Remaining lex tokens o\n";
      cout << r;
    }

    int max3(const int x, const int y, const int z) {
      int d = (x > y) ? x : y;
      d = (d > z) ? d : z;
      return(d);
    }

    ptok mktok(const ptok::ptok_type t, const string s = "") {
      ptok p;

      p.type = t;
      p.val = s;
      return(p);
    }

    void prptoks(const ptoks& p) {
      ptoks::const_iterator i;
      for (i = p.begin(); i != p.end(); i++)
	cout << (*i) << endl;
    }

    pares mknull(const pares& p) {
      Tokens t = GLTOKS(p);
      ptoks r;
      pares q(r, t);
      return(q);
    }


    // p may ref a temp...  Problem???
    pares mkpares(const ptok& p, const Tokens& t) {
      ptoks r;
      r.push_back(p);
      pares q(r, t);
      return(q);
    }

    int ptokn(const pares& p) {
      ptoks r = GPTOKS(p);
      return(r.size());
    }

    int depth(pares *p, pares *q=0, pares *r=0) {
      int d = ptokn(*p);
      if (d == 0)
	return(0);

      if (q != 0) {
	int e = ptokn(*q);
	if (e == 0)
	  return(0);
	d += e;

	if (r != 0) {
	  e = ptokn(*r);
	  if (e == 0)
	    return(0);
	  d += e;
	}
      }
      return(d);
    }

    ptok gptok(const pares& p) {
      ptoks r = GPTOKS(p);
      return(r.front());
    }


    // p : q  | q
    pares twoway(const pares& p, const pares& q) {
      ptoks pt = GPTOKS(p);
      ptoks qt = GPTOKS(q);
      Tokens lt = GLTOKS(q);

      qt.insert(qt.end(), pt.begin(), pt.end());
      pares res(pt, lt);
      return(res);
    }

    // p : r : q
    pares threeway(const pares& p, const pares& q, const pares& r) {
      ptoks pt = GPTOKS(p);
      ptoks qt = GPTOKS(q);
      ptoks rt = GPTOKS(r);
      Tokens lt = GLTOKS(r);

      pt.insert(pt.end(), rt.begin(), rt.end());
      pt.insert(pt.end(), qt.begin(), qt.end());
      pares res(pt, lt);
      return(res);
    }

    pares zap(const pares& p) {
      ptoks notoks;
      Tokens toks = GLTOKS(p);
      pares res(notoks, toks);
      return(res);
    }



    pares numid(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);


      if (ltoks.empty())
	return(null);

      Token t = ltoks.pop();
      if (t.type == Token::ID) {
	ptok p;

	if (t.datum == "id")
	  p = mktok(ptok::id);
	else if (t.datum == "resid")
	  p = mktok(ptok::resid);
	else
	  return(null);

	return(mkpares(p, ltoks));
      }

      return(null);
    }



    pares alphid(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);

      if (ltoks.empty())
	return(null);

      Token t = ltoks.pop();
      if (t.type == Token::ID) {
	ptok p;

	if (t.datum == "name")
	  p = mktok(ptok::name);
	else if (t.datum == "resname")
	  p = mktok(ptok::resname);
	else if (t.datum == "segid")
	  p = mktok(ptok::segid);
	else
	  return(null);

	return(mkpares(p, ltoks));
      }

      return(null);
    }


    pares number(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);

      if (ltoks.empty())
	return(null);

      Token t = ltoks.pop();
      if (t.type == Token::NUMERIC) {
	ptok p = mktok(ptok::num, t.datum);
	return(mkpares(p, ltoks));
      }

      
      return(null);
    }



    pares alpha(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);

      if (ltoks.empty())
	return(null);

      Token t = ltoks.pop();
      if (t.type == Token::STRING) {
	ptok p = mktok(ptok::alpha, t.datum);
	return(mkpares(p, ltoks));
      }

      
      return(null);
    }



    pares numeric(const pares p) {

      pares p1 = numid(p);
      pares p2 = number(p);

      if (ptokn(p1) >= ptokn(p2))
	return(p1);
      
      return(p2);
    }




    pares alphabetic(const pares p) {

      pares p1 = alphid(p);
      pares p2 = alpha(p);

      if (ptokn(p1) >= ptokn(p2))
	return(p1);
      
      return(p2);
    }


    pares literal(const pares p) {

      pares p1 = alphabetic(p);
      pares p2 = numeric(p);

      if (ptokn(p1) >= ptokn(p2))
	return(p1);
      
      return(p2);
    }

    pares unop(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);

      if (ltoks.empty())
	return(null);

      Token t = ltoks.pop();
      if (t.type == Token::OPERATOR) {
	ptok p;
	if (t.datum == "!")
	  p = mktok(ptok::lnot);
	else
	  return(null);
	return(mkpares(p, ltoks));
      }
      return(null);
    }


    pares relop(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);

      if (ltoks.empty())
	return(null);

      Token t = ltoks.pop();
      if (t.type == Token::OPERATOR) {
	ptok p;
	if (t.datum == "<")
	  p = mktok(ptok::lt);
	else if (t.datum == "<=")
	  p = mktok(ptok::lte);
	else if (t.datum == ">=")
	  p = mktok(ptok::gte);
	else if (t.datum == ">")
	  p = mktok(ptok::gt);
	else if (t.datum == "==")
	  p = mktok(ptok::eq);
	else if (t.datum == "!=")
	  p = mktok(ptok::neq);
	else
	  return(null);
	return(mkpares(p, ltoks));
      }
      return(null);
    }


    pares binop(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);

      if (ltoks.empty())
	return(null);

      Token t = ltoks.pop();
      if (t.type == Token::OPERATOR) {
	ptok p;
	if (t.datum == "&&")
	  p = mktok(ptok::land);
	else if (t.datum == "||")
	  p = mktok(ptok::lor);
	else
	  return(null);
	return(mkpares(p, ltoks));
      }
      return(null);
    }

    pares regexp(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);

      if (ltoks.empty())
	return(null);

      pares p1 = alphid(p);
      if (ptokn(p1) == 0)
	return(null);

      ltoks = GLTOKS(p1);
      Token t = ltoks.pop();
      if (t.type == Token::OPERATOR && t.datum == "=~") {
	ptoks notoks;
	pares p2 = alpha(pares(notoks, ltoks));
	if (ptokn(p2) > 0) {
	  pares q = mkpares(mktok(ptok::regex), ltoks);
	  return(threeway(p1, q, p2));
	}
      }
      return(null);
    }


    pares rexp(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);

      if (ltoks.empty())
	return(null);

      pares p1 = regexp(p);
      
      pares p2a = literal(p);
      pares p2b = relop(zap(p2a));
      pares p2c = literal(zap(p2b));

      if (ptokn(p1) > 0)
	return(p1);
      if (ptokn(p2a) >0 && ptokn(p2b) >0 && ptokn(p2c) >0)
	return(threeway(p2a, p2b, p2c));
      
      return(null);
    }


    pares lexp(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);

      if (ltoks.empty())
	return(null);

    
      pares p1 = rexp(zap(p));

      pares p2a = unop(zap(p));
      pares p2b = rexp(zap(p2a));

      pares p3a = rexp(zap(p));
      pares p3b = binop(zap(p3a));
      pares p3c = rexp(zap(p3b));

      int d1 = depth(&p1);
      int d2 = depth(&p2a, &p2b);
      int d3 = depth(&p3a, &p3b, &p3c);

      int m = max3(d1, d2, d3);
      
      if (m == d1)
	return(p1);
      else if (m == d2)
	return(twoway(p2a, p2b));

      return(threeway(p3a, p3b, p3c));

    }


    pares sexp(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);

      if (ltoks.empty())
	return(null);

      Token t = ltoks.pop();
      if (t.type == Token::LPAR) {
	ptoks notoks;
	pares q(notoks, ltoks);
	pares r = lexp(q);
	if (ptokn(r) > 0) {
	  ltoks = GLTOKS(r);
	  t = ltoks.pop();
	  if (t.type == Token::RPAR) {
	    ptoks toks = GPTOKS(r);
	    pares res(toks, ltoks);
	    return(res);
	  }
	}
      }
      
      return(lexp(p));
    }


    pares expr(const pares p) {
      Tokens ltoks = GLTOKS(p);
      pares null = mknull(p);

      if (ltoks.empty())
	return(null);

    
      pares p1 = sexp(zap(p));

      pares p2a = unop(zap(p));
      pares p2b = sexp(zap(p2a));

      pares p3a = sexp(zap(p));
      pares p3b = binop(zap(p3a));
      pares p3c = sexp(zap(p3b));

      int d1 = depth(&p1);
      int d2 = depth(&p2a, &p2b);
      int d3 = depth(&p3a, &p3b, &p3c);

      int m = max3(d1, d2, d3);
      
      if (m == d1)
	return(p1);
      else if (m == d2)
	return(twoway(p2a, p2b));

      return(threeway(p3a, p3b, p3c));

    }



    
  };

};


