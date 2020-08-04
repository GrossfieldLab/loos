/*
 *  Compute membrane property distribution about a protein
 *  Can compute chain molecular order parameter, tilt vector, or density
 *
 *  Alan Grossfield
 *  (adapted from some code by Josh Horn, adapted from some stuff I wrote)
 */

#include <boost/lexical_cast.hpp>
#include <sstream>
#include "loos.hpp"

//* Virtual base class for CalcProperty.  Needed because CalcProperty is
//  a template, and you can't instantiate a template, only a particular version
//  of it.
class CalcPropertyBase
{
public:
    virtual void normalize(uint frames)
        {
        }

    virtual void calc(const loos::AtomicGroup &group,
                      const uint xbin, const uint ybin)
        {
        }

    virtual void set(const uint xbin, const uint ybin, const double val)
        {
        }

    virtual void set(const uint xbin, const uint ybin, const loos::GCoord val)
        {
        }

    virtual void incr(const uint xbin, const uint ybin, const double val)
        {
        }

    virtual void incr(const uint xbin, const uint ybin, const loos::GCoord val)
        {
        }

    //* Need to use a pointer to val instead of a return val so
    //  that the template below can overload it.
    //  NOTE: if you're going to create a new CalcProperty subclass
    //        that returns something other than a double or a GCoord, you'll
    //        need to add new virtual version of all of the methods here
    virtual void get(const uint xbin, const uint ybin, double *val)
        {
        }

    virtual void get(const uint xbin, const uint ybin, loos::GCoord *val)
        {
        }

    virtual const uint get_norm(const uint xbin, const uint ybin)
        {
        uint i = 1;
        return(i);
        }

    virtual const std::string print(const uint xbin, const uint ybin)
        {
        std::string s(" ");
        return(s);
        }

};

template <class T> class CalcProperty : public CalcPropertyBase
{
protected:
    uint _xbins, _ybins;
public:
    CalcProperty( uint xbins, uint ybins ):
                _xbins(xbins),
                _ybins(ybins),
                _storage(xbins*ybins),
                _norm(xbins*ybins)
        {
        }

    // This will usually be ok, but will occassionally
    // get overridden by subclasses, eg ones that return
    // a density
    void normalize (uint frames)
        {
        T val;
        for (uint i=0; i< _xbins; ++i)
            {
            for (uint j=0; j< _ybins; ++j)
                {
                get(i, j, &val);
                uint norm = get_norm(i,j);
                if (norm)
                    {
                    set(i,j, val /get_norm(i,j));
                    }
                else
                    set(i,j, -999999999.);
                }
            }
        }

    // This will get overridden by subclasses
    void calc(const loos::AtomicGroup &group, const uint xbin, const uint ybin)
        {
        incr(xbin, ybin, static_cast<T>(1.0));
        }

    void set(const uint xbin, const uint ybin, const T val)
        {
        uint index= xbin * _ybins + ybin;
        _storage[index] = val;
        }

    void incr(const uint xbin, const uint ybin, const T val)
        {
        uint index= xbin * _ybins + ybin;
        _storage[index] += val;
        _norm[index]++;
        }

    void get(const uint xbin, const uint ybin, T *val)
        {
        uint index= xbin * _ybins + ybin;
        *val = _storage[index];
        }

    const uint get_norm(const uint xbin, const uint ybin)
        {
        uint index= xbin * _ybins + ybin;
        return(_norm[index]);
        }

    const std::string print(const uint xbin, const uint ybin)
        {
        T val;
        get(xbin, ybin, &val);
        return(boost::lexical_cast<std::string>(val));
        }

private:
    std::vector<T> _storage;
    std::vector<uint> _norm;

};

//* calculate the density distribution, in groups/Ang^2
class CalcDensity : public CalcProperty<double>
{
public:
    CalcDensity( uint xbins, uint ybins, double xwidth, double ywidth) :
            CalcProperty<double>(xbins, ybins)
        {
        _bin_area = xwidth * ywidth;
        }

    void normalize (uint frames)
        {
        double norm =  _bin_area * frames;
        double val;
        for (uint i=0; i< this->_xbins; ++i)
            {
            for (uint j=0; j< this->_ybins; ++j)
                {
                this->get(i,j, &val);
                this->set(i, j, val / norm);
                }
            }
        }

private:
    double _bin_area;
};

//* Calculate molecular order parameter
class CalcMolOrder : public CalcProperty<double>
{
public:
   CalcMolOrder (uint xbins, uint ybins) : CalcProperty<double>(xbins, ybins)
        {
        }

   void calc(const loos::AtomicGroup &group, const uint xbin, const uint ybin)
        {
        std::vector<loos::GCoord> axes = group.principalAxes();
        loos::GCoord ave = 0.5* (axes[1] + axes[2]);
        double cosine = ave.z() / ave.length();
        double val = fabs(1.5*cosine*cosine - 0.5);
        incr(xbin, ybin, val);
        }

};

//* Calculate height
class CalcHeight : public CalcProperty<double>
{
public:
   CalcHeight (uint xbins, uint ybins) : CalcProperty<double>(xbins, ybins)
        {
        }

   // This implicitly assumes the membrane center is z=0
   void calc(const loos::AtomicGroup &group, const uint xbin, const uint ybin)
        {
        loos::GCoord centroid = group.centroid();
        incr(xbin, ybin, centroid.z());
        }

};

//* Calculate the in-plane "orientation field" for the group
class CalcOrientVector : public CalcProperty<loos::GCoord>
{
public:
    CalcOrientVector(uint xbins, uint ybins) :
        CalcProperty<loos::GCoord>(xbins, ybins)
        {
        }

    void calc(const loos::AtomicGroup &group, const uint xbin, const uint ybin)
        {
        std::vector<loos::GCoord> axes = group.principalAxes();
        // Force a consistent sign convention on the principal axis
        // by insisting it point toward the center of the membrane.
        // So, if the molecule is in the +z leaflet, the axis must point
        // "downward"
        loos::GCoord centroid = group.centroid();
        if (axes[0].z()*centroid.z() > 0)
            {
            axes[0] *= -1;
            }
        axes[0].z() = 0.0;
        incr(xbin, ybin, axes[0]);
        }

    const std::string print(const uint xbin, const uint ybin)
        {
        std::stringstream s;
        loos::GCoord tmp;
        get(xbin, ybin, &tmp);
        s << tmp.x();
        s << "\t";
        s << tmp.y();
        return(s.str());
        }

    void normalize (uint frames)
        {
        loos::GCoord val;
        loos::GCoord zero(0., 0., 0.);
        for (uint i=0; i< _xbins; ++i)
            {
            for (uint j=0; j< _ybins; ++j)
                {
                get(i, j, &val);
                uint norm = get_norm(i,j);
                if (norm)
                    {
                    set(i,j, val /get_norm(i,j));
                    }
                else
                    set(i,j, zero);
                }
            }
        }
};

//* Calculate dipole moment
class CalcDipole : public CalcProperty<loos::GCoord>
{
public:
    CalcDipole(uint xbins, uint ybins) :
        CalcProperty<loos::GCoord>(xbins, ybins)
        {
        }

    void calc(const loos::AtomicGroup &group, const uint xbin, const uint ybin)
        {
        incr(xbin, ybin, group.dipoleMoment());
        }

    const std::string print(const uint xbin, const uint ybin)
        {
        std::stringstream s;
        loos::GCoord tmp;
        get(xbin, ybin, &tmp);
        s << tmp.x();
        s << "\t";
        s << tmp.y();
        s << "\t";
        s << tmp.z();
        return(s.str());
        }

    void normalize (uint frames)
        {
        loos::GCoord val;
        loos::GCoord zero(0., 0., 0.);
        for (uint i=0; i< _xbins; ++i)
            {
            for (uint j=0; j< _ybins; ++j)
                {
                get(i, j, &val);
                uint norm = get_norm(i,j);
                if (norm)
                    {
                    set(i,j, val /get_norm(i,j));
                    }
                else
                    set(i,j, zero);
                }
            }
        }

};
