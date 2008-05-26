#include <iostream>
#include <string>
#include <vector>

#include <TimeSeries.hpp>

using namespace std;

int main(int argc, char *argv[]) {
    vector<double> tmp;
    for (int i=0; i<20; i++) {
        tmp.push_back(i);
    }

    TimeSeries foo;
    cout << foo.size() << endl;

    TimeSeries foo2(tmp);
    TimeSeries foo3(foo2);

    for (int i=0; i<20; i++) {
        cout << i << "  " << foo2[i] << "  " << foo3[i] << endl;
    }

    foo2 *= 5.0;

    TimeSeries foo4 = foo2 - foo3;
    TimeSeries foo5 = 1.0 + 2.0*foo3 - foo4;

    cout << endl;
    for (int i=0; i<20; i++) {
        cout << i << "  " << foo2[i] << "  " << foo4[i] << "  "
             << foo5[i] << endl;
    }

    cout << foo5.average() << "  "
         << foo5.variance() << "  "
         << foo5.stdev() << "  "
         << foo5.sterr() << "  "
         << endl;

}
