#ifndef TQUATERNION_H
#define TQUATERNION_H
#include "LA.h"


class TQuaternion
{
protected:
long double q;
TVector Q;

public:
    TQuaternion();
    TQuaternion(long double l0, long double l1, long double l2, long double l3);
    TQuaternion(long double phi, const TVector& V);
    TQuaternion(const TQuaternion& rval);

    inline long double scal() const { return q; }
    inline TVector vect() const { return Q; };

    TQuaternion& operator = (const TQuaternion& rval);
    TQuaternion operator - (const TQuaternion& rval) const;
    TQuaternion operator + (const TQuaternion& rval) const;
    TQuaternion operator * (const TQuaternion& rval) const;
    TQuaternion operator * (const TVector& rval) const;
    TQuaternion operator ! () const;
    double& operator [] (const int i);

    TQuaternion& norm();

    TQuaternion conj() const;

    TMatrix rotateMatrix() const;

    static TQuaternion KrilAngles(long double phi, long double psi, long double theta);

};

#endif // TQUATERNION_H
