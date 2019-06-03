#pragma once//---------------------------------------------------------------------------

#include "model.h"
#include <tuple>
#include <fstream>

//---------------------------------------------------------------------------

struct Date
{
    int year, month, day;
    int hour, minute;
    long double seconds;
    bool operator <(const Date& x)
    {
        return std::tie(year, month, day, hour, minute, seconds) < std::tie(x.year, x.month, x.day, x.hour, x.minute, x.seconds);
    }
};

//---------------------------------------------------------------------------
class EarthSolarRotation : public TModel
//Earth rotation around Solar, rough calculation
{
protected:
    int sun_time{0}, days{0}, sun_time_8_20{0};
    static const long double startConditions[5][6];
    const long double mu = 132712.43994e+6, //km^3/c^2
                      J2000 = 2451544.5,   //Julian days
                      omega = 7.292115E-5, //rad/c
                      Re = 6371.3; //km
    long double gnomonHeight = 0.0, ro = 1.0, S = 0.0, Srad = 0.0,  Sg0 = 0.0, Tstar = 0.0, shadowMaxLen = 10.0; //Gnomon height, distanse to Solar, Current Star Time Place calculation, CSTP radians
    long double longtitude = 0.0/*0.645771823236*/, latitude = 0.0;//3.14159/2./*0.959931088595*/; //Gnomon position
    TVector Rsh, Rg, r0, ReStar, Re0, r;//shadow, gnomon, from Earth center to Earth surface, shadow from top of gnomon to surface, from Solar to Earth normalized
    TMatrix A;
    Date checkDay, startDates[5];
    Date d0 ;//start Dates
    bool daytime = false, night = false; //is it day or night
    std::ofstream file;
    std::ofstream fileLLS;
    std::ofstream sun_time_file;


    //Date start, finish, ephemerids;
public:
    EarthSolarRotation();
    EarthSolarRotation( long double t0, long double tk, TVector& V );
    EarthSolarRotation( Date dk, long double latitude, long double longtitude, bool flag);//, int date );
    long double checkLen(long double len);
    long double toJulianDate( Date date );
    void getMatrixA(long double phi, long double s);
    void getRight( const TVector& X, long double t, TVector& Y );
    void do_thing( const TVector& X, long double t );
    void finish();
    void sunlight();
    long double t_temp;
    long double d;
    int flag_l;
    bool once;
};
