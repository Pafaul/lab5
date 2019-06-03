//---------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include "custom.h"
#include <fstream>

#include "LA.h"

#define PI 3.14159265359

//---------------------------------------------------------------------------
const long double EarthSolarRotation::startConditions[5][6] = {
    {-2.566123740124270E7, 1.33935023154466E7, 5.805149372446711E7, -2.983549561177192E1, -4.846747552523134, -2.100585886567924},
    {-1.491959571635911E8, 1.457650675437370E6, 6.257355288057094E5, -5.907576993605990E-1, -2.743215365812734E1, -1.189141047828979E1},
    {-1.636143269850053E5, -1.3844477675798751E8, -6.001897492552029E7, 2.928899278768503E1, -7.429357814949981E-2, -3.187227439684210E-2},
    {1.469991745241464E8, -3.624128890102431E5, -1.574568651085924E5, -1.97050005690881E-1, 2.722474550315330E1, 1.180232151662914E1},
    {6.034600352945230E5, 1.360371485420887E8, 5.897533756750560E7, -3.028866166273934E1, 1.088711793301731E-1, 4.771257185837333E-2}
};

void printDate1(Date date)
{
    std::cout << "Chosen day: " << std::endl;
    std::cout << "Year: " << date.year << ", Month: " << date.month << ", Day: " << date.day << ", Time: 00:00:00" << std::endl;
}

EarthSolarRotation::EarthSolarRotation() : TModel()
{
    A.resize(3, 3);
    X0.resize(6);
    X0[0] = -2.566123740124270e+7L; //km
    X0[1] = 1.339350231544666e+8L;
    X0[2] = 5.805149372446711e+7L;
    X0[3] = -2.983549561177192e+1L; //km/c
    X0[4] = -4.846747552523134L;
    X0[5] = -2.100585886567924L;
}

long double EarthSolarRotation::checkLen(long double len)
{

    if (len > shadowMaxLen) return shadowMaxLen; else return len;
}

EarthSolarRotation::EarthSolarRotation( long double t0, long double tk, TVector& V ) : EarthSolarRotation()
{
    this->t0 = t0;
    this->t1 = tk;
    this->X0 = TVector(V);
}

EarthSolarRotation::EarthSolarRotation( Date dk, long double latitude, long double longtitude, bool flag)
{
    std::cout << "latitude : " << latitude << "; longtitude : " << longtitude << std::endl;
    long double temp;
    startDates[0].year = 2019;
    startDates[0].month = 1;
    startDates[0].day = 1;
    startDates[1].year = 2019;
    startDates[1].month = 3;
    startDates[1].day = 21;
    startDates[2].year = 2019;
    startDates[2].month = 6;
    startDates[2].day = 22;
    startDates[3].year = 2019;
    startDates[3].month = 9;
    startDates[3].day = 21;
    startDates[4].year = 2019;
    startDates[4].month = 12;
    startDates[4].day = 22;

    for (int i = 4; i >= 0; i--) {
        if ( startDates[i] < dk) {
            d0 = startDates[i];
            X0.resize(6);
            for (int j = 0; j < 6; j++)
                X0[j] = startConditions[i][j];
            break;
        }
    }
    d0 = startDates[0];
    checkDay = d0;
    dk.day += 1; dk.hour = 23; dk.minute = 59; dk.seconds = 59;
    Tstar = ((int) (toJulianDate(d0) - J2000) )/36525.;
    Sg0 = 24110.54841 + 8640184.812866*Tstar + 0.093104*pow(Tstar, 2) - 6.2*10E-6*pow(Tstar, 3); modf(Sg0, &temp);
    Sg0 = 2*M_PI/86400*(fmod(Sg0, 86400.)); std::cout << " Sg0 = " << Sg0 << std::endl;
    gnomonHeight = 1.0;
    Rsh = TVector(3);
    Rg = TVector(3);
    r0 = TVector(3);
    r = TVector(3);
    ReStar = TVector(3);
    Re0 = TVector(3);
    this->latitude = latitude; this->longtitude = longtitude;
    std::cout<<std::fixed;
    this->t0 = toJulianDate(d0)*86400.; std::cout << "T0: " << this->t0 << std::endl; printDate1(d0);
    if (flag == false)
       {
            this->t1 = toJulianDate(d0)*86400.; std::cout << "T1: " << this->t1 << std::endl; printDate1(dk);
       }
    else
        {
            this->t1 = toJulianDate(dk)*86400; std::cout << "T1: " << this->t1 << std::endl; printDate1(dk);
        }
    std::cout<<"s perehodom na letnee? "<<std::endl;
    std::cin>>flag_l;
    once = false;

    file.open("output.txt", std::ios_base::out);
    fileLLS.open("testLLS.txt", std::ios_base::out);
    sun_time_file.open("sun_time_file.txt");
    std::cout<<std::fixed;
    std::cout<<(t1-t0)/(3600.*24.)<<std::endl;
}

long double EarthSolarRotation::toJulianDate(Date date)
{
    int a = (14 - date.month)/12,
        M = date.month+12*a - 3,
        Y = date.year+4800 - a,
        JDN = date.day + ((153*M+2)/5) + 365*Y
            +(Y/4) - (Y/100) + (Y/400) - 32045;
    return JDN + (date.hour - 12)/24. + (date.minute)/1440. + date.seconds/86400.;
}

void EarthSolarRotation::getMatrixA(long double phi, long double s)
{
    A.resize(3, 3);
    A(0, 0) = - sin(phi) * cos(s);
    A(0, 1) = - sin(phi) * sin(s);
    A(0, 2) =   cos(phi);
    A(1, 0) =   cos(phi) * cos(s);
    A(1, 1) =   cos(phi) * sin(s);
    A(1, 2) = - sin(phi);
    A(2, 0) = - sin(s);
    A(2, 1) =   cos(s);
    A(2, 2) =   0;

}

void EarthSolarRotation::getRight( const TVector& X, long double t, TVector& Y )
{

    Y.resize(6);
    Y[0] = X[3];
    Y[1] = X[4];
    Y[2] = X[5];
    ro = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
    Y[3] = -mu*X[0]/pow(ro, 3.);
    Y[4] = -mu*X[1]/pow(ro, 3.);
    Y[5] = -mu*X[2]/pow(ro, 3.);
}

void EarthSolarRotation::do_thing(const TVector &X, long double t)
{
    if (once == false){
        //d = 1 + (toJulianDate(checkDay)*86400.)/(60.*60.*24) - (toJulianDate(d0)*86400.)/(60.*60.*24);
        d = 1 + (t)/(60.*60.*24) - (toJulianDate(d0)*86400.)/(60.*60.*24);
        once = true;
    }

    if(floor((t-t0)/(3600.*24.)) >= d)
    {
        if (toJulianDate(checkDay)*86400. < t)
        {
            S = Sg0 + omega*(t - t0) + longtitude;
            for(int i = 0; i < 3; i++) Re0[i] = X[i]/X.length();

            r[0] = Re*cos(latitude)*cos(S);
            r[1] = Re*cos(latitude)*sin(S);
            r[2] = Re*sin(latitude);

            r0 = r*(1/r.length());

            Rg = r0*gnomonHeight;
            ReStar = Re0*(-gnomonHeight/(Re0*r0));
            Rsh = Rg + ReStar;
            getMatrixA(latitude, S); Rsh = A * Rsh;

            long double alfa = acos(Re0*r/(Re0.length()*r.length()));

            if (alfa > M_PI/2.0) daytime = true; else daytime = false;

            if (daytime == true)
            {
                //если солнце над местным горизонтом - добавляем "светлую" минуту
                ++sun_time;
                t_temp = (t-toJulianDate(checkDay)*86400.)/(60.*60.) - (days-d)*24 + ceil(longtitude / 0.26);
                if ((days - d >= 90) && (days - d <= 300) && (flag_l == 1))
                {
                    t_temp += 1;
                }
                /*file << t_temp << " "
                 << checkLen(sqrt(pow(Rsh[0],2)+pow(Rsh[1],2))) << " "
                 << Rsh[0] << " " << Rsh[1] << " " << Rsh[2] << " "
                 << daytime << std::endl;*/
                //если время между 8 и 20 часами текущего дня, то добавляем светлую минуту в этом отрезке
                if ((t_temp > 8.0) && (t_temp < 20.0 ))
                {
                    ++sun_time_8_20;
                }
            }

            if (t > t0 + days * 86400.)
            {
                if(sun_time != 0)
                {
                    if(sun_time >= 720.0) sun_time = 720.0;
                    sun_time_file << days - d <<" "<<sun_time / 60.<<" "<<sun_time_8_20 / 60.<< std::endl;
                }
                ++days;
                sun_time = 0;
                sun_time_8_20 = 0;
            }
        }
    }
}

void EarthSolarRotation::sunlight()
{

}

void EarthSolarRotation::finish()
{
    file.close();
    fileLLS.close();
    sun_time_file.close();
}
