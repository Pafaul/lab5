/* Standard library	*/


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "iomanip"

#include "LA.h"

#define M_PI   3.14159265358979323846264338327950288
#define degToRad(angleInDegrees) ((angleInDegrees) * M_PI / 180.0)
#define radToDeg(angleInRadians) ((angleInRadians) * 180.0 / M_PI)

#include "integrator.h"
//#include "custom.h"
#include "custom.h"

void calcModel(TModel* model);
void printDate(Date date);
long double getString( std::string s );
Date createDate(void);

int main()
{
    bool flag = false;
    std::cout<<"daily evolution of gnomon shadow or early evolution of duration of daylight?(0|1)"<<std::endl;
    std::cin >> flag;
    Date startDate = createDate();
    TModel* model = new EarthSolarRotation( startDate, degToRad(getString("↕ (latitude)")), degToRad(getString("↔ (longtitde)")), flag);//, dateChoice());
    TIntegrator* Integrator = new TDormandPrinceIntegrator();
    Integrator->setPrecision(1E-16);
    Integrator->Run( model );
    calcModel(model);
    model->finish();
    delete model;
    delete Integrator;
	return 0;
}


void calcModel(TModel* model){

    std::ofstream file("Integr_res.txt");

        TMatrix Result = model->getResult();

        for (int i=0; i<Result.rowCount(); i++)
        {
            for (int j=0; j<Result.colCount(); j++)
            {
                file<<std::fixed;
                file << Result(i, j) << " ";
                std::cout<<std::fixed;
                //cout<<Result(i, j) << " ";
            }
            file << std::endl;
        }

        file.close();
}

void printDate(Date date)
{
    std::cout << "Year: " << date.year << ", Month: " << date.month << ", Day: " << date.day << ", Time: 00:00:00" << std::endl;
}

long double getString( std::string s )
{
    long double la = 0.0;
    std::cout << "Input " << s << " : ______\b\b\b\b\b\b"; std::cin >> la;
    return la;
}

Date createDate(void)
{
    Date date;
    std::cout << "Create date." << std::endl;
    std::cout << "Input year: "; std::cin >> date.year;
    std::cout << "Input month: "; std::cin >> date.month;
    std::cout << "Input day: "; std::cin >> date.day;
    std::cout << "Created date: " << std::endl; printDate(date);
    date.hour = 0; date.minute = 0; date.seconds = 0;
    return date;
}
