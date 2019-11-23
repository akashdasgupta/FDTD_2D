#include "Ez_point.h"

std::ostream& operator<< (std::ostream &out, const YeePoint &point)
{
    out << point.m_Ez << "," 
        << point.m_Hx << "," 
        << point.m_Hy << std::endl; 
    return out;
}

YeePoint operator+ (YeePoint point1, YeePoint point2)
{
    YeePoint new_point;
    
    new_point.SetEz(point1.GetEz() + point2.GetEz());
    new_point.SetHx(point1.GetHx() + point2.GetHx());
    new_point.SetHy(point1.GetHy() + point2.GetHy());

    return new_point;
}

// Setters:
YeePoint::YeePoint(){}
YeePoint::YeePoint(double setter[3])
{
    m_Ez = setter[0];
    m_Hx = setter[1];
    m_Hy = setter[2];
}


void YeePoint::SetEz(double value){m_Ez=value;}
void YeePoint::SetHx(double value){m_Hx=value;}
void YeePoint::SetHy(double value){m_Hy=value;}

