#include "Ez_point.h"

std::ostream& operator<< (std::ostream &out, const Point &point)
{
    out << point.m_Ez << "," 
        << point.m_Dz << "," 
        << point.m_Hx << "," 
        << point.m_Hy << std::endl; 
    return out;
}

Point operator+ (Point point1, Point point2)
{
    Point new_point;
    
    new_point.SetEz(point1.GetEz() + point2.GetEz());
    new_point.SetDz(point1.GetDz() + point2.GetDz());
    new_point.SetHx(point1.GetHx() + point2.GetHx());
    new_point.SetHy(point1.GetHy() + point2.GetHy());

    return new_point;
}

// Setters:
Point::Point(){}
Point::Point(double setter[4])
{
    m_Ez = setter[0];
    m_Dz = setter[1];
    m_Hx = setter[1];
    m_Hy = setter[3];
}


void Point::SetEz(double value){m_Ez=value;}
void Point::SetDz(double value){m_Dz=value;}
void Point::SetHx(double value){m_Hx=value;}
void Point::SetHy(double value){m_Hy=value;}

void Point::InjectEz(double value){m_Ez+=value;}
void Point::InjectDz(double value){m_Dz+=value;}
void Point::InjectHx(double value){m_Hx+=value;}
void Point::InjectHy(double value){m_Hy+=value;}

