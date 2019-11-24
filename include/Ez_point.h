#ifndef EZPOINT
#define EZPOINT
#include <iostream>


class Point
{
private:
    double m_Ez{0};
    double m_Dz{0};
    double m_Hx{0};
    double m_Hy{0};


public:
    Point();
    Point(double setter[4]);
    
    void SetEz(double value);
    void SetDz(double value);
    void SetHx(double value);
    void SetHy(double value);

    void InjectEz(double value);
    void InjectDz(double value);
    void InjectHx(double value);
    void InjectHy(double value);
    
    friend std::ostream& operator<< (std::ostream &out, const Point &point);
    
    double GetEz(){return m_Ez;}
    double GetDz(){return m_Dz;}
    double GetHx(){return m_Hx;}
    double GetHy(){return m_Hy;}

};

Point operator+ (Point point1, Point point2); 
#endif
