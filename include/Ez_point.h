#ifndef EZPOINT
#define EZPOINT
#include <iostream>


class YeePoint
{
private:
    double m_Ez{0};
    double m_Hx{0};
    double m_Hy{0};


public:
    YeePoint();
    YeePoint(double setter[3]);
    
    void SetEz(double value);
    void SetHx(double value);
    void SetHy(double value);

    
    friend std::ostream& operator<< (std::ostream &out, const YeePoint &point);
    
    double GetEz(){return m_Ez;}
    double GetHx(){return m_Hx;}
    double GetHy(){return m_Hy;}

};

YeePoint operator+ (YeePoint point1, YeePoint point2); 
#endif
