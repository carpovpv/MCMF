#include <stdio.h>


int main()
{
	double x, y, ot;
	double m1, m2, m_ot1, m_ot2;

	m1  = 10.0;
	m2 = -10.0;

	m_ot1 = m_ot2  = 0;


	while(scanf("%lf %lf %lf", &x, &y, &ot) == 3)
	{

		double r = y - (1.0 - x);
		if(r >= 0)
		{
			if( r <= m1)
			{
				m1 = r;
				m_ot1 = ot;
			}
		}
		else
		{
			if(r > m2)
			{
				m2 = r;
				m_ot2 = ot;
			}
		}

	}
	ot = (m_ot1 + m_ot2) / 2.0;
	printf("%.5f\n", m_ot2);

	return 0;
}
