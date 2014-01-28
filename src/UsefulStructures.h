/*
 * UsefulStructures.h
 *
 *  Created on: 18 бер. 2013
 *      Author: kopp
 */

#ifndef USEFULSTRUCTURES_H_
#define USEFULSTRUCTURES_H_
#include <vector>

struct Range
{
	double m_min;
	double m_max;
	size_t m_sampling;

	double getStep() const
	{
		return (m_sampling > 1) ? (m_max - m_min) / (m_sampling - 1) : 0.0;
	}
	void toVector(std::vector<double>& vec, double scal = 1.0) const
	{
		double step = getStep();
		for (size_t i = 0; i < m_sampling; ++i)
		{
			vec.push_back(scal * (m_min + step * i));
		}
	}
	bool good() const
	{
		if ((m_min <= m_max) && (m_sampling > 0))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

std::ostream& operator<<(std::ostream& out, const Range& range);

#endif /* USEFULSTRUCTURES_H_ */
