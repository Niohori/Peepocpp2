#ifndef VISION__H
#define VISION__H
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>



static std::vector<double> rotate_point(std::vector<double> target, const double& angle)
{
	double s = sin(angle);
	double c = cos(angle);
	// rotate point
	double xnew = target[0] * c + target[1] * s;
	double ynew = -target[0] * s + target[1] * c;
	return { xnew,ynew };
}


static double collision(const std::vector<double>& observer, const std::vector<double>& target, const double& direction,
	const double& edge_1, const double& edge_2, const double& view_radius)
{
	//first calculate coordinates target and observerin rotated axis (as to direction as referential)

	std::vector<double> observer_ = rotate_point(observer, direction);
	std::vector<double> target_ = rotate_point(target, direction);
	//calculate angle vs the x-axis
	if (target_[0] == 0.0) { target_[0] += 0.00000001; }//safety in case atan infinite
	double angle = atan(target_[1] / target_[0]);
	double max_a = std::max(edge_1, edge_2);
	double min_a = std::min(edge_1, edge_2);
	if ((angle <= max_a) && (angle >= min_a))
	{
		double distance = sqrt(pow((observer_[1] - target_[1]), 2.) + pow((observer_[0] - target_[0]), 2.));
		if (distance <= view_radius) { distance; }
	}
	return -1.0;
}


static bool collision(const std::vector<double>& observer, const std::vector<double>& target, const double& dist)
{
	double distance = sqrt(pow((observer[1] - target[1]), 2.) + pow((observer[0] - target[0]), 2.));
	if (distance <= dist) { return true; }
	return false;
}


#endif
