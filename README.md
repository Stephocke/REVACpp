# REVACpp
Relevance Estimation and Value Calibration in C++

```c++

#include "Revac.h"
#include <chrono>	// time


// sphere fitness function
double sphere(const ParaVector& paraVector)
{
	int i = 0;
	return std::accumulate(paraVector.data().begin(), paraVector.data().end(), 0.0,
		[&](double& sum, const double& val)
		{;
			double us = unscale0_1(val, paraVector.bounds()[i].min, paraVector.bounds()[i].max);
			i++;
			return sum + us * us;
		});
}

// rosenbrock fitness function
double rosenbrock(const ParaVector& paraVector)
{
	double obj = 0;
	for (int i = 0; i < paraVector.k() - 1; i++)
	{
		double xi = unscale0_1(paraVector.at(i), paraVector.bounds()[i].min, paraVector.bounds()[i].max);
		double xiplus = unscale0_1(paraVector.at(i + 1), paraVector.bounds()[i + 1].min, paraVector.bounds()[i + 1].max);

		obj += 100 * std::pow(xiplus - std::pow(xi, 2), 2) + std::pow(1 - xi, 2);
	}
	return obj;
}

int main()
{
	gen.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());

	int p = 50;   // number of parents
	int h = 5;	 // considered neighbors in mutation
	int k = 60;	 // number of parameter
	int s = 100; // pop size
	std::vector<Bound> bounds(k, Bound{ -10, 10 });

	ParaVector solution = REVAC(k, bounds.data(), rosenbrock, 20000, p, h, s, std::cout);

	std::cout << solution;
}
