# REVACpp
Relevance Estimation and Value Calibration in C++

```c++

#include "Revac.h"
#include <chrono>	// time


// sphere fitness function
double sphere(const ParaVector& paraVector)
{
	double obj = 0;
	for (int i = 0; i < paraVector.k(); i++)
	{
		double xi = paraVector.unscaled_at(i);
		obj += std::pow(xi, 2);
	}
	return obj;
}

// rosenbrock fitness function
double rosenbrock(const ParaVector& paraVector)
{
	double obj = 0;
	for (int i = 0; i < paraVector.k() - 1; i++)
	{
		double xi = paraVector.unscaled_at(i);
		double xiplus = paraVector.unscaled_at(i + 1);

		obj += 100 * std::pow(xiplus - std::pow(xi, 2), 2) + std::pow(1 - xi, 2);
	}
	return obj;
}

int main()
{
	revacGen.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());

	int p = 50;   // number of parents
	int h = 5;	 // considered neighbors in mutation
	int k = 60;	 // number of parameter
	int s = 100; // pop size
	std::vector<Bound> bounds(k, Bound{ -10, 10 });

	ParaVector solution = REVAC(k, bounds.data(), sphere, 20000, p, h, s, std::cout);

	std::cout << solution;
}
