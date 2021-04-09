#ifndef REVAC_H
#define REVAC_H

#include <vector>	// container
#include <math.h>	// some math :)
#include <numeric>	// accumulate
#include <random>	// default_random_engine / distributions
#include <iostream> // cout


std::default_random_engine gen;

double scale0_1(double val, double min, double max)
{
	return (val - min) / (max - min);
}

double unscale0_1(double val, double min, double max)
{
	return val * (max - min) + min;
}

struct Bound
{
	double min, max;
};

class ParaVector;
using Population = std::vector<ParaVector>;


class ParaVector
{
public:
	ParaVector() :m_age(0), m_fitness(INFINITY) {}

	int k()const { return int(m_data.size()); }
	int age()const { return m_age; }
	const Bound*bounds()const { return m_bounds.data(); }
	double at(int i)const { return m_data[i]; }
	double fitness()const { return m_fitness; }
	const std::vector<double>& data()const { return m_data; }

	void clear();
	void init(int k, const Bound*bounds, double(*objective)(const ParaVector&) = nullptr);
	void updateAge() { m_age++; }

	friend ParaVector uscanning(Population&pop, int p, double(*objective)(const ParaVector&));
	friend void mutation(Population&pop, ParaVector& offspring, int p, int k, int h, double(*objective)(const ParaVector&));

	friend std::ostream& operator<<(std::ostream&stream, const ParaVector& vec);
private:
	int m_age;
	double m_fitness;
	std::vector<double> m_data;
	std::vector<Bound> m_bounds;
	void evaluate(double(*objective)(const ParaVector&)) { m_fitness = objective(*this); }
};



// computes the average fitness of a population vector
double avgFitness(const Population& pop)
{
	int n = pop.size();
	double sum = std::accumulate(pop.begin(), pop.end(), 0.0, [&](double val, const ParaVector& vec)
		{
			return val + vec.fitness();
		});
	return sum / n;
}


// outputs the unscaled values of the m_data vector
std::ostream& operator<<(std::ostream&stream, const ParaVector& vec)
{
	int i = 0;
	for (auto& el : vec.data())
	{
		stream << unscale0_1(el, vec.bounds()[i].min, vec.bounds()[i].max) << "\t";
		i++;
	}
	stream << "\n";

	return stream;
}


// uniform scanning crossover 
// pop - population vector 
// p - number of parents concerned
// objective - objective function to evaluate a ParaVector
ParaVector uscanning(Population&pop, int p, double(*objective)(const ParaVector&))
{
	int k = pop[0].k();
	std::vector<int> idx(p);
	std::iota(idx.begin(), idx.end(), 0);
	std::shuffle(idx.begin(), idx.end(), gen);

	int pidx = -1;
	int pctr = 0;
	int Ei = int(std::ceil(1.0 / p * k));
	ParaVector offspring;
	offspring.init(k, pop[0].bounds());
	for (int i = 0; i < k; i++)
	{
		if (i%Ei == 0 && pctr < p)
		{
			pidx = idx[pctr];
			pctr++;
		}
		offspring.m_data[i] = pop[pidx].at(i);
	}

	offspring.evaluate(objective);

	return offspring;
}

// REVAC mutation operator
// pop - population
// offspring - created from crossover
// p - number of parents
// k - column or parameter i respectively
// h - REVAC parameter determines the range of the mutation interval
void mutation(Population&pop, ParaVector& offspring, int p, int k, int h, double(*objective)(const ParaVector&))
{
	std::vector<double> D(p, 0);

	for (int i = 0; i < p; i++)
	{
		D[i] = pop.data()[i].at(k);
	}

	std::sort(D.begin(), D.end());
	int lb = std::distance(D.begin(), std::lower_bound(D.begin(), D.end(), offspring.data()[k]));
	int ub = std::distance(D.begin(), std::upper_bound(D.begin(), D.end(), offspring.data()[k]));

	double mintervalBegin, mintervalEnd;
	if (lb - (h - 1) < 0)
	{
		mintervalBegin = std::max(offspring.bounds()[k].min,
			D[0] - std::max(0.001, (D[p - 1] - D[0]) / p * ((h - 1) - lb)));
	}
	else
	{
		mintervalBegin = D[lb - (h - 1)];
	}
	if (ub + h - 1 >= p)
	{
		mintervalEnd = std::min(offspring.bounds()[k].max,
			D[p - 1] + std::max(0.001, (D[p - 1] - D[0]) / p * (ub + h - 1 - (p - 1))));
	}
	else
	{
		mintervalEnd = D[ub + h - 1];
	}

	offspring.m_data[k] = std::uniform_real_distribution<double>(mintervalBegin, mintervalEnd)(gen);
	offspring.evaluate(objective);
}


// updates the age of all ParaVectors within the population 
// returns the index of the element comprising the highest age
int updateAge(Population& pop)
{
	int amaxIdx = -1;
	int amax = 0;
	// starting at the end to avoid removing the best individual at the start of revac
	for (int i = pop.size() - 1; i >= 0; i--)
	{
		pop[i].updateAge();
		if (amax < pop[i].age())
		{
			amax = pop[i].age();
			amaxIdx = i;
		}
	}

	return amaxIdx;
}


// Relevance estimation and value calibration
// k - number of parameters
// bounds - array of bounds associated to the k variables
// objective - objective function to evaluate a ParaVector
// iterations - number of iterations 
// p - number of parents
// h - numbers of neighbors considered during mutation 
// s - population size
// stream - output stream, e.g., a file or the console
ParaVector REVAC(int k, const Bound*bounds, double(*objective)(const ParaVector&), int iterations = 10000, int p = 50, int h = -1, int s = 100,std::ostream& stream = std::cout)
{
	// initialization
	Population pop(s, ParaVector());
	double bestFitness = INFINITY;
	ParaVector bestConfig;
	for (int i = 0; i < s; i++)
	{
		pop[i].init(k, bounds, objective);
	}
	h = h == -1 ? p / 10 : h;

	// run until termination
	for (int it = 0; it < iterations; it++)
	{
		std::sort(pop.begin(), pop.end(), [](const ParaVector& a, const ParaVector& b)
			{
				return a.fitness() < b.fitness();
			});

		if (pop[0].fitness() < bestFitness)
		{
			bestFitness = pop[0].fitness();
			bestConfig = pop[0];
		}

		if (it % 100 == 0)
			stream << "Iteration " << it << ": Avg=" << avgFitness(pop)
			<< "\tBestPop=" << pop[0].fitness()
			<< "\tBestOpt=" << bestFitness << "\n";

		ParaVector offspring = uscanning(pop, p, objective);
		for (int i = 0; i < k; i++)
			mutation(pop, offspring, p, i, h, objective);

		int idxOldest = updateAge(pop);
		pop[idxOldest] = offspring;
	}

	return bestConfig;
}


// clears and initializes the m_data vector with random values drawn from a uniform distribution
// k - number of parameters
// bounds - array comprising the bounds associated with the k variables
// objective - objective function to evaluate the initialized paraVector - if null fitness is INFINITY
void ParaVector::init(int k, const Bound*bounds, double(*objective)(const ParaVector&))
{
	clear();
	m_bounds.reserve(k);
	m_data.reserve(k);
	for (int i = 0; i < k; i++)
	{
		m_bounds.push_back(bounds[i]);
		double val = std::uniform_real_distribution<double>(bounds[i].min, bounds[i].max)(gen);
		m_data.push_back(scale0_1(val, bounds[i].min, bounds[i].max));
	}

	if (objective != nullptr)
	{
		evaluate(objective);
	}
}

// clears the m_data vector
void ParaVector::clear()
{
	m_fitness = INFINITY;
	m_age = 0;
	m_data.clear();
	m_bounds.clear();
}

#endif // !REVAC_H
