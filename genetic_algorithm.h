#ifndef GENETICALGORITHM__H
#define GENETICALGORITHM__H


#include <dlib/matrix.h>

#include <boost/range/join.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/assign.hpp>
#include <boost/variant.hpp>
//#include <boost/filesystem.hpp>

#include <nlohmann/json.hpp>

#include <stdlib.h>    
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <time.h>
#include <chrono>
#include <fstream>
#include <tuple>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <random>
//#include <variant>//only in C++17-> use boost instead

#include "peepo_network.h"
#include "combinations.h"

using json = nlohmann::json;
using std::setw;
using std::setprecision;
using std::left;
using std::fixed;
using std::right;
using std::scientific;





using namespace boost::assign;
using namespace boost::adaptors;
//using namespace boost::filesystem;
using namespace boost;



struct	Agent
{
	Agent(void);
	Agent(const double&, const PeepoNetwork&, const double& = 0.0, const double& = 0.0);
	Agent& operator=(const Agent&);
	Agent(const Agent&);
	double fitness;
	PeepoNetwork network;
	double mut_top;
	double mut_cpd;
	bool operator<(const Agent& );
	bool operator>(const Agent&);
};

std::vector<json> get_topologies(PeepoNetwork&, const unsigned&);
matrix<unsigned> get_adjency_map(const std::vector<std::vector<std::string> >&,
	const std::vector<std::string>&, const std::vector<std::string>&);
std::vector<std::vector<std::string> > adjency_to_edges(const matrix<unsigned>&,
	const std::vector<std::string>&, const std::vector<std::string>&);
bool valid_graph(matrix<unsigned>&,
	const std::vector<std::string>&, const std::vector<std::string>&);
std::vector<std::vector<Agent>> get_couples(std::vector<Agent >);
bool is_in_edges(const std::vector<std::vector<std::string> >&, const std::string&, const std::string&);




class GeneticAlgorithm {

private://members
	std::ofstream log;// ("logfile.txt");//TEMPORARY FOR DEBUGGING PURPOSES
// ----------------------------------------------------------------------------------------

	long long seed;// = std::chrono::system_clock::now().time_since_epoch().count();//seed
	double POPULATION_CONVERGENCE_RATIO = 1. / 5.;
	int NUMBER_OF_MATING_PARENTS = 2;
	double PI = 3.141596;

	double population_fitness_score;

private://functions
	std::vector<double> generate_rand_vec(const unsigned&, const double&);
	std::vector<Agent> copy_agents(const std::vector<Agent>&);


public://members
	bool fast;
	int npop;
	std::vector<Agent> population;
	std::vector<Agent> selected_parents;
	int number_of_parents;
	double convergence_ratio;
	double p_mut_top;
	double p_mut_cpd;
	std::vector<std::string> root_nodes;
	std::vector<std::string> leaf_nodes;
	json cardinality_map;
	PeepoNetwork peepo;
	unsigned max_removal;
	unsigned treshold;//controls the allowed distance of the parents(i.e.how much edges they differ from each other)
	Agent best_chromosome;
	std::map<double, std::vector<Agent> > last_generation;
	unsigned convergence_period;
	std::vector<double> convergence_history;
	double convergence_sensitivity;


public://functions
	GeneticAlgorithm();
	GeneticAlgorithm(const std::string& ,
		const bool& ,
		const unsigned& ,
		const double& ,
		const unsigned&  ,
		const double& ,
		const double& ,
		const unsigned& );
	~GeneticAlgorithm();
	std::vector<Agent> get_population(void);
	void initialize(const std::string&);
	void get_parents(void);
	Agent get_optimal_newtork(void);
	bool check_convergence(void);
	std::vector<Agent> evolve(std::vector<Agent>, double&, bool&);
	std::vector<Agent> cross_over(void);
	PeepoNetwork prune_or_grow(PeepoNetwork);
	PeepoNetwork mutate_cpds(PeepoNetwork);
	t_cpd_leaf ga_child_cpd(const std::vector<unsigned>&, const std::vector<double>&);
	t_cpd_leaf normalize_distribution(t_cpd_leaf&);
};
#endif
