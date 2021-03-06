#ifndef GENERATIVE_MODEL__H
#define GENERATIVE_MODE__H

#include "peepo_network.h"
#include "sensory_input.h"


class GenerativeModel {
private:
	
	

public:
	PeepoNetwork peepo_network;
	GenerativeModel();
	GenerativeModel(PeepoNetwork&, SensoryInput&);
	std::shared_ptr<SensoryInput> sensory_input;
	double process(void);
	std::map<std::string, std::vector<double> > predict(void);
	std::vector<double> error_(const std::vector<double>&, const std::vector<double>&);
	double error_size(const std::vector<double>&, const std::vector<double>&);
	double precision(const std::vector<double>&);
	void error_minimization(const std::string&, const double&, const std::vector<double>&, const std::vector<double>&);
	void hypothesis_update(const std::string&, const std::vector<double>&, const std::vector<double>&);
	std::map<std::string, unsigned >  get_root_values(void);
	std::vector<std::string> get_roots(void);
	std::vector<std::string >  get_leafs(void);
	//unsigned get_node_index(const std::string&);// NOT USED
	bool is_leaf(const std::string&);
	//bool is_root(const std::string&);//NOT USED
};

#endif

