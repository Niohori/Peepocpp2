#include "generative_model.h"


GenerativeModel::GenerativeModel() {;} 


GenerativeModel::GenerativeModel(PeepoNetwork& a_peepo, SensoryInput& a_sensory):
	peepo_network(a_peepo),
	sensory_input(a_sensory.clone())
{
	bool ok = peepo_network.to_bayesian_network();
};


double GenerativeModel::process(void)
{
	double total_prediction_error_size = 0.0;
	std::map<std::string, std::vector<double> > predic = predict();
	for (auto node:predic) {
		std::string node_name = node.first;
		if (is_leaf(node_name)) {
			std::vector<double> prediction = node.second;
			std::vector<double> observation = sensory_input->value(node_name);
			std::vector<double> prediction_error = error_(prediction, observation);
			double prediction_error_size = error_size(prediction, observation);
			double precision_ = precision(prediction);
			total_prediction_error_size += prediction_error_size;
			if (prediction_error_size > 0.1) {
				error_minimization(node_name, precision_, prediction_error, prediction);
			}
		}
	}
	return total_prediction_error_size;
};

std::map<std::string, std::vector<double> > GenerativeModel::predict(void)
{
	std::map<std::string, unsigned>   evidences = get_root_values();
	return peepo_network.get_inference(evidences);
}

std::vector<double> GenerativeModel::error_(const std::vector<double>& pred, const std::vector<double>& obs)
{
	std::vector<double> err;
	if (pred.size() != obs.size()) {
		std::cout << "Error in error estimation (prediction and observation have different size." << std::endl;
		std::cin.get();
	}
	for (int i = 0; i < pred.size(); i++) {
		err.push_back(obs[i] - pred[i]);
	}
	return err;
}


double GenerativeModel::error_size(const std::vector<double>& pred, const std::vector<double>& obs)
{
	if (pred.size() != obs.size()) {
		std::cout << "Error in error size estimation (prediction and observation have different size." << std::endl;
		std::cin.get();
	}
	double entropy = 0.f;
	for (int i = 0; i < pred.size(); i++) {
		entropy += obs[i] * log(pred[i] / obs[i]);
	}
	return entropy;
}

double GenerativeModel::precision(const std::vector<double>& pred)
{
	double entropy = 0.f;
	for (int i = 0; i < pred.size(); i++) {
		entropy -=pred[i] * log2f(pred[i] );
	}
	return entropy;
}

void GenerativeModel::error_minimization(const std::string& node_name, const double& precision, 
										const std::vector<double>& prediction_error, 
										const std::vector<double>& prediction)
{

	hypothesis_update(node_name,  prediction_error, prediction);
	//?? why both functions as error_minimization doesn't preform anything
}

void GenerativeModel::hypothesis_update(const std::string& node_name, 
										const std::vector<double>& prediction_error,
										const std::vector<double>& prediction)
{
	bool is_action = false;
	for(auto node: peepo_network.get_pro_nodes()){
		if (node == node_name) { is_action = true; break; }
	}
	if (is_action) {
		sensory_input->action(node_name, prediction);
	}
	else{
		std::map<std::string, unsigned>   evidence;
		std::vector<double> sum = prediction_error;
		std::transform(sum.begin(), sum.end(), prediction.begin(), sum.begin(), std::plus<double>());
		evidence[node_name] = std::max_element(sum.begin(), sum.end()) - sum.begin();
		std::map<std::string, std::vector<double> > inference = peepo_network.get_inference(evidence);
		for (auto root : get_roots()) {
			//old_hypo = self.bayesian_network.states[root_index].distribution.items()--> ??
			std::vector<double>  new_hypo = inference.at(root);
			peepo_network.add_cpds(root, new_hypo);
		}
	bool ok = peepo_network.to_bayesian_network();//we need to refresh the bayesian newtork with the new cpds of the root
	}
}


std::map<std::string, unsigned > GenerativeModel::get_root_values(void)
{
	std::map<std::string, unsigned > roots;
	for (auto node : get_roots())
	{
		std::vector<double> cpd = peepo_network.get_cpds_root(node);
		roots[node] = std::max_element(cpd.begin(), cpd.end()) - cpd.begin();
	}
	return roots;
}


std::vector<std::string > GenerativeModel::get_roots(void)
{
	return peepo_network.get_active_root_nodes();
}


std::vector<std::string > GenerativeModel::get_leafs(void)
{
	return peepo_network.get_leaf_nodes();
}



bool GenerativeModel::is_leaf(const std::string& node)
{
	for (auto leaf : get_leafs()) {
		if (leaf == node) { return true; }
	}
	return false;
}

/*
bool GenerativeModel::is_root(const std::string& node)
{
	for (auto root : get_roots()) {
		if (root == node) { return true; }
	}
	return false;
}*/


//unsigned GenerativeModel::get_node_index(const std::string&) { ; }