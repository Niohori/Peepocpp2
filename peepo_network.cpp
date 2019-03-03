#include "peepo_network.h"




/*std::vector<json> get_topologies(PeepoNetwork& peepo_network, unsigned max_removal)
{
	std::vector<std::vector<std::string>  > max_edges = fully_connected_network(peepo_network);
	if (max_removal >= max_edges.size()) { max_removal = max_edges.size() - 1; }
	std::vector<json> topologies;
	json a;
	for (unsigned comb = max_edges.size(); comb > max_removal; comb--) {
		Combinations<std::vector<std::string> > comb_obj(max_edges, comb);
		std::vector<std::vector<std::vector<std::string>>> combinations = comb_obj.get();
		for (auto edges : combinations) {
			a["egdes"] = edges;
			a["entropy"] = edges.size();
			topologies.push_back(a);
		}
	}
	return topologies;
}

std::vector<std::vector<std::string>  >   fully_connected_network(PeepoNetwork& peepo_network)
{
	std::vector<std::string> root_nodes = peepo_network.get_root_nodes();
	std::vector<std::string> leaf_nodes = peepo_network.get_leaf_nodes();
	peepo_network.edges.clear();
	for(unsigned root = 0; root < root_nodes.size(); root++){
		for (unsigned leaf = 0; leaf < leaf_nodes.size(); leaf++) {
			peepo_network.add_edge({ root_nodes[root],leaf_nodes[leaf] });
		}
	}
	std::vector<std::vector<std::string> > a = peepo_network.edges;
	return a;
}*/

void write_to_file(const std::string& name, const json&  peepo_network) {
	std::string my_path = std::experimental::filesystem::current_path().u8string()+ "\\resources";
	if (!std::experimental::filesystem::exists(my_path)) {std::experimental::filesystem::create_directory(my_path);};
	my_path += "\\" + name + ".json";
	std::ofstream o(my_path);
	o << std::setw(4) << std::boolalpha << peepo_network << std::endl;
	o.close();
}


PeepoNetwork::PeepoNetwork()
{
	//ctor;
	 /* initialize random seed: */
	srand(time(NULL));
	
}


PeepoNetwork::PeepoNetwork(const std::string& source):identification(source + ".json")
{
	//identification = source + ".json";
	//log =  std::ofstream("logfile.txt");
	from_json();
	/* initialize random seed: */
	srand(time(NULL));
}


PeepoNetwork:: ~PeepoNetwork()
{
 //destructor;
}


/*
PeepoNetwork& PeepoNetwork::operator=(const PeepoNetwork& apeepo)
{
	
	//*this = PeepoNetwork(apeepo);
	return *this;
}*/


PeepoNetwork::PeepoNetwork (const PeepoNetwork& apeepo)
{
	identification = apeepo.identification;
	description = apeepo.description;
	train_from = apeepo.train_from;
	train_data = apeepo.train_data;
	frozen = apeepo.frozen;
	bel_nodes = apeepo.bel_nodes;
	mem_nodes = apeepo.mem_nodes;
	lan_nodes = apeepo.lan_nodes;
	ext_nodes = apeepo.ext_nodes;
	pro_nodes = apeepo.pro_nodes;
	int_nodes = apeepo.int_nodes;
	edges = apeepo.edges;
	cpds = apeepo.cpds;
	assemble();
	omega_map = apeepo.omega_map;
	//directed_graph<bayes_node>::kernel_1a_c* bn = new directed_graph<bayes_node>::kernel_1a_c;
	//solution.reset();
	//bn = directed_graph<bayes_node>::kernel_1a_c(apeepo.bn);
	//directed_graph<bayes_node>::kernel_1a_c bn(apeepo.bn);

}

void PeepoNetwork::from_json(void)
{
		// read the JSON file
		std::ifstream ntw_(identification);
		//stream the file to the network object
		json ntw;
		ntw_ >> ntw;
		from_json(ntw);
}

void PeepoNetwork::from_json(const json& nw)
{
	//update class members
	if (identification == "") { identification = nw["header"]["identification"].get<std::string>(); }
	description = nw["header"]["description"].get<std::string>();
	train_from = nw["header"]["train_from"].get<std::string>();
	train_data = nw["header"]["train_from"].get<std::string>();
	frozen = nw["header"]["frozen"].get<bool>();
	date = nw["header"]["date"].get<std::string>();
	bel_nodes = nw["nodes"]["RON"]["BEL"].get< std::vector<json> >();
	mem_nodes = nw["nodes"]["RON"]["MEM"].get< std::vector<json> >();
	lan_nodes = nw["nodes"]["LAN"].get< std::vector<json> >();
	ext_nodes = nw["nodes"]["LEN"]["EXT"].get< std::vector<json> >();
	int_nodes = nw["nodes"]["LEN"]["INT"].get< std::vector<json> >();
	pro_nodes = nw["nodes"]["LEN"]["PRO"].get< std::vector<json> >();
	if (frozen) {
		edges = nw["edges"].get< std::vector<std::vector<std::string> > >();
		cpds = nw["cpds"].get<json>();
	}
	assemble();
}


void PeepoNetwork::from_json(const std::string& source)
{
	identification = source + ".json";
	from_json();
}

json PeepoNetwork::to_json(void)
{
	return network;
}


json PeepoNetwork::make_network(void)
{
	json nw;
	nw["header"]["identification"] = identification;
	nw["header"]["description"] = description;
	nw["header"]["train_from"] = train_from;
	nw["header"]["train_from"] = train_data;
	nw["header"]["frozen"] =  frozen ;
	nw["header"]["date"] = date; 
	nw["nodes"]["RON"]["BEL"] = bel_nodes;
	nw["nodes"]["RON"]["MEM"] = mem_nodes;
	nw["nodes"]["LAN"] = lan_nodes;
	nw["nodes"]["LEN"]["EXT"] = ext_nodes;
	nw["nodes"]["LEN"]["INT"] = int_nodes;
	nw["nodes"]["LEN"]["PRO"] = pro_nodes;
	nw["edges"] = edges;
	nw["cpds"] = cpds;
	return nw;
	}


void  PeepoNetwork::assemble(void)
{
	network.clear();
	network = make_network();
	make_cardinality_map();
	make_omega_map();
}


void  PeepoNetwork::disassemble(void)
{
	//network.clear();
	edges.clear();
	cpds.clear();
	omega_map.clear();
}


void PeepoNetwork::make_cardinality_map(void)
{
	auto range1 = boost::join(boost::join(bel_nodes, mem_nodes), lan_nodes);
	auto range2 = boost::join(range1, ext_nodes);
	auto range3 = boost::join(range2, int_nodes);
	auto range  = boost::join(range3, pro_nodes);
	cardinality_map.clear();
	for (auto it = boost::begin(range); it != boost::end(range); ++it){
		auto node = *it;
		cardinality_map.emplace(node["name"], node["card"]);
	}
}

json PeepoNetwork::get_cardinality_map(void)
{
	return cardinality_map;
}

void PeepoNetwork::make_omega_map(void)//TO DO : improve by use of existing cardinality_map
{
	omega_map.clear();
	auto range1 = boost::join(boost::join(bel_nodes, mem_nodes), lan_nodes);
	auto range2 = boost::join(range1, ext_nodes);
	auto range3 = boost::join(range2, int_nodes);
	auto range  = boost::join(range3, pro_nodes);
	for (auto it = boost::begin(range); it != boost::end(range); ++it) {
		auto node = *it;
		std::vector<std::string> incoming_nodes = get_incoming_edges(node["name"]);
		double max_omega = 2.0f*3.141596f;
		for (auto itt = boost::begin(incoming_nodes); itt != boost::end(incoming_nodes); ++itt) {
			max_omega *= cardinality_map[*itt];
		}
		std::vector<double> omega(node["card"],0.0 );
		for (unsigned k = 0; k < omega.size(); k++) { omega[k] =  static_cast <double> (std::rand()) / static_cast <double> (RAND_MAX);; }
		omega_map.emplace(node["name"], omega);
	}

}

json PeepoNetwork::get_omega_map(void)
{
	return omega_map;
}

void PeepoNetwork::add_omega_map(const std::string& node, const std::vector<double>& omega)
{
	omega_map.update({ { node,omega } });
}


std::vector<std::vector<std::string> >  PeepoNetwork::get_edges(void)
{ 
	return edges;
}



std::vector<std::string>  PeepoNetwork::get_nodes(void)
{
	std::vector<std::string> nodes;
	for (auto& node : cardinality_map.items()){nodes.push_back(node.key());}
	std::sort(nodes.begin(), nodes.end());
	return nodes;
}


std::vector<std::string> PeepoNetwork::get_root_nodes(void)
{
	std::vector<std::string> root_nodes;
	auto range = boost::join(bel_nodes, mem_nodes);
	for_each(range.begin(), range.end(), [&root_nodes](json node) {root_nodes.push_back(node["name"]); });
	std::sort(root_nodes.begin(), root_nodes.end());
	return root_nodes;
}


std::vector<std::string> PeepoNetwork::get_active_root_nodes(void)
{
	std::vector<std::string> active_root_nodes;
	for (auto leaf : get_leaf_nodes())
	{
		for (auto root : get_incoming_edges(leaf)) {
			active_root_nodes.push_back(root);
		}
	}
	std::sort(active_root_nodes.begin(), active_root_nodes.end());
	active_root_nodes.erase(std::unique(active_root_nodes.begin(), active_root_nodes.end()), active_root_nodes.end());
	return active_root_nodes;
}

std::vector<std::string> PeepoNetwork::get_active_nodes(void)
{
	std::vector<std::string> active_nodes = get_active_root_nodes();
	std::vector<std::string> leaf_nodes = get_leaf_nodes();
	active_nodes.insert(std::end(active_nodes), std::begin(leaf_nodes), std::end(leaf_nodes));
	return active_nodes;
}

std::vector<std::string> PeepoNetwork::get_leaf_nodes(void)
{
	std::vector<std::string> leaf_nodes;
	auto range2 = boost::join(pro_nodes, ext_nodes);
	auto range = boost::join(range2, int_nodes);
	for_each(range.begin(), range.end(), [&leaf_nodes](json node) {leaf_nodes.push_back(node["name"]); });
	std::sort(leaf_nodes.begin(), leaf_nodes.end());
	return leaf_nodes;
}


std::vector<std::string> PeepoNetwork::get_lan_nodes(void)
{
	std::vector<std::string> lan_;
	for_each(lan_nodes.begin(), lan_nodes.end(), [&lan_](json node) {lan_.push_back(node["name"]); });
	std::sort(lan_.begin(), lan_.end());
	return lan_;
}


std::vector<std::string> PeepoNetwork::get_bel_nodes(void)
{
	std::vector<std::string> bel_;
	for_each(bel_nodes.begin(), bel_nodes.end(), [&bel_](json node) {bel_.push_back(node["name"]); });
	std::sort(bel_.begin(),bel_.end());
	return bel_;
}


std::vector<std::string> PeepoNetwork::get_mem_nodes(void)
{
	std::vector<std::string> mem_;
	for_each(mem_nodes.begin(), mem_nodes.end(), [&mem_](json node) {mem_.push_back(node["name"]); });
	std::sort(mem_.begin(), mem_.end());
	return mem_;
}


std::vector<std::string> PeepoNetwork::get_ext_nodes(void)
{
	std::vector<std::string> ext_;
	for_each(ext_nodes.begin(), ext_nodes.end(), [&ext_](json node) {ext_.push_back(node["name"]); });
	std::sort(ext_.begin(), ext_.end());
	return ext_;
}


std::vector<std::string> PeepoNetwork::get_int_nodes(void)
{
	std::vector<std::string> int_;
	for_each(int_nodes.begin(), int_nodes.end(), [&int_](json node) {int_.push_back(node["name"]); });
	std::sort(int_.begin(), int_.end());
	return int_;
}


std::vector<std::string> PeepoNetwork::get_pro_nodes(void)
{
	std::vector<std::string> pro_;
	for_each(pro_nodes.begin(), pro_nodes.end(), [&pro_](json node) {pro_.push_back(node["name"]); });
	std::sort(pro_.begin(), pro_.end());
	return pro_;
}


void PeepoNetwork::add_belief_node(const std::string& node, const unsigned& cardinality)
{
	json a;
	a["name"] = node;
	a["card"] = cardinality;
	bel_nodes.push_back(a);
	cardinality_map[node] = cardinality;
}

void PeepoNetwork::add_memory_node(const std::string& node, const unsigned& cardinality)
{
	json a;
	a["name"] = node;
	a["card"] = cardinality;
	mem_nodes.push_back(a);
	cardinality_map[node] = cardinality;
}


void PeepoNetwork::remove_belief_node(const std::string& node)
{
	std::vector<json> belief_nodes;
	for (auto noode : bel_nodes) { if (noode["name"] != node) { belief_nodes.push_back(noode); } };
	bel_nodes = belief_nodes;
	remove_node(node);
}


void PeepoNetwork::remove_node(const std::string& node)
{
	cpds.erase(node);
	cardinality_map.erase(node);
	omega_map.erase(node);
	for (auto parent : get_incoming_edges(node)) { remove_edge({ parent,node }); };
	for (auto child : get_outgoing_edges(node)) { remove_edge({ node,child}); };
}


std::vector<std::string>  PeepoNetwork::get_incoming_edges( std::string node)
{
	std::vector<std::string> incoming_edges;
	for_each(edges.begin(), edges.end(), [node, &incoming_edges](std::vector<std::string> edge) {if (edge[1] == node) {incoming_edges.push_back(edge[0]); } });
	std::sort(incoming_edges.begin(), incoming_edges.end());
	return incoming_edges;
}


std::vector<std::string>  PeepoNetwork::get_outgoing_edges(const std::string& node)
{
	std::vector<std::string>outcoming_edges;
	for_each(edges.begin(), edges.end(), [node, &outcoming_edges](std::vector<std::string> edge) {if (edge[0] == node) { outcoming_edges.push_back(edge[1]); } });
	std::sort(outcoming_edges.begin(), outcoming_edges.end());
	return outcoming_edges;
}


void PeepoNetwork::set_edges(const std::vector<std::vector<std::string> >& edges_)
{
	edges = edges_;
}


void PeepoNetwork::add_edge(const std::vector<std::string>& edge)
{
	edges.push_back(edge);
}


void PeepoNetwork::remove_edge(const std::vector<std::string>& edge)
{
	edges.erase(std::remove_if(edges.begin(), edges.end(), [&edge](std::vector<std::string> edg) {return  edg[0] == edge[0] && edg[1] == edge[1]; }), edges.end());
}


void PeepoNetwork::set_cpds(const json& cpds_)
{
	cpds = cpds_;
}


json PeepoNetwork::get_cpds()
{
	return cpds;
}


t_cpd PeepoNetwork::get_cpds(const std::string& node)
{
	std::vector<std::string> root_nodes = get_root_nodes();
	bool is_root = false;
	for_each(root_nodes.begin(), root_nodes.end(), [&is_root, node](std::string node_) { if (node == node_) { is_root = true; }});
	if (is_root) {return cpds[node].get <t_cpd_root>(); }
	return cpds[node].get <t_cpd_leaf>();
}


t_cpd_root PeepoNetwork::get_cpds_root(const std::string& node)
{
	return cpds[node].get <t_cpd_root>();

}



void  PeepoNetwork::add_cpds(const std::string& node, const t_cpd_leaf& cpd)
{
	cpds.update({ { node,cpd } });
}


void  PeepoNetwork::add_cpds(const std::string node, const t_cpd_root cpd)
{
	cpds.update({ { node,cpd } });
}



bool PeepoNetwork::to_bayesian_network(void)
{
	assemble();
	/*   TO DELETE BECAME OBSOLENT AFTER ERROR FIXING IN GENETIC ALGO
	//clean up edges-> this should not be necessary -> flaw to find in GeneticAlgorithm->probably mutation yoels redundants egdes?
	std::sort(edges.begin(), edges.end());
	std::vector<std::vector<std::string>> edges_ = edges;
	for (size_t i = 0; i < edges_.size(); ++i)
	{
		edges_.erase(std::remove_if(edges_.begin() + i + 1, edges_.end(),
			[&edges_, &i](const std::vector<std::string>& _rhs)
		{
			return edges_[i] == _rhs;
		}),edges_.end());
	}
	edges = edges_;
	*/

	//create the necessary containers
	std::vector<std::string> root_nodes = get_active_root_nodes();
	std::vector<std::string> leaf_nodes = get_leaf_nodes();
	std::vector<std::string> all_nodes = get_active_nodes();
	unsigned number_of_nodes = all_nodes.size();
	unsigned number_of_edges = edges.size();
	

	// Create a dictionary of node names coupled to a number, as dLib respresent node by integers/unsigned
	for (unsigned node = 0; node < number_of_nodes; node++) { node_dic[all_nodes[node]] = node; };

	//get a dlib pointer
	bn.reset();
	bn = std::make_unique<directed_graph<bayes_node>::kernel_1a_c>();
	
	//inform how many the nodes  the network will have
	bn->set_number_of_nodes(number_of_nodes);

	//add the edges
	for (unsigned edge = 0; edge < number_of_edges; edge++) { bn->add_edge(node_dic[edges[edge][0]], node_dic[edges[edge][1]]); }
	//check if the graph is not broken in sevreal subgarphs
	if (!graph_is_connected(*bn)) {return false;}

	// Now we inform all the nodes in the network about their cardinality
	for (unsigned node = 0; node < number_of_nodes; node++) { set_node_num_values(*bn, node_dic[all_nodes[node]], cardinality_map[all_nodes[node]]); }

	//Add all the conditional probability information for each node
	// Each node's conditional probability is dependent on the state of its parents.  
	// To specify this state we need to use the assignment object.  This assignment 
	// object allows us to specify the state of each nodes parents. 
	// For parents nodes this object is empty use parents_state.clear() if previously assigned

	//assignment parent_state;
	assignment parent_state;

	//first the root nodes
	for (unsigned node = 0; node < root_nodes.size(); node++){
		unsigned cardinality = cardinality_map[root_nodes[node]];
		t_cpd_root probability = cpds[root_nodes[node]].get<t_cpd_root>();;
		parent_state.clear();
		for (unsigned state = 0; state < cardinality; state++) { set_node_probability(*bn, node_dic[root_nodes[node]], state, parent_state, probability[state]); }
	}

	//now the leaf nodes
	for (unsigned leaf = 0; leaf < leaf_nodes.size(); leaf++) {
		unsigned cardinality = cardinality_map[leaf_nodes[leaf]];
		t_cpd_leaf probability = cpds[leaf_nodes[leaf]].get<t_cpd_leaf>();;
		//get incoming nodes, their cardinality and make a matrix of posssible parents state combinations
		std::vector<std::string> in_nodes = get_incoming_edges(leaf_nodes[leaf]);
		std::vector<int> parents_states_cardinality;
		for (unsigned parent = 0; parent < in_nodes.size(); parent++) { parents_states_cardinality.push_back(cardinality_map[in_nodes[parent]]); }
		std::vector<std::vector<int> > parents_states_matrix = States::get_index_matrix(parents_states_cardinality);
		parent_state.clear();
		std::vector<unsigned> local_dict;//used as a conversion map for that specific leaf node
		for (auto parent : in_nodes) { 
			local_dict.push_back(node_dic[parent]); 
			parent_state.add(node_dic[parent], 0); 
		}
		for (unsigned states = 0; states < parents_states_matrix[0].size(); states++) {
			for (unsigned parent = 0; parent < parents_states_matrix.size(); parent++) { parent_state[local_dict[parent]] = parents_states_matrix[parent][states]; }
			for (unsigned leaf_state = 0; leaf_state < cardinality; leaf_state++) { set_node_probability(*bn, node_dic[leaf_nodes[leaf]], leaf_state, parent_state, probability[leaf_state][states]); }
		}
	}

	return true;
}

std::map<std::string, std::vector<double>> PeepoNetwork::get_inference(std::map<std::string, unsigned >& evidences)
{
	
	std::map<std::string, std::vector<double> > states;
	
	for (std::map<std::string, unsigned >::iterator it = evidences.begin(); it != evidences.end(); ++it)	{
		set_node_value(*bn, node_dic[it->first], it->second);
		set_node_as_evidence(*bn, node_dic[it->first]);
	}
	// create a solution to this bayesian network using the join tree algorithm

	join_tree_type join_tree;
	create_moral_graph(*bn, join_tree);
	create_join_tree(join_tree, join_tree);
	bayesian_network_join_tree solution(*bn, join_tree);

	std::vector<std::string> mynodes = get_nodes();
	for (unsigned node = 0; node < mynodes.size(); node++) {
		unsigned anode = node_dic[mynodes[node]];
		matrix<double> prob = solution.probability(anode);
		std::vector<double> result;
		for (auto p : prob) { result.push_back(p); }
		states[mynodes[node]] = result;
	}
	
	/*
	 // Use an approximate uinference i.e. the gibbs sampler object (!code works but  is incomplete!)
	bayesian_network_gibbs_sampler sampler;


	// To use this algorithm all we do is go into a loop for a certain number of times
	// and each time through we sample the bayesian network.  Then we count how 
	// many times a node has a certain state.  Then the probability of that node
	// having that state is just its count/total times through the loop. 

	// The following code illustrates the general procedure.

	std::vector <unsigned long> gibbs_count(number_of_nodes, 0);

	// The more times you let the loop run the more accurate the result will be.  Here we loop
	// 2000 times.
	const long rounds = 2000;
	for (long i = 0; i < rounds; ++i)
	{
		sampler.sample_graph(bn);
		for (unsigned node = 0; node < number_of_nodes; node++) {
			if (node_value(bn, node) == 1) { gibbs_count[node]++; }
		}

	}
	std::vector<double> gibbs_result(number_of_nodes, 0.0);
	std::cout << "Using the approximate Gibbs Sampler algorithm:\n";
	for (unsigned i = 0; i < number_of_nodes; i++)
	{
		gibbs_result[i] = (double)gibbs_count[i] / (double)rounds;
		std::cout << "node " << i << " -> " << gibbs_result[i] << std::endl;
	}

	*/
	return states;
}

unsigned PeepoNetwork::get_state(const matrix<double>& prob)
{
	return index_of_max(prob);
}


