#include "genetic_algorithm.h"
//using namespace dlib;

using namespace boost::assign;
using namespace boost::adaptors;
//using namespace boost::filesystem;
using namespace boost;

Agent::Agent(void)
{
	//constructor
}

Agent::Agent(const double& afitness, const PeepoNetwork& a_network, const double& atopp, const double& acpdp) :
	fitness(afitness),
	network(a_network),
	mut_top(atopp),
	mut_cpd(acpdp)
{}

Agent& Agent::operator=(const Agent& a)
{
	fitness = a.fitness;
	PeepoNetwork network(a.network);
	mut_top = a.mut_top;
	mut_cpd = a.mut_cpd;
	network.assemble();
	return *this;

}

Agent::Agent(const Agent& a) :
	fitness(a.fitness),
	network(a.network),
	mut_top(a.mut_top),
	mut_cpd(a.mut_cpd) {
	network.assemble();
}

bool Agent::operator<(const Agent& x)
{
	return fitness < x.fitness;
}



bool Agent::operator>(const Agent& x)
{
	return fitness > x.fitness;
}


std::vector<json> get_topologies(PeepoNetwork& peepo, const unsigned& npop)
{

	std::vector<std::string> root_nodes = peepo.get_root_nodes();
	std::vector<std::string> leaf_nodes = peepo.get_leaf_nodes();
    std::vector<json> topologies;
	std::vector<std::vector<std::string> > a_topology;
	matrix<unsigned > map = ones_matrix<unsigned>(leaf_nodes.size(), root_nodes.size());
	map  = ones_matrix<unsigned>(leaf_nodes.size(), root_nodes.size());
	json a;
	// initialize random seed: 
	long long seed = std::chrono::system_clock::now().time_since_epoch().count();//seed
	std::default_random_engine dre(seed);//engine
	std::uniform_real_distribution<double> di(0.0, 2.0*0.999999999999);//uniform distribution in interval [0,2)
	int count = 0;
	while (true) {
		for (int row = 0; row < map.nr(); row++) {
			for (int edge = 0; edge < map.nc(); edge++) {
				map(row,edge) = (unsigned)floor(di(dre));
			}
		}
		bool valid = valid_graph(map, root_nodes, leaf_nodes);
		if (valid) {
			a["egdes"] = adjency_to_edges(map, root_nodes, leaf_nodes);
			a["entropy"] = adjency_to_edges(map, root_nodes, leaf_nodes).size();
			topologies.push_back(a);
		}
		if (topologies.size() > npop) { break; }
		count++;
		if (count > 2000)//safety net in order to limit the search time
		{
			topologies.insert(topologies.end(), topologies.begin(), topologies.end());//taking duplicates of already fount topologies
		}
	} //while loop
	std::random_shuffle(topologies.begin(), topologies.end());
	topologies.resize(npop-1);
	map = ones_matrix<unsigned>(leaf_nodes.size(), root_nodes.size());//adding a fully connected network
	a["egdes"] = adjency_to_edges(map, root_nodes, leaf_nodes);
	a["entropy"] = adjency_to_edges(map, root_nodes, leaf_nodes).size();
	topologies.push_back(a);
	
	//topologies.clear();
		for(unsigned i= 0; i < npop; i++)
	{
		//topologies.push_back(a);
	}

	return topologies;
}

bool valid_graph(matrix<unsigned>& map_,
	const std::vector<std::string>& root_nodes,
	const std::vector<std::string>& leaf_nodes)
{
	//check for  a minimum of edgs

	if (sum(map_) < unsigned(map_.nc())*unsigned(map_.nr())/2) { return false; }
	
	//check for orphan leaf nodes
	if (prod(sum_cols(map_)) == 0) {  return false; }


	//check for unconnected root_nodes and adapt adjecnty map accordingly
	std::vector<unsigned> loose_roots;
	matrix<unsigned> sm_rows = sum_rows(map_);
	for (int root = 0; root < sm_rows.nc(); root++) {
		if (sm_rows(root) == 0) { loose_roots.push_back(root); }
	}

	if (loose_roots.size() != 0) {
		for (unsigned i = 0; i < loose_roots.size(); i++) {
			map_ = remove_col(map_, loose_roots[i]);
		}
	}
	//the next piece ofcode uses the bfs algo to check whether the network does not contains isolated garphs 
	//(i.e. a multigraph is excluded)
	//construct a normal adjancy matrix
	unsigned R = map_.nc();
	unsigned V = R + leaf_nodes.size();

	std::vector<std::vector<unsigned> > map(V, std::vector<unsigned>(V, 0));
	for (unsigned leaf = 0; leaf < leaf_nodes.size(); leaf++) {
		for (unsigned root = 0; root < R; root++) {
			if (map_(leaf,root) == 1) {
				map[root][leaf + R] = 1;
				map[leaf + R][root] = 1;//the true adgency matrix is symmetric}
			}

		}
	}
	//set diagonal to 1
	for (unsigned row = 0; row < R + leaf_nodes.size(); row++) { map[row][row] = 1; }
	//fill the rest
	std::vector<std::vector<unsigned> >   adj(V, std::vector<unsigned>(0, 0));
	for (unsigned row = 0; row < R + leaf_nodes.size(); row++) {
		for (unsigned col = 0; col < R + leaf_nodes.size(); col++) {
			if (map[row][col] == 1) { adj[row].push_back(col); }
		}
	}

	//loop through all leaf_nodes using the bfs algo
	for (unsigned leaf = R; leaf < R + leaf_nodes.size(); leaf++) {
		unsigned s = leaf;
		// Mark all the vertices as not visited
		std::vector<bool> visited(V, false);

		// Create a queue for BFS 
		std::list<unsigned> queue;

		// Mark the current node as visited and enqueue it 
		visited[s] = true;
		queue.push_back(s);
		// 'i' will be used to get all adjacent 
		// vertices of a vertex 
		std::vector<unsigned>::iterator i;

		while (!queue.empty())
		{
			// Dequeue a vertex from queue and print it 
			s = queue.front();
			queue.pop_front();

			// Get all adjacent vertices of the dequeued 
			// vertex s. If a adjacent has not been visited, 
			// then mark it visited and enqueue it 
			for (i = adj[s].begin(); i != adj[s].end(); ++i)
			{
				if (!visited[*i])
				{
					visited[*i] = true;
					queue.push_back(*i);
				}
			}
		}
		for (unsigned i = 0; i < visited.size(); i++) {
			if (visited[i] == false) { return false; }//an isolated subgraph ahs been encontered
		}
	}//end of leaf loop
	return true;
}


matrix<unsigned> get_adjency_map(const std::vector<std::vector<std::string> >& adges,
	const std::vector<std::string>& root_nodes,
	const std::vector<std::string>& leaf_nodes)
{
	matrix<unsigned> map;
	map = zeros_matrix<unsigned>(leaf_nodes.size(),root_nodes.size());
	for (unsigned row = 0; row < root_nodes.size(); row++) {
		for (unsigned col = 0; col < leaf_nodes.size(); col++) {
			for (auto edge : adges) {
				if (edge[1] == leaf_nodes[col] && edge[0] == root_nodes[row]) {
					map(col,row) = 1;
				}
			}
		}
	}
	return map;
}


std::vector<std::vector<std::string> > adjency_to_edges(const matrix<unsigned>& map,
	const std::vector<std::string>& root_nodes,
	const std::vector<std::string>& leaf_nodes)
{
	std::vector<std::vector<std::string> > edges;
	for (unsigned col = 0; col < root_nodes.size(); col++) {
		for (unsigned row = 0; row < leaf_nodes.size(); row++) {
			if (map(row,col) == 1) { edges.push_back({ root_nodes[col],leaf_nodes[row] }); }
		}
	}
	return edges;
}

std::vector<std::vector<Agent>> get_couples(std::vector<Agent > parents)
{
	std::vector<int > agents;
	for (int i = 0; i < parents.size(); i++) {
		agents.push_back(i);
	}
	std::vector<std::vector<Agent>> mating_couples;
	Combinations<int> comb_obj(agents, 2);
	std::vector<std::vector<int>> mating_indexes = comb_obj.get();
	for (int a = 0; a < mating_indexes.size(); a++)
	{
		mating_couples.push_back({(parents[mating_indexes[a][0]]),(parents[mating_indexes[a][1]]) });
	}

	return mating_couples;

}

bool is_in_edges(const std::vector<std::vector<std::string> >& edges, const std::string&  root, const std::string& leaf)
{
	bool stop = true;
	for (auto edge : edges) {
		if (edge[0] == root && edge[1] == leaf) { stop = false; break; }
	}
	return stop;
}


GeneticAlgorithm::GeneticAlgorithm()
{
	//ctor;
	number_of_parents = int(NUMBER_OF_MATING_PARENTS);
	convergence_ratio = int(npop * POPULATION_CONVERGENCE_RATIO);
	/* initialize random seed: */
	srand(time(NULL));
	seed = std::chrono::system_clock::now().time_since_epoch().count();//seed
}

GeneticAlgorithm::GeneticAlgorithm(
	const std::string& source,
	const bool& fast,
	const unsigned& convergence_period,
	const double& convergence_sensitivity_percent,
	const unsigned&  Npop,
	const double& p_mut_top,
	const double& p_mut_cpd,
	const unsigned& max_removal) :
	fast(fast), 
	convergence_period(convergence_period),
	convergence_sensitivity(convergence_sensitivity_percent/100.0),
	npop(Npop),
	p_mut_top(p_mut_top),
	p_mut_cpd(p_mut_cpd),
	max_removal(max_removal),
	convergence_ratio(double(Npop * POPULATION_CONVERGENCE_RATIO)),
	number_of_parents(NUMBER_OF_MATING_PARENTS)
{
	number_of_parents = int(NUMBER_OF_MATING_PARENTS);
	/* initialize random seed: */
	srand(time(NULL));
	seed = std::chrono::system_clock::now().time_since_epoch().count();//seed
	initialize(source);
}


GeneticAlgorithm:: ~GeneticAlgorithm()
{
	//destructor;
}

std::vector<double> GeneticAlgorithm::generate_rand_vec(const unsigned& dim, const double& coef)
{
	std::default_random_engine dre(seed);//engine
	std::uniform_real_distribution<double> di(0.0, coef);//distribution
	std::vector<double> my_rand(dim);
	std::generate(my_rand.begin(), my_rand.end(), [&] { return di(dre); });
	return my_rand;
}




std::vector<Agent> GeneticAlgorithm::get_population(void)
{
	best_chromosome = population[0];
	best_chromosome.fitness = 0.0;
	last_generation.clear();
	last_generation.emplace(0.0, population);
	return population;
}


void GeneticAlgorithm::initialize(const std::string& origin)
{
	peepo.from_json(origin);
	root_nodes = peepo.get_root_nodes();
	leaf_nodes = peepo.get_leaf_nodes();
	cardinality_map = peepo.get_cardinality_map();
	std::vector<json> topologies = get_topologies(peepo, npop);
	for (unsigned t = 0; t < topologies.size();t++) {
		peepo.edges = topologies[t]["egdes"].get<std::vector<std::vector<std::string> > >();
		
		for (auto node : peepo.get_nodes()) {
			std::vector<std::string> incoming_nodes = peepo.get_incoming_edges(node);
			t_cpd_root  root_cpd;
			t_cpd_leaf leaf_cpd;
			std::vector<double> my_omega;
			if (incoming_nodes.size() == 0) {
				for (unsigned n = 0; n < cardinality_map[node]; n++) {
					root_cpd.push_back(1.0 / double(cardinality_map[node]));
					my_omega.push_back(0.0);
				}
			}
			else {
				std::vector<unsigned> my_card_parents;
				unsigned prod = 1;
				for (auto nod : incoming_nodes) { my_card_parents.push_back(cardinality_map[nod]); prod *= cardinality_map[nod]; };
				double max_omega = 2.0*PI*prod;
				std::vector<double> my_omega = generate_rand_vec(cardinality_map[node], max_omega);
				leaf_cpd = ga_child_cpd(my_card_parents, my_omega);
			}
			if (root_cpd.size() == 0){
				peepo.add_cpds(node, leaf_cpd);
			}
			else {
				peepo.add_cpds(node, root_cpd);
			}
			peepo.add_omega_map(node, my_omega);
		}//end of nodes loop
		peepo.assemble();
		PeepoNetwork a_peepo(peepo);//copy constructor
		Agent agent(0.0, a_peepo);
		population.push_back(agent);
		peepo.disassemble();
	}//end of topologies loop
	if (convergence_ratio == 0.0) { convergence_ratio = 1.0; };
	}


void GeneticAlgorithm::get_parents(void)
{
	selected_parents.clear();
	std::sort(population.begin(), population.end(), [](Agent a, Agent b) {return a.fitness > b.fitness; });
	if (population[0].fitness >= best_chromosome.fitness) {
		best_chromosome = population[0];
	}
	//calculate the fitness score over the best self.convergence_ratio best performing chromosomes
	population_fitness_score = 0.0;
	for (unsigned i = 0; i < convergence_ratio; i++) {
		population_fitness_score += population[i].fitness / convergence_ratio;
	}
	
	if (fast) {
		for (int i = 0; i < number_of_parents; i++) {
			selected_parents.push_back(population[i]);
			selected_parents[i].fitness = 0.0;
			population[i].fitness = 0.0;
		}
	}
	else {
		std::vector<unsigned> pool;
		for (int index = 0; index < npop; index++) {
			int repeat = npop - index;
			for (int i = 0; i < repeat; i++) {
				pool.push_back(index);
			}
		}
		std::random_shuffle(pool.begin(), pool.end());
		for (int draw = 0; draw < number_of_parents; draw++) {
			unsigned pool_index = std::rand() % (pool.size() - 1);
			unsigned parent_index = pool[pool_index];
			population[parent_index].fitness = 0.0;
			selected_parents.push_back(population[parent_index]);
		}
	}
	std::random_shuffle(selected_parents.begin(), selected_parents.end());
}

Agent GeneticAlgorithm::get_optimal_newtork(void)
{
	return best_chromosome;
}

bool GeneticAlgorithm::check_convergence(void)
{
	if (convergence_history.size() <= convergence_period) {
		convergence_history.push_back(population_fitness_score);
	}
	else
	{
		convergence_history.push_back(population_fitness_score);
		convergence_history.erase(convergence_history.begin());
		double mean = 0.0;
		double std = 0.0;
		{
			unsigned len = convergence_history.size();
			for (auto fit : convergence_history) { [&mean, len](const double& fit) {mean += fit / len; }; };
			for (auto fit : convergence_history) { [&std, mean](const double& fit) {std += (fit - mean)*(fit - mean); }; }
			std = std::sqrt(std / len);
		}//scope bracket
		double acceptable_std = mean * convergence_sensitivity / 2.0;
		if (std < acceptable_std) { return true; }
	}
	double prev_fitness = last_generation.begin()->first;
	double actual_fitness = population_fitness_score;
	if (actual_fitness < prev_fitness) {
		actual_fitness = (prev_fitness + actual_fitness) / 2.0;
		std::vector<Agent> actual_parents = copy_agents(selected_parents);
		selected_parents = copy_agents(last_generation.begin()->second);
		selected_parents.insert(std::end(selected_parents), std::begin(actual_parents), std::end(actual_parents));
		std::random_shuffle(selected_parents.begin(),selected_parents.end());
		selected_parents.resize(npop);
	}
	last_generation.clear();
	last_generation.insert(std::pair<double, std::vector<Agent> >(actual_fitness, selected_parents));
	selected_parents.resize(number_of_parents);
	return false;
}


std::vector<Agent> GeneticAlgorithm::copy_agents(const std::vector<Agent>& in_agents)
{
	std::vector<Agent> agents;
	for (auto par : in_agents) { 
		agents.push_back(Agent(par)); }
	for (int a = 0; a < agents.size(); a++) {
		agents[a].network.assemble();
	}
	return agents;
}

std::vector<Agent> GeneticAlgorithm::evolve(std::vector<Agent> in_population, double& in_fitness, bool& converging)
{
	population = copy_agents(in_population);
	converging = false;
	if (population.size() == 0) { in_fitness = -1; return in_population;	}
	get_parents();
	converging = check_convergence();
	if (converging) { in_fitness = population_fitness_score; return in_population; }
	if(population_fitness_score < 0) { in_fitness = -1.0;  return in_population; }
	//Cross-over
	std::vector<Agent> selected_offsprings = cross_over();
	while (true) {//expand the population to npop; the added peepo which are clones, will be mutated
		int n_chrom = selected_parents.size() + selected_offsprings.size();
		if (n_chrom > npop) { break; }
		for (auto x : selected_parents) {
			Agent y(x.fitness, PeepoNetwork(x.network), 0.0,0.0);
			y.network.assemble();
			selected_offsprings.push_back(y);
		}
	}
	//Mutation
	std::vector<Agent> selected_offsprings_c;
	for (unsigned offspr = 0; offspr < selected_offsprings.size(); offspr++) {
		Agent offspring(selected_offsprings[offspr]);
		double m_top = offspring.mut_top;
		double m_cpd = offspring.mut_cpd;
		//check if treshold is met and mutate accordingly
		if (m_top < p_mut_top) {
			PeepoNetwork nk = prune_or_grow(offspring.network);
			offspring.network.edges = nk.edges;
			offspring.network.cpds = nk.cpds;
			offspring.network.omega_map = nk.omega_map;
			offspring.network.assemble();
		}
		
		if (m_cpd < p_mut_cpd) {
			PeepoNetwork nk = mutate_cpds(offspring.network);
			offspring.network.edges = nk.edges;
			offspring.network.cpds = nk.cpds;
			offspring.network.omega_map = nk.omega_map;
			offspring.network.assemble();
		}
		offspring.network.assemble();
		selected_offsprings_c.push_back(offspring);
	}
	selected_offsprings.clear();
	selected_offsprings = selected_offsprings_c;
	//collect parents and offsprings
	std::random_shuffle(selected_offsprings.begin(), selected_offsprings.end());
	population.clear();
	for (auto par : selected_parents) {	population.push_back(par);}
	for(auto off : selected_offsprings) { population.push_back(off); }
	//prune the population to npop
	population.resize(npop);
	in_fitness = population_fitness_score;
	return copy_agents(population);
}



std::vector<Agent> GeneticAlgorithm::cross_over(void)
{

	std::vector<Agent> selected_offsprings;

	std::vector<std::string> rt_nodes = selected_parents[0].network.get_root_nodes();
	std::vector<std::string> lf_nodes = selected_parents[0].network.get_leaf_nodes();
	std::vector<std::vector<Agent>> mating_couples = get_couples(selected_parents);
	if (mating_couples.size() > npop) {
		std::random_shuffle(mating_couples.begin(), mating_couples.end());
		mating_couples.resize(npop);
	}
	for (auto chrom : mating_couples) {
		
		matrix<unsigned> map_1 = get_adjency_map(chrom[0].network.get_edges(), chrom[0].network.get_root_nodes(), chrom[0].network.get_leaf_nodes());
		matrix<unsigned> map_2 = get_adjency_map(chrom[1].network.get_edges(), chrom[1].network.get_root_nodes(), chrom[1].network.get_leaf_nodes());
		matrix<int> diff = abs(matrix_cast<int>(map_1) - matrix_cast<int>(map_2));
		if (sum(diff) == 0 || sum(diff) > int(treshold)) {
		//-> there is no difference in topology between the parents, two offsprings, identical to the parents are created but they will have 100% chance to mutate,  both in topology as for cpd
		// or the difference between parents is too big.We assume cloning of the  parents(will be mutated!)
			Agent ag_1(chrom[0]);
			Agent ag_2(chrom[1]);
			ag_1.mut_cpd = 0.0;
			ag_1.mut_top = 0.0;
			ag_2.mut_cpd = 0.0;
			ag_2.mut_top = 0.0;
			selected_offsprings.push_back(ag_1);
			selected_offsprings.push_back(ag_2);
			continue;
		}
		//we now construct offsprings for all other cases:
		//if there are q positions with differnet values in the two adjency matrices, thean, we will   have 2 ^ q - 2 offsprings
		std::vector<std::vector<int>> indices;
		for (int row = 0; row < diff.nr(); row++) {
			for (int col = 0; col < diff.nc(); col++) {
				if (diff(row, col) == 1) {
				indices.push_back({ row,col });//containes the amtrix indiced with different values in map_1 & map_2
				}
			}
		}
		std::vector<std::vector<int> > combination;
		combination.push_back({ 1,0 });//this is the possible combinations when we have only 1 difference
		combination.push_back({ 0,1 });//this is the possible combinations when we have only 1 difference
		if (indices.size() > 1) {
			std::vector<int> combo(indices.size(), 2);
			combination = States::transpose(States::get_index_matrix(combo));
		}
		for (auto comb : combination) {
			matrix<unsigned> map = map_1;
			for (unsigned pos = 0; pos < indices.size(); pos++) {
				map(indices[pos][0], indices[pos][1]) = comb[pos];
			}
			//check if some leaf nodes get orphan, if yes we set skip to True and skip this combination
			matrix<unsigned> check_map = sum_cols(map);//TO CHECK!!or rows??
			bool skip = false;
			if (min(check_map) == 0) {
				skip = true;
			}
			//check if this combination is equal to one of the parents, if yes we set skip to True and skip this combination
			if (map == map_1 || map == map_2) {
				skip = true;
			}
			if (skip) { continue; }

			//for each map(combination) we create an offspring
			std::vector<std::vector<std::string> >  edges = adjency_to_edges(map, rt_nodes, lf_nodes);
			PeepoNetwork a_peepo(chrom[0].network);
			a_peepo.disassemble();
			a_peepo.edges = edges;

			for (auto node : a_peepo.get_nodes()) {
				std::vector<std::string> incoming_nodes = a_peepo.get_incoming_edges(node);
				t_cpd_root my_cpd_root;
				t_cpd_leaf my_cpd_leaf;
				std::vector<double> my_omega;
				if (incoming_nodes.size() == 0) {
					my_cpd_root = t_cpd_root(a_peepo.cardinality_map[node].get<unsigned>(), 1.0 / a_peepo.cardinality_map[node].get<unsigned>());
					my_omega.clear();
				}
				else {
					std::vector<unsigned> my_card_parents;
					unsigned prod_card = 1;
					for (auto nod : incoming_nodes) {
						my_card_parents.push_back(a_peepo.cardinality_map[nod].get<unsigned>());
						prod_card *= a_peepo.cardinality_map[nod].get<unsigned>();
					}
					double max_omega = 2.0*PI*prod_card;
					my_omega = generate_rand_vec(a_peepo.cardinality_map[node].get<unsigned>(), max_omega);
					my_cpd_leaf = ga_child_cpd(my_card_parents, my_omega);
				}
				if (my_cpd_root.size() == 0) {
					a_peepo.add_cpds(node, my_cpd_leaf);
				}
				else {
					a_peepo.add_cpds(node, my_cpd_root);
				}
				a_peepo.add_omega_map(node, my_omega);
				}
				a_peepo.assemble();

				//check whether the combination is a cloning of one parent, if YES mutation will occur
				double mut_top = double(std::rand() % 100) / 100.0;
				double mut_cpd = double(std::rand() % 100) / 100.0;
				if (comb[0] == comb[1]) {
				mut_top = 0.0;
				mut_cpd = 0.0;
				}
				Agent pp(0.0, a_peepo, mut_top, mut_cpd);
				selected_offsprings.push_back(pp);
		}//combination loop
	}//mating_couples loop
	return selected_offsprings;
}


PeepoNetwork GeneticAlgorithm::prune_or_grow(PeepoNetwork network)
{
	std::vector<std::string> nodes = network.get_nodes();
	std::vector<std::string> rt_nodes = network.get_root_nodes();
	std::vector<std::string> lf_nodes = network.get_leaf_nodes();
	std::vector<std::vector<std::string> > edges = network.get_edges();
	std::vector<std::string> incoming_edges;
	std::string leaf;
	std::string root;
	std::string  p_or_g = "prune";
	bool stop = false;
	while (true) {//loop while finding a acceptable pait of root & leafs to prrune or grow
		std::random_shuffle(lf_nodes.begin(), lf_nodes.end());
		leaf = lf_nodes[0];
		incoming_edges = network.get_incoming_edges(leaf);
		if (incoming_edges.size() == 1) {//no pruning allowed
			p_or_g = "grow";
			while (true) {//choose a root node not yet reached by the leaf
				std::random_shuffle(rt_nodes.begin(), rt_nodes.end());
				root = rt_nodes[0];
				stop = is_in_edges(edges, root, leaf);
				if (stop) { break;}
			}
		}
		if (stop) { break;}
		if (incoming_edges.size() == rt_nodes.size()) {//no grow allowed
			p_or_g = "prune";
			std::random_shuffle(rt_nodes.begin(), rt_nodes.end());
			root = rt_nodes[0];
			break;
		}
		//choosing randomly prune or grow
		int k = std::rand() % 100;
		if (k > 49) { p_or_g = "grow"; }
		if (p_or_g == "prune") {
			std::random_shuffle(incoming_edges.begin(), incoming_edges.end());
			root = incoming_edges[0];
			break;
		}
		else {
			while (true) {//choose a root node not yet reached by the leaf
				std::random_shuffle(rt_nodes.begin(), rt_nodes.end());
				root = rt_nodes[0];
				stop = is_in_edges(edges, root, leaf);
				if (stop) { break; }
			}
			break;
		}
	}

	if (p_or_g == "prune") {
		std::vector < std::vector<std::string>>  edges_;
		for (auto edge : edges) {
			if(!(edge[0] == root &&  edge[0] == leaf)){edges_.push_back(edge); }
		}
		edges = edges_;
	}
	else {
		edges.push_back({ root,leaf });
	}
	//json cpds = network.cpds;
	json omega_map = network.omega_map;
	network.set_edges(edges);
	incoming_edges = network.get_incoming_edges(leaf);
	std::vector<unsigned> parents_card;
	for (auto x : incoming_edges) {
		parents_card.push_back(cardinality_map[x]);
	}
	t_cpd_leaf new_cpd = ga_child_cpd(parents_card, omega_map[leaf]);
	network.add_cpds(leaf, new_cpd);
	//network.omega_map = omega_map;
	network.assemble();
	return network;
}


PeepoNetwork GeneticAlgorithm::mutate_cpds(PeepoNetwork network)
{
	std::vector<std::string> leafs = network.get_leaf_nodes();
	json omega_map = network.omega_map;
	json card_map = network.get_cardinality_map();
	double epsilon = 0.05 + double(std::rand() % 1000 / 1000.0*0.7);
	for (auto leaf : leafs) {
		std::vector<std::string> incoming_edges = network.get_incoming_edges(leaf);
		std::vector<unsigned> parents_card;
		for (auto x : incoming_edges) {
			parents_card.push_back(card_map[x]);
		}
		std::vector<double>  my_omega = omega_map[leaf].get<std::vector<double>>();
		for (unsigned i = 0; i < my_omega.size(); i++) {
			my_omega[i] += (0.5 - double(std::rand() % 1000) / 1000.0)*epsilon;
		}
		t_cpd_leaf my_cpd = ga_child_cpd(parents_card, my_omega);
		network.omega_map[leaf] = my_omega;
		network.cpds[leaf] = my_cpd;
	}
	network.assemble();
	return PeepoNetwork(network);
}



t_cpd_leaf GeneticAlgorithm::ga_child_cpd(const std::vector<unsigned>& card_parents, const std::vector<double>& omega)
{
	t_cpd_leaf pdf;
	double phase_shift = omega[0];
	unsigned n_comb = 1;
	for (auto card : card_parents) { n_comb *= card; };
	for (auto ang : omega)
	{
		t_cpd_root pdf_row;
		for (unsigned col = 0; col < n_comb; col++){pdf_row.push_back(sin(ang*(col + 1) + phase_shift) + 1.2);}
		pdf.push_back(pdf_row);
	}
	return normalize_distribution(pdf);
}


t_cpd_leaf GeneticAlgorithm::normalize_distribution(t_cpd_leaf& matrix_)
{
	t_cpd_leaf matrix = matrix_;
	std::vector<double> sum(matrix[0].size(), 0.0);
	for (unsigned column = 0; column < sum.size(); column++)
	{
		for (unsigned row = 0; row < matrix.size(); row++)
		{
			sum[column] += matrix[row][column];
		}
	}
	for (unsigned column = 0; column < sum.size(); column++)
	{
		for (unsigned row = 0; row < matrix.size(); row++)
		{
			matrix[row][column] /= sum[column];;
		}
	}
	return matrix;
}

