#include "organism.h"



SensoryInputPeepo::SensoryInputPeepo() { ; }


SensoryInputPeepo::SensoryInputPeepo(const Peepo& peep) :
	peepo(peep)
{
	;
}


void SensoryInputPeepo::action(const std::string& node, const  std::vector<double>& prediction)
{

	if (std::max_element(prediction.begin(), prediction.end()) - prediction.begin() == 0)
	{
		peepo.motor.at(get_direction(node)) = false;
	}
	else {
		peepo.motor.at(get_direction(node)) = true;
	}
}


std::vector<double> SensoryInputPeepo::value(const std::string& name)
{
	if (name.find(VISION) != std::string::npos){
		if (peepo.view.at(get_quadrant(name))) {
			return { 0.1,0.9 };
		}
		else {
			return {0.9, 0.1};
		}
	}
	if (name.find(MOTOR) != std::string::npos) {
		if (peepo.motor.at(get_quadrant(name))) {
			return { 0.1,0.9 };
		}
		else {
			return {0.9, 0.1};
		}
	}
	if (name.find(ID_ENNEMY) != std::string::npos) {
		if (peepo.is_an_ennemy) {
			return { 0.1,0.9 };
		}
		else {
			return {0.9, 0.1};
		}
	}
	if (name.find(ID_FOOD) != std::string::npos) {
		if (peepo.is_food) {
			return { 0.1,0.9 };
		}
		else {
			return {0.9, 0.1};
		}
	}
	return { 0.5,0.5 };
}


std::string SensoryInputPeepo::get_quadrant(const std::string& name)
{
	std::string quad;
	std::vector<std::string> quadrants = { "1","2","3","4","5","6" };
	for (auto qd : quadrants)
	{
		if (name.find(qd)) { quad = qd; break; }
	}
	return quad;
}



std::string SensoryInputPeepo::get_direction(const std::string& name)
{
	std::string direction;
	for (auto dir : { LEFT, RIGHT })
	{
		if (name.find(dir)) { direction  = dir; break; }
	}
	return direction;
}


Peepo::Peepo()
{ ; }


Peepo::Peepo(const std::string& name_, const PeepoNetwork& network_, const bool& graphical_,
	std::vector<double>& pos_, std::vector<Obstacle>& obstacles_) :
	name(name_),
	network(network_),
	graphical(graphical_),
	pos(pos_),
	obstacles(obstacles_),
	generative_model(GenerativeModel(network, SensoryInputPeepo(*this))),
	rotation(0.f),
	stomach(0),
	bang(0),
	is_an_ennemy(false),
	is_food(false)
{
	motor[LEFT] = false;
	motor[RIGHT] = false;
	view["1"] = false;
	view["2"] = false;
	view["3"] = false;
	view["4"] = false;
	view["5"] = false;
	view["6"] = false;
	for (double angle = -30.0; angle < 30.0; angle += 10.0) {
		sectors.push_back({ angle, angle + 10.0 });
	}
}

void Peepo::update() 
{
	generative_model.process();
	// Move peepo
	double factor1 = PEEPO_SPEED;// *delta_time;
	pos[0] += std::cos(rotation) * factor1;
	pos[1] += std::sin(rotation) * factor1;
	if (motor[LEFT]) {
		rotation -= 10. / 180.*PI;
		if (rotation < 0.0) {
			rotation = 2.0*PI;//?? 
		}
	}
	if (motor[RIGHT]) {
		rotation += 10. / 180.*PI;
		if (rotation> 2.0*PI) {
			rotation = 0.0;
		}
	}
	calculate_obstacles();
	if (graphical) {

		//TO DO but NECSSARY?
	}
	// Check collisions between the peepo and the screen
	if (pos[0] - SIZE_PEEPO <= 0.)
	{
		rotation = PI - rotation;
		pos[0] = SIZE_PEEPO + 0.1;
	}
	if (pos[0] + SIZE_PEEPO >= WIN_SIZE)
	{
		rotation = PI - rotation;
		pos[0] = WIN_SIZE - SIZE_PEEPO - 0.1;
	}
	if (pos[1] - SIZE_PEEPO <= 0.)
	{
		rotation = -rotation;
		pos[1] =  SIZE_PEEPO + 0.1;
	}
	if (pos[1] + SIZE_PEEPO >= WIN_SIZE)
	{
		rotation = -rotation;
		pos[1] = WIN_SIZE - SIZE_PEEPO - 0.1;
	}

	edge_right = { pos[0] + PEEPO_RADIUS*std::cos(rotation + 30.0 / 180.0*PI),
				   pos[1] + PEEPO_RADIUS*std::sin(rotation + 30.0 / 180.0*PI)};
	edge_left  = { pos[0] + PEEPO_RADIUS*std::cos(rotation - 30.0 / 180.0*PI),
				   pos[1] + PEEPO_RADIUS*std::sin(rotation - 30.0 / 180.0*PI)};
}

void Peepo::calculate_obstacles(void)
{
	is_an_ennemy = false;
	is_food = false;
	for (auto vw : view) {
		view[vw.first] = false;
	}
	//check if collision occured with food or ennemy
	int to_remove = -1;
	int count = 0;
	for (auto obst : obstacles) {
		if (collision(pos, { obst.x, obst.y }, SIZE_OBST+SIZE_PEEPO))
		{
			if (obst.type == "food") {
				to_remove = count;
				stomach++;
			}
			if (obst.type == "ennemy") {
				//to_remove = count;
				bang++;
			}
		}
		count++;
	}
	//remove eaten food
	if (to_remove >= 0) {
		obstacles.erase(obstacles.begin() + to_remove);
	}
	//observations
	double closest_distance = 10000.0;
	relevant_sector.index = 0;
	for (int index = 0; index < sectors.size(); index++) {
		auto sector = sectors[index];
		double lower_edge = sector[0];
		double upper_edge = sector[1];
		for (auto obstacle : obstacles) {
			double is_collision = collision(pos, { obstacle.x, obstacle.y }, 
											rotation, lower_edge, upper_edge, PEEPO_RADIUS);
			if (is_collision > 0.0) {
				if (is_collision < closest_distance) {
					closest_distance = is_collision;
					relevant_sector.index = index + 1;
					relevant_sector.xy = { obstacle.x,obstacle.y };
					relevant_sector.distance = closest_distance;
					if (obstacle.type == "food") {
						is_an_ennemy = false;
						is_food = true;
					}
					if (obstacle.type == "ennemy") {
						is_an_ennemy = true;
						is_food = false;
					}
				}
			}
		}
	}
	std::stringstream ss; //from <sstream>
	ss << relevant_sector.index;
	std::string only_true = ss.str();
	view.at(only_true) = true;
	double sight_angle = 0.0;
	if (only_true == "0") { sight_angle = rotation; }
	if (only_true == "1") { sight_angle = rotation - 25.0 / 180.0*PI; }
	if (only_true == "2") { sight_angle = rotation - 15.0 / 180.0*PI; }
	if (only_true == "3") { sight_angle = rotation -  5.0 / 180.0*PI; }
	if (only_true == "4") { sight_angle = rotation +  5.0 / 180.0*PI; }
	if (only_true == "5") { sight_angle = rotation + 15.0 / 180.0*PI; }
	if (only_true == "6") { sight_angle = rotation + 25.0 / 180.0*PI; }
	edge_middle = { pos[0] + relevant_sector.distance*std::cos(sight_angle),
					pos[1] + relevant_sector.distance*std::sin(sight_angle) };

}