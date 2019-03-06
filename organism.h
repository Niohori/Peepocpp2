#ifndef ORGANISM__H
#define ORGANISM__H


#include "peepo_network.h"
#include "sensory_input.h"
#include "generative_model.h"
#include "vision.h"
#include "SFML/Graphics.hpp"
#include <cmath>

const std::string LEFT = "left";

const std::string RIGHT = "right";

const std::string VISION = "vision";

const std::string MOTOR = "motor";

const std::string ID_ENNEMY = "ennemy";

const std::string ID_FOOD = "food";


// Define some global constants
const int WIN_SIZE = 900;
const double SIZE_OBST = 5.f;
const double SIZE_PEEPO = 5.f;
const double PEEPO_RADIUS = 50.0;
const double PI = 3.14159f;
const double PEEPO_SPEED = 2.f;


struct Obstacle
{
	Obstacle(const std::string& id, const std::string& type, const double& x, const double& y) :
		id(id),
		type(type),
		x(x),
		y(y) {
		if (type == "food") {
			shape.setRadius(SIZE_OBST - 3);
			shape.setOutlineThickness(3);
			shape.setOutlineColor(sf::Color::Black);
			shape.setFillColor(sf::Color::Green);
			shape.setOrigin(SIZE_OBST / 2, SIZE_OBST / 2);
		}
		if (type == "ennemy") {
			shape.setRadius(SIZE_OBST - 3);
			shape.setOutlineThickness(3);
			shape.setOutlineColor(sf::Color::Black);
			shape.setFillColor(sf::Color::Red);
			shape.setOrigin(SIZE_OBST / 2, SIZE_OBST / 2);
		}
	}
	std::string id;
	std::string type;
	double x;
	double y;
	sf::CircleShape shape;
};



class Peepo
{

private:
	//members
	//std::vector<Obstacle> obstacles;
	std::vector<std::vector<double> > sectors;
	struct RelevantSector {
		int index;
		std::vector<double> xy;
		double distance;
	};
	RelevantSector relevant_sector;

public:
	//members
	std::string name;
	PeepoNetwork network;
	bool graphical;
	std::shared_ptr<SensoryInput> sensory_input;

	std::vector<double> pos;
	std::vector<double> edge_right;
	std::vector<double> edge_left;
	std::vector<double> edge_middle;

	std::vector<Obstacle> obstacles;
	std::map<std::string, bool> motor;
	std::map<std::string, bool> view;
	bool is_an_ennemy;
	bool is_food;
	std::shared_ptr<GenerativeModel> generative_model;
	unsigned stomach;
	unsigned bang;
	double rotation;

	//TO COMPLETE




public:
	//functions
	Peepo();
	Peepo(const Peepo&);
	Peepo(const std::string&, const PeepoNetwork&, const bool&,  std::vector<double>&, std::vector<Obstacle>&);
	void assemble_obstacles(void);
	void update();
	void draw(void);
	//void make_image();//
	void calculate_obstacles(void);
};



class SensoryInputPeepo :public SensoryInput
{
private:
	//Peepo peepo;

public:
	Peepo peepo;
	SensoryInputPeepo();
	SensoryInputPeepo(Peepo&);
	std::shared_ptr<SensoryInput> clone() ;
	void action(const std::string&, const std::vector<double>&);
	std::map<std::string, bool> get_motor();
	std::vector<double> value(const std::string&);
	static std::string get_quadrant(const std::string&);
	static std::string get_direction(const std::string&);
};




#endif