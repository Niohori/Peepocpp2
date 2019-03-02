
////////////////////////////////////////////////////////////
// Headers
////////////////////////////////////////////////////////////
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>
#include "SFML/Graphics.hpp"
#include <cmath>
#include <ctime>
#include <algorithm>
#include <numeric>
#include <time.h>



#include "peepo_network.h"
#include "genetic_algorithm.h"
#include "organism.h"

#include "nlohmann/json.hpp"



using json = nlohmann::json;




std::vector<Obstacle>read_obstacles(const bool& graphical)
{
	std::vector<Obstacle> objects;
	// read a JSON file
	std::string my_path = std::experimental::filesystem::current_path().u8string() + "\\resources";
	if (!std::experimental::filesystem::exists(my_path)) { std::experimental::filesystem::create_directory(my_path); };
	my_path += "\\game_of_life_world.json";
	std::ifstream i(my_path);
	json js;
	i >> js;
	//fill container
	for (json::iterator it = js.begin(); it != js.end(); ++it) 
	{
		std::string id = it.key();
		std::string type = (*it)["type"].get<std::string>();
		double x = (*it)["x"].get<double>();
		double y = (*it)["y"].get<double>();
		Obstacle obst(id, type, x, y);
		obst.shape.setPosition(obst.x, obst.y);
		objects.push_back(obst);
	}
	return objects;
}



void generate_obstacles(const unsigned& num_food, const unsigned& num_ennemies)
{
	json js;
	for (unsigned i = 0; i < num_food; i++)
	{
		std::stringstream ss; //from <sstream>
		ss << i;
		std::string id  = "obj_"+ ss.str();
		js[id]["type"] = "food";
		js[id]["x"] = double(std::rand() % WIN_SIZE);
		js[id]["y"] = double(std::rand() % WIN_SIZE);
	}
	for (unsigned i = 0; i < num_ennemies; i++)
	{
		std::stringstream ss; //from <sstream>
		ss << i + num_food; 
		std::string id = "obj_" + ss.str();
		js[id]["type"] = "ennemy";
		js[id]["x"] = double(std::rand() % WIN_SIZE);
		js[id]["y"] = double(std::rand() % WIN_SIZE);
	}
	std::string my_path = std::experimental::filesystem::current_path().u8string() + "\\resources";
	if (!std::experimental::filesystem::exists(my_path)) { std::experimental::filesystem::create_directory(my_path); };
	my_path += "\\game_of_life_world.json";
	std::ofstream o(my_path);
	o << std::setw(4) << js << std::endl;
	o.close();
}


static class World {

public:
	//World(void) { ; };
	World(const bool& graphical_, sf::RenderWindow& window_) :
		graphical(graphical_),
		window(window_) {};
	bool graphical;
	sf::RenderWindow& window;

	void main_loop(const unsigned& max_age, const bool& verify) {

		bool isPlaying = true;
		std::vector<Obstacle> obstacles = read_obstacles(graphical);
		std::vector<double> pos = { WIN_SIZE / 2, WIN_SIZE / 2 };
		double peepo_angle = (std::rand() % 180)*PI / 180.f; // to be changed later
		// Create  grapphical peepo_
		sf::CircleShape peepo_;
		peepo_.setRadius(SIZE_PEEPO - 3);
		peepo_.setOutlineThickness(3);
		peepo_.setOutlineColor(sf::Color::Black);
		peepo_.setFillColor(sf::Color::White);
		peepo_.setOrigin(SIZE_PEEPO / 2, SIZE_PEEPO / 2);
		peepo_.setPosition(pos[0],pos[1]);


		
		// get network peepo
		PeepoNetwork network("best_life_game_network");
		//create organism
		std::string name = "archibald";
		Peepo peepo(name, network, true, pos, obstacles);
		//display loop
		while (window.isOpen())
		{
			// Handle events
			sf::Event event;
			
			while (window.pollEvent(event))
			{

				// Window closed or escape key pressed: exit
				if ((event.type == sf::Event::Closed) ||
					((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Escape)))
				{
					window.close();
					break;
				}

				// Window size changed, adjust view appropriately
				if (event.type == sf::Event::Resized)
				{
					sf::View view;
					view.setSize(WIN_SIZE, WIN_SIZE);
					view.setCenter(WIN_SIZE / 2.f, WIN_SIZE / 2.f);
					window.setView(view);
				}
			}

			if (isPlaying)
			{
				peepo.update();
				pos = peepo.pos;
				peepo_angle = peepo.rotation;
				// Move peepo
				peepo_.setPosition(pos[0], pos[1]);
			}

			// Clear the window and draw the elements
			window.clear(sf::Color(150, 150, 150));
			if (isPlaying)
			{
				window.draw(peepo_);

				/*sf::Vertex edge_right[] =
				{
					sf::Vertex(sf::Vector2f(pos[0], pos[1])),
					sf::Vertex(sf::Vector2f(peepo.edge_right[0], peepo.edge_right[1]))
				};
				window.draw(edge_right, 2, sf::Lines);

				sf::Vertex edge_left[] =
				{
					sf::Vertex(sf::Vector2f(pos[0], pos[1])),
					sf::Vertex(sf::Vector2f(peepo.edge_left[0], peepo.edge_left[1]))
				};
				window.draw(edge_left, 2, sf::Lines);

				sf::Vertex edge_middle[] =
				{
					sf::Vertex(sf::Vector2f(pos[0], pos[1])),
					sf::Vertex(sf::Vector2f(peepo.edge_middle[0], peepo.edge_middle[1]))
				};
				window.draw(edge_middle, 2, sf::Lines);*/

				for (auto obst : obstacles) {
					window.draw(obst.shape);
				}
			}
			// Display things on screen
			window.display();
		}
	}
};

static void verification(const bool& graphical)
{
	unsigned max_age = 2000;
	// Create the window of the application
	sf::RenderWindow window(sf::VideoMode(WIN_SIZE, WIN_SIZE, 32), "Game of life",
		sf::Style::Titlebar | sf::Style::Default);
	window.setVerticalSyncEnabled(false);
	World world(graphical, window);
	world.main_loop(max_age, true);
	std::cin.get();
}


int main()
{
	std::srand(static_cast<unsigned int>(std::time(NULL)));
	//create obstacle
	unsigned num_food = 200;
	unsigned num_ennemies = 200;
	//generate_obstacles(num_food, num_ennemies);
	verification(true);
	return 0;
}