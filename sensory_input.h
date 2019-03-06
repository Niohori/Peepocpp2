#ifndef SENSORY_INPUT__H
#define SENSORY_INPUT__H


class SensoryInput {
public:
	//virtual ~SensoryInput() = 0;
	virtual void action(const std::string& a, const std::vector<double>& ) = 0;
	virtual std::vector<double> value(const std::string& a) = 0;
	virtual std::shared_ptr<SensoryInput> clone()  = 0;
	virtual std::map<std::string, bool> get_motor() = 0;
};

#endif
