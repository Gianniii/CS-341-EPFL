
#include "lsystem.h"
#include <stack>
#include <memory>
#include <iostream>

/*
Provided utilities:

- Dice class (utils/misc.h)
	Produces random values uniformly distributed between 0 and 1
	Example:
		Dice d;
		double random_val = d.roll();

- write_string_to_file (utils/misc.h)
	Write string data into a text file.
	Example:
		write_string_to_file("ala ma kota!", "ala.txt");
*/

std::string LindenmayerSystemDeterministic::expandSymbol(unsigned char const& sym) {
	/*============================================================
		TODO 1.1
		For a given symbol in the sequence, what should it be replaced with after expansion?
		The rules are in this->rules, see lsystem.h for details.
	*/
	auto r = rules.find(sym);
	if(r == rules.end()){
		return {char(sym)}; // this constructs string from char 
	}else{
		return r->second;
	}

	/*
	You may find useful:
		map.find: Iterator to an element with key equivalent to key. If no such element is found, past-the-end (see end()) iterator is returned.
		http://en.cppreference.com/w/cpp/container/unordered_map/find
	============================================================
	*/
}

std::string LindenmayerSystem::expandOnce(std::string const& symbol_sequence) {
	/*============================================================
		TODO 1.2
		Perform one iteration of grammar expansion on `symbol_sequence`.
		Use the expandSymbol method
	*/
	std::string s = "";

	for (char c : symbol_sequence){
		s += expandSymbol(c);
	}

	return s;

	//============================================================
}

std::string LindenmayerSystem::expand(std::string const& initial, uint32_t num_iters) {
	/*============================================================
		TODO 1.3
		Perform `num_iters` iterations of grammar expansion (use expandOnce)
	*/
	std::string s = initial;
	for (uint32_t i = 0; i < num_iters; i++)
	{
		s = expandOnce(s);
	}

	return s;
	

	//============================================================
}

std::vector<Segment> LindenmayerSystem::draw(std::string const& symbols) {
	/*============================================================
		TODO 2.1
		Build line segments according to the sequence of symbols
		The initial position is (0, 0) and the initial direction is "up" (0, 1)
		Segment is std::pair<vec2, vec2>
	*/
	
	std::vector<Segment> drawing;
	vec2 curr_pos = vec2(0,0);
	vec2 curr_dir = vec2(0,1);
	double delta = deg2rad(rotation_angle_deg);
	float cos_d = cosf(delta);
	float sin_d = sinf(delta);
	float cos_dm = cosf(-delta);
	float sin_dm = sinf(-delta);
	mat2 rot_mat_p = mat2(cos_d, -sin_d, sin_d, cos_d);
	mat2 rot_mat_m = mat2(cos_dm, -sin_dm, sin_dm, cos_dm);
	std::stack<Segment> stack;
	Segment seg;
	vec2 next;

	for (char c : symbols)
		switch (c){
			
		case '+':{
			curr_dir = normalize(rot_mat_p * curr_dir);
			break;
		}

		case '-':{
			curr_dir = normalize(rot_mat_m * curr_dir);
			break;
		}

		case '[':{
			stack.push({curr_pos, curr_dir});
			break;
		}

		case ']':{
			seg = stack.top();
			stack.pop();
			curr_pos = seg.first;
			curr_dir = seg.second;
			break;
		}
		
		default:{
			next = curr_pos + curr_dir;
			drawing.push_back({curr_pos,next});
			curr_pos = next;
		}

	}
	return drawing;
	
	//============================================================
}

std::string LindenmayerSystemStochastic::expandSymbol(unsigned char const& sym) {
	/*============================================================
		TODO 4.1
		For a given symbol in the sequence, what should it be replaced with after expansion?
		(stochastic case)
		The rules are in this->rules, but now these are stochastic rules because this method belongs to the LindenmayerSystemStochastic class, see lsystem.h for details.
	*/
	
	
	auto r = rules.find(sym);
	if(r == rules.end()){
		return {char(sym)}; // this constructs string from char 
	}else{
		double d = dice.roll();
		double total_p = 0.0;
		std::vector<StochasticRule> s_rules = r->second;
		for(size_t i = 0; i < s_rules.size(); i++){
			double prob = s_rules[i].probability;
			if (d > total_p && d <= prob + total_p){
				return s_rules[i].expansion;
			} else
			{
				total_p += prob;
			}
						
		}
		return s_rules[0].expansion; //if dice.roll = 0
	}
	//============================================================
}

void LindenmayerSystemDeterministic::addRuleDeterministic(unsigned char sym, std::string const& expansion) {
	rules[sym] = expansion;
}

void LindenmayerSystemStochastic::addRuleStochastic(unsigned char sym, std::vector<StochasticRule> expansions_with_ps) {
	rules[sym] = expansions_with_ps;
}
