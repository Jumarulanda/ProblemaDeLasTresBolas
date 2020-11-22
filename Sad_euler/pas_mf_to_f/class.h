#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "f_friend.h"

using std::string;
using std::cin;
using std::endl;
using std::cout;

class test_class {
	public:
		test_class();

		void print_something(string);

	friend string f_friend :: modify_string(string (*)(string), string);

	private:
		string f_to_print(string);
};
