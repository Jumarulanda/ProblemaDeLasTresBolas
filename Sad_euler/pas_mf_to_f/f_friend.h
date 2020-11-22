#include <iostream>
#include <string>

using std::string;

namespace f_friend {
	string modify_string(string (*)(string), string);
}
