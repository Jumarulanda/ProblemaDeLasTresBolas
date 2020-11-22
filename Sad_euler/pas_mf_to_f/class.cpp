#include "class.h"

test_class :: test_class() {
	cout << "Se ha construido la clase, tío" << endl;
}

void test_class :: print_something(string argument) {
	string f = [] (string estri) -> string { return this -> f_to_print(estri); };
	string modi_arg = f_friend :: modify_string(f,argument);
	cout << modi_arg <<  endl;
	/* cout << f_to_print(argument) << endl; */
}

string test_class :: f_to_print(string argument) {
	argument += ": eeeh, argumentado";
	return argument;
}
