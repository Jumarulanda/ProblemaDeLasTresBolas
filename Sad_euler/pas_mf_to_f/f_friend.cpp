#include "f_friend.h"

f_friend ::  modfify_string(string (*f) (string), string modi) {
	string modifai = f(modi);

	return modifai;
}
