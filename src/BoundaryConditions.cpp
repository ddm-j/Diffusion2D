#include <BoundaryConditions.h>

static bool isLuaScript(const std::string& input) {
	// Detect if Passed string is a lua script file path or Lua code
	std::string extension = ".lua";
	if (input.length() >= extension.length()) {
		return input.substr(input.length() - extension.length()) == extension;
	}
	return false;
}

void BoundaryConditions::initializeLua() {
	// Initialize the Lua VM
	if (L == nullptr) { // Check if Lua state has already been created
		L = luaL_newstate();
		luaL_openlibs(L);
		std::cout << "Lua state initialized." << std::endl;
	}
	else {
		std::cout << "Lua state already initialized." << std::endl;
	}
}

void BoundaryConditions::loadLuaFunction(const std::string funcName, const std::string UDF, const bool is_file, int& ref) {
	// Load a lua function from code or file
	if (is_file) {
		luaL_dofile(L, UDF.c_str());
	}
	else {
		luaL_dostring(L, UDF.c_str());
	}

	// Push the Function onto the stack
	lua_getglobal(L, funcName.c_str());

	// Get the reference
	if (lua_isfunction(L, -1)) {
		ref = luaL_ref(L, LUA_REGISTRYINDEX);
	}
	else {
		throw std::invalid_argument(funcName + " is not a function.");
	}
}

double BoundaryConditions::callLuaFunction(int& ref, double arg) {
	// Push Function onto stack
	lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
	
	// Push Argument onto stack
	lua_pushnumber(L, arg);

	// Call the function
	if (lua_pcall(L, 1, 1, 0) != 0) {
		throw std::runtime_error(lua_tostring(L, -1));
	}

	double res = static_cast<float>(lua_tonumber(L, -1));
	lua_pop(L, 1);
	return res;
}

void BoundaryConditions::updateBCs(const double arg) {
	// Call all UDF BCs and update their values
	for (int index : udfBounds) {
		bcvs[index] = callLuaFunction(luaRefs[index], arg);
	}
}

void BoundaryConditions::addBC(char face, const int type, const double value) {
	// Add a UDF boundary condition (Dirichlet or Neumann)
	face = toupper(face);

	// Check the Input
	if (std::string("NSEW").find(face) == std::string::npos) {
		throw std::invalid_argument("Invalid boundary face specified. Should be either N, S, E, or W.");
	}
	if ((type < 0) || (type > 1)) {
		throw std::invalid_argument("Invalid BC type specified. A standard boundary is either 0-Dirichlet or 1-Neumann.");
	}

	// Add Input
	int face_index = bcMap[face];
	bcs[face_index] = type;
	bcvs[face_index] = value;
}

void BoundaryConditions::addBC(char face, const int type, const std::string funcName, const std::string UDF) {
	// Add a standard boundary condition (Dirichlet or Neumann)
	face = toupper(face);

	// Check the Input
	if (std::string("NSEW").find(face) == std::string::npos) {
		throw std::invalid_argument("Invalid boundary face specified. Should be either N, S, E, or W.");
	}
	if ((type < 2) || (type > 3)) {
		throw std::invalid_argument("Invalid BC type specified. A UDF boundary is either 2-Dirichlet or 3-Neumann.");
	}

	// Initialize the Lua State if Necessary
	initializeLua();

	// Load the Function
	has_UDFs = true;
	int face_index = bcMap[face];
	udfBounds.push_back(face_index);
	loadLuaFunction(funcName, UDF, isLuaScript(UDF), luaRefs[face_index]);

	std::cout << "Lua Function loaded: " << luaRefs[face_index] << std::endl;

	// Get the Initial Boundary Condition Value
	bcs[face_index] = type;
	bcvs[face_index] = callLuaFunction(luaRefs[face_index], 1.0);
}