all: carmfind.cpp
	g++ -std=c++17 carmfind.cpp -lgmpxx -lgmp -o carmfind
