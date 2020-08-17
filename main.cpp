#include <iostream>
#include "MeshGeneration.h"

int main(){

	std::cout << "Start Mesh Generation." << std::endl;

	MeshGeneration* meshGen = new MeshGeneration;
	meshGen->readInputFile(); // Read data from input file
	meshGen->calcMeshData(); // Generate mesh data
	meshGen->calcResisivityDistribution(); // Calculate resistivity distribution

	meshGen->outputMeshData(); // Output mesh data
	//meshGen->outputOCD(); // Output OCD file
	meshGen->outputVTK(); // Output VTK file
	delete meshGen;

	std::cout << "End Mesh Generation." << std::endl;

}
