#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "MeshGeneration.h"

double MeshGeneration::EPS = 1e-12;
double MeshGeneration::distanceConversion = 1000;

// Constructer
MeshGeneration::MeshGeneration():
	m_numX(0),
	m_numY(0),
	m_numZ(0),
	m_numElemTotal(0),
	m_numNodeTotal(0),
	m_CoordinatesX(NULL),
	m_CoordinatesY(NULL),
	m_CoordinatesZ(NULL),
	m_orderOfNumbering(XYZ),
	m_numLayers(0),
	m_numResistivityBlocks(0),
	m_elemGroupingZ(NULL),
	m_numElemGroupX(NULL),
	m_numElemGroupY(NULL),
	m_elemGroupingX(NULL),
	m_elemGroupingY(NULL),
	m_neighborElements(NULL),
	m_locXYZ(NULL),
	m_nodesOfElements(NULL),
	m_numResisivityBlockAccumulated(NULL),
	m_resistivityBlockID(NULL),
	m_XYZ2ElemID(NULL),
	m_initialResistivity(100.0),
	//m_airResistivity(1.0e12),
	m_airResistivity(-1.0),
	m_ResistivityValues(NULL),
	m_fixResistivityValues(NULL),
	m_numResisivityAnomalies(0),
	m_anomalyData(NULL),
	m_hasDivisionNumberRead(false),
	m_hasXCoordinatesRead(false),
	m_hasYCoordinatesRead(false),
	m_hasZCoordinatesRead(false),
	m_hasNumberingMethodRead(false),
	m_hasLayersRead(false),
	m_hasLayerRead(NULL),
	m_hasInitialResistivityRead(false),
	m_hasAirResistivityRead(false),
	m_hasAnomaliesRead(false),
	m_meshType(HEXA)
{}

// Destructer
MeshGeneration::~MeshGeneration(){

	if( m_CoordinatesX != NULL ){
		delete[] m_CoordinatesX;
	}

	if( m_CoordinatesY != NULL ){
		delete[] m_CoordinatesY;
	}

	if( m_CoordinatesZ != NULL ){
		delete[] m_CoordinatesZ;
	}

	if( m_elemGroupingZ != NULL ){
		delete[] m_elemGroupingZ;
	}
	
	if( m_numElemGroupX != NULL ){
		delete[] m_numElemGroupX;
	}

	if( m_numElemGroupY != NULL ){
		delete[] m_numElemGroupY;
	}

	if( m_elemGroupingX != NULL ){
		for( int i = 0 ; i < m_numLayers; ++i ){
			if( m_elemGroupingX[i] != NULL ){
				delete[] m_elemGroupingX[i];
			}
		}
		delete[] m_elemGroupingX;
	}

	if( m_elemGroupingY != NULL ){
		for( int i = 0 ; i < m_numLayers; ++i ){
			if( m_elemGroupingY[i] != NULL ){
				delete[] m_elemGroupingY[i];
			}
		}
		delete[] m_elemGroupingY;
	}

	if( m_neighborElements != NULL ){
		delete[] m_neighborElements;
	}

	if( m_locXYZ != NULL ){
		delete[] m_locXYZ;
	}

	if( m_nodesOfElements != NULL ){
		delete[] m_nodesOfElements;
	}

	if( m_numResisivityBlockAccumulated != NULL ){
		delete[] m_numResisivityBlockAccumulated;
	}

	if( m_resistivityBlockID != NULL ){
		delete[] m_resistivityBlockID;
	}

	if( m_XYZ2ElemID != NULL ){
		for( int ix = 0 ; ix < m_numX; ++ix ){
			for( int iy = 0 ; iy < m_numY; ++iy ){
				if( m_XYZ2ElemID[ix][iy] != NULL ){
					delete[] m_XYZ2ElemID[ix][iy];
				}
			}
			if( m_XYZ2ElemID[ix] != NULL ){
				delete[] m_XYZ2ElemID[ix];
			}
		}
		delete[] m_XYZ2ElemID;
	}

	if( m_ResistivityValues != NULL ){
		delete[] m_ResistivityValues;
	}

	if( m_fixResistivityValues != NULL){
		delete [] m_fixResistivityValues;
		m_fixResistivityValues = NULL;
	}

	if( m_hasLayerRead != NULL ){
		delete[] m_hasLayerRead;
	}

	if( m_anomalyData != NULL ){
		delete[] m_anomalyData;
	}

}

void MeshGeneration::readInputFile(){

	std::cout << "Read Input file." << std::endl;	
	
	// Read parameters from input file
	std::ifstream inputFile( "meshgen.inp", std::ios::in );
	if( inputFile.fail() )
	{
		std::cerr << "File open error : meshgen.inp !!" << std::endl;
		exit(1);
	}

	while(!inputFile.eof())
	{
		std::string line;
		inputFile >> line;

		//if( line.substr(0,16).compare("DIVISION_NUMBERS") == 0 ){
		if( line == "DIVISION_NUMBERS" ){
			readDivisionNumbers(inputFile);
		}
		//else if( line.substr(0,13).compare("X_COORDINATES") == 0 ){
		else if( line == "X_COORDINATES" ){
			readXCoordinates(inputFile);
		}
		//else if( line.substr(0,13).compare("Y_COORDINATES") == 0 ){
		else if( line == "Y_COORDINATES" ){
			readYCoordinates(inputFile);
		}
		//else if( line.substr(0,13).compare("Z_COORDINATES") == 0 ){
		else if( line == "Z_COORDINATES" ){
			readZCoordinates(inputFile);
		}
		//else if( line.substr(0,16).compare("NUMBERING_METHOD") == 0 ){
		else if( line == "NUMBERING_METHOD" ){
			readNumberingMethod(inputFile);
		}
		//else if( line.substr(0,6).compare("LAYERS") == 0 ){
		else if( line == "LAYERS" ){
			readLayers(inputFile);
		}
		//else if( line.substr(0,5).compare("LAYER") == 0 ){
		else if( line == "LAYER" ){
			int layerIDStart(0);
			int layerIDEnd(0);
			inputFile >> layerIDStart;
			inputFile >> layerIDEnd;
			if( layerIDStart == 1  ){
				std::cerr << "Resistivity block division cannot be specified for layer one, in which resistivity value is uniform regardless of the input." << std::endl;
				exit(1);
			}else if( layerIDStart < 2 ){
				std::cerr << "IDs of layers must be larger than 1 . : StartID = " << layerIDStart << ", EndID = " << layerIDEnd << std::endl;
				exit(1);
			}else if( layerIDEnd > m_numLayers ){
				std::cerr << "IDs of layers must less than or equal to total number of layer. StartID = " << layerIDStart << ", EndID = " << layerIDEnd << std::endl;
				exit(1);
			}else{
				readLayer( inputFile, layerIDStart, layerIDEnd );
			}
		}
		//else if( line.substr(0,19).compare("INITIAL_RESISTIVITY") == 0 ){
		else if( line == "INITIAL_RESISTIVITY" ){
			readInitialResistivity(inputFile);
		}
		//else if( line.substr(0,23).compare("AIR_RESISTIVITY") == 0 ){
		else if( line == "AIR_RESISTIVITY" ){
			readAirResistivity(inputFile);
		}
		//else if( line.substr(0,9).compare("ANOMALIES") == 0 ){
		else if( line == "ANOMALIES" ){
			readAnomalies(inputFile);
		}
		//else if( line.substr(0,4).compare("HEXA") == 0 ){
		else if( line == "HEXA" ){
			m_meshType = HEXA;
		}
		//else if( line.substr(0,5).compare("TETRA") == 0 ){
		else if( line == "TETRA" ){
			m_meshType = TETRA;
		}
		//else if( line.substr(0,6).compare("TETRA2") == 0 ){
		else if( line == "TETRA2" ){
			m_meshType = TETRA2;
		}
		//else if( line.substr(0,3).compare("END") == 0 ){
		else if( line == "END" ){
			std::cout << "End of the data." << std::endl;
			break;
		}
		else{
			std::cerr << "Improper data : " << line << std::endl;
			exit(1);
		}

	}
	inputFile.close();

	if( m_meshType == HEXA ){
		std::cout << "Type of mesh : Hexahedron" << std::endl; 
	}else if( m_meshType == TETRA ){
		std::cout << "Type of mesh : Tetrahedron" << std::endl; 
	}else if( m_meshType == TETRA2 ){
		std::cout << "Type of mesh : Tetrahedron2" << std::endl; 
	}

	// Check whether data of every layer has been read
	for( int i = 1; i < m_numLayers; ++i ){
		if( m_hasLayerRead[i] == false ){
			std::cerr << "Data of layer " << i+1 << " has not been set." << std::endl;	
			exit(1);
		}
	}

	// Total numbers of elements
	m_numElemTotal = m_numX * m_numY * m_numZ; 

	// Total numbers of nodes
	m_numNodeTotal = ( m_numX + 1 ) * ( m_numY + 1 ) * ( m_numZ + 1 ); 
}

// Read division numbers.
void MeshGeneration::readDivisionNumbers(std::ifstream& infile){

	if( m_hasDivisionNumberRead == true ){
		std::cerr << "Division numbers has already been read." << std::endl;	
		exit(1);
	}
	else{
		std::cout << "Read division numbers." << std::endl;	
	}

	infile >> m_numX;
	infile >> m_numY;
	infile >> m_numZ;

	std::cout << m_numX << " " << m_numY << " " << m_numZ << std::endl;

	m_hasDivisionNumberRead = true;

}

//// Read edge lengths
//void MeshGeneration::readEdgeLengths(std::ifstream& infile){
//
//	if( m_hasEdgeLengthsRead == true ){
//		std::cerr << "Edge lengths has already been read." << std::endl;	
//		exit(1);
//	}
//	else{
//		std::cout << "Read edge lengths." << std::endl;
//	}
//
//	m_CoordinatesX = new double[ m_numX + 1 ];// lenghs of edges parallel to the X direction.
//	m_CoordinatesY = new double[ m_numY + 1 ];// lenghs of edges parallel to the Y direction.
//	m_CoordinatesZ = new double[ m_numZ + 1 ];// lenghs of edges parallel to the Z direction.
//
//	m_CoordinatesX[0] = 0.0;
//	for( int i = 1; i < m_numX + 1; ++i ){
//		infile >> m_CoordinatesX[i];
//		std::cout << m_CoordinatesX[i] << std::endl;
//	}
//
//	m_CoordinatesY[0] = 0.0;
//	for( int i = 1; i < m_numY + 1; ++i ){
//		infile >> m_CoordinatesY[i];
//		std::cout << m_CoordinatesY[i] << std::endl;
//	}
//
//	m_CoordinatesZ[0] = 0.0;
//	for( int i = 1; i < m_numZ + 1; ++i ){
//		infile >> m_CoordinatesZ[i];
//		std::cout << m_CoordinatesZ[i] << std::endl;
//	}
//
//	m_hasEdgeLengthsRead = true;
//
//};

// Read X coordinates
void MeshGeneration::readXCoordinates(std::ifstream& infile){

	if( m_hasXCoordinatesRead == true ){
		std::cerr << "X coordinates have already been read." << std::endl;	
		exit(1);
	}
	else{
		std::cout << "Read X coordinates." << std::endl;
	}

	m_CoordinatesX = new double[ m_numX + 1 ];// X coorinates

	for( int i = 0; i < m_numX + 1; ++i ){
		infile >> m_CoordinatesX[i];
		std::cout << m_CoordinatesX[i] << " ";
	}
	std::cout << std::endl;

	m_hasXCoordinatesRead = true;

};

// Read Y coordinates
void MeshGeneration::readYCoordinates(std::ifstream& infile){

	if( m_hasYCoordinatesRead == true ){
		std::cerr << "Y coordinates have already been read." << std::endl;	
		exit(1);
	}
	else{
		std::cout << "Read Y coordinates." << std::endl;
	}

	m_CoordinatesY = new double[ m_numY + 1 ];// Y coorinates

	for( int i = 0; i < m_numY + 1; ++i ){
		infile >> m_CoordinatesY[i];
		std::cout << m_CoordinatesY[i] << " " ;
	}
	std::cout << std::endl;

	m_hasYCoordinatesRead = true;

};

// Read Z coordinates
void MeshGeneration::readZCoordinates(std::ifstream& infile){

	if( m_hasZCoordinatesRead == true ){
		std::cerr << "Z coordinates have already been read." << std::endl;	
		exit(1);
	}
	else{
		std::cout << "Read Z coordinates." << std::endl;
	}

	m_CoordinatesZ = new double[ m_numZ + 1 ];// Z coorinates

	for( int i = 0; i < m_numZ + 1; ++i ){
		infile >> m_CoordinatesZ[i];
		std::cout << m_CoordinatesZ[i] << " ";
	}
	std::cout << std::endl;

	m_hasZCoordinatesRead = true;

};

// Read numbering method
void MeshGeneration::readNumberingMethod(std::ifstream& infile){

	if( m_hasNumberingMethodRead == true ){
		std::cerr << "Numbering method has already been read." << std::endl;	
		exit(1);
	}
	else{
		std::cout << "Read numbering method." << std::endl;	
	}

	infile >> m_orderOfNumbering;

	if( m_orderOfNumbering != 1 ){
		std::cerr << "Paramter specifing order of numbering must be one !!" << std::endl;
		exit(1);
	}

	std::cout << m_orderOfNumbering << std::endl;

	if( m_orderOfNumbering < 1 || m_orderOfNumbering > 3 ){
		std::cerr << "Wrong input for paramter specifing order of numbering : " << m_orderOfNumbering << std::endl;
		exit(1);
	}

	m_hasNumberingMethodRead = true;
}

//// Read data of layers
//void MeshGeneration::readEarthLayers(std::ifstream& infile){
//
//	if( m_hasEarthLayersRead == true ){
//		std::cerr << "Data of layers of the earth has already been read." << std::endl;	
//		exit(1);
//	}
//	else{
//		std::cout << "Read data of layers of the earth." << std::endl;	
//	}
//
//	infile >> m_numEarthLayers;
//
//	std::cout << m_numEarthLayers << std::endl;
//
//	m_elemGroupingZ = new int[ m_numEarthLayers + 1 ];
//	m_elemGroupingZ[0] = 0;
//	for( int i = 1; i < m_numEarthLayers + 1; ++i ){
//		infile >> m_elemGroupingZ[i];
//		std::cout << m_elemGroupingZ[i] << std::endl;
//	}
//
//	m_numElemGroupX = new int[m_numEarthLayers];
//	m_numElemGroupY = new int[m_numEarthLayers];
//	m_elemGroupingX = new int*[m_numEarthLayers];
//	m_elemGroupingY = new int*[m_numEarthLayers];
//
//	m_hasEarthLayersRead = true;
//}
//
//// Read number of air layers
//void MeshGeneration::readAirLayer(std::ifstream& infile){
//
//	if( m_hasAirLayerRead == true ){
//		std::cerr << "Number of air layers has already been read." << std::endl;	
//		exit(1);
//	}
//	else{
//		std::cout << "Read number of air layers." << std::endl;	
//	}
//
//	infile >> m_numElemOfAir;
//
//	std::cout << m_numElemOfAir << std::endl;
//
//	m_hasAirLayerRead = true;
//}

// Read data of layers
void MeshGeneration::readLayers(std::ifstream& infile){

	if( m_hasLayersRead == true ){
		std::cerr << "Data of layers of the earth has already been read." << std::endl;	
		exit(1);
	}
	else{
		std::cout << "Read data of layers of the earth." << std::endl;	
	}

	infile >> m_numLayers;

	std::cout << m_numLayers << " ";

	m_hasLayerRead = new bool[ m_numLayers ];
	for( int i = 0; i < m_numLayers; ++i ){
		m_hasLayerRead[i] = false;
	}

	m_elemGroupingZ = new int[ m_numLayers + 1 ];
	m_elemGroupingZ[0] = 0;
	for( int i = 1; i < m_numLayers + 1; ++i ){
		infile >> m_elemGroupingZ[i];
		std::cout << m_elemGroupingZ[i] << " " ;
	}
	std::cout << std::endl;

	m_numElemGroupX = new int[m_numLayers];
	m_numElemGroupY = new int[m_numLayers];
	m_elemGroupingX = new int*[m_numLayers];
	m_elemGroupingY = new int*[m_numLayers];

	// The air layer
	m_numElemGroupX[0] = 1;
	m_elemGroupingX[0] = new int[2];
	m_elemGroupingX[0][0] = 0;
	m_elemGroupingX[0][1] = m_numX;

	m_numElemGroupY[0] = 1;
	m_elemGroupingY[0] = new int[2];
	m_elemGroupingY[0][0] = 0;
	m_elemGroupingY[0][1] = m_numY;
	m_hasLayerRead[0] = true;

	//Data check
	if( m_elemGroupingZ[m_numLayers] != m_numZ ){
		std::cerr << "The Last number of the data of layers must be equal to total division number of Z direction." << std::endl;
		exit(1);
	}

	m_hasLayersRead = true;
}

// Read data of each layer of the earth
void MeshGeneration::readLayer(std::ifstream& infile, const int layerIDStart, const int layerIDEnd){

	if( layerIDStart < 2 || layerIDEnd > m_numLayers ){
		std::cerr << "IDs of layers must be and more than 1 and less than or equal to "
			<< m_numLayers << ". : StartID = " << layerIDStart << ", EndID = " << layerIDEnd << std::endl;
		exit(1);
	}

	const int layerIDZeroStart = layerIDStart - 1;
	const int layerIDZeroEnd = layerIDEnd - 1;

	std::cout << "Read data of the layers from " << layerIDStart << " to " << layerIDEnd << "." << std::endl;
	for( int i = layerIDZeroStart; i < layerIDZeroEnd; ++i ){
		if( m_hasLayerRead[i] == true ){
			std::cerr << "Data of layer " << i+1 << " has already been read." << std::endl;
			exit(1);
		}
	}

	// Read data of X direction
	int nGroupX(0);
	infile >> nGroupX;
	int* elemGroupingXWork = NULL; 
	if( nGroupX < 0 ){

		nGroupX = m_numX;
		elemGroupingXWork = new int[nGroupX+1];

		elemGroupingXWork[0] = 0;
		for( int i = 1; i < nGroupX + 1; ++i ){
			elemGroupingXWork[i] = i;
		}

	}else{

		elemGroupingXWork = new int[nGroupX+1];

		elemGroupingXWork[0] = 0;
		for( int i = 1; i < nGroupX + 1; ++i ){
			infile >> elemGroupingXWork[i];
		}

	}
	// Read data of Y direction
	int nGroupY(0);
	infile >> nGroupY;
	int* elemGroupingYWork = NULL; 
	if( nGroupY < 0 ){

		nGroupY = m_numY;
		elemGroupingYWork = new int[nGroupY+1];

		elemGroupingYWork[0] = 0;
		for( int i = 1; i < nGroupY + 1; ++i ){
			elemGroupingYWork[i] = i;
		}

	}else{

		elemGroupingYWork = new int[nGroupY+1];

		elemGroupingYWork[0] = 0;
		for( int i = 1; i < nGroupY + 1; ++i ){
			infile >> elemGroupingYWork[i];
		}

	}

	for( int iLayer = layerIDZeroStart; iLayer <= layerIDZeroEnd; ++iLayer ){

		std::cout << "Layer# " << iLayer+1 << std::endl;

		// X direction
		std::cout << "X direction" << std::endl;
		std::cout << nGroupX << " ";

		m_numElemGroupX[iLayer] = nGroupX;
		m_elemGroupingX[iLayer] = new int[nGroupX+1];

		m_elemGroupingX[iLayer][0] = 0;
		for( int i = 1; i < m_numElemGroupX[iLayer] + 1; ++i ){
			m_elemGroupingX[iLayer][i] = elemGroupingXWork[i]; 
		}

		for( int i = 1; i < m_numElemGroupX[iLayer] + 1; ++i ){
			std::cout << m_elemGroupingX[iLayer][i] << " ";
		}

		std::cout << std::endl;

		// Y direction
		std::cout << "Y direction" << std::endl;
		std::cout << nGroupY << " ";

		m_numElemGroupY[iLayer] = nGroupY;
		m_elemGroupingY[iLayer] = new int[nGroupY+1];

		m_elemGroupingY[iLayer][0] = 0;
		for( int i = 1; i < m_numElemGroupY[iLayer] + 1; ++i ){
			m_elemGroupingY[iLayer][i] = elemGroupingYWork[i]; 
		}

		for( int i = 1; i < m_numElemGroupY[iLayer] + 1; ++i ){
			std::cout << m_elemGroupingY[iLayer][i] << " ";
		}

		std::cout << std::endl;

		m_hasLayerRead[iLayer] = true;

	}

	delete [] elemGroupingXWork;
	delete [] elemGroupingYWork;

}

// Read initial resistivity
void MeshGeneration::readInitialResistivity(std::ifstream& infile){

	if( m_hasInitialResistivityRead == true ){
		std::cerr << "Initial resistivity has already been read." << std::endl;	
		exit(1);
	}
	else{
		std::cout << "Read initial resistivity." << std::endl;	
	}

	infile >> m_initialResistivity;

	if( m_initialResistivity < 0 ){
		std::cerr << "Initial resistivity value is less than zero !" << std::endl;
		exit(1);			
	}

	std::cout << m_initialResistivity << std::endl;

	m_hasInitialResistivityRead = true;

}

// Read resistivity of the air
void MeshGeneration::readAirResistivity(std::ifstream& infile){

	if( m_hasAirResistivityRead == true ){
		std::cerr << "Air resistivity has already been read." << std::endl;	
		exit(1);
	}
	else{
		std::cout << "Read resistivity of the air." << std::endl;	
	}

	infile >> m_airResistivity;

	//if( m_airResistivity < 0 ){
	//	std::cerr << "Resistivity value of the air is less than zero !" << std::endl;
	//	exit(1);			
	//}
	
	std::cout << m_airResistivity << std::endl;

	m_hasAirResistivityRead = true;

}

// Read data of resistivity anomalies
void MeshGeneration::readAnomalies(std::ifstream& infile){

	if( m_hasAnomaliesRead == true ){
		std::cerr << "Data of resistivity anomalies has already been read." << std::endl;	
		exit(1);
	}
	else{
		std::cout << "Read data of resistivity anomalies." << std::endl;	
	}

	infile >> m_numResisivityAnomalies;

	std::cout << m_numResisivityAnomalies << std::endl;

	m_anomalyData = new struct DataOfAnomaly[m_numResisivityAnomalies];

	for( int i = 0; i < m_numResisivityAnomalies; ++i ){

		infile >> m_anomalyData[i].XStart; 
		if( m_anomalyData[i].XStart < m_CoordinatesX[0] - EPS || 
			m_anomalyData[i].XStart > m_CoordinatesX[m_numX] + EPS ){
			std::cerr << "Lower limit of X coordinate of anomaly " << i << " is out of the range from " 
			          << m_CoordinatesX[0] << " to " << m_CoordinatesX[m_numX]
			          << "." << std::endl;
			exit(1);			
		}

		infile >> m_anomalyData[i].XEnd; 
		if( m_anomalyData[i].XEnd < m_CoordinatesX[0] - EPS || 
			m_anomalyData[i].XEnd > m_CoordinatesX[m_numX] + EPS ){
			std::cerr << "Upper limit of X coordinate of anomaly " << i << " is out of the range from " 
			          << m_CoordinatesX[0] << " to " << m_CoordinatesX[m_numX]
			          << "." << std::endl;
			exit(1);			
		}
		if( m_anomalyData[i].XEnd <= m_anomalyData[i].XStart ){
			std::cerr << "Upper limit of X coordinate of anomaly " << i << " must be larger than lower limit of it." << std::endl;
			exit(1);
		}

		infile >> m_anomalyData[i].YStart; 
		if( m_anomalyData[i].YStart < m_CoordinatesY[0] - EPS || 
			m_anomalyData[i].YStart > m_CoordinatesY[m_numY] + EPS ){
			std::cerr << "Lower limit of Y coordinate of anomaly " << i << " is out of the range from " 
			          << m_CoordinatesY[0] << " to " << m_CoordinatesY[m_numY]
			          << "." << std::endl;
			exit(1);			
		}

		infile >> m_anomalyData[i].YEnd; 
		if( m_anomalyData[i].YEnd < m_CoordinatesY[0] - EPS || 
			m_anomalyData[i].YEnd > m_CoordinatesY[m_numY] + EPS ){
			std::cerr << "Upper limit of Y coordinate of anomaly " << i << " is out of the range from " 
			          << m_CoordinatesY[0] << " to " << m_CoordinatesY[m_numY]
			          << "." << std::endl;
			exit(1);			
		}
		if( m_anomalyData[i].YEnd <= m_anomalyData[i].YStart ){
			std::cerr << "Upper limit of Y coordinate of anomaly " << i << " must be larger than lower limit of it." << std::endl;
			exit(1);
		}

		infile >> m_anomalyData[i].ZStart; 
		if( m_anomalyData[i].ZStart < m_CoordinatesZ[ m_elemGroupingZ[1] ] - EPS || 
			m_anomalyData[i].ZStart > m_CoordinatesZ[ m_numZ ] + EPS ){
			std::cerr << "Lower limit of Z coordinate of anomaly " << i << " is out of the range from " 
			          << m_CoordinatesZ[ m_elemGroupingZ[1] ] << " to " << m_CoordinatesZ[ m_numZ ]
			          << "." << std::endl;
			exit(1);			
		}

		infile >> m_anomalyData[i].ZEnd; 
		if( m_anomalyData[i].ZEnd < m_CoordinatesZ[ m_elemGroupingZ[1] ] - EPS || 
			m_anomalyData[i].ZEnd > m_CoordinatesZ[ m_numZ ] + EPS ){
			std::cerr << "Upper limit of Z coordinate of anomaly " << i << " is out of the range from " 
			          << m_CoordinatesZ[ m_elemGroupingZ[1] ] << " to " << m_CoordinatesZ[ m_numZ ]
			          << "." << std::endl;
			exit(1);			
		}
		else if( m_anomalyData[i].ZEnd <= m_anomalyData[i].ZStart ){
			std::cerr << "Upper limit of Z coordinate of anomaly " << i << " must be larger than lower limit of it." << std::endl;
			exit(1);
		}

		infile >> m_anomalyData[i].resistivityValue; 
		if( m_anomalyData[i].resistivityValue < 0 ){
			std::cerr << "Resistivity value of anomaly " << i << " is less than zero !" << std::endl;
			exit(1);			
		}
	
		int ibuf(0);
		infile >> ibuf;
		m_anomalyData[i].fixResistivityValue = ( ibuf == 1 ? true : false );

		std::cout << m_anomalyData[i].XStart << " " << m_anomalyData[i].XEnd << " "
		          << m_anomalyData[i].YStart << " " << m_anomalyData[i].YEnd << " "
				  << m_anomalyData[i].ZStart << " " << m_anomalyData[i].ZEnd << " "
				  << m_anomalyData[i].resistivityValue;

		if( m_anomalyData[i].fixResistivityValue ){
			std::cout << " Fix"  << std::endl;
		}else{
			std::cout << " Free"  << std::endl;
		}

	}

	m_hasAnomaliesRead = true;
}

// Calculate ID of neighbor element for tetra mesh
int MeshGeneration::calcNeighborElementIDOfTetraMesh( const int iElem, const int iSubElem, const int iFace ) const{

	if( m_meshType != TETRA ){
		std::cerr << "calcNeighborElementIDOfTetraMesh can be use only for tetra mesh !!" << std::endl;
		exit(1);
	}

	if( m_orderOfNumbering != 1 ){
		std::cerr << "Paramter specifing order of numbering must be one !!" << std::endl;
		exit(1);
	}

	switch( iSubElem ){
		case 0:
			switch( iFace ){
				case 0:
					return iElem * 6 + 2;
					break;
				case 1:
					return iElem * 6 + 1;
					break;
				case 2:
					return m_neighborElements[iElem*6    ] < 0 ? -1 : m_neighborElements[iElem*6    ] * 6 + 3;
					break;
				case 3:
					return m_neighborElements[iElem*6 + 4] < 0 ? -1 : m_neighborElements[iElem*6 + 4] * 6 + 4;
					break;
				default:
					std::cerr << "Unknown neighbor face ID !! iFace = " << iFace << std::endl;
					exit(1);
					break;
			}
			break;
		case 1:
			switch( iFace ){
				case 0:
					return iElem * 6; 
					break;
				case 1:
					return iElem * 6 + 3;
					break;
				case 2:
					return m_neighborElements[iElem*6 + 3] < 0 ? -1 : m_neighborElements[iElem*6 + 3] * 6 + 2;
					break;
				case 3:
					return m_neighborElements[iElem*6 + 4] < 0 ? -1 : m_neighborElements[iElem*6 + 4] * 6 + 5;
					break;
				default:
					std::cerr << "Unknown neighbor face ID !! iFace = " << iFace << std::endl;
					exit(1);
					break;
			}
			break;
		case 2:
			switch( iFace ){
				case 0:
					return m_neighborElements[iElem*6    ] < 0 ? -1 : m_neighborElements[iElem*6    ] * 6 + 5;
					break;
				case 1:
					return iElem * 6 + 4;
					break;
				case 2:
					return iElem * 6 + 0;
					break;
				case 3:
					return m_neighborElements[iElem*6 + 2] < 0 ? -1 : m_neighborElements[iElem*6 + 2] * 6 + 1;
					break;
				default:
					std::cerr << "Unknown neighbor face ID !! iFace = " << iFace << std::endl;
					exit(1);
					break;
			}
			break;
		case 3:
			switch( iFace ){
				case 0:
					return m_neighborElements[iElem*6 + 3] < 0 ? -1 : m_neighborElements[iElem*6 + 3] * 6 + 4;
					break;
				case 1:
					return iElem * 6 + 5;
					break;
				case 2:
					return m_neighborElements[iElem*6 + 1] < 0 ? -1 : m_neighborElements[iElem*6 + 1] * 6;
					break;
				case 3:
					return iElem * 6 + 1;
					break;
				default:
					std::cerr << "Unknown neighbor face ID !! iFace = " << iFace << std::endl;
					exit(1);
					break;
			}
			break;
		case 4:
			switch( iFace ){
				case 0:
					return iElem * 6 + 5;
					break;
				case 1:
					return m_neighborElements[iElem*6 + 2] < 0 ? -1 : m_neighborElements[iElem*6 + 2] * 6 + 3;
					break;
				case 2:
					return iElem * 6 + 2;
					break;
				case 3:
					return m_neighborElements[iElem*6 + 5] < 0 ? -1 : m_neighborElements[iElem*6 + 5] * 6;
					break;
				default:
					std::cerr << "Unknown neighbor face ID !! iFace = " << iFace << std::endl;
					exit(1);
					break;
			}
			break;
		case 5:
			switch( iFace ){
				case 0:
					return m_neighborElements[iElem*6 + 1] < 0 ? -1 : m_neighborElements[iElem*6 + 1] * 6 + 2;
					break;
				case 1:
					return iElem * 6 + 4;
					break;
				case 2:
					return iElem * 6 + 3;
					break;
				case 3:
					return m_neighborElements[iElem*6 + 5] < 0 ? -1 : m_neighborElements[iElem*6 + 5] * 6 + 1;
					break;
				default:
					std::cerr << "Unknown neighbor face ID !! iFace = " << iFace << std::endl;
					exit(1);
					break;
			}
			break;
		default:
			std::cerr << "Unknown sub-element ID!! iSubElem = " << iSubElem << std::endl;
			exit(1);
			break;
	}

	return -1;
	
}

// Calculate ID of neighbor element for tetra mesh
int MeshGeneration::calcNodeIDOfTetraMesh( const int iElem, const int iSubElem, const int iNode ) const{

	if( m_meshType != TETRA ){
		std::cerr << "calcNodeIDOfTetraMesh can be use only for tetra mesh !!" << std::endl;
		exit(1);
	}

	if( m_orderOfNumbering != 1 ){
		std::cerr << "Paramter specifing order of numbering must be one !!" << std::endl;
		exit(1);
	}

	switch( iSubElem ){
		case 0:
			switch( iNode ){
				case 0:
					return m_nodesOfElements[iElem*8 + 3];
					break;
				case 1:
					return m_nodesOfElements[iElem*8    ];
					break;
				case 2:
					return m_nodesOfElements[iElem*8 + 1];
					break;
				case 3:
					return m_nodesOfElements[iElem*8 + 7];
					break;
				default:
					std::cerr << "Unknown neighbor node ID !! iNode = " << iNode << std::endl;
					exit(1);
					break;
			}
			break;
		case 1:
			switch( iNode ){
				case 0:
					return m_nodesOfElements[iElem*8 + 2];
					break;
				case 1:
					return m_nodesOfElements[iElem*8 + 3];
					break;
				case 2:
					return m_nodesOfElements[iElem*8 + 1];
					break;
				case 3:
					return m_nodesOfElements[iElem*8 + 7];
					break;
				default:
					std::cerr << "Unknown neighbor node ID !! iNode = " << iNode << std::endl;
					exit(1);
					break;
			}
			break;
		case 2:
			switch( iNode ){
				case 0:
					return m_nodesOfElements[iElem*8 + 1];
					break;
				case 1:
					return m_nodesOfElements[iElem*8    ];
					break;
				case 2:
					return m_nodesOfElements[iElem*8 + 4];
					break;
				case 3:
					return m_nodesOfElements[iElem*8 + 7];
					break;
				default:
					std::cerr << "Unknown neighbor node ID !! iNode = " << iNode << std::endl;
					exit(1);
					break;
			}
			break;
		case 3:
			switch( iNode ){
				case 0:
					return m_nodesOfElements[iElem*8 + 1];
					break;
				case 1:
					return m_nodesOfElements[iElem*8 + 2];
					break;
				case 2:
					return m_nodesOfElements[iElem*8 + 7];
					break;
				case 3:
					return m_nodesOfElements[iElem*8 + 6];
					break;
				default:
					std::cerr << "Unknown neighbor node ID !! iNode = " << iNode << std::endl;
					exit(1);
					break;
			}
			break;
		case 4:
			switch( iNode ){
				case 0:
					return m_nodesOfElements[iElem*8 + 4];
					break;
				case 1:
					return m_nodesOfElements[iElem*8 + 7];
					break;
				case 2:
					return m_nodesOfElements[iElem*8 + 5];
					break;
				case 3:
					return m_nodesOfElements[iElem*8 + 1];
					break;
				default:
					std::cerr << "Unknown neighbor node ID !! iNode = " << iNode << std::endl;
					exit(1);
					break;
			}
			break;
		case 5:
			switch( iNode ){
				case 0:
					return m_nodesOfElements[iElem*8 + 7];
					break;
				case 1:
					return m_nodesOfElements[iElem*8 + 6];
					break;
				case 2:
					return m_nodesOfElements[iElem*8 + 5];
					break;
				case 3:
					return m_nodesOfElements[iElem*8 + 1];
					break;
				default:
					std::cerr << "Unknown neighbor node ID !! iNode = " << iNode << std::endl;
					exit(1);
					break;
			}
			break;
		default:
			std::cerr << "Unknown sub-element ID!! iSubElem = " << iSubElem << std::endl;
			exit(1);
			break;
	}

	return -1;

}

// Calculate type of sub tetra mesh
int MeshGeneration::typeOfSubTetraMesh2( const int ix, const int iy, const int iz ) const{

	if( m_meshType != TETRA2 ){
		std::cerr << "typeOfSubTetraMesh can be use only for tetra2 mesh !!" << std::endl;
		exit(1);
	}

	if( m_orderOfNumbering != 1 ){
		std::cerr << "Paramter specifing order of numbering must be one !!" << std::endl;
		exit(1);
	}

	if( ( ix + iy + iz ) % 2 == 0 ){
		return 0;
	}else{
		return 1;
	}

}

// Calculate ID of neighbor element for tetra mesh2
int MeshGeneration::calcNeighborElementIDOfTetraMesh2( const int ix, const int iy, const int iz, const int iSubElem, const int iFace ) const{

	if( m_meshType != TETRA2 ){
		std::cerr << "typeOfSubTetraMesh can be use only for tetra2 mesh !!" << std::endl;
		exit(1);
	}

	if( m_orderOfNumbering != 1 ){
		std::cerr << "Paramter specifing order of numbering must be one !!" << std::endl;
		exit(1);
	}

	const int type = typeOfSubTetraMesh2( ix, iy, iz );

	const int offset = m_XYZ2ElemID[ix][iy][iz] * 5;

	int offsetXMinus = -9;
	int offsetXPlus  = -9;
	int offsetYMinus = -9;
	int offsetYPlus  = -9;
	int offsetZMinus = -9;
	int offsetZPlus  = -9;

	// Outer boundary
	if( ix == 0 ){
		if( type == 0 ){
			if( ( iSubElem == 0 && iFace == 3 ) || ( iSubElem == 2 && iFace == 2 ) ){
				return -1;
			}
		}else{
			if( ( iSubElem == 0 && iFace == 0 ) || ( iSubElem == 3 && iFace == 3 ) ){
				return -1;
			}
		}
	}else{
		offsetXMinus = m_XYZ2ElemID[ix-1][iy][iz] * 5;
	}
	
	if( ix == m_numX - 1 ){
		if( type == 0 ){
			if( ( iSubElem == 1 && iFace == 0 ) || ( iSubElem == 3 && iFace == 2 ) ){
				return -1;
			}
		}else{
			if( ( iSubElem == 1 && iFace == 1 ) || ( iSubElem == 2 && iFace == 3 ) ){
				return -1;
			}
		}
	}else{
		offsetXPlus = m_XYZ2ElemID[ix+1][iy][iz] * 5;
	}

	if( iy == 0 ){
		if( type == 0 ){
			if( ( iSubElem == 1 && iFace == 1 ) || ( iSubElem == 2 && iFace == 3 ) ){
				return -1;
			}
		}else{
			if( ( iSubElem == 0 && iFace == 3 ) || ( iSubElem == 2 && iFace == 2 ) ){
				return -1;
			}
		}
	}else{
		offsetYMinus = m_XYZ2ElemID[ix][iy-1][iz] * 5; 
	}
	
	if( iy == m_numY - 1 ){
		if( type == 0 ){
			if( ( iSubElem == 0 && iFace == 0 ) || ( iSubElem == 3 && iFace == 3 ) ){
				return -1;
			}
		}else{
			if( ( iSubElem == 1 && iFace == 0 ) || ( iSubElem == 3 && iFace == 2 ) ){
				return -1;
			}
		}
	}else{
		offsetYPlus = m_XYZ2ElemID[ix][iy+1][iz] * 5; 
	}

	if( iz == 0 ){
		if( ( iSubElem == 0 && iFace == 2 ) || ( iSubElem == 1 && iFace == 2 ) ){
			return -1;
		}
	}else{
		offsetZMinus = m_XYZ2ElemID[ix][iy][iz-1] * 5; 
	}
	
	if( iz == m_numZ - 1 ){
		if( ( iSubElem == 2 && iFace == 0 ) || ( iSubElem == 3 && iFace == 0 ) ){
			return -1;
		}
	}else{
		offsetZPlus = m_XYZ2ElemID[ix][iy][iz+1] * 5; 
	}
	
	switch( iSubElem ){
		case 0:
			switch( iFace ){
				case 0:
					if( type == 0 ){
						return offsetYPlus;
					}else{
						return offsetXMinus + 1;
					}
					break;
				case 1:
					return offset + 4;
					break;
				case 2:
					if( type == 0 ){
						return offsetZMinus + 3;
					}else{
						return offsetZMinus + 2;
					}
					break;
				case 3:
					if( type == 0 ){
						return offsetXMinus + 1;
					}else{
						return offsetYMinus;
					}
					break;
				default:
					std::cerr << "Wrong face ID !! : iFace = " << iFace << std::endl;
					exit(1);
					break;
			}
			break;
		case 1:
			switch( iFace ){
				case 0:
					if( type == 0 ){
						return offsetXPlus;
					}else{
						return offsetYPlus + 1;
					}
					break;
				case 1:
					if( type == 0 ){
						return offsetYMinus + 1;
					}else{
						return offsetXPlus;
					}
					break;
				case 2:
					if( type == 0 ){
						return offsetZMinus + 2;
					}else{
						return offsetZMinus + 3;
					}
					break;
				case 3:
					return offset + 4;
					break;
				default:
					std::cerr << "Wrong face ID !! : iFace = " << iFace << std::endl;
					exit(1);
					break;
			}
			break;
		case 2:
			switch( iFace ){
				case 0:
					if( type == 0 ){
						return offsetZPlus;
					}else{
						return offsetZPlus + 1;
					}
					break;
				case 1:
					return offset + 4;
					break;
				case 2:
					if( type == 0 ){
						return offsetXMinus + 2;
					}else{
						return offsetYMinus + 3;
					}
					break;
				case 3:
					if( type == 0 ){
						return offsetYMinus + 3;
					}else{
						return offsetXPlus + 2;
					}
					break;
				default:
					std::cerr << "Wrong face ID !! : iFace = " << iFace << std::endl;
					exit(1);
					break;
			}
			break;
		case 3:
			switch( iFace ){
				case 0:
					if( type == 0 ){
						return offsetZPlus + 1;
					}else{
						return offsetZPlus;
					}
					break;
				case 1:
					return offset + 4;
					break;
				case 2:
					if( type == 0 ){
						return offsetXPlus + 3;
					}else{
						return offsetYPlus + 2;
					}
					break;
				case 3:
					if( type == 0 ){
						return offsetYPlus + 2;
					}else{
						return offsetXMinus + 3;
					}
					break;
				default:
					std::cerr << "Wrong face ID !! : iFace = " << iFace << std::endl;
					exit(1);
					break;
			}
			break;
		case 4:
			switch( iFace ){
				case 0:
					return offset + 3;
					break;
				case 1:
					return offset + 2;
					break;
				case 2:
					return offset + 1;
					break;
				case 3:
					return offset + 0;
					break;
				default:
					std::cerr << "Wrong face ID !! : iFace = " << iFace << std::endl;
					exit(1);
					break;
			}
			break;
		default:
			std::cerr << "Wrong sub element ID !! : iSubElem = " << iSubElem << std::endl;
			exit(1);
			break;
	}

	return -9;

}

// Calculate ID of neighbor element for tetra mesh2
int MeshGeneration::calcNodeIDOfTetraMesh2( const int ix, const int iy, const int iz, const int iSubElem, const int iNode ) const{

	if( m_meshType != TETRA2 ){
		std::cerr << "typeOfSubTetraMesh can be use only for tetra2 mesh !!" << std::endl;
		exit(1);
	}

	if( m_orderOfNumbering != 1 ){
		std::cerr << "Paramter specifing order of numbering must be one !!" << std::endl;
		exit(1);
	}

	const int iElem = m_XYZ2ElemID[ix][iy][iz];

	const int type = typeOfSubTetraMesh2( ix, iy, iz );

	switch( iSubElem ){
		case 0:
			switch( iNode ){
				case 0:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8    ];
					}else{
						return m_nodesOfElements[iElem*8 + 1];
					}
					break;
				case 1:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 3];
					}else{
						return m_nodesOfElements[iElem*8    ];
					}
					break;
				case 2:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 7];
					}else{
						return m_nodesOfElements[iElem*8 + 4];
					}
					break;
				case 3:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 2];
					}else{
						return m_nodesOfElements[iElem*8 + 3];
					}
					break;
				default:
					std::cerr << "Wrong node ID !! : iNode = " << iNode << std::endl;
					exit(1);
					break;
			}
			break;
		case 1:
			switch( iNode ){
				case 0:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8    ];
					}else{
						return m_nodesOfElements[iElem*8 + 1];
					}
					break;
				case 1:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 2];
					}else{
						return m_nodesOfElements[iElem*8 + 3];
					}
					break;
				case 2:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 5];
					}else{
						return m_nodesOfElements[iElem*8 + 6];
					}
					break;
				case 3:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 1];
					}else{
						return m_nodesOfElements[iElem*8 + 2];
					}
					break;
				default:
					std::cerr << "Wrong node ID !! : iNode = " << iNode << std::endl;
					exit(1);
					break;
			}
			break;
		case 2:
			switch( iNode ){
				case 0:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8    ];
					}else{
						return m_nodesOfElements[iElem*8 + 1];
					}
					break;
				case 1:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 4];
					}else{
						return m_nodesOfElements[iElem*8 + 5];
					}
					break;
				case 2:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 5];
					}else{
						return m_nodesOfElements[iElem*8 + 6];
					}
					break;
				case 3:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 7];
					}else{
						return m_nodesOfElements[iElem*8 + 4];
					}
					break;
				default:
					std::cerr << "Wrong node ID !! : iNode = " << iNode << std::endl;
					exit(1);
					break;
			}
			break;
		case 3:
			switch( iNode ){
				case 0:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 2];
					}else{
						return m_nodesOfElements[iElem*8 + 3];
					}
					break;
				case 1:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 6];
					}else{
						return m_nodesOfElements[iElem*8 + 7];
					}
					break;
				case 2:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 7];
					}else{
						return m_nodesOfElements[iElem*8 + 4];
					}
					break;
				case 3:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 5];
					}else{
						return m_nodesOfElements[iElem*8 + 6];
					}
					break;
				default:
					std::cerr << "Wrong node ID !! : iNode = " << iNode << std::endl;
					exit(1);
					break;
			}
			break;
		case 4:
			switch( iNode ){
				case 0:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8    ];
					}else{
						return m_nodesOfElements[iElem*8 + 1];
					}
					break;
				case 1:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 2];
					}else{
						return m_nodesOfElements[iElem*8 + 3];
					}
					break;
				case 2:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 7];
					}else{
						return m_nodesOfElements[iElem*8 + 4];
					}
					break;
				case 3:
					if( type == 0 ){
						return m_nodesOfElements[iElem*8 + 5];
					}else{
						return m_nodesOfElements[iElem*8 + 6];
					}
					break;
				default:
					std::cerr << "Wrong node ID !! : iNode = " << iNode << std::endl;
					exit(1);
					break;
			}
			break;
		default:
			std::cerr << "Wrong sub element ID !! : iSubElem = " << iSubElem << std::endl;
			exit(1);
			break;
	}


}

// Caluculate mesh data
void MeshGeneration::calcMeshData(){

	std::cout << "Generating mesh data." << std::endl;

	int num1(0);
	int num2(0);
	int num3(0);
	switch( m_orderOfNumbering ){
		case XYZ:
			std::cout << "Numbering elements. X => Y => Z " << std::endl;
			num1 = m_numX;
			num2 = m_numY;
			num3 = m_numZ;
			break;
		case YZX:
			std::cout << "Numbering elements. Y => Z => X " << std::endl;
			num1 = m_numY;
			num2 = m_numZ;
			num3 = m_numX;
			break;
		case ZXY:
			std::cout << "Numbering elements. Z => X => Y" << std::endl;
			num1 = m_numZ;
			num2 = m_numX;
			num3 = m_numY;
			break;
		default:
			std::cerr << "Wrong input for paramter specifing order of numbering : " << m_orderOfNumbering << std::endl;
			exit(1);
			break;
	}

	const int numNode1 = num1 + 1;
	const int numNode2 = num2 + 1;
	const int numNode3 = num3 + 1;

	m_locXYZ = new double[ m_numNodeTotal * 3 ]; // Location of nodes

	int iNode(0);
	for( int i3 = 0 ; i3 < numNode3; ++i3 ){
		for( int i2 = 0 ; i2 < numNode2; ++i2 ){
			for( int i1 = 0 ; i1 < numNode1; ++i1 ){

				int iX(0);
				int iY(0);
				int iZ(0);
				switch( m_orderOfNumbering ){
					case XYZ:
						iX = i1;
						iY = i2;
						iZ = i3;
						break;
					case YZX:
						iX = i3;
						iY = i1;
						iZ = i2;
						break;
					case ZXY:
						iX = i2;
						iY = i3;
						iZ = i1;
						break;
				}

				// Location of nodes
				m_locXYZ[ iNode*3     ] = m_CoordinatesX[iX];
				m_locXYZ[ iNode*3 + 1 ] = m_CoordinatesY[iY];
				m_locXYZ[ iNode*3 + 2 ] = m_CoordinatesZ[iZ];

				++iNode;
			}
		}
	}

	// Caluculate accumulated numbers of resisivity blocks 
	calcNumResisivityBlockAccumulated();

	// Number of nodes of each plane
	const int numNodeOfPlane = numNode1 * numNode2;

	m_neighborElements = new int[ m_numElemTotal * 6 ]; // IDs of neighbor elements
	m_nodesOfElements = new int[ m_numElemTotal * 8 ];// IDs of nodes belonging to each element
	m_resistivityBlockID = new int[ m_numElemTotal ];// IDs of resistivity blocks
	m_XYZ2ElemID = new int**[ m_numX ];// Three dimensional array containing element ID
	for( int ix = 0 ; ix < m_numX; ++ix ){
		m_XYZ2ElemID[ix] = new int*[ m_numY ];
		for( int iy = 0 ; iy < m_numY; ++iy ){
			m_XYZ2ElemID[ix][iy] = new int[ m_numZ ];
			for( int iz = 0 ; iz < m_numZ; ++iz ){
				m_XYZ2ElemID[ix][iy][iz] = 0;// Initialization
			}
		}
	}

	// Caluculate neighbor elements and locations of center points.
	int iElem(0);
	for( int i3 = 0 ; i3 < num3; ++i3 ){
		for( int i2 = 0 ; i2 < num2; ++i2 ){
			for( int i1 = 0 ; i1 < num1; ++i1 ){

				// IDs of neighbor elements
				m_neighborElements[iElem*6    ] = iElem - 1;
				m_neighborElements[iElem*6 + 1] = iElem + 1;
				m_neighborElements[iElem*6 + 2] = iElem - num1;
				m_neighborElements[iElem*6 + 3] = iElem + num1;
				m_neighborElements[iElem*6 + 4] = iElem - num1 * num2;
				m_neighborElements[iElem*6 + 5] = iElem + num1 * num2;

				// IDs of neighbor elements are setted to be -1 at boundaries.
				if( i1 == 0 ) {
					m_neighborElements[iElem*6    ] = -1;
				}
				if( i1 == num1 - 1 ) {
					m_neighborElements[iElem*6 + 1] = -1;
				}
				if( i2 == 0 ) {
					m_neighborElements[iElem*6 + 2] = -1;
				}
				if( i2 == num2 - 1 ) {
					m_neighborElements[iElem*6 + 3] = -1;
				}
				if( i3 == 0 ) {
					m_neighborElements[iElem*6 + 4] = -1;
				}
				if( i3 == num3 - 1 ) {
					m_neighborElements[iElem*6 + 5] = -1;
				}
				
				int iX(0);
				int iY(0);
				int iZ(0);
				switch( m_orderOfNumbering ){
					case XYZ:
						iX = i1;
						iY = i2;
						iZ = i3;
						break;
					case YZX:
						iX = i3;
						iY = i1;
						iZ = i2;
						break;
					case ZXY:
						iX = i2;
						iY = i3;
						iZ = i1;
						break;
				}

				//// Edge length of Element
				//m_CoordinatesXOfElements[iElem] = m_CoordinatesX[iX+1] - m_CoordinatesX[iX];
				//m_CoordinatesYOfElements[iElem] = m_CoordinatesY[iY+1] - m_CoordinatesY[iY];
				//m_CoordinatesZOfElements[iElem] = m_CoordinatesZ[iZ+1] - m_CoordinatesZ[iZ];

				//// Location of center of Element
				//m_xOfCenterPoint[iElem] = m_CoordinatesXOfElements[iElem] / 2.0 + m_CoordinatesX[iX];
				//m_yOfCenterPoint[iElem] = m_CoordinatesYOfElements[iElem] / 2.0 + m_CoordinatesY[iY];
				//m_zOfCenterPoint[iElem] = m_CoordinatesZOfElements[iElem] / 2.0 + m_CoordinatesZ[iZ];

				const int node1 = i3 * numNodeOfPlane + i2 * numNode1 + i1;
				const int node2 = node1            + 1;
				const int node3 = node1 + numNode1 + 1;
				const int node4 = node1 + numNode1;

				// IDs of nodes belonging to each element
				m_nodesOfElements[ iElem*8     ] = node1;
				m_nodesOfElements[ iElem*8 + 1 ] = node2;
				m_nodesOfElements[ iElem*8 + 2 ] = node3;
				m_nodesOfElements[ iElem*8 + 3 ] = node4;
				m_nodesOfElements[ iElem*8 + 4 ] = node1 + numNodeOfPlane;
				m_nodesOfElements[ iElem*8 + 5 ] = node2 + numNodeOfPlane;
				m_nodesOfElements[ iElem*8 + 6 ] = node3 + numNodeOfPlane;
				m_nodesOfElements[ iElem*8 + 7 ] = node4 + numNodeOfPlane;

				// IDs of resistivity blocks
				//m_resistivityBlockID[iElem] = calcResisivityBlockID( iX+1, iY+1, iZ+1 );
				m_resistivityBlockID[iElem] = calcResisivityBlockID( iX, iY, iZ );

				// Three dimensional array containing element ID
				m_XYZ2ElemID[iX][iY][iZ] = iElem;

				++iElem;
			}
		}
	}
	
}

// Caluculate resistivity distribution
void MeshGeneration::calcResisivityDistribution(){

	// Initialize resistivity values
	m_ResistivityValues = new double[m_numResistivityBlocks];
	m_ResistivityValues[0] = m_airResistivity;
	for( int i = 1; i < m_numResistivityBlocks; ++i ){
		m_ResistivityValues[i] = m_initialResistivity;
	}
	
	m_fixResistivityValues = new bool[m_numResistivityBlocks];
	m_fixResistivityValues[0] = true;
	for( int i = 1; i < m_numResistivityBlocks; ++i ){
		m_fixResistivityValues[i] = false;
	}

	for( int iLayer = 1; iLayer < m_numLayers; ++iLayer ){
		for( int ix = 0; ix < m_numElemGroupX[iLayer]; ++ix ){
			for( int iy = 0; iy < m_numElemGroupY[iLayer]; ++iy ){

				const double blkX1 = m_CoordinatesX[ m_elemGroupingX[iLayer][ix]   ];
				const double blkX2 = m_CoordinatesX[ m_elemGroupingX[iLayer][ix+1] ];
				const double blkY1 = m_CoordinatesY[ m_elemGroupingY[iLayer][iy]   ];  
				const double blkY2 = m_CoordinatesY[ m_elemGroupingY[iLayer][iy+1] ];
				const double blkZ1 = m_CoordinatesZ[ m_elemGroupingZ[iLayer  ]     ];
				const double blkZ2 = m_CoordinatesZ[ m_elemGroupingZ[iLayer+1]     ];
				bool includeAnomaly(false);

				for( int iano = 0; iano < m_numResisivityAnomalies; ++iano ){
					const double x1 = m_anomalyData[iano].XStart;
					const double x2 = m_anomalyData[iano].XEnd;
					const double y1 = m_anomalyData[iano].YStart;
					const double y2 = m_anomalyData[iano].YEnd;
					const double z1 = m_anomalyData[iano].ZStart;
					const double z2 = m_anomalyData[iano].ZEnd;

					if( x1 > blkX2 - EPS || x2 < blkX1 + EPS ||
						y1 > blkY2 - EPS || y2 < blkY1 + EPS ||
						z1 > blkZ2 - EPS || z2 < blkZ1 + EPS ){
						continue;
					}

					const int blkID = calcResisivityBlockID( m_elemGroupingX[iLayer][ix], m_elemGroupingY[iLayer][iy], m_elemGroupingZ[iLayer] );
					m_ResistivityValues[blkID] = m_anomalyData[iano].resistivityValue;
					m_fixResistivityValues[blkID] = m_anomalyData[iano].fixResistivityValue;

					std::cout << "Block " << blkID << " includes anomaly " << iano << ". Its resistivity value is set to " <<
						m_anomalyData[iano].resistivityValue << " [Ohm-m] " <<
						": X from " << blkX1 << " [km] to " << blkX2 << " [km] " <<
						", Y from " << blkY1 << " [km] to " << blkY2 << " [km] " <<
						", Z from " << blkZ1 << " [km] to " << blkZ2 << " [km] . "; 
					if( m_fixResistivityValues[blkID] ){
						std::cout << " Its resistivity value is fixed." << std::endl;
					}else{
						std::cout << " Its resistivity value is modifiable." << std::endl;
					}

					if( includeAnomaly ){
						std::cerr << "Block " << blkID << " includes more than one anomalies. Its resistivity value is overwritten with " <<
						m_anomalyData[iano].resistivityValue << " [Ohm-m] " <<
						": X from " << blkX1 << " [km] to " << blkX2 << " [km] " <<
						", Y from " << blkY1 << " [km] to " << blkY2 << " [km] " <<
						", Z from " << blkZ1 << " [km] to " << blkZ2 << " [km] . "; 

						if( m_fixResistivityValues[blkID] ){
							std::cerr << " Its resistivity value is fixed." << std::endl;
						}else{
							std::cerr << " Its resistivity value is modifiable." << std::endl;
						}

					}else{
						includeAnomaly = true;
					}

				}

			}
		}
	}

}

// Caluculate accumulated numbers of resisivity blocks 
void MeshGeneration::calcNumResisivityBlockAccumulated(){
	
	m_numResisivityBlockAccumulated = new int[m_numLayers+1];

	m_numResisivityBlockAccumulated[0] = 0;
	int icount(0);
	for( int iLayer = 0; iLayer < m_numLayers; ++iLayer ){
		icount += m_numElemGroupX[iLayer] * m_numElemGroupY[iLayer];
		m_numResisivityBlockAccumulated[iLayer + 1] = icount;
	}

	m_numResistivityBlocks = m_numResisivityBlockAccumulated[m_numLayers];
}

// Caluculate ID of resisivity block
int MeshGeneration::calcResisivityBlockID( const int ix, const int iy, const int iz ){

	//if( iz <= m_elemGroupingZ[1] ){ // corresponds to the air
	if( iz < m_elemGroupingZ[1] ){ // corresponds to the air
		return 0;
	}

	int iLayer(1);
	for( ; iLayer < m_numLayers + 1; ++iLayer ){
		//if( iz > m_elemGroupingZ[iLayer] && iz <= m_elemGroupingZ[iLayer+1] ){
		if( iz >= m_elemGroupingZ[iLayer] && iz < m_elemGroupingZ[iLayer+1] ){
			break;
		}
	}

	int iGx(0);	
	for( ; iGx < m_numElemGroupX[iLayer]; ++iGx ){
		//if( ix > m_elemGroupingX[iLayer][iGx] && ix <= m_elemGroupingX[iLayer][iGx+1] ){
		if( ix >= m_elemGroupingX[iLayer][iGx] && ix < m_elemGroupingX[iLayer][iGx+1] ){
			break;
		}
	}

	int iGy(0);	
	for( ; iGy < m_numElemGroupY[iLayer]; ++iGy ){
		//if( iy > m_elemGroupingY[iLayer][iGy] && iy <= m_elemGroupingY[iLayer][iGy+1] ){
		if( iy >= m_elemGroupingY[iLayer][iGy] && iy < m_elemGroupingY[iLayer][iGy+1] ){
			break;
		}
	}

	//int	blkID = m_numResisivityBlockAccumulated[iLayer] + m_numElemGroupY[iLayer]*iGy + iGx;
	int	blkID = m_numResisivityBlockAccumulated[iLayer] + m_numElemGroupX[iLayer]*iGy + iGx;

	return blkID;
}


// Output mesh data and resistivity model data
void MeshGeneration::outputMeshData(){

	std::cout << "Outputing mesh data." << std::endl;

	// Post processing
	FILE *fp;

	//------------------------
	//--- Output mesh data ---
	//------------------------
	//if( (fp = fopen("model_0.dat", "w")) == NULL ) {
	if( (fp = fopen("mesh.dat", "w")) == NULL ) {
		printf("File open error : model_0.dat !! \n");
		exit(1);
	}

	if( m_meshType == HEXA ){// Hexa mesh

		fprintf(fp, "%s\n", "HEXA" );

		fprintf(fp, "%10d%10d%10d%10d\n",m_numX, m_numY, m_numZ, m_elemGroupingZ[1] );
	
		for( int iNode = 0; iNode < m_numNodeTotal; ++iNode ){
			fprintf(fp, "%10d%20f%20f%20f\n",
				iNode,
				m_locXYZ[iNode*3    ] * distanceConversion,			
				m_locXYZ[iNode*3 + 1] * distanceConversion,
				m_locXYZ[iNode*3 + 2] * distanceConversion );
		}

		for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){
			//fprintf(fp, "%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n",
			fprintf(fp, "%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n",
				iElem,
				m_neighborElements[iElem*6    ],
				m_neighborElements[iElem*6 + 1],
				m_neighborElements[iElem*6 + 2],
				m_neighborElements[iElem*6 + 3],
				m_neighborElements[iElem*6 + 4],
				m_neighborElements[iElem*6 + 5],
				m_nodesOfElements[iElem*8    ],
				m_nodesOfElements[iElem*8 + 1],
				m_nodesOfElements[iElem*8 + 2],
				m_nodesOfElements[iElem*8 + 3],
				m_nodesOfElements[iElem*8 + 4],
				m_nodesOfElements[iElem*8 + 5],
				m_nodesOfElements[iElem*8 + 6],
				m_nodesOfElements[iElem*8 + 7] );
		}

		//--------------------------------------------------
		// Elements and nodes belonging to boundary planes
		//--------------------------------------------------
		// Y-Z Plane ( Minus Side )
		fprintf(fp, "%10d\n", m_numY * m_numZ );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int iy = 0; iy < m_numY; ++iy ){
				const int iElem = m_XYZ2ElemID[0][iy][iz];
				int node1 = 0;
				int node2 = 0;
				int node3 = 0;
				int node4 = 0;
				switch( m_orderOfNumbering ){
					case XYZ:
						node1 = m_nodesOfElements[ iElem*8 + 7 ];
						node2 = m_nodesOfElements[ iElem*8 + 4 ];
						node3 = m_nodesOfElements[ iElem*8     ];
						node4 = m_nodesOfElements[ iElem*8 + 3 ];
						break;
					case YZX:
						node1 = m_nodesOfElements[ iElem*8 + 2 ];
						node2 = m_nodesOfElements[ iElem*8 + 3 ];
						node3 = m_nodesOfElements[ iElem*8     ];
						node4 = m_nodesOfElements[ iElem*8 + 1 ];
						break;
					case ZXY:					
						node1 = m_nodesOfElements[ iElem*8 + 5 ];
						node2 = m_nodesOfElements[ iElem*8 + 1 ];
						node3 = m_nodesOfElements[ iElem*8     ];
						node4 = m_nodesOfElements[ iElem*8 + 4 ];
						break;
				}
				fprintf(fp, "%10d%10d%10d%10d%10d\n", iElem, node1, node2, node3, node4 );
			}
		}

		// Y-Z Plane ( Plus Side )
		fprintf(fp, "%10d\n", m_numY * m_numZ );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int iy = 0; iy < m_numY; ++iy ){
				const int iElem = m_XYZ2ElemID[m_numX-1][iy][iz];
				int node1 = 0;
				int node2 = 0;
				int node3 = 0;
				int node4 = 0;
				switch( m_orderOfNumbering ){
					case XYZ:
						node1 = m_nodesOfElements[ iElem*8 + 6 ];
						node2 = m_nodesOfElements[ iElem*8 + 5 ];
						node3 = m_nodesOfElements[ iElem*8 + 1 ];
						node4 = m_nodesOfElements[ iElem*8 + 2 ];
						break;
					case YZX:
						node1 = m_nodesOfElements[ iElem*8 + 6 ];
						node2 = m_nodesOfElements[ iElem*8 + 7 ];
						node3 = m_nodesOfElements[ iElem*8 + 4 ];
						node4 = m_nodesOfElements[ iElem*8 + 5 ];
						break;
					case ZXY:
						node1 = m_nodesOfElements[ iElem*8 + 6 ];
						node2 = m_nodesOfElements[ iElem*8 + 2 ];
						node3 = m_nodesOfElements[ iElem*8 + 3 ];
						node4 = m_nodesOfElements[ iElem*8 + 7 ];
						break;
				}
				fprintf(fp, "%10d%10d%10d%10d%10d\n", iElem, node1, node2, node3, node4 );
			}
		}

		// Z-X Plane ( Minus Side )
		fprintf(fp, "%10d\n", m_numZ * m_numX );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int ix = 0; ix < m_numX; ++ix ){
			//for( int ix = m_numX - 1; ix >= 0; --ix ){
				const int iElem = m_XYZ2ElemID[ix][0][iz];
				int node1 = 0;
				int node2 = 0;
				int node3 = 0;
				int node4 = 0;
				switch( m_orderOfNumbering ){
					case XYZ:
						node1 = m_nodesOfElements[ iElem*8 + 5 ];
						node2 = m_nodesOfElements[ iElem*8 + 4 ];
						node3 = m_nodesOfElements[ iElem*8     ];
						node4 = m_nodesOfElements[ iElem*8 + 1 ];
						break;
					case YZX:
						node1 = m_nodesOfElements[ iElem*8 + 7 ];
						node2 = m_nodesOfElements[ iElem*8 + 3 ];
						node3 = m_nodesOfElements[ iElem*8     ];
						node4 = m_nodesOfElements[ iElem*8 + 4 ];
						break;
					case ZXY:
						node1 = m_nodesOfElements[ iElem*8 + 2 ];
						node2 = m_nodesOfElements[ iElem*8 + 1 ];
						node3 = m_nodesOfElements[ iElem*8     ];
						node4 = m_nodesOfElements[ iElem*8 + 3 ];
						break;
				}
				fprintf(fp, "%10d%10d%10d%10d%10d\n", iElem, node1, node2, node3, node4 );
			}
		}

		// Z-X Plane ( Plus Side )
		fprintf(fp, "%10d\n", m_numZ * m_numX );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iElem = m_XYZ2ElemID[ix][m_numY-1][iz];
				int node1 = 0;
				int node2 = 0;
				int node3 = 0;
				int node4 = 0;
				switch( m_orderOfNumbering ){
					case XYZ:
						node1 = m_nodesOfElements[ iElem*8 + 6 ];
						node2 = m_nodesOfElements[ iElem*8 + 7 ];
						node3 = m_nodesOfElements[ iElem*8 + 3 ];
						node4 = m_nodesOfElements[ iElem*8 + 2 ];
						break;
					case YZX:
						node1 = m_nodesOfElements[ iElem*8 + 6 ];
						node2 = m_nodesOfElements[ iElem*8 + 2 ];
						node3 = m_nodesOfElements[ iElem*8 + 1 ];
						node4 = m_nodesOfElements[ iElem*8 + 5 ];
						break;
					case ZXY:												
						node1 = m_nodesOfElements[ iElem*8 + 6 ];
						node2 = m_nodesOfElements[ iElem*8 + 5 ];
						node3 = m_nodesOfElements[ iElem*8 + 4 ];
						node4 = m_nodesOfElements[ iElem*8 + 7 ];
						break;
				}
				fprintf(fp, "%10d%10d%10d%10d%10d\n", iElem, node1, node2, node3, node4 );
			}
		}

		// X-Y Plane ( Minus Side ) => Top Boundary
		fprintf(fp, "%10d\n", m_numX * m_numY );
		for( int iy = 0; iy < m_numY; ++iy ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iElem = m_XYZ2ElemID[ix][iy][0];
				int node1 = 0;
				int node2 = 0;
				int node3 = 0;
				int node4 = 0;
				switch( m_orderOfNumbering ){
					case XYZ:
						node1 = m_nodesOfElements[ iElem*8 + 2 ];
						node2 = m_nodesOfElements[ iElem*8 + 3 ];
						node3 = m_nodesOfElements[ iElem*8     ];
						node4 = m_nodesOfElements[ iElem*8 + 1 ];
						break;
					case YZX:
						node1 = m_nodesOfElements[ iElem*8 + 5 ];					
						node2 = m_nodesOfElements[ iElem*8 + 1 ];
						node3 = m_nodesOfElements[ iElem*8     ];
						node4 = m_nodesOfElements[ iElem*8 + 4 ];
						break;
					case ZXY:
						node1 = m_nodesOfElements[ iElem*8 + 7 ];
						node2 = m_nodesOfElements[ iElem*8 + 4 ];
						node3 = m_nodesOfElements[ iElem*8     ];
						node4 = m_nodesOfElements[ iElem*8 + 3 ];
						break;
				}
				fprintf(fp, "%10d%10d%10d%10d%10d\n", iElem, node1, node2, node3, node4 );
			}
		}

		// X-Y Plane ( Plus Side ) => Bottom Boundary
		fprintf(fp, "%10d\n", m_numX * m_numY );
		for( int iy = 0; iy < m_numY; ++iy ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iElem = m_XYZ2ElemID[ix][iy][m_numZ-1];
				int node1 = 0;
				int node2 = 0;
				int node3 = 0;
				int node4 = 0;
				switch( m_orderOfNumbering ){
					case XYZ:					
						node1 = m_nodesOfElements[ iElem*8 + 6 ];
						node2 = m_nodesOfElements[ iElem*8 + 7 ];
						node3 = m_nodesOfElements[ iElem*8 + 4 ];
						node4 = m_nodesOfElements[ iElem*8 + 5 ];
						break;
					case YZX:
						node1 = m_nodesOfElements[ iElem*8 + 6 ];
						node2 = m_nodesOfElements[ iElem*8 + 2 ];
						node3 = m_nodesOfElements[ iElem*8 + 3 ];
						node4 = m_nodesOfElements[ iElem*8 + 7 ];
						break;
					case ZXY:
						node1 = m_nodesOfElements[ iElem*8 + 6 ];
						node2 = m_nodesOfElements[ iElem*8 + 5 ];
						node3 = m_nodesOfElements[ iElem*8 + 1 ];
						node4 = m_nodesOfElements[ iElem*8 + 2 ];
						break;
				}
				fprintf(fp, "%10d%10d%10d%10d%10d\n", iElem, node1, node2, node3, node4 );
			}
		}

		fclose(fp);

		//------------------------------
		//--- resistivity model data ---
		//------------------------------
		//if( (fp = fopen("resistivity_model0.dat", "w")) == NULL ) {
		if( (fp = fopen("resistivity_block_iter0.dat", "w")) == NULL ) {
			printf("File open error : resistivity_block_iter0.dat !! \n");
			exit(1);
		}

		fprintf(fp, "%10d%10d\n",m_numX*m_numY*m_numZ, m_numResistivityBlocks );

		for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){
			fprintf(fp, "%10d%10d\n", iElem, m_resistivityBlockID[iElem] );
		}

#ifdef _OLD
		fprintf(fp, "%10d%5s%15e%10d\n", 0, "     ", m_ResistivityValues[0] , 1 );
#else
		fprintf(fp, "%10d%5s%15e%15e%15e%15e%10d\n", 0, "     ", m_ResistivityValues[0], 1.0e-20, 1.0e+20, 1.0, 1 );
#endif
		for( int iBlk = 1; iBlk < m_numResistivityBlocks; ++iBlk ){
			//fprintf(fp, "%10d%5s%15e\n", iBlk, "     ", m_ResistivityValues[iBlk] );
			//fprintf(fp, "%10d%5s%15e%10d\n", iBlk, "     ", m_ResistivityValues[iBlk] , 0 );
#ifdef _OLD
			fprintf(fp, "%10d%5s%15e%10d\n", iBlk, "     ", m_ResistivityValues[iBlk] , m_fixResistivityValues[iBlk] ? 1 : 0 );
#else
			fprintf(fp, "%10d%5s%15e%15e%15e%15e%10d\n", iBlk, "     ", m_ResistivityValues[iBlk], 1.0e-20, 1.0e+20, 1.0, m_fixResistivityValues[iBlk] ? 1 : 0 );
#endif
		}

		fclose(fp);

	}else if( m_meshType == TETRA ){// Tetra mesh

		if( m_orderOfNumbering != 1 ){
			std::cerr << "Paramter specifing order of numbering must be one !!" << std::endl;
			exit(1);
		}

		fprintf(fp, "%s\n", "TETRA" );

		fprintf(fp, "%10d\n", m_numNodeTotal );
	
		for( int iNode = 0; iNode < m_numNodeTotal; ++iNode ){
			fprintf(fp, "%10d%20f%20f%20f\n",
				iNode,
				m_locXYZ[iNode*3    ] * distanceConversion,			
				m_locXYZ[iNode*3 + 1] * distanceConversion,
				m_locXYZ[iNode*3 + 2] * distanceConversion );
		}

		fprintf(fp, "%10d\n", m_numElemTotal * 6 );

		for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){

			for( int isub = 0; isub < 6 ; ++isub ){

				fprintf(fp, "%10d%10d%10d%10d%10d%10d%10d%10d%10d\n",
					iElem * 6 + isub,
					calcNeighborElementIDOfTetraMesh( iElem, isub, 0 ),
					calcNeighborElementIDOfTetraMesh( iElem, isub, 1 ),
					calcNeighborElementIDOfTetraMesh( iElem, isub, 2 ),
					calcNeighborElementIDOfTetraMesh( iElem, isub, 3 ),
					calcNodeIDOfTetraMesh( iElem, isub, 0 ),
					calcNodeIDOfTetraMesh( iElem, isub, 1 ),
					calcNodeIDOfTetraMesh( iElem, isub, 2 ),
					calcNodeIDOfTetraMesh( iElem, isub, 3 ) );

			}

		}

		//--------------------------------------------------
		// Elements and faces belonging to boundary planes
		//--------------------------------------------------
		// Y-Z Plane ( Minus Side )
		fprintf(fp, "%10d\n", m_numY * m_numZ * 2 );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int iy = 0; iy < m_numY; ++iy ){
				const int iElem = m_XYZ2ElemID[0][iy][iz];
				fprintf(fp, "%10d%10d\n", iElem * 6    , 2 );
				fprintf(fp, "%10d%10d\n", iElem * 6 + 2, 0 );
			}
		}

		// Y-Z Plane ( Plus Side )
		fprintf(fp, "%10d\n", m_numY * m_numZ * 2 );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int iy = 0; iy < m_numY; ++iy ){
				const int iElem = m_XYZ2ElemID[m_numX-1][iy][iz];
				fprintf(fp, "%10d%10d\n", iElem * 6 + 3, 2 );
				fprintf(fp, "%10d%10d\n", iElem * 6 + 5, 0 );
			}
		}

		// Z-X Plane ( Minus Side )
		fprintf(fp, "%10d\n", m_numZ * m_numX * 2 );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iElem = m_XYZ2ElemID[ix][0][iz];
				fprintf(fp, "%10d%10d\n", iElem * 6 + 2, 3 );
				fprintf(fp, "%10d%10d\n", iElem * 6 + 4, 1 );
			}
		}

		// Z-X Plane ( Plus Side )
		fprintf(fp, "%10d\n", m_numZ * m_numX * 2 );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iElem = m_XYZ2ElemID[ix][m_numY-1][iz];
				fprintf(fp, "%10d%10d\n", iElem * 6 + 1, 2 );
				fprintf(fp, "%10d%10d\n", iElem * 6 + 3, 0 );
			}
		}

		// X-Y Plane ( Minus Side ) => Top Boundary
		fprintf(fp, "%10d\n", m_numX * m_numY * 2 );
		for( int iy = 0; iy < m_numY; ++iy ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iElem = m_XYZ2ElemID[ix][iy][0];
				fprintf(fp, "%10d%10d\n", iElem * 6    , 3 );
				fprintf(fp, "%10d%10d\n", iElem * 6 + 1, 3 );
			}
		}

		// X-Y Plane ( Plus Side ) => Bottom Boundary
		fprintf(fp, "%10d\n", m_numX * m_numY * 2 );
		for( int iy = 0; iy < m_numY; ++iy ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iElem = m_XYZ2ElemID[ix][iy][m_numZ-1];
				fprintf(fp, "%10d%10d\n", iElem * 6 + 4, 3 );
				fprintf(fp, "%10d%10d\n", iElem * 6 + 5, 3 );
			}
		}

		// Land surface
		fprintf(fp, "%10d\n", m_numX * m_numY * 2 );
		for( int iy = 0; iy < m_numY; ++iy ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iElem = m_XYZ2ElemID[ix][iy][m_elemGroupingZ[1]];
				fprintf(fp, "%10d%10d\n", iElem * 6    , 3 );
				fprintf(fp, "%10d%10d\n", iElem * 6 + 1, 3 );
			}
		}

		fclose(fp);

		//------------------------------
		//--- resistivity model data ---
		//------------------------------
		if( (fp = fopen("resistivity_block_iter0.dat", "w")) == NULL ) {
			printf("File open error : resistivity_block_iter0.dat !! \n");
			exit(1);
		}

		fprintf(fp, "%10d%10d\n",m_numElemTotal * 6, m_numResistivityBlocks );

		for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){
			for( int i = 0; i < 6; ++i ){
				fprintf(fp, "%10d%10d\n", iElem * 6 + i, m_resistivityBlockID[iElem] );
			}
		}

		fprintf(fp, "%10d%5s%15e%10d\n", 0, "     ", m_ResistivityValues[0] , 1 );
		for( int iBlk = 1; iBlk < m_numResistivityBlocks; ++iBlk ){
			fprintf(fp, "%10d%5s%15e%10d\n", iBlk, "     ", m_ResistivityValues[iBlk] , m_fixResistivityValues[iBlk] ? 1 : 0 );
		}

		fclose(fp);

	}else if( m_meshType == TETRA2 ){// Another tetra mesh

		if( m_orderOfNumbering != 1 ){
			std::cerr << "Paramter specifing order of numbering must be one !!" << std::endl;
			exit(1);
		}

		fprintf(fp, "%s\n", "TETRA" );

		fprintf(fp, "%10d\n", m_numNodeTotal );
	
		for( int iNode = 0; iNode < m_numNodeTotal; ++iNode ){
			fprintf(fp, "%10d%20f%20f%20f\n",
				iNode,
				m_locXYZ[iNode*3    ] * distanceConversion,			
				m_locXYZ[iNode*3 + 1] * distanceConversion,
				m_locXYZ[iNode*3 + 2] * distanceConversion );
		}

		// Total numbers of elements
		fprintf(fp, "%10d\n", m_numElemTotal * 5 );

		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int iy = 0; iy < m_numY; ++iy ){
				for( int ix = 0; ix < m_numX; ++ix ){
					const int offset = m_XYZ2ElemID[ix][iy][iz] * 5;
					for( int isub = 0; isub < 5; ++isub ){
						fprintf(fp, "%10d%10d%10d%10d%10d%10d%10d%10d%10d\n",
							offset + isub,
							calcNeighborElementIDOfTetraMesh2( ix, iy, iz, isub, 0 ),
							calcNeighborElementIDOfTetraMesh2( ix, iy, iz, isub, 1 ),
							calcNeighborElementIDOfTetraMesh2( ix, iy, iz, isub, 2 ),
							calcNeighborElementIDOfTetraMesh2( ix, iy, iz, isub, 3 ),
							calcNodeIDOfTetraMesh2( ix, iy, iz, isub, 0 ),
							calcNodeIDOfTetraMesh2( ix, iy, iz, isub, 1 ),
							calcNodeIDOfTetraMesh2( ix, iy, iz, isub, 2 ),
							calcNodeIDOfTetraMesh2( ix, iy, iz, isub, 3 ) );
					}
				}
			}
		}

		//--------------------------------------------------
		// Elements and faces belonging to boundary planes
		//--------------------------------------------------
		// Y-Z Plane ( Minus Side )
		fprintf(fp, "%10d\n", m_numY * m_numZ * 2 );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int iy = 0; iy < m_numY; ++iy ){
				const int ix = 0;
				const int offset = m_XYZ2ElemID[ix][iy][iz] * 5;
				if( typeOfSubTetraMesh2( ix, iy, iz ) == 0 ){
					fprintf(fp, "%10d%10d\n", offset    , 3 );
					fprintf(fp, "%10d%10d\n", offset + 2, 2 );
				}else{
					fprintf(fp, "%10d%10d\n", offset    , 0 );
					fprintf(fp, "%10d%10d\n", offset + 3, 3 );
				}
			}
		}

		// Y-Z Plane ( Plus Side )
		fprintf(fp, "%10d\n", m_numY * m_numZ * 2 );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int iy = 0; iy < m_numY; ++iy ){
				const int ix = m_numX-1;
				const int offset = m_XYZ2ElemID[ix][iy][iz] * 5;
				if( typeOfSubTetraMesh2( ix, iy, iz ) == 0 ){
					fprintf(fp, "%10d%10d\n", offset + 1, 0 );
					fprintf(fp, "%10d%10d\n", offset + 3, 2 );
				}else{
					fprintf(fp, "%10d%10d\n", offset + 1, 1 );
					fprintf(fp, "%10d%10d\n", offset + 2, 3 );
				}
			}
		}

		// Z-X Plane ( Minus Side )
		fprintf(fp, "%10d\n", m_numZ * m_numX * 2 );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iy = 0;
				const int offset = m_XYZ2ElemID[ix][iy][iz] * 5;
				if( typeOfSubTetraMesh2( ix, iy, iz ) == 0 ){
					fprintf(fp, "%10d%10d\n", offset + 1, 1 );
					fprintf(fp, "%10d%10d\n", offset + 2, 3 );
				}else{
					fprintf(fp, "%10d%10d\n", offset    , 3 );
					fprintf(fp, "%10d%10d\n", offset + 2, 2 );
				}
			}
		}

		// Z-X Plane ( Plus Side )
		fprintf(fp, "%10d\n", m_numZ * m_numX * 2 );
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iy = m_numY-1;
				const int offset = m_XYZ2ElemID[ix][iy][iz] * 5;
				if( typeOfSubTetraMesh2( ix, iy, iz ) == 0 ){
					fprintf(fp, "%10d%10d\n", offset    , 0 );
					fprintf(fp, "%10d%10d\n", offset + 3, 3 );
				}else{
					fprintf(fp, "%10d%10d\n", offset + 1, 0 );
					fprintf(fp, "%10d%10d\n", offset + 3, 2 );
				}
			}
		}

		// X-Y Plane ( Minus Side ) => Top Boundary
		fprintf(fp, "%10d\n", m_numX * m_numY * 2 );
		for( int iy = 0; iy < m_numY; ++iy ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iz = 0;
				const int offset = m_XYZ2ElemID[ix][iy][iz] * 5;
				fprintf(fp, "%10d%10d\n", offset    , 2 );
				fprintf(fp, "%10d%10d\n", offset + 1, 2 );
			}
		}

		// X-Y Plane ( Plus Side ) => Bottom Boundary
		fprintf(fp, "%10d\n", m_numX * m_numY * 2 );
		for( int iy = 0; iy < m_numY; ++iy ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iz = m_numZ-1;
				const int offset = m_XYZ2ElemID[ix][iy][iz] * 5;
				fprintf(fp, "%10d%10d\n", offset + 2, 0 );
				fprintf(fp, "%10d%10d\n", offset + 3, 0 );
			}
		}

		// Land surface
		fprintf(fp, "%10d\n", m_numX * m_numY * 2 );
		for( int iy = 0; iy < m_numY; ++iy ){
			for( int ix = 0; ix < m_numX; ++ix ){
				const int iz = m_elemGroupingZ[1];
				const int offset = m_XYZ2ElemID[ix][iy][iz] * 5;
				fprintf(fp, "%10d%10d\n", offset    , 2 );
				fprintf(fp, "%10d%10d\n", offset + 1, 2 );
			}
		}

		fclose(fp);

		//------------------------------
		//--- resistivity model data ---
		//------------------------------
		if( (fp = fopen("resistivity_block_iter0.dat", "w")) == NULL ) {
			printf("File open error : resistivity_block_iter0.dat !! \n");
			exit(1);
		}

		fprintf(fp, "%10d%10d\n",m_numElemTotal * 5, m_numResistivityBlocks );

		for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){
			for( int i = 0; i < 5; ++i ){
				fprintf(fp, "%10d%10d\n", iElem * 5 + i, m_resistivityBlockID[iElem] );
			}
		}

		fprintf(fp, "%10d%5s%15e%10d\n", 0, "     ", m_ResistivityValues[0] , 1 );
		for( int iBlk = 1; iBlk < m_numResistivityBlocks; ++iBlk ){
			fprintf(fp, "%10d%5s%15e%10d\n", iBlk, "     ", m_ResistivityValues[iBlk] , m_fixResistivityValues[iBlk] ? 1 : 0 );
		}

		fclose(fp);
		
	}else{
		std::cerr << "Improper mesh type !" << std::endl;
		exit(1);
	}

}


// Output VTK file
void MeshGeneration::outputVTK(){

	std::cout << "Output mesh data to a VTK file." << std::endl;

	if( m_meshType == TETRA || m_meshType == TETRA2 ){// Tetra mesh
		if( m_orderOfNumbering != 1 ){
			std::cerr << "Paramter specifing order of numbering must be one !!" << std::endl;
			exit(1);
		}
	}

#ifndef _BINARY_OUT_VTK

	// VTK file
	std::ofstream vtkFile( "MeshData.vtk", std::ios::out );

	vtkFile << "# vtk DataFile Version 2.0" << std::endl;
	vtkFile << "MeshData" << std::endl;
	vtkFile << "ASCII" << std::endl;
	vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
	vtkFile << "POINTS " << m_numNodeTotal << " float" << std::endl;

	for( int iNode = 0; iNode < m_numNodeTotal; ++iNode ){
		vtkFile << m_locXYZ[iNode*3  ]*distanceConversion << " "
			    << m_locXYZ[iNode*3+1]*distanceConversion << " "
				<< m_locXYZ[iNode*3+2]*distanceConversion << std::endl;
	}
	

	if( m_meshType == HEXA ){// Hexa mesh

		vtkFile << "CELLS " << m_numElemTotal << " " << m_numElemTotal*9 << std::endl;
		for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){
			vtkFile << 8 << " "
					<< m_nodesOfElements[iElem*8    ] << " " 
					<< m_nodesOfElements[iElem*8 + 1] << " " 
					<< m_nodesOfElements[iElem*8 + 2] << " " 
					<< m_nodesOfElements[iElem*8 + 3] << " " 
					<< m_nodesOfElements[iElem*8 + 4] << " " 
					<< m_nodesOfElements[iElem*8 + 5] << " " 
					<< m_nodesOfElements[iElem*8 + 6] << " " 
					<< m_nodesOfElements[iElem*8 + 7] << std::endl;
		}

		vtkFile << "CELL_TYPES " << m_numElemTotal << std::endl;
		for( int iElem = 0 ; iElem < m_numElemTotal; ++iElem ){
			vtkFile << "12" << std::endl;
		}

		vtkFile << "CELL_DATA " << m_numElemTotal << std::endl;
		vtkFile << "SCALARS BlockID int" <<  std::endl;
		vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < m_numElemTotal; ++iElem ){
			vtkFile << m_resistivityBlockID[iElem] << std::endl;
		}

		vtkFile << "SCALARS Resistivity[Ohm-m] float" <<  std::endl;
		vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < m_numElemTotal; ++iElem ){
			vtkFile << m_ResistivityValues[ m_resistivityBlockID[iElem] ] << std::endl;
		}

		vtkFile << "SCALARS ElemID int" <<  std::endl;
		vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < m_numElemTotal; ++iElem ){
			vtkFile << iElem << std::endl;
		}

	}else if( m_meshType == TETRA ){// Tetra mesh

		const int numElem = m_numElemTotal * 6;

		vtkFile << "CELLS " << numElem << " " << numElem * 5 << std::endl;
		for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){
			for( int isub = 0; isub < 6; ++isub ){
				vtkFile << 4 << " "
						<< calcNodeIDOfTetraMesh( iElem, isub, 0 ) << " " 
						<< calcNodeIDOfTetraMesh( iElem, isub, 1 ) << " " 
						<< calcNodeIDOfTetraMesh( iElem, isub, 2 ) << " " 
						<< calcNodeIDOfTetraMesh( iElem, isub, 3 ) << std::endl;
			}
		}

		vtkFile << "CELL_TYPES " << numElem << std::endl;
		for( int iElem = 0 ; iElem < numElem; ++iElem ){
			vtkFile << "10" << std::endl;
		}

		vtkFile << "CELL_DATA " << numElem << std::endl;
		vtkFile << "SCALARS BlockID int" <<  std::endl;
		vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < m_numElemTotal; ++iElem ){
			for( int isub = 0; isub < 6; ++isub ){
				vtkFile << m_resistivityBlockID[iElem] << std::endl;
			}
		}

		vtkFile << "SCALARS Resistivity[Ohm-m] float" <<  std::endl;
		vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < m_numElemTotal; ++iElem ){
			for( int isub = 0; isub < 6; ++isub ){
				vtkFile << m_ResistivityValues[ m_resistivityBlockID[iElem] ] << std::endl;
			}
		}

		vtkFile << "SCALARS ElemID int" <<  std::endl;
		vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < numElem; ++iElem ){
			vtkFile << iElem << std::endl;
		}

	}else if( m_meshType == TETRA2 ){// Tetra mesh2

		const int numElem = m_numElemTotal * 5;

		vtkFile << "CELLS " << numElem << " " << numElem * 5 << std::endl;
		for( int iz = 0; iz < m_numZ; ++iz ){
			for( int iy = 0; iy < m_numY; ++iy ){
				for( int ix = 0; ix < m_numX; ++ix ){
					for( int isub = 0; isub < 5; ++isub ){

						vtkFile << 4 << " "
								<< calcNodeIDOfTetraMesh2( ix, iy, iz, isub, 0 ) << " " 
								<< calcNodeIDOfTetraMesh2( ix, iy, iz, isub, 1 ) << " " 
								<< calcNodeIDOfTetraMesh2( ix, iy, iz, isub, 2 ) << " " 
								<< calcNodeIDOfTetraMesh2( ix, iy, iz, isub, 3 ) << std::endl;

					}
				}
			}
		}

		vtkFile << "CELL_TYPES " << numElem << std::endl;
		for( int iElem = 0 ; iElem < numElem; ++iElem ){
			vtkFile << "10" << std::endl;
		}

		vtkFile << "CELL_DATA " << numElem << std::endl;
		vtkFile << "SCALARS BlockID int" <<  std::endl;
		vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < m_numElemTotal; ++iElem ){
			for( int isub = 0; isub < 5; ++isub ){
				vtkFile << m_resistivityBlockID[iElem] << std::endl;
			}
		}

		vtkFile << "SCALARS Resistivity[Ohm-m] float" <<  std::endl;
		vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < m_numElemTotal; ++iElem ){
			for( int isub = 0; isub < 5; ++isub ){
				vtkFile << m_ResistivityValues[ m_resistivityBlockID[iElem] ] << std::endl;
			}
		}

		vtkFile << "SCALARS ElemID int" <<  std::endl;
		vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < numElem; ++iElem ){
			vtkFile << iElem << std::endl;
		}

	}else{
		std::cerr << "Improper mesh type !" << std::endl;
		exit(1);
	}

	vtkFile << "POINT_DATA " << m_numNodeTotal << std::endl;
	vtkFile << "SCALARS NodeID int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( int iNode = 0 ; iNode < m_numNodeTotal; ++iNode ){
		vtkFile << iNode << std::endl;
	}

	vtkFile.close();

#else

	// VTK file
	std::ofstream vtkFile( "MeshData.vtk", std::ios::out | std::ios::binary );

	char buf[80];
	sprintf( buf, "# vtk DataFile Version 2.0\n" );
	vtkFile.write( buf, strlen(buf) );
	sprintf( buf, "MeshData\n" );
	vtkFile.write( buf, strlen(buf) );
	sprintf( buf, "BINARY\n" );
	vtkFile.write( buf, strlen(buf) );
	sprintf( buf, "DATASET UNSTRUCTURED_GRID\n" );
	vtkFile.write( buf, strlen(buf) );

	sprintf( buf, "POINTS %d double\n", m_numNodeTotal );
	vtkFile.write( buf, strlen(buf) );
	vtkFile.write( (char*) m_locXYZ, sizeof( double ) * m_numNodeTotal * 3 );

	//sprintf( buf, "CELLS %d %d\n", m_numElemTotal, m_numElemTotal * 9 );
	//vtkFile.write( buf, strlen(buf) );
	//vtkFile.write( (char*) m_nodesOfElements, sizeof( int ) * m_numElemTotal * 9 );

	//sprintf( buf, "CELL_TYPES %d\n", m_numElemTotal );
	//vtkFile.write( buf, strlen(buf) );
	//int* dummy = new int[m_numElemTotal];
	//for( int iElem = 0 ; iElem < m_numElemTotal; ++iElem ){
	//	dummy[iElem] = 12;
	//}
	//vtkFile.write( (char*) dummy, sizeof( double ) * m_numElemTotal );

	//delete [] dummy;

	vtkFile.close();

#endif

}
