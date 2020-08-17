#ifndef DBLDEF_MESH_GENERATION
#define DBLDEF_MESH_GENERATION

#include <fstream>

struct DataOfAnomaly{
	double XStart;
	double XEnd;
	double YStart;
	double YEnd;
	double ZStart;
	double ZEnd;
	double resistivityValue;
	bool fixResistivityValue;
};

class MeshGeneration{
	public:
		// Constructer
		explicit MeshGeneration();

		// Destructer
		~MeshGeneration();
		
		// Read input file
		void readInputFile();

		// Caluculate mesh data
		void calcMeshData();

		// Caluculate resistivity distribution
		void calcResisivityDistribution();

		// Caluculate ID of resisivity block
		int calcResisivityBlockID( const int ix, const int iy, const int iz );

		// Caluculate accumulated numbers of resisivity blocks 
		void calcNumResisivityBlockAccumulated();

		// Output mesh data
		void outputMeshData();

		// Output VTK file
		void outputVTK();

private:
		enum numbering{
			XYZ=1,
			YZX=2,
			ZXY=3
		};

		enum MeshType{
			HEXA = 0,
			TETRA,
			TETRA2
		};

		static double EPS;

		static double distanceConversion;


		// Division number (X direction)
		int m_numX;

		// Division number (Y direction)
		int m_numY;

		// Division number (Z direction)
		int m_numZ;

		// Total numbers of elements
		int m_numElemTotal; 

		// Total numbers of nodes
		int m_numNodeTotal; 

		// X coordinates
		double* m_CoordinatesX;

		// Y coordinates
		double* m_CoordinatesY;

		// Z coordinates
		double* m_CoordinatesZ;

		// Order of numbering
		int m_orderOfNumbering;

		// Total number of layers
		int m_numLayers;

		// Total number of resistivity blocks
		int m_numResistivityBlocks;

		// Element grouping of each layer (Z direction, Accumulated values))
		int* m_elemGroupingZ;

		// Number of element group of each layer (X direction)
		int* m_numElemGroupX;

		// Number of element group of each layer (Y direction)
		int* m_numElemGroupY;

		// Element grouping of each layer (X direction, Accumulated values))
		int** m_elemGroupingX;

		// Element grouping of each layer (Y direction, Accumulated values))
		int** m_elemGroupingY;

		// Array of neighbor elements
		int* m_neighborElements;

		// Location of nodes
		double* m_locXYZ;
		
		// IDs of nodes belonging to each element
		int* m_nodesOfElements;

		// Caluculate accumulated numbers of resisivity blocks 
		int* m_numResisivityBlockAccumulated;

		// ID of resisitivity blocks
		int* m_resistivityBlockID;

		// Three dimensional array containing element ID
		int*** m_XYZ2ElemID;

		// Initial resistivity value
		double m_initialResistivity;

		// Resistivity value of the air
		double m_airResistivity;

		// Resistivity value including anomalies
		double* m_ResistivityValues;

		// Flag specifing fix resistivity value or not
		bool* m_fixResistivityValues;

		// Number of resistivity anomalies
		int m_numResisivityAnomalies;

		struct DataOfAnomaly *m_anomalyData;

		// Whether or not "DIVISION NUMBERS" has already read
		bool m_hasDivisionNumberRead;

		// Whether or not "X COORDINATES" has already read
		bool m_hasXCoordinatesRead;

		// Whether or not "Y COORDINATES" has already read
		bool m_hasYCoordinatesRead;

		// Whether or not "Z COORDINATES" has already read
		bool m_hasZCoordinatesRead;

		// Whether or not "NUMBERING METHOD" has already read
		bool m_hasNumberingMethodRead;

		// Whether or not "LAYERS" has already read
		bool m_hasLayersRead;

		// Whether or not "LAYER" of each layer has already read
		bool* m_hasLayerRead;

		// Whether or not "INITIAL RESISTIVITY" has already read
		bool m_hasInitialResistivityRead;

		// Whether or not "AIR RESISTIVITY" has already read
		bool m_hasAirResistivityRead;
	
		// Whether or not "ANOMALIES" has already read
		bool m_hasAnomaliesRead;

		// Type of mesh
		int m_meshType;

		// Read division numbers from input file
		void readDivisionNumbers(std::ifstream& infile);

		// Read edge lengths
		void readEdgeLengths(std::ifstream& infile);

		// Read X coordinates
		void readXCoordinates(std::ifstream& infile);

		// Read Y coordinates
		void readYCoordinates(std::ifstream& infile);

		// Read Z coordinates
		void readZCoordinates(std::ifstream& infile);

		// Read numbering method
		void readNumberingMethod(std::ifstream& infile);

		// Read data of layers
		void readLayers(std::ifstream& infile);

		// Read data of each layer
		//void readLayer(std::ifstream& infile, const int layerID);
		void readLayer(std::ifstream& infile, const int layerIDStart, const int layerIDEnd);

		// Read initial resistivity
		void readInitialResistivity(std::ifstream& infile);

		// Read resistivity of the air
		void readAirResistivity(std::ifstream& infile);

		// Read data of resistivity anomalies
		void readAnomalies(std::ifstream& infile);

		// Calculate ID of neighbor element for tetra mesh
		int calcNeighborElementIDOfTetraMesh( const int iElem, const int iSubElem, const int iFace ) const;

		// Calculate ID of neighbor element for tetra mesh
		int calcNodeIDOfTetraMesh( const int iElem, const int iSubElem, const int iNode ) const;

		// Calculate type of sub tetra mesh
		int typeOfSubTetraMesh2( const int ix, const int iy, const int iz ) const;

		// Calculate ID of neighbor element for tetra mesh2
		int calcNeighborElementIDOfTetraMesh2( const int ix, const int iy, const int iz, const int iSubElem, const int iFace ) const;

		// Calculate ID of neighbor element for tetra mesh2
		int calcNodeIDOfTetraMesh2(  const int ix, const int iy, const int iz, const int iSubElem, const int iNode ) const;

		//// Calculate elemdnt ID and face ID of boundary plane for tetra mesh2
		//int calcElemAndFacdIDOfBoundaryForTetraMesh2( const int iPlane, const int ix, const int iy, const int iz, const int iSubElem, const int iNode ) const;

};

#endif

