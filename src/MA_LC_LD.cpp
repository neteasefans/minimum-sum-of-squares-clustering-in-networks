#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <time.h>
#include <iostream>
#include<iomanip>
#include<io.h>
using namespace std;


struct Edge {
	bool tabu; 					// true - if edge is tabu
	int firstEndpointIndex;
	int secondEndpointIndex;
	double length;
};

struct Interior {
	int edgeIndex;
	double length; 				// length  belongs to (0;length of current edge); length - distance between firstEndpoint of edge and interior point 
};

struct Vertex {
	bool tabu; 					// true - if vertex is tabu
	int clusterIndex;
	double distanceBetweenVertexAndPrototype;
	double bottleneckPoint;
	double secondMinimumDistance;
	int clusterIndexOfSecondMinimumDistance;
};

struct IndexesArray {
	int * items;
	int size;
};

struct Prototype {
	double clusterObjective;
								// (vertexIndex = -1 and interior.edgeIndex != -1) or (vertexIndex != -1 and interior.edgeIndex = -1) 
	int vertexIndex;
	Interior interior;
	IndexesArray assignedVertices;
};

struct ShortestPathToInteriorPoint {
	IndexesArray pathThroughFirstEndpoint;
	IndexesArray pathThroughSecondEndpoint;
	IndexesArray pathThroughDependOn_x;
	double sum_pathThroughFirstEndpoint;
	double sum_pathThroughSecondEndpoint;
};

struct BestOneClusterSolution {
	int edgeIndex;
	double e;
	double objective;
};

struct Clustering {
	Prototype* prototypeArray;
	double objective;
};

struct Interval {
	double endpoint;
	bool isLeftEndpoint;
	int vertex;
};


struct IntervalArray {
	Interval * items;
	int size;
};

struct Webber_Points
{
	Interior interior;
	int vertex_number;
	bool isVertex;
};
//for tabu search
typedef struct Tabu
{
	int tabuTenure;
	int *tabuList;
}Tabu;
Tabu Tb;
double **dis;				//preserve the minimum and the second minimum distances between each vertex and the prototypes
double objective;
int *cluLen;
int **clu;
int *mark;
int *sol;
int *flag_add;

//for lc_ld
int real_len_s;
Webber_Points *wp_points;
const double MAX_DOUBLE = std::numeric_limits<double>::max();

int verticesNumber;
int edgesNumber;
int clustersNumber;
int counter;
int timeLimit;

Edge* edgeArray = 0;
Vertex* vertexArray = 0;
Clustering bestSolution;
Clustering currentSolution;
Clustering tmpSolution;
Clustering tmpSolution_update;
Clustering *population;
double ** shortestPathMatrix = 0;
int** isEdgeExist = 0;
ShortestPathToInteriorPoint shPath;
double* segmentEnds = 0;
IndexesArray potentialAssignedVert;
IntervalArray intervalArray;
time_t startTime, endTime, time_1;
double best_time_per;
#define NUM_LC_LD 20
#define PopNum 5
#define DEVIATION 0.00001
#define DEBUG
#define NUM1 100
static int compare(const void* a, const void* b) {
	if (*(double *)b < *(double *)a) { return 1; }
	if (*(double *)b > *(double *)a) { return -1; }
	return 0;
}

void exitOnError(const char * error, const char * what)
{
	printf("%s %s \n", error, what);
	exit(1);
}

/*find the lengths of the shortest paths between all pairs of vertices*/
void Floyd_harshall_algorithm()
{
	for (int k = 1; k <= verticesNumber; k++) {
		for (int i = 1; i <= verticesNumber; i++) {
			for (int j = 1; j <= verticesNumber; j++) {
				if (shortestPathMatrix[i][k] != MAX_DOUBLE && shortestPathMatrix[k][j] != MAX_DOUBLE && shortestPathMatrix[i][j] > shortestPathMatrix[i][k] + shortestPathMatrix[k][j]) {
					shortestPathMatrix[i][j] = shortestPathMatrix[i][k] + shortestPathMatrix[k][j];
				}
			}
		}
	}
}

void readInstanceAndInitializeVariables(char name[])
{
	ifstream f(name);
	char buffer[256];
	if (!f.is_open())
	{
		exitOnError("Error opening file: ", name);
	}

	do
	{
		f.getline(buffer, 250);
	} while (buffer[0] == 'c' || buffer[0] == '\0');
	sscanf(buffer, "%d %d %d", &verticesNumber, &edgesNumber, &clustersNumber);

	if (verticesNumber < 1 || edgesNumber < 1 || clustersNumber < 1) {
		exitOnError("Wrong number of vertices or/and number of edges or/and number of clusters", "");
	}
	if (verticesNumber < clustersNumber) {
		exitOnError("Number of clusters should be less or equal to number of vertices", "");
	}
	edgeArray = new Edge[edgesNumber + 1];
	vertexArray = new Vertex[verticesNumber + 1];
	bestSolution.prototypeArray = new Prototype[clustersNumber + 1];
	bestSolution.objective = MAX_DOUBLE;
	currentSolution.prototypeArray = new Prototype[clustersNumber + 1];
	currentSolution.objective = MAX_DOUBLE;
	tmpSolution.prototypeArray = new Prototype[clustersNumber + 1];
	tmpSolution_update.prototypeArray = new Prototype[clustersNumber + 1];
	tmpSolution_update.objective = MAX_DOUBLE;
	tmpSolution.objective = MAX_DOUBLE;
	shortestPathMatrix = new double*[verticesNumber + 1];
	isEdgeExist = new int *[verticesNumber + 1];
	segmentEnds = new double[2 * verticesNumber];
	shPath.pathThroughDependOn_x.items = new int[verticesNumber];
	shPath.pathThroughFirstEndpoint.items = new int[verticesNumber];
	shPath.pathThroughSecondEndpoint.items = new int[verticesNumber];
	potentialAssignedVert.items = new int[verticesNumber];
	intervalArray.items = new Interval[2 * verticesNumber];
	shPath.pathThroughDependOn_x.size = 0;
	shPath.pathThroughFirstEndpoint.size = 0;
	shPath.pathThroughSecondEndpoint.size = 0;
	shPath.sum_pathThroughFirstEndpoint = 0;
	shPath.sum_pathThroughSecondEndpoint = 0;
	potentialAssignedVert.size = 0;
	intervalArray.size = 0;
	int allo_len = verticesNumber + clustersNumber*NUM_LC_LD + 1;
	Tb.tabuList = new int[allo_len];
	dis = new double*[verticesNumber];	
	cluLen = new int[clustersNumber + 2];
	clu = new int*[clustersNumber + 2];
	mark = new int[verticesNumber + 1];
	sol = new int[clustersNumber + 1];
	for (int i = 0; i < verticesNumber; i++)
		dis[i] = new double[2];
	for (int i = 0; i <= clustersNumber; i++)
		clu[i] = new int[verticesNumber];
	wp_points = new Webber_Points[verticesNumber + clustersNumber*NUM_LC_LD];
	flag_add = new int[verticesNumber + clustersNumber*NUM_LC_LD];
	population = new Clustering[PopNum + 1];
	for (int i = 1; i <= PopNum; i++)
		population[i].prototypeArray = new Prototype[clustersNumber + 1];
	for (int i = 1; i <= PopNum; i++)
		for (int j = 1; j <= clustersNumber; j++)
			population[i].prototypeArray[j].assignedVertices.items = new int[verticesNumber];
	for (int i = 1; i <= verticesNumber; i++) {
		shortestPathMatrix[i] = new double[verticesNumber + 1];
		vertexArray[i].tabu = false;
		vertexArray[i].distanceBetweenVertexAndPrototype = MAX_DOUBLE;
		isEdgeExist[i] = new int[verticesNumber + 1];
		for (int j = 1; j <= verticesNumber; j++) {
			isEdgeExist[i][j] = -1;
		}
		if (i <= clustersNumber) {
			bestSolution.prototypeArray[i].interior.edgeIndex = -1;
			bestSolution.prototypeArray[i].interior.length = 0;
			bestSolution.prototypeArray[i].vertexIndex = -1;
			bestSolution.prototypeArray[i].assignedVertices.items = new int[verticesNumber];
			bestSolution.prototypeArray[i].assignedVertices.size = 0;

			currentSolution.prototypeArray[i].interior.edgeIndex = -1;
			currentSolution.prototypeArray[i].interior.length = 0;
			currentSolution.prototypeArray[i].vertexIndex = -1;
			currentSolution.prototypeArray[i].assignedVertices.items = new int[verticesNumber];
			currentSolution.prototypeArray[i].assignedVertices.size = 0;

			tmpSolution.prototypeArray[i].interior.edgeIndex = -1;
			tmpSolution.prototypeArray[i].interior.length = 0;
			tmpSolution.prototypeArray[i].vertexIndex = -1;
			tmpSolution.prototypeArray[i].assignedVertices.items = new int[verticesNumber];
			tmpSolution.prototypeArray[i].assignedVertices.size = 0;

			tmpSolution_update.prototypeArray[i].interior.edgeIndex = -1;
			tmpSolution_update.prototypeArray[i].interior.length = 0;
			tmpSolution_update.prototypeArray[i].vertexIndex = -1;
			tmpSolution_update.prototypeArray[i].assignedVertices.items = new int[verticesNumber];
			tmpSolution_update.prototypeArray[i].assignedVertices.size = 0;

			for (int j = 0; j < verticesNumber; j++) {
				bestSolution.prototypeArray[i].assignedVertices.items[j] = -1;
				currentSolution.prototypeArray[i].assignedVertices.items[j] = -1;
				tmpSolution.prototypeArray[i].assignedVertices.items[j] = -1;
				tmpSolution_update.prototypeArray[i].assignedVertices.items[j] = -1;
			}
		}
		for (int j = 1; j <= verticesNumber; j++) {
			if (i != j) {
				shortestPathMatrix[i][j] = MAX_DOUBLE;
			}
			else {
				shortestPathMatrix[i][j] = 0;
			}
		}
	}

	int edgeIndex = 1;
	while (!f.eof())
	{
		f.getline(buffer, 250);
		if (buffer[0] == 'c' || buffer[0] == '\0') continue;
		int result = sscanf(buffer, "%d %d %lf", &edgeArray[edgeIndex].firstEndpointIndex, &edgeArray[edgeIndex].secondEndpointIndex, &edgeArray[edgeIndex].length);
		if (result != 3) {
			exitOnError("Wrong line: ", buffer);
		}
		if (isEdgeExist[edgeArray[edgeIndex].firstEndpointIndex][edgeArray[edgeIndex].secondEndpointIndex] == -1) {
			isEdgeExist[edgeArray[edgeIndex].firstEndpointIndex][edgeArray[edgeIndex].secondEndpointIndex] = edgeIndex;
			isEdgeExist[edgeArray[edgeIndex].secondEndpointIndex][edgeArray[edgeIndex].firstEndpointIndex] = edgeIndex;
			if (edgeArray[edgeIndex].firstEndpointIndex > edgeArray[edgeIndex].secondEndpointIndex) {
				int tmp = edgeArray[edgeIndex].firstEndpointIndex;
				edgeArray[edgeIndex].firstEndpointIndex = edgeArray[edgeIndex].secondEndpointIndex;
				edgeArray[edgeIndex].secondEndpointIndex = tmp;
			}
			edgeArray[edgeIndex].tabu = false;
			shortestPathMatrix[edgeArray[edgeIndex].firstEndpointIndex][edgeArray[edgeIndex].secondEndpointIndex] = edgeArray[edgeIndex].length;
			shortestPathMatrix[edgeArray[edgeIndex].secondEndpointIndex][edgeArray[edgeIndex].firstEndpointIndex] = edgeArray[edgeIndex].length;
			edgeIndex++;
		}
		else {
			int curEdgeIndex = isEdgeExist[edgeArray[edgeIndex].firstEndpointIndex][edgeArray[edgeIndex].secondEndpointIndex];
			if (edgeArray[edgeIndex].length < edgeArray[curEdgeIndex].length) {
				shortestPathMatrix[edgeArray[curEdgeIndex].firstEndpointIndex][edgeArray[curEdgeIndex].secondEndpointIndex] = edgeArray[edgeIndex].length;
				shortestPathMatrix[edgeArray[curEdgeIndex].secondEndpointIndex][edgeArray[curEdgeIndex].firstEndpointIndex] = edgeArray[edgeIndex].length;
			}
			edgesNumber--;
		}
	}

	if (edgeIndex != (edgesNumber + 1)) {
		exitOnError("Wrong number of edges", "");
	}
	double density = (double)edgesNumber / (verticesNumber * (verticesNumber - 1) / 2);
	printf("|E| = %d |V| = %d density = %6.5f  \n", edgesNumber, verticesNumber, density);
	f.close();
}

std::random_device rd;
std::mt19937_64 generator(rd());
std::uniform_int_distribution<int>* vertexDiscreteDistribution = NULL;
std::uniform_int_distribution<int>* edgeDiscreteDistribution = NULL;
std::uniform_real_distribution<double> continuousDistributionWithout0and1(0.00000000000000000000000000000000001, 1);
std::uniform_real_distribution<double> continuousDistributionWith0and1(0, 1.00000000000000000000000000000000001);
std::uniform_int_distribution<int> coin(0, 1);

void generate_prototype_vertex_with_max_dist(int clusterIndex, Prototype* prototypeArray);


void generate_prototype_inside_edge(int clusterIndex, Prototype* prototypeArray) {
	int currentEdge;
	do {
		currentEdge = (*edgeDiscreteDistribution)(generator);
	} while (edgeArray[currentEdge].tabu || vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu || vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu);

	double x = continuousDistributionWithout0and1(generator);
	prototypeArray[clusterIndex].interior.edgeIndex = currentEdge;
	prototypeArray[clusterIndex].interior.length = x * edgeArray[currentEdge].length;
	prototypeArray[clusterIndex].vertexIndex = -1;
	edgeArray[currentEdge].tabu = true;
	vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu = true;
	vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu = true;
}

void generate_prototype_on_vertex(int clusterIndex, Prototype* prototypeArray) {
	int currentVertex;
	do {
		currentVertex = (*vertexDiscreteDistribution)(generator);
	} while (vertexArray[currentVertex].tabu);

	vertexArray[currentVertex].tabu = true;
	prototypeArray[clusterIndex].vertexIndex = currentVertex;
	prototypeArray[clusterIndex].interior.edgeIndex = -1;
}

void generate_prototype_on_vertex_or_edge_50_50(int clusterIndex, Prototype* prototypeArray) {
	if (coin(generator)) {
		generate_prototype_on_vertex(clusterIndex, prototypeArray);
	}
	else {
		generate_prototype_inside_edge(clusterIndex, prototypeArray);
	}
}

void generate_prototype_on_sorted_edges(Prototype* prototypeArray) { 		// it is possible to make it more effective
	int* indexes = new int[edgesNumber + 1];
	int i, j;
	for (i = 1; i <= edgesNumber; i++) {
		indexes[i] = i;
	}
	for (i = 1; i <= edgesNumber; i++) { 									// sort indexes by increasing edge length
		for (j = i; j <= edgesNumber; j++) {
			if (edgeArray[indexes[i]].length > edgeArray[indexes[j]].length) {
				int tmp;
				tmp = indexes[j];
				indexes[j] = indexes[i];
				indexes[i] = tmp;
			}
		}
	}
	int startIndex = 1;
	int currentEdge = 0;
	for (i = 1; i <= clustersNumber; i++) {
		int counter = 0;
		for (j = startIndex; (edgeArray[indexes[j]].length == edgeArray[indexes[startIndex]].length) && j <= edgesNumber; j++)
		{
			if (!edgeArray[indexes[j]].tabu && !vertexArray[edgeArray[indexes[j]].firstEndpointIndex].tabu && !vertexArray[edgeArray[indexes[j]].secondEndpointIndex].tabu) {
				counter++;
			}
		}
		if (counter == 0) {
			startIndex = j;
			i--;
			continue;
		}
		if (startIndex != (j - 1)) {
			std::uniform_int_distribution<int> discreteDistribution(startIndex, j - 1);
			do {
				currentEdge = indexes[discreteDistribution(generator)];
			} while (edgeArray[currentEdge].tabu || vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu || vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu);
		}
		else {
			currentEdge = indexes[startIndex];
		}

		double x = continuousDistributionWithout0and1(generator);
		prototypeArray[i].interior.edgeIndex = currentEdge;
		prototypeArray[i].interior.length = x * edgeArray[currentEdge].length;
		prototypeArray[i].vertexIndex = -1;
		edgeArray[currentEdge].tabu = true;
		vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu = true;
		vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu = true;
		if (counter == 1) {
			startIndex = j;
		}
	}
	delete [] indexes;
}

void generate_prototype(char generatorName, int clusterIndex, Prototype* prototypeArray) {
	switch (generatorName) {
	case 'v': {
		generate_prototype_on_vertex(clusterIndex, prototypeArray);
		break;
	}
	case 'e': {
		generate_prototype_inside_edge(clusterIndex, prototypeArray);
		break;
	}
	case 'r': {
		generate_prototype_on_vertex_or_edge_50_50(clusterIndex, prototypeArray);
		break;
	}
	default: {
		exitOnError("Wrong generator parameter: ", &generatorName);
	}
	}
}

void generate_prototype_vertex_with_max_dist(int clusterIndex, Prototype* prototypeArray) {
	int indexOfVertexWithMaxDist = 1;
	for (int vertexIndex = 2; vertexIndex <= verticesNumber; vertexIndex++) {
		if (vertexArray[vertexIndex].distanceBetweenVertexAndPrototype > vertexArray[indexOfVertexWithMaxDist].distanceBetweenVertexAndPrototype) {
			indexOfVertexWithMaxDist = vertexIndex;
		}
	}
	prototypeArray[clusterIndex].interior.edgeIndex = -1;
	prototypeArray[clusterIndex].vertexIndex = indexOfVertexWithMaxDist;
	vertexArray[indexOfVertexWithMaxDist].distanceBetweenVertexAndPrototype = 0;
}

void generate_prototype_for_remove_degeneracy(char generatorName, int clusterIndex, Prototype* prototypeArray) {
	if (generatorName == 'm') {
		generate_prototype_vertex_with_max_dist(clusterIndex, prototypeArray);
	}
	else {
		generate_prototype(generatorName, clusterIndex, prototypeArray);
	}
}

/* return true if prototypes are changed*/
bool allocate(Clustering& currentSolution) {
	Prototype* prototypeArray = currentSolution.prototypeArray;
	currentSolution.objective = 0;
	bool functionResult = false;
	double minDistance;
	double currentDistance;
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		prototypeArray[cluster].clusterObjective = 0;
		prototypeArray[cluster].assignedVertices.size = 0;
	}
	for (int currentVertex = 1; currentVertex <= verticesNumber; currentVertex++) {
		minDistance = MAX_DOUBLE;
		int nearestPrototypeIndex = -1;
		// find the nearest prototype
		for (int cluster = 1; cluster <= clustersNumber; cluster++) {
			if (prototypeArray[cluster].vertexIndex != -1) {
				currentDistance = shortestPathMatrix[currentVertex][prototypeArray[cluster].vertexIndex];
			}
			else {
				int edgeIndex = prototypeArray[cluster].interior.edgeIndex;
				double lengthFromFirstendpoint = prototypeArray[cluster].interior.length;

				currentDistance = shortestPathMatrix[currentVertex][edgeArray[edgeIndex].firstEndpointIndex] + lengthFromFirstendpoint;
				double tmpDistance = shortestPathMatrix[currentVertex][edgeArray[edgeIndex].secondEndpointIndex] + edgeArray[edgeIndex].length - lengthFromFirstendpoint;
				if (tmpDistance < currentDistance) {
					currentDistance = tmpDistance;
				}
			}
			if (minDistance > currentDistance) {
				minDistance = currentDistance;
				nearestPrototypeIndex = cluster;
			}
		}
		if (prototypeArray[nearestPrototypeIndex].assignedVertices.items[prototypeArray[nearestPrototypeIndex].assignedVertices.size] != currentVertex)
		{
			functionResult = true;
		}
		prototypeArray[nearestPrototypeIndex].assignedVertices.items[prototypeArray[nearestPrototypeIndex].assignedVertices.size++] = currentVertex;
		vertexArray[currentVertex].clusterIndex = nearestPrototypeIndex;
		vertexArray[currentVertex].distanceBetweenVertexAndPrototype = minDistance;
		prototypeArray[nearestPrototypeIndex].clusterObjective += pow(minDistance, 2);
	}
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		prototypeArray[cluster].assignedVertices.items[prototypeArray[cluster].assignedVertices.size] = -1;
		currentSolution.objective += prototypeArray[cluster].clusterObjective;
	}
	return functionResult;
}

inline double computeDistance(Prototype& currentPrototype, int currentVertex) {
	if (currentPrototype.vertexIndex != -1) {
		return shortestPathMatrix[currentVertex][currentPrototype.vertexIndex];
	}
	else {
		int edgeIndex = currentPrototype.interior.edgeIndex;
		double lengthFromFirstendpoint = currentPrototype.interior.length;
		return min(shortestPathMatrix[currentVertex][edgeArray[edgeIndex].firstEndpointIndex] + lengthFromFirstendpoint,
			shortestPathMatrix[currentVertex][edgeArray[edgeIndex].secondEndpointIndex] + edgeArray[edgeIndex].length - lengthFromFirstendpoint);
	}
}

inline void updateTabuInformation(Prototype& currentPrototype) {
	if (currentPrototype.vertexIndex != -1) {
		vertexArray[currentPrototype.vertexIndex].tabu = true;
	}
	else {
		int edgeIndex = currentPrototype.interior.edgeIndex;
		edgeArray[edgeIndex].tabu = true;
		vertexArray[edgeArray[edgeIndex].firstEndpointIndex].tabu = true;
		vertexArray[edgeArray[edgeIndex].secondEndpointIndex].tabu = true;
	}
}

inline void clearTabuInformation(Prototype& currentPrototype) {
	if (currentPrototype.vertexIndex != -1) {
		vertexArray[currentPrototype.vertexIndex].tabu = false;
	}
	else {
		int edgeIndex = currentPrototype.interior.edgeIndex;
		edgeArray[edgeIndex].tabu = false;
		vertexArray[edgeArray[edgeIndex].firstEndpointIndex].tabu = false;
		vertexArray[edgeArray[edgeIndex].secondEndpointIndex].tabu = false;
	}
}

inline void clearTabuInformationForNotRemaining(bool* remainingPrototypes, Prototype* prototypeArray) {
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		if (!remainingPrototypes[cluster]) {
			clearTabuInformation(prototypeArray[cluster]);
		}
	}
}

inline void clearAllTabuInformation() {
	for (int i = 1; i <= verticesNumber; i++) {
		vertexArray[i].tabu = false;
	}
	for (int i = 1; i <= edgesNumber; i++) {
		edgeArray[i].tabu = false;
	}
}

void shaking(bool* remainingPrototypes, int numberOfDifferencesBetweenSolutions, Prototype* prototypeArray, char generatorName) {
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_int_distribution<int> discreteDistribution(1, clustersNumber);
	int currentNumberOfDiff;
	int currentCluster;
	if (numberOfDifferencesBetweenSolutions > (clustersNumber / 2)) {
		for (int cluster = 1; cluster <= clustersNumber; cluster++) {
			remainingPrototypes[cluster] = false;
		}
		currentNumberOfDiff = clustersNumber;
		while (currentNumberOfDiff > numberOfDifferencesBetweenSolutions) {
			currentCluster = discreteDistribution(generator);
			if (!remainingPrototypes[currentCluster]) {
				remainingPrototypes[currentCluster] = true;
				currentNumberOfDiff--;
			}
		}
	}
	else {
		for (int cluster = 1; cluster <= clustersNumber; cluster++) {
			remainingPrototypes[cluster] = true;
		}
		currentNumberOfDiff = 0;
		while (currentNumberOfDiff < numberOfDifferencesBetweenSolutions) {
			currentCluster = discreteDistribution(generator);
			if (remainingPrototypes[currentCluster]) {
				remainingPrototypes[currentCluster] = false;
				currentNumberOfDiff++;
			}
		}

	}
	clearTabuInformationForNotRemaining(remainingPrototypes, prototypeArray);
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		if (!remainingPrototypes[cluster]) {
			//generate prototype
			generate_prototype(generatorName, cluster, prototypeArray);
		}
	}
}

int remove_degeneracy(bool* remainingPrototypes, Prototype* prototypeArray, char generatorName) {
	int counterDegeneracy = 0;
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		if (prototypeArray[cluster].assignedVertices.size == 0) {
			clearTabuInformation(prototypeArray[cluster]);
			++counterDegeneracy;
		}
	}
	if (counterDegeneracy == 0) { return 0; }
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		if (prototypeArray[cluster].assignedVertices.size == 0) {
			remainingPrototypes[cluster] = false;
			//generate prototype
			generate_prototype_for_remove_degeneracy(generatorName, cluster, prototypeArray);
		}
	}
	return counterDegeneracy;
}

void printEdge(int index) {
	cout << "(" << edgeArray[index].firstEndpointIndex << "; " << edgeArray[index].secondEndpointIndex << ")";
}

void printPrototypes(Prototype* prototypeArray) {
	for (int i = 1; i <= clustersNumber; i++) {
		cout << "Prototype " << i << " ";
		if (prototypeArray[i].vertexIndex != -1) {
			cout << "vertex " << prototypeArray[i].vertexIndex;
		}
		else {
			cout << "edge ";
			printEdge(prototypeArray[i].interior.edgeIndex);
			cout << " x=" << (prototypeArray[i].interior.length / edgeArray[prototypeArray[i].interior.edgeIndex].length);
		}
		cout << '\n';
	}
}

double computeObjectiveFunctionForClusterAccordingToShPath(Prototype& currentPrototype, int edgeIndex, double e) {
	double sum = 0;
	if (e > 0 && e < edgeArray[edgeIndex].length) {
		for (int i = 0; i < shPath.pathThroughFirstEndpoint.size; i++) {
			sum += pow((shortestPathMatrix[shPath.pathThroughFirstEndpoint.items[i]][edgeArray[edgeIndex].firstEndpointIndex] + e), 2);
		}
		for (int i = 0; i < shPath.pathThroughDependOn_x.size; i++) {
			sum += pow((shortestPathMatrix[shPath.pathThroughDependOn_x.items[i]][edgeArray[edgeIndex].firstEndpointIndex] + e), 2);
		}
		e = edgeArray[edgeIndex].length - e;
		for (int i = 0; i < shPath.pathThroughSecondEndpoint.size; i++) {
			sum += pow((shortestPathMatrix[shPath.pathThroughSecondEndpoint.items[i]][edgeArray[edgeIndex].secondEndpointIndex] + e), 2);
		}
	}
	else {
		int vertexPr;
		if (e <= 0) {
			vertexPr = edgeArray[edgeIndex].firstEndpointIndex;
		}
		else {
			vertexPr = edgeArray[edgeIndex].secondEndpointIndex;
		}
		for (int index = 0; index < currentPrototype.assignedVertices.size; index++) {
			sum += pow(shortestPathMatrix[vertexPr][currentPrototype.assignedVertices.items[index]], 2);
		}
	}
	return sum;
}

double computeObjectiveForClusterVertices(IndexesArray& assignedVertices, int edgeIndex, double e) {
	double sum = 0;
	double dist = 0;
	int vertexIndex;
	if (e > 0 && e < edgeArray[edgeIndex].length) {
		for (int index = 0; index < assignedVertices.size; index++) {
			vertexIndex = assignedVertices.items[index];
			dist = min(shortestPathMatrix[vertexIndex][edgeArray[edgeIndex].firstEndpointIndex] + e, shortestPathMatrix[vertexIndex][edgeArray[edgeIndex].secondEndpointIndex] + edgeArray[edgeIndex].length - e);
			sum += pow(dist, 2);
		}
	}
	else {
		if (e <= 0) {
			vertexIndex = edgeArray[edgeIndex].firstEndpointIndex;
		}
		else {
			vertexIndex = edgeArray[edgeIndex].secondEndpointIndex;
		}
		for (int index = 0; index < assignedVertices.size; index++) {
			sum += pow(shortestPathMatrix[vertexIndex][assignedVertices.items[index]], 2);
		}
	}
	return sum;
}

double computeObjectiveForCluster(Prototype& currentPrototype) {
	double sum = 0;
	double dist = 0;
	int vertexIndex;
	if (currentPrototype.interior.edgeIndex != -1) {
		for (int index = 0; index < currentPrototype.assignedVertices.size; index++) {
			vertexIndex = currentPrototype.assignedVertices.items[index];
			dist = min(shortestPathMatrix[vertexIndex][edgeArray[currentPrototype.interior.edgeIndex].firstEndpointIndex] + currentPrototype.interior.length, shortestPathMatrix[vertexIndex][edgeArray[currentPrototype.interior.edgeIndex].secondEndpointIndex] + edgeArray[currentPrototype.interior.edgeIndex].length - currentPrototype.interior.length);
			sum += pow(dist, 2);
		}
	}
	else {
		for (int index = 0; index < currentPrototype.assignedVertices.size; index++) {
			sum += pow(shortestPathMatrix[currentPrototype.vertexIndex][currentPrototype.assignedVertices.items[index]], 2);
		}
	}
	return sum;
}

void computeBottleneckPointsAndSegmentEnds(IndexesArray assignedVertices, int currentEdge, double intervalPoint1, double intervalPoint2, int& sizeOfSegmentEndsArray) {
	shPath.pathThroughDependOn_x.size = 0;
	shPath.pathThroughFirstEndpoint.size = 0;
	shPath.pathThroughSecondEndpoint.size = 0;
	shPath.sum_pathThroughFirstEndpoint = 0;
	shPath.sum_pathThroughSecondEndpoint = 0;
	for (int currentIndex = 0; currentIndex < assignedVertices.size; currentIndex++) {
		int currentVertex = assignedVertices.items[currentIndex];
		vertexArray[currentVertex].bottleneckPoint = 0.5*(edgeArray[currentEdge].length + shortestPathMatrix[currentVertex][edgeArray[currentEdge].secondEndpointIndex] - shortestPathMatrix[currentVertex][edgeArray[currentEdge].firstEndpointIndex]);
		if (vertexArray[currentVertex].bottleneckPoint <= intervalPoint1) {
			shPath.pathThroughSecondEndpoint.items[shPath.pathThroughSecondEndpoint.size++] = currentVertex;
			shPath.sum_pathThroughSecondEndpoint += shortestPathMatrix[currentVertex][edgeArray[currentEdge].secondEndpointIndex];
		}
		else if (vertexArray[currentVertex].bottleneckPoint >= intervalPoint2) {
			shPath.pathThroughFirstEndpoint.items[shPath.pathThroughFirstEndpoint.size++] = currentVertex;
			shPath.sum_pathThroughFirstEndpoint += shortestPathMatrix[currentVertex][edgeArray[currentEdge].firstEndpointIndex];
		}
		else {
			shPath.pathThroughDependOn_x.items[shPath.pathThroughDependOn_x.size++] = currentVertex;
			segmentEnds[sizeOfSegmentEndsArray++] = vertexArray[currentVertex].bottleneckPoint;
			shPath.sum_pathThroughFirstEndpoint += shortestPathMatrix[currentVertex][edgeArray[currentEdge].firstEndpointIndex];
		}
	}
}

bool solve_1_prototype_for_edge(int currentEdge, Prototype& currentPrototype, BestOneClusterSolution& bestSolution, double intervalPoint1, double intervalPoint2) {
	if (intervalPoint1 > intervalPoint2 || intervalPoint1 > edgeArray[currentEdge].length || intervalPoint2 < 0) {
		return true;
	}
	else if (intervalPoint1 == intervalPoint2) {
		double result = computeObjectiveForClusterVertices(currentPrototype.assignedVertices, currentEdge, intervalPoint1);
		if (result < bestSolution.objective) {
			bestSolution.edgeIndex = currentEdge;
			bestSolution.e = intervalPoint1;
			bestSolution.objective = result;
		}
		return true;
	}

	// compute bottleneck points
	segmentEnds[0] = intervalPoint1;
	segmentEnds[1] = intervalPoint2;
	int sizeOfSegmentEndsArray = 2;
	computeBottleneckPointsAndSegmentEnds(currentPrototype.assignedVertices, currentEdge, intervalPoint1, intervalPoint2, sizeOfSegmentEndsArray);

	// sort segmentEnds by incresing order
	qsort(segmentEnds, sizeOfSegmentEndsArray, sizeof(double), compare);
	double e = 0;
	// find all local minimums
	for (int segmentIndex = 0; segmentIndex < (sizeOfSegmentEndsArray - 1); segmentIndex++) {
		if (segmentEnds[segmentIndex] == segmentEnds[segmentIndex + 1]) { continue; }
		if (segmentIndex != 0) {
			for (int index = 0; index < shPath.pathThroughDependOn_x.size; index++) {
				int currentVertex = shPath.pathThroughDependOn_x.items[index];
				if (vertexArray[currentVertex].bottleneckPoint <= segmentEnds[segmentIndex]) {
					shPath.pathThroughSecondEndpoint.items[shPath.pathThroughSecondEndpoint.size++] = currentVertex;
					if (shPath.pathThroughDependOn_x.size > 1) {
						shPath.pathThroughDependOn_x.items[index] = shPath.pathThroughDependOn_x.items[shPath.pathThroughDependOn_x.size - 1];
					}
					--shPath.pathThroughDependOn_x.size;
					shPath.sum_pathThroughFirstEndpoint -= shortestPathMatrix[currentVertex][edgeArray[currentEdge].firstEndpointIndex];
					shPath.sum_pathThroughSecondEndpoint += shortestPathMatrix[currentVertex][edgeArray[currentEdge].secondEndpointIndex];
					index--;
				}
			}
		}
		double leftSum = shPath.sum_pathThroughFirstEndpoint + segmentEnds[segmentIndex] * (shPath.pathThroughFirstEndpoint.size + shPath.pathThroughDependOn_x.size);
		double rightSum = shPath.sum_pathThroughSecondEndpoint + (edgeArray[currentEdge].length - segmentEnds[segmentIndex]) * (shPath.pathThroughSecondEndpoint.size);
		if (leftSum < rightSum) {
			if (segmentEnds[segmentIndex] + (rightSum - leftSum) / (double)currentPrototype.assignedVertices.size > segmentEnds[segmentIndex + 1]) {
				if ((segmentIndex + 2) != sizeOfSegmentEndsArray) {
					continue;
				}
				e = segmentEnds[segmentIndex + 1];
			}
			else { e = segmentEnds[segmentIndex] + (rightSum - leftSum) / (double)currentPrototype.assignedVertices.size; }
		}
		else {
			if (segmentIndex != 0) { continue; }
			e = 0;
		}
		double result = computeObjectiveFunctionForClusterAccordingToShPath(currentPrototype, currentEdge, e);
		if (result < bestSolution.objective) {
			bestSolution.edgeIndex = currentEdge;
			bestSolution.e = e;
			bestSolution.objective = result;
		}
	}
	return true;
}

// change prototype
void changePrototypeToBestOneClusterSolution(Prototype& currentPrototype, BestOneClusterSolution& bestSolution) {
	currentPrototype.clusterObjective = bestSolution.objective;
	if (bestSolution.e > 0 && bestSolution.e < edgeArray[bestSolution.edgeIndex].length) {
		currentPrototype.vertexIndex = -1;
		currentPrototype.interior.edgeIndex = bestSolution.edgeIndex;
		currentPrototype.interior.length = bestSolution.e;
	}
	else if (bestSolution.e == 0) {
		currentPrototype.vertexIndex = edgeArray[bestSolution.edgeIndex].firstEndpointIndex;
		currentPrototype.interior.edgeIndex = -1;
	}
	else {
		currentPrototype.vertexIndex = edgeArray[bestSolution.edgeIndex].secondEndpointIndex;
		currentPrototype.interior.edgeIndex = -1;
	}
}

void solve_1_prototype_for_every_cluster(Clustering& currentSolution) {
	BestOneClusterSolution bestSolution;
	Prototype* prototypeArray = currentSolution.prototypeArray;
	currentSolution.objective = 0;
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		clearTabuInformation(prototypeArray[cluster]);
		if (prototypeArray[cluster].assignedVertices.size > 1) {
			bestSolution.objective = prototypeArray[cluster].clusterObjective;
			int clusterSize = prototypeArray[cluster].assignedVertices.size;
			for (int i = 0; i < clusterSize; i++) {
				for (int j = 0; j < clusterSize; j++) {
					int currentEdge = isEdgeExist[prototypeArray[cluster].assignedVertices.items[i]][prototypeArray[cluster].assignedVertices.items[j]];
					if (currentEdge != -1) {
						solve_1_prototype_for_edge(currentEdge, prototypeArray[cluster], bestSolution, 0.0, edgeArray[currentEdge].length);
					}
				}
			}
			if (bestSolution.objective < prototypeArray[cluster].clusterObjective) {
				changePrototypeToBestOneClusterSolution(prototypeArray[cluster], bestSolution);
			}
		}
		else if (prototypeArray[cluster].assignedVertices.size == 1) {
			prototypeArray[cluster].vertexIndex = prototypeArray[cluster].assignedVertices.items[0];
			prototypeArray[cluster].interior.edgeIndex = -1;
			prototypeArray[cluster].clusterObjective = 0;
		}
		else {
			cout << "Error! solve_1_prototype_for_every_cluster (852 line)\n";
		}
		updateTabuInformation(prototypeArray[cluster]);
		double tmpResult = computeObjectiveForCluster(prototypeArray[cluster]);
		if (abs(tmpResult - prototypeArray[cluster].clusterObjective) > 0.0001) {
			cout << "Error! solve_1_prototype_for_every_cluster (857 line)\n";
		}
		currentSolution.objective += prototypeArray[cluster].clusterObjective;
	}
}

void copyPrototypesAndUpdateTabu(Prototype* from, Prototype* to, bool* remainingPrototypes) {
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		if (!remainingPrototypes[cluster]) {
			to[cluster].vertexIndex = from[cluster].vertexIndex;
			to[cluster].interior.edgeIndex = from[cluster].interior.edgeIndex;
			to[cluster].interior.length = from[cluster].interior.length;
			updateTabuInformation(to[cluster]);
		}
	}
}

void copyAllInformation(Clustering& from, Clustering& to) {
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		to.prototypeArray[cluster].vertexIndex = from.prototypeArray[cluster].vertexIndex;
		to.prototypeArray[cluster].interior.edgeIndex = from.prototypeArray[cluster].interior.edgeIndex;
		to.prototypeArray[cluster].interior.length = from.prototypeArray[cluster].interior.length;
	}
	to.objective = from.objective;
}

void k_means(Clustering& currentSolution, bool* remainingPrototypes, char remove_degeneracyGen) {
	bool isDifferent = true;
	counter = 0;
	do {
		isDifferent = allocate(currentSolution);
		if (isDifferent || counter == 0) {
			if (remove_degeneracy(remainingPrototypes, currentSolution.prototypeArray, remove_degeneracyGen)) {
				isDifferent = true;
				continue;
			}
			solve_1_prototype_for_every_cluster(currentSolution);
			isDifferent = true;
			//cout << "current.objective=" << currentSolution.objective << endl;
		}
		counter++;
		time(&endTime);
	} while (isDifferent && counter < 1 && difftime(endTime, startTime) < timeLimit);
	if (isDifferent) {
		for (int currentVertex = 1; currentVertex <= verticesNumber; currentVertex++) {
			vertexArray[currentVertex].distanceBetweenVertexAndPrototype = computeDistance(currentSolution.prototypeArray[vertexArray[currentVertex].clusterIndex], currentVertex);
		}
	}
}

void printClusters(Clustering& solution)
{
	Prototype* prototypeArray = solution.prototypeArray;
	int* objectsClusters = new int[verticesNumber + 1];
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		for (int vertexIndex = 0; vertexIndex < prototypeArray[cluster].assignedVertices.size; ++vertexIndex)
			objectsClusters[prototypeArray[cluster].assignedVertices.items[vertexIndex]] = cluster;
	}
	if (verticesNumber > 0)
	{
		cout << objectsClusters[1];
		for (int vertexIndex = 2; vertexIndex <= verticesNumber; ++vertexIndex)
			cout << ", " << objectsClusters[vertexIndex];
		cout << '\n';
	}
	delete [] objectsClusters;
	//delete [] prototypeArray;
}

template<typename T>
void bubble_sort(T *arr, int len)
{
	int i, j;
	T temp;
	for (j = 0; j < len - 1; j++)
	{
		for (i = 0; i < len - 1 - j; i++)
		{
			if (arr[i] > arr[i + 1])
			{
				temp = arr[i];
				arr[i] = arr[i + 1];
				arr[i + 1] = temp;
			}
		}
	}
}

//build tje dis, cluLen and clu arrays 
void build_data(double **dis, int *cluLen, int **clu, const int num_p, const int *sol)
{
	int i, loop, index, index_first;
	double dis_min;	
	for (int i = 0; i < verticesNumber; i++)
	{
		loop = 0;
		index_first = -1;
		while (loop < 2)
		{
			dis_min = DBL_MAX;
			for (int j = 0; j < num_p; j++)
			{
				int proto = sol[j];
				if (proto < verticesNumber)
				{
					//cout << "shortpath=" << shortestPathMatrix[i + 1][proto + 1];
					if (shortestPathMatrix[i + 1][proto + 1] < dis_min && j != index_first)		
					{
						dis_min = shortestPathMatrix[i + 1][proto + 1];
						index = j;
					}
				}
				else
				{
					int curEdge = wp_points[proto].interior.edgeIndex;
					double length = wp_points[proto].interior.length;
					double left = shortestPathMatrix[i + 1][edgeArray[curEdge].firstEndpointIndex] + length;
					double right = shortestPathMatrix[i + 1][edgeArray[curEdge].secondEndpointIndex] + edgeArray[curEdge].length - length;
					if (left > right)
						left = right;
					if (left < dis_min && j != index_first)
					{
						dis_min = left;
						index = j;
					}
				}
			}
			if (loop == 0)
				mark[i] = index;
			index_first = index;
			dis[i][loop] = dis_min;
			loop++;
		}
	}	
	for (i = 0; i < num_p; i++)
		cluLen[i] = 0;
	for (i = 0; i < verticesNumber; i++)
		clu[mark[i]][cluLen[mark[i]]++] = i;
}

int drop_func(int *sol, int iter, int *tbList, int*cluLen, int**clu, double **dis, int &num_p, double &objective)
{
	int i, j, k, vec, drop_node, tail, drop_node_real;
	double delta, delta1, delta_dif;
	int num_best, tabu_num_best;
	int best_v[NUM1], tabu_best_v[NUM1];
	double delta_best, tabu_delta_best;
	delta_best = DBL_MAX;
	tabu_delta_best = DBL_MAX;
	num_best = 0;
	tabu_num_best = 0;
	for (i = 0; i < num_p; i++)
	{
		delta = 0;
		delta1 = 0;
		for (j = 0; j < cluLen[i]; j++)
		{
			vec = clu[i][j];
			delta += dis[vec][0] * dis[vec][0];
		}
		for (k = 0; k < cluLen[i]; k++)
		{
			vec = clu[i][k];
			delta1 += dis[vec][1] * dis[vec][1];
		}
		delta_dif = delta1 - delta;
		if (tbList[sol[i]] <= iter)
		{
			if (delta_dif < delta_best)
			{
				best_v[0] = i;
				delta_best = delta_dif;
				num_best = 1;
			}
			else if (delta_dif == delta_best &&num_best < NUM1)
			{
				best_v[num_best] = i;
				num_best++;
			}
		}
		else
		{
			if (delta_dif < tabu_delta_best)
			{
				tabu_best_v[0] = i;
				tabu_delta_best = delta_dif;
				tabu_num_best = 1;
			}
			else if (delta_dif == tabu_delta_best && tabu_num_best < NUM1)
			{
				tabu_best_v[tabu_num_best] = i;
				tabu_num_best++;
			}
		}
	}
	//aspiration criterion
	if ((tabu_num_best>0 && tabu_delta_best < delta_best &&objective + tabu_delta_best < bestSolution.objective &&num_p - 1 == clustersNumber)
		|| (tabu_num_best > 0 && num_best == 0))
	{
		drop_node = tabu_best_v[rand() % tabu_num_best];
		drop_node_real = sol[drop_node];
		num_p--;
		tail = sol[num_p];
		sol[drop_node] = tail;
		objective += tabu_delta_best;
	}
	else
	{
		drop_node = best_v[rand() % num_best];
		drop_node_real = sol[drop_node];
		num_p--;
		tail = sol[num_p];
		sol[drop_node] = tail;
		objective += delta_best;
	}
	//cout << "drop_cluster="<<drop_node<<", drop_node_real=" << drop_node_real << " ";	
	return sol[drop_node];
}

int add_func(int *sol, int iter, int *tbList, double **dis, int &num_p, double &objective)
{
	int i, j, add_node;
	int num_best;
	int best_v[NUM1];
	double delta, delta_best;
	for (i = 0; i < verticesNumber + NUM_LC_LD*clustersNumber; i++)
		flag_add[i] = 0;												//the set of vertices that can execute the add move 
	for (i = 0; i < num_p; i++)
	{
		if (sol[i]<verticesNumber)
			flag_add[sol[i]] = 1;
		else
		{
			int curEdge = wp_points[sol[i]].interior.edgeIndex;
			flag_add[edgeArray[curEdge].firstEndpointIndex - 1] = 1;
			flag_add[edgeArray[curEdge].secondEndpointIndex - 1] = 1;
			flag_add[sol[i]] = 1;
		}
	}
	delta_best = DBL_MAX;
	num_best = 0;
	for (i = 0; i < real_len_s; i++)
	{
		delta = 0;
		if (flag_add[i] == 0 && tbList[i] <= iter)
		{
			for (j = 0; j < verticesNumber; j++)
			{
				if (i < verticesNumber)
				{
					if (shortestPathMatrix[j + 1][i + 1] < dis[j][0])
						delta += shortestPathMatrix[j + 1][i + 1] * shortestPathMatrix[j + 1][i + 1] - dis[j][0] * dis[j][0];
				}
				else
				{
					int curEdge = wp_points[i].interior.edgeIndex;
					double length = wp_points[i].interior.length;
					double left = shortestPathMatrix[j + 1][edgeArray[curEdge].firstEndpointIndex] + length;
					double right = shortestPathMatrix[j + 1][edgeArray[curEdge].secondEndpointIndex] + edgeArray[curEdge].length - length;
					if (left > right)
						left = right;
					if (left < dis[j][0])
						delta += left*left - dis[j][0] * dis[j][0];
				}
			}
			if (delta < delta_best)
			{
				best_v[0] = i;
				delta_best = delta;
				num_best = 1;
				//cout << "i=" << i << endl;
			}
			else if (delta == delta_best &&num_best < NUM1)
			{
				best_v[num_best] = i;
				num_best++;
			}
		}
	}
	//getchar();
	add_node = best_v[rand() % num_best];
	sol[num_p] = add_node;
	num_p++;
	objective += delta_best;
	tbList[add_node] = iter + rand() % 10 + Tb.tabuTenure;
	return add_node;
}

void copy_information(int *sol, double objective, Clustering &to)
{
	for (int i = 0; i < clustersNumber; i++)
	{
		if (sol[i] < verticesNumber)
		{
			to.prototypeArray[i + 1].vertexIndex = sol[i] + 1;
			to.prototypeArray[i + 1].interior.edgeIndex = -1;
		}
		else
		{
			to.prototypeArray[i + 1].vertexIndex = -1;
			to.prototypeArray[i + 1].interior.edgeIndex = wp_points[sol[i]].interior.edgeIndex;
			to.prototypeArray[i + 1].interior.length = wp_points[sol[i]].interior.length;
		}
	}
	to.objective = objective;
}

double tabu_search(Clustering &currentSol)
{
	int i, iter, dial, num_p, wp_isExist, wp_number;
	int drop_node, add_node;
	double local_minima;
	local_minima = currentSol.objective;
	if (real_len_s >= verticesNumber + (NUM_LC_LD - 1)*clustersNumber)
		real_len_s = verticesNumber;
	for (i = 1; i <= clustersNumber; i++)
	{
		if (currentSol.prototypeArray[i].vertexIndex != -1)
			sol[i - 1] = currentSol.prototypeArray[i].vertexIndex - 1;
		else
		{
			wp_isExist = 0;
			for (int j = verticesNumber; j < real_len_s; j++)
			{
				if (wp_points[j].interior.edgeIndex == currentSol.prototypeArray[i].interior.edgeIndex && fabs(wp_points[j].interior.length - currentSol.prototypeArray[i].interior.length) < DEVIATION)
				{
					wp_isExist = 1;
					wp_number = j;
					break;
				}
			}
			if (!wp_isExist && real_len_s<verticesNumber + NUM_LC_LD*clustersNumber)
			{
				sol[i - 1] = real_len_s;
				wp_points[real_len_s].isVertex = false;
				wp_points[real_len_s].vertex_number = real_len_s;
				wp_points[real_len_s].interior.edgeIndex = currentSol.prototypeArray[i].interior.edgeIndex;
				wp_points[real_len_s].interior.length = currentSol.prototypeArray[i].interior.length;
				real_len_s++;
			}
			else if (wp_isExist)
			{
				sol[i - 1] = wp_number;
			}
		}
	}
	//cout << "real_len_s=" << real_len_s <<",objective="<<currentSol.objective<< endl;
	//printPrototypes(currentSol.prototypeArray);
	for (i = 0; i < verticesNumber + NUM_LC_LD *clustersNumber; i++)
		Tb.tabuList[i] = 0;
	iter = 0; num_p = clustersNumber;
	objective = currentSol.objective;
	while (iter <50)
	{
		dial = rand() % 2;
		//first drop move£¬followed by add move
		if (dial == 0)
		{
			build_data(dis, cluLen, clu, num_p, sol);
			drop_node = drop_func(sol, iter, Tb.tabuList, cluLen, clu, dis, num_p, objective);

			build_data(dis, cluLen, clu, num_p, sol);
			add_node = add_func(sol, iter, Tb.tabuList, dis, num_p, objective);
		}
		//first add move£¬followed by drop move
		else
		{
			build_data(dis, cluLen, clu, num_p, sol);
			add_node = add_func(sol, iter, Tb.tabuList, dis, num_p, objective);

			build_data(dis, cluLen, clu, num_p, sol);
			drop_node = drop_func(sol, iter, Tb.tabuList, cluLen, clu, dis, num_p, objective);
		}
		iter++;
		if (bestSolution.objective - objective > DEVIATION)						//feasible solution		
		{
			copy_information(sol, objective, bestSolution);
			best_time_per = difftime(time(&time_1), startTime) * 1000;
		}
		if (objective < local_minima)
			copy_information(sol, objective, currentSol);
	}
	return currentSol.objective;
}

//tabu search as the discrete local search£¬ and k-means as the continous local search
void reformulation_local_search(char *argv[], Clustering &currentSolution)
{
	int flag_stop = 0;
	int iter = 0;
	double xc, xd;
	bool* remainingPrototypes = new bool[clustersNumber + 1];
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		remainingPrototypes[cluster] = false;
	}
	tabu_search(currentSolution);
	//cout << "here 2 " << endl;
	copyAllInformation(currentSolution, tmpSolution_update);
	while (!flag_stop &&iter<5)
	{
		//step 2
		//local_continous(argv, currentSolution, 0);
		k_means(currentSolution, remainingPrototypes, argv[5][0]);
		//cout << "here 3 " << endl;
		if (bestSolution.objective -currentSolution.objective > DEVIATION)
		{
			copyAllInformation(currentSolution, bestSolution);
			best_time_per = difftime(time(&time_1), startTime) * 1000;
		}
		xc = currentSolution.objective;
		xd = tabu_search(currentSolution);
		if (currentSolution.objective<tmpSolution_update.objective)
			copyAllInformation(currentSolution, tmpSolution_update);
		//printf("xc=%2f,xd=%2f\n", xc, xd);
		if (fabs(xd - xc) < DEVIATION)
			flag_stop = 1;
		iter++;
	}
	delete [] remainingPrototypes;	
}

void initial_population(char *argv[])
{
	int i = 1, j, flag_dif;
	bool* remainingPrototypes = new bool[clustersNumber + 1];
	while (i <= PopNum &&difftime(endTime,startTime)<timeLimit)
	{
		flag_dif = 0;
		clearAllTabuInformation();
		for (j = 1; j <= clustersNumber; j++)
			generate_prototype_on_vertex_or_edge_50_50(j, population[i].prototypeArray);
		//printPrototypes(population[i].prototypeArray);
		for (int cluster = 1; cluster <= clustersNumber; cluster++) {
			remainingPrototypes[cluster] = false;
		}
		allocate(population[i]);
		//cout << "here 1" << endl;
		reformulation_local_search(argv, population[i]);
		
		//cout << "i=" << i << ",objective=" << population[i].objective << ",bestSolution.objective=" << bestSolution.objective << endl;
		printf("i=%d,objective=%2f,bestSolution.objective=%2f\n", i, population[i].objective, bestSolution.objective);
		for (j = 1; j < i; j++)
		{
			if (fabs(population[i].objective - population[j].objective)<DEVIATION)
			{
				flag_dif = 1;
				break;
			}
		}
		//cout << "pop" << " " << i << ",objective=" << population[i].objective << endl;
		if (flag_dif == 0)
			i++;
		time(&endTime);
	}
	delete [] remainingPrototypes;
}

//crossover: common shared prototypes
void cross_over(int *vertex_proto, int *edge_proto, int *vertex_proto_2, int *edge_proto_2, double *interior_edge, double *interior_edge_2)
{
	int vertex_proto_len = 0;
	int edge_proto_len = 0;
	int vertex_proto_len_2 = 0;
	int edge_proto_len_2 = 0;
	int par1, par2;
	int index1 = 0, index2 = 0, len_ch = 1;
	std::uniform_int_distribution<int> parent(1, PopNum);
	par1 = parent(rd);
	par2 = parent(rd);
	while (par2 == par1)
		par2 = parent(rd);
	for (int i = 1; i <= clustersNumber; i++)
	{
		if (population[par1].prototypeArray[i].vertexIndex != -1)
			vertex_proto[vertex_proto_len++] = population[par1].prototypeArray[i].vertexIndex;
		else
		{
			edge_proto[edge_proto_len++] = population[par1].prototypeArray[i].interior.edgeIndex;
			interior_edge[population[par1].prototypeArray[i].interior.edgeIndex] = population[par1].prototypeArray[i].interior.length;
		}

		if (population[par2].prototypeArray[i].vertexIndex != -1)
			vertex_proto_2[vertex_proto_len_2++] = population[par2].prototypeArray[i].vertexIndex;
		else
		{
			edge_proto_2[edge_proto_len_2++] = population[par2].prototypeArray[i].interior.edgeIndex;
			interior_edge_2[population[par2].prototypeArray[i].interior.edgeIndex] = population[par2].prototypeArray[i].interior.length;
		}
	}
	bubble_sort(vertex_proto, vertex_proto_len);
	bubble_sort(vertex_proto_2, vertex_proto_len_2);
	bubble_sort(edge_proto, edge_proto_len);
	bubble_sort(edge_proto_2, edge_proto_len_2);
	for (int i = 1; i <= clustersNumber; i++)
	{
		currentSolution.prototypeArray[i].vertexIndex = -1;
		currentSolution.prototypeArray[i].interior.edgeIndex = -1;
	}
	while (index1 < vertex_proto_len && index2 < vertex_proto_len_2)
	{
		if (vertex_proto[index1] > vertex_proto_2[index2])
			index2++;
		else if (vertex_proto[index1] == vertex_proto_2[index2])
		{
			currentSolution.prototypeArray[len_ch].vertexIndex = vertex_proto[index1];
			currentSolution.prototypeArray[len_ch].interior.edgeIndex = -1;
			len_ch++;
			index2++;
			index1++;
		}
		else
			index1++;
	}
	index1 = 0; index2 = 0;
	while (index1 < edge_proto_len && index2 < edge_proto_len_2)
	{
		if (edge_proto[index1] > edge_proto_2[index2])
			index2++;
		else if (edge_proto[index1] == edge_proto_2[index2] && fabs(interior_edge[edge_proto[index1]] - interior_edge_2[edge_proto_2[index2]]) < DEVIATION)
			//else if (edge_proto[index1] == edge_proto_2[index2])
		{
			currentSolution.prototypeArray[len_ch].interior.edgeIndex = edge_proto[index1];
			currentSolution.prototypeArray[len_ch].interior.length = interior_edge[edge_proto[index1]];
			currentSolution.prototypeArray[len_ch].vertexIndex = -1;
			len_ch++;
			index2++;
			index1++;
		}
		else
			index1++;
	}
#ifdef DEBUG1
	cout << "len_ch=" << len_ch - 1 << endl;
	getchar();
#endif
	clearAllTabuInformation();
	for (int i = 1; i <= len_ch - 1; i++)
	{	updateTabuInformation(currentSolution.prototypeArray[i]);		
	}
	while (len_ch <= clustersNumber)
	{
		generate_prototype_on_vertex_or_edge_50_50(len_ch, currentSolution.prototypeArray);
		len_ch++;
	}
	allocate(currentSolution);
#ifdef DEBUG1
	cout << "currentSolution(child):" << endl;
	printPrototypes(currentSolution.prototypeArray);
	getchar();
#endif
}

void mutation(Clustering &currentSolution, char *argv[])
{
	//cout << "currentSolution.objective=" << currentSolution.objective << endl;
	//int k = rand() % 5 + 1;
	int k = 0.1*clustersNumber + rand() % 5;
	bool* remainingPrototypes = new bool[clustersNumber + 1];
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		remainingPrototypes[cluster] = false;
	}
	shaking(remainingPrototypes, k, currentSolution.prototypeArray, argv[4][0]);
	allocate(currentSolution);
	delete [] remainingPrototypes;
	//cout << "after mutation, currentSolution.objective=" << currentSolution.objective << endl;
}

void update_population(Clustering currentSolution)
{
	int i, index, flag_dif;
	double ff_max = -MAX_DOUBLE;
	for (i = 1; i <= PopNum; i++)
	{
		if (population[i].objective > ff_max)
		{
			ff_max = population[i].objective;
			index = i;
		}
	}
	flag_dif = 0;
	for (i = 1; i <= PopNum; i++)
	{
		if (fabs(currentSolution.objective - population[i].objective) < DEVIATION)							
			flag_dif = 1;
	}
	if (currentSolution.objective < ff_max && flag_dif == 0)
		copyAllInformation(currentSolution, population[index]);
}

void compute_similarity(int *vertex_proto, int *edge_proto, int *vertex_proto_2, int *edge_proto_2, double *interior_edge, double *interior_edge_2)
{
	int len_same = 0, index1, index2;
	double similar;
	int vertex_proto_len = 0;
	int edge_proto_len = 0;
	int vertex_proto_len_2 = 0;
	int edge_proto_len_2 = 0;
	for (int i = 1; i <= PopNum; i++)
	{
		for (int j = i + 1; j <= PopNum; j++)
		{
			vertex_proto_len = 0;
			edge_proto_len = 0;
			vertex_proto_len_2 = 0;
			edge_proto_len_2 = 0;
			for (int k = 1; k <= clustersNumber; k++)
			{
				if (population[i].prototypeArray[k].vertexIndex != -1)
					vertex_proto[vertex_proto_len++] = population[i].prototypeArray[k].vertexIndex;
				else
				{
					edge_proto[edge_proto_len++] = population[i].prototypeArray[k].interior.edgeIndex;
					interior_edge[population[i].prototypeArray[k].interior.edgeIndex] = population[i].prototypeArray[k].interior.length;
				}

				if (population[j].prototypeArray[k].vertexIndex != -1)
					vertex_proto_2[vertex_proto_len_2++] = population[j].prototypeArray[k].vertexIndex;
				else
				{
					edge_proto_2[edge_proto_len_2++] = population[j].prototypeArray[k].interior.edgeIndex;
					interior_edge_2[population[j].prototypeArray[k].interior.edgeIndex] = population[j].prototypeArray[k].interior.length;
				}
			}
			bubble_sort(vertex_proto, vertex_proto_len);
			bubble_sort(vertex_proto_2, vertex_proto_len_2);
			bubble_sort(edge_proto, edge_proto_len);
			bubble_sort(edge_proto_2, edge_proto_len_2);
			index1 = 0;
			index2 = 0;
			while (index1 < vertex_proto_len && index2 < vertex_proto_len_2)
			{
				if (vertex_proto[index1] > vertex_proto_2[index2])
					index2++;
				else if (vertex_proto[index1] == vertex_proto_2[index2])
				{
					len_same++;
					index2++;
					index1++;
				}
				else
					index1++;
			}
			index1 = 0; index2 = 0;
			while (index1 < edge_proto_len && index2 < edge_proto_len_2)
			{
				if (edge_proto[index1] > edge_proto_2[index2])
					index2++;
				else if (edge_proto[index1] == edge_proto_2[index2] && fabs(interior_edge[edge_proto[index1]] - interior_edge_2[edge_proto_2[index2]]) < DEVIATION)
					//else if (edge_proto[index1] == edge_proto_2[index2])
				{
					len_same++;
					index2++;
					index1++;
				}
				else
					index1++;
			}
		}
	}
	similar = 1.0*len_same / (clustersNumber * PopNum*(PopNum - 1) / 2);
	cout << "similar=" << similar << endl;
}

//reformulation local search as the local search within the memetic frawework
double memetic(char *argv[])
{
	int iter = 0;
	real_len_s = verticesNumber;
	bestSolution.objective = DBL_MAX;
	int *vertex_proto = new int[clustersNumber + 1];
	int *edge_proto = new int[clustersNumber + 1];
	int *vertex_proto_2 = new int[clustersNumber + 1];
	int *edge_proto_2 = new int[clustersNumber + 1];
	double *interior_edge_length = new double[edgesNumber + 1];
	double *interior_edge_length_2 = new double[edgesNumber + 1];
	vertexDiscreteDistribution = new std::uniform_int_distribution<int>(1, verticesNumber);
	edgeDiscreteDistribution = new std::uniform_int_distribution<int>(1, edgesNumber);
	time(&startTime);
	initial_population(argv);
	//while (iter < 100)
	while (difftime(endTime, startTime)<timeLimit)
	{
		cross_over(vertex_proto, edge_proto, vertex_proto_2, edge_proto_2, interior_edge_length, interior_edge_length_2);
		mutation(currentSolution, argv);
		reformulation_local_search(argv, currentSolution);
		update_population(tmpSolution_update);
		printf("generations=%d,objective=%.2f,ObjBest=%.2f\n", iter++, currentSolution.objective, bestSolution.objective);
		time(&endTime);
	}	
	compute_similarity(vertex_proto, edge_proto, vertex_proto_2, edge_proto_2, interior_edge_length, interior_edge_length_2);
	delete [] vertex_proto;
	delete [] edge_proto;
	delete [] vertex_proto_2;
	delete [] edge_proto_2;
	delete [] interior_edge_length;
	delete [] interior_edge_length_2;
	return bestSolution.objective;
}
/*
int main(int argc, char * argv[])
{
	srand(unsigned(time(NULL)));
	if (argc < 6)
	{
		exitOnError("Usage: program.x <filename>", "");
	}
	argv[1] = "F:\\MSSC\\MSSC_net\\dataset\\pmed21.txt";
	timeLimit = atoi(argv[2]);
	Tb.tabuTenure = 10;
	readInstanceAndInitializeVariables(argv[1]);
	Floyd_harshall_algorithm();
	vertexDiscreteDistribution = new std::uniform_int_distribution<int>(1, verticesNumber);
	edgeDiscreteDistribution = new std::uniform_int_distribution<int>(1, edgesNumber);
	bestSolution.objective = MAX_DOUBLE;
	memetic(argv);
	printf("obj_mem=%2f\n", bestSolution.objective);
	getchar();
	return 0;
}*/

int main(int argc, char *argv[])
{
	srand(unsigned(time(NULL)));
	if (argc < 6)
	{
		exitOnError("Usage: program.x <filename>", "");
	}
	int instance_num = 0;
	int instance_total_num = 40;
	char **argv_ins = new char*[instance_total_num];				//the number of instances
	for (int i = 0; i < instance_total_num; i++)
		argv_ins[i] = new char[100];
	double best_value = DBL_MAX;
	double avg_value = 0;
	double avg_time = 0;

	ofstream resultsFile;
	ofstream valuesFile;
	resultsFile.open("E:\\resultados_MA_LC_LD_2.txt", ofstream::app);
	valuesFile.open("E:\\solucoes_MA_LC_LD_2.txt", ofstream::app);
	const char *to_search = "F:\\MSSC\\MSSC_net\\dataset\\*.txt";    
	char *path = new char[100];
	long handle;                                                    
	struct _finddata_t fileinfo;                                 
	handle = _findfirst(to_search, &fileinfo);                      
	if (-1 == handle)
		return -1;
	strcpy(argv_ins[instance_num++], fileinfo.name);
	while (!_findnext(handle, &fileinfo))                              //looks up the "to_search" directory to find all the instances
	{
		strcpy(argv_ins[instance_num++], fileinfo.name);
		cout << "name=" << fileinfo.name << endl;
	}
	_findclose(handle);                                         	  //close handle
	int runs = 10;
	timeLimit = atoi(argv[2]);
	Tb.tabuTenure = 20;	
	char *file_name = new char[100];
	double temp_best;
	for (int i = 0; i < instance_num; i++)
	{
		strcpy(path, "F:\\MSSC\\MSSC_net\\dataset\\");
		strcpy(file_name, strcat(path, argv_ins[i]));
		readInstanceAndInitializeVariables(file_name);
		Floyd_harshall_algorithm();
		valuesFile << file_name << ":";
		avg_value = 0;
		best_value = DBL_MAX;
		for (int run = 0; run < runs; run++)
		{
			memetic(argv);
			temp_best = bestSolution.objective;
			avg_value += temp_best;
			avg_time += best_time_per;
			if (temp_best < best_value)
				best_value = temp_best;
			valuesFile << setprecision(2) << setiosflags(ios::fixed) << temp_best << ";";
		}
		resultsFile << file_name << ":bestV=" << setprecision(2) << setiosflags(ios::fixed) << best_value;
		resultsFile << ",avgV=" << setprecision(2) << setiosflags(ios::fixed) << avg_value / runs;
		resultsFile << ",avg_time=" << setprecision(2) << setiosflags(ios::fixed) << avg_time / runs << endl;
		valuesFile << endl;
	}
	getchar();
}