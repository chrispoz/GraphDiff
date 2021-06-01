#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <functional>
#include <thread>
#include "JSONBoostGraphMap.hpp"
#include "json_boost_graphviz_utils.hpp"
#include "json_boost_graphviz_enums.hpp"

using namespace std;
using boost::function;
using std::vector;
using std::map;
using std::string;

int main(int argc, char **argv)
{
	//
	//	Get the JSON strings to compare.
	if(argc < 7)
	{
		std::cout << "Requires 6 arguments - terminating" << std::endl;
		return 0;
	}

	//	First pair of files - JSON input files to diff.
	string json1 = argv[1];
	string json2 = argv[2];

	//	2nd pair of files - graphviz output file names.
	string outputFile1 = argv[3];
	string outputFile2 = argv[4];

	//	3rd pair of files - image output file names.
	string outputGraphicsFile1 = argv[5];
	string outputGraphicsFile2 = argv[6];

	//	Declare input streams and buffers to read files.
	ifstream if1(json1.c_str());
	string json1Buffer((std::istreambuf_iterator<char>(if1)),
		std::istreambuf_iterator<char>());

	ifstream if2(json2.c_str());
	string json2Buffer((std::istreambuf_iterator<char>(if2)),
		std::istreambuf_iterator<char>());

	//	Load the json files into the json graph mapper.
	DiGraph newDg;
	JSONBoostGraphMapping<DiGraph, NamePropMap, GraphPropMap> graphMap(newDg);
	if(graphMap.loadJSONToDAG(json1Buffer.c_str(), SCHEME_NODE_SET_TYPES, SET_1_NOT_2))
	{
		graphMap.createOutput(outputFile1, outputGraphicsFile1);
	}

	DiGraph newDg2;
	JSONBoostGraphMapping<DiGraph, NamePropMap, GraphPropMap> graphMap2(newDg2);
	if(graphMap2.loadJSONToDAG(json2Buffer.c_str(), SCHEME_NODE_SET_TYPES, SET_2_NOT_1))
	{
		graphMap2.createOutput(outputFile2, outputGraphicsFile2);
	}

	//	Create an output graph and a graph property writer to create the output 
	//	graphviz and image files.
	DiGraph gOut;
	GraphPropertyWriter<DiGraph, GraphPropMap> outGraph(gOut);

	//	Diff the graph maps and populate the outGraph with diff results.
	if(graphMap2.diff2(graphMap, outGraph))
	{
		//	Define the output files
		string outputFileDiff = outputFile1 + ".diff.gv";
		string outputGraphicsFileDiff = outputGraphicsFile1 + "diff.png";

		//	Create the output.
		outGraph.createOutput(outputFileDiff, outputGraphicsFileDiff);
	}
	else
	{
		std::cout << "The graphs passed all equivalency tests." << std::endl;
	}  

/*
	DiGraph newDg;
	JSONBoostGraphMapping<DiGraph, NamePropMap, GraphPropMap> graphMap(newDg);
	if(graphMap.loadJSON(json1Buffer.c_str(), SCHEME_NODE_SET_TYPES, SET_1))
	{
		graphMap.createOutput(outputFile1, outputGraphicsFile1);
	}

	DiGraph newDg2;
	JSONBoostGraphMapping<DiGraph, NamePropMap, GraphPropMap> graphMap2(newDg2);
	if(graphMap2.loadJSON(json2Buffer.c_str(), SCHEME_NODE_SET_TYPES, SET_2))
	{
		graphMap2.createOutput(outputFile2, outputGraphicsFile2);
	}

	//	Create an output graph and a graph property writer to create the output 
	//	graphviz and image files.
	DiGraph gOut;
	//GraphPropertyWriter<DiGraph, NamePropMap, GraphPropMap> outGraph(gOut);
	GraphPropertyWriter<DiGraph, GraphPropMap> outGraph(gOut);

	//	Diff the graph maps and populate the outGraph with diff results.
	if(graphMap2.diff2(graphMap, outGraph))
	{
		//	Define the output files
		string outputFileDiff = outputFile1 + ".diff.gv";
		string outputGraphicsFileDiff = outputGraphicsFile1 + "diff.png";

		//	Create the output.
		outGraph.createOutput(outputFileDiff, outputGraphicsFileDiff);
	}
	else
	{
		std::cout << "The graphs passed all equivalency tests." << std::endl;
	}  
*/


	return 0;
}



/*
	//constructBidirectionalGraph("bidirectional.gv", "bidirectional.png");

	DoubleArray da(3);
	da[0] = 5.76;
	da[1] = 23.04;
	da[2] = 96.04;
	IntArray di(4);
	di[0] = 1;
	di[1] = 3;
	di[2] = 16;
	di[3] = 10;
	DoubleArrayP daP(3);
	daP[0] = 15.76;
	daP[1] = 123.04;
	daP[2] = 2396.04;

	//	Use a template function on the array.
	//double sum = array_sum<DoubleArray>(da, 3);
	array_traits<DoubleArray>::value_type dSum = sum<DoubleArray>(da, 3);
	array_traits<DoubleArrayP>::value_type dpSum = sum<DoubleArrayP>(daP, 3);
	array_traits<IntArray>::value_type iSum = sum<IntArray>(di, 4);
	std::cout << "Sum = " << dSum << std::endl;
	std::cout << "Sum = " << iSum << std::endl;
	std::cout << "Sum = " << dpSum << std::endl;


class DoubleArrayP
{
public:
	typedef double value_type;
	DoubleArrayP(int size = 10)
	{
		m_data = new double[size];
		for(int i = 0; i < size; ++i)
		{
			m_data[i] = 0.0;
		}
	}
	~DoubleArrayP()
	{
		delete [] m_data;
	}
	double& operator[](int i)
	{
		return m_data[i];
	}
private:
	double *m_data;
};

class DoubleArray
{
public:
	DoubleArray(int size = 10)
	{
		m_data.assign(size, 0.0);
	}
	typedef double value_type;
	double& operator[](int i)
	{
		return m_data[i];
	}

private:
	std::vector<double> m_data;
};

class IntArray
{
public:
	IntArray(int size = 10)
	{
		m_data.assign(size, 0);
	}
	typedef int value_type;
	int& operator[](int i)
	{
		return m_data[i];
	}

private:
	std::vector<int> m_data;
};

//	Main traits definition.
template <typename Array>
struct array_traits
{
	typedef typename Array::value_type value_type;
};


template <>
struct array_traits<double *>
{
	typedef double value_type;
};


template <typename Array>
typename array_traits<Array>::value_type sum(Array& v, int n)
{
	typename array_traits<Array>::value_type total = 0;
	for(int i = 0; i < n; i++)
	{
		total += v[i];
	}

	return total;
}

rapidjson::Value::ConstMemberIterator mapJSONToGraph(rapidjson::Document& doc,
	rapidjson::Value::ConstMemberIterator iterIn, 
	vector< JSONToBoostMap >& parentListMap,
	map<JSONNode, int>& objNodeMap,
	int& vLevel,
	DiGraphWrapper& graph);


//	Given a digraph with properties, an output viz file path and
//	output png file path, create the viz file and graphics file from the graph.
void outputGraph(DiGraphWrapper& wrapper, const std::string& outVizFile,
				 const std::string& outPngFile);

void constructBidirectionalGraph(const std::string& outFileViz,
	const std::string& outFilePng);
void writeGraphOutputFile(const string& dotInFile, const string& pngFileOut);

rapidjson::Value::ConstMemberIterator  mapJSONToGraph(rapidjson::Document& doc,
	rapidjson::Value::ConstMemberIterator iterIn,
	vector< JSONToBoostMap >& parentListMap,
	map<JSONNode, int>& objNodeMap,
	int& vLevel,
	DiGraphWrapper& graph)
{
	if(iterIn == doc.MemberEnd())
	{
		return iterIn;
	}

	bool bIsRootNode = (iterIn == doc.MemberBegin());

	//	See if it's already in the map.
	map<JSONNode, int>::const_iterator mapIter =
		objNodeMap.find(iterIn);

	bool bAlreadyExists = (mapIter != objNodeMap.end());

	string type = kTypeNames[iterIn->value.GetType()] ;

	if(type == "Object")
	{
		//  For an object, it is important that we know that a member
		//	has already been visited.  We can do this by keeping a list
		//	of the visited nodes and reference them by their member pointer.
		if(!bAlreadyExists)
		{
			std::cout << "Type - " << type << " name: " << iterIn->name.GetString() << endl;
			objNodeMap.insert(std::make_pair(iterIn, objNodeMap.size()));

			//	Create the vertex and set the nodeNames and nodeIndices properties.
			Vertex v = boost::add_vertex(graph.g_);
			graph.nodeNames[v] = iterIn->name.GetString();
			graph.nodeIndices[v] = objNodeMap.size();

			if(bIsRootNode)
			{
				parentListMap.push_back(JSONToBoostMap(JSONNode(iterIn), v));

			}
			else
			{
				parentListMap.push_back(JSONToBoostMap(JSONNode(iterIn), v));
				Vertex parent = parentListMap[vLevel - 1].second();
				boost::add_edge(parent, v, graph.g_);
			}

			if(iterIn->value.MemberBegin() != iterIn->value.MemberEnd())
			{
				vLevel++;
				return iterIn->value.MemberBegin();
			}
		}
		else
		{
			//	Go through the children and find one that doesn't already exist
			rapidjson::Value::ConstMemberIterator childIter = iterIn->value.MemberBegin();

			bool bFound = false;
			while(childIter != iterIn->value.MemberEnd())
			{
				mapIter = objNodeMap.find(childIter);

				if(mapIter == objNodeMap.end())
				{
					//	NOTE -- Don't add to map - function will call itself next
					//	with this iterator.
					return childIter;
				}
				++childIter;
			}

			//	If we get here, we did not find any children that haven't been set.
			return iterIn;
		}
	}
	else if(type == "Array")
	{
		if(mapIter == objNodeMap.end())
		{
			objNodeMap.insert(std::make_pair(JSONNode(iterIn), objNodeMap.size()));
		}
		//	Get the size of the array.
		const rapidjson::SizeType size = iterIn->value.Size();
		std::cout << "Type - " << type << " name: " << iterIn->name.GetString() << endl;
		std::cout << "\tElements:\t[";

		rapidjson::Value::ConstValueIterator valIter = iterIn->value.Begin();
		while(valIter != iterIn->value.End())
		{
			objNodeMap.insert(std::make_pair(valIter, objNodeMap.size()));

			Vertex v = boost::add_vertex(graph.g_);

			//	Add parent edge to element
			Vertex parent = parentListMap[vLevel - 1].second();
			boost::add_edge(parent, v, graph.g_);

			string val = valIter->GetString();
			graph.nodeNames[v] = val;
			graph.nodeIndices[v] = objNodeMap.size();

			auto tmpIter = valIter;
			tmpIter++;
			if(tmpIter != iterIn->value.End())
			{
				std::cout << val << ",";
			}
			else
			{
				std::cout << val;
			}
			valIter++;
		}
		std::cout << "]" << endl;
		return iterIn;
	}
	else
	{
		if(mapIter == objNodeMap.end())
		{
			objNodeMap.insert(std::make_pair(iterIn, objNodeMap.size()));

			Vertex v = boost::add_vertex(graph.g_);
			string str = iterIn->name.GetString();
			str += "/";
			str += iterIn->value.GetString();
			graph.nodeNames[v] = str;
			graph.nodeIndices[v] = objNodeMap.size();

			Vertex parent = parentListMap[vLevel - 1].second();
			boost::add_edge(parent, v, graph.g_);

		}
		std::cout << "Type - " << type << " name: " << iterIn->name.GetString() << " value: " << iterIn->value.GetString() << endl;
	}

	return ++iterIn;
}

//	Graph will take as inputs the number of
//	vertices, the number of edges, 
void constructBidirectionalGraph(const std::string& outFileViz,
	const std::string& outFilePng)
{
	//	Start with 5 vertices A-E are the names.
	enum labels {A = 0, B, C, D, E, F, G};
	static char *labelNames [] =
	{
		"A",
		"B",
		"C",
		"D",
		"E",
		"F",
		"G"
	};

	//	Define the graph.
	const int nVertices = (G - A) + 1;
	DiGraph graph(nVertices);

	typedef std::pair<int, int> IntEdge;
	IntEdge edge_array [] =
	{
		IntEdge(A, B),
		IntEdge(A, E),
		IntEdge(B, C),
		IntEdge(C, D),
		IntEdge(B, D),
		IntEdge(E, F),
		IntEdge(E, G)
	};

	const int numEdges = sizeof(edge_array)/sizeof(edge_array[0]);

	IndexMap index = boost::get(boost::vertex_index, graph);
	typedef boost::graph_traits<DiGraph>::vertex_iterator vertex_iter;
	std::pair<vertex_iter, vertex_iter> vp;

	for(vp = boost::vertices(graph); vp.first != vp.second; ++vp.first)
	{
		Vertex v = *vp.first;

		std::cout << "Index = " << index[v] << std::endl;
	}

	//	Construct the edges
	for(int i = 0; i < numEdges; ++i)
	{
		boost::add_edge(edge_array[i].first, edge_array[i].second, graph);
	}

	//	Now write the graph out to graphviz format.  
	std::ofstream outputViz(outFileViz, std::ofstream::out);

	boost::write_graphviz(outputViz, graph, boost::make_label_writer(labelNames));

	//	Close out the graphviz file.
	outputViz.close();

	//	Now that we have it, write it to file.
	writeGraphOutputFile(outFileViz, outFilePng);
}

void outputGraph(DiGraphWrapper& wrapper, const std::string& outVizFile,
				 const std::string& outPngFile)
{
	//	Now write the graph out to graphviz format.  
	std::ofstream outputViz(outVizFile, std::ofstream::out);

	boost::write_graphviz(outputViz, wrapper.g_, boost::make_label_writer(wrapper.nodeNames));

	//	Close out the graphviz file.
	outputViz.close();

	//	Now that we have it, write it to file.
	writeGraphOutputFile(outVizFile, outPngFile);
}
*/