#ifndef _JSON_BOOST_GRAPH_MAP_HPP_
#define _JSON_BOOST_GRAPH_MAP_HPP_

//	Shettep MSFT
#pragma warning(disable : 4018)

#include <boost/graph/directed_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graphviz.hpp>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <graphviz/cgraph.h>
#include <graphviz/gvc.h>
#include <graphviz/cdt.h>
#include <iostream>
#include <utility>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "../rapidjson/include/document.h"

#include "json_boost_graphviz_enums.hpp"
#include "json_boost_graphviz_utils.hpp"

typedef GraphPropertyWriter<DiGraph, GraphPropMap> GraphOutputGenerator;
typedef std::map<std::string, GraphOutputGenerator > OutputGraphMap;


union NodeData
{
	rapidjson::Value::ConstMemberIterator member_;
	rapidjson::Value::ConstValueIterator value_;
};

struct JSONNode
{
	bool operator<(const JSONNode& other) const
	{
		return (data_.member_ < other.data_.member_) ||
			(data_.value_ < other.data_.value_);
	}
	JSONNode(rapidjson::Value::ConstMemberIterator p)
	{
		data_.member_ = p;
	}
	JSONNode(rapidjson::Value::ConstValueIterator p)
	{
		data_.value_ = p;
	}
	NodeData data_;
};

class ParentNodeMap
{
public:
	//	Construct parent nodes with the node itself and the
	//	graph so we can insert sibling edges.
	ParentNodeMap(JSONToBoostNodeMap parent, DiGraph& graph) :
		parentNode(parent),
		graph_(graph)
	{
	}

	bool addChildNode(BoostNode childNode)
	{
		//	Create an edge from the last sibling.
		if(childNodes.size() > 0)
		{
			//	Create an edge from the sibling to the current node
			boost::add_edge(childNodes[childNodes.size() - 1], childNode, graph_);
		}

		//	Add the child node to the list
		childNodes.push_back(childNode);

		return true;
	}

	JSONNode& first()
	{
		return parentNode.first();
	}

	Vertex& second()
	{
		return parentNode.second();
	}

	ParentNodeMap& operator=(const ParentNodeMap& nodeMap)
	{
		parentNode = nodeMap.parentNode;
		graph_ = nodeMap.graph_;
		childNodes = nodeMap.childNodes;
		return *this;
	}
	JSONToBoostNodeMap parentNode;
	std::vector<BoostNode> childNodes;
	DiGraph& graph_;
};

template <class Graph, class NameMap, class Index2AttributeMap>
class JSONBoostGraphMapping
{
public:
	JSONBoostGraphMapping(Graph& g):
		graphDataSet_(g)
	{
	}

	virtual ~JSONBoostGraphMapping()
	{
	}

	bool loadJSONToDAG(const std::string& jsonString, MapColorSchemes colorScheme, int colorValue)
	{
		if(document_.Parse<0>(jsonString.c_str()).HasParseError())
		{
			std::cout << document_.GetParseError() << " - terminating." << std::endl;
			return false;
 		}

		rapidjson::Value::ConstMemberIterator docIter = document_.MemberBegin();
		rapidjson::Value::ConstMemberIterator docEnd = document_.MemberEnd();

		std::map<JSONNode, int> objNodeMap;


		//	Store the parent node data.
		std::vector< ParentNodeMap > parentListMap;
		//	Define variable to track nesting
		int vNestLevel = 0;

		while(docIter != docEnd)
		{
			auto tmpIter = docIter;

			docIter = mapJSONToDAG(document_, docIter,  
				parentListMap, objNodeMap, vNestLevel, colorScheme, colorValue);
			if(tmpIter == docIter)
			{
				//	Now find the correct parent.
				if(vNestLevel > 0)
				{
					--vNestLevel;

					//	Remove the last entry.
					parentListMap.erase(parentListMap.begin() + vNestLevel);
				}
				else
				{
				}

				//	Get next iterator - depends on what level we're on.
				//	When level is 0, we have completely visited the subgraph and can advance
				//	to the next sibling.
				docIter = (vNestLevel > 0) ? parentListMap[vNestLevel - 1].first().data_.member_ : ++docIter;
			}
		}

		return true;
	}

	bool loadJSON(const std::string& jsonString, MapColorSchemes colorScheme, int colorValue)
	{
		if(document_.Parse<0>(jsonString.c_str()).HasParseError())
		{
			std::cout << document_.GetParseError() << " - terminating." << std::endl;
			return false;
 		}

		rapidjson::Value::ConstMemberIterator docIter = document_.MemberBegin();
		rapidjson::Value::ConstMemberIterator docEnd = document_.MemberEnd();

		std::map<JSONNode, int> objNodeMap;
		std::vector< JSONToBoostMap > parentListMap;

		//	Define variable to track nesting
		int vNestLevel = 0;

		while(docIter != docEnd)
		{
			auto tmpIter = docIter;

			docIter = mapJSONToGraph(document_, docIter,  
				parentListMap, objNodeMap, vNestLevel, colorScheme, colorValue);
			if(tmpIter == docIter)
			{
				//	Now find the correct parent.
				if(vNestLevel > 0)
				{
					--vNestLevel;

					//	Remove the last entry.
					parentListMap.erase(parentListMap.begin() + vNestLevel);
				}

				//	Get next iterator - depends on what level we're on.
				//	When level is 0, we have completely visited the subgraph and can advance
				//	to the next sibling.
				docIter = (vNestLevel > 0) ? parentListMap[vNestLevel - 1].first().data_.member_ : ++docIter;
			}
		}

		return true;
	}

	const Graph& getGraph() const;
	const rapidjson::Document& getDocument() const;

	//	Given a digraph with properties, an output viz file path and
	//	output png file path, create the viz file and graphics file from the graph.
	void createOutput(const std::string& outVizFile,
			const std::string& outPngFile)
	{
		graphDataSet_.createOutput(outVizFile, outPngFile);
	}

	//	New diff function:  support notion of sibling nodes in graph.
	//	We want the LHS diff nodes G1\G2 to be on the LHS
	//	and RHS diff nodes G2\G1 to be on the RHS.
	bool diff2(JSONBoostGraphMapping& graphIn, 
		GraphPropertyWriter<Graph, Index2AttributeMap>& objectData)
	{
		auto vertexIter1 = boost::vertices(graphDataSet_.g_);
		auto vertexIter2 = boost::vertices(graphIn.graphDataSet_.g_);
		
		//	Use an ordered list visitor to generate the
		//	ordered vertex lists so we can visit each
		//	vertex in each graph in a standard loop in the same order as a full dfs traversal.
		std::vector<Vertex> orderedList1;
		std::vector<int> objectIndexList1;
		dfs_ordered_list_visitor<Vertex, Index2AttributeMap> orderedListVisitor1(orderedList1, objectIndexList1,
			graphDataSet_.propertyMap_);
		boost::depth_first_search(graphDataSet_.g_, boost::visitor(orderedListVisitor1));

		std::vector<Vertex> orderedList2;
		std::vector<int> objectIndexList2;
		dfs_ordered_list_visitor<Vertex, Index2AttributeMap> orderedListVisitor2(orderedList2, objectIndexList2,
			graphIn.graphDataSet_.propertyMap_);
		boost::depth_first_search(graphIn.graphDataSet_.g_, boost::visitor(orderedListVisitor2));

		std::cout << "Input graphs are being validated." << std::endl;

		std::set<Vertex> set1;
		std::set<int> indexSet1;
		std::set<GraphPropertyStruct> propSet1;
		dfs_set_visitor<Vertex> setVisitor1(set1, indexSet1, propSet1);
		boost::depth_first_search(graphDataSet_.g_, boost::visitor(setVisitor1));
		
		std::set<Vertex> set2;
		std::set<int> indexSet2;
		std::set<GraphPropertyStruct> propSet2;

		dfs_set_visitor<Vertex> setVisitor2(set2, indexSet2, propSet2);
		boost::depth_first_search(graphIn.graphDataSet_.g_, boost::visitor(setVisitor2));

		//	First check is to compare the size of the lists to the size of the sets.
		//	This ensures each vertex in each graph is uniquely defined.
		if((set1.size() == orderedList1.size()) &&
			(set2.size() == orderedList2.size()))
		{
			std::cout << "Graphs are properly defined.  Continuing. " << std::endl;

			//	Get the set intersection of the object properties from P1->P2.
			std::vector<GraphPropertyStruct> p12IntersectionSet(std::min(set1.size(), set2.size()));
			std::vector<GraphPropertyStruct>::iterator gps12IntersectionInter = std::set_intersection(propSet1.begin(), propSet1.end(),
				propSet2.begin(), propSet2.end(), p12IntersectionSet.begin(), GraphPropertyStructComparator());
			std::cout << "-------- INTERSECTION ------------------" << std::endl;

			//	Resize and sort the output set.
			p12IntersectionSet.resize(gps12IntersectionInter - p12IntersectionSet.begin());
			sort(p12IntersectionSet.begin(), p12IntersectionSet.end(), GraphPropertyStructIndexComparator());
			for(int i = 0; i < p12IntersectionSet.size(); ++i)
			{
				if(p12IntersectionSet[i].dataType_ == rapidjson::Type::kStringType)
				{
					std::cout << "\tIndex: " << p12IntersectionSet[i].index_ << ", Name: " << p12IntersectionSet[i].name_ << " has value \t" << p12IntersectionSet[i].value_ << std::endl;
				}
				else
				{
					std::cout << "Object Index: " << p12IntersectionSet[i].index_ << ", Name\t"  << p12IntersectionSet[i].name_ << std::endl;
				}
			}

			//	Get the set intersection of the object properties from P2->P1.
			std::vector<GraphPropertyStruct> p21IntersectionSet(std::min(set1.size(), set2.size()));
			std::vector<GraphPropertyStruct>::iterator gps21IntersectionInter = std::set_intersection(propSet2.begin(), propSet2.end(),
				propSet1.begin(), propSet1.end(), p21IntersectionSet.begin(), GraphPropertyStructComparator());
			p21IntersectionSet.resize(gps21IntersectionInter - p21IntersectionSet.begin());
			sort(p21IntersectionSet.begin(), p21IntersectionSet.end(), GraphPropertyStructIndexComparator());
			std::cout << "-------- INTERSECTION ------------------" << std::endl;

			//	Map the indexes of intersection from graph2 to graph1
			//	This allows us to use graph1 indexes when we create
			//	the output graph.
			//	We know the sizes must be the same.
			map<int,int> intersectionMapping;
			for(int i = 0; i < p12IntersectionSet.size(); ++i)
			{
				intersectionMapping.insert(std::make_pair(p21IntersectionSet[i].index_, p12IntersectionSet[i].index_));
			}

			//	Get the set difference of the object properties.
			std::vector<GraphPropertyStruct> p12DiffSet(std::max(set1.size(), set2.size()));
			std::vector<GraphPropertyStruct>::iterator gpsDiff12Iter = std::set_difference(propSet1.begin(), propSet1.end(),
				propSet2.begin(), propSet2.end(), p12DiffSet.begin(), GraphPropertyStructComparator());
			p12DiffSet.resize(gpsDiff12Iter - p12DiffSet.begin());
			std::sort(p12DiffSet.begin(), p12DiffSet.end(), GraphPropertyStructIndexComparator());

			std::cout << "-------- SET DIFF 1->2 ------------------" << std::endl;

			//	Resize the output set.
			for(int i = 0; i < p12DiffSet.size(); ++i)
			{
				if(p12DiffSet[i].dataType_ == rapidjson::Type::kStringType)
				{
					std::cout << "\tIndex: " << p12DiffSet[i].index_ <<  ", Name: " << p12DiffSet[i].name_ << " has value \t" << p12DiffSet[i].value_ << std::endl;
				}
				else
				{
					std::cout << "Object Index: " << p12DiffSet[i].index_ << ", Name\t"  << p12DiffSet[i].name_ << std::endl;
				}
			}

			std::vector<GraphPropertyStruct> p21DiffSet(std::max(set1.size(), set2.size()));
			std::vector<GraphPropertyStruct>::iterator gpsDiff21Iter = std::set_difference(propSet2.begin(), propSet2.end(),
				propSet1.begin(), propSet1.end(), p21DiffSet.begin(), GraphPropertyStructComparator());
			std::cout << "-------- SET DIFF 2->1 ------------------" << std::endl;
			p21DiffSet.resize(gpsDiff21Iter - p21DiffSet.begin());
			std::sort(p21DiffSet.begin(), p21DiffSet.end(), GraphPropertyStructIndexComparator());
			for(int i = 0; i < p21DiffSet.size(); ++i)
			{
				if(p21DiffSet[i].dataType_ == rapidjson::Type::kStringType)
				{
					std::cout << "\tIndex: " << p21DiffSet[i].index_ << ", Name: " << p21DiffSet[i].name_ << " has value \t" << p21DiffSet[i].value_ << std::endl;
				}
				else
				{
					std::cout << "Object Index: " << p21DiffSet[i].index_ << ", Name \t"  << p21DiffSet[i].name_ << std::endl;
				}
			}

			//	Set up output graph with name and properties map to graph 1 (this)
			//	and then update the properties in the graph to indicate the match
			//	status based on the data being in the intersection set or the p1->p2 difference set.
			//	Then add edges from from p2->p1 difference set.
			objectData.g_ = graphDataSet_.g_;
			objectData.propertyMap_ = graphDataSet_.propertyMap_;

			//	Go through the intersection index list for g1 and update the status to indicate a match.
			for(int i = 0; i < p12IntersectionSet.size(); ++i)
			{
				GraphPropertyStruct& gpsInt = p12IntersectionSet[i];
				int graph1Index = gpsInt.index_;

				Vertex curVertex = orderedList1[graph1Index];

				GraphPropertyStruct& gps = boost::get(boost::vertex_index2, objectData.g_, curVertex);
				gps.addUpdateNodeMappingType(MapColorSchemes::SCHEME_NODE_SET_TYPES, SetTypes::SET_INTERSECTION);
			}

			//	Go through the diff set and update the matching data
			for(int i = 0; i < p12DiffSet.size(); ++i)
			{
				GraphPropertyStruct& gps12Diff = p12DiffSet[i];
				int graph1Index = gps12Diff.index_;

				Vertex graph1Vertex = orderedList1[graph1Index];
				GraphPropertyStruct& gps = boost::get(boost::vertex_index2, objectData.g_, graph1Vertex);

				//	Update the node 
				gps.addUpdateNodeMappingType(MapColorSchemes::SCHEME_NODE_SET_TYPES, SetTypes::SET_2_NOT_1);
			}

			//	Go through the p21 diff set and add each vertex. 
			//	Also for each vertex, determine the parent.
			//	How to search for the parent:
			//	1) There should be exactly one in edge to the vertex - use
			//		in_edges to determine the parent (as an index into graph2).
			//	2) Take the index and see if there is a mapping back to graph1.
			//	   If YES --> get the vertex for the parent in graph1 and add the edge from that parent
			//	   If NO  --> Iterate through all of the vertex parents and repeat step2
			//	for each vertex, then .  Then create
			//	an edge from the parent to the new vertex.
			for(int i = 0; i < p21DiffSet.size(); ++i)
			{
				GraphPropertyStruct& gps21Diff = p21DiffSet[i];

				//	Index of this vertex from graph2 ordered list
				int graph2Index = gps21Diff.index_;

				//	Get the vertex.
				Vertex graph2Vertex = orderedList2[graph2Index];

				//	Create a new vertex and set the properties.
				Vertex v = boost::add_vertex(objectData.g_);
				objectData.propertyMap_ = boost::get(boost::vertex_index2, objectData.g_);

				GraphPropertyStruct& gps = boost::get(boost::vertex_index2, objectData.g_, v);
				gps.name_ = gps21Diff.name_;
				gps.value_ = gps21Diff.value_;
				gps.dataType_ = gps21Diff.dataType_;
				gps.addUpdateNodeMappingType(MapColorSchemes::SCHEME_NODE_SET_TYPES, SetTypes::SET_1_NOT_2);
				orderedList1.push_back(v);
				gps.index_ = orderedList1.size();
				objectData.propertyMap_[v] = gps;

				bool bFoundParent = false;
				//	See if the node has in_edges from its graph
				boost::graph_traits<DiGraph>::in_edge_iterator outI, outEnd;
				int edgeCount = 0;

				do
				{
					//	Now find the equivalent parent in this list if there is one.
					for(boost::tie(outI, outEnd) = boost::in_edges(graph2Vertex, graphIn.graphDataSet_.g_);
						outI != outEnd; ++outI)
					{
						Vertex src = boost::source(*outI, graphIn.graphDataSet_.g_);

						//	Get the property for the source
						GraphPropertyStruct& gps2 = boost::get(boost::vertex_index2, graphIn.graphDataSet_.g_, src);

						//	We now have the index of the parent node in its graph.  
						int src2Index = gps2.index_;

						//	Do a search of the lookup mappings to see if there is a common parent.
						map<int,int>::iterator indexIter = intersectionMapping.find(src2Index);

						if(indexIter != intersectionMapping.end())
						{
							int src1Index = indexIter->second;

							Vertex parent = orderedList1[src1Index];
							boost::add_edge(parent, v, objectData.g_);
							bFoundParent = true;
						}
						else
						{
							//	Node not found in graph 1.  Check graph2's parent node then.
							graph2Vertex = src;
							break;
						}
					}
				}while(!bFoundParent);
			}
		}
		else
		{
			std::cout << "Graphs are not properly defined.  Exiting." << std::endl;
			return false;
		}
		return true;
	}

private:

	GraphPropertyWriter<Graph, Index2AttributeMap> graphDataSet_;
	rapidjson::Document document_;

	rapidjson::Value::ConstMemberIterator mapJSONToDAG(rapidjson::Document& doc,
		rapidjson::Value::ConstMemberIterator iterIn,
		std::vector< ParentNodeMap >& parentListMap,
		std::map<JSONNode, int>& objNodeMap,
		int& vNestLevel, MapColorSchemes colorScheme, int colorValue)
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
			//	of the visited nodes.
			if(!bAlreadyExists)
			{
				std::cout << "Adding vertex " << objNodeMap.size() << ": " << iterIn->name.GetString() << std::endl;
				objNodeMap.insert(std::make_pair(iterIn, objNodeMap.size()));

				//	Create the vertex and set the nodeNames and nodeIndices properties.
				Vertex v = boost::add_vertex(graphDataSet_.g_);
				GraphPropertyStruct& gps = boost::get(boost::vertex_index2, graphDataSet_.g_, v);
				gps.name_ = iterIn->name.GetString();
				gps.index_ = objNodeMap.size();
				gps.shapeType_ = MapShapeTypes::SHAPE_OBJECT;
				gps.dataType_ = iterIn->value.GetType();
				gps.rank_ = vNestLevel;
				gps.addUpdateNodeMappingType(colorScheme, colorValue);
				if(bIsRootNode)
				{
					//	Root node has no incoming parent edges 
					ParentNodeMap parentNode(JSONToBoostMap(JSONNode(iterIn), v), graphDataSet_.g_);
					parentListMap.push_back(parentNode);
				}
				else
				{
					Vertex parent = vNestLevel > 0 ? parentListMap[vNestLevel - 1].second() : NULL;
					ParentNodeMap parentNode(JSONToBoostMap(JSONNode(iterIn), v), graphDataSet_.g_);
					parentListMap.push_back(parentNode);
					parentListMap[vNestLevel - 1].addChildNode(v);

					//	Create an edge from the parent to this vertex.
					boost::add_edge(parent, v, graphDataSet_.g_);
					std::cout << "Adding edge e(" << graphDataSet_.propertyMap_[parent].name_ << "," << graphDataSet_.propertyMap_[v].name_ << ")" << std::endl;
				}

				if(iterIn->value.MemberBegin() != iterIn->value.MemberEnd())
				{
					vNestLevel++;
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
				Vertex parent = vNestLevel > 0 ? parentListMap[vNestLevel - 1].second() : NULL;

				std::cout << "Adding vertex " << objNodeMap.size() << ": " << iterIn->name.GetString() << std::endl;
				objNodeMap.insert(std::make_pair(JSONNode(iterIn), objNodeMap.size()));
				Vertex v = boost::add_vertex(graphDataSet_.g_);

				GraphPropertyStruct& gps = boost::get(boost::vertex_index2, graphDataSet_.g_, v);
				gps.index_ = objNodeMap.size();
				gps.shapeType_ = MapShapeTypes::SHAPE_ARRAY;
				gps.dataType_ = iterIn->value.GetType();
				gps.name_ = iterIn->name.GetString();
				gps.rank_ = vNestLevel;
				gps.addUpdateNodeMappingType(colorScheme, colorValue);

				boost::add_edge(parent, v, graphDataSet_.g_);

				ParentNodeMap parentNode(JSONToBoostMap(JSONNode(iterIn), v), graphDataSet_.g_);
				parentListMap.push_back(parentNode);
				parentListMap[vNestLevel - 1].addChildNode(v);

				std::cout << "Adding edge e(" << graphDataSet_.propertyMap_[parent].name_ << "," << graphDataSet_.propertyMap_[v].name_ << ")" << std::endl;

				vNestLevel++;
			}

			//	Get the size of the array.
			const rapidjson::SizeType size = iterIn->value.Size();

			rapidjson::Value::ConstValueIterator valIter = iterIn->value.Begin();
			while(valIter != iterIn->value.End())
			{
				Vertex parent = vNestLevel > 0 ? parentListMap[vNestLevel - 1].second() : NULL;

				std::cout << "Adding vertex " << objNodeMap.size() << ": " << valIter->GetString() << std::endl;
				string val = valIter->GetString();

				objNodeMap.insert(std::make_pair(valIter, objNodeMap.size()));

				Vertex v = boost::add_vertex(graphDataSet_.g_);

				GraphPropertyStruct& gps = boost::get(boost::vertex_index2, graphDataSet_.g_, v);
				gps.index_ = objNodeMap.size();
				gps.shapeType_ = MapShapeTypes::SHAPE_VALUE;
				gps.dataType_ = valIter->GetType();
				gps.name_ = val;
				gps.value_ = val;
				gps.rank_ = vNestLevel;

				gps.addUpdateNodeMappingType(colorScheme, colorValue);

				//	Add parent edge to element
				boost::add_edge(parent, v, graphDataSet_.g_);
				ParentNodeMap parentNode(JSONToBoostMap(JSONNode(iterIn), v), graphDataSet_.g_);
				parentListMap.push_back(parentNode);

				parentListMap[vNestLevel - 1].addChildNode(v);

				std::cout << "Adding edge e(" << graphDataSet_.propertyMap_[parent].name_ << "," << graphDataSet_.propertyMap_[v].name_ << ")" << std::endl;

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
			return iterIn;
		}
		else
		{
			if(mapIter == objNodeMap.end())
			{
				std::cout << "Adding vertex " << objNodeMap.size() << ": " << iterIn->name.GetString() << std::endl;
				objNodeMap.insert(std::make_pair(iterIn, objNodeMap.size()));

				Vertex v = boost::add_vertex(graphDataSet_.g_);
				GraphPropertyStruct& gps = boost::get(boost::vertex_index2, graphDataSet_.g_, v);
				gps.index_ = objNodeMap.size();
				gps.shapeType_ = MapShapeTypes::SHAPE_VALUE;
				gps.dataType_ = iterIn->value.GetType();
				gps.name_ = iterIn->name.GetString();
				gps.value_ = iterIn->value.GetString();
				gps.rank_ = vNestLevel;
				gps.addUpdateNodeMappingType(colorScheme, colorValue);

				Vertex parent = vNestLevel > 0 ? parentListMap[vNestLevel - 1].second() : NULL;
				boost::add_edge(parent, v, graphDataSet_.g_);

				parentListMap[vNestLevel - 1].addChildNode(v);

				std::cout << "Adding edge e(" << graphDataSet_.propertyMap_[parent].name_ << "," << graphDataSet_.propertyMap_[v].name_ << ")" << std::endl;
			}
		}

		return ++iterIn;
	}

	//	Takes an input JSON document and maps to a DFS ordered directed acyclic graph.
	//	Does not create edges between sibling nodes.
	//	Sets the fill color to a value for a given color scheme.
	rapidjson::Value::ConstMemberIterator mapJSONToGraph(rapidjson::Document& doc,
		rapidjson::Value::ConstMemberIterator iterIn,
		std::vector< JSONToBoostMap >& parentListMap,
		std::map<JSONNode, int>& objNodeMap,
		int& vNestLevel, MapColorSchemes colorScheme, int colorValue)
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
				std::cout << "Adding vertex " << objNodeMap.size() << ": " << iterIn->name.GetString() << std::endl;
				objNodeMap.insert(std::make_pair(iterIn, objNodeMap.size()));

				//	Create the vertex and set the nodeNames and nodeIndices properties.
				Vertex v = boost::add_vertex(graphDataSet_.g_);
				GraphPropertyStruct& gps = boost::get(boost::vertex_index2, graphDataSet_.g_, v);
				gps.name_ = iterIn->name.GetString();
				gps.index_ = objNodeMap.size();
				gps.shapeType_ = MapShapeTypes::SHAPE_OBJECT;
				gps.dataType_ = iterIn->value.GetType();
				gps.addUpdateNodeMappingType(colorScheme, colorValue);
				if(bIsRootNode)
				{
					parentListMap.push_back(JSONToBoostMap(JSONNode(iterIn), v));
				}
				else
				{
					parentListMap.push_back(JSONToBoostMap(JSONNode(iterIn), v));
					Vertex parent = parentListMap[vNestLevel - 1].second();
					boost::add_edge(parent, v, graphDataSet_.g_);
					std::cout << "Adding edge e(" << graphDataSet_.propertyMap_[parent].name_ << "," << graphDataSet_.propertyMap_[v].name_ << ")" << std::endl;
				}

				if(iterIn->value.MemberBegin() != iterIn->value.MemberEnd())
				{
					vNestLevel++;
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
				std::cout << "Adding vertex " << objNodeMap.size() << ": " << iterIn->name.GetString() << std::endl;
				objNodeMap.insert(std::make_pair(JSONNode(iterIn), objNodeMap.size()));
				Vertex v = boost::add_vertex(graphDataSet_.g_);

				GraphPropertyStruct& gps = boost::get(boost::vertex_index2, graphDataSet_.g_, v);
				gps.index_ = objNodeMap.size();
				gps.shapeType_ = MapShapeTypes::SHAPE_ARRAY;
				gps.dataType_ = iterIn->value.GetType();
				gps.name_ = iterIn->name.GetString();
				gps.addUpdateNodeMappingType(colorScheme, colorValue);

				Vertex parent = parentListMap[vNestLevel - 1].second();
				boost::add_edge(parent, v, graphDataSet_.g_);

				parentListMap.push_back(JSONToBoostMap(JSONNode(iterIn), v));

				std::cout << "Adding edge e(" << graphDataSet_.propertyMap_[parent].name_ << "," << graphDataSet_.propertyMap_[v].name_ << ")" << std::endl;

				vNestLevel++;
			}

			//	Get the size of the array.
			const rapidjson::SizeType size = iterIn->value.Size();

			rapidjson::Value::ConstValueIterator valIter = iterIn->value.Begin();
			while(valIter != iterIn->value.End())
			{
				std::cout << "Adding vertex " << objNodeMap.size() << ": " << valIter->GetString() << std::endl;
				string val = valIter->GetString();

				objNodeMap.insert(std::make_pair(valIter, objNodeMap.size()));

				Vertex v = boost::add_vertex(graphDataSet_.g_);

				GraphPropertyStruct& gps = boost::get(boost::vertex_index2, graphDataSet_.g_, v);
				gps.index_ = objNodeMap.size();
				gps.shapeType_ = MapShapeTypes::SHAPE_VALUE;
				gps.dataType_ = valIter->GetType();
				gps.name_ = val;
				gps.value_ = val;
				gps.addUpdateNodeMappingType(colorScheme, colorValue);

				//	Add parent edge to element
				Vertex parent = parentListMap[vNestLevel - 1].second();
				boost::add_edge(parent, v, graphDataSet_.g_);

				std::cout << "Adding edge e(" << graphDataSet_.propertyMap_[parent].name_ << "," << graphDataSet_.propertyMap_[v].name_ << ")" << std::endl;

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
			return iterIn;
		}
		else
		{
			if(mapIter == objNodeMap.end())
			{
				std::cout << "Adding vertex " << objNodeMap.size() << ": " << iterIn->name.GetString() << std::endl;
				objNodeMap.insert(std::make_pair(iterIn, objNodeMap.size()));

				Vertex v = boost::add_vertex(graphDataSet_.g_);
				GraphPropertyStruct& gps = boost::get(boost::vertex_index2, graphDataSet_.g_, v);
				gps.index_ = objNodeMap.size();
				gps.shapeType_ = MapShapeTypes::SHAPE_VALUE;
				gps.dataType_ = iterIn->value.GetType();
				gps.name_ = iterIn->name.GetString();
				gps.value_ = iterIn->value.GetString();
				gps.addUpdateNodeMappingType(colorScheme, colorValue);

				Vertex parent = parentListMap[vNestLevel - 1].second();
				boost::add_edge(parent, v, graphDataSet_.g_);
				std::cout << "Adding edge e(" << graphDataSet_.propertyMap_[parent].name_ << "," << graphDataSet_.propertyMap_[v].name_ << ")" << std::endl;
			}
		}

		return ++iterIn;
	}

	void writeGraphOutputFile(const std::string& dotInFile, const std::string& pngFileOut)
	{
		graphDataSet_.writeGraphOutputFile(dotInFile, pngFileOut);
	}
};

template <class G, class N, class I>
bool diff(JSONBoostGraphMapping<G,N,I>& graph1, JSONBoostGraphMapping<G,N,I>& graph2, 
	OutputGraphMap& outputGraphMap)
{
	return false;
}



#endif