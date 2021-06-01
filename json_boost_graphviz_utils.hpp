#ifndef _JSON_BOOST_GRAPHVIZ_UTILS_H_
#define _JSON_BOOST_GRAPHVIZ_UTILS_H_

#include <vector>
#include <iosfwd>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graphviz.hpp>

#include "../rapidjson/include/document.h"

namespace
{
	static const std::string  emptySet = "<&#946;>";
}

struct GraphPropertyStruct
{
	GraphPropertyStruct():
		index_(0),
		name_(""),
		dataType_(JSONDataType::kNullType),
		shapeType_(MapShapeTypes::SHAPE_VALUE),
		sizeType_(ShapeSizeTypes::SHAPE_NORMAL),
		rank_(0)
	{
	}
	GraphPropertyStruct(const GraphPropertyStruct& gps):
		index_(gps.index_),
		name_(gps.name_),
		value_(gps.value_),
		shapeType_(gps.shapeType_),
		sizeType_(gps.sizeType_),
		dataType_(gps.dataType_),
		nodeTypeMap_(gps.nodeTypeMap_),
		rank_(gps.rank_)
	{
	}
	GraphPropertyStruct(int ind, const std::string& name,
		JSONDataType dt,
		MapShapeTypes st, 
		ShapeSizeTypes sizeType):
		index_(ind),
		dataType_(dt),
		name_(name),
		shapeType_(st),
		sizeType_(sizeType),
		rank_(0)
	{
	}

	//	Ensure < operator is weak order compliant.
	bool operator<(const GraphPropertyStruct& rhs) const
	{
		//	The index cannot determine equality in
		//	this context
//			if(index_ < rhs.index_)
//			{
//				return true;
//			}
//			else if(index_ > rhs.index_)
//			{
//				return false;
//			}
//			else
		if(name_ < rhs.name_)
		{
			return true;
		}
		else if(name_ > rhs.name_)
		{
			return false;
		}
		else
		{
			if(value_ < rhs.value_)
			{
				return true;
			}
			else if(value_ > rhs.value_)
			{
				return false;
			}
			else
			{
				if(rank_ < rhs.rank_)
				{
					return true;
				}
				else if(rank_ > rhs.rank_)
				{
					return false;
				}
				else
				{
					if(dataType_ < rhs.dataType_)
					{
						return true;
					}
					if(dataType_ > rhs.dataType_)
					{
						return false;
					}
				}
			}
		}
		return false;
	}

	GraphPropertyStruct& operator=(const GraphPropertyStruct& rhs)
	{
		name_ = rhs.name_;
		value_ = rhs.value_;
		index_ = rhs.index_;
		shapeType_ = rhs.shapeType_;
		sizeType_ = rhs.sizeType_;
		dataType_ = rhs.dataType_;
		nodeTypeMap_ = rhs.nodeTypeMap_;
		rank_ = rhs.rank_;
		return *this;
	}

	friend std::ostream& operator<<(std::ostream& os, const GraphPropertyStruct& gps)
	{
		os << gps.index_;
		return os;
	}

	bool addUpdateNodeMappingType(MapColorSchemes scheme, int value)
	{
		std::map<MapColorSchemes, int>::iterator mapIter = nodeTypeMap_.find(scheme);
		bool bFound = mapIter != nodeTypeMap_.end();

		if(!bFound)
		{
			auto retVal = nodeTypeMap_.insert(std::make_pair(scheme, value));	
			if(!retVal.second)
			{
				std::cout << "Error inserting into node type map" << std::endl;
			}
		}
		else
		{
			mapIter->second = value;
		}
		return bFound;
	}

	int index_;
	int rank_;
	std::string name_;
	MapShapeTypes shapeType_;
	ShapeSizeTypes sizeType_;

	//MapColorTypes colorType_;

	//	Color schemes used to display node characteristics
	std::map<MapColorSchemes, int> nodeTypeMap_;
	JSONDataType  dataType_;
	std::string value_;
};

bool operator==(const GraphPropertyStruct& lhs, const GraphPropertyStruct& rhs)
{
	return ((lhs.name_ == rhs.name_) &&
		(lhs.value_ == rhs.value_) &&
//		(lhs.colorType_ == rhs.colorType_)  &&
		(lhs.dataType_ == rhs.dataType_) );
}

struct GraphPropertyStructComparator
{
	GraphPropertyStructComparator()
	{
	}
	GraphPropertyStructComparator(GraphPropertyStruct& gps):
		gps_(gps)
	{
	}
	bool operator()(const GraphPropertyStruct& gps1,
					const GraphPropertyStruct& gps2) const
	{
		return gps1 < gps2;
	}

	bool operator==(const GraphPropertyStruct& gps2) const
	{
		return (gps_ == gps2);
	}

	GraphPropertyStruct gps_;
};

struct GraphPropertyStructIndexComparator
{
	GraphPropertyStructIndexComparator()
	{
	}
	bool operator()(const GraphPropertyStruct& gps1,
					const GraphPropertyStruct& gps2) const
	{
		return gps1.index_ < gps2.index_;
	}
};

typedef rapidjson::Value::ConstMemberIterator MEMBER_ITERATOR;
typedef rapidjson::Value::ConstValueIterator VALUE_ITERATOR;

typedef boost::property<boost::vertex_name_t, std::string, boost::property<boost::vertex_index2_t, GraphPropertyStruct> > VertexProperty;
typedef boost::property<boost::edge_name_t, std::string, boost::property<boost::edge_index_t, GraphPropertyStruct> > EdgeProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexProperty> DiGraphOld;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexProperty, EdgeProperty> DiGraph;

//	Define the vertex and edge descriptor
typedef boost::graph_traits<DiGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<DiGraph>::edge_descriptor Edge;

//	Define the name map and property map
typedef boost::property_map<DiGraph, boost::vertex_index_t>::type IndexMap;
typedef boost::property_map<DiGraph, boost::vertex_index2_t>::type GraphPropMap;
typedef boost::property_map<DiGraph, boost::vertex_name_t>::type NamePropMap;

typedef boost::graph_traits <DiGraph>::in_edge_iterator in_edge_iterator;
typedef boost::graph_traits<DiGraph>::out_edge_iterator out_edge_iterator;
namespace boost
{

  template <typename Graph, typename VertexWriter>
  inline void
  write_graphviz_yo(std::ostream& out, const Graph& g, VertexWriter vw
                 BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,boost::vertex_list_graph_tag))
  {
    boost::default_writer dw;
    boost::default_writer gw;
	boost::write_graphviz_yo(out, g, vw, dw, gw);
  }

  template <typename Graph, typename VertexPropertiesWriter,
            typename EdgePropertiesWriter, typename GraphPropertiesWriter>
  inline void
  write_graphviz_yo(std::ostream& out, const Graph& g,
                 VertexPropertiesWriter vpw,
                 EdgePropertiesWriter epw,
                 GraphPropertiesWriter gpw
                 BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,vertex_list_graph_tag))
  { 
	  write_graphviz_yo(out, g, vpw, epw, gpw, get(vertex_index, g)); 
  }

template <typename Graph, typename VertexPropertiesWriter,
        typename EdgePropertiesWriter, typename GraphPropertiesWriter,
        typename VertexID>
  inline void
  write_graphviz_yo
    (std::ostream& out, const Graph& g,
     VertexPropertiesWriter vpw,
     EdgePropertiesWriter epw,
     GraphPropertiesWriter gpw,
     VertexID vertex_id
     BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,vertex_list_graph_tag))
  {
    BOOST_CONCEPT_ASSERT((EdgeListGraphConcept<Graph>));

    typedef typename graph_traits<Graph>::directed_category cat_type;
    typedef graphviz_io_traits<cat_type> Traits;
    std::string name = "G";
    out << Traits::name() << " " << escape_dot_string(name) << " {" << std::endl;
	out << "rankdir=TB;" << std::endl;

    gpw(out); //print graph properties

	std::map<int, vector<Vertex> > rankVertexInfo;

    typename graph_traits<Graph>::vertex_iterator i, end;

    for(boost::tie(i,end) = vertices(g); i != end; ++i) 
	{
		Vertex vertex = get(vertex_id, *i);
		const GraphPropertyStruct& gps = boost::get(boost::vertex_index2, g, *i);
		std::cout << "Index: " << gps.index_ << " vertex: " << *i << " name: " << gps.name_ << " rank: " << gps.rank_ << std::endl;
		int rank = gps.rank_;
		std::map<int, vector<Vertex> >::iterator mapIter = rankVertexInfo.find(rank);
		if(mapIter == rankVertexInfo.end())
		{
			vector<Vertex> vertexVector;
			vertexVector.push_back(vertex);
			rankVertexInfo.insert(std::map<int, vector<Vertex> >::value_type(rank, vertexVector));
		}
		else
		{
			vector<Vertex>& rankVector = mapIter->second;
			rankVector.push_back(vertex);
		}

      out << escape_dot_string(vertex);
      vpw(out, *i); //print vertex attributes
      out << ";" << std::endl;
    }
    typename graph_traits<Graph>::edge_iterator ei, edge_end;
    for(boost::tie(ei, edge_end) = edges(g); ei != edge_end; ++ei) 
	{
      out << escape_dot_string(get(vertex_id, source(*ei, g))) << Traits::delimiter() << escape_dot_string(get(vertex_id, target(*ei, g))) << " ";
      epw(out, *ei); //print edge attributes
      out << ";" << std::endl;
    }

	//	Now add in the rank info.

	map<int, vector<Vertex> >::iterator rankMapIter =  rankVertexInfo.begin();

	while(rankMapIter != rankVertexInfo.end())
	{
		int rank = rankMapIter->first;

		vector<Vertex>& nodes = rankMapIter->second;

		out << "{ rank=same ";
		for(int j = 0; j < nodes.size(); ++j)
		{
			const GraphPropertyStruct& gps = boost::get(boost::vertex_index2, g, nodes[j]);
			if(j < nodes.size() - 1)
			{
				out << gps.index_ - 1 << ",";
			}
			else
			{
				out << gps.index_ - 1;
			}
		}
		out << "}" << std::endl;
		++rankMapIter;
	}
/*
graph 
{ 
		rankdir=LR;
		a -- { b c d }; b -- { c e }; c -- { e f }; d -- { f g }; e -- h; 
		f -- { h i j g }; g -- k; h -- { o l }; i -- { l m j }; j -- { m n k }; 
		k -- { n r }; l -- { o m }; m -- { o p n }; n -- { q r }; 
		o -- { s p }; p -- { s t q }; q -- { t r }; r -- t; s -- z; t -- z; 
		{ rank=same, b, c, d }
		{ rank=same, e, f, g }
		{ rank=same, h, i, j, k }
		{ rank=same, l, m, n }
		{ rank=same, o, p, q, r }
		{ rank=same, s, t }
}

*/
	
    out << "}" << std::endl;
  }

}

//template <class Graph, class NameMap, class Index2PropertyMap>
template <class Graph, class Index2PropertyMap>
struct GraphPropertyWriter
{
	GraphPropertyWriter(Graph& g):
		g_(g),
		propertyMap_()
	{
		propertyMap_ = boost::get(boost::vertex_index2, g_);
	}
	Graph& getGraph() const
	{
		return g_;
	}
	void createOutput(const std::string& outVizFile,
			const std::string& outPngFile)
	{
		//	Now write the graph out to graphviz format.  
		std::ofstream outputViz(outVizFile, std::ofstream::out);

//		boost::write_graphviz(outputViz, g_, name_graph_prop_writer<NameMap,Index2PropertyMap>(nameMap_,
//			propertyMap_));

		//	Override the write_graphviz with our function.
//		boost::write_graphviz(outputViz, g_, graph_prop_writer<Index2PropertyMap>(propertyMap_));
//		int x = 1;
		write_graphviz_yo(outputViz, g_, graph_prop_writer<Index2PropertyMap>(propertyMap_));

		//	Close out the graphviz file.
		outputViz.close();

		//	Now that we have it, write it to file.
		writeGraphOutputFile(outVizFile, outPngFile);
	}
	
	void writeGraphOutputFile(const std::string& dotInFile, const std::string& pngFileOut)
	{
		std::ifstream inFS(dotInFile.c_str());
		inFS.seekg(0, inFS.end);
		int bufLen = static_cast<int>(inFS.tellg());
		inFS.seekg(0, inFS.beg);

		std::string buffer(bufLen, 0);

		inFS.read(&buffer[0], bufLen);
		inFS.close();

		Agraph_t *G = NULL;
		if(bufLen > 0)
		{
			GVC_t *gvc = gvContext();

			G = agmemread(&buffer[0]);
			gvLayout (gvc, G, "dot");

			int rc = gvRenderFilename (gvc, G, "png", pngFileOut.c_str());
			if(rc != 0)
			{
				std::cout << "Error - rc = " << rc << std::endl;
			}
			gvFreeLayout(gvc, G);    

			agclose (G);    
			gvFreeContext(gvc);
		}

	}

	Graph& g_;
	Index2PropertyMap propertyMap_;
};

template <typename T1, typename T2> class MapVertexToVertex
{
public:
	MapVertexToVertex(T1 t1, T2 t2) :
		mapValue_(std::make_pair(t1, t2))
	{
	}

	T1& first()
	{
		return mapValue_.first;
	}
	T2& second()
	{
		return mapValue_.second;
	}
protected:
	std::pair<T1, T2> mapValue_;
};

struct JSONNode;
typedef MapVertexToVertex<JSONNode, Vertex> JSONToBoostMap;

typedef Vertex BoostNode;
typedef MapVertexToVertex<JSONNode, Vertex> JSONToBoostNodeMap;

template <typename VertexNameMap>
class dfs_set_visitor : public boost::default_dfs_visitor
{
public:
	dfs_set_visitor(std::set<Vertex>& vSet,
		std::set<int>& objIndexSet,
		std::set<GraphPropertyStruct>& objPropertySet):
		orderedOutputSet_(vSet),
		objectIndexSet_(objIndexSet),
		objectPropertySet_(objPropertySet)
	{
	}

	dfs_set_visitor()
	{
	}

	template<typename Vertex, typename Graph>
	void discover_vertex(Vertex v, const Graph& g)
	{
		orderedOutputSet_.insert(v);

		//	If of type object or array, then add to output list.
		auto gps = boost::get(boost::vertex_index2, g, v);
		objectPropertySet_.insert(gps);

		if((rapidjson::kObjectType == gps.dataType_) ||
			(rapidjson::kArrayType == gps.dataType_))
		{
			objectIndexSet_.insert(orderedOutputSet_.size() - 1);
		}
	}

	template<typename Vertex, typename Graph>
	void finish_vertex(Vertex u, const Graph& g)
	{
	}
private:
	std::set<Vertex>& orderedOutputSet_;
	std::set<int>& objectIndexSet_;
	std::set<GraphPropertyStruct>& objectPropertySet_;
};

template <typename VertexNameMap, typename PropertyMap>
class dfs_ordered_list_visitor : public boost::default_dfs_visitor
{
public:
	dfs_ordered_list_visitor(std::vector<Vertex>& vList,
		std::vector<int>& objectIndexList, 
		PropertyMap& propMap):
		orderedOutputList_(vList),
		objectIndexList_(objectIndexList),
		propertyMap_(propMap)
	{
	}

	dfs_ordered_list_visitor()
	{
	}

	template<typename Vertex, typename Graph>
	void discover_vertex(Vertex v, const Graph& g)
	{
		propertyMap_[v].index_ = orderedOutputList_.size();
		orderedOutputList_.push_back(v);

		//	If of type object or array, then add to output list.
		auto gps = boost::get(boost::vertex_index2, g, v);

		if((rapidjson::kObjectType == gps.dataType_) ||
			(rapidjson::kArrayType == gps.dataType_))
		{
			objectIndexList_.push_back(orderedOutputList_.size() - 1);

		}
	}
	template<typename Vertex, typename Graph>
	void finish_vertex(Vertex u, const Graph& g)
	{
	}
private:
	std::vector<Vertex>& orderedOutputList_;
	std::vector<int>& objectIndexList_;
	PropertyMap propertyMap_;
};

template <typename VertexNameMap>
class dfs_digraph_compare_visitor : public boost::default_dfs_visitor
{
public:
	dfs_digraph_compare_visitor(Vertex v):
		v_(v)
	{
	}

	template <class Graph>
	std::pair<Graph, Graph> compare(const Vertex& v,
		const Graph& g)
	{
		Graph leftSubgraph;
		Graph rightSubgraph;

		//	This function will return a graph that is a subgraph of
		//	of v when (v != v_).

		if(v != v_)
		{
			//	Starting with v, return a subgraph of v_, and v
			
		}
		return std::make_pair(leftSubgraph, rightSubgraph);
	}

	dfs_digraph_compare_visitor()
	{
	}

	template<typename Vertex, typename Graph>
	void discover_vertex(Vertex v, const Graph& g)
	{
		bIsSame_ = (v == v_);
		std::cout << "Visiting vertex " << v << std::endl;
		std::cout << "IsSame = " << bIsSame_;
	}
	template<typename Vertex, typename Graph>
	void finish_vertex(Vertex u, const Graph& g)
	{
		//	Works.
		//std::cout << "Finishing vertex " << u << std::endl;
	}
	Vertex& v()
	{
		return v_;
	}
private:
	//	May not be valid.  The comparator needs to be able
	//	to store a copy so it can
	Vertex v_;
	bool bIsSame_;
};

template <class Name>
class label_writer 
  {
  public:
    label_writer(Name _name) : name(_name) {}
    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const 
	{
      out << "[label=" << escape_dot_string(get(name, v)) << "]";
    }
  private:
    Name name;
  };

template <class Name, class GraphItemPropStruct>
class indexed_label_writer : public boost::label_writer<Name>
  {
  public:
	indexed_label_writer(Name n, GraphItemPropStruct& gps):
		gps_(gps),
		label_writer(gps.name_),
	{
	}
    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const 
	{
		std::stringstream tmp;
	
		GraphProperyStruct& gps = boost::get(gps_, v);

		if (!gps.value_.empty())
		{
			tmp << gps.index_ << ": " << gps.name_ << "/" << gps.value_;
		}
		else
		{
			stringstream tmpName;
			tmpName << gps.index_ << ": " << gps.name_ << "/{" << "&#8709" << "}";
			//tmpName << gps.name_ <<  "/{" << "&#8709" << "}";
			tmp <<  boost::escape_dot_string(tmpName.str());

		}
		out << "[";
		out << "label=" << boost::escape_dot_string(tmp.str());

		out << ",";
		out << "shape="<< shapeValues[gps.shapeType_];
		out << ",";
		out << "width=" << ShapeSizeStrings[gps.sizeType_] << ",height=" << ShapeSizeStrings[gps.sizeType_];
		out << ",";

		std::map<MapColorSchemes, int>::iterator iter = gps.nodeTypeMap_.begin();
		stringstream colorStream;
		while(iter != gps.nodeTypeMap_.end())
		{
			MapColorSchemes scheme = iter->first;
			int value = iter->second;

			if(scheme == SCHEME_NODE_SET_TYPES)
			{
				colorStream << "#";
				if(value == 2)
				{
					colorStream << "00" << "ff" << "00" << "7f";
				}
				else if(value == 3)
				{
					colorStream << "7f" << "00" << "00" << "ff";
				}
				else if(value == 4)
				{
					colorStream << "00" << "00" << "7f" << "ff";
				}
			}
			++iter;
		}

		string colorString = colorStream.str();
		if(colorString.length() > 0)
		{
			out << "fillcolor=\"" << colorString << "\"";
			out << ",";
			out << "color=black";
			out << ",";
			out << "style=filled";
		}
		out << "]";
	}
  private:
	GraphItemPropStruct& gps_;
};

template <class Name, class GraphItemPropStruct>
class name_graph_prop_writer : public boost::label_writer<Name>
{
  public:
		name_graph_prop_writer(Name n, GraphItemPropStruct gps):
			label_writer(n),
			name_(n),
			gps_(gps)
		{
		}
    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const 
	{
		std::stringstream tmp;
		GraphPropertyStruct& gps = boost::get(gps_, v);

		if (!gps.value_.empty())
		{
			//tmp << gps.value_;
			tmp << gps.name_ << "/" << gps.value_;
		}
		else
		{
			stringstream tmpName;
			tmpName << gps.name_ ;
			tmp <<  tmpName.str();
		}
		out << "[";
		out << "label=" << boost::escape_dot_string(tmp.str());
		out << ",";
		out << "shape="<< shapeValues[gps.shapeType_];
		out << ",";
		out << "width=" << ShapeSizeStrings[gps.sizeType_] << ",height=" << ShapeSizeStrings[gps.sizeType_];
		out << ",";
		std::map<MapColorSchemes, int>::iterator iter = gps.nodeTypeMap_.begin();
		stringstream colorStream;
		if(iter == gps.nodeTypeMap_.end())
		{
			bool bBad = true;
		}
		while(iter != gps.nodeTypeMap_.end())
		{
			MapColorSchemes scheme = iter->first;
			int value = iter->second;
			if(scheme == SCHEME_NODE_SET_TYPES)
			{
				colorStream << "#";
				if(value == 2)
				{
					colorStream << "00" << "ef" << "00" << "7f";
				}
				else if(value == 3)
				{
					colorStream << "fe" << "00" << "00" << "5f";
				}
				else if(value == 4)
				{
					colorStream << "00" << "00" << "ee" << "5f";
				}

			}
			++iter;
		}
		string colorString = colorStream.str();
		if(colorString.length() > 0)
		{
			out << "fillcolor=\"" << colorString << "\"";
			out << ",";
			out << "color=black";
			out << ",";
			out << "style=filled";
		}
		out << "]";
	
    }
  private:
	  Name name_;
	  GraphItemPropStruct gps_;
};

template <class GraphItemPropStruct>
class graph_prop_writer : public boost::label_writer<GraphItemPropStruct>
{
  public:
		graph_prop_writer(const GraphItemPropStruct& gps):
			label_writer(gps),
			gps_(gps)
		{
		}
    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const 
	{
		std::stringstream tmp;
		GraphPropertyStruct& gps = boost::get(gps_, v);

		if (!gps.value_.empty())
		{
			//tmp << gps.value_;
			tmp << gps.name_ << "/" << gps.value_;
		}
		else
		{
			stringstream tmpName;
			tmpName << gps.name_ ;
			tmp <<  tmpName.str();
		}
		out << "[";
		out << "label=" << boost::escape_dot_string(tmp.str());
		out << ",";
		out << "shape="<< shapeValues[gps.shapeType_];
		out << ",";
		out << "width=" << ShapeSizeStrings[gps.sizeType_] << ",height=" << ShapeSizeStrings[gps.sizeType_];
		out << ",";
		std::map<MapColorSchemes, int>::iterator iter = gps.nodeTypeMap_.begin();
		stringstream colorStream;
		while(iter != gps.nodeTypeMap_.end())
		{
			MapColorSchemes scheme = iter->first;
			int value = iter->second;
			if(scheme == SCHEME_NODE_SET_TYPES)
			{
				colorStream << "#";
				if(value == 2)
				{
					colorStream << "00" << "ef" << "00" << "7f";
				}
				else if(value == 3)
				{
					colorStream << "fe" << "00" << "00" << "5f";
				}
				else if(value == 4)
				{
					colorStream << "00" << "00" << "ee" << "5f";
				}

			}
			++iter;
		}
		string colorString = colorStream.str();
		if(colorString.length() > 0)
		{
			out << "fillcolor=\"" << colorString << "\"";
			out << ",";
			out << "color=black";
			out << ",";
			out << "style=filled";
		}
		out << "]";
	
    }
private:
	GraphItemPropStruct gps_;
};

template <class Name, class Index, class ShapeType, class ColorType>
class name_index_shape_writer : public boost::label_writer<Name>
{
  public:
	name_index_shape_writer(Name n, Index ind, ShapeType st, ColorType ct):
		label_writer(n),
		index(ind),
		myName(n),
		shapeType(st),
		colorType(ct)
	{
	}
    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const 
	{
		std::stringstream tmp;
		tmp << boost::get(index, v) << "-" << boost::get(myName, v);
		out << "[";
		out << "color=" << colorValues[boost::get(shapeType, v)];
		out << ",";
		out << "label=" << boost::escape_dot_string(tmp.str());
		out << ",";
		out << "shape="<< shapeValues[boost::get(colorType, v)];
		out << "]";
	
    }
  private:
	  Index index;
	  Name myName;
	  ShapeType shapeType;
	  ColorType colorType;
};
#endif