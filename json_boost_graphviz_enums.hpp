#ifndef _JSON_BOOST_GRAPHVIZ_ENUMS_H_
#define _JSON_BOOST_GRAPHVIZ_ENUMS_H_

#include "../rapidjson/include/rapidjson.h"

//	Map the parse type to a display type -->
//	Object --> Ellipse
//  Primitive --> Rectangle
enum MapShapeTypes
{
	SHAPE_NONE = 0,
	SHAPE_OBJECT = 1,
	SHAPE_ARRAY = 2,
	SHAPE_VALUE = 3,
	SHAPE_FIELD_COMPARE = 4
};

static const std::string shapeValues [] =
{
	"none",
	"component",
	"box",
	"ellipse",
	"doublecircle"
};

enum ShapeSizeTypes
{
	SHAPE_VERY_SMALL = 0,
	SHAPE_SMALL = 1,
	SHAPE_NORMAL = 2,
	SHAPE_LARGE = 3,
	SHAPE_VERY_LARGE = 4
};

enum MatchTypes
{
	NO_MATCH = 0,
	MATCH = 1
};

enum SetTypes
{
	SET_EMPTY = 1,				//	Empty set
	SET_INTERSECTION = 2,		//	Set intersection
	SET_UNION = 3,				//	Union of two sets
	SET_1_NOT_2 = 3,			//	In set 1 but not in set 2
	SET_2_NOT_1 = 4				//	In set 2 but not in set 1
};

static const std::string SetTypeStrings[] =
{
	"0",
	"1",
	"2",
	"3"
};

static const std::string ShapeSizeStrings [] =
{
	"0.1",
	"0.2",
	"1",
	"2",
	"4"
};

enum MapColorSchemes
{
	SCHEME_NODE_TYPES = 1,			//	Parent or Leaf
	SCHEME_NODE_VALUE_TYPES = 2,	//  If leaf, what is the data type of the value (string, int, etc.)
	SCHEME_NODE_SET_TYPES = 3,		//  Node in both graphs; Node in Graph1 but not in Graph2; Node in Graph2 but not in Graph1
	SCHEMA_NODE_EQUIV_TYPES = 4    //	Nodes are equivalent; nodes are not equivalent; indeterminate
};

enum MapColorTypes
{
	COLOR_NONE = 0,
	COLOR_GREEN = 1,
	COLOR_RED = 2,
	COLOR_BLUE = 3,
	COLOR_WHITE = 4,
	COLOR_BLACK = 5,
	COLOR_GREY = 6,
	COLOR_YELLOW = 7,
	COLOR_GREEN_YELLOW = 8,
	COLOR_LIGHT_RED = 9,
	COLOR_LIGHT_BLUE = 10,
	COLOR_ROYAL_BLUE = 11
};

static const std::string colorValues [] =
{
	"antiquewhite",
	"green",
	"red",
	"deepskyblue",
	"white",
	"black",
	"grey",
	"yellow",
	"greenyellow",
	"firebrick1",
	"aquamarine",
	"royalblue3"
};

typedef rapidjson::Type JSONDataType;

static std::string kTypeNames[] = 
{ 
	"Null", 
	"False", 
	"True", 
	"Object", 
	"Array", 
	"String", 
	"Number" 
};


#endif