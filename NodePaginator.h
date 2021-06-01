//
//
//  Node paginator will provide the geometry model for rendering acyclic graphs
//	in pixel space
//  Inputs:  Page width in pixel units (W)
//			 Page height in pixel units (H)
//			 X Margin (w)
//			 Y Margin (h)
//			 Node X Margin (wn)
//			 Node Y Margin (hn)
//			 Directed Acyclic Graph (dag)
//	The font size of each node may be different as may the rendering styles (e.g. bold, italic).
//	Regardless we need to provide a flexible way to calculate the pixel width and height values
//  given a string length, a font spec, a style (italic, bold,  etc.), padding for each node in a graph.
//  NameConfig nameConfig = new NameConfig(name, fontName, styleFlags, fontSize)
//  ValueConfig valueConfig = new ValueConfig(valueType, *pFn(*) valueRenderFn)
//	void valueRenderFunction(int &nodeWidth, int& nodeHeight);

//	Each page will contain a certain node width and height.  The more nested the graph, the
//	larger its height will be.  The more sibling or array nodes, the larger the width will be.

//	A node is an entity with a name property and an optional value property.  The width of
//	a node will depend on the maximum of the rendering of the name or the value or some
//	representation of the value (e.g. bitmap/icon).  The name
//	and value may have different styles so it is the net width of rendering the name or the
//	value.  
//  W(node) = max(W(name), W(value)) for any node
//  W(name) = strlen(name)*fontFactorName + padName
//  W(value) = strlen(value)*fontFactorValue + padValue

//	For a graph of height "h", there are layout options to provide more meaningful
//	graph renderings.  
//	1)  Equal spacing -- page will be subdivided in X and Y directions and nodes
//	    will have proportional widths, heights,and spacing.
//  2)  Optimized -- will render graph on a portion of the page provided the rendering
//					 will be legible.
//  3)  Max depth -- page will be rendered with maximum graph level as specified.
//  4)  Filter on name -- regular expression to match a name
//  5)  Filter on value -- regular expression to match a value

//	Geometry model
//	Given a page of width W and height H, we want to provide a 2-D rendering model
//	for an arbitrary directed acyclic graph.  To render the model, it will be necessary
//  to expect that for any given level in the graph, no horizontal nodes will overlap, 
//  and also that no two nodes overlap vertically.

//  EQUAL SPACE MODEL
//	The padding between each node will be constant, and will be calculated
//  given the width and height of each cell in each vertical row of the graph.
//	Examples:  N = 1, assume we have a node width in pixels of w(node), and height h(node).
//  The padding in the X direction would follow the relation w + 2*padding = W
//  so that padding = (W - w)/2 
//  If N = 2, then w1 + w2 + padding*3 = W
//  so that padding = (W - w1 - w2)/3
//  And so on with general relation for any given vertical level
//  padding(x) = (W - sum(wi))/(N + 1)
//	Noting that each node may be of different width, and also height (less common though)
//  Similarly in the Y direction, with M = number of nodes.
//	padding(y) = (H - sum(hi))/(M + 1)
//
//	OPTIMIZED MODEL
//	This model will provide a meaningful rendering of the graph based on the number of
//	nodes in the vertical and horizontal levels, but will maintain the aspect ratio of
//	each node.  This differs from equal space model as equal space model will always
//	occupy the entire page.  This model will simply calculate the widths based on the
//	node layouts and provide standard options for rendering:
//	LEFT, XCENTER, RIGHT, TOP, YCENTER, BOTTOM
//
//	MAX DEPTH
//	This model will render a maximum vertical depth specified as an input parameter.
//
//	FILTER ON NAME
//	This model accepts a regex input and will provide a rendering with enlarged/highlighted 
//	nodes where filter matches regex.
//
//	FILTER ON VALUE
//	This model accepts a regex input and will provide a rendering with enlarged/highlighted
//	nodes where filter matches regex.
//
//	

