//	Bayesian graph project
//
//
//	This project will provide a Bayesian model using Boost Graph
//	library.  It is educational in nature and will leverage the
//	rich feature set of Boost to perform some complex calculations
//	on Bayesian networks.

//	Bayes theorem:
//	P(A|B) = (1/P(B))*P(A)*P(B|A)
//	The probability of A happening given that B has in fact occurred
//	is equal to the product of probability of A happening P(A) times
//  the probability of B given that A has occurred P(B|A) divided by
//  the probability of B

//	Methodology -- construct a graph that will represent the pdf for a measurable
//  chain of events.
//  Application of the function on any non-leaf node will provide the cumulative
//  probability at that node as the aggregration/integration of the cdf, given evidence
//  provided at that node.