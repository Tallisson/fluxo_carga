#ifndef JACOBIAN_H_
#define JACOBIAN_H_

#include "graph.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

#include <armadillo>

using namespace arma;

namespace ns3
{

class Jacobian : public Object
{
public:
	Jacobian();
	virtual ~Jacobian();

	static TypeId GetTypeId();

	mat GetMatrix () const;
	void SetMatrix (uint64_t numLines, uint64_t numCols);

	void SetJ1(uint64_t numRows, uint64_t numCols);
	mat GetJ1 () const;

	void SetJ2(uint64_t numRows, uint64_t numCols);
	mat GetJ2 () const;

	void SetJ3(uint64_t numRows, uint64_t numCols);
	mat GetJ3 () const;

	void SetJ4(uint64_t numRows, uint64_t numCols);
	mat GetJ4 () const;

	mat CalcJac(Ptr<Graph> graph);
	vec SolveSys (vec b);

	void Zeros (void);
	void Zeros (uint32_t numRows, uint32_t numCols);

	mat GetJqv(void);

private:
	mat m_matrix;
	mat m_j1;
	mat m_j2;
	mat m_j3;
	mat m_j4;

	void CalcDPk(Ptr<Graph> graph);
	void CalcDQk(Ptr<Graph> graph);
	void CalcJDQ(Ptr<Graph> graph);
};

}

#endif  // JACOBIAN_H_
