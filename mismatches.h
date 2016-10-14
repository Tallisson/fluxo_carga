#ifndef MISMATCHES_H_
#define MISMATCHES_H_

#include "graph.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

#include <armadillo>

using namespace arma;

namespace ns3
{

class Mismatch : public Object
{
private:
	vec m_mis;

	void CalcPkB(Ptr<Graph> graph);
	void CalcQkB(Ptr<Graph> graph);
public:
	Mismatch();
	virtual ~Mismatch();

	static TypeId GetTypeId();

	vec CalcMismatches(Ptr<Graph> graph);

	vec GetMis(void);
};

}

#endif  // MISMATCHES_H_
