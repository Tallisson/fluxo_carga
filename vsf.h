#ifndef VSF_H_
#define VSF_H_

#include "v-control.h"
#include "graph.h"
#include "jacobian.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

#include <armadillo>

namespace ns3
{

class Vsf : public VControl
{
public:
	Vsf();

	virtual ~Vsf();

	static TypeId GetTypeId(void);

	virtual bool DoControl (mat jac, Ptr<Graph> graph);

	static const double LIMIAR = 0.01;
};

}

#endif /* VSF_H_ */
