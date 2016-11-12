#ifndef MAX_Q_H_
#define MAX_Q_H_

#include "v-control.h"
#include "graph.h"
#include "jacobian.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

namespace ns3
{

class MaxQ : public VControl
{
public:
	MaxQ();

	virtual ~MaxQ();

	static TypeId GetTypeId(void);

	virtual bool DoControl (mat jqv, Ptr<Graph> graph);

	//uint32_t MaxV (Ptr<Graph> graph, vec q, Ptr<Bus> modBus);

  static const double LIMIAR = 0.01;
};

}

#endif /* MaxQ_H_ */
