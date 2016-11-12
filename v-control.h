#ifndef V_CONTROL_H_
#define V_CONTROL_H_

#include "graph.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

#include <armadillo>

namespace ns3
{

class VControl : public Object
{
public:
	static TypeId GetTypeId(void);

	Ptr<Bus> MaxDsv (Ptr<Graph> graph);

	uint32_t MaxV (Ptr<Graph> graph, arma::vec vt, Ptr<Bus> modBus);

	virtual bool DoControl (arma::mat jqv, Ptr<Graph> graph) = 0;


};

}

#endif /* V_CONTROL_H_ */
