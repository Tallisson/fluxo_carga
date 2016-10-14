#ifndef CONTROL_H_
#define CONTROL_H_

#include "graph.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

namespace ns3
{

class Control : public Object
{
public:
	/*Control();
	virtual ~Control();*/

	static TypeId GetTypeId(void);

	virtual bool DoControl (Ptr<Graph> graph) = 0;
	virtual bool DoRestore (Ptr<Graph> graph) = 0;
};

}

#endif /* CONTROL_H_ */
