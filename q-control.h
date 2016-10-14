#ifndef QCONTROL_H_
#define QCONTROL_H_

#include "control.h"
#include "graph.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

namespace ns3
{

class QControl : public Control
{
private:
	void UpdateOrd(Ptr<Graph> graph);

public:
	QControl();

	virtual ~QControl();

	static TypeId GetTypeId(void);

	virtual bool DoControl (Ptr<Graph> graph);
	virtual bool DoRestore (Ptr<Graph> graph);
};

}

#endif /* QCONTROL_H_ */
