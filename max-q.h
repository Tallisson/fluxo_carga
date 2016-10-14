#ifndef MAX_Q_H_
#define MAX_Q_H_

#include "control.h"
#include "graph.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

namespace ns3
{

class MaxQ : public Control
{
private:
	void UpdateOrd(Ptr<Graph> graph);

	Ptr<Bus> MaxDsv(Ptr<Graph> graph);

public:
	MaxQ();

	virtual ~MaxQ();

	static TypeId GetTypeId(void);

	virtual bool DoControl (Ptr<Graph> graph);
	virtual bool DoRestore (Ptr<Graph> graph);

  static const double LIMIAR = 0.01;

  Ptr<Jacobian> m_jac;
};

}

#endif /* MaxQ_H_ */
