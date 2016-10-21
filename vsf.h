#ifndef VSF_H_
#define VSF_H_

#include "control.h"
#include "graph.h"
#include "jacobian.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

namespace ns3
{

class Vsf : public Control
{
private:
	void UpdateOrd(Ptr<Graph> graph);

	Ptr<Jacobian> m_jac;
public:
	Vsf();

	virtual ~Vsf();

	static TypeId GetTypeId(void);

	virtual bool DoControl (Ptr<Graph> graph);
	virtual bool DoRestore (Ptr<Graph> graph);

	Ptr<Bus> MaxDsv (Ptr<Graph> graph);

	void SetJqv(Ptr<Jacobian> jac);

	static const double LIMIAR = 0.01;
};

}

#endif /* VSF_H_ */
