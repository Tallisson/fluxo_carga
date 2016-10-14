#ifndef LOAD_FLOW_H_
#define LOAD_FLOW_H_

#include "graph.h"
#include "jacobian.h"
#include "mismatches.h"
#include "utils.h"
#include "control.h"
#include "q-control.h"
#include "report.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

#include <string>

namespace ns3
{

class LoadFlow : public Object
{
private:
	Ptr<Graph> m_graph;
	Ptr<Control> m_qControl;
	Ptr<Mismatch> m_mismatches;
	Ptr<Jacobian> m_jac;
	Ptr<Report> m_report;

	Sts_t m_sts;

	double m_precision;
	double m_totalL;
	uint32_t m_iter;
	uint32_t m_maxIter;

	vec m_losses;

	void InitX0(void);

	void InitJ(void);

	void CalcLosses(void);

	bool m_verbose;

	std::string m_dir;
	std::string m_file;
public:
	LoadFlow(void);
	virtual ~LoadFlow(void);

	static TypeId GetTypeId(void);

	Ptr<Graph> GetGraph(void) const;
	//void SetGraph(Ptr<Graph> graph);

	Ptr<Control> GetQControl(void) const;
	void SetQControl(Ptr<Control> qControl);

	void Prepare(std::string cdf);

	void Execute();

	void SetPrecision(double precision);

	void SetDir(std::string dir);

	void SetFile(std::string file);
};

}

#endif /* LOAD_FLOW */