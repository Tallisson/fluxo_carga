#ifndef GRAPH_H_
#define GRAPH_H_

#include "branch.h"
#include "bus.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

#include <string>
#include <map>

namespace ns3 {

class Graph : public Object
{
private:
	uint32_t m_numBus;
	uint32_t m_numBranch;

	uint32_t m_numPQ;
	uint32_t m_numPV;
	uint32_t m_numSlack;
	uint32_t m_posSlack;
	//uint32_t m_noSlack;

	std::map<uint32_t, Ptr<Bus> > m_buses;
	std::vector<Ptr<Branch> > m_branches;

public:
	Graph();
	virtual ~Graph();

	static TypeId GetTypeId (void);

	void AddBus(Ptr<Bus> bus);
	void Assoc(Ptr<Branch> branch);

	uint32_t GetPosSlack(void) const;
	void SetPosSlack (uint32_t posSlack);

	std::map<uint32_t, Ptr<Bus> > GetBuses() const;
	Ptr<Bus> GetBus(uint32_t idBus);

	std::vector<Ptr<Branch> > GetBranches() const;

	uint32_t GetNumBus() const;
	uint32_t GetNumBranch() const;

	uint32_t GetNumPQ() const;
	void SetNumPQ (uint32_t numPQ);

	void SetNumPV (uint32_t numPV);
	uint32_t GetNumPV() const;
};

}

#endif  // GRAPH_H_
