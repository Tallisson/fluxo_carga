#ifndef IO_UTILS_H_
#define IO_UTILS_H_

#include "ns3/object.h"
#include "ns3/type-id.h"
#include "ns3/graph.h"

#include <vector>
#include <string>

namespace ns3
{

struct data
{
public:
	std::string m_label;
	double m_value;
};

typedef data Data_t;

class IOUtils : public Object {
public:
	IOUtils();
	virtual ~IOUtils();

	static TypeId GetTypeId();

	static std::string Print(std::vector<Data_t> values, uint8_t numSpace = 8);

	static std::string Format(double v, uint32_t p = 5, uint32_t maxS = 10);

	static std::string Format (std::string v, uint32_t maxS = 9);

	static std::string WriteCdf(Ptr<Graph> graph, double getBase);

	static std::string WriteCnz(Ptr<Graph> graph, double getBase);

	static void WriteFile (std::string filename, uint32_t choose, Ptr<Graph> graph, double base = 100);

	static const std::string SEP;
};

}

#endif /* BRANCH_H_ */
