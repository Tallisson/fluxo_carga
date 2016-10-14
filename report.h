#ifndef REPORT_H_
#define REPORT_H_

#include "graph.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

#include <map>
#include <string>
#include <vector>

namespace ns3
{

struct TBus {
public:
	DBus_t m_data;
	double m_crt;
};

typedef TBus TBus_t;

class Report : public Object {
public:
	Report(void);

	~Report(void);

	static TypeId GetTypeId(void);

	void StoreData(Ptr<Graph> graph, double base);

	void StoreL(double l);

	void Save(void);

	void SetFileLTotal(std::string fileLTotal);

	void SetFileCrt(std::string crt);

	void SetFileQg(std::string fileQg);

	void SetFileV(std::string fileV);

	void SetFileState(std::string fileState);
private:
	std::string m_fileV;

	std::string m_fileQg;

	std::string m_fileState;

	std::string m_fileL_km;

	std::string m_fileLTotal;

	std::string m_fileCrt;

	std::map<uint32_t, std::vector<TBus_t> > m_states;

	std::vector<double> m_l;
};

}

#endif /*REPORT_H_*/
