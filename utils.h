#ifndef UTILS_H_
#define UTILS_H_

#include "bus.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <boost/cstdint.hpp>
#include <boost/regex.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/container/vector.hpp>
#include <string>
#include <vector>
#include <map>

namespace ns3 {

struct sts
{
public:
  bool m_hasTape;
  double m_baseMVA;

  uint32_t m_nb;
  uint32_t m_npv;
  uint32_t m_npq;
  uint32_t m_nc;
  uint32_t m_ntap;
  uint32_t m_posSlack;

  std::vector<DBus_t> buses;
  std::vector<DBranch_t> branches;

  uint32_t BusFromNex (uint32_t nex);
};
typedef sts Sts_t;


/*
 *
 */
class Utils : public Object
{
private:
	Sts_t m_s;
public:
	Utils();
	virtual ~Utils();

  Sts_t Read(std::string filename);

  static const double STEP = 0.00625;
  static const double TAPMIN_MIN = 0.75;
  static const double TAPMAX_MAX = 1.25;

  static const double TAPMIN = 0.90;
  static const double TAPMAX = 1.10;

  static const double VMIN_MIN = 0.65;
  static const double VMAX_MAX = 1.35;

  static const double VMIN = 0.90;
  static const double VMAX = 1.10;

  static TypeId GetTypeId(void);

  Sts_t GetSts(void) const;

  std::string Format (double v, uint32_t p, uint32_t maxS);
};

}
#endif /* UTILS_H_ */
