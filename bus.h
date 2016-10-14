#ifndef BUS_H_
#define BUS_H_

#include "branch.h"

#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

#include <string>
#include <vector>

namespace ns3
{

/*
 * Struct to store buses variables
 */
struct bus
{
public:
  uint32_t m_nex;
  uint32_t m_nin;
  std::string m_nome;
  double m_area;
  double m_zona;
  uint32_t m_tipo;
  double m_v;
  double m_ang;
  double m_pc;
  double m_qc;
  double m_pg;
  double m_qg;
  double m_base_kV;
  double m_vg_o;
  double m_qgmax;
  double m_qgmin;
  double m_vmax;
  double m_vmin;
  double m_gsh;
  double m_bsh;
  uint32_t m_ctrl_rem;
  uint32_t m_ordPV;
  uint32_t m_posPV;
  uint32_t m_ordPQ;
  uint32_t m_posPQ;
  uint32_t m_ord;
};
typedef bus DBus_t;

class Bus : public Object {
public:
  enum Type
  {
    SLACK = 3,
    GENERATION = 2,
    LOAD = 0,
		LOSS_CONTROL_REACT = 4
  };

  enum TapType
	{
  	IMP = 0,
		TAP = 1
	};

  enum Violation
	{
  	NORMAL = 0,
  	MIN_VOLTAGE_VIOLATION = 1,
		MAX_VOLTAGE_VIOLATION = 2
	};

  static const double MIN_VOLTAGE_GR = 1.00;
  static const double MIN_VOLTAGE_IEEE = 1.10;
  static const double MIN_VOLTAGE_ONS = 0.95;
  static const double MAX_VOLTAGE_ONS = 1.05;
  static const double MAX_VOLTAGE_IEEE = 0.90;
  static const double INCREMENT_STEP = 0.01;

  Bus();
  virtual ~Bus();

  static TypeId GetTypeId();

  DBus_t GetBus(void) const;
  void SetBus(DBus_t bus);

  Type GetType(void);
  void SetType(enum Type);

  void AddBranch(Ptr<Branch> branch, Ptr<Bus> neighbor);

  std::vector<Ptr<Branch> > GetBranches() const;
  std::vector<Ptr<Bus> > GetNeighbors() const;

  double CalcPg(void);
  double CalcQg (void);

  void SetTap(TapType tap);
  TapType GetTap(void) const;

  void Print(void);

  void SetCrt(double crt);
  double GetCrt(void) const;

  double CalcDsv(void);

  double GetDsv(void);

  Bus::Violation GetStatus (void);
private:
	double m_aCalc;
	double m_vCalc;
	double m_qgCalc;
	double m_pgCalc;
	double m_dsv;

	TapType m_tap;

	DBus_t m_bus;
	std::vector<Ptr<Branch> > m_branches;
	std::vector<Ptr<Bus> > m_neighbors;

	Type m_type;

	double m_crt;

	Violation m_status;
};
}
#endif /* BUS_H_ */
