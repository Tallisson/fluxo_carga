#ifndef BRANCH_H_
#define BRANCH_H_

#include "ns3/double.h"
#include "ns3/ptr.h"
#include "ns3/object.h"
#include "ns3/type-id.h"

#include <map>
#include <string>

namespace ns3
{

/*
 * Struct to store branches variables
 */
struct branch
{
public:
  uint32_t m_ni;
  uint32_t m_nf;
  uint32_t m_area;
  uint32_t m_zona;
  uint32_t m_circ;
  uint32_t m_tipo;
  uint32_t m_ordtap;
  uint32_t m_postap;
  double m_r;
  double m_x;
  double m_bsh;
  double m_line_rating1;
  double m_line_rating2;
  double m_line_rating3;
  double m_t_ctrl;
  double m_side;
  double m_tap;
  double m_def;
  double m_tapmin;
  double m_tapmax;
  double m_passo;
  double m_ctrl_min;
  double m_ctrl_max;
  double m_g;
  double m_b;
};
typedef branch DBranch_t;

class Branch : public Object {
private:
	DBranch_t m_branch;

	double m_p_km;
	double m_p_mk;

	double m_q_km;
	double m_q_mk;

	double m_p_km_L;
	double m_p_mk_L;
	double m_q_km_L;
	double m_q_mk_L;

	double m_l;
public:
	Branch();
	virtual ~Branch();

	static TypeId GetTypeId();

	void SetBranch(DBranch_t branch);
	DBranch_t GetBranch() const;

	/*double CalcPkm(DoubleValue vK, DoubleValue vM,
			DoubleValue aK, DoubleValue aM, double theta_km);
	double CalcPmk(DoubleValue vK, DoubleValue vM,
			DoubleValue aK, DoubleValue aM, double theta_km);

	double CalcQkm(DoubleValue vK, DoubleValue vM,
			DoubleValue aK, DoubleValue aM, double theta_km);
	double CalcQmk(DoubleValue vK, DoubleValue vM,
			DoubleValue aK, DoubleValue aM, double theta_km);*/

	double CalcPkmL(DoubleValue vK, DoubleValue vM,
			DoubleValue aK, DoubleValue aM);
	double CalcPmkL(DoubleValue vK, DoubleValue vM,
			DoubleValue aK, DoubleValue aM);

	double CalcQkmL(DoubleValue vK, DoubleValue vM,
			DoubleValue aK, DoubleValue aM);
	double CalcQmkL(DoubleValue vK, DoubleValue vM,
			DoubleValue aK, DoubleValue aM);

	double CalcL (DoubleValue vK, DoubleValue vM,
			DoubleValue aK, DoubleValue aM);

	void Print(void);
};

}

#endif /* BRANCH_H_ */
