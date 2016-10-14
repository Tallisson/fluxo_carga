/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 */

#include "branch.h"

#include "ns3/log.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"

#include <math.h>

namespace ns3
{

NS_LOG_COMPONENT_DEFINE ("Branch");
NS_OBJECT_ENSURE_REGISTERED (Branch);

Branch::Branch ()
{
}

Branch::~Branch ()
{
}

TypeId
Branch::GetTypeId (void)
{
	static TypeId tid = TypeId ("ns3::Branch")
		      .SetParent<Object> ()
		      .AddConstructor<Branch> ()
					.AddAttribute("Qkm",
	            "Reactive Power km",
	            DoubleValue (0), /* TODO */
	            MakeDoubleAccessor (&Branch::m_p_km),
	            MakeDoubleChecker<double> ())
					.AddAttribute("Qmk",
					    "Reactive Power mk",
					    DoubleValue (0), /* TODO */
					    MakeDoubleAccessor (&Branch::m_p_mk),
					    MakeDoubleChecker<double> ())
					.AddAttribute("Pkm",
					    "Active Power km",
							DoubleValue (0), /* TODO */
							MakeDoubleAccessor (&Branch::m_p_km),
							MakeDoubleChecker<double> ())
					.AddAttribute("Pmk",
							"Active Power mk",
							DoubleValue (0), /* TODO */
							MakeDoubleAccessor (&Branch::m_p_mk),
							MakeDoubleChecker<double> ())
	;

	return tid;
}

/*double
Branch::CalcPkm(DoubleValue vK, DoubleValue vM,
								DoubleValue aK, DoubleValue aM, double theta_km)
{
	m_p_km = m_branch.m_g * 1 / std::pow (m_branch.m_tap, 2) * std::pow (vK.Get (), 2) -
						1/m_branch.m_tap * vK.Get () * vM.Get () *
						( m_branch.m_g * cos (theta_km) + m_branch.m_b * sin (theta_km) );

	return m_p_km;
}

double
Branch::CalcPmk(DoubleValue vK, DoubleValue vM,
								DoubleValue aK, DoubleValue aM, double theta_km)
{
	m_p_mk = m_branch.m_g * std::pow (vM.Get (), 2) -
						1/m_branch.m_tap * vK.Get () * vM.Get () *
						( m_branch.m_g * cos (theta_km) - m_branch.m_b * sin (theta_km) );

	return m_p_mk;
}

double
Branch::CalcQkm(DoubleValue vK, DoubleValue vM,
								DoubleValue aK, DoubleValue aM, double theta_km)
{
	m_q_km = -( m_branch.m_b * ( 1 / std::pow (m_branch.m_tap, 2) ) + m_branch.bsh ) *
							std::pow (vK.Get (), 2) +
							1/m_branch.m_tap * vK.Get () * vM.Get () *
							( m_branch.m_b * cos (theta_km) - m_branch.m_g * sin (theta_km) );

	return m_q_km;
}

double
Branch::CalcQmk(DoubleValue vK, DoubleValue vM,
								DoubleValue aK, DoubleValue aM, double theta_km)
{
	m_q_mk = -( m_branch.m_b + m_branch.bsh ) * std::pow (vM.Get (), 2) +
								1/m_branch.m_tap * vK.Get () * vM.Get () *
								( m_branch.m_b * cos (theta_km) + m_branch.m_g * sin (theta_km) );

	return m_q_mk;
}*/

void
Branch::SetBranch(DBranch_t branch)
{
	m_branch = branch;
}

DBranch_t
Branch::GetBranch() const
{
	return m_branch;
}

double
Branch::CalcPkmL(DoubleValue vK, DoubleValue vM,
		DoubleValue aK, DoubleValue aM)
{
	// s.branch.g(i)*s.bus.v(k)^2 - s.bus.v(k)*s.bus.v(m)*(s.branch.g(i)*cos(akm)+s.branch.b(i)*sin(akm));
	double theta_km = aK.Get () - aM.Get ();
	m_p_km_L = GetBranch ().m_g * pow(vK.Get (), 2) -
			   vK.Get () * vM.Get () *
			   ( GetBranch ().m_g * cos (theta_km) + GetBranch ().m_b * sin (theta_km ) );

	return m_p_km_L;
}

double
Branch::CalcPmkL(DoubleValue vK, DoubleValue vM,
		DoubleValue aK, DoubleValue aM)
{
	// s.branch.g(i)*s.bus.v(m)^2 - s.bus.v(k)*s.bus.v(m)*(s.branch.g(i)*cos(akm)-s.branch.b(i)*sin(akm));
	double theta_km = aK.Get () - aM.Get ();
	m_p_mk_L = GetBranch ().m_g * pow(vM.Get (), 2) -
			   vK.Get () * vM.Get () *
			   (GetBranch ().m_g * cos(theta_km) - GetBranch ().m_b * sin(theta_km) );

	return m_p_mk_L;
}

double
Branch::CalcQkmL(DoubleValue vK, DoubleValue vM,
		DoubleValue aK, DoubleValue aM)
{
	// -(s.branch.b(i)+s.branch.bsh(i))*s.bus.v(k)^2 + s.bus.v(k)*s.bus.v(m)*(s.branch.b(i)*cos(akm)-s.branch.g(i)*sin(akm))
	double theta_km = aK.Get () - aM.Get ();
	m_q_km_L = -(GetBranch ().m_b + GetBranch ().m_bsh) * pow(vK.Get (), 2) +
				vK.Get () * vM.Get () *
				( GetBranch ().m_b * cos(theta_km) - GetBranch ().m_g * sin (theta_km) );

	return m_q_km_L;
}

double
Branch::CalcQmkL(DoubleValue vK, DoubleValue vM,
		DoubleValue aK, DoubleValue aM)
{
	// -(s.branch.b(i)+s.branch.bsh(i))*s.bus.v(m)^2 + s.bus.v(k)*s.bus.v(m)*(s.branch.b(i)*cos(akm)+s.branch.g(i)*sin(akm));
	double theta_km = aK.Get () - aM.Get ();
	m_q_mk_L = -(GetBranch ().m_b + GetBranch ().m_bsh) *
				pow(vM.Get (), 2) + vK.Get () * vM.Get () *
				( GetBranch ().m_b * cos (theta_km) + GetBranch ().m_g * sin (theta_km) );

	return m_q_mk_L;
}

double
Branch::CalcL (DoubleValue vK, DoubleValue vM,
		DoubleValue aK, DoubleValue aM)
{
	// s.branch.g(i)*(s.bus.v(k)^2 + s.bus.v(m)^2 - 2*s.bus.v(k)*s.bus.v(m)*cos(akm));
	double theta_km = aK.Get () - aM.Get ();
	m_l = GetBranch ().m_g *
		  ( pow(vK.Get (), 2) + pow(vM.Get (), 2) - 2 * vK.Get () * vM.Get () * cos (theta_km) );

	return m_l;
}

void
Branch::Print(void)
{
	std::cout << "Branch ( " << m_branch.m_ni << ", " << m_branch.m_nf << ") => "
			 	 	 	<< m_branch.m_b << "\t" << m_branch.m_bsh << "\t" << m_branch.m_g
						<< "\t" << m_branch.m_tap << "\t" << m_branch.m_def << std::endl;
}

}
