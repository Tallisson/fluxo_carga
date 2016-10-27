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

#include "bus.h"

#include "ns3/log.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"
#include "ns3/boolean.h"
#include "ns3/io-utils.h"

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <sstream>

namespace ns3
{

NS_LOG_COMPONENT_DEFINE ("Bus");
NS_OBJECT_ENSURE_REGISTERED (Bus);

Bus::Bus() :
	m_tap (Bus::IMP), m_crt (0.0)
{
	m_dsv = 0;
}

Bus::~Bus()
{
}

TypeId
Bus::GetTypeId (void)
{
	static TypeId tid = TypeId ("ns3::Bus")
	      .SetParent<Object> ()
	      .AddConstructor<Bus> ()
				.AddAttribute("ACalc",
            "Calculated Angle",
            DoubleValue (0), /* TODO */
            MakeDoubleAccessor (&Bus::m_aCalc),
            MakeDoubleChecker<double> ())
				.AddAttribute("VCalc",
				    "Calculated Voltage",
						DoubleValue (0), /* TODO */
						MakeDoubleAccessor (&Bus::m_vCalc),
						MakeDoubleChecker<double> ())
				.AddAttribute("QCalc",
						"Calculated Q",
						DoubleValue (0), /* TODO */
						MakeDoubleAccessor (&Bus::m_qgCalc),
						MakeDoubleChecker<double> ())
				.AddAttribute("PCalc",
						"Calculated P",
						DoubleValue (0), /* TODO */
						MakeDoubleAccessor (&Bus::m_pgCalc),
						MakeDoubleChecker<double> ())
				;

	return tid;
}

void
Bus::SetBus(DBus_t bus)
{
	m_bus = bus;
}

DBus_t
Bus::GetBus() const
{
	return m_bus;
}

void
Bus::SetType(Bus::Type type)
{
	m_type = type;
}

Bus::Type
Bus::GetType(void)
{
	return m_type;
}

void
Bus::AddBranch(Ptr<Branch> branch, Ptr<Bus> neighbor)
{
	m_branches.push_back (branch);
	m_neighbors.push_back (neighbor);
}

std::vector<Ptr<Branch> >
Bus::GetBranches () const
{
	return m_branches;
}

double
Bus::CalcPg (void)
{
	double pgK = 0;
	for (uint32_t i = 0; i < m_branches.size (); i++)
		{
			Ptr<Branch> branch = m_branches.at (i);
			Ptr<Bus> busM = m_neighbors.at (i);
			DBranch_t dataBranch = branch->GetBranch ();

			DoubleValue aM, vM;
			busM->GetAttribute ("ACalc", aM);
			busM->GetAttribute ("VCalc", vM);
	    double theta_km = m_aCalc - aM.Get ();
	    double vK = m_vCalc;

	    if (m_bus.m_nin == dataBranch.m_nf)
	    	{
	        pgK += ( dataBranch.m_g * (1 / pow(dataBranch.m_tap, 2)) *
	        			 pow(vK, 2) - (1 / dataBranch.m_tap) * vK * vM.Get () *
								 ( dataBranch.m_g * cos (theta_km - dataBranch.m_def) +
								 dataBranch.m_b * sin (theta_km - dataBranch.m_def) ) );
	    	}
	    else
	    	{
	        pgK += ( dataBranch.m_g * pow(vK, 2) - (1 / dataBranch.m_tap) *
	        			 vK * vM.Get () * ( dataBranch.m_g * cos(theta_km + dataBranch.m_def) +
	        			 dataBranch.m_b * sin(theta_km + dataBranch.m_def) ) );
	    	}
		}
	pgK += - GetBus ().m_gsh * pow(m_vCalc, 2) + GetBus ().m_pc;
	m_pgCalc = pgK;

	return m_pgCalc;
}

double
Bus::CalcQg (void)
{
	double qgK = 0;
	for (uint32_t i = 0; i < m_branches.size (); i++)
		{
			Ptr<Branch> branch = m_branches.at (i);
			Ptr<Bus> busM = m_neighbors.at (i);
			DBranch_t dataBranch = branch->GetBranch ();

			DoubleValue aM, vM;
			busM->GetAttribute ("ACalc", aM);
			busM->GetAttribute ("VCalc", vM);
	    double theta_km = m_aCalc - aM.Get ();
	    double vK = m_vCalc;

	    if (m_bus.m_nin == dataBranch.m_ni)
	    	{
	        qgK += ( -(dataBranch.m_b * (1 / pow(dataBranch.m_tap, 2)) + dataBranch.m_bsh) *
	        				pow(vK, 2) + (1 / dataBranch.m_tap) * vK * vM.Get () *
									( dataBranch.m_b * cos (theta_km - dataBranch.m_def) -
									dataBranch.m_g * sin (theta_km - dataBranch.m_def) ) );
	    	}
	    else
	    	{
	        qgK += ( -dataBranch.m_b * pow(vK, 2) + ( 1 / dataBranch.m_tap) * vK *vM.Get () *
	        			 ( dataBranch.m_b * cos (theta_km + dataBranch.m_def) - dataBranch.m_g *
	        				sin (theta_km + dataBranch.m_def) ) );
	    	}
		}
	qgK += - GetBus ().m_bsh * pow(m_vCalc, 2) + GetBus ().m_qc;
	m_qgCalc = qgK;

	return m_qgCalc;
}

std::vector<Ptr<Bus> >
Bus::GetNeighbors (void) const
{
	return m_neighbors;
}

void
Bus::SetTap(Bus::TapType tap)
{
	m_tap = tap;
}

Bus::TapType
Bus::GetTap(void) const
{
	return m_tap;
}

void
Bus::Print(void)
{
	std::cout << IOUtils::Format(m_vCalc) << IOUtils::Format(m_aCalc) << std::endl;
}

void
Bus::SetCrt(double crt)
{
	m_crt = crt;
}

double
Bus::GetCrt(void) const
{
	return m_crt;
}

double
Bus::CalcDsv (void)
{
	DoubleValue v;
	GetAttribute("VCalc", v);
	m_dsv = 0;
	if(v.Get () < Bus::MIN_VOLTAGE_ONS)
		{
			m_dsv = Bus::MIN_VOLTAGE_ONS - v.Get ();
			m_status = MIN_VOLTAGE_VIOLATION;
		}

	if (v.Get () > Bus::MAX_VOLTAGE_ONS)
		{
			m_dsv = v.Get () - Bus::MAX_VOLTAGE_ONS;
			m_status = MAX_VOLTAGE_VIOLATION;
		}

	return m_dsv;
}

double
Bus::GetDsv(void)
{
	return m_dsv;
}

Bus::Violation
Bus::GetStatus (void)
{
	return m_status;
}

}
