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

#include "mismatches.h"

#include "graph.h"

#include <ns3/log.h>

#include <math.h>
#include <iostream>

namespace ns3
{

NS_LOG_COMPONENT_DEFINE ("Mismatch");
NS_OBJECT_ENSURE_REGISTERED (Mismatch);

Mismatch::Mismatch ()
{
}

Mismatch::~Mismatch ()
{
}

TypeId
Mismatch::GetTypeId ()
{
	static TypeId tid = TypeId ("ns3::Mismatch")
					.SetParent<Object> ()
					.AddConstructor<Mismatch> ()
	;

	return tid;
}

void
Mismatch::CalcPkB (Ptr<Graph> graph)
{
	// Balanço de potência ativa:
	for (uint32_t i = 0; i < graph->GetNumBus (); i++)
		{
			Ptr<Bus> busK = graph->GetBus (i+1);
			if (busK->GetType () == Bus::SLACK)
				{
					continue;
				}
			uint32_t k = busK->GetBus ().m_ord;
			m_mis (k) = 0;

			std::vector<Ptr<Branch> > branches = busK->GetBranches();
			std::vector<Ptr<Bus> > neighbors = busK->GetNeighbors ();
			for (uint32_t j = 0; j < branches.size (); j++)
				{
					DBranch_t dataBranch = branches.at (j)->GetBranch ();
					Ptr<Bus> busM = neighbors.at (j);

					DoubleValue vK, vM, aK, aM;
					busK->GetAttribute ("VCalc", vK);
					busM->GetAttribute ("VCalc", vM);
					busK->GetAttribute ("ACalc", aK);
					busM->GetAttribute ("ACalc", aM);

	        double theta_km = aK.Get () - aM.Get ();

	        //if (busK->GetBus ().m_nin == dataBranch.m_ni)
	        if (dataBranch.m_tipo == 1 && busK->GetTap () == Bus::TAP)
	        	{
            	/*mis(k-1) =
            			(s.branch.g(km)*(1/s.branch.tap(km)^2)*s.bus.v(k)^2 -
            			(1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
            			(s.branch.g(km)*cos(akm - s.branch.def(km))+s.branch.b(km)*
            			sin(akm - s.branch.def(km)))) + mis(k-1);
            	*/
	            m_mis(k) = (dataBranch.m_g * 1 / pow (dataBranch.m_tap, 2) * pow(vK.Get (), 2) -
	            					 (1 / dataBranch.m_tap) * vK.Get () * vM.Get () *
												 ( dataBranch.m_g * cos(theta_km - dataBranch.m_def) +
													dataBranch.m_b * sin(theta_km - dataBranch.m_def) ) ) + m_mis(k);
	        	}
	        else
	        	{
	        		/*
	        		 * mis(k-1) = (s.branch.g(km)*s.bus.v(k)^2 - (1/s.branch.tap(km))*
	        		 * s.bus.v(k)*s.bus.v(m)*(s.branch.g(km)*cos(akm + s.branch.def(km))+
	        		 * s.branch.b(km)*sin(akm + s.branch.def(km)))) + mis(k-1);
	        		 */
	            m_mis(k) = (dataBranch.m_g * pow(vK.Get (), 2) -
											 	 (1 / dataBranch.m_tap) * vK.Get () * vM.Get () *
												 ( dataBranch.m_g * cos(theta_km + dataBranch.m_def) +
													dataBranch.m_b * sin(theta_km + dataBranch.m_def) ) ) + m_mis(k);
	        	}
				}

				DoubleValue vK, pgK;
				busK->GetAttribute ("VCalc", vK);
				busK->GetAttribute("PCalc", pgK);
				/*
				 * s.bus.pg(k) + s.bus.gsh(k)*s.bus.v(k)^2 - s.bus.pc(k) - mis(k-1);
				 */
				m_mis(k) = pgK.Get () + busK->GetBus ().m_gsh *
									 pow(vK.Get (), 2) - busK->GetBus ().m_pc - m_mis(k);
		}
}

void
Mismatch::CalcQkB (Ptr<Graph> graph)
{
	// Balanço de potência reativa:
	for (uint32_t i = 0; i < graph->GetNumBus (); i++)
		{
			Ptr<Bus> busK = graph->GetBus (i + 1);
			if (busK->GetType() == Bus::LOSS_CONTROL_REACT ||
						busK->GetType () == Bus::LOAD)
				{
					uint32_t index = graph->GetNumBus () - 1 + busK->GetBus ().m_ordPQ;
					m_mis (index) = 0;

					std::vector<Ptr<Branch> > branches = busK->GetBranches();
					std::vector<Ptr<Bus> > neighbors = busK->GetNeighbors ();
					for (uint32_t j = 0; j < branches.size (); j++)
						{
							DBranch_t dataBranch = branches.at (j)->GetBranch ();
							Ptr<Bus> busM = neighbors.at (j);

							DoubleValue vK, vM, aK, aM;
							busK->GetAttribute ("VCalc", vK);
							busM->GetAttribute ("VCalc", vM);
							busK->GetAttribute ("ACalc", aK);
							busM->GetAttribute ("ACalc", aM);
							double theta_km = aK.Get () - aM.Get ();

							//if (busK->GetBus ().m_nin == dataBranch.m_ni)
							if (dataBranch.m_tipo == 1 && busK->GetTap () == Bus::TAP)
								{
									/*
									 * mis(s.nb-1+s.bus.ordPQ(k)) =
									 * (-(s.branch.b(km)*(1/(s.branch.tap(km)^2))+s.branch.bsh(km))*
									 * s.bus.v(k)^2 + (1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
									 * (s.branch.b(km)*cos(akm - s.branch.def(km))-s.branch.g(km)*
									 * sin(akm - s.branch.def(km)))) + mis(s.nb-1+s.bus.ordPQ(k));
									 */
									m_mis(index) = (-(dataBranch.m_b * (1 / pow(dataBranch.m_tap, 2)) + dataBranch.m_bsh ) *
																	pow (vK.Get (), 2) + (1/dataBranch.m_tap) * vK.Get () * vM.Get () *
																	(dataBranch.m_b * cos(theta_km - dataBranch.m_def) - dataBranch.m_g *
																	sin(theta_km - dataBranch.m_def) ) ) + m_mis(index);
								}
							else
								{
									/*
									 * mis(s.nb-1+s.bus.ordPQ(k)) =
									 * (-s.branch.b(km)*s.bus.v(k)^2 + (1/s.branch.tap(km))*s.bus.v(k)*
									 * s.bus.v(m)*(s.branch.b(km)*cos(akm + s.branch.def(km))-s.branch.g(km)*
									 * sin(akm + s.branch.def(km)))) + mis(s.nb-1+s.bus.ordPQ(k));
									 */
									m_mis(index) = (-(dataBranch.m_b) * pow (vK.Get (), 2) +
																	(1/dataBranch.m_tap) * vK.Get () * vM.Get () *
																	(dataBranch.m_b * cos(theta_km + dataBranch.m_def) - dataBranch.m_g *
																	sin(theta_km + dataBranch.m_def) ) ) + m_mis(index);
								}
						}

					DoubleValue vK, qgK;
					busK->GetAttribute ("VCalc", vK);
					busK->GetAttribute("QCalc", qgK);
					/*
					 * mis(s.nb-1+s.bus.ordPQ(k)) =
					 * s.bus.qg(k) + s.bus.bsh(k)*s.bus.v(k)^2 - s.bus.qc(k) - mis(s.nb-1+s.bus.ordPQ(k));
					*/
					m_mis(index) = qgK.Get () + busK->GetBus ().m_bsh * pow (vK.Get (), 2) -
												 busK->GetBus ().m_qc - m_mis(index);
				}
		}
}

vec
Mismatch::CalcMismatches(Ptr<Graph> graph)
{
	m_mis = zeros<vec>(graph->GetNumPQ () * 2 + graph->GetNumPV ());
	CalcPkB (graph);
	CalcQkB (graph);

	return GetMis ();
}

vec
Mismatch::GetMis(void)
{
	return m_mis;
}

}
