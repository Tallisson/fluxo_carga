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

#include "q-control.h"
#include "graph.h"

#include <ns3/double.h>
#include <ns3/log.h>
#include <ns3/uinteger.h>
#include <ns3/assert.h>

#include <map>
#include <math.h>
#include <vector>

namespace ns3
{

QControl::QControl ()
{
}

QControl::~QControl ()
{
}

TypeId
QControl::GetTypeId(void)
{
	static TypeId tid = TypeId ("ns3::QControl")
		    .SetParent<Object> ()
		    .AddConstructor<QControl> ()
	;

	return tid;
}

bool
QControl::DoControl (Ptr<Graph> graph)
{
	// Cálculo da geração de potência reativa e teste de violação dos limites:
	bool update = false;

	for (uint32_t i = 0; i < graph->GetNumBus (); i++)
		{
			Ptr<Bus> busK = graph->GetBus (i + 1);

	    if (busK->GetType () == Bus::SLACK ||
	    			busK->GetType () == Bus::GENERATION)
	    	{
	    		double qg = 0;
	    		std::vector<Ptr<Branch> > branches = busK->GetBranches ();
	    		std::vector<Ptr<Bus> > neighbors = busK->GetNeighbors ();

	        for (uint32_t i = 0; i < branches.size (); i++)
	        	{
	        		DBranch_t dataBranch = branches.at (i)->GetBranch ();
	        		Ptr<Bus> busM = neighbors.at (i);

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
	            		 * s.bus.qg(k) =
	            		 * (-(s.branch.b(km)*(1/(s.branch.tap(km)^2))+s.branch.bsh(km))*
	            		 * s.bus.v(k)^2 + (1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m) *
	            		 * (s.branch.b(km)*cos(akm - s.branch.def(km))-s.branch.g(km)*sin(akm - s.branch.def(km)))) +
	            		 * s.bus.qg(k);
	            		 */
	            		double a = -(dataBranch.m_b * (1 / pow(dataBranch.m_tap, 2)) + dataBranch.m_bsh);
	            		double b = pow(vK.Get (), 2);
	            		double c = (1 / dataBranch.m_tap) * vK.Get () * vM.Get ();
	            		double d = dataBranch.m_b * cos(theta_km - dataBranch.m_def) - dataBranch.m_g *
											 	 	 	 sin(theta_km - dataBranch.m_def);
	                qg = (a * b + c * d) + qg;
	            	}
	            else
	            	{
	            		/*
	            		 * s.bus.qg(k) = (-s.branch.b(km)*s.bus.v(k)^2 + (1/s.branch.tap(km)) *
	            		 * s.bus.v(k)*s.bus.v(m)*(s.branch.b(km)*cos(akm + s.branch.def(km)) -
	            		 * s.branch.g(km)*sin(akm + s.branch.def(km)))) + s.bus.qg(k);
	            		 */
	                qg = (-(dataBranch.m_b) * pow(vK.Get (), 2) +
	                		 (1 / dataBranch.m_tap) * vK.Get () * vM.Get () *
											 (dataBranch.m_b * cos (theta_km + dataBranch.m_def) - dataBranch.m_g *
											 sin(theta_km + dataBranch.m_def)) ) + qg;
	            	}
	        	}

	        DoubleValue vK, aK, qgK;
	        busK->GetAttribute ("VCalc", vK);
	        busK->GetAttribute ("ACalc", aK);
	        // s.bus.qg(k) = - s.bus.bsh(k)*s.bus.v(k)^2 + s.bus.qc(k) + s.bus.qg(k);
	        qg = - busK->GetBus ().m_bsh * pow(vK.Get (), 2) + busK->GetBus ().m_qc + qg;
	        /*std::cout << "Bus" << busK->GetBus ().m_nin << ", Qg = " << qg << ", qgMin = " << busK->GetBus ().m_qgmin <<
	        						 ", qgMax = " << busK->GetBus ().m_qgmax << ", ThetaK = " << aK.Get () << std::endl;*/
	        if (busK->GetBus ().m_nin == 2)
						{
							std::cout << "Qg = " << qg << std::endl;
						}
	        if (qg < busK->GetBus ().m_qgmin)
	        	{
	            qg = busK->GetBus ().m_qgmin;
	            busK->SetType (Bus::LOSS_CONTROL_REACT);
	            update = true;
	        	}
	        else if (qg > busK->GetBus ().m_qgmax)
						{
	            qg = busK->GetBus ().m_qgmax;
	            busK->SetType (Bus::LOSS_CONTROL_REACT);
	            update = true;
						}

	        qgK.Set(qg);
	        busK->SetAttribute ("QCalc", qgK);
	    	}
		}
	if (update)
		{
			UpdateOrd (graph);
		}

	return update;
}

bool
QControl::DoRestore (Ptr<Graph> graph)
{
	// Verificando quais barras recuperaram controle de reativo:
	bool update = false;
	for (uint32_t i = 0; i < graph->GetNumBus (); i++)
		{
			Ptr<Bus> busK = graph->GetBus (i + 1);

			if (busK->GetType () == Bus::LOSS_CONTROL_REACT)
				{
					DoubleValue vK, qgK;
					busK->GetAttribute("VCalc", vK);
					busK->GetAttribute("QCalc", qgK);
					if (qgK.Get () == busK->GetBus ().m_qgmin &&
								vK.Get () <= busK->GetBus ().m_v)
						{
							std::cout << "Min = " << vK.Get () << ", " << busK->GetBus ().m_v << std::endl;
							busK->SetAttribute("QCalc", DoubleValue (0));
							busK->SetType (Bus::GENERATION);
							update = true;
						}

					if (qgK.Get () == busK->GetBus ().m_qgmax &&
								vK.Get () >= busK->GetBus ().m_v)
						{
							std::cout << "Max = " << vK.Get () << ", " << busK->GetBus ().m_v << std::endl;
							busK->SetAttribute("QCalc", DoubleValue (0));
							busK->SetType (Bus::GENERATION);
							update = true;
						}
				}
		}

	if (update)
		{
			UpdateOrd (graph);
		}

	return update;
}

void
QControl::UpdateOrd(Ptr<Graph> graph)
{
	uint32_t cont = 0;
	uint32_t cont1 = 0;
	uint32_t cont2 = 0;

	uint32_t npv = 0;
	uint32_t npq = 0;

	for (uint32_t i = 0; i < graph->GetNumBus (); i++)
		{
			Ptr<Bus> bus = graph->GetBus (i+1);
			DBus_t oldBus = bus->GetBus ();
			DBus_t newBus = oldBus;
	    if (bus->GetType () == Bus::GENERATION)
	    	{
	    		newBus.m_ord = cont++;
	    		newBus.m_ordPV = cont1++;
	    		newBus.m_posPQ = i;
	    		bus->SetBus (newBus);
	        npv++;
	    	}

	    if (bus->GetType () == 1 || bus->GetType () == Bus::LOAD ||
	    					bus->GetType () == Bus::LOSS_CONTROL_REACT)
	    	{
	    		newBus.m_ord = cont++;
	        newBus.m_ordPQ = cont2++;
	        newBus.m_posPQ = i;
	        bus->SetBus (newBus);
	        npq++;
	    	}

	    graph->SetNumPQ (npq);
	    graph->SetNumPV (npv);
		}
}

}
