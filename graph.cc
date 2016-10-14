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

#include "graph.h"

#include "ns3/log.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"

namespace ns3
{

NS_LOG_COMPONENT_DEFINE ("Graph");
NS_OBJECT_ENSURE_REGISTERED (Graph);

Graph::Graph () :
	m_numBus(0), m_numBranch(0),
	m_numPQ (0), m_numPV (0),
	m_numSlack (0)
{
}

Graph::~Graph ()
{
	m_buses.clear ();
	m_numBus = 0;
	m_numBranch = 0;
	m_numPQ = 0;
	m_numPV = 0;
	m_numSlack = 0;
}

TypeId
Graph::GetTypeId (void)
{
	static TypeId tid = TypeId ("ns3::Graph")
	      .SetParent<Object> ()
	      .AddConstructor<Graph> ()
				.AddAttribute("NumBus",
						"Number of Buses",
						UintegerValue (0), /* TODO */
						MakeDoubleAccessor (&Graph::m_numBus),
						MakeDoubleChecker<uint32_t> ())
				.AddAttribute("NumBranch",
						"Number of Branches",
						UintegerValue (0), /* TODO */
						MakeDoubleAccessor (&Graph::m_numBranch),
						MakeDoubleChecker<uint32_t> ())
				.AddAttribute("NumPQ",
						"Number of PQ buses",
						UintegerValue (0), /* TODO */
						MakeDoubleAccessor (&Graph::m_numPQ),
						MakeDoubleChecker<uint32_t> ())
				.AddAttribute("NumPV",
						"Number of PQ buses",
						UintegerValue (0), /* TODO */
						MakeDoubleAccessor (&Graph::m_numPV),
						MakeDoubleChecker<uint32_t> ())
				.AddAttribute("NumSlack",
						"Number of Slack buses",
						UintegerValue (1), /* TODO */
						MakeDoubleAccessor (&Graph::m_numSlack),
						MakeDoubleChecker<uint32_t> ())
	;

	return tid;
}

void
Graph::AddBus(Ptr<Bus> bus)
{
	uint32_t id = bus->GetBus ().m_nin;
	m_buses.insert (std::pair<uint32_t, Ptr<Bus> > (id, bus));
	m_numBus++;

	if (bus->GetBus ().m_tipo == Bus::SLACK)
		{
			bus->SetType (Bus::SLACK);
			bus->SetAttribute ("VCalc", DoubleValue (bus->GetBus ().m_v));
			bus->SetAttribute ("ACalc", DoubleValue (bus->GetBus ().m_ang));
			m_numSlack++;
		}

	if (bus->GetBus ().m_tipo == Bus::GENERATION)
		{
			bus->SetType (Bus::GENERATION);
			bus->SetAttribute ("VCalc", DoubleValue (bus->GetBus ().m_v));
			bus->SetAttribute ("PCalc", DoubleValue (bus->GetBus ().m_pg));
			m_numPV++;
		}

	if (bus->GetBus ().m_tipo == Bus::LOAD
			|| bus->GetBus ().m_tipo == Bus::LOSS_CONTROL_REACT)
		{
			if (bus->GetBus ().m_tipo == 0)
				{
					bus->SetType (Bus::LOAD);
				}
			else if(bus->GetBus ().m_tipo == 4)
				{
					bus->SetType (Bus::LOSS_CONTROL_REACT);
				}
			bus->SetAttribute ("PCalc", DoubleValue (bus->GetBus ().m_pg));
			bus->SetAttribute ("QCalc", DoubleValue (bus->GetBus ().m_qg));
			m_numPQ++;
		}

}

void
Graph::Assoc(Ptr<Branch> branch)
{
	uint32_t idK = branch->GetBranch ().m_nf;
	uint32_t idM = branch->GetBranch ().m_ni;

	std::map<uint32_t, Ptr<Bus> >::iterator itK, itM;
	Ptr<Bus> busK = m_buses.find(idK)->second;
	Ptr<Bus> busM = m_buses.find(idM)->second;

	busK->AddBranch (branch, busM);
	busM->AddBranch (branch, busK);

	m_branches.push_back (branch);
	m_numBranch++;
}

std::map<uint32_t, Ptr<Bus> >
Graph::GetBuses() const
{
	return m_buses;
}

uint32_t
Graph::GetNumBus() const
{
	return m_numBus;
}

uint32_t
Graph::GetNumBranch() const
{
	return m_numBranch;
}

Ptr<Bus>
Graph::GetBus(uint32_t idBus)
{
	std::map<uint32_t, Ptr<Bus> >::iterator it = m_buses.find (idBus);
	Ptr<Bus> bus = 0;
	if (it != m_buses.end ())
		{
			bus = it->second;
		}

	return bus;
}

std::vector<Ptr<Branch> >
Graph::GetBranches() const
{
	return m_branches;
}

uint32_t
Graph::GetPosSlack(void) const
{
	return m_posSlack;
}

void
Graph::SetPosSlack (uint32_t posSlack)
{
	m_posSlack = posSlack;
}

uint32_t
Graph::GetNumPQ () const
{
	return m_numPQ;
}

void
Graph::SetNumPQ (uint32_t numPQ)
{
	m_numPQ = numPQ;
}

uint32_t
Graph::GetNumPV () const
{
	return m_numPV;
}

void
Graph::SetNumPV (uint32_t numPV)
{
	m_numPV = numPV;
}

}
