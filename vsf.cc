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

#include "vsf.h"
#include "graph.h"

#include <ns3/double.h>
#include <ns3/log.h>
#include <ns3/uinteger.h>
#include <ns3/assert.h>

#include <armadillo>
#include <map>
#include <math.h>
#include <vector>

using namespace arma;

namespace ns3
{

NS_LOG_COMPONENT_DEFINE ("Vsf");
NS_OBJECT_ENSURE_REGISTERED (Vsf);

Vsf::Vsf ()
{
}

Vsf::~Vsf ()
{
}

TypeId
Vsf::GetTypeId(void)
{
	static TypeId tid = TypeId ("ns3::Vsf")
		    .SetParent<Object> ()
		    .AddConstructor<Vsf> ()
	;

	return tid;
}

bool
Vsf::DoControl (mat jqv, Ptr<Graph> graph)
{
	bool control = false;
	vec deltaV = zeros<vec> (jqv.n_cols);
	for (uint32_t i = 1; i <= graph->GetNumBus (); i++)
		{
			Ptr<Bus> bus = graph->GetBus (i);
			if (bus->GetType () != Bus::LOAD &&
						bus->GetType () != Bus::LOSS_CONTROL_REACT)
				{
					continue;
				}
			double dsv = bus->CalcDsv();
			if (dsv != 0)
				{
					deltaV (i - 1) = dsv;
					control = true;
				}
		}

	if (control == true)
		{
			Ptr<Bus> maxVlt = MaxDsv (graph);
			mat invJqv = inv (jqv);

			subview_col<double> s = invJqv.col(maxVlt->GetBus ().m_nin - 1);
			vec vsf = zeros<vec> (s.n_elem);
			for(uint32_t i = 0; i < s.n_elem; i++)
				{
					vsf (i) = s (i);
				}

			uint32_t idBus = MaxV (graph, vsf, maxVlt);
			DoubleValue t;
			maxVlt->GetAttribute ("VCalc", t);
			std::cout << "Max Vlt = " << maxVlt->GetBus ().m_nin << " = " << t.Get () << ", Max = " << vsf (idBus - 1) << std::endl;
			vec aux = zeros<vec> (vsf.n_elem);
			aux (idBus - 1) = vsf (idBus - 1);
			vec deltaVIjt = inv (jqv) * aux;

			double m_alpha = 1;
			double value = fabs (deltaVIjt (idBus - 1));
			while (m_alpha * value < LIMIAR)
				{
					m_alpha++;
				}
			std::cout << "Variação de Tensão => " << value << ", alpha = " << m_alpha << "\n";

			Ptr<Bus> bus = graph->GetBus (idBus);
			value = deltaVIjt (idBus - 1);
			if (maxVlt->GetStatus () == Bus::MIN_VOLTAGE_VIOLATION && value < 0)
				{
					value = fabs (value);
				}

			if (maxVlt->GetStatus () == Bus::MAX_VOLTAGE_VIOLATION && value > 0)
				{
					value *= -1;
				}
			DoubleValue v;
			bus->GetAttribute("VCalc", v);
			double newValue = v.Get () + (value * m_alpha);
			if (newValue < Bus::MIN_VOLTAGE_GR)
				{
					newValue = Bus::MIN_VOLTAGE_GR;
				}
			if (newValue > Bus::MAX_VOLTAGE_ONS)
				{
					newValue = Bus::MAX_VOLTAGE_ONS;
				}

			std::cout << "Value = " << v.Get () << " Bus = " << bus->GetType () << ", Voltage + value = " << (newValue) << std::endl;
			bus->SetAttribute ("VCalc", DoubleValue (newValue));
		}

	return control;
}

}
