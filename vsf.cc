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

Ptr<Bus>
Vsf::MaxDsv (Ptr<Graph> graph)
{
  Ptr<Bus> max = NULL;
  for (uint32_t i = 1; i <= graph->GetNumBus (); i++)
    {
      Ptr<Bus> bus = graph->GetBus (i);
      if (bus->GetType () != Bus::LOAD)
        {
          continue;
        }
      if (max == NULL)
        {
          if (bus->GetDsv () != 0)
            {
              max = bus;
            }
        }
      else if (fabs (bus->GetDsv ()) > fabs (max->GetDsv ()))
        {
          max = bus;
        }
    }

  return max;
}

bool
Vsf::DoControl (Ptr<Graph> graph)
{
	bool control = false;

  std::vector<double> volts;
  std::vector<uint32_t> index;
  for (uint32_t i = 1; i <= graph->GetNumBus (); i++)
    {
      Ptr<Bus> bus = graph->GetBus (i);
      if (bus->GetType () != Bus::LOAD)
        {
          continue;
        }
      double dsv = bus->CalcDsv();
      if (dsv != 0)
        {
      		control = true;
          volts.push_back (dsv);
          index.push_back (bus->GetBus ().m_nin);
        }
    }

	if (volts.size () > 0)
		{
			Ptr<Bus> maxVlt = MaxDsv (graph);

			mat jqv = m_jac->GetJqv ();
      mat invJqv = inv (jqv);

			vec deltaV = zeros<vec> (jqv.n_cols);
			for (uint32_t i = 0; i < volts.size (); i++)
				{
					uint32_t ind = index.at (i) - 1;
					deltaV (ind) = volts.at (i);
				}

			vec deltaQ = (jqv * deltaV);

			NS_LOG_INFO ("Delta V: " << endl << deltaV);
			NS_LOG_INFO ("Delta Q: " << endl << deltaQ);
			double maxQ = 0;
			uint32_t maxElem = 0;
			for (uint32_t i = 1; i <= deltaQ.n_elem; i++)
				{
					Ptr<Bus> crtBus = graph->GetBus (i);

					if (crtBus->GetType () == Bus::LOAD)
						{
							continue;
						}

					if (crtBus->GetType () != Bus::GENERATION)
						{
							continue;
						}

					DoubleValue v;
					crtBus->GetAttribute ("VCalc", v);
					if (v.Get () == Bus::MAX_VOLTAGE_ONS && maxVlt->GetStatus () == Bus::MIN_VOLTAGE_VIOLATION)
						{
							continue;
						}

					if (v.Get () == Bus::MIN_VOLTAGE_GR && maxVlt->GetStatus () == Bus::MAX_VOLTAGE_VIOLATION)
						{
							continue;
						}

          uint32_t idX = i - 1;
          uint32_t idY = maxVlt->GetBus ().m_nin - 1;

          double value = invJqv (idX, idY);
          crtBus->SetCrt (value);
					std::cout << "Bus " << i << "=> value = " << value << "\n";
					if (value == 0)
						{
							continue;
						}
					if ( maxElem == 0 || maxQ < value )
						{
							maxQ = value;
							maxElem = i;
						}
				}
			if (maxElem == 0)
				{
					return false;
				}

			std::cout << "Max elemento: " << (maxElem-1) << ", max value = " << maxQ << std::endl;
			vec auxQ = zeros<vec> (deltaQ.n_elem);
			auxQ(maxElem -1) = maxQ;
			vec deltaVIjt = inv (jqv) * auxQ;
			std::cout << "AuxQ: " << endl << auxQ << std::endl;
			std::cout << "Delta V injected: " << std::endl << deltaVIjt << std::endl;

			double m_alpha = 1;

			double value = fabs (deltaVIjt (maxElem));
			std::cout << "Variação de Tensão => " << deltaVIjt << "\n";
			while (m_alpha * value < LIMIAR)
				{
					m_alpha++;
				}

			Ptr<Bus> bus = graph->GetBus (maxElem);
			std::cout << "Regulando Tensões " << std::endl;
			value = deltaVIjt (maxElem - 1);
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
			std::cout << "Incremento na barra " << (maxElem) << " = " << (value * m_alpha) << std::endl;
			std::cout << "Variação de potência reativa na barra " << (maxElem) << " => " << maxQ << std::endl;
			std::cout << "ALPHA = " << m_alpha << std::endl;

			std::cout << "Voltage + value = " << (newValue) << std::endl;

			if ( newValue > Bus::MAX_VOLTAGE_ONS )
				{
					newValue = Bus::MAX_VOLTAGE_ONS;
				}
			if ( newValue < Bus::MIN_VOLTAGE_ONS )
				{
					newValue = Bus::MIN_VOLTAGE_ONS;
				}
			bus->SetAttribute ("VCalc", DoubleValue (newValue));
		}

	return control;
}

bool
Vsf::DoRestore (Ptr<Graph> graph)
{
	return false;
}

void
Vsf::SetJac (Ptr<Jacobian> jac)
{
	m_jac = jac;
}

}
