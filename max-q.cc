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

#include "max-q.h"
#include "graph.h"

#include <ns3/double.h>
#include <ns3/log.h>
#include <ns3/uinteger.h>
#include <ns3/assert.h>

#include <armadillo>
#include <iostream>
#include <map>
#include <math.h>
#include <vector>

using namespace arma;

namespace ns3
{

NS_LOG_COMPONENT_DEFINE ("MaxQ");
NS_OBJECT_ENSURE_REGISTERED (MaxQ);

MaxQ::MaxQ ()
{
}

MaxQ::~MaxQ ()
{
}

TypeId
MaxQ::GetTypeId(void)
{
	static TypeId tid = TypeId ("ns3::MaxQ")
		    .SetParent<Object> ()
		    .AddConstructor<MaxQ> ()
	;

	return tid;
}

Ptr<Bus>
MaxQ::MaxDsv (Ptr<Graph> graph)
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
MaxQ::DoControl (Ptr<Graph> graph)
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
          index.push_back (i - 1);
        }
    }

  if (volts.size () > 0)
    {
  		Ptr<Bus> maxVlt = MaxDsv (graph);

      mat jqv = m_jac->GetJqv ();
      vec deltaV = zeros<vec> (jqv.n_cols);
      for (uint32_t i = 0; i < volts.size (); i++)
        {
          uint32_t ind = index.at (i);
          deltaV(ind) = volts.at (i);
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

          DoubleValue v;
          crtBus->GetAttribute ("VCalc", v);
          if (v == Bus::MAX_VOLTAGE_ONS && maxVlt->GetStatus () == Bus::MIN_VOLTAGE_VIOLATION)
            {
              continue;
            }

          if (v == Bus::MIN_VOLTAGE_GR && maxVlt->GetStatus () == Bus::MAX_VOLTAGE_VIOLATION)
            {
              continue;
            }

          double value = deltaQ (i-1);
          crtBus->SetCrt (deltaQ (i-1));
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
			std::cout << "Variação de Tensão => " << deltaVIjt () << "\n";
			if (value == 0)
				{
					continue;
				}
			while (m_alpha * value < LIMIAR)
				{
					m_alpha++;
				}

      Ptr<Bus> bus = graph->GetBus (maxElem);
      std::cout << "Regulando Tensões " << std::endl;
      double value = deltaVIjt (maxElem - 1);
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
MaxQ::DoRestore (Ptr<Graph> graph)
{
	return false;
}

}