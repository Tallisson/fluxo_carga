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

#include "v-control.h"

#include <ns3/assert.h>
#include <ns3/double.h>
#include <ns3/log.h>

namespace ns3
{

Ptr<Bus>
VControl::MaxDsv (Ptr<Graph> graph)
{
  Ptr<Bus> maxBus = graph->GetBus (1);
  double maxDsv = maxBus->GetDsv ();
  for (uint32_t i = 1; i <= graph->GetNumBus (); i++)
    {
      Ptr<Bus> bus = graph->GetBus (i);
      if (bus->GetType () != Bus::LOAD &&
      			bus->GetType () != Bus::LOSS_CONTROL_REACT)
        {
          continue;
        }

			if (bus->GetDsv () > maxDsv)
				{
					maxBus = bus;
					maxDsv = bus->GetDsv ();
				}

    }

  return maxBus;
}

uint32_t
VControl::MaxV (Ptr<Graph> graph, arma::vec vt, Ptr<Bus> modBus)
{
	double maxQ = vt (0);
	uint32_t maxBus = 1;
	for (uint32_t i = 0; i < vt.n_elem; ++i)
		{
			Ptr<Bus> crtBus = graph->GetBus (i+1);
			if (crtBus->GetType () == Bus::LOAD ||
		          		crtBus->GetType () == Bus::LOSS_CONTROL_REACT)
		  	{
					continue;
		  	}

			  crtBus->SetCrt (vt (i));
				DoubleValue v;
			  crtBus->GetAttribute ("VCalc", v);

				if (v.Get () == Bus::MAX_VOLTAGE_ONS && modBus->GetStatus () == Bus::MIN_VOLTAGE_VIOLATION)
					{
						continue;
					}

				if (v.Get () == Bus::MIN_VOLTAGE_GR && modBus->GetStatus () == Bus::MAX_VOLTAGE_VIOLATION)
					{
						continue;
					}

				if (vt (i) > maxQ)
					{
						maxQ = vt (i);
						maxBus = i + 1;
					}
		}

	return maxBus;
}


TypeId
VControl::GetTypeId(void)
{
	static TypeId tid = TypeId ("ns3::VControl")
	    .SetParent<Object> ()
	;

	return tid;
}

}
