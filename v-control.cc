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

TypeId
VControl::GetTypeId(void)
{
	static TypeId tid = TypeId ("ns3::VControl")
	    .SetParent<Object> ()
	;

	return tid;
}

}
