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

#include "io-utils.h"

#include "ns3/log.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>

namespace ns3
{

const std::string SEP = " ";

IOUtils::IOUtils()
{
}

IOUtils::~IOUtils()
{
}

TypeId
IOUtils::GetTypeId()
{
	static TypeId tid = TypeId ("ns3::IOUtils")
			      .SetParent<Object> ()
			      .AddConstructor<IOUtils> ()
	;

	return tid;
}

std::string
IOUtils::Print(std::vector<Data_t> values, uint8_t precision)
{
	std::ostringstream osV, osL;

	osV.unsetf (std::ios::floatfield);
	osV.precision (precision);
	osV.setf ( std::ios::fixed, std:: ios::floatfield ); // floatfield set to fixed

	for (uint32_t i = 0; i < values.size (); i++) {
		std::string label = values.at (i).m_label;
		osL << label;
		uint32_t size_s = label.length ();
		if (size_s < precision)
			{
				uint32_t space = precision - size_s;
				for (uint32_t j = 0; j < space; j++)
					{
						osL << " ";
					}
				if (i < values.size () - 1)
					{
						osL << "|";
					}
			}
	}
	osL << std::endl;

	for (uint32_t i = 0; i < values.size (); i++)
		{
			osV << values.at (i).m_value;
			if (i < values.size () - 1)
				{
					osV << "|";
				}
		}
	std::cout << osV.str ();
	osL << osV.str ();
	return osL.str ();
}

std::string
IOUtils::Format (double v, uint32_t p, uint32_t maxS)
{
  /*std::string sep = " ";
  std::ostringstream os;
  os.unsetf (std::ios::floatfield);
  os.precision (p);
  os.setf ( std::ios::fixed, std:: ios::floatfield ); // floatfield set to fixed
  os << v << sep;

  std::string s;
  for (uint32_t i = 0; i < (maxS - os.str ().size ()); i++)
    {
      s.append (sep);
    }
  std::ostringstream t;
  t << (s + os.str ());

  return ( t.str () );*/

	std::ostringstream os;
	std::ostringstream osE;
	if (p > 0)
		{
			os.unsetf (std::ios::floatfield);
			os.precision (p);
			os.setf ( std::ios::fixed, std:: ios::floatfield ); // floatfield set to fixed
		}
	else
		{
			if (v == 0)
				{
					os.unsetf (std::ios::floatfield);
					os.precision (1);
					os.setf ( std::ios::fixed, std:: ios::floatfield ); // floatfield set to fixed
				}
		}

	os << v;
	uint32_t l = os.str ().size ();
	if (l > maxS)
		{
			uint32_t prec = 1;
			uint32_t s;
			std::string str;
			do {
					std::ostringstream aux;
					aux.unsetf (std::ios::floatfield);
					aux.precision (prec++);
					aux.setf ( std::ios::fixed, std:: ios::floatfield ); // floatfield set to fixed
					aux << v;
					str = aux.str ();
					s = str.size ();
			} while (s < maxS);
			osE << str;
		}
	else
		{
			osE << os.str ();
		}

	std::string s;
	int32_t size = (maxS - os.str ().size ());
	for (int32_t i = 0; i < size; i++)
		{
			s.append (" ");
		}

	std::ostringstream t;
	t << (s + osE.str ());

	return ( t.str () );
}

std::string
IOUtils::Format (std::string v, uint32_t maxS)
{
	/*std::string sep = " ";
	std::ostringstream os;
	os << v;

	std::string s;
  for (uint32_t i = 0; i < (maxS - os.str ().size ()); i++)
    {
      s.append (sep);
    }
  std::ostringstream t;
  t << (s + os.str ());

  return ( t.str () );*/

	std::ostringstream os;
	os << v;

	std::string aux;
	int32_t size = (maxS - os.str ().size ());
	for (int32_t i = 0; i < size; i++)
		{
			aux.append (" ");
		}
	std::ostringstream t;
	t << (aux + os.str ());

	return t.str ();
}


std::string
IOUtils::WriteCdf (Ptr<Graph> graph, double getBase)
{
	std::ostringstream otp;
	otp << "BUS DATA FOLLOWS                            " << graph->GetNumBus () << " ITEMS" << std::endl;

	for (uint32_t i = 1; i <= graph->GetNumBus (); i++)
		{
			Ptr<Bus> bus = graph->GetBus (i);

			DBus_t data = bus->GetBus ();

			std::ostringstream os;
			std::ostringstream n;
			n << data.m_nex;
			std::string nex = n.str ();
			for (uint32_t i = 0; i < (4 - nex.size ()); i++)
				{
					os << " ";
				}
			os << nex;

			DoubleValue vC, aC, pG, qG;
			bus->GetAttribute ("VCalc", vC);
			bus->GetAttribute ("ACalc", aC);
			bus->GetAttribute ("PCalc", pG);
			bus->GetAttribute ("QCalc", qG);

			std::string v = Format (vC.Get (), 3, 6);
			std::string a = Format (aC.Get () * 180 / M_PI, 2, 6);
			std::string pc = Format (data.m_pc * getBase, 0, 8);
			std::string qc = Format (data.m_qc * getBase, 0, 8);
			std::string pg = Format (pG.Get () * getBase, 0, 8);
			std::string qg = Format (qG.Get () * getBase, 0, 7);

			std::string base_kV = Format (data.m_base_kV, 0, 7);
		  std::string vg_o;
		  if (data.m_vg_o != 0.65)
			  {
			  	vg_o = Format (data.m_vg_o, 0, 6);
			  }
			else
				{
					vg_o = Format (0, 0, 6);
				}
		  std::string max_l;
		  std::string min_l;
		  double max = data.m_qgmax * getBase;
		  if (max == 999)
		  	{
		  		max = 0;
		  	}
		  double min = data.m_qgmin * getBase;
		  if (min == -999)
		  	{
		  		min = 0;
		  	}
		  switch (data.m_tipo)
		  {
		  	case 3:
		  	case 2:
		  		max_l = Format (max, 0, 7);
		  		min_l = Format (min, 0, 7);
		  		break;
		  	case 1:
		  		max_l = Format (data.m_vmax, 0, 7);
					min_l = Format (data.m_vmin, 0, 7);
					break;
		  	default:
		  		max_l = Format (0, 0, 7);
		  		min_l = Format (0, 0, 7);
		  		break;
		  }

			std::string gsh = Format (data.m_gsh, 0, 7);
			std::string bsh = Format (data.m_bsh, 0, 7);
			std::string ctrl_rem = Format (data.m_ctrl_rem, 0, 4);

			os << " " << data.m_nome << " " << Format(data.m_area, 0, 1) << " " << Format(data.m_zona, 0, 2)
				 << " " << Format(data.m_tipo, 0, 2) << " " << v << " " << a << " " << pc << " " << qc
				 << " " << pg << " " << qg <<  " " << base_kV << " " << vg_o << " " << max_l << " "
				 << min_l << " " << gsh << " " << bsh << " " << ctrl_rem << " " << Format(nex, 4);

			otp << os.str () << std::endl;
		}

	otp << "-999" << std::endl
		 	 								<< "BRANCH DATA FOLLOWS                         " << graph->GetNumBranch() << " ITEMS" << std::endl;

	std::vector< Ptr<Branch> > branches = graph->GetBranches ();
	for (uint32_t t = 0; t < branches.size (); t++)
		{
			Ptr<Branch> branch = branches.at (t);
			DBranch_t data = branch->GetBranch ();
			std::ostringstream os;
			std::ostringstream i;
			i << data.m_ni;
			std::string ni = i.str ();
			for (uint32_t i = 0; i < (4 - ni.size ()); i++)
				{
					os << " ";
				}
			os << ni << " ";

			std::ostringstream f;
			f << data.m_nf;
			std::string nf = f.str ();
			for (uint32_t i = 0; i < (4 - nf.size ()); i++)
				{
					os << " ";
				}
			os << nf << " ";

			std::string area = Format (data.m_area, 0, 2);
			std::string zona = Format (data.m_zona, 0, 2);

			std::string r = Format(data.m_r, 0, 9);
			std::string x = Format(data.m_x, 0, 9);
			std::string bsh = Format(data.m_bsh * 2, 0, 9);
			std::string tap;
			if (data.m_tap == 1)
				{
					tap = Format(0, 0, 7);
				}
			else
				{
					tap = Format(data.m_tap, 0, 7);
				}

			std::string line_rating1 = Format(data.m_line_rating1, 0, 5);
			std::string line_rating2 = Format(data.m_line_rating2, 0, 5);
			std::string line_rating3 = Format(data.m_line_rating3, 0, 5);
			std::string t_ctrl = Format(data.m_t_ctrl, 0, 4);
			std::string def = Format(data.m_def, 0, 7);
			std::string passo = Format(data.m_passo, 0, 6);
			std::string ctrl_min = Format (data.m_ctrl_min, 0, 7);
			std::string ctrl_max = Format (data.m_ctrl_max, 0, 6);

			std::string tapmin = Format (data.m_tapmin, 0, 6);
			std::string tapmax = Format (data.m_tapmax, 0, 6);
			std::string c = Format (t + 1, 0, 5);

			os << area << " " << zona << " " << data.m_circ << " " << data.m_tipo << " " << r << " "
				 << x << " " << bsh << " " << line_rating1 << " " << line_rating2 << " " << line_rating3
				 << " " << t_ctrl << " " << data.m_side << " " << tap << " " << def << " " << tapmin
				 << " " << tapmax << " " << passo << " " << ctrl_min << " " << ctrl_max << " " << c;

		 	otp << os.str () << std::endl;
		}
  otp << "-999" << std::endl;

  return otp.str ();
}

std::string
IOUtils::WriteCnz (Ptr<Graph> graph, double getBase)
{
	std::ostringstream otp;
	otp << "BUS DATA FOLLOWS                            " << graph->GetNumBus () << " ITEMS" << std::endl;

	for (uint32_t i = 1; i <= graph->GetNumBus (); i++)
		{
			Ptr<Bus> bus = graph->GetBus (i);

			DBus_t data = bus->GetBus ();

			std::ostringstream os;
			std::ostringstream n;
			n << data.m_nex;
			std::string nex = n.str ();
			for (uint32_t i = 0; i < (4 - nex.size ()); i++)
				{
					os << " ";
				}
			os << nex;

			double area = (data.m_area == 1 ? 0 : data.m_area);

			DoubleValue vC, aC, pG, qG;
			bus->GetAttribute ("VCalc", vC);
			bus->GetAttribute ("ACalc", aC);
			bus->GetAttribute ("PCalc", pG);
			bus->GetAttribute ("QCalc", qG);

			std::string v = Format (vC.Get (), 4, 6);
			std::string a = Format (aC.Get () * 180 / M_PI, 4, 6);
			std::string pc = Format (data.m_pc * getBase, 6, 8);
			std::string qc = Format (data.m_qc * getBase, 6, 8);
			std::string pg = Format (pG.Get () * getBase, 6, 8);
			std::string qg = Format (qG.Get () * getBase, 6, 7);

			std::string base_kV = Format (data.m_base_kV, 5, 7);
			std::string vg_o;
			if (data.m_vg_o != 0.65)
				{
					vg_o = Format (data.m_vg_o, 5, 6);;
				}
			else
				{
					vg_o = Format (0, 5, 6);
				}
			std::string max_l;
			std::string min_l;
			double max = data.m_qgmax * getBase;
			if (max == 999)
				{
					max = 0;
				}
			double min = data.m_qgmin * getBase;
			if (min == -999)
				{
					min = 0;
				}
			switch (data.m_tipo)
			{
				case 3:
				case 2:
					max_l = Format (max, 5, 7);
					min_l = Format (min, 5, 7);
					break;
				case 1:
					max_l = Format (data.m_vmax, 5, 7);
					min_l = Format (data.m_vmin, 5, 7);
					break;
				default:
					max_l = Format (0, 5, 7);
					min_l = Format (0, 5, 7);
					break;
			}

			std::string gsh = Format (data.m_gsh, 5, 7);
			std::string bsh = Format (data.m_bsh, 5, 7);
			std::string ctrl_rem = Format (data.m_ctrl_rem, 0, 4);

			os << " " << data.m_nome << " " << area << " " << Format (data.m_zona, 0, 2) << " "
					<< Format (data.m_tipo, 0, 2) << " " << v << " " << a << " " << pc << " " << qc
					<< " " << pg << " " << qg <<  " " << base_kV << " " << vg_o << " " << max_l
					<< " " << min_l << " " << gsh << " " << bsh << " " << ctrl_rem << " " << Format(nex, 4);
			otp << os.str () << std::endl;
		}
	otp << "-999" << std::endl
								<< "BRANCH DATA FOLLOWS                         " << graph->GetNumBranch () << " ITEMS" << std::endl;


	std::vector< Ptr<Branch> > conn = graph->GetBranches ();
	for (uint32_t t = 0; t < conn.size (); t++)
		{
			Ptr<Branch> branch = conn.at (t);
			DBranch_t data = branch->GetBranch ();
			std::ostringstream os;
			std::ostringstream i;
			i << data.m_ni;
			std::string ni = i.str ();
			for (uint32_t i = 0; i < (4 - ni.size ()); i++)
				{
					os << " ";
				}
			os << ni << " ";

			std::ostringstream f;
			f << data.m_nf;
			std::string nf = f.str ();
			for (uint32_t i = 0; i < (4 - nf.size ()); i++)
				{
					os << " ";
				}
			os << nf << " ";

			std::string area = Format (data.m_area, 0, 2);
			std::string zona = Format (data.m_zona, 0, 2);

			std::string r = Format(data.m_r, 8, 9);
			std::string x = Format(data.m_x, 8, 9);
			std::string bsh = Format(data.m_bsh * 2, 8, 9);
			std::string tap;
			if (data.m_tap == 1)
				{
					tap = Format(0, 5, 7);
				}
			else
				{
					tap = Format(data.m_tap, 5, 7);
				}

			std::string line_rating1 = Format(data.m_line_rating1, 0, 5);
			std::string line_rating2 = Format(data.m_line_rating2, 0, 5);
			std::string line_rating3 = Format(data.m_line_rating3, 0, 5);
			std::string t_ctrl = Format(data.m_t_ctrl, 0, 4);
			std::string def = Format(data.m_def, 5, 7);
			std::string passo = Format(data.m_passo, 4, 6);
			std::string ctrl_min = Format (data.m_ctrl_min, 4, 6);
			std::string ctrl_max = Format (data.m_ctrl_max, 4, 6);

			std::string tapmin = Format (data.m_tapmin, 4, 6);
			std::string tapmax = Format (data.m_tapmax, 4, 6);
			std::string c = Format (t + 1, 0, 5);

			os << area << " " << zona << " " << data.m_circ << " " << data.m_tipo << " " << r << " "
				 << x << " " << bsh << " " << line_rating1 << " " << line_rating2 << " "
				 << line_rating3 << " " << t_ctrl << " " << data.m_side << " " << tap << " "
				 << def << " " << tapmin << " " << tapmax << " " << passo << " " << ctrl_min << " "
				 << ctrl_max << " " << c
			;
			otp << os.str () << std::endl;
		}
	otp << "-999" << std::endl;

	return otp.str ();
}

void
IOUtils::WriteFile (std::string filename, uint32_t choose, Ptr<Graph> graph, double base)
{
	std::ofstream myfile;
  myfile.open (filename.c_str ());
  if (myfile.is_open ())
  	{
  		std::string data;
		  switch (choose)
		  	{
		  		case 0:
		  			data = WriteCdf (graph, base);
		  			break;
		  		case 1:
		  			data = WriteCnz (graph, base);
		  			break;
		  	}

		  myfile << data;
      myfile.close ();
  	}
}


}
