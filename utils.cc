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

#include "utils.h"

#include <ns3/log.h>
#include <ns3/assert.h>

#include "branch.h"
#include "bus.h"

#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <sstream>

using namespace boost;
using namespace std;

namespace ns3
{

NS_LOG_COMPONENT_DEFINE ("Utils");

uint32_t
Sts_t::BusFromNex (uint32_t nex)
{
  uint32_t sb = buses.size ();
  for (uint32_t i = 0; i < sb; i++)
    {
      DBus_t bus = buses.at (i);
      if (bus.m_nex == nex)
        {
          return bus.m_nin;
        }
    }

  return 0;
}

Utils::Utils ()
{
}

Utils::~Utils ()
{
}

TypeId
Utils::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::Utils")
    .SetParent<Object> ()
    .AddConstructor<Utils> ()
  ;

  return tid;
}

Sts_t
Utils::GetSts (void) const
{
  return m_s;
}

Sts_t
Utils::Read (string filename)
{
  ifstream file (filename.c_str ());
  if (!file)
    {
      return m_s;
    }

  m_s.m_posSlack = 0;
  m_s.m_hasTape = false;
  // Get baseMVA.
  boost::iostreams::filtering_istream in;
  in.push(file);

  // Get base MVA
  string titleCdf;
  std::getline (in, titleCdf);
  if (titleCdf.substr (0, 4).compare ("TAPE") == 0)
    {
      std::getline (in, titleCdf);
      m_s.m_hasTape = true;
    }

  m_s.m_baseMVA = atof (titleCdf.substr (31, 5).c_str ());
  NS_ASSERT_MSG ((m_s.m_baseMVA != 0), "Base MVA is equal to zero");

  // Search string 'BUS DATA FOLLOWS'.
  int i = 0;
  while (true)
    {
      string line;
      std::getline (in, line);
      if (line.substr (0,16).compare ("BUS DATA FOLLOWS") == 0)
        {
          break;
        }
      i++;
    }
  // Leitura e condicionamento dos dados de barra.
  uint32_t ibus = 0;
  while (true)
    {
      string line;
      std::getline (in, line);
      if (line.substr (0, 4).compare ("-999") == 0)
        {
          break;
        }

      ibus++;
      // Número externo da barra:
      DBus_t bus;
      bus.m_nex = atoi (line.substr(0,4).c_str ());

      // Número interno da barra:
      bus.m_nin = ibus;
      // Identificação alfanumérica da barra:
      bus.m_nome = line.substr (5,12);
      // Área:
      bus.m_area = atoi (line.substr (18,2).c_str ());
      // Zona para cálculo de perdas:
      bus.m_zona = atof (line.substr (21,2).c_str ());

      // Tipo da barra:
      bus.m_tipo = atoi (line.substr (24,2).c_str ()); // De acordo com a Figura 3.
      if (bus.m_tipo == 3)
        {
          m_s.m_posSlack = ibus;
        }

      if (bus.m_tipo < 0 || bus.m_tipo > 3)
        {
          bus.m_tipo = 0;
          NS_LOG_WARN ("Tipo da barra definido fora dos padrões CDF, considerou-se que esta barra seja do tipo (PQ)");
        }

      //// Magnitude de tensão:
      bus.m_v = atof (line.substr (27,5).c_str ());
      if (bus.m_v < VMIN_MIN)
        {
          bus.m_v = 1;
          NS_LOG_WARN ("Magnitude de tensão definida abaixo do valor mínimo para divergência automática do caso.");
        }
      else if (bus.m_v > VMAX_MAX)
        {
          bus.m_v = 1;
          NS_LOG_WARN ("Magnitude de tensão definida acima do valor máximo para divergência automática do caso.");
        }

      // Ângulo de fase de tensão:
      bus.m_ang = atof (line.substr (33,6).c_str ())* (M_PI / 180.0);
      // Carga ativa:
      bus.m_pc = atof (line.substr (40,8).c_str ()) / m_s.m_baseMVA;
      // Carga reativa:
      bus.m_qc = atof (line.substr (49,8).c_str ()) / m_s.m_baseMVA;
      // Geração ativa:
      bus.m_pg = atof (line.substr (58,8).c_str ()) / m_s.m_baseMVA;
      // Geração reativa:
      bus.m_qg = atof (line.substr (67,7).c_str ()) / m_s.m_baseMVA;
      // Vbase:
      bus.m_base_kV = atof (line.substr (76,6).c_str ());
      // Magnitude de tensão desejada para barras controladas por geração de reativos ou transformadores:
      bus.m_vg_o = atof (line.substr (84,5).c_str ());

      if (bus.m_tipo == 0)
        {
          if (bus.m_vg_o < VMIN_MIN)
            {
              bus.m_vg_o = 1;
              NS_LOG_INFO ("Magnitude da tensão desejada definida abaixo do valor mínimo para divergência automática do caso.");
            }
          else if (bus.m_vg_o > VMAX_MAX)
            {
              bus.m_vg_o = 1;
              NS_LOG_INFO ("Magnitude da tensão desejada definida acima do valor máximo para divergência automática do caso.");
            }
        }

      if (bus.m_tipo == 3 || bus.m_tipo == 2)
        {
          // Geração reativa máxima:
          bus.m_qgmax = atof (line.substr (90,7).c_str ()) / m_s.m_baseMVA;
          // Geração reativa mínima:
          bus.m_qgmin = atof (line.substr (98,7).c_str ()) / m_s.m_baseMVA;

          if (bus.m_qgmax == 0 && bus.m_qgmin == 0)
          	{
          		bus.m_qgmax = 999.9999;
          		bus.m_qgmin = -999.9999;
          	}
          // Magnitude de tensão máxima:
          bus.m_vmax = VMAX;
          // Magnitude de tensão mínima:
          bus.m_vmin = VMIN;
        }
      else if(bus.m_tipo == 1)
        {
          // Geração reativa máxima:
          bus.m_qgmax = 0.0;
          // Geração reativa mínima:
          bus.m_qgmin = 0.0;

          // Magnitude de tensão máxima:
          bus.m_vmax = atof (line.substr (90,7).c_str ());
          if (bus.m_vmax == 0)
            {
              bus.m_vmax = VMAX;
              NS_LOG_WARN ("Limite máximo da magnitude de tensão não definido.");
            }
          else if (bus.m_vmax > VMAX)
            {
              bus.m_vmax = VMAX;
              NS_LOG_WARN ("Limite máximo da magnitude de tensão definido acima do valor máximo padrão");
            }
          else if (bus.m_vmax < VMIN)
            {
              bus.m_vmax = VMAX;
              NS_LOG_WARN ("Limite máximo da magnitude de tensão definido abaixo do valor mínimo padrão");
            }
          // Magnitude de tensão mínima:
          bus.m_vmin = atof (line.substr (98,7).c_str ());
          if (bus.m_vmin == 0)
            {
              bus.m_vmin = VMIN;
              NS_LOG_WARN ("Limite mínimo da magnitude de tensão não definido.");
            }
          else if (bus.m_vmin < VMIN)
            {
              bus.m_vmin = VMIN;
              NS_LOG_INFO ("Limite mínimo da magnitude de tensão definido abaixo do valor mínimo padrão.");
            }
          else if (bus.m_vmin > VMAX)
            {
              bus.m_vmin = VMIN;
              NS_LOG_INFO ("Limite mínimo da magnitude de tensão definido acima do valor máximo padrão.");
            }
        }
      else
        {
          // Geração reativa máxima:
          bus.m_qgmax = 0.0;
          // Geração reativa mínima:
          bus.m_qgmin = 0.0;

          // Magnitude de tensão máxima:
          bus.m_vmax = VMAX;
          // Magnitude de tensão mínima:
          bus.m_vmin = VMIN;
        }

      // Condutância do shunt de barra:
      bus.m_gsh = atof (line.substr (106,7).c_str ());
      // Susceptância do shunt de barra:
      //bus.bsh = atof (line.substr (114,8).c_str()) / s.m_baseMVA;
      bus.m_bsh = atof (line.substr (114,7).c_str());
      // Barra remota controlada pela geração de reativos na barra atual:
      if (bus.m_tipo == 3 || bus.m_tipo == 2)
        {
          bus.m_ctrl_rem = atof (line.substr (123,4).c_str ());
          if (bus.m_ctrl_rem == 0)
            {
              bus.m_ctrl_rem = bus.m_nex;
            }
        }
      else
        {
          bus.m_ctrl_rem = 0;
        }
      // Add bus
      m_s.buses.push_back (bus);
    }
  // Criação de dados de barra adicionais.
  m_s.m_nb = ibus;
  uint32_t cont1 = 0;
  uint32_t cont2 = 0;
  uint32_t noSlack = 0;
  m_s.m_npv = 0;
  m_s.m_npq = 0;
  for (uint32_t ibus = 0; ibus < m_s.m_nb; ibus++)
    {
      DBus_t bus = m_s.buses.at (ibus);
      if (bus.m_tipo == 3)
        {
          continue;
        }

      if (bus.m_tipo == 2)
        {
          bus.m_ordPV = cont1++;
          bus.m_posPV = ibus;
          bus.m_ord = noSlack++;
          m_s.m_npv++;
        }
      else if (bus.m_tipo == 1 || bus.m_tipo == 0)
        {
          bus.m_ordPQ = cont2++;
          bus.m_posPQ = ibus;
          bus.m_ord = noSlack++;
          m_s.m_npq++;
        }

      m_s.buses.at (ibus) = bus;
    }

  if (m_s.m_posSlack == 0)
    {
      DBus_t bus = m_s.buses.at (0);
      bus.m_tipo = 3;
      m_s.m_posSlack = 1;
      m_s.buses.at (0) = bus;
      NS_LOG_INFO ("Não foi encontrada nenhuma barra do tipo SLACK. A primeira barra do tipo PV dos dados de barra foi transformada em barra SLACK");
    }

  // Busca pela string 'BRANCH DATA FOLLOWS'.
  while (true)
    {
      string line; std::getline (in, line);
      if (line.substr (0,19).compare ("BRANCH DATA FOLLOWS") == 0)
        {
          break;
        }
    }
  //// Leitura dos dados de circuito.
  uint32_t ibranch = 0;
  while (true)
    {
      string line; std::getline (in, line);
      if (line.substr (0,4).compare ("-999") == 0)
        {
          break;
        }

      ibranch++;
      DBranch_t branch;
      // Da barra (ou barra do tap para trafos):
      branch.m_ni = atoi (line.substr (0,4).c_str ());
      // Para barra (ou barra da impedância para trafos):
      branch.m_nf = atoi (line.substr (5,4).c_str ());
      // Área:
      branch.m_area = atoi (line.substr (10,2).c_str ());
      // Zona para cálculo de perdas:
      branch.m_zona = atoi (line.substr (13,2).c_str ());
      // Quantidade de circuitos paralelos:
      branch.m_circ = atoi (line.substr (16, 1).c_str ());
      // Tipo do ramo:
      branch.m_tipo = atoi (line.substr(18, 1).c_str ());
      if (branch.m_tipo < 0 || branch.m_tipo > 4)
        {
          branch.m_tipo = 0;
          NS_LOG_WARN ("Tipo do ramo definido fora dos padrões CDF; considerou-se que este ramo seja do tipo '0'.");
        }
      // Resistência série:
      branch.m_r = atof (line.substr(20, 8).c_str ());
      // Reatância série:
      branch.m_x = atof (line.substr (29,9).c_str ());

      NS_ASSERT_MSG (branch.m_x != 0, "O ramo possui reatância igual a zero.");
      // Susceptância shunt de linha:
      branch.m_bsh = atof (line.substr (40,10).c_str () ) / 2; // De acordo com o texto do artigo.

      // Line ratings:
      branch.m_line_rating1 = atof (line.substr (50,5).c_str ());
      branch.m_line_rating2 = atof(line.substr (56,5).c_str ());
      branch.m_line_rating3 = atof(line.substr (62,5).c_str ());

      // Terminal com magnitude de tensão controlada por tap:
      branch.m_t_ctrl = atof(line.substr (68,4).c_str ());
      if (branch.m_t_ctrl == 0 && branch.m_tipo == 2)
        {
          branch.m_t_ctrl = branch.m_nf;
          NS_LOG_INFO ("Barra controlada pelo trafo do ramo não definida.");
        }

      // Localização da barra controlada pelo transformador:
      branch.m_side = atof(line.substr (73, 1).c_str ());
      //**************************************************//
      // Implementar testes de verificação do valor lido! //
      //**************************************************//

      // Tap:
      branch.m_tap = atof(line.substr (75, 6).c_str ());
      if (branch.m_tipo == 0)
        {
          branch.m_tap = 1.0;
        }
      else if (branch.m_tap == 0)
        {
          branch.m_tap = 1.0;
          NS_LOG_WARN ("Tap do trafo do ramo não definido.");
        }
      // Ângulo de defasamento:
      branch.m_def = atof (line.substr (83,6).c_str () )* (M_PI/180);
      if (branch.m_tipo == 2 || branch.m_tipo == 3)
        {
          // Tap mínimo:
          branch.m_tapmin = atof(line.substr (90,6).c_str ());
          if (branch.m_tapmin == 0)
            {
              branch.m_tapmin = TAPMIN;
              NS_LOG_WARN ("Limite mínimo do tap do trafo não definido.");
            }
          else if (branch.m_tapmin < TAPMIN_MIN || branch.m_tapmin > 1)
            {
              branch.m_tapmin = TAPMIN_MIN;
              NS_LOG_WARN ("Limite mínimo do tap do trafo do ramo definido fora dos valores operacionais.");
            }
          // Tap máximo:
          branch.m_tapmax = atof(line.substr (97,6).c_str ());
          if (branch.m_tapmax == 0)
            {
              branch.m_tapmax = TAPMAX;
              NS_LOG_WARN ("Limite máximo do tap do trafo não definido.");
            }
          else if (branch.m_tapmax > TAPMAX_MAX || branch.m_tapmax < 1)
            {
              branch.m_tapmax = TAPMAX;
              NS_LOG_WARN ("Limite máximo do tap do trafo definido fora dos valores operacionais.");
            }
          // Tamanho do passo entre posições intermediárias do tap:
          branch.m_passo = atof(line.substr (104,6).c_str ());
          if (branch.m_passo == 0)
            {
              branch.m_passo = STEP;
              NS_LOG_WARN ("Tamanho do passo entre posições intermediárias do tap do trafo não definido.");
            }
          // Limite mínimo da grandeza (magnitude de tensão, fluxo de potência reativa ou fluxo de potência ativa) controlada:
          branch.m_ctrl_min = atof(line.substr (112,6).c_str ());
          // Limite máximo da grandeza (magnitude de tensão, fluxo de potência reativa ou fluxo de potência ativa) controlada:
          branch.m_ctrl_max = atof(line.substr (119,4).c_str ());
        }
      else
        {
          // Tap mínimo:
          branch.m_tapmin = 0;
          // Tap máximo:
          branch.m_tapmax = 0;
          // Tamanho do passo entre posições intermediárias do tap:
          branch.m_passo = 0;
          // Limite mínimo da grandeza (magnitude de tensão, fluxo de potência reativa ou fluxo de potência ativa) controlada:
          branch.m_ctrl_min = 0;
          // Limite máximo da grandeza (magnitude de tensão, fluxo de potência reativa ou fluxo de potência ativa) controlada:
          branch.m_ctrl_max = 0;
        }
      m_s.branches.push_back (branch);
    }
  // Criação de dados de circuito adicionais.
  m_s.m_nc = ibranch;
  cont1 = 0;
  m_s.m_ntap = 0;
  for (uint32_t ibranch=0; ibranch < m_s.m_nc; ibranch++)
    {
      DBranch_t branch = m_s.branches.at (ibranch);
      branch.m_g =  branch.m_r / (pow (branch.m_r, 2) + pow (branch.m_x, 2));
      branch.m_b = -branch.m_x / (pow (branch.m_r, 2) + pow (branch.m_x, 2));

      if (branch.m_tipo == 2)
        {
          branch.m_ordtap = cont1++;
          branch.m_postap = ibranch;
          m_s.m_ntap++;
        }

      m_s.branches.at (ibranch) = branch;
    }

  // Busca pela string 'LOSS ZONES FOLLOW'.
  while (true)
    {
      string line; std::getline (in, line);
      if (line.substr (0,17).compare ("LOSS ZONES FOLLOW") == 0)
        {
          break;
        }
    }

  // Leitura dos dados de zonas para o cálculo de perdas.
  while (true)
    {
      string line; std::getline (in, line);
      if (line.substr (0,3).compare ("-99") == 0)
        {
          break;
        }
    }

  // Busca pela string 'INTERCHANGE DATA FOLLOWS'.
  while (true)
    {
      string line; std::getline (in, line);
      if (line.substr (0,24).compare ("INTERCHANGE DATA FOLLOWS") == 0)
        {
          break;
        }
    }

  // Leitura dos dados de intercâmbio de potência ativa entre áreas.
  while (true)
    {
      string line; std::getline (in, line);
      if (line.substr (0,2).compare ("-9") == 0)
        {
          break;
        }
    }

  // Busca pela string 'TIE LINES FOLLOW'.
  while (true)
    {
      string line; std::getline (in, line);
      if (line.substr (0,16).compare ("TIE LINES FOLLOW") == 0)
        {
          break;
        }
    }

  // Leitura dos dados de 'tie lines'.
  while (true)
    {
      string line; std::getline (in, line);
      if (line.substr (0,4).compare ("-999") == 0)
        {
          break;
        }
    }

  // Busca pela string 'END OF DATA'.
  while (true)
    {
      string line; std::getline (in, line);
      if (line.substr (0,11).compare ("END OF DATA") == 0)
        {
          break;
        }
    }
  // Fecha o arquivo de entrada.
  if (file.is_open ())
    {
      file.close ();
    }

  return m_s;
}

std::string
Utils::Format (double v, uint32_t p, uint32_t maxS)
{
  std::string sep = " ";
  ostringstream os;
  os.unsetf (std::ios::floatfield);
  os.precision (p);
  os.setf ( std::ios::fixed, std:: ios::floatfield ); // floatfield set to fixed
  os << v << sep;

  std::string s;
  for (uint32_t i = 0; i < (maxS - os.str ().size ()); i++)
    {
      s.append (sep);
    }
  ostringstream t;
  t << (s + os.str ());

  return ( t.str () );
}

}
