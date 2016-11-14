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

#include "load-flow.h"

#include "q-control.h"

#include "ns3/log.h"
#include "ns3/double.h"
#include "ns3/uinteger.h"
#include "ns3/io-utils.h"

#include <iostream>
#include <armadillo>
#include <map>

#include <sstream>

using namespace arma;

namespace ns3 {

NS_LOG_COMPONENT_DEFINE("LoadFlow");
NS_OBJECT_ENSURE_REGISTERED(LoadFlow);

LoadFlow::LoadFlow(void) :
	m_precision(0.00001), m_iter(0), m_maxIter(10)
{
	m_qControl = CreateObject<QControl> ();
	m_graph = CreateObject<Graph> ();
	m_mismatches = CreateObject<Mismatch> ();
	m_jac = CreateObject<Jacobian> ();
	m_verbose = true;
	m_dir = "/home/thiago/pwns-ns-3/resultados-fluxo";

	m_report = CreateObject<Report> ();
	m_report->SetFileCrt ("/home/thiago/pwns-ns-3/resultados-fluxo/crt.txt");
	m_report->SetFileLTotal ("/home/thiago/pwns-ns-3/resultados-fluxo/ltotal.txt");
	m_report->SetFileState ("/home/thiago/pwns-ns-3/resultados-fluxo/state.txt");
	m_report->SetFileV ("/home/thiago/pwns-ns-3/resultados-fluxo/v.txt");
	m_report->SetFileQg ("/home/thiago/pwns-ns-3/resultados-fluxo/qg.txt");
}

LoadFlow::~LoadFlow(void)
{
	m_graph = NULL;
	m_qControl = NULL;
	m_mismatches = NULL;
	m_jac = NULL;
	m_report = NULL;
}

TypeId
LoadFlow::GetTypeId(void)
{
	static TypeId tid =
			TypeId("ns3::LoadFlow").SetParent<Object>().AddConstructor<LoadFlow>();

	return tid;
}

Ptr<Graph>
LoadFlow::GetGraph(void) const
{
	return m_graph;
}

/*void
 LoadFlow::SetGraph(Ptr<Graph> graph)
 {
 m_graph = graph;
 }*/

Ptr<Control>
LoadFlow::GetQControl(void) const
{
	return m_qControl;
}

void
LoadFlow::SetQControl(Ptr<Control> qControl)
{
	m_qControl = qControl;
}

void
LoadFlow::InitX0(void)
{
	// Passo 2 - Arbitrar um ponto inicial:
	for (uint32_t i = 0; i < m_graph->GetNumBus(); i++)
		{
			Ptr<Bus> bus = m_graph->GetBus(i + 1);
			if (bus->GetType() != Bus::SLACK) {
				double angSlack =
						m_graph->GetBus(m_graph->GetPosSlack())->GetBus().m_ang;
				bus->SetAttribute("ACalc", DoubleValue(angSlack));
			}

			if (bus->GetType() == Bus::LOAD) {
				bus->SetAttribute("VCalc", DoubleValue(1.0));
			}
		}
}

void
LoadFlow::Prepare(std::string cdf)
{
	// Definir topologia completa da rede.
	Ptr<Utils> utils = CreateObject<Utils>();
	m_sts = utils->Read(cdf);

	m_graph->SetPosSlack(m_sts.m_posSlack);
	for (uint32_t i = 0; i < m_sts.buses.size(); i++)
		{
			DBus_t busData = m_sts.buses.at(i);
			Ptr<Bus> bus = CreateObject<Bus>();
			bus->SetBus(busData);
			m_graph->AddBus(bus);
		}

	for (uint32_t i = 0; i < m_sts.branches.size(); i++)
		{
			DBranch_t dataBranch = m_sts.branches.at(i);
			Ptr<Branch> branch = CreateObject<Branch>();
			branch->SetBranch(dataBranch);
			m_graph->Assoc(branch);

			if (dataBranch.m_tipo == 0)
				{
					dataBranch.m_tap = 1.0;
					branch->SetBranch(dataBranch);

					Ptr<Bus> busNi = m_graph->GetBus(dataBranch.m_ni);
					Ptr<Bus> busNf = m_graph->GetBus(dataBranch.m_nf);
					busNi->SetTap(Bus::IMP);
					busNf->SetTap(Bus::IMP);
				}
			else if (dataBranch.m_tipo == 1)
				{
					Ptr<Bus> busNi = m_graph->GetBus(dataBranch.m_ni);
					Ptr<Bus> busNf = m_graph->GetBus(dataBranch.m_nf);
					busNi->SetTap(Bus::TAP);
					busNf->SetTap(Bus::IMP);
				}
		}
}

void
LoadFlow::InitJ(void)
{
	uint32_t nPQ = m_graph->GetNumPQ ();
	uint32_t nPV = m_graph->GetNumPV ();
	uint32_t size = m_graph->GetNumPQ () * 2 + m_graph->GetNumPV();
	uint32_t numB = m_graph->GetNumBus ();
	m_jac->SetMatrix(size, size);

	m_jac->SetJ1 ((nPQ + nPV), (nPQ + nPV));
	m_jac->SetJ2 ((nPQ + nPV), numB);
	m_jac->SetJ3 (numB, (nPQ + nPV));
	m_jac->SetJ4 (numB, numB);

}

void
LoadFlow::Reset(void)
{
	uint32_t size = m_graph->GetNumPQ() * 2 + m_graph->GetNumPV();
	m_jac->Zeros(size, size);
	m_b = zeros<vec>(size);
}

void
LoadFlow::Execute()
{
	/*
	 * Calcular o vetor dos mismatches, a matriz Jacobiana e resolver
	 * o sistema de equações lineares:
	 */
	InitX0 ();
	for(uint32_t i = 1; i <= m_graph->GetNumBus(); i++)
		{
			Ptr<Bus> bus = m_graph->GetBus(i);
			bus->Print();
		}

	InitJ ();

	m_b = m_mismatches->CalcMismatches(m_graph);

	bool execute = false;
	uint32_t nextIter, nextCrt;
	nextCrt = 0;
	do//while (execute)
		{
			nextIter = 0;
			m_iter = 0;
			while (nextIter == 0)
				{
					mat m = m_jac->CalcJac (m_graph);
					vec dx = m_jac->SolveSys (m_b);
					//std::cout << "dx = " << std::endl << dx << std::endl;

					//std::cout << "Dx: \n" << dx << std::endl;
					// Atualizar 'a' e 'V':

					for (uint32_t i = 0; i < m_graph->GetNumBus(); i++)
						{
							Ptr<Bus> bus = m_graph->GetBus(i + 1);
							if (bus->GetType() != Bus::SLACK)
								{
									DoubleValue ang;
									bus->GetAttribute("ACalc", ang);
									double angD = ang.Get();

									angD += dx(bus->GetBus().m_ord);
									//std::cout << "Ângulo Anterior " << ang.Get () << ", Atualização = " << dx (bus->GetBus ().m_ord) << " Novo Ângulo = " << angD << std::endl;
									bus->SetAttribute("ACalc", DoubleValue(angD));
								}

							if (bus->GetType() != Bus::SLACK && bus->GetType() != Bus::GENERATION)
								{
									DoubleValue v;
									bus->GetAttribute("VCalc", v);
									double vD = v.Get();

									uint32_t ind = m_graph->GetNumBus() - 1 + bus->GetBus().m_ordPQ;
									vD += dx(ind);
									//std::cout << "Mag. Tensão Anterior " << v.Get () << ", Atualização = " << dx (ind) << ", Nova Mag. Tensão = " << vD << std::endl;
									bus->SetAttribute("VCalc", DoubleValue(vD));
								}
						}
						bool crt = false;
						if (m_qControl != NULL)
							{
								crt = m_qControl->GetObject<QControl> ()->DoRestore (m_graph);
								if (crt == true)
									{
										InitJ ();
									}
							}


						m_iter++;

						// Qlim
						std::cout << "Iter " << m_iter << std::endl;
						if (m_qControl != NULL)
							{
								crt = m_qControl->GetObject<QControl> ()->DoControl (m_graph);
								if (crt == true)
									{
										InitJ ();
									}
							}

						/*
						 * Cálculo do vetor dos mismatches com os valores corrigidos de
						 * 'a' e 'V':
						 */
						m_b = m_mismatches->CalcMismatches(m_graph);
						//std::cout << "Erros: \n" << b;
						std::cout << "+++++++++++++++++++++++++++++++++++++++++\n";
						for (uint32_t i = 0; i < m_graph->GetNumBus(); i++)
							{
								Ptr<Bus> bus = m_graph->GetBus(i + 1);
								//bus->Print();
							}
						// Teste de convergência:

						double maxB = max(abs(m_b));

						if (maxB <= m_precision)
							{
								nextIter = 1;
							}
						else
							{
								m_jac->Zeros();
								nextIter = 0;

								// Critério de saída do laço:
								if (m_iter == m_maxIter)
									{
										nextIter = 2;
									}
							}
							if(nextIter != 0)
								{
									CalcLosses();
									m_report->StoreData (m_graph, m_sts.m_baseMVA);
									m_report->StoreL (m_totalL);
								}
				}

			if(m_vControl != NULL)
				{
					execute = false;
					if (nextCrt < 20)
						{
							execute = m_vControl->DoControl (m_jac->GetJqv (), m_graph);
						}
				}
			if (execute == true)
				{
					nextCrt++;
				}
		} while (execute);

	if (nextIter == 1)
		{
			std::cout << "O método de Newton-Raphson convergiu em " << m_iter
					<< " iterações" << std::endl;
		}
	else
		{
			std::cout << "O número máximo de iterações foi atingido e o método de Newton-Raphson não convergiu..."
								<< std::endl;
		}

	for (uint32_t i = 0; i < m_graph->GetNumBus(); i++)
		{
			Ptr<Bus> bus = m_graph->GetBus(i + 1);
			DoubleValue v;
			bus->GetAttribute("VCalc", v);
			if ((bus->GetType () == Bus::LOAD || bus->GetType () == 4) && (v.Get () < 0.95 || v.Get () > 1.05))
				{
					bus->Print();
				}
		}

	std::string filename = m_dir + "/" + m_file;
	IOUtils::WriteFile (filename, 0, m_graph, m_sts.m_baseMVA);
	m_report->Save (m_graph);
}

void
LoadFlow::CalcLosses(void)
{
	std::vector<double> p_km, p_mk, q_km, q_mk, pL, losses;
	m_totalL = 0;
	std::vector<Ptr<Branch> > branches = m_graph->GetBranches ();
	for (uint32_t i = 0; i < m_graph->GetNumBranch (); i++)
		{
			Ptr<Branch> branch = branches.at (i);
			Ptr<Bus> busK = m_graph->GetBus (branch->GetBranch ().m_ni);
			Ptr<Bus> busM = m_graph->GetBus (branch->GetBranch ().m_nf);

			DoubleValue vK, vM, aK, aM;
			busK->GetAttribute("VCalc", vK);
			busK->GetAttribute("ACalc", aK);
			busM->GetAttribute("VCalc", vM);
			busM->GetAttribute("ACalc", aM);

			p_km.push_back (branch->CalcPkmL (vK, vM, aK, aM));
			p_mk.push_back (branch->CalcPmkL (vK, vM, aK, aM));
			q_km.push_back (branch->CalcQkmL (vK, vM, aK, aM));
			q_mk.push_back (branch->CalcQmkL (vK, vM, aK, aM));

			double l = branch->CalcL (vK, vM, aK, aM);
			losses.push_back (l);
			m_totalL += l;
		}

	m_totalL *= m_sts.m_baseMVA;
	if (m_verbose)
		{
			/*std::cout << IOUtils::Format ("Pkm") << IOUtils::Format ("Pmk")
								<< IOUtils::Format ("Perdas") << IOUtils::Format ("Qkm")
								<< IOUtils::Format ("Qmk") << std::endl;
			for (uint32_t i = 0; i < p_km.size (); i++)
				{
					std::cout << IOUtils::Format (p_km.at (i)) << IOUtils::Format (p_mk.at (i))
										<< IOUtils::Format (losses.at (i)) << IOUtils::Format (q_km.at (i))
										<< IOUtils::Format (q_mk.at (i)) << std::endl;
				}*/
		}

	NS_LOG_INFO("Perda Total: " << m_totalL);
	//std::cout << "Perda Total: " << (m_totalL) << std::endl;
}

void
LoadFlow::SetPrecision(double precision)
{
	m_precision = precision;
}

void
LoadFlow::SetDir (std::string dir)
{
	m_dir = dir;
}

void
LoadFlow::SetFile (std::string file)
{
	m_file = file;
}

void
LoadFlow::SetVControl (Ptr<VControl> vControl)
{
	m_vControl = vControl;
}

}
