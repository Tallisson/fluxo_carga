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

#include "jacobian.h"

#include "ns3/log.h"
#include "ns3/uinteger.h"

#include <iostream>

namespace ns3
{

NS_LOG_COMPONENT_DEFINE ("Jacobian");
NS_OBJECT_ENSURE_REGISTERED (Jacobian);

Jacobian::Jacobian ()
{

}

Jacobian::~Jacobian ()
{
}


TypeId
Jacobian::GetTypeId (void)
{
	static TypeId tid = TypeId ("ns3::Jacobian")
		      .SetParent<Object> ()
		      .AddConstructor<Jacobian> ()
	;

	return tid;
}

/*void
Jacobian::SetGraph (Ptr<Graph> graph)
{
	m_graph = graph;

	UintegerValue nPQ, nPV;
	m_graph->GetAttribute("NumPQ", nPQ);
	m_graph->GetAttribute("NumPQ", nPV);

	uint64_t num = nPQ.Get () * 2 + nPV.Get ();
	SetMatrix (num, num);
}*/

mat
Jacobian::GetMatrix () const
{
	return m_matrix;
}

void
Jacobian::SetMatrix (uint64_t numRows, uint64_t numCols)
{
	m_matrix = zeros(numRows, numCols);
}

mat
Jacobian::GetJ1 () const
{
	return m_j1;
}
void
Jacobian::SetJ1 (uint64_t numRows, uint64_t numCols)
{
	m_j1 = zeros(numRows, numCols);
}

mat
Jacobian::GetJ2 () const
{
	return m_j2;
}
void
Jacobian::SetJ2 (uint64_t numRows, uint64_t numCols)
{
	m_j2 = zeros(numRows, numCols);
}

mat
Jacobian::GetJ3 () const
{
	return m_j3;
}
void
Jacobian::SetJ3 (uint64_t numRows, uint64_t numCols)
{
	m_j3 = zeros(numRows, numCols);
}

mat
Jacobian::GetJ4 () const
{
	return m_j1;
}
void
Jacobian::SetJ4 (uint64_t numRows, uint64_t numCols)
{
	m_j4 = zeros(numRows, numCols);
}

void
Jacobian::CalcDPk (Ptr<Graph> graph)
{
	for (uint32_t i = 0; i < graph->GetNumBus (); i++)
		{
			Ptr<Bus> busK = graph->GetBus (i + 1);
			if (busK->GetType () != Bus::SLACK)
				{
					std::vector<Ptr<Branch> > branches = busK->GetBranches ();
					std::vector<Ptr<Bus> > neighbors = busK->GetNeighbors ();
					for (uint32_t j = 0; j < branches.size (); j++)
						{
								DBranch_t dataBranch = branches.at (j)->GetBranch ();
								Ptr<Bus> busM = neighbors.at (j);

								DoubleValue vK, vM, aK, aM;
								busK->GetAttribute ("VCalc", vK);
								busM->GetAttribute ("VCalc", vM);
								busK->GetAttribute ("ACalc", aK);
								busM->GetAttribute ("ACalc", aM);

								double theta_km = aK.Get () - aM.Get ();
								double tap = dataBranch.m_tap;

								/*
								 * dPk em relação a 'ak'
								 * Jac(k-1, k-1) = -(1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
								 * (s.branch.g(km)*sin(akm)-s.branch.b(km)*cos(akm)) + Jac(k-1, k-1)
								 */
								uint32_t k = busK->GetBus ().m_ord;
								//if (busM->GetType () != Bus::SLACK)
									{
										m_matrix (k, k) = - (1 / tap) * vK.Get () * vM.Get () *
																			( dataBranch.m_g * sin (theta_km) - dataBranch.m_b * cos (theta_km) ) +
																			m_matrix (k, k);
									}

								/*
								 * dPk em relação a 'am' (exceto quando m for a barra slack).
								 *
								 * Jac(k-1, m-1) = (1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
								 * (s.branch.g(km)*sin(akm)-s.branch.b(km)*cos(akm)) + Jac(k-1, m-1);
								 */
								if (busM->GetType () != Bus::SLACK)
									{
										uint32_t m = busM->GetBus ().m_ord;
										m_matrix (k, m) = (1 / tap) * vK.Get () * vM.Get () *
																					( dataBranch.m_g * sin (theta_km) - dataBranch.m_b * cos (theta_km) ) +
																					m_matrix (k, m);
									}


								/*
								 * dPk em relação a 'vk'.
								 * I:
								 * Jac(k-1, s.nb-1+s.bus.ordPQ(k)) =
								 * -2*s.branch.g(km)*(1/s.branch.tap(km))*s.bus.v(k) +
								 * (1/s.branch.tap(km))*s.bus.v(m)*(s.branch.g(km)*cos(akm)+s.branch.b(km)*sin(akm)) +
								 * Jac(k-1, s.nb-1+s.bus.ordPQ(k));
								 *
								 * II:
								 * Jac(k-1, s.nb-1+s.bus.ordPQ(k)) =
								 * -2*s.branch.g(km)*s.bus.v(k) +
								 * (1/s.branch.tap(km))*s.bus.v(m)*(s.branch.g(km)*cos(akm)+s.branch.b(km)*sin(akm)) +
								 * Jac(k-1, s.nb-1+s.bus.ordPQ(k));
								 */
								if (busK->GetType () == Bus::LOAD || busK->GetType () == Bus::LOSS_CONTROL_REACT)
									{
										uint32_t index = graph->GetNumBus () - 1 + busK->GetBus ().m_ordPQ;

										if (dataBranch.m_tipo == 1 && busK->GetTap () == Bus::TAP)
											{
												m_matrix (k, index) = -2 * dataBranch.m_g * (1 / tap) * vK.Get () +
																						(1 / tap) * vM.Get () *
																						(dataBranch.m_g * cos (theta_km) + dataBranch.m_b * sin(theta_km)) +
																						m_matrix (k, index);
											}
										else
											{
												m_matrix (k, index) = -2 * dataBranch.m_g * vK.Get () +
																							(1 / tap) * vM.Get () *
																							(dataBranch.m_g * cos (theta_km) + dataBranch.m_b * sin(theta_km)) +
																							m_matrix (k, index);
											}
									}

								/*
								 * dPk em relação a 'vm'.
								 * Jac(k-1, s.nb-1+s.bus.ordPQ(m)) =
								 * (1/s.branch.tap(km))*s.bus.v(k)*
								 * (s.branch.g(km)*cos(akm)+s.branch.b(km)*sin(akm)) + Jac(k-1, s.nb-1+s.bus.ordPQ(m));
								 */
								if (busM->GetType () == Bus::LOAD ||
											busM->GetType () == Bus::LOSS_CONTROL_REACT)
									{
										uint32_t index = graph->GetNumBus () - 1 + busM->GetBus ().m_ordPQ;
										m_matrix (k, index) = (1 / tap) * vK.Get () *
																				( dataBranch.m_g * cos (theta_km) + dataBranch.m_b * sin (theta_km) ) +
																				m_matrix (k, index);
									}
						}
				}
		}
}

void
Jacobian::CalcDQk (Ptr<Graph> graph)
{
	for (uint32_t i = 0; i < graph->GetNumBus (); i++)
		{
			Ptr<Bus> busK = graph->GetBus (i + 1);
			if (busK->GetType () == Bus::LOAD ||
						busK->GetType () == Bus::LOSS_CONTROL_REACT)
				{
					uint32_t indexK = graph->GetNumBus () - 1 + busK->GetBus ().m_ordPQ;
					std::vector<Ptr<Branch> > branches = busK->GetBranches ();
					std::vector<Ptr<Bus> > neighbors = busK->GetNeighbors ();

					for (uint32_t j = 0; j < branches.size (); j++)
						{
							Ptr<Branch> branch = branches.at (j);
							DBranch_t dataBranch = branch->GetBranch ();
							Ptr<Bus> busM = neighbors.at (j);

							DoubleValue vK, vM, aK, aM;
							busK->GetAttribute ("VCalc", vK);
							busM->GetAttribute ("VCalc", vM);
							busK->GetAttribute ("ACalc", aK);
							busM->GetAttribute ("ACalc", aM);

							double theta_km = aK.Get () - aM.Get ();
							double tap = dataBranch.m_tap;

							uint32_t k = busK->GetBus ().m_ord;
							/*
							 * dQk em relaçao a 'ak'.
							 * Jac(s.nb-1+s.bus.ordPQ(k), k-1) =
							 * (1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
							 * (s.branch.b(km)*sin(akm)+s.branch.g(km)*cos(akm)) + Jac(s.nb-1+s.bus.ordPQ(k), k-1)
							 */
							if (busK->GetType () != Bus::SLACK)
								{
									m_matrix (indexK, k) = (1 / tap) * vK.Get () * vM.Get () *
																				( dataBranch.m_b * sin (theta_km) + dataBranch.m_g * cos (theta_km) ) +
																				m_matrix (indexK, k);
								}

								/*
								 * dQk em relaçao a 'am' (exceto quando m for a barra slack).
								 * Jac(s.nb-1+s.bus.ordPQ(k), m-1) =
								 * -(1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
								 * (s.branch.g(km)*cos(akm)+s.branch.b(km)*sin(akm)) + Jac(s.nb-1+s.bus.ordPQ(k), m-1)
								 */
							if (busM->GetType () != Bus::SLACK)
								{
									uint32_t m = busM->GetBus ().m_ord;
									m_matrix (indexK, m) = -(1 / tap) * vK.Get () * vM.Get () *
																				( dataBranch.m_g * cos (theta_km) + dataBranch.m_b * sin (theta_km) ) +
																				m_matrix(indexK, m);
								}

							/*
							 * dQk em relaçao a 'vk'
							 * I:
							 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k)) =
							 * 2*((1/s.branch.tap(km)^2)*s.branch.b(km)+s.branch.bsh(km))*s.bus.v(k) -
							 * (1/s.branch.tap(km))*s.bus.v(m)*(s.branch.b(km)*cos(akm)-s.branch.g(km)*sin(akm)) +
							 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k))
							 *
							 * II:
							 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k)) =
							 * 2*(s.branch.b(km)+s.branch.bsh(km))*s.bus.v(k) -
							 * (1/s.branch.tap(km))*s.bus.v(m)*(s.branch.b(km)*cos(akm)-s.branch.g(km)*sin(akm)) +
							 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k));
							 */
							/*if (busM->GetType () == Bus::LOAD ||
										busM->GetType () == Bus::LOSS_CONTROL_REACT)*/
								{
									if (dataBranch.m_tipo == 1 && busK->GetTap () == Bus::TAP)
										{
											m_matrix (indexK, indexK) = 2 * ( pow((1 / tap), 2) * dataBranch.m_b + dataBranch.m_bsh ) *
																			vK.Get () - (1 / tap) * vM.Get () *
																			(dataBranch.m_b * cos (theta_km) - dataBranch.m_g * sin(theta_km)) +
																			m_matrix (indexK, indexK);
										}
									else
										{
											m_matrix (indexK, indexK) = 2 * ( dataBranch.m_b + dataBranch.m_bsh ) *
																			vK.Get () - (1 / tap) * vM.Get () *
																			(dataBranch.m_b * cos (theta_km) - dataBranch.m_g * sin(theta_km)) +
																			m_matrix (indexK, indexK);
										}
								}

								/* dQk em relacao a 'vm'.
								 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(m)) = -(1/s.branch.tap(km))*s.bus.v(k)*
								 * (s.branch.b(km)*cos(akm)-s.branch.g(km)*sin(akm)) + Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(m))
								 */
							if (busM->GetType () == Bus::LOAD ||
										busM->GetType () == Bus::LOSS_CONTROL_REACT)
								{
									uint32_t indexM = graph->GetNumBus () - 1 + busM->GetBus ().m_ordPQ;
									m_matrix(indexK, indexM) = -(1 / tap) * vK.Get () *
																	( dataBranch.m_b * cos (theta_km) - dataBranch.m_g * sin (theta_km) ) +
																	m_matrix(indexK, indexM);
								}
						 }
					/*
					 * dQk em relaçao a 'vk' (continuação)
					 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k)) =
					 * 2*s.bus.bsh(k)*s.bus.v(k) + Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k))
					 */
					DoubleValue vK;
					busK->GetAttribute("VCalc", vK);
					m_matrix (indexK, indexK) =
							(2 * busK->GetBus ().m_bsh * vK.Get ()) + m_matrix (indexK, indexK);
				}
		}
}

mat
Jacobian::CalcJac (Ptr<Graph> graph)
{
	CalcDPk (graph);
	CalcDQk (graph);

	return GetMatrix ();
}

void
Jacobian::Zeros (void)
{
	uint32_t size = m_matrix.n_cols;
	m_matrix = zeros(size, size);
}

void
Jacobian::Zeros (uint32_t numRows, uint32_t numCols)
{
	m_matrix.clear ();
	m_matrix = zeros(numRows, numCols);
}

vec
Jacobian::SolveSys (vec b)
{
	return vec(inv (m_matrix) * -b);
}

mat
Jacobian::GetJqv (void)
{
  mat j1; // = *(m_lf->GetJ1 ()->GetData ()) * -1;
  mat j2; //= *(m_lf->GetJ2 ()->GetData ()) * -1;
  mat j3; //= *(m_lf->GetJ3 ()->GetData ()) * -1;
  mat j4; // = *(m_lf->GetJ4 ()->GetData ()) * -1;

  mat jqv = j4 - (j3 * inv (j1) * j2);

	return inv (m_matrix);
}

}

