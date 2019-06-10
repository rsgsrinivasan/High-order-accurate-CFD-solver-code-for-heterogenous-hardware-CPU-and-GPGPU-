//#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <string>
#include "functions.h"
#include <openmpi-x86_64/mpi.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using std::cout;
using std::cin;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::istringstream;

void roe_average()
{

	double sqr_b1, sqr_i, sqr_b3, sqr_c0, sqr_c2, sqr_d4, sqr_d5;
	sqr_i = sqrt(FLOW[i].rho);
	sqr_b1 = sqrt(FLOW[node[i].n_n[1]].rho);
	sqr_b3 = sqrt(FLOW[node[i].n_n[3]].rho);
	sqr_c0 = sqrt(FLOW[node[i].n_n[0]].rho);
	sqr_c2 = sqrt(FLOW[node[i].n_n[2]].rho);
	sqr_d4 = sqrt(FLOW[node[i].n_n[4]].rho);
	sqr_d5 = sqrt(FLOW[node[i].n_n[5]].rho);


	ROE_AVER[i].roe_rho_ip = sqr_b1 * sqr_i;
	ROE_AVER[i].roe_u_ip = (sqr_b1*FLOW[node[i].n_n[1]].u + sqr_i * FLOW[i].u) / (sqr_b1 + sqr_i);
	ROE_AVER[i].roe_v_ip = (sqr_b1*FLOW[node[i].n_n[1]].v + sqr_i * FLOW[i].v) / (sqr_b1 + sqr_i);
	ROE_AVER[i].roe_w_ip = (sqr_b1*FLOW[node[i].n_n[1]].w + sqr_i * FLOW[i].w) / (sqr_b1 + sqr_i);
	ROE_AVER[i].roe_h_ip = (sqr_b1*(((FLOW[node[i].n_n[1]].p / 0.4) + 0.5*FLOW[node[i].n_n[1]].rho * (FLOW[node[i].n_n[1]].u * FLOW[node[i].n_n[1]].u + FLOW[node[i].n_n[1]].v * FLOW[node[i].n_n[1]].v + FLOW[node[i].n_n[1]].w * FLOW[node[i].n_n[1]].w) + FLOW[node[i].n_n[1]].p) / FLOW[node[i].n_n[1]].rho) + \
		sqr_i*(((FLOW[i].p / 0.4) + 0.5*FLOW[i].rho * (FLOW[i].u * FLOW[i].u + FLOW[i].v * FLOW[i].v + FLOW[i].w * FLOW[i].w) + FLOW[i].p) / FLOW[i].rho)) / (sqr_b1 + sqr_i);

	ROE_AVER[i].roe_rho_im = sqr_b3 * sqr_i;
	ROE_AVER[i].roe_u_im = (sqr_b3*FLOW[node[i].n_n[3]].u + sqr_i * FLOW[i].u) / (sqr_b3 + sqr_i);
	ROE_AVER[i].roe_v_im = (sqr_b3*FLOW[node[i].n_n[3]].v + sqr_i * FLOW[i].v) / (sqr_b3 + sqr_i);
	ROE_AVER[i].roe_w_im = (sqr_b3*FLOW[node[i].n_n[3]].w + sqr_i * FLOW[i].w) / (sqr_b3 + sqr_i);
	ROE_AVER[i].roe_h_im = (sqr_b3*(((FLOW[node[i].n_n[3]].p / 0.4) + 0.5*FLOW[node[i].n_n[3]].rho * (FLOW[node[i].n_n[3]].u * FLOW[node[i].n_n[3]].u + FLOW[node[i].n_n[3]].v * FLOW[node[i].n_n[3]].v + FLOW[node[i].n_n[3]].w * FLOW[node[i].n_n[3]].w) + FLOW[node[i].n_n[3]].p) / FLOW[node[i].n_n[3]].rho) + \
		sqr_i*(((FLOW[i].p / 0.4) + 0.5*FLOW[i].rho * (FLOW[i].u * FLOW[i].u + FLOW[i].v * FLOW[i].v + FLOW[i].w * FLOW[i].w) + FLOW[i].p) / FLOW[i].rho)) / (sqr_b3 + sqr_i);

	ROE_AVER[i].roe_rho_jp = sqr_c0 * sqr_i;
	ROE_AVER[i].roe_u_jp = (sqr_c0*FLOW[node[i].n_n[0]].u + sqr_i * FLOW[i].u) / (sqr_c0 + sqr_i);
	ROE_AVER[i].roe_v_jp = (sqr_c0*FLOW[node[i].n_n[0]].v + sqr_i * FLOW[i].v) / (sqr_c0 + sqr_i);
	ROE_AVER[i].roe_w_jp = (sqr_c0*FLOW[node[i].n_n[0]].w + sqr_i * FLOW[i].w) / (sqr_c0 + sqr_i);
	ROE_AVER[i].roe_h_jp = (sqr_c0*(((FLOW[node[i].n_n[0]].p / 0.4) + 0.5*FLOW[node[i].n_n[0]].rho * (FLOW[node[i].n_n[0]].u * FLOW[node[i].n_n[0]].u + FLOW[node[i].n_n[0]].v * FLOW[node[i].n_n[0]].v + FLOW[node[i].n_n[0]].w * FLOW[node[i].n_n[0]].w) + FLOW[node[i].n_n[0]].p) / FLOW[node[i].n_n[0]].rho) + \
		sqr_i*(((FLOW[i].p / 0.4) + 0.5*FLOW[i].rho * (FLOW[i].u * FLOW[i].u + FLOW[i].v * FLOW[i].v + FLOW[i].w * FLOW[i].w) + FLOW[i].p) / FLOW[i].rho)) / (sqr_c0 + sqr_i);

	ROE_AVER[i].roe_rho_jm = sqr_c2 * sqr_i;
	ROE_AVER[i].roe_u_jm = (sqr_c2*FLOW[node[i].n_n[2]].u + sqr_i * FLOW[i].u) / (sqr_c2 + sqr_i);
	ROE_AVER[i].roe_v_jm = (sqr_c2*FLOW[node[i].n_n[2]].v + sqr_i * FLOW[i].v) / (sqr_c2 + sqr_i);
	ROE_AVER[i].roe_w_jm = (sqr_c2*FLOW[node[i].n_n[2]].w + sqr_i * FLOW[i].w) / (sqr_c2 + sqr_i);
	ROE_AVER[i].roe_h_jm = (sqr_c2*(((FLOW[node[i].n_n[2]].p / 0.4) + 0.5*FLOW[node[i].n_n[2]].rho * (FLOW[node[i].n_n[2]].u * FLOW[node[i].n_n[2]].u + FLOW[node[i].n_n[2]].v * FLOW[node[i].n_n[2]].v + FLOW[node[i].n_n[2]].w * FLOW[node[i].n_n[2]].w) + FLOW[node[i].n_n[2]].p) / FLOW[node[i].n_n[2]].rho) + \
		sqr_i*(((FLOW[i].p / 0.4) + 0.5*FLOW[i].rho * (FLOW[i].u * FLOW[i].u + FLOW[i].v * FLOW[i].v + FLOW[i].w * FLOW[i].w) + FLOW[i].p) / FLOW[i].rho)) / (sqr_c2 + sqr_i);

	ROE_AVER[i].roe_rho_kp = sqr_d4 * sqr_i;
	ROE_AVER[i].roe_u_kp = (sqr_d4*FLOW[node[i].n_n[4]].u + sqr_i * FLOW[i].u) / (sqr_d4 + sqr_i);
	ROE_AVER[i].roe_v_kp = (sqr_d4*FLOW[node[i].n_n[4]].v + sqr_i * FLOW[i].v) / (sqr_d4 + sqr_i);
	ROE_AVER[i].roe_w_kp = (sqr_d4*FLOW[node[i].n_n[4]].w + sqr_i * FLOW[i].w) / (sqr_d4 + sqr_i);
	ROE_AVER[i].roe_h_kp = (sqr_d4*(((FLOW[node[i].n_n[4]].p / 0.4) + 0.5*FLOW[node[i].n_n[4]].rho * (FLOW[node[i].n_n[4]].u * FLOW[node[i].n_n[4]].u + FLOW[node[i].n_n[4]].v * FLOW[node[i].n_n[4]].v + FLOW[node[i].n_n[4]].w * FLOW[node[i].n_n[4]].w) + FLOW[node[i].n_n[4]].p) / FLOW[node[i].n_n[4]].rho) + \
		sqr_i*(((FLOW[i].p / 0.4) + 0.5*FLOW[i].rho * (FLOW[i].u * FLOW[i].u + FLOW[i].v * FLOW[i].v + FLOW[i].w * FLOW[i].w) + FLOW[i].p) / FLOW[i].rho)) / (sqr_d4 + sqr_i);

	ROE_AVER[i].roe_rho_km = sqr_d5 * sqr_i;
	ROE_AVER[i].roe_u_km = (sqr_d5*FLOW[node[i].n_n[5]].u + sqr_i * FLOW[i].u) / (sqr_d5 + sqr_i);
	ROE_AVER[i].roe_v_km = (sqr_d5*FLOW[node[i].n_n[5]].v + sqr_i * FLOW[i].v) / (sqr_d5 + sqr_i);
	ROE_AVER[i].roe_w_km = (sqr_d5*FLOW[node[i].n_n[5]].w + sqr_i * FLOW[i].w) / (sqr_d5 + sqr_i);
	ROE_AVER[i].roe_h_km = (sqr_d5*(((FLOW[node[i].n_n[5]].p / 0.4) + 0.5*FLOW[node[i].n_n[5]].rho * (FLOW[node[i].n_n[5]].u * FLOW[node[i].n_n[5]].u + FLOW[node[i].n_n[5]].v * FLOW[node[i].n_n[5]].v + FLOW[node[i].n_n[5]].w * FLOW[node[i].n_n[5]].w) + FLOW[node[i].n_n[5]].p) / FLOW[node[i].n_n[5]].rho) + \
		sqr_i*(((FLOW[i].p / 0.4) + 0.5*FLOW[i].rho * (FLOW[i].u * FLOW[i].u + FLOW[i].v * FLOW[i].v + FLOW[i].w * FLOW[i].w) + FLOW[i].p) / FLOW[i].rho)) / (sqr_d5 + sqr_i);

	if (i <= sd_node && j == 0)
	{
		diver[i][0].u = FLOW[i].u;
		diver[i][0].v = FLOW[i].v;
		diver[i][0].e = FLOW[i].e;
	}

}


