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

void weno_solver_0()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(ROE_AVER[i].roe_h_jp - 0.5*(ROE_AVER[i].roe_u_jp * ROE_AVER[i].roe_u_jp + ROE_AVER[i].roe_v_jp * ROE_AVER[i].roe_v_jp + ROE_AVER[i].roe_w_jp * ROE_AVER[i].roe_w_jp)));

	//	if (hk == 0)
	{
		zeta_xi = metric[i].eta_xjp;
		zeta_yi = metric[i].eta_yjp;
		zeta_zi = metric[i].eta_zjp;
	}

	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * ROE_AVER[i].roe_u_jp + ky_bar * ROE_AVER[i].roe_v_jp + kz_bar * ROE_AVER[i].roe_w_jp;
	phi_sq = 0.5*0.4*(ROE_AVER[i].roe_u_jp * ROE_AVER[i].roe_u_jp + ROE_AVER[i].roe_v_jp * ROE_AVER[i].roe_v_jp + ROE_AVER[i].roe_w_jp * ROE_AVER[i].roe_w_jp);
	alpha = ROE_AVER[i].roe_rho_jp / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[i].roe_rho_jp * roe_a[i]);

	r_eigen_Qjp[0][0] = kx_bar;
	r_eigen_Qjp[0][1] = ky_bar;
	r_eigen_Qjp[0][2] = kz_bar;
	r_eigen_Qjp[0][3] = alpha;
	r_eigen_Qjp[0][4] = alpha;

	r_eigen_Qjp[1][0] = kx_bar * ROE_AVER[i].roe_u_jp;
	r_eigen_Qjp[1][1] = ky_bar * ROE_AVER[i].roe_u_jp - kz_bar * ROE_AVER[i].roe_rho_jp;
	r_eigen_Qjp[1][2] = kz_bar * ROE_AVER[i].roe_u_jp + ky_bar * ROE_AVER[i].roe_rho_jp;
	r_eigen_Qjp[1][3] = alpha * (ROE_AVER[i].roe_u_jp + kx_bar * roe_a[i]);
	r_eigen_Qjp[1][4] = alpha * (ROE_AVER[i].roe_u_jp - kx_bar * roe_a[i]);

	r_eigen_Qjp[2][0] = kx_bar * ROE_AVER[i].roe_v_jp + kz_bar * ROE_AVER[i].roe_rho_jp;
	r_eigen_Qjp[2][1] = ky_bar * ROE_AVER[i].roe_v_jp;
	r_eigen_Qjp[2][2] = kz_bar * ROE_AVER[i].roe_v_jp - kx_bar * ROE_AVER[i].roe_rho_jp;
	r_eigen_Qjp[2][3] = alpha * (ROE_AVER[i].roe_v_jp + ky_bar * roe_a[i]);
	r_eigen_Qjp[2][4] = alpha * (ROE_AVER[i].roe_v_jp - ky_bar * roe_a[i]);

	r_eigen_Qjp[3][0] = kx_bar * ROE_AVER[i].roe_w_jp - ky_bar * ROE_AVER[i].roe_rho_jp;
	r_eigen_Qjp[3][1] = ky_bar * ROE_AVER[i].roe_w_jp + kx_bar * ROE_AVER[i].roe_rho_jp;
	r_eigen_Qjp[3][2] = kz_bar * ROE_AVER[i].roe_w_jp;
	r_eigen_Qjp[3][3] = alpha * (ROE_AVER[i].roe_w_jp + kz_bar * roe_a[i]);
	r_eigen_Qjp[3][4] = alpha * (ROE_AVER[i].roe_w_jp - kz_bar * roe_a[i]);

	r_eigen_Qjp[4][0] = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_jp * (kz_bar*ROE_AVER[i].roe_v_jp - ky_bar * ROE_AVER[i].roe_w_jp);
	r_eigen_Qjp[4][1] = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_jp * (kx_bar*ROE_AVER[i].roe_w_jp - kz_bar * ROE_AVER[i].roe_u_jp);
	r_eigen_Qjp[4][2] = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_jp * (ky_bar*ROE_AVER[i].roe_u_jp - kx_bar * ROE_AVER[i].roe_v_jp);
	r_eigen_Qjp[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qjp[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qjp[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*ROE_AVER[i].roe_v_jp - ky_bar * ROE_AVER[i].roe_w_jp) / ROE_AVER[i].roe_rho_jp);
	l_eigen_Qjp[0][1] = kx_bar * 0.4*ROE_AVER[i].roe_u_jp / (roe_a[i] * roe_a[i]);
	l_eigen_Qjp[0][2] = (kx_bar*0.4*ROE_AVER[i].roe_v_jp / (roe_a[i] * roe_a[i])) + (kz_bar / ROE_AVER[i].roe_rho_jp);
	l_eigen_Qjp[0][3] = (kx_bar*0.4*ROE_AVER[i].roe_w_jp / (roe_a[i] * roe_a[i])) - (ky_bar / ROE_AVER[i].roe_rho_jp);
	l_eigen_Qjp[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjp[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*ROE_AVER[i].roe_w_jp - kz_bar * ROE_AVER[i].roe_u_jp) / ROE_AVER[i].roe_rho_jp);
	l_eigen_Qjp[1][1] = (ky_bar*0.4*ROE_AVER[i].roe_u_jp / (roe_a[i] * roe_a[i])) - kz_bar / ROE_AVER[i].roe_rho_jp;
	l_eigen_Qjp[1][2] = ky_bar * 0.4*ROE_AVER[i].roe_v_jp / (roe_a[i] * roe_a[i]);
	l_eigen_Qjp[1][3] = (ky_bar*0.4*ROE_AVER[i].roe_w_jp / (roe_a[i] * roe_a[i])) + (kx_bar / ROE_AVER[i].roe_rho_jp);
	l_eigen_Qjp[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjp[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*ROE_AVER[i].roe_u_jp - kx_bar * ROE_AVER[i].roe_v_jp) / ROE_AVER[i].roe_rho_jp);
	l_eigen_Qjp[2][1] = kz_bar * 0.4*ROE_AVER[i].roe_u_jp / (roe_a[i] * roe_a[i]) + (ky_bar / ROE_AVER[i].roe_rho_jp);
	l_eigen_Qjp[2][2] = kz_bar * 0.4*ROE_AVER[i].roe_v_jp / (roe_a[i] * roe_a[i]) - (kx_bar / ROE_AVER[i].roe_rho_jp);
	l_eigen_Qjp[2][3] = kz_bar * 0.4*ROE_AVER[i].roe_w_jp / (roe_a[i] * roe_a[i]);
	l_eigen_Qjp[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjp[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qjp[3][1] = -beta * (0.4*ROE_AVER[i].roe_u_jp - kx_bar * roe_a[i]);
	l_eigen_Qjp[3][2] = -beta * (0.4*ROE_AVER[i].roe_v_jp - ky_bar * roe_a[i]);
	l_eigen_Qjp[3][3] = -beta * (0.4*ROE_AVER[i].roe_w_jp - kz_bar * roe_a[i]);
	l_eigen_Qjp[3][4] = beta * 0.4;

	l_eigen_Qjp[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qjp[4][1] = -beta * (0.4*ROE_AVER[i].roe_u_jp + kx_bar * roe_a[i]);
	l_eigen_Qjp[4][2] = -beta * (0.4*ROE_AVER[i].roe_v_jp + ky_bar * roe_a[i]);
	l_eigen_Qjp[4][3] = -beta * (0.4*ROE_AVER[i].roe_w_jp + kz_bar * roe_a[i]);
	l_eigen_Qjp[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		//	if (hk == 0)
		{
			Qj_iplus[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k] = l_eigen_Qjp[m][0] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][0].ip + l_eigen_Qjp[m][1] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][1].ip + l_eigen_Qjp[m][2] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][2].ip + l_eigen_Qjp[m][3] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][3].ip + l_eigen_Qjp[m][4] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][4].ip;
			Qj_iplus[node[node[i].n_n[0]].n_n[0]][k] = l_eigen_Qjp[m][0] * Q[node[node[i].n_n[0]].n_n[0]][0].ip + l_eigen_Qjp[m][1] * Q[node[node[i].n_n[0]].n_n[0]][1].ip + l_eigen_Qjp[m][2] * Q[node[node[i].n_n[0]].n_n[0]][2].ip + l_eigen_Qjp[m][3] * Q[node[node[i].n_n[0]].n_n[0]][3].ip + l_eigen_Qjp[m][4] * Q[node[node[i].n_n[0]].n_n[0]][4].ip;
			Qj_iplus[node[i].n_n[0]][k] = l_eigen_Qjp[m][0] * Q[node[i].n_n[0]][0].ip + l_eigen_Qjp[m][1] * Q[node[i].n_n[0]][1].ip + l_eigen_Qjp[m][2] * Q[node[i].n_n[0]][2].ip + l_eigen_Qjp[m][3] * Q[node[i].n_n[0]][3].ip + l_eigen_Qjp[m][4] * Q[node[i].n_n[0]][4].ip;
			Qj_iplus[i][k] = l_eigen_Qjp[m][0] * Q[i][0].ip + l_eigen_Qjp[m][1] * Q[i][1].ip + l_eigen_Qjp[m][2] * Q[i][2].ip + l_eigen_Qjp[m][3] * Q[i][3].ip + l_eigen_Qjp[m][4] * Q[i][4].ip;
			Qj_iplus[node[i].n_n[2]][k] = l_eigen_Qjp[m][0] * Q[node[i].n_n[2]][0].ip + l_eigen_Qjp[m][1] * Q[node[i].n_n[2]][1].ip + l_eigen_Qjp[m][2] * Q[node[i].n_n[2]][2].ip + l_eigen_Qjp[m][3] * Q[node[i].n_n[2]][3].ip + l_eigen_Qjp[m][4] * Q[node[i].n_n[2]][4].ip;
			Qj_iplus[node[node[i].n_n[2]].n_n[2]][k] = l_eigen_Qjp[m][0] * Q[node[node[i].n_n[2]].n_n[2]][0].ip + l_eigen_Qjp[m][1] * Q[node[node[i].n_n[2]].n_n[2]][1].ip + l_eigen_Qjp[m][2] * Q[node[node[i].n_n[2]].n_n[2]][2].ip + l_eigen_Qjp[m][3] * Q[node[node[i].n_n[2]].n_n[2]][3].ip + l_eigen_Qjp[m][4] * Q[node[node[i].n_n[2]].n_n[2]][4].ip;
			Qj_iplus[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k] = l_eigen_Qjp[m][0] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][0].ip + l_eigen_Qjp[m][1] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][1].ip + l_eigen_Qjp[m][2] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][2].ip + l_eigen_Qjp[m][3] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][3].ip + l_eigen_Qjp[m][4] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][4].ip;
		}

		m++;
	}
}


void weno_solver_1()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(ROE_AVER[i].roe_h_ip - 0.5*(ROE_AVER[i].roe_u_ip * ROE_AVER[i].roe_u_ip + ROE_AVER[i].roe_v_ip * ROE_AVER[i].roe_v_ip + ROE_AVER[i].roe_w_ip * ROE_AVER[i].roe_w_ip)));

	zeta_xi = metric[i].zeta_xip;
	zeta_yi = metric[i].zeta_yip;
	zeta_zi = metric[i].zeta_zip;

	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * ROE_AVER[i].roe_u_ip + ky_bar * ROE_AVER[i].roe_v_ip + kz_bar * ROE_AVER[i].roe_w_ip;
	phi_sq = 0.5*0.4*(ROE_AVER[i].roe_u_ip * ROE_AVER[i].roe_u_ip + ROE_AVER[i].roe_v_ip * ROE_AVER[i].roe_v_ip + ROE_AVER[i].roe_w_ip * ROE_AVER[i].roe_w_ip);
	alpha = ROE_AVER[i].roe_rho_ip / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[i].roe_rho_ip * roe_a[i]);

	r_eigen_Qip[0][0] = kx_bar;
	r_eigen_Qip[0][1] = ky_bar;
	r_eigen_Qip[0][2] = kz_bar;
	r_eigen_Qip[0][3] = alpha;
	r_eigen_Qip[0][4] = alpha;

	r_eigen_Qip[1][0] = kx_bar * ROE_AVER[i].roe_u_ip;
	r_eigen_Qip[1][1] = ky_bar * ROE_AVER[i].roe_u_ip - kz_bar * ROE_AVER[i].roe_rho_ip;
	r_eigen_Qip[1][2] = kz_bar * ROE_AVER[i].roe_u_ip + ky_bar * ROE_AVER[i].roe_rho_ip;
	r_eigen_Qip[1][3] = alpha * (ROE_AVER[i].roe_u_ip + kx_bar * roe_a[i]);
	r_eigen_Qip[1][4] = alpha * (ROE_AVER[i].roe_u_ip - kx_bar * roe_a[i]);

	r_eigen_Qip[2][0] = kx_bar * ROE_AVER[i].roe_v_ip + kz_bar * ROE_AVER[i].roe_rho_ip;
	r_eigen_Qip[2][1] = ky_bar * ROE_AVER[i].roe_v_ip;
	r_eigen_Qip[2][2] = kz_bar * ROE_AVER[i].roe_v_ip - kx_bar * ROE_AVER[i].roe_rho_ip;
	r_eigen_Qip[2][3] = alpha * (ROE_AVER[i].roe_v_ip + ky_bar * roe_a[i]);
	r_eigen_Qip[2][4] = alpha * (ROE_AVER[i].roe_v_ip - ky_bar * roe_a[i]);

	r_eigen_Qip[3][0] = kx_bar * ROE_AVER[i].roe_w_ip - ky_bar * ROE_AVER[i].roe_rho_ip;
	r_eigen_Qip[3][1] = ky_bar * ROE_AVER[i].roe_w_ip + kx_bar * ROE_AVER[i].roe_rho_ip;
	r_eigen_Qip[3][2] = kz_bar * ROE_AVER[i].roe_w_ip;
	r_eigen_Qip[3][3] = alpha * (ROE_AVER[i].roe_w_ip + kz_bar * roe_a[i]);
	r_eigen_Qip[3][4] = alpha * (ROE_AVER[i].roe_w_ip - kz_bar * roe_a[i]);

	r_eigen_Qip[4][0] = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_ip * (kz_bar*ROE_AVER[i].roe_v_ip - ky_bar * ROE_AVER[i].roe_w_ip);
	r_eigen_Qip[4][1] = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_ip * (kx_bar*ROE_AVER[i].roe_w_ip - kz_bar * ROE_AVER[i].roe_u_ip);
	r_eigen_Qip[4][2] = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_ip * (ky_bar*ROE_AVER[i].roe_u_ip - kx_bar * ROE_AVER[i].roe_v_ip);
	r_eigen_Qip[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qip[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qip[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*ROE_AVER[i].roe_v_ip - ky_bar * ROE_AVER[i].roe_w_ip) / ROE_AVER[i].roe_rho_ip);
	l_eigen_Qip[0][1] = kx_bar * 0.4*ROE_AVER[i].roe_u_ip / (roe_a[i] * roe_a[i]);
	l_eigen_Qip[0][2] = (kx_bar*0.4*ROE_AVER[i].roe_v_ip / (roe_a[i] * roe_a[i])) + (kz_bar / ROE_AVER[i].roe_rho_ip);
	l_eigen_Qip[0][3] = (kx_bar*0.4*ROE_AVER[i].roe_w_ip / (roe_a[i] * roe_a[i])) - (ky_bar / ROE_AVER[i].roe_rho_ip);
	l_eigen_Qip[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qip[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*ROE_AVER[i].roe_w_ip - kz_bar * ROE_AVER[i].roe_u_ip) / ROE_AVER[i].roe_rho_ip);
	l_eigen_Qip[1][1] = (ky_bar*0.4*ROE_AVER[i].roe_u_ip / (roe_a[i] * roe_a[i])) - kz_bar / ROE_AVER[i].roe_rho_ip;
	l_eigen_Qip[1][2] = ky_bar * 0.4*ROE_AVER[i].roe_v_ip / (roe_a[i] * roe_a[i]);
	l_eigen_Qip[1][3] = (ky_bar*0.4*ROE_AVER[i].roe_w_ip / (roe_a[i] * roe_a[i])) + (kx_bar / ROE_AVER[i].roe_rho_ip);
	l_eigen_Qip[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qip[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*ROE_AVER[i].roe_u_ip - kx_bar * ROE_AVER[i].roe_v_ip) / ROE_AVER[i].roe_rho_ip);
	l_eigen_Qip[2][1] = kz_bar * 0.4*ROE_AVER[i].roe_u_ip / (roe_a[i] * roe_a[i]) + (ky_bar / ROE_AVER[i].roe_rho_ip);
	l_eigen_Qip[2][2] = kz_bar * 0.4*ROE_AVER[i].roe_v_ip / (roe_a[i] * roe_a[i]) - (kx_bar / ROE_AVER[i].roe_rho_ip);
	l_eigen_Qip[2][3] = kz_bar * 0.4*ROE_AVER[i].roe_w_ip / (roe_a[i] * roe_a[i]);
	l_eigen_Qip[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qip[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qip[3][1] = -beta * (0.4*ROE_AVER[i].roe_u_ip - kx_bar * roe_a[i]);
	l_eigen_Qip[3][2] = -beta * (0.4*ROE_AVER[i].roe_v_ip - ky_bar * roe_a[i]);
	l_eigen_Qip[3][3] = -beta * (0.4*ROE_AVER[i].roe_w_ip - kz_bar * roe_a[i]);
	l_eigen_Qip[3][4] = beta * 0.4;

	l_eigen_Qip[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qip[4][1] = -beta * (0.4*ROE_AVER[i].roe_u_ip + kx_bar * roe_a[i]);
	l_eigen_Qip[4][2] = -beta * (0.4*ROE_AVER[i].roe_v_ip + ky_bar * roe_a[i]);
	l_eigen_Qip[4][3] = -beta * (0.4*ROE_AVER[i].roe_w_ip + kz_bar * roe_a[i]);
	l_eigen_Qip[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		Qi_iplus[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k] = l_eigen_Qip[m][0] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][0].ip + l_eigen_Qip[m][1] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][1].ip + l_eigen_Qip[m][2] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][2].ip + l_eigen_Qip[m][3] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][3].ip + l_eigen_Qip[m][4] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][4].ip;
		Qi_iplus[node[node[i].n_n[1]].n_n[1]][k] = l_eigen_Qip[m][0] * Q[node[node[i].n_n[1]].n_n[1]][0].ip + l_eigen_Qip[m][1] * Q[node[node[i].n_n[1]].n_n[1]][1].ip + l_eigen_Qip[m][2] * Q[node[node[i].n_n[1]].n_n[1]][2].ip + l_eigen_Qip[m][3] * Q[node[node[i].n_n[1]].n_n[1]][3].ip + l_eigen_Qip[m][4] * Q[node[node[i].n_n[1]].n_n[1]][4].ip;
		Qi_iplus[node[i].n_n[1]][k] = l_eigen_Qip[m][0] * Q[node[i].n_n[1]][0].ip + l_eigen_Qip[m][1] * Q[node[i].n_n[1]][1].ip + l_eigen_Qip[m][2] * Q[node[i].n_n[1]][2].ip + l_eigen_Qip[m][3] * Q[node[i].n_n[1]][3].ip + l_eigen_Qip[m][4] * Q[node[i].n_n[1]][4].ip;
		Qi_iplus[i][k] = l_eigen_Qip[m][0] * Q[i][0].ip + l_eigen_Qip[m][1] * Q[i][1].ip + l_eigen_Qip[m][2] * Q[i][2].ip + l_eigen_Qip[m][3] * Q[i][3].ip + l_eigen_Qip[m][4] * Q[i][4].ip;
		Qi_iplus[node[i].n_n[3]][k] = l_eigen_Qip[m][0] * Q[node[i].n_n[3]][0].ip + l_eigen_Qip[m][1] * Q[node[i].n_n[3]][1].ip + l_eigen_Qip[m][2] * Q[node[i].n_n[3]][2].ip + l_eigen_Qip[m][3] * Q[node[i].n_n[3]][3].ip + l_eigen_Qip[m][4] * Q[node[i].n_n[3]][4].ip;
		Qi_iplus[node[node[i].n_n[3]].n_n[3]][k] = l_eigen_Qip[m][0] * Q[node[node[i].n_n[3]].n_n[3]][0].ip + l_eigen_Qip[m][1] * Q[node[node[i].n_n[3]].n_n[3]][1].ip + l_eigen_Qip[m][2] * Q[node[node[i].n_n[3]].n_n[3]][2].ip + l_eigen_Qip[m][3] * Q[node[node[i].n_n[3]].n_n[3]][3].ip + l_eigen_Qip[m][4] * Q[node[node[i].n_n[3]].n_n[3]][4].ip;
		Qi_iplus[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k] = l_eigen_Qip[m][0] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][0].ip + l_eigen_Qip[m][1] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][1].ip + l_eigen_Qip[m][2] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][2].ip + l_eigen_Qip[m][3] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][3].ip + l_eigen_Qip[m][4] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][4].ip;


		m++;
	}
}

void weno_solver_2()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(ROE_AVER[i].roe_h_jm - 0.5*(ROE_AVER[i].roe_u_jm * ROE_AVER[i].roe_u_jm + ROE_AVER[i].roe_v_jm * ROE_AVER[i].roe_v_jm + ROE_AVER[i].roe_w_jm * ROE_AVER[i].roe_w_jm)));

	zeta_xi = metric[i].eta_xjm;
	zeta_yi = metric[i].eta_yjm;
	zeta_zi = metric[i].eta_zjm;


	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * ROE_AVER[i].roe_u_jm + ky_bar * ROE_AVER[i].roe_v_jm + kz_bar * ROE_AVER[i].roe_w_jm;
	phi_sq = 0.5*0.4*(ROE_AVER[i].roe_u_jm * ROE_AVER[i].roe_u_jm + ROE_AVER[i].roe_v_jm * ROE_AVER[i].roe_v_jm + ROE_AVER[i].roe_w_jm * ROE_AVER[i].roe_w_jm);
	alpha = ROE_AVER[i].roe_rho_jm / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[i].roe_rho_jm * roe_a[i]);

	r_eigen_Qjm[0][0] = kx_bar;
	r_eigen_Qjm[0][1] = ky_bar;
	r_eigen_Qjm[0][2] = kz_bar;
	r_eigen_Qjm[0][3] = alpha;
	r_eigen_Qjm[0][4] = alpha;

	r_eigen_Qjm[1][0] = kx_bar * ROE_AVER[i].roe_u_jm;
	r_eigen_Qjm[1][1] = ky_bar * ROE_AVER[i].roe_u_jm - kz_bar * ROE_AVER[i].roe_rho_jm;
	r_eigen_Qjm[1][2] = kz_bar * ROE_AVER[i].roe_u_jm + ky_bar * ROE_AVER[i].roe_rho_jm;
	r_eigen_Qjm[1][3] = alpha * (ROE_AVER[i].roe_u_jm + kx_bar * roe_a[i]);
	r_eigen_Qjm[1][4] = alpha * (ROE_AVER[i].roe_u_jm - kx_bar * roe_a[i]);

	r_eigen_Qjm[2][0] = kx_bar * ROE_AVER[i].roe_v_jm + kz_bar * ROE_AVER[i].roe_rho_jm;
	r_eigen_Qjm[2][1] = ky_bar * ROE_AVER[i].roe_v_jm;
	r_eigen_Qjm[2][2] = kz_bar * ROE_AVER[i].roe_v_jm - kx_bar * ROE_AVER[i].roe_rho_jm;
	r_eigen_Qjm[2][3] = alpha * (ROE_AVER[i].roe_v_jm + ky_bar * roe_a[i]);
	r_eigen_Qjm[2][4] = alpha * (ROE_AVER[i].roe_v_jm - ky_bar * roe_a[i]);

	r_eigen_Qjm[3][0] = kx_bar * ROE_AVER[i].roe_w_jm - ky_bar * ROE_AVER[i].roe_rho_jm;
	r_eigen_Qjm[3][1] = ky_bar * ROE_AVER[i].roe_w_jm + kx_bar * ROE_AVER[i].roe_rho_jm;
	r_eigen_Qjm[3][2] = kz_bar * ROE_AVER[i].roe_w_jm;
	r_eigen_Qjm[3][3] = alpha * (ROE_AVER[i].roe_w_jm + kz_bar * roe_a[i]);
	r_eigen_Qjm[3][4] = alpha * (ROE_AVER[i].roe_w_jm - kz_bar * roe_a[i]);

	r_eigen_Qjm[4][0] = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_jm * (kz_bar*ROE_AVER[i].roe_v_jm - ky_bar * ROE_AVER[i].roe_w_jm);
	r_eigen_Qjm[4][1] = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_jm * (kx_bar*ROE_AVER[i].roe_w_jm - kz_bar * ROE_AVER[i].roe_u_jm);
	r_eigen_Qjm[4][2] = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_jm * (ky_bar*ROE_AVER[i].roe_u_jm - kx_bar * ROE_AVER[i].roe_v_jm);
	r_eigen_Qjm[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qjm[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qjm[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*ROE_AVER[i].roe_v_jm - ky_bar * ROE_AVER[i].roe_w_jm) / ROE_AVER[i].roe_rho_jm);
	l_eigen_Qjm[0][1] = kx_bar * 0.4*ROE_AVER[i].roe_u_jm / (roe_a[i] * roe_a[i]);
	l_eigen_Qjm[0][2] = (kx_bar*0.4*ROE_AVER[i].roe_v_jm / (roe_a[i] * roe_a[i])) + (kz_bar / ROE_AVER[i].roe_rho_jm);
	l_eigen_Qjm[0][3] = (kx_bar*0.4*ROE_AVER[i].roe_w_jm / (roe_a[i] * roe_a[i])) - (ky_bar / ROE_AVER[i].roe_rho_jm);
	l_eigen_Qjm[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjm[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*ROE_AVER[i].roe_w_jm - kz_bar * ROE_AVER[i].roe_u_jm) / ROE_AVER[i].roe_rho_jm);
	l_eigen_Qjm[1][1] = (ky_bar*0.4*ROE_AVER[i].roe_u_jm / (roe_a[i] * roe_a[i])) - kz_bar / ROE_AVER[i].roe_rho_jm;
	l_eigen_Qjm[1][2] = ky_bar * 0.4*ROE_AVER[i].roe_v_jm / (roe_a[i] * roe_a[i]);
	l_eigen_Qjm[1][3] = (ky_bar*0.4*ROE_AVER[i].roe_w_jm / (roe_a[i] * roe_a[i])) + (kx_bar / ROE_AVER[i].roe_rho_jm);
	l_eigen_Qjm[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjm[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*ROE_AVER[i].roe_u_jm - kx_bar * ROE_AVER[i].roe_v_jm) / ROE_AVER[i].roe_rho_jm);
	l_eigen_Qjm[2][1] = kz_bar * 0.4*ROE_AVER[i].roe_u_jm / (roe_a[i] * roe_a[i]) + (ky_bar / ROE_AVER[i].roe_rho_jm);
	l_eigen_Qjm[2][2] = kz_bar * 0.4*ROE_AVER[i].roe_v_jm / (roe_a[i] * roe_a[i]) - (kx_bar / ROE_AVER[i].roe_rho_jm);
	l_eigen_Qjm[2][3] = kz_bar * 0.4*ROE_AVER[i].roe_w_jm / (roe_a[i] * roe_a[i]);
	l_eigen_Qjm[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjm[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qjm[3][1] = -beta * (0.4*ROE_AVER[i].roe_u_jm - kx_bar * roe_a[i]);
	l_eigen_Qjm[3][2] = -beta * (0.4*ROE_AVER[i].roe_v_jm - ky_bar * roe_a[i]);
	l_eigen_Qjm[3][3] = -beta * (0.4*ROE_AVER[i].roe_w_jm - kz_bar * roe_a[i]);
	l_eigen_Qjm[3][4] = beta * 0.4;

	l_eigen_Qjm[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qjm[4][1] = -beta * (0.4*ROE_AVER[i].roe_u_jm + kx_bar * roe_a[i]);
	l_eigen_Qjm[4][2] = -beta * (0.4*ROE_AVER[i].roe_v_jm + ky_bar * roe_a[i]);
	l_eigen_Qjm[4][3] = -beta * (0.4*ROE_AVER[i].roe_w_jm + kz_bar * roe_a[i]);
	l_eigen_Qjm[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		//	if (hk == 2)
		{
			Qj_iminus_n[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k] = l_eigen_Qjm[m][0] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][0].ip + l_eigen_Qjm[m][1] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][1].ip + l_eigen_Qjm[m][2] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][2].ip + l_eigen_Qjm[m][3] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][3].ip + l_eigen_Qjm[m][4] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][4].ip;
			Qj_iminus_n[node[node[i].n_n[0]].n_n[0]][k] = l_eigen_Qjm[m][0] * Q[node[node[i].n_n[0]].n_n[0]][0].ip + l_eigen_Qjm[m][1] * Q[node[node[i].n_n[0]].n_n[0]][1].ip + l_eigen_Qjm[m][2] * Q[node[node[i].n_n[0]].n_n[0]][2].ip + l_eigen_Qjm[m][3] * Q[node[node[i].n_n[0]].n_n[0]][3].ip + l_eigen_Qjm[m][4] * Q[node[node[i].n_n[0]].n_n[0]][4].ip;
			Qj_iminus_n[node[i].n_n[0]][k] = l_eigen_Qjm[m][0] * Q[node[i].n_n[0]][0].ip + l_eigen_Qjm[m][1] * Q[node[i].n_n[0]][1].ip + l_eigen_Qjm[m][2] * Q[node[i].n_n[0]][2].ip + l_eigen_Qjm[m][3] * Q[node[i].n_n[0]][3].ip + l_eigen_Qjm[m][4] * Q[node[i].n_n[0]][4].ip;
			Qj_iminus_n[i][k] = l_eigen_Qjm[m][0] * Q[i][0].ip + l_eigen_Qjm[m][1] * Q[i][1].ip + l_eigen_Qjm[m][2] * Q[i][2].ip + l_eigen_Qjm[m][3] * Q[i][3].ip + l_eigen_Qjm[m][4] * Q[i][4].ip;
			Qj_iminus_n[node[i].n_n[2]][k] = l_eigen_Qjm[m][0] * Q[node[i].n_n[2]][0].ip + l_eigen_Qjm[m][1] * Q[node[i].n_n[2]][1].ip + l_eigen_Qjm[m][2] * Q[node[i].n_n[2]][2].ip + l_eigen_Qjm[m][3] * Q[node[i].n_n[2]][3].ip + l_eigen_Qjm[m][4] * Q[node[i].n_n[2]][4].ip;
			Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k] = l_eigen_Qjm[m][0] * Q[node[node[i].n_n[2]].n_n[2]][0].ip + l_eigen_Qjm[m][1] * Q[node[node[i].n_n[2]].n_n[2]][1].ip + l_eigen_Qjm[m][2] * Q[node[node[i].n_n[2]].n_n[2]][2].ip + l_eigen_Qjm[m][3] * Q[node[node[i].n_n[2]].n_n[2]][3].ip + l_eigen_Qjm[m][4] * Q[node[node[i].n_n[2]].n_n[2]][4].ip;
			Qj_iminus_n[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k] = l_eigen_Qjm[m][0] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][0].ip + l_eigen_Qjm[m][1] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][1].ip + l_eigen_Qjm[m][2] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][2].ip + l_eigen_Qjm[m][3] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][3].ip + l_eigen_Qjm[m][4] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][4].ip;
		}

		m++;
	}
}

void weno_solver_3()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(ROE_AVER[i].roe_h_im - 0.5*(ROE_AVER[i].roe_u_im * ROE_AVER[i].roe_u_im + ROE_AVER[i].roe_v_im * ROE_AVER[i].roe_v_im + ROE_AVER[i].roe_w_im * ROE_AVER[i].roe_w_im)));

	//	if(hk == 3)
	{
		zeta_xi = metric[i].zeta_xim;
		zeta_yi = metric[i].zeta_yim;
		zeta_zi = metric[i].zeta_zim;
	}

	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * ROE_AVER[i].roe_u_im + ky_bar * ROE_AVER[i].roe_v_im + kz_bar * ROE_AVER[i].roe_w_im;
	phi_sq = 0.5*0.4*(ROE_AVER[i].roe_u_im * ROE_AVER[i].roe_u_im + ROE_AVER[i].roe_v_im * ROE_AVER[i].roe_v_im + ROE_AVER[i].roe_w_im * ROE_AVER[i].roe_w_im);
	alpha = ROE_AVER[i].roe_rho_im / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[i].roe_rho_im * roe_a[i]);

	r_eigen_Qim[0][0] = kx_bar;
	r_eigen_Qim[0][1] = ky_bar;
	r_eigen_Qim[0][2] = kz_bar;
	r_eigen_Qim[0][3] = alpha;
	r_eigen_Qim[0][4] = alpha;

	r_eigen_Qim[1][0] = kx_bar * ROE_AVER[i].roe_u_im;
	r_eigen_Qim[1][1] = ky_bar * ROE_AVER[i].roe_u_im - kz_bar * ROE_AVER[i].roe_rho_im;
	r_eigen_Qim[1][2] = kz_bar * ROE_AVER[i].roe_u_im + ky_bar * ROE_AVER[i].roe_rho_im;
	r_eigen_Qim[1][3] = alpha * (ROE_AVER[i].roe_u_im + kx_bar * roe_a[i]);
	r_eigen_Qim[1][4] = alpha * (ROE_AVER[i].roe_u_im - kx_bar * roe_a[i]);

	r_eigen_Qim[2][0] = kx_bar * ROE_AVER[i].roe_v_im + kz_bar * ROE_AVER[i].roe_rho_im;
	r_eigen_Qim[2][1] = ky_bar * ROE_AVER[i].roe_v_im;
	r_eigen_Qim[2][2] = kz_bar * ROE_AVER[i].roe_v_im - kx_bar * ROE_AVER[i].roe_rho_im;
	r_eigen_Qim[2][3] = alpha * (ROE_AVER[i].roe_v_im + ky_bar * roe_a[i]);
	r_eigen_Qim[2][4] = alpha * (ROE_AVER[i].roe_v_im - ky_bar * roe_a[i]);

	r_eigen_Qim[3][0] = kx_bar * ROE_AVER[i].roe_w_im - ky_bar * ROE_AVER[i].roe_rho_im;
	r_eigen_Qim[3][1] = ky_bar * ROE_AVER[i].roe_w_im + kx_bar * ROE_AVER[i].roe_rho_im;
	r_eigen_Qim[3][2] = kz_bar * ROE_AVER[i].roe_w_im;
	r_eigen_Qim[3][3] = alpha * (ROE_AVER[i].roe_w_im + kz_bar * roe_a[i]);
	r_eigen_Qim[3][4] = alpha * (ROE_AVER[i].roe_w_im - kz_bar * roe_a[i]);

	r_eigen_Qim[4][0] = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_im * (kz_bar*ROE_AVER[i].roe_v_im - ky_bar * ROE_AVER[i].roe_w_im);
	r_eigen_Qim[4][1] = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_im * (kx_bar*ROE_AVER[i].roe_w_im - kz_bar * ROE_AVER[i].roe_u_im);
	r_eigen_Qim[4][2] = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_im * (ky_bar*ROE_AVER[i].roe_u_im - kx_bar * ROE_AVER[i].roe_v_im);
	r_eigen_Qim[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qim[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qim[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*ROE_AVER[i].roe_v_im - ky_bar * ROE_AVER[i].roe_w_im) / ROE_AVER[i].roe_rho_im);
	l_eigen_Qim[0][1] = kx_bar * 0.4*ROE_AVER[i].roe_u_im / (roe_a[i] * roe_a[i]);
	l_eigen_Qim[0][2] = (kx_bar*0.4*ROE_AVER[i].roe_v_im / (roe_a[i] * roe_a[i])) + (kz_bar / ROE_AVER[i].roe_rho_im);
	l_eigen_Qim[0][3] = (kx_bar*0.4*ROE_AVER[i].roe_w_im / (roe_a[i] * roe_a[i])) - (ky_bar / ROE_AVER[i].roe_rho_im);
	l_eigen_Qim[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qim[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*ROE_AVER[i].roe_w_im - kz_bar * ROE_AVER[i].roe_u_im) / ROE_AVER[i].roe_rho_im);
	l_eigen_Qim[1][1] = (ky_bar*0.4*ROE_AVER[i].roe_u_im / (roe_a[i] * roe_a[i])) - kz_bar / ROE_AVER[i].roe_rho_im;
	l_eigen_Qim[1][2] = ky_bar * 0.4*ROE_AVER[i].roe_v_im / (roe_a[i] * roe_a[i]);
	l_eigen_Qim[1][3] = (ky_bar*0.4*ROE_AVER[i].roe_w_im / (roe_a[i] * roe_a[i])) + (kx_bar / ROE_AVER[i].roe_rho_im);
	l_eigen_Qim[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qim[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*ROE_AVER[i].roe_u_im - kx_bar * ROE_AVER[i].roe_v_im) / ROE_AVER[i].roe_rho_im);
	l_eigen_Qim[2][1] = kz_bar * 0.4*ROE_AVER[i].roe_u_im / (roe_a[i] * roe_a[i]) + (ky_bar / ROE_AVER[i].roe_rho_im);
	l_eigen_Qim[2][2] = kz_bar * 0.4*ROE_AVER[i].roe_v_im / (roe_a[i] * roe_a[i]) - (kx_bar / ROE_AVER[i].roe_rho_im);
	l_eigen_Qim[2][3] = kz_bar * 0.4*ROE_AVER[i].roe_w_im / (roe_a[i] * roe_a[i]);
	l_eigen_Qim[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qim[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qim[3][1] = -beta * (0.4*ROE_AVER[i].roe_u_im - kx_bar * roe_a[i]);
	l_eigen_Qim[3][2] = -beta * (0.4*ROE_AVER[i].roe_v_im - ky_bar * roe_a[i]);
	l_eigen_Qim[3][3] = -beta * (0.4*ROE_AVER[i].roe_w_im - kz_bar * roe_a[i]);
	l_eigen_Qim[3][4] = beta * 0.4;

	l_eigen_Qim[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qim[4][1] = -beta * (0.4*ROE_AVER[i].roe_u_im + kx_bar * roe_a[i]);
	l_eigen_Qim[4][2] = -beta * (0.4*ROE_AVER[i].roe_v_im + ky_bar * roe_a[i]);
	l_eigen_Qim[4][3] = -beta * (0.4*ROE_AVER[i].roe_w_im + kz_bar * roe_a[i]);
	l_eigen_Qim[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		//	if(hk == 3 )
		{
			Qi_iminus_n[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k] = l_eigen_Qim[m][0] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][0].ip + l_eigen_Qim[m][1] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][1].ip + l_eigen_Qim[m][2] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][2].ip + l_eigen_Qim[m][3] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][3].ip + l_eigen_Qim[m][4] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][4].ip;
			Qi_iminus_n[node[node[i].n_n[1]].n_n[1]][k] = l_eigen_Qim[m][0] * Q[node[node[i].n_n[1]].n_n[1]][0].ip + l_eigen_Qim[m][1] * Q[node[node[i].n_n[1]].n_n[1]][1].ip + l_eigen_Qim[m][2] * Q[node[node[i].n_n[1]].n_n[1]][2].ip + l_eigen_Qim[m][3] * Q[node[node[i].n_n[1]].n_n[1]][3].ip + l_eigen_Qim[m][4] * Q[node[node[i].n_n[1]].n_n[1]][4].ip;
			Qi_iminus_n[node[i].n_n[1]][k] = l_eigen_Qim[m][0] * Q[node[i].n_n[1]][0].ip + l_eigen_Qim[m][1] * Q[node[i].n_n[1]][1].ip + l_eigen_Qim[m][2] * Q[node[i].n_n[1]][2].ip + l_eigen_Qim[m][3] * Q[node[i].n_n[1]][3].ip + l_eigen_Qim[m][4] * Q[node[i].n_n[1]][4].ip;
			Qi_iminus_n[i][k] = l_eigen_Qim[m][0] * Q[i][0].ip + l_eigen_Qim[m][1] * Q[i][1].ip + l_eigen_Qim[m][2] * Q[i][2].ip + l_eigen_Qim[m][3] * Q[i][3].ip + l_eigen_Qim[m][4] * Q[i][4].ip;
			Qi_iminus_n[node[i].n_n[3]][k] = l_eigen_Qim[m][0] * Q[node[i].n_n[3]][0].ip + l_eigen_Qim[m][1] * Q[node[i].n_n[3]][1].ip + l_eigen_Qim[m][2] * Q[node[i].n_n[3]][2].ip + l_eigen_Qim[m][3] * Q[node[i].n_n[3]][3].ip + l_eigen_Qim[m][4] * Q[node[i].n_n[3]][4].ip;
			Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k] = l_eigen_Qim[m][0] * Q[node[node[i].n_n[3]].n_n[3]][0].ip + l_eigen_Qim[m][1] * Q[node[node[i].n_n[3]].n_n[3]][1].ip + l_eigen_Qim[m][2] * Q[node[node[i].n_n[3]].n_n[3]][2].ip + l_eigen_Qim[m][3] * Q[node[node[i].n_n[3]].n_n[3]][3].ip + l_eigen_Qim[m][4] * Q[node[node[i].n_n[3]].n_n[3]][4].ip;
			Qi_iminus_n[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k] = l_eigen_Qim[m][0] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][0].ip + l_eigen_Qim[m][1] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][1].ip + l_eigen_Qim[m][2] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][2].ip + l_eigen_Qim[m][3] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][3].ip + l_eigen_Qim[m][4] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][4].ip;
		}
		m++;
	}
}

void weno_solver_4()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(ROE_AVER[i].roe_h_kp - 0.5*(ROE_AVER[i].roe_u_kp * ROE_AVER[i].roe_u_kp + ROE_AVER[i].roe_v_kp * ROE_AVER[i].roe_v_kp + ROE_AVER[i].roe_w_kp * ROE_AVER[i].roe_w_kp)));

	//	if (hk == 4)
	{
		zeta_xi = metric[i].xi_xkp;
		zeta_yi = metric[i].xi_ykp;
		zeta_zi = metric[i].xi_zkp;
	}

	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * ROE_AVER[i].roe_u_kp + ky_bar * ROE_AVER[i].roe_v_kp + kz_bar * ROE_AVER[i].roe_w_kp;
	phi_sq = 0.5*0.4*(ROE_AVER[i].roe_u_kp * ROE_AVER[i].roe_u_kp + ROE_AVER[i].roe_v_kp * ROE_AVER[i].roe_v_kp + ROE_AVER[i].roe_w_kp * ROE_AVER[i].roe_w_kp);
	alpha = ROE_AVER[i].roe_rho_kp / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[i].roe_rho_kp * roe_a[i]);

	r_eigen_Qkp[0][0] = kx_bar;
	r_eigen_Qkp[0][1] = ky_bar;
	r_eigen_Qkp[0][2] = kz_bar;
	r_eigen_Qkp[0][3] = alpha;
	r_eigen_Qkp[0][4] = alpha;

	r_eigen_Qkp[1][0] = kx_bar * ROE_AVER[i].roe_u_kp;
	r_eigen_Qkp[1][1] = ky_bar * ROE_AVER[i].roe_u_kp - kz_bar * ROE_AVER[i].roe_rho_kp;
	r_eigen_Qkp[1][2] = kz_bar * ROE_AVER[i].roe_u_kp + ky_bar * ROE_AVER[i].roe_rho_kp;
	r_eigen_Qkp[1][3] = alpha * (ROE_AVER[i].roe_u_kp + kx_bar * roe_a[i]);
	r_eigen_Qkp[1][4] = alpha * (ROE_AVER[i].roe_u_kp - kx_bar * roe_a[i]);

	r_eigen_Qkp[2][0] = kx_bar * ROE_AVER[i].roe_v_kp + kz_bar * ROE_AVER[i].roe_rho_kp;
	r_eigen_Qkp[2][1] = ky_bar * ROE_AVER[i].roe_v_kp;
	r_eigen_Qkp[2][2] = kz_bar * ROE_AVER[i].roe_v_kp - kx_bar * ROE_AVER[i].roe_rho_kp;
	r_eigen_Qkp[2][3] = alpha * (ROE_AVER[i].roe_v_kp + ky_bar * roe_a[i]);
	r_eigen_Qkp[2][4] = alpha * (ROE_AVER[i].roe_v_kp - ky_bar * roe_a[i]);

	r_eigen_Qkp[3][0] = kx_bar * ROE_AVER[i].roe_w_kp - ky_bar * ROE_AVER[i].roe_rho_kp;
	r_eigen_Qkp[3][1] = ky_bar * ROE_AVER[i].roe_w_kp + kx_bar * ROE_AVER[i].roe_rho_kp;
	r_eigen_Qkp[3][2] = kz_bar * ROE_AVER[i].roe_w_kp;
	r_eigen_Qkp[3][3] = alpha * (ROE_AVER[i].roe_w_kp + kz_bar * roe_a[i]);
	r_eigen_Qkp[3][4] = alpha * (ROE_AVER[i].roe_w_kp - kz_bar * roe_a[i]);

	r_eigen_Qkp[4][0] = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_kp * (kz_bar*ROE_AVER[i].roe_v_kp - ky_bar * ROE_AVER[i].roe_w_kp);
	r_eigen_Qkp[4][1] = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_kp * (kx_bar*ROE_AVER[i].roe_w_kp - kz_bar * ROE_AVER[i].roe_u_kp);
	r_eigen_Qkp[4][2] = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_kp * (ky_bar*ROE_AVER[i].roe_u_kp - kx_bar * ROE_AVER[i].roe_v_kp);
	r_eigen_Qkp[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qkp[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qkp[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*ROE_AVER[i].roe_v_kp - ky_bar * ROE_AVER[i].roe_w_kp) / ROE_AVER[i].roe_rho_kp);
	l_eigen_Qkp[0][1] = kx_bar * 0.4*ROE_AVER[i].roe_u_kp / (roe_a[i] * roe_a[i]);
	l_eigen_Qkp[0][2] = (kx_bar*0.4*ROE_AVER[i].roe_v_kp / (roe_a[i] * roe_a[i])) + (kz_bar / ROE_AVER[i].roe_rho_kp);
	l_eigen_Qkp[0][3] = (kx_bar*0.4*ROE_AVER[i].roe_w_kp / (roe_a[i] * roe_a[i])) - (ky_bar / ROE_AVER[i].roe_rho_kp);
	l_eigen_Qkp[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkp[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*ROE_AVER[i].roe_w_kp - kz_bar * ROE_AVER[i].roe_u_kp) / ROE_AVER[i].roe_rho_kp);
	l_eigen_Qkp[1][1] = (ky_bar*0.4*ROE_AVER[i].roe_u_kp / (roe_a[i] * roe_a[i])) - kz_bar / ROE_AVER[i].roe_rho_kp;
	l_eigen_Qkp[1][2] = ky_bar * 0.4*ROE_AVER[i].roe_v_kp / (roe_a[i] * roe_a[i]);
	l_eigen_Qkp[1][3] = (ky_bar*0.4*ROE_AVER[i].roe_w_kp / (roe_a[i] * roe_a[i])) + (kx_bar / ROE_AVER[i].roe_rho_kp);
	l_eigen_Qkp[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkp[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*ROE_AVER[i].roe_u_kp - kx_bar * ROE_AVER[i].roe_v_kp) / ROE_AVER[i].roe_rho_kp);
	l_eigen_Qkp[2][1] = kz_bar * 0.4*ROE_AVER[i].roe_u_kp / (roe_a[i] * roe_a[i]) + (ky_bar / ROE_AVER[i].roe_rho_kp);
	l_eigen_Qkp[2][2] = kz_bar * 0.4*ROE_AVER[i].roe_v_kp / (roe_a[i] * roe_a[i]) - (kx_bar / ROE_AVER[i].roe_rho_kp);
	l_eigen_Qkp[2][3] = kz_bar * 0.4*ROE_AVER[i].roe_w_kp / (roe_a[i] * roe_a[i]);
	l_eigen_Qkp[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkp[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qkp[3][1] = -beta * (0.4*ROE_AVER[i].roe_u_kp - kx_bar * roe_a[i]);
	l_eigen_Qkp[3][2] = -beta * (0.4*ROE_AVER[i].roe_v_kp - ky_bar * roe_a[i]);
	l_eigen_Qkp[3][3] = -beta * (0.4*ROE_AVER[i].roe_w_kp - kz_bar * roe_a[i]);
	l_eigen_Qkp[3][4] = beta * 0.4;

	l_eigen_Qkp[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qkp[4][1] = -beta * (0.4*ROE_AVER[i].roe_u_kp + kx_bar * roe_a[i]);
	l_eigen_Qkp[4][2] = -beta * (0.4*ROE_AVER[i].roe_v_kp + ky_bar * roe_a[i]);
	l_eigen_Qkp[4][3] = -beta * (0.4*ROE_AVER[i].roe_w_kp + kz_bar * roe_a[i]);
	l_eigen_Qkp[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		//	if (hk == 4)
		{
			Qk_iplus[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k] = l_eigen_Qkp[m][0] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][0].ip + l_eigen_Qkp[m][1] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][1].ip + l_eigen_Qkp[m][2] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][2].ip + l_eigen_Qkp[m][3] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][3].ip + l_eigen_Qkp[m][4] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][4].ip;
			Qk_iplus[node[node[i].n_n[4]].n_n[4]][k] = l_eigen_Qkp[m][0] * Q[node[node[i].n_n[4]].n_n[4]][0].ip + l_eigen_Qkp[m][1] * Q[node[node[i].n_n[4]].n_n[4]][1].ip + l_eigen_Qkp[m][2] * Q[node[node[i].n_n[4]].n_n[4]][2].ip + l_eigen_Qkp[m][3] * Q[node[node[i].n_n[4]].n_n[4]][3].ip + l_eigen_Qkp[m][4] * Q[node[node[i].n_n[4]].n_n[4]][4].ip;
			Qk_iplus[node[i].n_n[4]][k] = l_eigen_Qkp[m][0] * Q[node[i].n_n[4]][0].ip + l_eigen_Qkp[m][1] * Q[node[i].n_n[4]][1].ip + l_eigen_Qkp[m][2] * Q[node[i].n_n[4]][2].ip + l_eigen_Qkp[m][3] * Q[node[i].n_n[4]][3].ip + l_eigen_Qkp[m][4] * Q[node[i].n_n[4]][4].ip;
			Qk_iplus[i][k] = l_eigen_Qkp[m][0] * Q[i][0].ip + l_eigen_Qkp[m][1] * Q[i][1].ip + l_eigen_Qkp[m][2] * Q[i][2].ip + l_eigen_Qkp[m][3] * Q[i][3].ip + l_eigen_Qkp[m][4] * Q[i][4].ip;
			Qk_iplus[node[i].n_n[5]][k] = l_eigen_Qkp[m][0] * Q[node[i].n_n[5]][0].ip + l_eigen_Qkp[m][1] * Q[node[i].n_n[5]][1].ip + l_eigen_Qkp[m][2] * Q[node[i].n_n[5]][2].ip + l_eigen_Qkp[m][3] * Q[node[i].n_n[5]][3].ip + l_eigen_Qkp[m][4] * Q[node[i].n_n[5]][4].ip;
			Qk_iplus[node[node[i].n_n[5]].n_n[5]][k] = l_eigen_Qkp[m][0] * Q[node[node[i].n_n[5]].n_n[5]][0].ip + l_eigen_Qkp[m][1] * Q[node[node[i].n_n[5]].n_n[5]][1].ip + l_eigen_Qkp[m][2] * Q[node[node[i].n_n[5]].n_n[5]][2].ip + l_eigen_Qkp[m][3] * Q[node[node[i].n_n[5]].n_n[5]][3].ip + l_eigen_Qkp[m][4] * Q[node[node[i].n_n[5]].n_n[5]][4].ip;
			Qk_iplus[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k] = l_eigen_Qkp[m][0] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][0].ip + l_eigen_Qkp[m][1] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][1].ip + l_eigen_Qkp[m][2] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][2].ip + l_eigen_Qkp[m][3] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][3].ip + l_eigen_Qkp[m][4] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][4].ip;
		}

		m++;
	}
}


void weno_solver_5()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(ROE_AVER[i].roe_h_km - 0.5*(ROE_AVER[i].roe_u_km * ROE_AVER[i].roe_u_km + ROE_AVER[i].roe_v_km * ROE_AVER[i].roe_v_km + ROE_AVER[i].roe_w_km * ROE_AVER[i].roe_w_km)));

	//	if(hk == 5)
	{
		zeta_xi = metric[i].xi_xkm;
		zeta_yi = metric[i].xi_ykm;
		zeta_zi = metric[i].xi_zkm;
	}

	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * ROE_AVER[i].roe_u_km + ky_bar * ROE_AVER[i].roe_v_km + kz_bar * ROE_AVER[i].roe_w_km;
	phi_sq = 0.5*0.4*(ROE_AVER[i].roe_u_km * ROE_AVER[i].roe_u_km + ROE_AVER[i].roe_v_km * ROE_AVER[i].roe_v_km + ROE_AVER[i].roe_w_km * ROE_AVER[i].roe_w_km);
	alpha = ROE_AVER[i].roe_rho_km / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[i].roe_rho_km * roe_a[i]);

	r_eigen_Qkm[0][0] = kx_bar;
	r_eigen_Qkm[0][1] = ky_bar;
	r_eigen_Qkm[0][2] = kz_bar;
	r_eigen_Qkm[0][3] = alpha;
	r_eigen_Qkm[0][4] = alpha;

	r_eigen_Qkm[1][0] = kx_bar * ROE_AVER[i].roe_u_km;
	r_eigen_Qkm[1][1] = ky_bar * ROE_AVER[i].roe_u_km - kz_bar * ROE_AVER[i].roe_rho_km;
	r_eigen_Qkm[1][2] = kz_bar * ROE_AVER[i].roe_u_km + ky_bar * ROE_AVER[i].roe_rho_km;
	r_eigen_Qkm[1][3] = alpha * (ROE_AVER[i].roe_u_km + kx_bar * roe_a[i]);
	r_eigen_Qkm[1][4] = alpha * (ROE_AVER[i].roe_u_km - kx_bar * roe_a[i]);

	r_eigen_Qkm[2][0] = kx_bar * ROE_AVER[i].roe_v_km + kz_bar * ROE_AVER[i].roe_rho_km;
	r_eigen_Qkm[2][1] = ky_bar * ROE_AVER[i].roe_v_km;
	r_eigen_Qkm[2][2] = kz_bar * ROE_AVER[i].roe_v_km - kx_bar * ROE_AVER[i].roe_rho_km;
	r_eigen_Qkm[2][3] = alpha * (ROE_AVER[i].roe_v_km + ky_bar * roe_a[i]);
	r_eigen_Qkm[2][4] = alpha * (ROE_AVER[i].roe_v_km - ky_bar * roe_a[i]);

	r_eigen_Qkm[3][0] = kx_bar * ROE_AVER[i].roe_w_km - ky_bar * ROE_AVER[i].roe_rho_km;
	r_eigen_Qkm[3][1] = ky_bar * ROE_AVER[i].roe_w_km + kx_bar * ROE_AVER[i].roe_rho_km;
	r_eigen_Qkm[3][2] = kz_bar * ROE_AVER[i].roe_w_km;
	r_eigen_Qkm[3][3] = alpha * (ROE_AVER[i].roe_w_km + kz_bar * roe_a[i]);
	r_eigen_Qkm[3][4] = alpha * (ROE_AVER[i].roe_w_km - kz_bar * roe_a[i]);

	r_eigen_Qkm[4][0] = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_km * (kz_bar*ROE_AVER[i].roe_v_km - ky_bar * ROE_AVER[i].roe_w_km);
	r_eigen_Qkm[4][1] = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_km * (kx_bar*ROE_AVER[i].roe_w_km - kz_bar * ROE_AVER[i].roe_u_km);
	r_eigen_Qkm[4][2] = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[i].roe_rho_km * (ky_bar*ROE_AVER[i].roe_u_km - kx_bar * ROE_AVER[i].roe_v_km);
	r_eigen_Qkm[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qkm[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qkm[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*ROE_AVER[i].roe_v_km - ky_bar * ROE_AVER[i].roe_w_km) / ROE_AVER[i].roe_rho_km);
	l_eigen_Qkm[0][1] = kx_bar * 0.4*ROE_AVER[i].roe_u_km / (roe_a[i] * roe_a[i]);
	l_eigen_Qkm[0][2] = (kx_bar*0.4*ROE_AVER[i].roe_v_km / (roe_a[i] * roe_a[i])) + (kz_bar / ROE_AVER[i].roe_rho_km);
	l_eigen_Qkm[0][3] = (kx_bar*0.4*ROE_AVER[i].roe_w_km / (roe_a[i] * roe_a[i])) - (ky_bar / ROE_AVER[i].roe_rho_km);
	l_eigen_Qkm[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkm[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*ROE_AVER[i].roe_w_km - kz_bar * ROE_AVER[i].roe_u_km) / ROE_AVER[i].roe_rho_km);
	l_eigen_Qkm[1][1] = (ky_bar*0.4*ROE_AVER[i].roe_u_km / (roe_a[i] * roe_a[i])) - kz_bar / ROE_AVER[i].roe_rho_km;
	l_eigen_Qkm[1][2] = ky_bar * 0.4*ROE_AVER[i].roe_v_km / (roe_a[i] * roe_a[i]);
	l_eigen_Qkm[1][3] = (ky_bar*0.4*ROE_AVER[i].roe_w_km / (roe_a[i] * roe_a[i])) + (kx_bar / ROE_AVER[i].roe_rho_km);
	l_eigen_Qkm[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkm[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*ROE_AVER[i].roe_u_km - kx_bar * ROE_AVER[i].roe_v_km) / ROE_AVER[i].roe_rho_km);
	l_eigen_Qkm[2][1] = kz_bar * 0.4*ROE_AVER[i].roe_u_km / (roe_a[i] * roe_a[i]) + (ky_bar / ROE_AVER[i].roe_rho_km);
	l_eigen_Qkm[2][2] = kz_bar * 0.4*ROE_AVER[i].roe_v_km / (roe_a[i] * roe_a[i]) - (kx_bar / ROE_AVER[i].roe_rho_km);
	l_eigen_Qkm[2][3] = kz_bar * 0.4*ROE_AVER[i].roe_w_km / (roe_a[i] * roe_a[i]);
	l_eigen_Qkm[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkm[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qkm[3][1] = -beta * (0.4*ROE_AVER[i].roe_u_km - kx_bar * roe_a[i]);
	l_eigen_Qkm[3][2] = -beta * (0.4*ROE_AVER[i].roe_v_km - ky_bar * roe_a[i]);
	l_eigen_Qkm[3][3] = -beta * (0.4*ROE_AVER[i].roe_w_km - kz_bar * roe_a[i]);
	l_eigen_Qkm[3][4] = beta * 0.4;

	l_eigen_Qkm[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qkm[4][1] = -beta * (0.4*ROE_AVER[i].roe_u_km + kx_bar * roe_a[i]);
	l_eigen_Qkm[4][2] = -beta * (0.4*ROE_AVER[i].roe_v_km + ky_bar * roe_a[i]);
	l_eigen_Qkm[4][3] = -beta * (0.4*ROE_AVER[i].roe_w_km + kz_bar * roe_a[i]);
	l_eigen_Qkm[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		//	if (hk == 5)
		{
			Qk_iminus_n[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k] = l_eigen_Qkm[m][0] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][0].ip + l_eigen_Qkm[m][1] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][1].ip + l_eigen_Qkm[m][2] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][2].ip + l_eigen_Qkm[m][3] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][3].ip + l_eigen_Qkm[m][4] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][4].ip;
			Qk_iminus_n[node[node[i].n_n[4]].n_n[4]][k] = l_eigen_Qkm[m][0] * Q[node[node[i].n_n[4]].n_n[4]][0].ip + l_eigen_Qkm[m][1] * Q[node[node[i].n_n[4]].n_n[4]][1].ip + l_eigen_Qkm[m][2] * Q[node[node[i].n_n[4]].n_n[4]][2].ip + l_eigen_Qkm[m][3] * Q[node[node[i].n_n[4]].n_n[4]][3].ip + l_eigen_Qkm[m][4] * Q[node[node[i].n_n[4]].n_n[4]][4].ip;
			Qk_iminus_n[node[i].n_n[4]][k] = l_eigen_Qkm[m][0] * Q[node[i].n_n[4]][0].ip + l_eigen_Qkm[m][1] * Q[node[i].n_n[4]][1].ip + l_eigen_Qkm[m][2] * Q[node[i].n_n[4]][2].ip + l_eigen_Qkm[m][3] * Q[node[i].n_n[4]][3].ip + l_eigen_Qkm[m][4] * Q[node[i].n_n[4]][4].ip;
			Qk_iminus_n[i][k] = l_eigen_Qkm[m][0] * Q[i][0].ip + l_eigen_Qkm[m][1] * Q[i][1].ip + l_eigen_Qkm[m][2] * Q[i][2].ip + l_eigen_Qkm[m][3] * Q[i][3].ip + l_eigen_Qkm[m][4] * Q[i][4].ip;
			Qk_iminus_n[node[i].n_n[5]][k] = l_eigen_Qkm[m][0] * Q[node[i].n_n[5]][0].ip + l_eigen_Qkm[m][1] * Q[node[i].n_n[5]][1].ip + l_eigen_Qkm[m][2] * Q[node[i].n_n[5]][2].ip + l_eigen_Qkm[m][3] * Q[node[i].n_n[5]][3].ip + l_eigen_Qkm[m][4] * Q[node[i].n_n[5]][4].ip;
			Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k] = l_eigen_Qkm[m][0] * Q[node[node[i].n_n[5]].n_n[5]][0].ip + l_eigen_Qkm[m][1] * Q[node[node[i].n_n[5]].n_n[5]][1].ip + l_eigen_Qkm[m][2] * Q[node[node[i].n_n[5]].n_n[5]][2].ip + l_eigen_Qkm[m][3] * Q[node[node[i].n_n[5]].n_n[5]][3].ip + l_eigen_Qkm[m][4] * Q[node[node[i].n_n[5]].n_n[5]][4].ip;
			Qk_iminus_n[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k] = l_eigen_Qkm[m][0] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][0].ip + l_eigen_Qkm[m][1] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][1].ip + l_eigen_Qkm[m][2] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][2].ip + l_eigen_Qkm[m][3] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][3].ip + l_eigen_Qkm[m][4] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][4].ip;
		}

		m++;
	}
}

