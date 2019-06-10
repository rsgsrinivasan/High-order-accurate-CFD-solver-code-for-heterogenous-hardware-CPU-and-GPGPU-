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

void intialise(int RK)
{
	int l, k, d_temp, d_temp6;
	double Cd, Cdr1, Cdr2, Cdr3, Cdc1, Cdc2, Cdc3, l_determinant, l_angle, l_angle2, Vt, Vn, Vt2, Vn2, Vn_sign;
	Cd = 0.07;
	Cdr1 = 0.0123; 
	Cdr2 = 0.0279;
	Cdr3 = 0.0264;
	Cdc1 = 0.0121;
	Cdc2 = 0.0138;
	Cdc3 = 0.0;
	
	if (initial == 0)
	{
		for (i=1; i<=sd_node; i++)
		{				
			FLOW[i].u = 1.0;
			FLOW[i].v = 0.0;
			FLOW[i].w = 0.0;
			FLOW[i].p = 1.0/(1.4*Mach*Mach);
			FLOW[i].t = 1.0;
			FLOW[i].rho = (1.4*Mach*Mach)*(FLOW[i].p/FLOW[i].t);
			FLOW[i].e = FLOW[i].p/(0.4*FLOW[i].rho);
			FLOW[i].a = sqrt(1.4*FLOW[i].p/FLOW[i].rho);																
		}		
	}
	
/****************************************************OUTLET********************************************************/		
	for (i=0; i<out_node; i++)
	{	
		if (node[outlet_node[i]].loc > 0 && node[outlet_node[i]].loc <= 6)
		{	
			if (node[outlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[outlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[outlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[outlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[outlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[outlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			FLOW[outlet_node[i]].u= FLOW[node[outlet_node[i]].n_n[k]].u;	
			FLOW[outlet_node[i]].v= FLOW[node[outlet_node[i]].n_n[k]].v;	
			FLOW[outlet_node[i]].w= FLOW[node[outlet_node[i]].n_n[k]].w;	
			FLOW[outlet_node[i]].p= FLOW[node[outlet_node[i]].n_n[k]].p;
		/* 	if (back_pressure != 0.0)
			{
				if (node[outlet_node[i]].y < 412.0 || node[outlet_node[i]].y > 767.0)
				{
					FLOW[outlet_node[i]].p = FLOW[node[outlet_node[i]].n_n[k]].p;	
				}
				if (node[outlet_node[i]].y > 412.0 && node[outlet_node[i]].y < 767.0)
				{
					FLOW[outlet_node[i]].p = back_pressure*(1.0/(1.4*Mach*Mach));	
				}
			}	 */			
			FLOW[outlet_node[i]].t= FLOW[node[outlet_node[i]].n_n[k]].t;	
			FLOW[outlet_node[i]].rho= (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
			FLOW[outlet_node[i]].a = sqrt(1.4*FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].rho); 

			FLOW[node[outlet_node[i]].n_n[l]].u = FLOW[outlet_node[i]].u;
			FLOW[node[outlet_node[i]].n_n[l]].v = FLOW[outlet_node[i]].v;
			FLOW[node[outlet_node[i]].n_n[l]].w = FLOW[outlet_node[i]].w;
			FLOW[node[outlet_node[i]].n_n[l]].p = FLOW[outlet_node[i]].p;
			FLOW[node[outlet_node[i]].n_n[l]].rho = FLOW[outlet_node[i]].rho;
			FLOW[node[outlet_node[i]].n_n[l]].t = FLOW[outlet_node[i]].t;	
			FLOW[node[outlet_node[i]].n_n[l]].a = FLOW[outlet_node[i]].a;
			FLOW[node[outlet_node[i]].n_n[l]].e = FLOW[outlet_node[i]].e;
			
			FLOW[node[node[outlet_node[i]].n_n[l]].n_n[l]].u = FLOW[outlet_node[i]].u;
			FLOW[node[node[outlet_node[i]].n_n[l]].n_n[l]].v = FLOW[outlet_node[i]].v;
			FLOW[node[node[outlet_node[i]].n_n[l]].n_n[l]].w = FLOW[outlet_node[i]].w;
			FLOW[node[node[outlet_node[i]].n_n[l]].n_n[l]].p = FLOW[outlet_node[i]].p;
			FLOW[node[node[outlet_node[i]].n_n[l]].n_n[l]].rho = FLOW[outlet_node[i]].rho;
			FLOW[node[node[outlet_node[i]].n_n[l]].n_n[l]].t = FLOW[outlet_node[i]].t;	
			FLOW[node[node[outlet_node[i]].n_n[l]].n_n[l]].a = FLOW[outlet_node[i]].a;
			FLOW[node[node[outlet_node[i]].n_n[l]].n_n[l]].e = FLOW[outlet_node[i]].e;
	
			FLOW[node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].u = FLOW[outlet_node[i]].u;
			FLOW[node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].v = FLOW[outlet_node[i]].v;
			FLOW[node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].w = FLOW[outlet_node[i]].w;
			FLOW[node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[outlet_node[i]].p;
			FLOW[node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[outlet_node[i]].rho;
			FLOW[node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[outlet_node[i]].t;
			FLOW[node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[outlet_node[i]].a;
			FLOW[node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[outlet_node[i]].e;
			
		}
		else if (node[outlet_node[i]].loc == 7 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);

		}
		else if (node[outlet_node[i]].loc == 8 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);			
		}
		else if (node[outlet_node[i]].loc == 9 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);			
		}
		else if (node[outlet_node[i]].loc == 10 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		
		else if (node[outlet_node[i]].loc == 11 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 12 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);			
		}
		else if (node[outlet_node[i]].loc == 13 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);		
		}
		else if (node[outlet_node[i]].loc == 14 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);	
		}
		else if (node[outlet_node[i]].loc == 15 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[5]].n_n[0]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[5]].n_n[0]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[5]].n_n[0]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[5]].n_n[0]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[5]].n_n[0]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 16 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[0]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[0]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[0]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[0]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[0]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 17 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[4]].n_n[0]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[4]].n_n[0]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[4]].n_n[0]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[4]].n_n[0]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[4]].n_n[0]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 18 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[0]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[0]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[0]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[0]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[0]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 19 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[5]].n_n[2]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[5]].n_n[2]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[5]].n_n[2]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[5]].n_n[2]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[5]].n_n[2]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 20 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[2]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[2]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[2]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[2]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[2]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 21 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[4]].n_n[2]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[4]].n_n[2]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[4]].n_n[2]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[4]].n_n[2]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[4]].n_n[2]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 22 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[2]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[2]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[2]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[2]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[2]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 23 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[5]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[5]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[5]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[5]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[5]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 24 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[4]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[4]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[4]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[4]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[1]].n_n[4]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 25 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[4]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[4]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[4]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[4]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[4]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
		else if (node[outlet_node[i]].loc == 26 )
		{
			FLOW[outlet_node[i]].u = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[5]].u;
			FLOW[outlet_node[i]].v = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[5]].v;
			FLOW[outlet_node[i]].w = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[5]].w;
			FLOW[outlet_node[i]].p = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[5]].p;
			FLOW[outlet_node[i]].t = FLOW[node[node[outlet_node[i]].n_n[3]].n_n[5]].t;
			FLOW[outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[outlet_node[i]].p/FLOW[outlet_node[i]].t);
			FLOW[outlet_node[i]].e = FLOW[outlet_node[i]].p/(0.4*FLOW[outlet_node[i]].rho);
		}
	}

/****************************************************INLET********************************************************/
	for (i=0; i< inl_node; i++)
	{	
		if (node[inlet_node[i]].loc > 0 && node[inlet_node[i]].loc <= 6)
		{	
			if (node[inlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[inlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[inlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[inlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[inlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[inlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}
		
			FLOW[inlet_node[i]].u = 1.0;		
			FLOW[inlet_node[i]].v = 0.0;	
			FLOW[inlet_node[i]].w = 0.0;				
			FLOW[inlet_node[i]].p = 1.0/(1.4*Mach*Mach);		
			FLOW[inlet_node[i]].t = 1.0;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
			FLOW[inlet_node[i]].a = sqrt(1.4*FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].rho); 

			FLOW[node[inlet_node[i]].n_n[l]].u = FLOW[inlet_node[i]].u;
			FLOW[node[inlet_node[i]].n_n[l]].v = FLOW[inlet_node[i]].v;
			FLOW[node[inlet_node[i]].n_n[l]].w = FLOW[inlet_node[i]].w;
			FLOW[node[inlet_node[i]].n_n[l]].p = FLOW[inlet_node[i]].p;
			FLOW[node[inlet_node[i]].n_n[l]].rho = FLOW[inlet_node[i]].rho;
			FLOW[node[inlet_node[i]].n_n[l]].t = FLOW[inlet_node[i]].t;	
			FLOW[node[inlet_node[i]].n_n[l]].a = FLOW[inlet_node[i]].a;
			FLOW[node[inlet_node[i]].n_n[l]].e = FLOW[inlet_node[i]].e;
			
			FLOW[node[node[inlet_node[i]].n_n[l]].n_n[l]].u = FLOW[inlet_node[i]].u;
			FLOW[node[node[inlet_node[i]].n_n[l]].n_n[l]].v = FLOW[inlet_node[i]].v;
			FLOW[node[node[inlet_node[i]].n_n[l]].n_n[l]].w = FLOW[inlet_node[i]].w;
			FLOW[node[node[inlet_node[i]].n_n[l]].n_n[l]].p = FLOW[inlet_node[i]].p;
			FLOW[node[node[inlet_node[i]].n_n[l]].n_n[l]].rho = FLOW[inlet_node[i]].rho;
			FLOW[node[node[inlet_node[i]].n_n[l]].n_n[l]].t = FLOW[inlet_node[i]].t;	
			FLOW[node[node[inlet_node[i]].n_n[l]].n_n[l]].a = FLOW[inlet_node[i]].a;
			FLOW[node[node[inlet_node[i]].n_n[l]].n_n[l]].e = FLOW[inlet_node[i]].e;
	
			FLOW[node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].u = FLOW[inlet_node[i]].u;
			FLOW[node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].v = FLOW[inlet_node[i]].v;
			FLOW[node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].w = FLOW[inlet_node[i]].w;
			FLOW[node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[inlet_node[i]].p;
			FLOW[node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[inlet_node[i]].rho;
			FLOW[node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[inlet_node[i]].t;
			FLOW[node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[inlet_node[i]].a;
			FLOW[node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[inlet_node[i]].e;
		
		}
		else if (node[inlet_node[i]].loc == 7 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);

		}
		else if (node[inlet_node[i]].loc == 8 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);			
		}
		else if (node[inlet_node[i]].loc == 9 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);			
		}
		else if (node[inlet_node[i]].loc == 10 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		
		else if (node[inlet_node[i]].loc == 11 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 12 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);			
		}
		else if (node[inlet_node[i]].loc == 13 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);		
		}
		else if (node[inlet_node[i]].loc == 14 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);	
		}
		else if (node[inlet_node[i]].loc == 15 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[5]].n_n[0]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[5]].n_n[0]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[5]].n_n[0]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[5]].n_n[0]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[5]].n_n[0]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 16 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[0]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[0]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[0]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[0]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[0]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 17 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[4]].n_n[0]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[4]].n_n[0]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[4]].n_n[0]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[4]].n_n[0]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[4]].n_n[0]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 18 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[0]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[0]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[0]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[0]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[0]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 19 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[5]].n_n[2]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[5]].n_n[2]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[5]].n_n[2]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[5]].n_n[2]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[5]].n_n[2]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 20 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[2]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[2]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[2]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[2]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[2]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 21 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[4]].n_n[2]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[4]].n_n[2]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[4]].n_n[2]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[4]].n_n[2]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[4]].n_n[2]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 22 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[2]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[2]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[2]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[2]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[2]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 23 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[5]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[5]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[5]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[5]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[5]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 24 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[4]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[4]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[4]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[4]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[1]].n_n[4]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 25 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[4]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[4]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[4]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[4]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[4]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
		else if (node[inlet_node[i]].loc == 26 )
		{
			FLOW[inlet_node[i]].u = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[5]].u;
			FLOW[inlet_node[i]].v = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[5]].v;
			FLOW[inlet_node[i]].w = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[5]].w;
			FLOW[inlet_node[i]].p = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[5]].p;
			FLOW[inlet_node[i]].t = FLOW[node[node[inlet_node[i]].n_n[3]].n_n[5]].t;
			FLOW[inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[inlet_node[i]].p/FLOW[inlet_node[i]].t);
			FLOW[inlet_node[i]].e = FLOW[inlet_node[i]].p/(0.4*FLOW[inlet_node[i]].rho);
		}
	}
	
	/****************************************************BOUNDARY********************************************************/
	for (i=0; i< bou_node; i++)
	{	
		if (node[boundary_node[i]].loc > 0 && node[boundary_node[i]].loc <= 6)
		{	
			if (node[boundary_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[boundary_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[boundary_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[boundary_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[boundary_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[boundary_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			FLOW[boundary_node[i]].u = FLOW[node[boundary_node[i]].n_n[k]].u;
			FLOW[boundary_node[i]].v = FLOW[node[boundary_node[i]].n_n[k]].v;
			FLOW[boundary_node[i]].w = FLOW[node[boundary_node[i]].n_n[k]].w;
			FLOW[boundary_node[i]].p = FLOW[node[boundary_node[i]].n_n[k]].p;
			FLOW[boundary_node[i]].rho = FLOW[node[boundary_node[i]].n_n[k]].rho;
			FLOW[boundary_node[i]].t = FLOW[node[boundary_node[i]].n_n[k]].t;	
			FLOW[boundary_node[i]].a = FLOW[node[boundary_node[i]].n_n[k]].a;
			FLOW[boundary_node[i]].e = FLOW[node[boundary_node[i]].n_n[k]].e;
			//mu[RK][boundary_node[i]] = mu[RK][boundary_node[i]];
			
			FLOW[node[boundary_node[i]].n_n[l]].u = FLOW[node[boundary_node[i]].n_n[k]].u;
			FLOW[node[boundary_node[i]].n_n[l]].v = FLOW[node[boundary_node[i]].n_n[k]].v;
			FLOW[node[boundary_node[i]].n_n[l]].w = FLOW[node[boundary_node[i]].n_n[k]].w;
			FLOW[node[boundary_node[i]].n_n[l]].p = FLOW[node[boundary_node[i]].n_n[k]].p;
			FLOW[node[boundary_node[i]].n_n[l]].rho = FLOW[node[boundary_node[i]].n_n[k]].rho;
			FLOW[node[boundary_node[i]].n_n[l]].t = FLOW[node[boundary_node[i]].n_n[k]].t;	
			FLOW[node[boundary_node[i]].n_n[l]].a = FLOW[node[boundary_node[i]].n_n[k]].a;
			FLOW[node[boundary_node[i]].n_n[l]].e = FLOW[node[boundary_node[i]].n_n[k]].e;
			//mFLOW[node[boundary_node[i]].n_n[l]] = mFLOW[boundary_node[i]];
			
			FLOW[node[node[boundary_node[i]].n_n[l]].n_n[l]].u = FLOW[node[node[boundary_node[i]].n_n[k]].n_n[k]].u;
			FLOW[node[node[boundary_node[i]].n_n[l]].n_n[l]].v = FLOW[node[node[boundary_node[i]].n_n[k]].n_n[k]].v;
			FLOW[node[node[boundary_node[i]].n_n[l]].n_n[l]].w = FLOW[node[node[boundary_node[i]].n_n[k]].n_n[k]].w;
			FLOW[node[node[boundary_node[i]].n_n[l]].n_n[l]].p = FLOW[node[node[boundary_node[i]].n_n[k]].n_n[k]].p;
			FLOW[node[node[boundary_node[i]].n_n[l]].n_n[l]].rho = FLOW[node[node[boundary_node[i]].n_n[k]].n_n[k]].rho;
			FLOW[node[node[boundary_node[i]].n_n[l]].n_n[l]].t = FLOW[node[node[boundary_node[i]].n_n[k]].n_n[k]].t;	
			FLOW[node[node[boundary_node[i]].n_n[l]].n_n[l]].a = FLOW[node[node[boundary_node[i]].n_n[k]].n_n[k]].a;
			FLOW[node[node[boundary_node[i]].n_n[l]].n_n[l]].e = FLOW[node[node[boundary_node[i]].n_n[k]].n_n[k]].e;
			//mFLOW[node[node[boundary_node[i]].n_n[l]].n_n[l]] = mFLOW[boundary_node[i]];
	
			FLOW[node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].u = FLOW[node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].u;
			FLOW[node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].v = FLOW[node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].v;
			FLOW[node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].w = FLOW[node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].w;
			FLOW[node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].p;
			FLOW[node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].rho;
			FLOW[node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].t;	
			FLOW[node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].a;
			FLOW[node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].e;
			//mFLOW[node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = mFLOW[boundary_node[i]];
		}
		else if (node[boundary_node[i]].loc == 7 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);

		}
		else if (node[boundary_node[i]].loc == 8 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);			
		}
		else if (node[boundary_node[i]].loc == 9 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);			
		}
		else if (node[boundary_node[i]].loc == 10 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		
		else if (node[boundary_node[i]].loc == 11 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 12 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);			
		}
		else if (node[boundary_node[i]].loc == 13 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);		
		}
		else if (node[boundary_node[i]].loc == 14 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);	
		}
		else if (node[boundary_node[i]].loc == 15 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[5]].n_n[0]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[5]].n_n[0]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[5]].n_n[0]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[5]].n_n[0]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[5]].n_n[0]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 16 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[0]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[0]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[0]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[0]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[0]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 17 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[4]].n_n[0]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[4]].n_n[0]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[4]].n_n[0]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[4]].n_n[0]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[4]].n_n[0]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 18 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[0]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[0]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[0]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[0]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[0]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 19 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[5]].n_n[2]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[5]].n_n[2]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[5]].n_n[2]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[5]].n_n[2]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[5]].n_n[2]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 20 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[2]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[2]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[2]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[2]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[2]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 21 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[4]].n_n[2]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[4]].n_n[2]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[4]].n_n[2]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[4]].n_n[2]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[4]].n_n[2]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 22 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[2]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[2]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[2]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[2]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[2]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 23 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[5]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[5]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[5]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[5]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[5]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 24 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[4]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[4]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[4]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[4]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[1]].n_n[4]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 25 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[4]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[4]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[4]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[4]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[4]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}
		else if (node[boundary_node[i]].loc == 26 )
		{
			FLOW[boundary_node[i]].u = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[5]].u;
			FLOW[boundary_node[i]].v = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[5]].v;
			FLOW[boundary_node[i]].w = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[5]].w;
			FLOW[boundary_node[i]].p = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[5]].p;
			FLOW[boundary_node[i]].t = FLOW[node[node[boundary_node[i]].n_n[3]].n_n[5]].t;
			FLOW[boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[boundary_node[i]].p/FLOW[boundary_node[i]].t);
			FLOW[boundary_node[i]].e = FLOW[boundary_node[i]].p/(0.4*FLOW[boundary_node[i]].rho);
		}	
	}
	
	/*****************************************************************************************************************************************************************/
	
	for (i=0; i< wal_node; i++)
	{	
		if (node[wall_node[i]].loc > 0 && node[wall_node[i]].loc <= 6)
		{	
			if (node[wall_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[wall_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[wall_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[wall_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[wall_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[wall_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[wall_node[i]].n_n[k]].p;
			FLOW[wall_node[i]].t = FLOW[node[wall_node[i]].n_n[k]].t;	
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);		
			FLOW[wall_node[i]].a = sqrt(1.4*FLOW[wall_node[i]].p/FLOW[wall_node[i]].rho);	
		
			FLOW[node[wall_node[i]].n_n[l]].u = (-1.0)*FLOW[node[wall_node[i]].n_n[k]].u;
			FLOW[node[wall_node[i]].n_n[l]].v = (-1.0)*FLOW[node[wall_node[i]].n_n[k]].v;
			FLOW[node[wall_node[i]].n_n[l]].w = (-1.0)*FLOW[node[wall_node[i]].n_n[k]].w;
			FLOW[node[wall_node[i]].n_n[l]].p = FLOW[node[wall_node[i]].n_n[k]].p;
			FLOW[node[wall_node[i]].n_n[l]].rho = FLOW[node[wall_node[i]].n_n[k]].rho;
			FLOW[node[wall_node[i]].n_n[l]].t = FLOW[node[wall_node[i]].n_n[k]].t;
			FLOW[node[wall_node[i]].n_n[l]].a = FLOW[node[wall_node[i]].n_n[k]].a;
			FLOW[node[wall_node[i]].n_n[l]].e = FLOW[node[wall_node[i]].n_n[k]].e;
			//mFLOW[node[wall_node[i]].n_n[l]] = mFLOW[wall_node[i]];
		
			FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[node[wall_node[i]].n_n[k]].n_n[k]].u;
			FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[node[wall_node[i]].n_n[k]].n_n[k]].v;
			FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[node[wall_node[i]].n_n[k]].n_n[k]].w;
			FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].p = FLOW[node[node[wall_node[i]].n_n[k]].n_n[k]].p;
			FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].rho = FLOW[node[node[wall_node[i]].n_n[k]].n_n[k]].rho;
			FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].t = FLOW[node[node[wall_node[i]].n_n[k]].n_n[k]].t;
			FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].a = FLOW[node[node[wall_node[i]].n_n[k]].n_n[k]].a;
			FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].e = FLOW[node[node[wall_node[i]].n_n[k]].n_n[k]].e;
			//mFLOW[node[node[wall_node[i]].n_n[l]].n_n[l]] = mFLOW[wall_node[i]];
		
			FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].u;
			FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].v;
			FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].w;
			FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].p;
			FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].rho;
			FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].t;
			FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].a;
			FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].e;
			//mFLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mFLOW[wall_node[i]];
		
		}
		else if (node[wall_node[i]].loc == 7 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[2]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[2]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);

		}
		else if (node[wall_node[i]].loc == 8 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[2]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[2]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);			
		}
		else if (node[wall_node[i]].loc == 9 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[2]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[2]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);			
		}
		else if (node[wall_node[i]].loc == 10 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[2]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[2]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		
		else if (node[wall_node[i]].loc == 11 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[0]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[0]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 12 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[0]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[0]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);			
		}
		else if (node[wall_node[i]].loc == 13 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[0]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[0]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);		
		}
		else if (node[wall_node[i]].loc == 14 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[0]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[0]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);	
		}
		else if (node[wall_node[i]].loc == 15 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[5]].n_n[0]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[5]].n_n[0]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 16 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[1]].n_n[0]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[1]].n_n[0]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 17 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[4]].n_n[0]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[4]].n_n[0]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 18 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[3]].n_n[0]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[3]].n_n[0]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 19 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[5]].n_n[2]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[5]].n_n[2]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 20 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[1]].n_n[2]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[1]].n_n[2]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 21 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[4]].n_n[2]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[4]].n_n[2]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 22 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[3]].n_n[2]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[3]].n_n[2]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 23 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[1]].n_n[5]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[1]].n_n[5]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 24 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[1]].n_n[4]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[1]].n_n[4]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 25 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[3]].n_n[4]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[3]].n_n[4]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 26 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[node[wall_node[i]].n_n[3]].n_n[5]].p;
			FLOW[wall_node[i]].t = FLOW[node[node[wall_node[i]].n_n[3]].n_n[5]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}	
		else if (node[wall_node[i]].loc == 29 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[singular[wall_node[i]].n_n[3]].n_n[4]].p;
			FLOW[wall_node[i]].t = FLOW[node[singular[wall_node[i]].n_n[3]].n_n[4]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		else if (node[wall_node[i]].loc == 30 )
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[node[singular[wall_node[i]].n_n[3]].n_n[5]].p;
			FLOW[wall_node[i]].t = FLOW[node[singular[wall_node[i]].n_n[3]].n_n[5]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);
		}
		
		if (node[wall_node[i]].ID == 10 && node[wall_node[i]].loc != 29 && node[wall_node[i]].loc != 30)
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[singular[wall_node[i]].n_n[3]].p;
			FLOW[wall_node[i]].t = FLOW[singular[wall_node[i]].n_n[3]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);		
			
			for (h=0;h<3;h++)
			{
				if (h==0)
				{
					l = 2; /********no element on SOUTH**************/
					k = 0;
				}
				if (h==1) 
				{
					l = 0; /********no element on NORTH**************/
					k = 2;
				}
				if (h==2)
				{
					l = 1; /********no element on EAST**************/
					k = 3;										
				}
				
				FLOW[node[wall_node[i]].n_n[l]].u = (-1.0)*FLOW[singular[wall_node[i]].n_n[k]].u;
				FLOW[node[wall_node[i]].n_n[l]].v = (-1.0)*FLOW[singular[wall_node[i]].n_n[k]].v;
				FLOW[node[wall_node[i]].n_n[l]].w = (-1.0)*FLOW[singular[wall_node[i]].n_n[k]].w;
				FLOW[node[wall_node[i]].n_n[l]].p = FLOW[singular[wall_node[i]].n_n[k]].p;
				FLOW[node[wall_node[i]].n_n[l]].rho = FLOW[singular[wall_node[i]].n_n[k]].rho;
				FLOW[node[wall_node[i]].n_n[l]].t = FLOW[singular[wall_node[i]].n_n[k]].t;	
				FLOW[node[wall_node[i]].n_n[l]].a = FLOW[wall_node[i]].a;
				FLOW[node[wall_node[i]].n_n[l]].e = FLOW[singular[wall_node[i]].n_n[k]].e;
				//mu[RK][node[wall_node[i]].n_n[l]] = mu[RK][wall_node[i]];
			
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].u;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].v;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].w;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].p = FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].p;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].rho = FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].rho;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].t = FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].t;	
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].a = FLOW[wall_node[i]].a;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].e = FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].e;
				//mu[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];
	
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].u;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].v;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].w;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].p;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].rho;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].t;	
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[wall_node[i]].a;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].e;
				//mu[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];	
			}
		}
		
		if (node[wall_node[i]].loc >= 33 && node[wall_node[i]].loc <= 56)
		{
			FLOW[wall_node[i]].u = 0.0;
			FLOW[wall_node[i]].v = 0.0;
			FLOW[wall_node[i]].w = 0.0;
			FLOW[wall_node[i]].p = FLOW[singular[wall_node[i]].n_n[0]].p;
			FLOW[wall_node[i]].t = FLOW[singular[wall_node[i]].n_n[0]].t;
			FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
			FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);		
			
			if ((node[wall_node[i]].loc >= 33 && node[wall_node[i]].loc <= 36) || (node[wall_node[i]].loc >= 41 && node[wall_node[i]].loc <= 44) || (node[wall_node[i]].loc >= 49 && node[wall_node[i]].loc <= 52))
			{
				d_temp6 = 2;
			}
			if ((node[wall_node[i]].loc >= 37 && node[wall_node[i]].loc <= 40) || (node[wall_node[i]].loc >= 45 && node[wall_node[i]].loc <= 48))
			{
				d_temp6 = 3;
			}
			if (node[wall_node[i]].loc >= 53 && node[wall_node[i]].loc <= 56)
			{
				d_temp6 = 1;
			}
		
			for (d_temp = 0; d_temp<d_temp6; d_temp++)
			{
				if (d_temp == 0 && ((node[wall_node[i]].loc >= 33 && node[wall_node[i]].loc <= 36) || (node[wall_node[i]].loc >= 37 && node[wall_node[i]].loc <= 40)))
				{
					l = 2; 
					k = 0;
				}
				if (d_temp == 0 && ((node[wall_node[i]].loc >= 41 && node[wall_node[i]].loc <= 44) || (node[wall_node[i]].loc >= 45 && node[wall_node[i]].loc <= 48)))
				{
					l = 0; 
					k = 2;
				}
				if (d_temp == 0 && (node[wall_node[i]].loc == 49 || node[wall_node[i]].loc == 52))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 0 && (node[wall_node[i]].loc == 50 || node[wall_node[i]].loc == 51))
				{
					l = 4; 
					k = 5;
				}
				if (d_temp == 1 && (node[wall_node[i]].loc == 33 || node[wall_node[i]].loc == 43 || node[wall_node[i]].loc == 51 || node[wall_node[i]].loc == 52 || node[wall_node[i]].loc == 39 || node[wall_node[i]].loc == 40 || node[wall_node[i]].loc == 47 || node[wall_node[i]].loc == 48))
				{
					l = 3; 
					k = 1;
				}
				if (d_temp == 1 && (node[wall_node[i]].loc == 35 || node[wall_node[i]].loc == 41 || node[wall_node[i]].loc == 49 || node[wall_node[i]].loc == 50 || node[wall_node[i]].loc == 37 || node[wall_node[i]].loc == 38 || node[wall_node[i]].loc == 45 || node[wall_node[i]].loc == 46  ))
				{
					l = 1; 
					k = 3;
				}
				if (d_temp == 1 && (node[wall_node[i]].loc == 34 || node[wall_node[i]].loc == 44))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 1 && (node[wall_node[i]].loc == 36 || node[wall_node[i]].loc == 42))
				{
					l = 4; 
					k = 5;
				}
				if (d_temp == 2 && (node[wall_node[i]].loc == 37 || node[wall_node[i]].loc == 40 || node[wall_node[i]].loc == 48 || node[wall_node[i]].loc == 45))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 2 && (node[wall_node[i]].loc == 38 || node[wall_node[i]].loc == 39 || node[wall_node[i]].loc == 46 || node[wall_node[i]].loc == 47))
				{
					l = 4; 
					k = 5;
				}
				
				FLOW[wall_node[i]].u = 0.0;
				FLOW[wall_node[i]].v = 0.0;
				FLOW[wall_node[i]].w = 0.0;
				FLOW[wall_node[i]].p = FLOW[singular[wall_node[i]].n_n[k]].p;
				FLOW[wall_node[i]].t = FLOW[singular[wall_node[i]].n_n[k]].t;
				FLOW[wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[wall_node[i]].p/FLOW[wall_node[i]].t);
				FLOW[wall_node[i]].e = FLOW[wall_node[i]].p/(0.4*FLOW[wall_node[i]].rho);	
				
				
				FLOW[node[wall_node[i]].n_n[l]].u = (-1.0)*FLOW[singular[wall_node[i]].n_n[k]].u;
				FLOW[node[wall_node[i]].n_n[l]].v = (-1.0)*FLOW[singular[wall_node[i]].n_n[k]].v;
				FLOW[node[wall_node[i]].n_n[l]].w = (-1.0)*FLOW[singular[wall_node[i]].n_n[k]].w;
				FLOW[node[wall_node[i]].n_n[l]].p = FLOW[singular[wall_node[i]].n_n[k]].p;
				FLOW[node[wall_node[i]].n_n[l]].rho = FLOW[singular[wall_node[i]].n_n[k]].rho;
				FLOW[node[wall_node[i]].n_n[l]].t = FLOW[singular[wall_node[i]].n_n[k]].t;	
				FLOW[node[wall_node[i]].n_n[l]].a = FLOW[wall_node[i]].a;
				FLOW[node[wall_node[i]].n_n[l]].e = FLOW[singular[wall_node[i]].n_n[k]].e;
				//mFLOW[node[wall_node[i]].n_n[l]] = mFLOW[wall_node[i]];
			
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].u;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].v;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].w;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].p = FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].p;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].rho = FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].rho;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].t = FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].t;	
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].a = FLOW[wall_node[i]].a;
				FLOW[node[node[wall_node[i]].n_n[l]].n_n[l]].e = FLOW[node[singular[wall_node[i]].n_n[k]].n_n[k]].e;
				//mFLOW[node[node[wall_node[i]].n_n[l]].n_n[l]] = mFLOW[wall_node[i]];
	
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].u;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].v;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].w;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].p;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].rho;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].t;	
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[wall_node[i]].a;
				FLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]].e;
				//mFLOW[node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mFLOW[wall_node[i]];	
				
			}
		}	
		
	}
	
	/****************************************************SD_OUTLET********************************************************/		
	for (i=0; i<sd_out_node; i++)
	{	
		if (node[sd_outlet_node[i]].loc > 0 && node[sd_outlet_node[i]].loc <= 6)
		{	
			if (node[sd_outlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_outlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_outlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_outlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_outlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_outlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			FLOW[sd_outlet_node[i]].u= FLOW[node[sd_outlet_node[i]].n_n[k]].u;	
			FLOW[sd_outlet_node[i]].v= FLOW[node[sd_outlet_node[i]].n_n[k]].v;	
			FLOW[sd_outlet_node[i]].w= FLOW[node[sd_outlet_node[i]].n_n[k]].w;	
			FLOW[sd_outlet_node[i]].p= FLOW[node[sd_outlet_node[i]].n_n[k]].p;
		/* 	if (back_pressure != 0.0)
			{
				if (node[sd_outlet_node[i]].y < 412.0 || node[sd_outlet_node[i]].y > 767.0)
				{
					FLOW[sd_outlet_node[i]].p = FLOW[node[sd_outlet_node[i]].n_n[k]].p;	
				}
				if (node[sd_outlet_node[i]].y > 412.0 && node[sd_outlet_node[i]].y < 767.0)
				{
					FLOW[sd_outlet_node[i]].p = back_pressure*(1.0/(1.4*Mach*Mach));	
				}
			}	 */			
			FLOW[sd_outlet_node[i]].t= FLOW[node[sd_outlet_node[i]].n_n[k]].t;	
			FLOW[sd_outlet_node[i]].rho= (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
			FLOW[sd_outlet_node[i]].a = sqrt(1.4*FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].rho); 

			FLOW[node[sd_outlet_node[i]].n_n[l]].u = FLOW[sd_outlet_node[i]].u;
			FLOW[node[sd_outlet_node[i]].n_n[l]].v = FLOW[sd_outlet_node[i]].v;
			FLOW[node[sd_outlet_node[i]].n_n[l]].w = FLOW[sd_outlet_node[i]].w;
			FLOW[node[sd_outlet_node[i]].n_n[l]].p = FLOW[sd_outlet_node[i]].p;
			FLOW[node[sd_outlet_node[i]].n_n[l]].rho = FLOW[sd_outlet_node[i]].rho;
			FLOW[node[sd_outlet_node[i]].n_n[l]].t = FLOW[sd_outlet_node[i]].t;	
			FLOW[node[sd_outlet_node[i]].n_n[l]].a = FLOW[sd_outlet_node[i]].a;
			FLOW[node[sd_outlet_node[i]].n_n[l]].e = FLOW[sd_outlet_node[i]].e;
			
			FLOW[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].u = FLOW[sd_outlet_node[i]].u;
			FLOW[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].v = FLOW[sd_outlet_node[i]].v;
			FLOW[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].w = FLOW[sd_outlet_node[i]].w;
			FLOW[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].p = FLOW[sd_outlet_node[i]].p;
			FLOW[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].rho = FLOW[sd_outlet_node[i]].rho;
			FLOW[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].t = FLOW[sd_outlet_node[i]].t;	
			FLOW[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].a = FLOW[sd_outlet_node[i]].a;
			FLOW[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].e = FLOW[sd_outlet_node[i]].e;
	
			FLOW[node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].u = FLOW[sd_outlet_node[i]].u;
			FLOW[node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].v = FLOW[sd_outlet_node[i]].v;
			FLOW[node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].w = FLOW[sd_outlet_node[i]].w;
			FLOW[node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[sd_outlet_node[i]].p;
			FLOW[node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[sd_outlet_node[i]].rho;
			FLOW[node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[sd_outlet_node[i]].t;
			FLOW[node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[sd_outlet_node[i]].a;
			FLOW[node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[sd_outlet_node[i]].e;
			
		}
		else if (node[sd_outlet_node[i]].loc == 7 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);

		}
		else if (node[sd_outlet_node[i]].loc == 8 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);			
		}
		else if (node[sd_outlet_node[i]].loc == 9 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);			
		}
		else if (node[sd_outlet_node[i]].loc == 10 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		
		else if (node[sd_outlet_node[i]].loc == 11 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 12 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);			
		}
		else if (node[sd_outlet_node[i]].loc == 13 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);		
		}
		else if (node[sd_outlet_node[i]].loc == 14 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);	
		}
		else if (node[sd_outlet_node[i]].loc == 15 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[5]].n_n[0]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[5]].n_n[0]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[5]].n_n[0]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[5]].n_n[0]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[5]].n_n[0]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 16 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[0]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[0]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[0]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[0]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[0]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 17 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[4]].n_n[0]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[4]].n_n[0]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[4]].n_n[0]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[4]].n_n[0]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[4]].n_n[0]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 18 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[0]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[0]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[0]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[0]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[0]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 19 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[5]].n_n[2]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[5]].n_n[2]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[5]].n_n[2]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[5]].n_n[2]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[5]].n_n[2]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 20 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[2]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[2]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[2]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[2]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[2]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 21 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[4]].n_n[2]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[4]].n_n[2]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[4]].n_n[2]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[4]].n_n[2]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[4]].n_n[2]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 22 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[2]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[2]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[2]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[2]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[2]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 23 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 24 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 25 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
		else if (node[sd_outlet_node[i]].loc == 26 )
		{
			FLOW[sd_outlet_node[i]].u = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].u;
			FLOW[sd_outlet_node[i]].v = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].v;
			FLOW[sd_outlet_node[i]].w = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].w;
			FLOW[sd_outlet_node[i]].p = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].p;
			FLOW[sd_outlet_node[i]].t = FLOW[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].t;
			FLOW[sd_outlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_outlet_node[i]].p/FLOW[sd_outlet_node[i]].t);
			FLOW[sd_outlet_node[i]].e = FLOW[sd_outlet_node[i]].p/(0.4*FLOW[sd_outlet_node[i]].rho);
		}
	}

/****************************************************INLET********************************************************/
	for (i=0; i< sd_inl_node; i++)
	{	
		if (node[sd_inlet_node[i]].loc > 0 && node[sd_inlet_node[i]].loc <= 6)
		{	
			if (node[sd_inlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_inlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_inlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_inlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_inlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_inlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}
		
			FLOW[sd_inlet_node[i]].u = 1.0;		
			FLOW[sd_inlet_node[i]].v = 0.0;	
			FLOW[sd_inlet_node[i]].w = 0.0;				
			FLOW[sd_inlet_node[i]].p = 1.0/(1.4*Mach*Mach);		
			FLOW[sd_inlet_node[i]].t = 1.0;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
			FLOW[sd_inlet_node[i]].a = sqrt(1.4*FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].rho); 

			FLOW[node[sd_inlet_node[i]].n_n[l]].u = FLOW[sd_inlet_node[i]].u;
			FLOW[node[sd_inlet_node[i]].n_n[l]].v = FLOW[sd_inlet_node[i]].v;
			FLOW[node[sd_inlet_node[i]].n_n[l]].w = FLOW[sd_inlet_node[i]].w;
			FLOW[node[sd_inlet_node[i]].n_n[l]].p = FLOW[sd_inlet_node[i]].p;
			FLOW[node[sd_inlet_node[i]].n_n[l]].rho = FLOW[sd_inlet_node[i]].rho;
			FLOW[node[sd_inlet_node[i]].n_n[l]].t = FLOW[sd_inlet_node[i]].t;	
			FLOW[node[sd_inlet_node[i]].n_n[l]].a = FLOW[sd_inlet_node[i]].a;
			FLOW[node[sd_inlet_node[i]].n_n[l]].e = FLOW[sd_inlet_node[i]].e;
			
			FLOW[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].u = FLOW[sd_inlet_node[i]].u;
			FLOW[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].v = FLOW[sd_inlet_node[i]].v;
			FLOW[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].w = FLOW[sd_inlet_node[i]].w;
			FLOW[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].p = FLOW[sd_inlet_node[i]].p;
			FLOW[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].rho = FLOW[sd_inlet_node[i]].rho;
			FLOW[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].t = FLOW[sd_inlet_node[i]].t;	
			FLOW[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].a = FLOW[sd_inlet_node[i]].a;
			FLOW[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].e = FLOW[sd_inlet_node[i]].e;
	
			FLOW[node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].u = FLOW[sd_inlet_node[i]].u;
			FLOW[node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].v = FLOW[sd_inlet_node[i]].v;
			FLOW[node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].w = FLOW[sd_inlet_node[i]].w;
			FLOW[node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[sd_inlet_node[i]].p;
			FLOW[node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[sd_inlet_node[i]].rho;
			FLOW[node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[sd_inlet_node[i]].t;
			FLOW[node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[sd_inlet_node[i]].a;
			FLOW[node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[sd_inlet_node[i]].e;
		
		}
		else if (node[sd_inlet_node[i]].loc == 7 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);

		}
		else if (node[sd_inlet_node[i]].loc == 8 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);			
		}
		else if (node[sd_inlet_node[i]].loc == 9 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);			
		}
		else if (node[sd_inlet_node[i]].loc == 10 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		
		else if (node[sd_inlet_node[i]].loc == 11 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 12 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);			
		}
		else if (node[sd_inlet_node[i]].loc == 13 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);		
		}
		else if (node[sd_inlet_node[i]].loc == 14 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);	
		}
		else if (node[sd_inlet_node[i]].loc == 15 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[5]].n_n[0]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[5]].n_n[0]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[5]].n_n[0]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[5]].n_n[0]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[5]].n_n[0]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 16 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[0]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[0]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[0]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[0]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[0]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 17 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[4]].n_n[0]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[4]].n_n[0]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[4]].n_n[0]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[4]].n_n[0]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[4]].n_n[0]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 18 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[0]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[0]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[0]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[0]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[0]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 19 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[5]].n_n[2]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[5]].n_n[2]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[5]].n_n[2]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[5]].n_n[2]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[5]].n_n[2]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 20 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[2]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[2]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[2]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[2]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[2]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 21 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[4]].n_n[2]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[4]].n_n[2]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[4]].n_n[2]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[4]].n_n[2]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[4]].n_n[2]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 22 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[2]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[2]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[2]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[2]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[2]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 23 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 24 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 25 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
		else if (node[sd_inlet_node[i]].loc == 26 )
		{
			FLOW[sd_inlet_node[i]].u = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].u;
			FLOW[sd_inlet_node[i]].v = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].v;
			FLOW[sd_inlet_node[i]].w = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].w;
			FLOW[sd_inlet_node[i]].p = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].p;
			FLOW[sd_inlet_node[i]].t = FLOW[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].t;
			FLOW[sd_inlet_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_inlet_node[i]].p/FLOW[sd_inlet_node[i]].t);
			FLOW[sd_inlet_node[i]].e = FLOW[sd_inlet_node[i]].p/(0.4*FLOW[sd_inlet_node[i]].rho);
		}
	}
	
	/****************************************************BOUNDARY********************************************************/
	for (i=0; i< sd_bou_node; i++)
	{	
		if (node[sd_boundary_node[i]].loc > 0 && node[sd_boundary_node[i]].loc <= 6)
		{	
			if (node[sd_boundary_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_boundary_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_boundary_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_boundary_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_boundary_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_boundary_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			FLOW[sd_boundary_node[i]].u = FLOW[node[sd_boundary_node[i]].n_n[k]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[sd_boundary_node[i]].n_n[k]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[sd_boundary_node[i]].n_n[k]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[sd_boundary_node[i]].n_n[k]].p;
			FLOW[sd_boundary_node[i]].rho = FLOW[node[sd_boundary_node[i]].n_n[k]].rho;
			FLOW[sd_boundary_node[i]].t = FLOW[node[sd_boundary_node[i]].n_n[k]].t;	
			FLOW[sd_boundary_node[i]].a = FLOW[node[sd_boundary_node[i]].n_n[k]].a;
			FLOW[sd_boundary_node[i]].e = FLOW[node[sd_boundary_node[i]].n_n[k]].e;
			//mu[RK][sd_boundary_node[i]] = mu[RK][sd_boundary_node[i]];
			
			FLOW[node[sd_boundary_node[i]].n_n[l]].u = FLOW[node[sd_boundary_node[i]].n_n[k]].u;
			FLOW[node[sd_boundary_node[i]].n_n[l]].v = FLOW[node[sd_boundary_node[i]].n_n[k]].v;
			FLOW[node[sd_boundary_node[i]].n_n[l]].w = FLOW[node[sd_boundary_node[i]].n_n[k]].w;
			FLOW[node[sd_boundary_node[i]].n_n[l]].p = FLOW[node[sd_boundary_node[i]].n_n[k]].p;
			FLOW[node[sd_boundary_node[i]].n_n[l]].rho = FLOW[node[sd_boundary_node[i]].n_n[k]].rho;
			FLOW[node[sd_boundary_node[i]].n_n[l]].t = FLOW[node[sd_boundary_node[i]].n_n[k]].t;	
			FLOW[node[sd_boundary_node[i]].n_n[l]].a = FLOW[node[sd_boundary_node[i]].n_n[k]].a;
			FLOW[node[sd_boundary_node[i]].n_n[l]].e = FLOW[node[sd_boundary_node[i]].n_n[k]].e;
			//mFLOW[node[sd_boundary_node[i]].n_n[l]] = mFLOW[sd_boundary_node[i]];
			
			FLOW[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].u = FLOW[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].u;
			FLOW[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].v = FLOW[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].v;
			FLOW[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].w = FLOW[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].w;
			FLOW[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].p = FLOW[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].p;
			FLOW[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].rho = FLOW[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].rho;
			FLOW[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].t = FLOW[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].t;	
			FLOW[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].a = FLOW[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].a;
			FLOW[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].e = FLOW[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].e;
			//mFLOW[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = mFLOW[sd_boundary_node[i]];
	
			FLOW[node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].u = FLOW[node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].u;
			FLOW[node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].v = FLOW[node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].v;
			FLOW[node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].w = FLOW[node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].w;
			FLOW[node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].p;
			FLOW[node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].rho;
			FLOW[node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].t;	
			FLOW[node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].a;
			FLOW[node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]].e;
			//mFLOW[node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = mFLOW[sd_boundary_node[i]];
		}
		else if (node[sd_boundary_node[i]].loc == 7 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);

		}
		else if (node[sd_boundary_node[i]].loc == 8 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);			
		}
		else if (node[sd_boundary_node[i]].loc == 9 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);			
		}
		else if (node[sd_boundary_node[i]].loc == 10 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		
		else if (node[sd_boundary_node[i]].loc == 11 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 12 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);			
		}
		else if (node[sd_boundary_node[i]].loc == 13 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);		
		}
		else if (node[sd_boundary_node[i]].loc == 14 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);	
		}
		else if (node[sd_boundary_node[i]].loc == 15 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[5]].n_n[0]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[5]].n_n[0]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[5]].n_n[0]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[5]].n_n[0]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[5]].n_n[0]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 16 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[0]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[0]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[0]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[0]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[0]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 17 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[4]].n_n[0]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[4]].n_n[0]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[4]].n_n[0]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[4]].n_n[0]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[4]].n_n[0]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 18 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[0]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[0]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[0]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[0]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[0]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 19 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[5]].n_n[2]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[5]].n_n[2]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[5]].n_n[2]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[5]].n_n[2]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[5]].n_n[2]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 20 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[2]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[2]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[2]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[2]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[2]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 21 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[4]].n_n[2]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[4]].n_n[2]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[4]].n_n[2]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[4]].n_n[2]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[4]].n_n[2]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 22 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[2]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[2]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[2]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[2]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[2]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 23 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 24 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 25 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}
		else if (node[sd_boundary_node[i]].loc == 26 )
		{
			FLOW[sd_boundary_node[i]].u = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].u;
			FLOW[sd_boundary_node[i]].v = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].v;
			FLOW[sd_boundary_node[i]].w = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].w;
			FLOW[sd_boundary_node[i]].p = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].p;
			FLOW[sd_boundary_node[i]].t = FLOW[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].t;
			FLOW[sd_boundary_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_boundary_node[i]].p/FLOW[sd_boundary_node[i]].t);
			FLOW[sd_boundary_node[i]].e = FLOW[sd_boundary_node[i]].p/(0.4*FLOW[sd_boundary_node[i]].rho);
		}	
	}
	
	/*****************************************************************************************************************************************************************/
	
	for (i=0; i< sd_wal_node; i++)
	{	
		if (node[sd_wall_node[i]].loc > 0 && node[sd_wall_node[i]].loc <= 6)
		{	
			if (node[sd_wall_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_wall_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_wall_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_wall_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_wall_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_wall_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[sd_wall_node[i]].n_n[k]].p;
			//FLOW[sd_wall_node[i]]= 1.0/(1.4*Mach*Mach);
			FLOW[sd_wall_node[i]].t = FLOW[node[sd_wall_node[i]].n_n[k]].t;	
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);		
			FLOW[sd_wall_node[i]].a = sqrt(1.4*FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].rho);	
		
			FLOW[node[sd_wall_node[i]].n_n[l]].u = (-1.0)*FLOW[node[sd_wall_node[i]].n_n[k]].u;
			FLOW[node[sd_wall_node[i]].n_n[l]].v = (-1.0)*FLOW[node[sd_wall_node[i]].n_n[k]].v;
			FLOW[node[sd_wall_node[i]].n_n[l]].w = (-1.0)*FLOW[node[sd_wall_node[i]].n_n[k]].w;
			FLOW[node[sd_wall_node[i]].n_n[l]].p = FLOW[node[sd_wall_node[i]].n_n[k]].p;
			FLOW[node[sd_wall_node[i]].n_n[l]].rho = FLOW[node[sd_wall_node[i]].n_n[k]].rho;
			FLOW[node[sd_wall_node[i]].n_n[l]].t = FLOW[node[sd_wall_node[i]].n_n[k]].t;
			FLOW[node[sd_wall_node[i]].n_n[l]].a = FLOW[node[sd_wall_node[i]].n_n[k]].a;
			FLOW[node[sd_wall_node[i]].n_n[l]].e = FLOW[node[sd_wall_node[i]].n_n[k]].e;
			//mFLOW[node[sd_wall_node[i]].n_n[l]] = mFLOW[sd_wall_node[i]];
		
			FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].u;
			FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].v;
			FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].w;
			FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].p = FLOW[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].p;
			FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].rho = FLOW[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].rho;
			FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].t = FLOW[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].t;
			FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].a = FLOW[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].a;
			FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].e = FLOW[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].e;
			//mFLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = mFLOW[sd_wall_node[i]];
		
			FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].u;
			FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].v;
			FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].w;
			FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].p;
			FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].rho;
			FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].t;
			FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].a;
			FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].e;
			//mFLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mFLOW[sd_wall_node[i]];
			
		}
		else if (node[sd_wall_node[i]].loc == 7 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[2]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[2]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);

		}
		else if (node[sd_wall_node[i]].loc == 8 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[2]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[2]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);			
		}
		else if (node[sd_wall_node[i]].loc == 9 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[2]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[2]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);			
		}
		else if (node[sd_wall_node[i]].loc == 10 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[2]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[2]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		
		else if (node[sd_wall_node[i]].loc == 11 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[0]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[0]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 12 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[0]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[0]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);			
		}
		else if (node[sd_wall_node[i]].loc == 13 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[0]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[0]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);		
		}
		else if (node[sd_wall_node[i]].loc == 14 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[0]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[0]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);	
		}
		else if (node[sd_wall_node[i]].loc == 15 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[5]].n_n[0]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[5]].n_n[0]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 16 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[1]].n_n[0]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[1]].n_n[0]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 17 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[4]].n_n[0]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[4]].n_n[0]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 18 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[3]].n_n[0]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[3]].n_n[0]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 19 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[5]].n_n[2]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[5]].n_n[2]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 20 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[1]].n_n[2]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[1]].n_n[2]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 21 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[4]].n_n[2]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[4]].n_n[2]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 22 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[3]].n_n[2]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[3]].n_n[2]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 23 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 24 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 25 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 26 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}	
		else if (node[sd_wall_node[i]].loc == 29 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[singular[sd_wall_node[i]].n_n[3]].n_n[4]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[singular[sd_wall_node[i]].n_n[3]].n_n[4]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		else if (node[sd_wall_node[i]].loc == 30 )
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[node[singular[sd_wall_node[i]].n_n[3]].n_n[5]].p;
			FLOW[sd_wall_node[i]].t = FLOW[node[singular[sd_wall_node[i]].n_n[3]].n_n[5]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);
		}
		
		if (node[sd_wall_node[i]].ID == 10 && node[sd_wall_node[i]].loc != 29 && node[sd_wall_node[i]].loc != 30)
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[singular[sd_wall_node[i]].n_n[3]].p;
			FLOW[sd_wall_node[i]].t = FLOW[singular[sd_wall_node[i]].n_n[3]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);		
			
			for (h=0;h<3;h++)
			{
				if (h==0)
				{
					l = 2; /********no element on SOUTH**************/
					k = 0;
				}
				if (h==1) 
				{
					l = 0; /********no element on NORTH**************/
					k = 2;
				}
				if (h==2)
				{
					l = 1; /********no element on EAST**************/
					k = 3;										
				}
				
				FLOW[node[sd_wall_node[i]].n_n[l]].u = (-1.0)*FLOW[singular[sd_wall_node[i]].n_n[k]].u;
				FLOW[node[sd_wall_node[i]].n_n[l]].v = (-1.0)*FLOW[singular[sd_wall_node[i]].n_n[k]].v;
				FLOW[node[sd_wall_node[i]].n_n[l]].w = (-1.0)*FLOW[singular[sd_wall_node[i]].n_n[k]].w;
				FLOW[node[sd_wall_node[i]].n_n[l]].p = FLOW[singular[sd_wall_node[i]].n_n[k]].p;
				FLOW[node[sd_wall_node[i]].n_n[l]].rho = FLOW[singular[sd_wall_node[i]].n_n[k]].rho;
				FLOW[node[sd_wall_node[i]].n_n[l]].t = FLOW[singular[sd_wall_node[i]].n_n[k]].t;	
				FLOW[node[sd_wall_node[i]].n_n[l]].a = FLOW[sd_wall_node[i]].a;
				FLOW[node[sd_wall_node[i]].n_n[l]].e = FLOW[singular[sd_wall_node[i]].n_n[k]].e;
				//mu[RK][node[sd_wall_node[i]].n_n[l]] = mu[RK][sd_wall_node[i]];
			
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].u;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].v;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].w;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].p = FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].p;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].rho = FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].rho;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].t = FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].t;	
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].a = FLOW[sd_wall_node[i]].a;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].e = FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].e;
				//mu[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];
	
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].u;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].v;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].w;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].p;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].rho;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].t;	
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[sd_wall_node[i]].a;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].e;
				//mu[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];	
			}
		}
		
		if (node[sd_wall_node[i]].loc >= 33 && node[sd_wall_node[i]].loc <= 56)
		{
			FLOW[sd_wall_node[i]].u = 0.0;
			FLOW[sd_wall_node[i]].v = 0.0;
			FLOW[sd_wall_node[i]].w = 0.0;
			FLOW[sd_wall_node[i]].p = FLOW[singular[sd_wall_node[i]].n_n[0]].p;
			FLOW[sd_wall_node[i]].t = FLOW[singular[sd_wall_node[i]].n_n[0]].t;
			FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
			FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);		
			
			if ((node[sd_wall_node[i]].loc >= 33 && node[sd_wall_node[i]].loc <= 36) || (node[sd_wall_node[i]].loc >= 41 && node[sd_wall_node[i]].loc <= 44) || (node[sd_wall_node[i]].loc >= 49 && node[sd_wall_node[i]].loc <= 52))
			{
				d_temp6 = 2;
			}
			if ((node[sd_wall_node[i]].loc >= 37 && node[sd_wall_node[i]].loc <= 40) || (node[sd_wall_node[i]].loc >= 45 && node[sd_wall_node[i]].loc <= 48))
			{
				d_temp6 = 3;
			}
			if (node[sd_wall_node[i]].loc >= 53 && node[sd_wall_node[i]].loc <= 56)
			{
				d_temp6 = 1;
			}
		
			for (d_temp = 0; d_temp<d_temp6; d_temp++)
			{
				if (d_temp == 0 && ((node[sd_wall_node[i]].loc >= 33 && node[sd_wall_node[i]].loc <= 36) || (node[sd_wall_node[i]].loc >= 37 && node[sd_wall_node[i]].loc <= 40)))
				{
					l = 2; 
					k = 0;
				}
				if (d_temp == 0 && ((node[sd_wall_node[i]].loc >= 41 && node[sd_wall_node[i]].loc <= 44) || (node[sd_wall_node[i]].loc >= 45 && node[sd_wall_node[i]].loc <= 48)))
				{
					l = 0; 
					k = 2;
				}
				if (d_temp == 0 && (node[sd_wall_node[i]].loc == 49 || node[sd_wall_node[i]].loc == 52))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 0 && (node[sd_wall_node[i]].loc == 50 || node[sd_wall_node[i]].loc == 51))
				{
					l = 4; 
					k = 5;
				}
				if (d_temp == 1 && (node[sd_wall_node[i]].loc == 33 || node[sd_wall_node[i]].loc == 43 || node[sd_wall_node[i]].loc == 51 || node[sd_wall_node[i]].loc == 52 || node[sd_wall_node[i]].loc == 39 || node[sd_wall_node[i]].loc == 40 || node[sd_wall_node[i]].loc == 47 || node[sd_wall_node[i]].loc == 48))
				{
					l = 3; 
					k = 1;
				}
				if (d_temp == 1 && (node[sd_wall_node[i]].loc == 35 || node[sd_wall_node[i]].loc == 41 || node[sd_wall_node[i]].loc == 49 || node[sd_wall_node[i]].loc == 50 || node[sd_wall_node[i]].loc == 37 || node[sd_wall_node[i]].loc == 38 || node[sd_wall_node[i]].loc == 45 || node[sd_wall_node[i]].loc == 46  ))
				{
					l = 1; 
					k = 3;
				}
				if (d_temp == 1 && (node[sd_wall_node[i]].loc == 34 || node[sd_wall_node[i]].loc == 44))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 1 && (node[sd_wall_node[i]].loc == 36 || node[sd_wall_node[i]].loc == 42))
				{
					l = 4; 
					k = 5;
				}
				if (d_temp == 2 && (node[sd_wall_node[i]].loc == 37 || node[sd_wall_node[i]].loc == 40 || node[sd_wall_node[i]].loc == 48 || node[sd_wall_node[i]].loc == 45))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 2 && (node[sd_wall_node[i]].loc == 38 || node[sd_wall_node[i]].loc == 39 || node[sd_wall_node[i]].loc == 46 || node[sd_wall_node[i]].loc == 47))
				{
					l = 4; 
					k = 5;
				}
				
				FLOW[sd_wall_node[i]].u = 0.0;
				FLOW[sd_wall_node[i]].v = 0.0;
				FLOW[sd_wall_node[i]].w = 0.0;
				FLOW[sd_wall_node[i]].p = FLOW[singular[sd_wall_node[i]].n_n[k]].p;
				FLOW[sd_wall_node[i]].t = FLOW[singular[sd_wall_node[i]].n_n[k]].t;
				FLOW[sd_wall_node[i]].rho = (1.4*Mach*Mach)*(FLOW[sd_wall_node[i]].p/FLOW[sd_wall_node[i]].t);
				FLOW[sd_wall_node[i]].e = FLOW[sd_wall_node[i]].p/(0.4*FLOW[sd_wall_node[i]].rho);	
				
				
				FLOW[node[sd_wall_node[i]].n_n[l]].u = (-1.0)*FLOW[singular[sd_wall_node[i]].n_n[k]].u;
				FLOW[node[sd_wall_node[i]].n_n[l]].v = (-1.0)*FLOW[singular[sd_wall_node[i]].n_n[k]].v;
				FLOW[node[sd_wall_node[i]].n_n[l]].w = (-1.0)*FLOW[singular[sd_wall_node[i]].n_n[k]].w;
				FLOW[node[sd_wall_node[i]].n_n[l]].p = FLOW[singular[sd_wall_node[i]].n_n[k]].p;
				FLOW[node[sd_wall_node[i]].n_n[l]].rho = FLOW[singular[sd_wall_node[i]].n_n[k]].rho;
				FLOW[node[sd_wall_node[i]].n_n[l]].t = FLOW[singular[sd_wall_node[i]].n_n[k]].t;	
				FLOW[node[sd_wall_node[i]].n_n[l]].a = FLOW[sd_wall_node[i]].a;
				FLOW[node[sd_wall_node[i]].n_n[l]].e = FLOW[singular[sd_wall_node[i]].n_n[k]].e;
				//mFLOW[node[sd_wall_node[i]].n_n[l]] = mFLOW[sd_wall_node[i]];
			
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].u;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].v;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].w;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].p = FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].p;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].rho = FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].rho;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].t = FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].t;	
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].a = FLOW[sd_wall_node[i]].a;
				FLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].e = FLOW[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].e;
				//mFLOW[node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = mFLOW[sd_wall_node[i]];
	
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].u = (-1.0)*FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].u;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].v = (-1.0)*FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].v;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].w = (-1.0)*FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].w;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].p = FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].p;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].rho = FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].rho;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].t = FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].t;	
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].a = FLOW[sd_wall_node[i]].a;
				FLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]].e = FLOW[node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]].e;
				//mFLOW[node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mFLOW[sd_wall_node[i]];	
				
			}
		}	
		
	}

}

