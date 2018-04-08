# include <iostream>
# include <fstream>
# include <string.h>
# include <algorithm>
# include <cmath>
#include <fstream>

#include <stdlib.h>
#include <stdio.h>






//using namespace std;






void navierstokes_semi_explicit(double Re, double conver_crit, double L1, double L2, int imax, int jmax, int cn, int cs, int cw, int ce, double u_vel_north, double u_vel_south, double u_vel_west, double u_vel_east, double v_vel_north, double v_vel_south, double v_vel_west, double v_vel_east, int pn, int ps, int pw, int pe, double p_north, double p_south, double p_west, double p_east)
{
	//inputs
	double rho = 1;
	double U = std::max(std::max(std::max(u_vel_north, u_vel_south), u_vel_east), u_vel_west);
	double V = std::max(std::max(std::max(v_vel_north, v_vel_south), v_vel_east), v_vel_west);
	double L = 1;
	double mu = 1 / Re;


	int i = 0, j = 0;

	double Dx = L1 / (imax - 2);
	double Dy = L2 / (jmax - 2);
	double DV = Dx*Dy;



	// Coordinates of face centers
	double *x = new double[imax];
	x[0] = 0;
	double C1 = (L1 / (imax - 2));
	for (i = 1; i < imax - 1; i++)
		x[i] = C1*i;
	double *y = new double[jmax];
	y[0] = 0;
	C1 = (L2 / (jmax - 2));
	for (j = 1; j < jmax - 1; j++)
		y[j] = C1*j;



	// Coordinates of Cell Centers; NEEDED FOR PLOTTING
	double *xc = new double[imax];
	xc[0] = 0;
	xc[imax - 1] = L1;
	for (i = 1; i < imax - 1; i++)
		xc[i] = (x[i] + x[i - 1]) / 2.0;
	double *yc = new double[jmax];
	yc[0] = 0;
	yc[jmax - 1] = L2;
	for (j = 1; j < jmax - 1; j++)
		yc[j] = (y[j] + y[j - 1]) / 2.0;

	//time step

	double Dt_advec = 1 / ((std::abs(U) / Dx) + (std::abs(V) / Dy));
	double Dt_cond = (1 / (2 * mu)) / ((1 / (Dx*Dx)) + (1 / (Dy*Dy)));
	double Dt = 0.5*std::min(Dt_advec, Dt_cond);


	double ** u;
	u = new double*[jmax];
	for (j = 0; j < jmax; j++)
		u[j] = new double[imax];
	double ** v;
	v = new double*[jmax];
	for (j = 0; j < jmax; j++)
		v[j] = new double[imax];
	double ** u_old;
	u_old = new double*[jmax];
	for (j = 0; j < jmax; j++)
		u_old[j] = new double[imax];
	double ** v_old;
	v_old = new double*[jmax];
	for (j = 0; j < jmax; j++)
		v_old[j] = new double[imax];
	double ** u_star;
	u_star = new double*[jmax];
	for (j = 0; j < jmax; j++)
		u_star[j] = new double[imax];
	double ** v_star;
	v_star = new double*[jmax];
	for (j = 0; j < jmax; j++)
		v_star[j] = new double[imax];




	double ** Au;
	Au = new double*[jmax];
	for (j = 0; j < jmax; j++)
		Au[j] = new double[imax];
	double ** Av;
	Av = new double*[jmax];
	for (j = 0; j < jmax; j++)
		Av[j] = new double[imax];
	double ** Du;
	Du = new double*[jmax];
	for (j = 0; j < jmax; j++)
		Du[j] = new double[imax];
	double ** Dv;
	Dv = new double*[jmax];
	for (j = 0; j < jmax; j++)
		Dv[j] = new double[imax];
	double ** Su;
	Su = new double*[jmax];
	for (j = 0; j < jmax; j++)
		Su[j] = new double[imax];
	double ** Sv;
	Sv = new double*[jmax];
	for (j = 0; j < jmax; j++)
		Sv[j] = new double[imax];
	double ** Sm_star;
	Sm_star = new double*[jmax];
	for (j = 0; j < jmax; j++)
		Sm_star[j] = new double[imax];
	double ** Sm_prime;
	Sm_prime = new double*[jmax];
	for (j = 0; j < jmax; j++)
		Sm_prime[j] = new double[imax];


	double ** p;
	p = new double*[jmax];
	for (j = 0; j < jmax; j++)
		p[j] = new double[imax];
	double ** p_star;
	p_star = new double*[jmax];
	for (j = 0; j < jmax; j++)
		p_star[j] = new double[imax];
	double ** p_prime;
	p_prime = new double*[jmax];
	for (j = 0; j < jmax; j++)
		p_prime[j] = new double[imax];
	p[j] = new double[imax];



	double ** mux;
	mux = new double*[jmax];
	for (j = 0; j < jmax; j++)
		mux[j] = new double[imax];
	double ** muy;
	muy = new double*[jmax];
	for (j = 0; j < jmax; j++)
		muy[j] = new double[imax];
	double ** mvx;
	mvx = new double*[jmax];
	for (j = 0; j < jmax; j++)
		mvx[j] = new double[imax];
	double ** mvy;
	mvy = new double*[jmax];
	for (j = 0; j < jmax; j++)
		mvy[j] = new double[imax];

	double ** aux;
	aux = new double*[jmax];
	for (j = 0; j < jmax; j++)
		aux[j] = new double[imax];
	double ** auy;
	auy = new double*[jmax];
	for (j = 0; j < jmax; j++)
		auy[j] = new double[imax];
	double ** avx;
	avx = new double*[jmax];
	for (j = 0; j < jmax; j++)
		avx[j] = new double[imax];
	double ** avy;
	avy = new double*[jmax];
	for (j = 0; j < jmax; j++)
		avy[j] = new double[imax];

	double ** dux;
	dux = new double*[jmax];
	for (j = 0; j < jmax; j++)
		dux[j] = new double[imax];
	double ** duy;
	duy = new double*[jmax];
	for (j = 0; j < jmax; j++)
		duy[j] = new double[imax];
	double ** dvx;
	dvx = new double*[jmax];
	for (j = 0; j < jmax; j++)
		dvx[j] = new double[imax];
	double ** dvy;
	dvy = new double*[jmax];
	for (j = 0; j < jmax; j++)
		dvy[j] = new double[imax];
	double ** my;
	my = new double*[jmax];
	for (j = 0; j < jmax; j++)
		my[j] = new double[imax];
	double ** mx;
	mx = new double*[jmax];
	for (j = 0; j < jmax; j++)
		mx[j] = new double[imax];
	double **my_star;
	my_star = new double*[jmax];
	for (j = 0; j < jmax; j++)
		my_star[j] = new double[imax];
	double ** mx_star;
	mx_star = new double*[jmax];
	for (j = 0; j < jmax; j++)
		mx_star[j] = new double[imax];
	double ** my_prime;
	my_prime = new double*[jmax];
	for (j = 0; j < jmax; j++)
		my_prime[j] = new double[imax];
	double ** mx_prime;
	mx_prime = new double*[jmax];
	for (j = 0; j < jmax; j++)
		mx_prime[j] = new double[imax];




	for (j = 0; j < jmax; j++)
		for (i = 0; i < imax; i++)

		{
			u[j][i] = 0;
			v[j][i] = 0;
			u_old[j][i] = 0;
			v_old[j][i] = 0;
			u_star[j][i] = 0;
			v_star[j][i] = 0;
			Au[j][i] = 0;
			Av[j][i] = 0;
			Du[j][i] = 0;
			Dv[j][i] = 0;
			Su[j][i] = 0;
			Sv[j][i] = 0;
			Sm_star[j][i] = 0;
			Sm_prime[j][i] = 0;
			mux[j][i] = 0;
			muy[j][i] = 0;
			mvx[j][i] = 0;
			mvy[j][i] = 0;
			aux[j][i] = 0;
			auy[j][i] = 0;
			avx[j][i] = 0;
			avy[j][i] = 0;
			dux[j][i] = 0;
			duy[j][i] = 0;
			dvx[j][i] = 0;
			dvy[j][i] = 0;
			p[j][i] = 0;
			p_star[j][i] = 0;
			p_prime[j][i] = 0;
			mx[j][i] = 0;
			mx_star[j][i] = 0;
			mx_prime[j][i] = 0;
			my[j][i] = 0;
			my_prime[j][i] = 0;
			my_star[j][i] = 0;
		}




	// velocity and pressure correction BC


	for (j = 0; j < jmax; j++)
	{
		u[j][0] = u_vel_west;
		u[j][imax - 2] = u_vel_east;
		v[j][0] = v_vel_west;
		v[j][imax - 1] = v_vel_east;
		p[j][0] = p_west;
		p[j][imax - 1] = p_east;
		p_star[j][0] = p_west;
		p_star[j][imax - 1] = p_east;
		p_prime[j][0] = p_prime[j][1];
		p_prime[j][imax - 1] = p_prime[j][imax - 2];
		u_star[j][0] = u_vel_west;
		u_star[j][imax - 2] = u_vel_east;
		v_star[j][0] = v_vel_west;
		v_star[j][imax - 1] = v_vel_east;
		u_old[j][0] = u_vel_west;
		u_old[j][imax - 2] = u_vel_east;
		v_old[j][0] = v_vel_west;
		v_old[j][imax - 1] = v_vel_east;

	}



	for (i = 0; i < imax; i++)
	{
		u[0][i] = u_vel_south;
		u[jmax - 1][i] = u_vel_north;
		v[0][i] = v_vel_south;
		v[jmax - 2][i] = v_vel_north;
		p[0][i] = p_south;
		p[jmax - 1][i] = p_north;
		p_star[0][i] = p_south;
		p_star[jmax - 1][i] = p_north;
		p_prime[jmax - 1][i] = p_prime[jmax - 2][i];
		p_prime[0][i] = p_prime[1][i];
		u_star[0][i] = u_vel_south;
		u_star[jmax - 1][i] = u_vel_north;
		v_star[0][i] = v_vel_south;
		v_star[jmax - 2][i] = v_vel_north;
		u_old[0][i] = u_vel_south;
		u_old[jmax - 1][i] = u_vel_north;
		v_old[0][i] = v_vel_south;
		v_old[jmax - 2][i] = v_vel_north;



	}

	for (j = 0; j < jmax; j++)
	{
		if (pw == 1)
			p[j][0] = p[j][1] - ((Dx / 2)*p_west);
		if (pe == 1)
			p[j][imax - 1] = p[j][imax - 2] + ((Dx / 2)*p_east);

	}



	for (i = 0; i < imax; i++)
	{
		if (pn == 1)
			p[jmax - 1][i] = p[jmax - 2][i] + ((Dy / 2)*p_north);
		if (ps == 1)
			p[0][i] = p[1][i] - ((Dy / 2)*p_south);

	}

	double Max = 0;
	double unsteadiness = 1;
	int n = 0;

	while (unsteadiness > conver_crit)
	{
		n = n + 1;
		for (j = 0; j < jmax; j++)

		{
			for (i = 0; i < imax; i++)
			{
				u_old[j][i] = u[j][i];
				v_old[j][i] = v[j][i];
			}
		}
		//u velocity and mass flux in x direction
		for (j = 1; j < jmax - 1; j++)
		{
			for (i = 0; i < imax - 2; i++)
			{
				mux[j][i] = rho*(u[j][i] + u[j][i + 1]) / 2;
				aux[j][i] = std::max(mux[j][i], 0.0)*u[j][i] + std::min(mux[j][i], 0.0)*u[j][i + 1];
				dux[j][i] = mu*(u[j][i + 1] - u[j][i]) / Dx;
			}

		}
		for (j = 0; j < jmax - 1; j++)
		{
			for (i = 1; i < imax - 2; i++)
			{
				muy[j][i] = rho*(v[j][i] + v[j][i + 1]) / 2;
				auy[j][i] = std::max(muy[j][i], 0.0)*u[j][i] + std::min(muy[j][i], 0.0)*u[j + 1][i];
				if (j == 0 || j == jmax - 2)
					duy[j][i] = mu*(u[j + 1][i] - u[j][i]) / (Dy / 2);
				else
					duy[j][i] = mu*(u[j + 1][i] - u[j][i]) / Dy;
			}
		}
		for (j = 1; j < jmax - 1; j++)
		{
			for (i = 1; i < imax - 2; i++)
			{

				Au[j][i] = (aux[j][i] - aux[j][i - 1])*Dy + (auy[j][i] - auy[j - 1][i])*Dx;
				Du[j][i] = (dux[j][i] - dux[j][i - 1])*Dy + (duy[j][i] - duy[j - 1][i])*Dx;
				Su[j][i] = (p[j][i] - p[j][i + 1])*Dy;
				u_star[j][i] = u[j][i] + (Dt / (rho*DV))*(Du[j][i] - Au[j][i] + Su[j][i]);
				mx_star[j][i] = rho*u_star[j][i];
			}
		}


		//v velocity and mass flux in y direction
		for (j = 1; j < jmax - 2; j++)
		{
			for (i = 0; i < imax - 1; i++)
			{
				mvx[j][i] = rho*(u[j][i] + u[j + 1][i]) / 2;
				avx[j][i] = std::max(mvx[j][i], 0.0)*v[j][i] + std::min(mvx[j][i], 0.0)*v[j][i + 1];
				if (i == 0 || i == imax - 2)
					dvx[j][i] = mu*(v[j][i + 1] - v[j][i]) / (Dx / 2);
				else
					dvx[j][i] = mu*(v[j][i + 1] - v[j][i]) / Dx;
			}
		}


		for (j = 0; j < jmax - 2; j++)
		{
			for (i = 1; i < imax - 1; i++)
			{
				mvy[j][i] = rho*(v[j][i] + v[j + 1][i]) / 2;
				avy[j][i] = std::max(mvy[j][i], 0.0)*v[j][i] + std::min(mvy[j][i], 0.0)*v[j + 1][i];
				dvy[j][i] = mu*(v[j + 1][i] - v[j][i]) / Dy;

			}
		}

		for (j = 1; j < jmax - 2; j++)
		{
			for (i = 1; i < imax - 1; i++)
			{


				Av[j][i] = (avx[j][i] - avx[j][i - 1])*Dy + (avy[j][i] - avy[j - 1][i])*Dx;
				Dv[j][i] = (dvx[j][i] - dvx[j][i - 1])*Dy + (dvy[j][i] - dvy[j - 1][i])*Dx;
				Sv[j][i] = (p[j][i] - p[j + 1][i])*Dx;
				v_star[j][i] = v[j][i] + (Dt / (rho*DV))*(Dv[j][i] - Av[j][i] + Sv[j][i]);
				my_star[j][i] = rho*v_star[j][i];
			}
		}
		for (j = 1; j < jmax - 1; j++)
		{
			if (cw == 1)
			{
				u[j][0] = u_star[j][1];
				v[j][0] = v_star[j][1];
			}
			if (ce == 1)
			{
				u[j][imax - 2] = u_star[j][imax - 3];
				v[j][imax - 1] = v_star[j][imax - 2];
			}
			mx_star[j][0] = rho*u[j][0];
			mx_star[j][imax - 2] = rho*u[j][imax - 2];
		}
		for (i = 1; i < imax - 1; i++)
		{
			if (cs == 1)
			{
				v[0][i] = v_star[1][i];
				u[0][i] = u_star[1][j];
			}
			if (cn == 1)
			{
				v[jmax - 2][i] = v_star[jmax - 3][i];
				u[jmax - 1][i] = u_star[jmax - 2][i];
			}
			my_star[0][i] = rho*v[0][i];
			my_star[jmax - 2][i] = rho*v[jmax - 2][i];
		}
		//mass source prediction
		for (j = 1; j < jmax - 1; j++)
		{
			for (i = 1; i < imax - 1; i++)
			{
				Sm_star[j][i] = (mx_star[j][i] - mx_star[j][i - 1])*Dy + (my_star[j][i] - my_star[j - 1][i])*Dx;
			}
		}

		//corrector step

		Max = 0;
		for (j = 0; j < jmax; j++)
		{
			for (i = 0; i < imax; i++)
			{
				p_prime[j][i] = 0;
				if (Max < std::abs(Sm_star[j][i]))
					Max = std::abs(Sm_star[j][i]);
			}
		}

		while (Max > 0.00000001)
		{
			double ap = Dt;
			for (j = 1; j < jmax - 1; j++)

			{
				for (i = 1; i < imax - 1; i++)
				{
					double delmx = (i == 1) ? (Dx / 2) : Dx;
					double delmy = (j == 1) ? (Dy / 2) : Dy;
					ap = Dt*((1 / Dx) + (1 / delmx))*Dy + ((1 / Dy) + (1 / delmy))*Dx;
					Sm_prime[j][i] = (-Dt*(((p_prime[j][i + 1] - p_prime[j][i]) / Dx) - ((p_prime[j][i] - p_prime[j][i - 1]) / delmx))*Dy) + (-Dt*(((p_prime[j + 1][i] - p_prime[j][i]) / Dy) - ((p_prime[j][i] - p_prime[j - 1][i]) / delmy))*Dx);
					p_prime[j][i] -= (1 / ap)*(Sm_star[j][i] + Sm_prime[j][i]);
				}
			}
			//bc
			for (j = 0; j < jmax; j++)
			{
				if (pw == 1)
					p_prime[j][0] = p_prime[j][1] - ((Dx / 2)*p_west);
				if (pe == 1)
					p_prime[j][imax - 1] = p_prime[j][imax - 2] + ((Dx / 2)*p_east);

			}



			for (i = 0; i < imax; i++)
			{
				if (pn == 1)
					p_prime[jmax - 1][i] = p_prime[jmax - 2][i] + ((Dy / 2)*p_north);
				if (ps == 1)
					p_prime[0][i] = p_prime[1][i] - ((Dy / 2)*p_south);

			}



			for (i = 0; i < imax - 1; i++)
			{
				for (j = 1; j < jmax - 1; j++)
				{
					double delx = (i == 0) || (i == imax - 2) ? (Dx / 2) : Dx;
					mx_prime[j][i] = (-Dt / delx)*(p_prime[j][i + 1] - p_prime[j][i]);
					mx_star[j][i] += mx_prime[j][i];
					u_star[j][i] = mx_star[j][i] / rho;
				}
			}
			for (j = 0; j < jmax - 1; j++)

			{
				for (i = 1; i < imax - 1; i++)
				{
					double dely = (j == 0) || (j == jmax - 2) ? (Dy / 2) : Dy;
					my_prime[j][i] = (-Dt / dely)*(p_prime[j + 1][i] - p_prime[j][i]);
					my_star[j][i] += my_prime[j][i];
					v_star[j][i] = my_star[j][i] / rho;
				}
			}
			for (j = 1; j < jmax - 1; j++)
			{
				if (cw == 1)
				{
					u[j][0] = u_star[j][1];
					v[j][0] = v_star[j][1];
				}
				if (ce == 1)
				{
					u[j][imax - 2] = u_star[j][imax - 3];
					v[j][imax - 1] = v_star[j][imax - 2];
				}
				mx_star[j][0] = rho*u[j][0];
				mx_star[j][imax - 2] = rho*u[j][imax - 2];
			}
			for (i = 1; i < imax - 1; i++)
			{
				if (cs == 1)
				{
					v[0][i] = v_star[1][i];
					u[0][i] = u_star[1][j];
				}
				if (cn == 1)
				{
					v[jmax - 2][i] = v_star[jmax - 3][i];
					u[jmax - 1][i] = u_star[jmax - 2][i];
				}
				my_star[0][i] = rho*v[0][i];
				my_star[jmax - 2][i] = rho*v[jmax - 2][i];
			}
			for (j = 1; j < jmax - 1; j++)

			{
				for (i = 1; i < imax - 1; i++)
				{
					Sm_star[j][i] = (mx_star[j][i] - mx_star[j][i - 1])*Dy + (my_star[j][i] - my_star[j - 1][i])*Dx;

				}
			}
			for (j = 0; j < jmax; j++)
			{
				for (i = 0; i < imax; i++)
				{

					p_star[j][i] += p_prime[j][i];
				}
			}
			Max = 0;
			for (j = 0; j < jmax; j++)
			{
				for (i = 0; i < imax; i++)
				{

					if (Max < std::abs(Sm_star[j][i]))
						Max = std::abs(Sm_star[j][i]);

				}
			}
		}

		for (j = 1; j < jmax - 1; j++)
			for (i = 1; i < imax - 2; i++)
				u[j][i] = u_star[j][i];
		for (j = 1; j < jmax - 2; j++)
			for (i = 1; i < imax - 1; i++)
				v[j][i] = v_star[j][i];
		unsteadiness = 0;
		for (j = 0; j < jmax; j++)
		{
			for (i = 0; i < imax; i++)
			{
				p[j][i] = p_star[j][i];
				//u[j][i] = u_star[j][i];
				//v[j][i] = v_star[j][i];
				mx[j][i] = mx_star[j][i];
				my[j][i] = my_star[j][i];
				if (unsteadiness < std::abs((u[j][i] - u_old[j][i]) / Dt))
					unsteadiness = std::abs((u[j][i] - u_old[j][i]) / Dt);
				if (unsteadiness < std::abs((v[j][i] - v_old[j][i]) / Dt))
					unsteadiness = std::abs((v[j][i] - v_old[j][i]) / Dt);

			}
		}
		unsteadiness = unsteadiness;
		//std::cout << '\n' << "time step " << n << " unsteadiness " << unsteadiness << '\n';
		printf("\n time step %d unsteadiness %f \n ", n, unsteadiness);

	}
	//graphical output





	//Vector plot

	double ** Vel1;
	Vel1 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		Vel1[j] = new double[imax];
	double ** Vel2;
	Vel2 = new double*[jmax];
	for (j = 0; j < jmax; j++)
		Vel2[j] = new double[imax];
	for (j = 0; j < jmax; j++)
	{
		for (i = 0; i < imax; i++)
		{
			if (i == 0)
			{
				Vel1[j][i] = u_vel_west;
				Vel2[j][i] = v_vel_west;
			}
			else
				if (i == imax - 1)
				{
					Vel1[j][i] = u_vel_east;
					Vel2[j][i] = v_vel_east;
				}
				else
					if (j == 0)
					{
						Vel1[j][i] = u_vel_south;
						Vel2[j][i] = v_vel_south;
					}
					else
						if (j == jmax - 1)
						{
							Vel1[j][i] = u_vel_north;
							Vel2[j][i] = v_vel_north;
						}
						else
						{   //0u0 1 1u 2 2u 3 3u 4 4u 5 5u6
							Vel1[j][i] = (u[j][i] + u[j][i - 1]) / 2;
							Vel2[j][i] = (v[j][i] + v[j - 1][i]) / 2;
						}

		}
	}
	std::ofstream out6;
	out6.open("vector_data.vtk");
	out6 << "# vtk DataFile Version 3.0" << '\n' << "vtk output" << '\n' << "ASCII" << '\n' << "DATASET STRUCTURED_GRID" << '\n';
	out6 << "DIMENSIONS " << jmax << " " << imax << " 1" << '\n';
	int no = imax*jmax;
	out6 << "POINTS " << no << " float" << '\n';
	for (i = 0; i < imax; i++)
	{
		for (j = 0; j < jmax; j++)
		{
			out6 << (float)xc[i] << " " << (float)yc[j] << " 0" << '\n';

		}
	}
	no = (imax - 1)*(jmax - 1);
	out6 << "CELL_DATA " << no << '\n';
	no = imax*jmax;
	out6 << "POINT_DATA " << no << '\n';
	out6 << "SCALARS Pressure float" << '\n' << "LOOKUP_TABLE default" << '\n';
	for (i = 0; i < imax; i++)
	{
		for (j = 0; j < jmax; j++)
		{
			out6 << (float)p[j][i] << '\n';

		}
	}
	out6 << "VECTORS Velocity float" << '\n';
	for (i = 0; i < imax; i++)
	{
		for (j = 0; j < jmax; j++)
		{
			out6 << Vel1[j][i] << '\n' << Vel2[j][i] << '\n' << "0" << '\n';

		}
	}

	out6.close();


}






int main()

{
	double u_vel_west =0, u_vel_east = 0, u_vel_south = 0, u_vel_north = 1;
	double v_vel_west = 0, v_vel_east = 0, v_vel_south = 0, v_vel_north = 0;
	double  p_north = 0, p_south = 0, p_west = 0, p_east = 0;
	int cn = 0; int ce = 0; int cw = 0; int cs = 0;
	int pn = 1; int ps = 1; int pw = 1; int pe = 1;
	double L1 = 1, L2 = 1;
	double Re = 100;
	int imax = 12, jmax = 10;
	double conver_crit = 0.001;
	navierstokes_semi_explicit(Re,conver_crit,L1,L2,imax,jmax,cn,cs,cw,ce, u_vel_north, u_vel_south, u_vel_west, u_vel_east , v_vel_north, v_vel_south,v_vel_west , v_vel_east,pn,ps,pw,pe,p_north,p_south,p_west,p_east);
	


}
