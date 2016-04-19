#include "Twin.h"

void integrator(double&, double&, double&, double&);
double rk4(double&, double&, double&);
double derivs(const double&, const double&);
void dvec(Twin&, const std::vector<double>&);

int main() {

	Twin t("Runge-Kutta Method");

	double y, xi, xf, dx, xout;

	while (true) {

		t.println("Enter the initial value of the dependent variable.");
		t.println("Example: 1.2");
		t.display();

		std::cin >> y;
		std::cin.ignore(1000, '\n');
		std::cin.clear();

		t.println("y:");
		t.println();
		t.println(y);
		t.println();
		t.println("Enter 'x' to save and continue.");
		t.display();

		char temp;
		std::cin >> temp;
		std::cin.ignore(1000, '\n');
		std::cin.clear();

		if (temp == 'x' || temp == 'X') break;

	}

	while (true) {

		t.println("Enter the initial value of the independent variable.");
		t.println("Example: 3.4");
		t.display();

		std::cin >> xi;
		std::cin.ignore(1000, '\n');
		std::cin.clear();

		t.println("xi:");
		t.println();
		t.println(xi);
		t.println();
		t.println("Enter 'x' to save and continue.");
		t.display();

		char temp;
		std::cin >> temp;
		std::cin.ignore(1000, '\n');
		std::cin.clear();

		if (temp == 'x' || temp == 'X') break;

	}

	while (true) {

		t.println("Enter the final value of the independent variable.");
		t.println("Example: 6.8");
		t.display();

		std::cin >> xf;
		std::cin.ignore(1000, '\n');
		std::cin.clear();

		t.println("xf:");
		t.println();
		t.println(xf);
		t.println();
		t.println("Enter 'x' to save and continue.");
		t.display();

		char temp;
		std::cin >> temp;
		std::cin.ignore(1000, '\n');
		std::cin.clear();

		if (temp == 'x' || temp == 'X') break;

	}

	while (true) {

		t.println("Enter the calculation step size.");
		t.println("Example: 1.0");
		t.display();

		std::cin >> dx;
		std::cin.ignore(1000, '\n');
		std::cin.clear();

		t.println("dx:");
		t.println();
		t.println(dx);
		t.println();
		t.println("Enter 'x' to save and continue.");
		t.display();

		char temp;
		std::cin >> temp;
		std::cin.ignore(1000, '\n');
		std::cin.clear();

		if (temp == 'x' || temp == 'X') break;

	}

	while (true) {

		t.println("Enter the output interval.");
		t.println("Example: 2.0");
		t.display();

		std::cin >> xout;
		std::cin.ignore(1000, '\n');
		std::cin.clear();

		t.println("xout:");
		t.println();
		t.println(xout);
		t.println();
		t.println("Enter 'x' to save and continue.");
		t.display();

		char temp;
		std::cin >> temp;
		std::cin.ignore(1000, '\n');
		std::cin.clear();

		if (temp == 'x' || temp == 'X') break;

	}

	double x = xi;
	std::vector<double> xpm = { x }, ypm = { y };

	while(x < xf) {

		double xend = x + xout;
		if (xend > xf) xend = xf;

		integrator(x, y, dx, xend);
		xpm.push_back(x);
		ypm.push_back(y);

	}

	t.println("Vector xpm:");
	dvec(t, xpm);
	t.display();

	std::cin.ignore(1000, '\n');
	std::cin.clear();

	t.println("Vector ypm:");
	dvec(t, ypm);
	t.display();


	std::cout << "Press enter to exit. . .";
	std::cin.ignore(1000, '\n');
	std::cin.clear();

	return EXIT_SUCCESS;
}

void integrator(double& x, double& y, double& h, double& xend) {

//	while (true) {

		if (xend - x < h) h = xend - x;
		y = rk4(x, y, h);
//		if (x >= xend) break;

//	}

}

double rk4(double& x, double& y, double& h) {

	double k1 = derivs(x, y);
	double ym = y + k1 * h / 2;
	double k2 = derivs(x + h / 2, ym);
	ym = y + k2 * h / 2;
	double k3 = derivs(x + h / 2, ym);
	double ye = y + k3 * h;
	double k4 = derivs(x + h, ye);
	double slope = (k1 + 2 * (k2 + k3) + k4) / 6;
	x += h;
	return y + slope * h;

}

double derivs(const double& x, const double& y) {

	return .026 * (1 - y / 12000)* y;

}

void dvec(Twin& t, const std::vector<double>& vec) {

	for (auto&& i : vec) {
		t.println("[ ");
		t.print(i);
		t.print(" ]");
	}

	t.println();
}