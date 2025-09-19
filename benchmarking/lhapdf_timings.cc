#include <LHAPDF/LHAPDF.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include <string>

double zeta_of_y(double y, double a) {
    return y + a * (1.0 - std::exp(-y));
}

double y_of_zeta(double zeta, double a) {
    double y = zeta;
    if (a != 0.0) {
        for (int iter = 0; iter < 100; ++iter) {
            double x = std::exp(-y);
            double diff = zeta - y - a * (1.0 - x);
            if (std::abs(diff) < 1e-12) break;
            double deriv = -1.0 - a * x;
            y -= diff / deriv;
        }
    }
    return y;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: ./lhapdf_timings PDFsetName outputdir\n";
        return 1;
    }
    std::string pdfname = argv[1];
    int member = 0;
    std::string outputdir = argv[2];

    // Output file
    std::ofstream fout(outputdir + "/" + pdfname + ".lhapdf.dat");
    fout << std::setprecision(10) << std::scientific;

    // Grid parameters (set to match Fortran defaults, or parse from command line)
    int nxQ = 5000; // or your desired value
    double ymax = 12.0;
    double Qinit = std::sqrt(2.0);
    double Qmax = 1e4;
    double grid_a = 9.0, gridQ_a = 3.0;

    int nz = static_cast<int>(std::round(std::sqrt(4.0 * nxQ)));
    int nQ = static_cast<int>(std::round(std::sqrt(0.25 * nxQ))) - 1;

    double zmax = zeta_of_y(ymax, grid_a);
    double zQmax = zeta_of_y(std::log(Qmax / Qinit), gridQ_a);

    // LHAPDF setup
    auto pdf = LHAPDF::mkPDF(pdfname, member);

    // Header
    for(int i = 0; i < 15; i++)
        fout << "#\n";
    fout << "# y=ln1/x Q pdf(-5:5)\n";

    // Grid evaluation
    std::vector<double> pdfval(13);
    for (int iQ = 0; iQ <= nQ; ++iQ) {
        for (int iz = 1; iz <= nz; ++iz) {
            double zeta = iz * zmax / nz;
            double y = y_of_zeta(zeta, grid_a);
            double zQ = (iQ + zeta / zmax) * zQmax / nQ;
            double Q = std::max(Qinit, std::min(Qmax, Qinit * std::exp(y_of_zeta(zQ, gridQ_a))));
            double x = std::exp(-y);

            pdf->xfxQ(x, Q, pdfval);

            fout << std::setw(20) << y << std::setw(20) << Q;
            for (int k = -5; k <= 5; ++k)
                fout << std::setw(20) << pdfval[k + 6];
            fout << "\n";
        }
        fout << "\n";
    }

    // Interpolation timing
    int nrep_interp = 1000000; // or your desired value
    int ny = static_cast<int>(std::round(std::sqrt(4.0 * nrep_interp)));
    int nQ_interp = static_cast<int>(std::round(std::sqrt(0.25 * nrep_interp))) - 1;

    std::vector<double> yvals(ny), xvals(ny), Qvals(nQ_interp);
    for (int i = 0; i < ny; ++i) {
        yvals[i] = (i + 1) * ymax / ny;
        xvals[i] = std::exp(-yvals[i]);
    }
    for (int i = 0; i < nQ_interp; ++i)
        Qvals[i] = Qinit + (i + 1) * (Qmax - Qinit) / nQ_interp;

    auto time_interp_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nQ_interp; ++j) {
            pdf->xfxQ(xvals[i], Qvals[j], pdfval);
        }
    }
    auto time_end = std::chrono::high_resolution_clock::now();
    double interp_ns = std::chrono::duration<double>(time_end - time_interp_start).count() / ny / nQ_interp * 1e9;

    fout << "# Interpolation time = " << std::fixed << std::setprecision(5)
         << interp_ns << " ns, nrep_interp = " << nrep_interp << "\n";

    for(int i = 0; i < 10; i++)
        fout << "#\n";

        fout.close();
    return 0;
}