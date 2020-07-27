// Fit charges on a surface of spherical shells to reproduce the electrostatic
// potential due to solvent at points on another surface

#include <loos.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

typedef Math::Matrix<double, Math::RowMajor> rMatrix;
typedef boost::variate_generator<base_generator_type&, boost::uniform_real<> >
    RNG;

string fullHelpMessage(void) {
    string msg =
        "\nSYNOPSIS\n"
        "\tFit charges on a surface of spherical shells to reproduce the\n"
        "\telectrostatic potential due to solvent at points on another surface."
        "\n\nDESCRIPTION\n"
        "\tGiven an atom selection, this tool will generate a grid of points\n"
        "\tequidistant on a surface defined by a union of spherical shells "
        "centered\n\ton the coordinates of those atoms in the first frame of "
        "the trajectory.\n\tA grid of equidistant points on a spherical shell "
        "with user-specified\n\tradius and point density is generated using a "
        "using a Fibonacci\n\tspiral [1]. For each atom, these points are "
        "rotated to a random\n\torientation using quaternions [2,3] generated "
        "using a random point on\n\tthe surface of a 4D hypersphere [4] and "
        "then translated to the\n\tcoordinates of the atom. Then, grid points "
        "closer than the specified\n\tradius to another atom are removed.\n\n"
        "\tTwo such surfaces are generated at different radii. This tool will\n"
        "\tcompute the classical electrostatic potential (ESP) at each atom in "
        "the\n\tselection and at each grid point for the closer radius due to "
        "solvent\n\tmolecules for which at least one atom is closer than a "
        "cutoff distance\n\tbut no atoms are closer than the farther radius to "
        "any atom in the\n\tselection. Then, point charges are placed at grid "
        "points for the farther\n\tradius, and their magnitude is fit by "
        "linear regression to reproduce the\n\tsolvent ESP averaged over the "
        "over the trajectory. The difference\n\tbetween the ESP at the atoms "
        "in the selection due to these point charges\n\tand that due to the "
        "solvent informs on the quality of the fit.\n\n"
        "\tThe first line of the output is the command-line input. The second "
        "line\n\tcontains the number of grid points on the ESP surface, the "
        "number of fit point\n\tcharges, the sum of the fit point charges, and "
        "the sum of square residuals from\n\tthe charge fit. The third line is "
        "the ESP at the atoms in the selection due to\n\tsolvent averaged over "
        "the trajectory. The third line is the ESP at the\n\tcoordinates of "
        "the atoms in the selection in the first frame due to the fit\n\tpoint "
        "charges. Then, for each grid point for the closer radius, the\n\t"
        "output contains the coordinates, the ESP due to solvent averaged over "
        "\n\tthe trajectory, and the ESP due to the fit point charges. "
        "Finally, the output contains the coordinates and magnitude\n\t of "
        "each fit point charge.\n"
        "\nEXAMPLES\n"
        "\nREFERENCES\n"
        "\t[1] Swinbank R & Purser RJ (2006) Q J R Meteorol Soc 132, 1769\n"
        "\t[2] Vesely FJ (1982) J Comp Phys 47, 291\n"
        "\t[3] Frenkel D & Smit B (2002) Understanding Molecular Simulation, "
        "2e, 49\n"
        "\t[4] Marsaglia G (1972) Ann Math Stat 43, 645\n";
    return msg;
}

class ToolOptions : public opts::OptionsPackage {

public:
    ToolOptions() {}

    void addGeneric(po::options_description &o) {
        o.add_options()
            ("grid_atoms,a", po::value<string>(&atom_selection)->default_value(
                "!(resname == 'WAT' || resname =~ '[+-]$')"), "Selection "
                "string for atoms used to generate the grid surface.")
            ("solvent,w", po::value<string>(&solvent_selection)->default_value(
                "resname == 'WAT' || resname =~ '[+-]$'"), "Selection string "
                "for solvent atoms.")
            ("cutoff,C", po::value<double>(&cutoff)->default_value(1000.0),
                "Distance in angstroms from any atom in the selection beyond "
                "which solvent molecules will be excluded.")
            ("fit_radius,R", po::value<string>(&fit_radius_str)->default_value(
                "10.0"), "Radius in angstroms of the spherical shells centered "
                "on each atom whose union will produce the surface containing "
                "the fit point charges.")
            ("fit_area,A", po::value<string>(&fit_area_str)->default_value(
                "10.0"), "Approximate area per point in square angstroms for "
                "the surface containing the fit point charges.")
            ("esp_radius,S", po::value<double>(&esp_radius)->default_value(3.0),
                "Radius in angstroms of the spherical shells centered on each "
                "atom whose union will produce the surface on which the ESP "
                "will be evaluated.")
            ("esp_area,B", po::value<double>(&esp_area)->default_value(1.0),
                "Approximate area per point in square angstroms for the "
                "surface on which the ESP will be evaluated.")
            ("weight,W", po::value<double>(&weight)->default_value(
                0.01 / 18.2223 / 18.2223), "Weight of harmonic restraint that "
                "restrains fit charges to zero.")
            ("close_solvent,c", po::value<bool>(&close_solvent)->default_value(
                false), "Include solvent closer than the closest fitting "
                "surface but farther than the ESP surface in the fitting "
                "target for the closest fitting surface.")
            ("esp_file,f", po::value<string>(&esp_file)->default_value(""),
                "Name of the output file to which ESP will be written.")
            ("seed", po::value<uint>(&seed)->default_value(0), "Seed for random"
                " number generator. 0 indicates to use the current time.")
            ("solvent_size", po::value<uint>(&solvent_size)->default_value(3),
                "Largest number of atoms in a solvent molecule.");
    }

    string print() const {
        ostringstream oss;
        oss << boost::format("grid_atoms=%s, solvent=%s, cutoff=%f, "
            "fit_radius=%s, fit_area=%s, esp_radius=%f, esp_area=%f, "
            "weight=%f, close_solvent=%d, esp_file=%s, seed=%d, "
            "solvent_size=%d") % atom_selection % solvent_selection % cutoff %
            fit_radius_str % fit_area_str % esp_radius % esp_area % weight %
            close_solvent % esp_file % seed % solvent_size;
        return (oss.str());
    }

    string atom_selection;
    string solvent_selection;
    double cutoff;
    string fit_radius_str;
    string fit_area_str;
    double esp_radius;
    double esp_area;
    double weight;
    bool close_solvent;
    string esp_file;
    uint seed;
    uint solvent_size;

}; // ToolOptions

// Get a vector of doubles from a string of numbers separated by spaces
//     input: string containing numbers separated by spaces
vector<double> stringToDoubles(const string& input) {

    istringstream input_stream(input);
    return vector<double>{istream_iterator<double>(input_stream),
        istream_iterator<double>()};

}

// Generate equidistant points on a spherical shell using a Fibonacci spiral
// Swinbank R & Purser RJ (2006) Q J R Meteorol Soc 132, 1769
//               N: number of points
//     sphere_grid: flattened vector of grid coordinates of size 3 * N
void generateEquidistantSphereGrid(const uint N, vector<double>& sphere_grid) {

    double z; // z = cos(theta), where theta is the polar angle from the +z axis
    double sin_theta; // sin(theta)
    double phi; // azimuthal angle from the +x axis

    // Interval between z coordinates of points
    const double z_interval = 2.0 / N;
    // Interval between phi coordinates of points. Equal to
    // 2 * pi * (1 - 1 / golden_ratio), where golden_ratio is (1 - sqrt(5)) / 2
    const double golden_angle = M_PI * (3 - sqrt(5));

    // Generate the grid
    z = 1.0 + z_interval * 0.5;
    phi = -golden_angle * 0.5;
    for (uint i = 0; i < N; ++i) {
        // Even distribution of z values from +1 to -1
        // Equivalent to z = 1 - 2 / N * (i + 0.5)
        z -= z_interval;
        // Rotate around the z axis by goldenAngle radians
        // Equivalent to phi = goldenAngle * (i + 0.5)
        phi += golden_angle;
        // z = cos(theta) and sin^2(theta) + cos^2(theta) = 1
        sin_theta = sqrt(1 - z * z);
        // Convert spherical coordinates to Cartesian coordinates
        sphere_grid.push_back(sin_theta * cos(phi));
        sphere_grid.push_back(sin_theta * sin(phi));
        sphere_grid.push_back(z);
    }

} // generateEquidistantSphereGrid()

// Generate a random rotation using a quaternion (q1, q2, q3, q0)
//      uniform_pm_one: function to produce uniform random numbers on (-1, 1)
//     rotation_matrix: random rotation matrix in three dimensions
void generateRandomRotation(RNG& uniform_pm_one,
    vector<double>& rotation_matrix) {

    double r1, r2; // Random numbers uniform on (-1, 1)
    double s1 = 2.0, s2 = 2.0; // Sums of squares of two random numbers
    double q1, q2, q3, q0; // Quaternion components
    // Intermediate quantities for construction of rotation matrix
    double q1_2, q2_2, q3_2, q0_2, q1_q2, q2_q3, q3_q1, q1_q0, q2_q0, q3_q0;

    // Generate a random point on a 4D hypersphere efficiently
    // Marsaglia G (1972) Ann Math Stat 43, 645
    while (s1 >= 1.0) {
        q1 = uniform_pm_one();
        q2 = uniform_pm_one();
        s1 = q1 * q1 + q2 * q2;
    }
    while (s2 >= 1.0) {
        r1 = uniform_pm_one();
        r2 = uniform_pm_one();
        s2 = r1 * r1 + r2 * r2;
    }
    q3 = r1 * sqrt((1.0 - s1) / s2);
    q0 = q3 * r2 / r1;

    // Construct the rotation matrix in terms of the quaternion
    // Quaternion components in terms of Euler angles (phi, theta, psi) are
    // q1 = sin(theta / 2) * cos((phi - psi) / 2)
    // q2 = sin(theta / 2) * sin((phi - psi) / 2)
    // q3 = cos(theta / 2) * sin((phi + psi) / 2)
    // q0 = cos(theta / 2) * cos((phi + psi) / 2)
    // Vesely FJ (1982) J Comp Phys 47, 291
    // Frenkel D & Smit B (2002) Understanding Molecular Simulation, 2e, 49
    q1_2 = q1 * q1;
    q2_2 = q2 * q2;
    q3_2 = q3 * q3;
    q0_2 = q0 * q0;
    q1_q2 = 2.0 * q1 * q2;
    q2_q3 = 2.0 * q2 * q3;
    q3_q1 = 2.0 * q3 * q1;
    q1_q0 = 2.0 * q1 * q0;
    q2_q0 = 2.0 * q2 * q0;
    q3_q0 = 2.0 * q3 * q0;
    rotation_matrix[0] = q1_2 - q2_2 - q3_2 + q0_2;
    rotation_matrix[1] = q1_q2 - q3_q0;
    rotation_matrix[2] = q3_q1 + q2_q0;
    rotation_matrix[3] = q1_q2 + q3_q0;
    rotation_matrix[4] = q2_2 - q3_2 - q1_2 + q0_2;
    rotation_matrix[5] = q2_q3 - q1_q0;
    rotation_matrix[6] = q3_q1 - q2_q0;
    rotation_matrix[7] = q2_q3 + q1_q0;
    rotation_matrix[8] = q3_2 - q1_2 - q2_2 + q0_2;

} // generateRandomRotation()

// Generate grid of equidistant points on a surface defined by a union of
// spherical shells centered on atoms in a selection
//       uniform_pm_one: function to produce uniform random numbers on (-1, 1)
//     grid_atom_coords: coordinates used to center the spherical shells
//               N_atom: number of atoms in the selection
//               radius: radius of the spherical shells
//              N_shell: number of points per spherical shell
//         surface_grid: flattened vector of the surface grid coordinates
void generateSurfaceGrid(RNG& uniform_pm_one,
    const vector<GCoord>& grid_atom_coords, const uint N_atom,
    const double radius, const uint N_shell, vector<GCoord>& surface_grid) {

    vector<double> sphere_grid; // Grid on a spherical shell
    vector<double> rotation_matrix(9); // Rotation matrix
    double atom_x, atom_y, atom_z; // Coordinates of an atom
    double sphere_x, sphere_y, sphere_z; // Coordinates of point in sphere_grid
    double grid_x, grid_y, grid_z; // Coordinates of point in surface_grid
    double dx, dy, dz; // Components of a vector from an atom to a grid point
    const double R_2 = radius * radius; // Radius squared
    bool record_grid_point; // Whether to include a grid point in surface_grid
    uint jj; // Index into sphere_grid

    // Generate grid of equidistant points on a spherical shell
    sphere_grid.reserve(3 * N_shell);
    generateEquidistantSphereGrid(N_shell, sphere_grid);

    // For each atom, rotate sphere_grid to a random orientation, expand to
    // radius, and translate to the coordinates of the atom. Then, delete grid
    // points within a distance of R from any other atom.
    for (uint i = 0; i < N_atom; ++i) {

        atom_x = grid_atom_coords[i][0];
        atom_y = grid_atom_coords[i][1];
        atom_z = grid_atom_coords[i][2];

        // Generate random rotation matrix in three dimensions
        generateRandomRotation(uniform_pm_one, rotation_matrix);

        for (uint j = 0; j < N_shell; ++j) {

            jj = 3 * j;
            sphere_x = sphere_grid[jj];
            sphere_y = sphere_grid[jj + 1];
            sphere_z = sphere_grid[jj + 2];

            // Apply random rotation, expand to radius, and translate
            // grid = (rotation_matrix @ sphere) * radius + atom
            grid_x = (rotation_matrix[0] * sphere_x + rotation_matrix[1] * 
                sphere_y + rotation_matrix[2] * sphere_z) * radius + atom_x;
            grid_y = (rotation_matrix[3] * sphere_x + rotation_matrix[4] * 
                sphere_y + rotation_matrix[5] * sphere_z) * radius + atom_y;
            grid_z = (rotation_matrix[6] * sphere_x + rotation_matrix[7] * 
                sphere_y + rotation_matrix[8] * sphere_z) * radius + atom_z;

            // Check whether grid point is closer than radius to any other atom
            record_grid_point = true;
            for (uint k = 0; k < N_atom; ++k) {

                // Don't check current atom
                if (k == i) continue;

                // Get components of vector from atom to grid point
                dx = grid_x - grid_atom_coords[k][0];
                dy = grid_y - grid_atom_coords[k][1];
                dz = grid_z - grid_atom_coords[k][2];

                // Check whether grid point is closer than R
                if (dx * dx + dy * dy + dz * dz < R_2) {
                    record_grid_point = false;
                    break;
                }

            } // Loop over atoms, k

            // Record grid point if it is farther than radius from other atoms
            if (record_grid_point) {
                surface_grid.push_back(GCoord(grid_x, grid_y, grid_z));
            }

        } // Loop over grid points, j

    } // Loop over atoms, i

} // generateSurfaceGrid()

// Calculate electrostatic potential (ESP) at points in esp_grid and at atoms in
// grid_atoms due to solvent molecules in which no atom is within a low cutoff
// but at least one atom is within a high cutoff of any atom in grid_atoms. The
// low cutoff is given by fit_radius[i], and the high cutoff is given by
// fit_radius[i - 1] for i > 1 and by high_cutoff for i == 0.
//     solvent_molecules: molecules in the solvent that generate the ESP
//            N_solvent: number of solvent molecules
//         solvent_size: largest number of atoms in a solvent molecule
//           grid_atoms: atoms at which the ESP will be evaluated
//               N_atom: number of atoms in grid_atoms
//             esp_grid: grid points at which the ESP will be evaluated
//         N_esp_points: number of grid points
//     solvent_boundary: boundaries on solvent distances to grid atoms
//             N_region: number of solvent-containing regions
//          solvent_esp: ESP due to solvent molecules
//        close_solvent: include solvent closer than the closest surface
void calculateFrameESP(vector<AtomicGroup>& solvent_molecules,
    const uint N_solvent, const uint solvent_size,
    const AtomicGroup& grid_atoms, const uint N_atom,
    const vector<GCoord>& esp_grid, const uint N_esp_points,
    const vector<double>& solvent_boundary, const uint N_region,
    rMatrix& solvent_esp, const bool close_solvent) {

    double R2_min; // Closest squared distance between solvent and grid atoms
    double solvent_charge; // Charge of a solvent atom
    GCoord solvent_coords; // Coordinates of a solvent atom
    uint surface; // Index of fitting surface

    // Matrix of squared distances between solvent atoms and grid atoms
    rMatrix R2(solvent_size, N_atom);

    // Loop over solvent molecules
    for (uint i = 0; i < N_solvent; ++i) {

        // Compute squared distances between all solvent atoms and all grid
        // atoms and find the minimum of such squared distances
        R2_min = solvent_boundary[0];

        // Loop over solvent atoms
        for (uint j = 0; j < solvent_molecules[i].size(); ++j) {

            solvent_coords = solvent_molecules[i][j]->coords();

            // Loop over grid atoms
            for (uint k = 0; k < N_atom; ++k) {
                R2(j,k) = solvent_coords.distance2(grid_atoms[k]->coords());
                if (R2(j,k) < R2_min)
                    R2_min = R2(j,k);
            }
        }

        // Identify to which region this solvent molecule belongs

        // Exclude solvent molecules farther than the cutoff
        if (R2_min >= solvent_boundary[0])
            continue;

        // Loop over solvent boundaries
        for (uint j = 1; j <= N_region; ++j) {
            if (R2_min >= solvent_boundary[j]) {
                surface = j - 1;
                break;
            }
        }

        // Solvent molecule closer than the closest fitting surface 
        if (R2_min < solvent_boundary[N_region]) {
            if (close_solvent && R2_min >= solvent_boundary[N_region + 1])
                surface = N_region - 1; // Include in ESP for closest surface
            else
                continue; // Exclude
        }

        // Calculate ESP due to this solvent molecule

        // Loop over solvent atoms
        for (uint j = 0; j < solvent_molecules[i].size(); ++j) {

            solvent_charge = solvent_molecules[i][j]->charge();
            solvent_coords = solvent_molecules[i][j]->coords();

            // Loop over grid points in esp_grid
            for (uint k = 0; k < N_esp_points; ++k)
                solvent_esp(surface, k) += solvent_charge /
                    solvent_coords.distance(esp_grid[k]);

            // Loop over atoms in grid_atoms
            for (uint k = 0; k < N_atom; ++k)
                solvent_esp(surface, k + N_esp_points) += solvent_charge /
                    sqrt(R2(j,k));

        }

    } // Loop over solvent molecules, i

} // calculateFrameESP()

// Write the solvent ESP at esp_grid and at grid atoms to an output file
//        esp_file: name of output file
//           frame: number of frames processed so far
//           index: index of current frame
//     solvent_esp: ESP due to solvent molecules
void writeESP(const string& esp_file, const uint frame, const uint index,
    const rMatrix& solvent_esp) {

    ofstream esp_out(esp_file.c_str(), ios::app); // Stream for output file
    double esp_sum; // Sum of ESP for different regions

    esp_out << index;

    // Loop over columns, i.e. points in esp_grid and grid_atoms
    for (uint i = 0; i < solvent_esp.cols(); ++i) {

        esp_sum = 0;

        // Loop over rows, i.e. fitting surfaces
        for (uint j = 0; j < solvent_esp.rows(); ++j)
            esp_sum += solvent_esp(j, i);

        // Print esp due to solvent in all regions averaged over the number of
        // frames processed so far
        esp_out << " " << boost::format("%12.8f") % (esp_sum / frame);

    }

    esp_out << "\n";

} // writeESP()

// Fit the magnitude of points charges located at grid points in fit_grid to
// reproduce the solvent ESP at grid points in esp_grid using linear regression
// from the LAPACK routine DGELS
//      one_over_R: fitting matrix containing reciprocal distance
//      target_esp: fitting target containing solvent ESP
//           N_esp: number of grid points at which ESP was evaluated
//     fit_charges: magnitudes of point charges to reproduce solvent ESP
//           N_fit: number of point charges
//          weight: weight of harmonic restaint of fit charges to zero
void fitChargeGrid(const rMatrix& one_over_R, const vector<double>& target_esp,
    const uint N_esp, vector<double>& fit_charges, const uint N_fit,
    const double weight) {

    // Number of rows in fitting matrix
    uint N_row;
    if (weight == 0)
        N_row = N_esp;
    else
        N_row = N_esp + N_fit;

    // Arrays passed to Fortran routines must be in column-major order
    Math::Matrix<double, Math::ColMajor> target(N_row, 1); // Fitting target
    Math::Matrix<double, Math::ColMajor> matrix(N_row, N_fit); // Fitting matrix

    // Copy target_esp to target and copy oneOverR to matrix because they will
    // be overwritten by DGELS
    for (uint i = 0; i < N_esp; ++i) {
        target[i] = target_esp[i];
        for (uint j = 0; j < N_fit; ++j)
            matrix(i, j) = one_over_R(i, j);
    }

    // Add harmonic restraints to restrain magnitude of fit charges to zero
    if (weight > 0)
        for (uint i = 0; i < N_fit; ++i)
            matrix(i + N_esp, i) = weight;

    // Set up arguments to DGELS
    char TRANS = 'N';
    int M = N_row;
    int N = N_fit;
    int NRHS = 1;
    int LDA = M;
    int LDB = M;
    double PRE_WORK;
    int LWORK = -1;
    int INFO;

    // First call to DGELS to get value of LWORK
    dgels_(&TRANS, &M, &N, &NRHS, matrix.get(), &LDA, target.get(), &LDB,
        &PRE_WORK, &LWORK, &INFO);
    if (INFO != 0)
        throw(NumericalError("DGELS failed to estimate value of LWORK."), INFO);
    LWORK = (int) PRE_WORK;
    double* WORK = new double[LWORK+1];

    // Second call to DGELS to solve linear regression
    dgels_(&TRANS, &M, &N, &NRHS, matrix.get(), &LDA, target.get(), &LDB, WORK,
        &LWORK, &INFO);
    if (INFO != 0)
        throw(NumericalError("DGELS failed to solve linear regression."), INFO);

    // Now target contains the solution vector (magnitudes of fit charges) in
    // the first N_fit elements, and the sum of square residuals (RSS) is the
    // sum of the squares of elements N_fit + 1 to N_esp
    double sum_Q = 0;
    for (uint i = 0; i < N_fit; ++i) {
        sum_Q += target[i];
        fit_charges.push_back(target[i]);
    }
    double RSS = 0;
    for (uint i = N_fit; i < N_row; ++i)
        RSS += target[i] * target[i];
    cout << boost::format("# %4d %4d %12.8f %14.8e\n") % N_esp % N_fit % sum_Q %
        RSS;

} // fitChargeGrid()

int main(int argc, char *argv[]) {

    // Command-line input
    string header = invocationHeader(argc, argv);

    // Set up tool options
    opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage());
    opts::BasicSelection *sopts = new opts::BasicSelection("all");
    opts::TrajectoryWithFrameIndices* tropts = 
        new opts::TrajectoryWithFrameIndices;
    ToolOptions *topts = new ToolOptions;
    opts::AggregateOptions options;
    options.add(bopts).add(sopts).add(tropts).add(topts);
    if (!options.parse(argc, argv))
        exit(-1);

    // Assign tool options to variables
    const double cutoff = topts->cutoff;
    const vector<double> fit_radius = stringToDoubles(topts->fit_radius_str);
    vector<double> fit_area = stringToDoubles(topts->fit_area_str);
    const double esp_radius = topts->esp_radius;
    const double esp_area = topts->esp_area;
    const double weight = topts->weight;
    const bool close_solvent = topts->close_solvent;
    const string esp_file = topts->esp_file;
    const uint seed = topts->seed;
    const uint solvent_size = topts->solvent_size;

    // Print command-line input
    cout << "# " << header << "\n";

    // Do some error checking on tool options

    // Check that some quantities are non-negative
    if (fit_radius[0] <= 0) {
        cout << boost::format("Error: fit_radius (%f) must be greater than "
            "zero.") % fit_radius[0];
        throw(LOOSError());
    }
    for (uint i = 0; i < fit_area.size(); ++i) {
        if (fit_area[i] <= 0) {
            cout << boost::format("Error: fit_area (%f) must be greater than "
                "zero.") % fit_area[i];
            throw(LOOSError());
        }
    }
    if (esp_radius <= 0) {
        cout << boost::format("Error: esp_radius (%f) must be greater than "
            "zero.") % esp_radius;
        throw(LOOSError());
    }
    if (esp_area <= 0) {
        cout << boost::format("Error: esp_area (%f) must be greater than "
            "zero.") % esp_area;
        throw(LOOSError());
    }
    if (weight < 0) {
        cout << boost::format("Error: weight (%f) must be greater than or "
            "equal to zero.") % weight;
        throw(LOOSError());
    }

    // Check that numbers in fit radius are given in descending order
    const uint N_fit_surface = fit_radius.size();
    for (uint i = 1; i < N_fit_surface; ++i) {
        if (fit_radius[i] > fit_radius[i - 1]) {
            cout << boost::format(" Error: fit radius (%s) must be given in "
                "descending order.\n\t%f is smaller than %f") %
                topts->fit_radius_str % fit_radius[i - 1] % fit_radius[i];
            throw(LOOSError());
        }
    }

    // Check that cutoff is greater than the largest number in fit_radius
    if (cutoff <= fit_radius[0]) {
        cout << boost::format("Error: cutoff (%f) must be larger than "
            "fit_radius (%f).") % cutoff % fit_radius[0];
        throw(LOOSError());
    }

    // If the size of fit_area is greater than one, check that it is the same
    // size as fit_radius
    if (fit_area.size() > 1 && fit_area.size() != N_fit_surface) {
        cout << boost::format(" Error: fit radius (%s) and fit area (%s)\n\t"
            "must be the same size or the size of fit_area must be equal to 1.")
            % topts->fit_radius_str % topts->fit_area_str;
        throw(LOOSError());
    }

    // If the size of fit_radius is greater than one but the size of fit_area is
    // one, then copy fit_area so that it is the same size as fit_radius
    if (N_fit_surface > 1 && fit_area.size() == 1)
        for (uint i = 1; i < N_fit_surface; ++i)
            fit_area.push_back(fit_area[0]);

    // Set seed for random number generator
    if (seed == 0)
        randomSeedRNG();
    else
        rng_singleton().seed(seed);

    // Set up random number generator to produce random numbers on (-1, 1)
    base_generator_type& rng_type = rng_singleton();
    boost::uniform_real<> range_pm_one(-1.0, 1.0);
    RNG uniform_pm_one(rng_type, range_pm_one);

    // Build LOOS system and generate atom selections
    AtomicGroup model = tropts->model;
    pTraj traj = tropts->trajectory;
    vector<uint> indices = tropts->frameList();
    AtomicGroup grid_atoms = selectAtoms(model, topts->atom_selection);
    AtomicGroup solvent_atoms = selectAtoms(model, topts->solvent_selection);

    // Number of frames in trajectory
    const uint N_frame = indices.size();

    // Number of atoms in grid_atoms selection
    const uint N_atom = grid_atoms.size();

    // Split solvent_atoms into molecules for cutoff determination
    vector<AtomicGroup> solvent_molecules = solvent_atoms.splitByResidue();

    // Number of molecules in solvent
    const uint N_solvent = solvent_molecules.size();

    // Number of points per spherical shell for the surface on which
    // electrostatic potential (ESP) will be evaluated
    const uint N_shell_esp = (uint) lround(4.0 * M_PI * esp_radius * 
        esp_radius / esp_area);

    // Number of points per spherical shell for the surfaces on which point
    // charges will be fit
    vector<uint> N_shell_fit(N_fit_surface);
    for (uint i = 0; i < N_fit_surface; ++i) 
        N_shell_fit[i] = (uint) lround(4.0 * M_PI * fit_radius[i] *
            fit_radius[i] / fit_area[i]);

    // Get coordinates of grid_atoms in the first frame
    traj->readFrame(indices[0]);
    traj->updateGroupCoords(grid_atoms);
    vector<GCoord> grid_atom_coords;
    grid_atom_coords.reserve(N_atom);
    for (uint i; i < N_atom; ++i)
        grid_atom_coords.push_back(grid_atoms[i]->coords());

    // Generate surface grid for calculating ESP using the coordinates of
    // grid_atoms in the first frame of the trajectory
    vector<GCoord> esp_grid;
    esp_grid.reserve(N_atom * N_shell_esp);
    generateSurfaceGrid(uniform_pm_one, grid_atom_coords, N_atom, esp_radius,
        N_shell_esp, esp_grid);

    // Actual number of grid points in the ESP surface grid
    const uint N_esp_points = esp_grid.size();

    // Generate one or more surface grids for fitting point charges using the
    // coordinates of grid_atoms in the first frame of the trajectory
    vector<vector<GCoord> > fit_grid(N_fit_surface);
    vector<uint> N_fit_points(N_fit_surface);
    for (uint i = 0; i < N_fit_surface; ++i) {
        fit_grid[i].reserve(N_atom * N_shell_fit[i]);
        generateSurfaceGrid(uniform_pm_one, grid_atom_coords, N_atom,
            fit_radius[i], N_shell_fit[i], fit_grid[i]);
        N_fit_points[i] = fit_grid[i].size();

        // Check whether the charge fit later will be underdetermined
        if (weight == 0 && N_esp_points < N_fit_points[i]) {
            cout << boost::format("Error: the charge fit at radius %5.2f will "
                "underdetermined because the number\nof charges to be fit "
                "(%4d) exceeds the number of grid points for which\nsolvent "
                "ESP will be calculated (%4d). Try to increase the area per "
                "point\nfor the fit surface using the '-A' argument or "
                "decrease the area per point\nfor the ESP surface using the "
                "'-B' argument.") % fit_radius[i] % N_fit_points[i] %
                N_esp_points;
            throw(LOOSError());
        }
    }

    // Vector of size N_fit_surface that contains vectors of size N_esp_points +
    // N_atom containing the ESP averaged over the trajectory at grid_atoms and
    // at points in esp_grid due to solvent molecules between two cutoff
    // distances from grid_atoms. The cutoff distances for the first vector are
    // fit_radius[0] and topts->cutoff, and the cutoff distances for subsequent
    // vectors are fit_radius[i] and fit_radius[i - 1]
    rMatrix solvent_esp(N_fit_surface, N_esp_points + N_atom);

    // Vector of boundaries for solvent-containing regions. Stored as squared
    // distances for easier comparison to decide to what region a solvent
    // molecule belongs.
    vector<double> solvent_boundary(N_fit_surface + 2);
    solvent_boundary[0] = cutoff * cutoff;
    for (uint i = 0; i < N_fit_surface; ++i)
        solvent_boundary[i + 1] = fit_radius[i] * fit_radius[i];
    solvent_boundary[N_fit_surface + 1] = esp_radius * esp_radius;

    // For first frame, the frame was already read and the grid_atoms coords
    // were updated. Just update solvent_atoms coords and calculate ESP.
    traj->updateGroupCoords(solvent_atoms);
    calculateFrameESP(solvent_molecules, N_solvent, solvent_size, grid_atoms,
        N_atom, esp_grid, N_esp_points, solvent_boundary, N_fit_surface,
        solvent_esp, close_solvent);

    // If a filename was given, print the ESP grid to an output file
    if (esp_file != "") {
        ofstream esp_out(esp_file.c_str());
        esp_out << "# " << header << "\n";
        esp_out.close();
        writeESP(esp_file, 1, indices[0], solvent_esp);
    }

    // Loop over trajectory and accumulate solvent ESP
    for (uint i = 1; i < N_frame; ++i) {
        traj->readFrame(indices[i]);
        traj->updateGroupCoords(grid_atoms);
        traj->updateGroupCoords(solvent_atoms);
        calculateFrameESP(solvent_molecules, N_solvent, solvent_size,
            grid_atoms, N_atom, esp_grid, N_esp_points, solvent_boundary,
            N_fit_surface, solvent_esp, close_solvent);
        if (esp_file != "")
            writeESP(esp_file, i + 1, indices[i], solvent_esp);
    }

    // Reciprocal distance between grid points in esp_grid and fit_grid
    vector<rMatrix> one_over_R;

    // Loop over fitting surfaces
    for (uint i = 0; i < N_fit_surface; ++i) {

        // Divide solvent_esp by N_frame to get ESP averaged over the trajectory
        for (uint j = 0; j < N_esp_points + N_atom; ++j)
            solvent_esp(i, j) /= N_frame;

        // Compute reciprocal distance between points in esp_grid and fit_grid
        one_over_R.push_back(rMatrix(N_esp_points, N_fit_points[i]));
        for (uint j = 0; j < N_esp_points; ++j)
            for (uint k = 0; k < N_fit_points[i]; ++k)
                one_over_R[i](j, k) = pow(esp_grid[j].distance2(fit_grid[i][k]),
                    -0.5);

    }

    // Fit magnitudes of point charges on farthest fitting surface to reproduce
    // ESP due to solvent farther than that surface on esp_grid
    vector<double> fitting_target;
    fitting_target.reserve(N_esp_points);
    for (uint i = 0; i < N_esp_points; ++i)
        fitting_target.push_back(solvent_esp(0, i));
    vector<vector<double> > fit_charges(N_fit_surface);
    fit_charges[0].reserve(N_fit_points[0]);
    fitChargeGrid(one_over_R[0], fitting_target, N_esp_points, fit_charges[0],
        N_fit_points[0], weight);

    // Calculate ESP at points in esp_grid and at the coordinates of grid_atoms
    // in the first frame due to fit charges on the farthest fitting surface
    rMatrix fit_esp(N_fit_surface, N_esp_points + N_atom);
    for (uint i = 0; i < N_fit_points[0]; ++i) {
        // Loop over grid points in esp_grid
        for (uint j = 0; j < N_esp_points; ++j)
            fit_esp(0, j) += fit_charges[0][i] * one_over_R[0](j, i);
        // Loop over atom coordinates in grid_atom_coords
        for (uint j = 0; j < N_atom; ++j)
            fit_esp(0, j + N_esp_points) += fit_charges[0][i] /
                fit_grid[0][i].distance(grid_atom_coords[j]);
    }

    // Perform charge fits for other fitting surfaces
    for (uint i = 1; i < N_fit_surface; ++i) {

        // Adjust fitting target to include the ESP from solvent in the current
        // region (solvent_esp[i]) and the error in the fit for the previous
        // fitting surface (fitting_target[i - 1] - fit_esp[i - 1])
        for (uint j = 0; j < N_esp_points; ++j)
            fitting_target[j] += solvent_esp(i, j) - fit_esp(i - 1, j);

        // Perform charge fit
        fit_charges[i].reserve(N_fit_points[i]);
        fitChargeGrid(one_over_R[i], fitting_target, N_esp_points,
            fit_charges[i], N_fit_points[i], weight);

        // Calculate ESP at points in esp_grid and at the coordinates of
        // grid_atoms in the first frame due to fit charges on this surface
        for (uint j = 0; j < N_fit_points[i]; ++j) {
            // Loop over grid points in esp_grid
            for (uint k = 0; k < N_esp_points; ++k)
                fit_esp(i, k) += fit_charges[i][j] * one_over_R[i](k, j);
            // Loop over atom coordinates in grid_atom_coords
            for (uint k = 0; k < N_atom; ++k)
                fit_esp(i, k + N_esp_points) += fit_charges[i][j] /
                    fit_grid[i][j].distance(grid_atom_coords[k]);
        }

    }

    // Print ESP at grid_atoms due to solvent averaged over the trajectory
    double solvent_sum;
    cout << "#";
    for (uint i = N_esp_points; i < N_esp_points + N_atom; ++i) {
        solvent_sum = 0;
        for (uint j = 0; j < N_fit_surface; ++j)
            solvent_sum += solvent_esp(j, i);
        cout << boost::format(" %12.8f") % solvent_sum;
    }
    cout << "\n";

    // Print ESP at coordinates of grid_atoms in the first frame due to fit
    // point charges on all fitting surfaces
    double fit_sum;
    cout << "#";
    for (uint i = N_esp_points; i < N_esp_points + N_atom; ++i) {
        fit_sum = 0;
        for (uint j = 0; j < N_fit_surface; ++j)
            fit_sum += fit_esp(j, i);
        cout << boost::format(" %12.8f") % fit_sum;
    }
    cout << "\n";

    // Print coordinates of grid points in esp_grid, the ESP due to solvent
    // averaged over the trajectory, and the ESP due to fit point charges
    for (uint i = 0; i < N_esp_points; ++i) {
        solvent_sum = 0;
        fit_sum = 0;
        for (uint j = 0; j < N_fit_surface; ++j) {
            solvent_sum += solvent_esp(j, i);
            fit_sum += fit_esp(j, i);
        }
        cout << boost::format("# %12.8f %12.8f %12.8f %12.8f %12.8f\n") %
            esp_grid[i][0] % esp_grid[i][1] % esp_grid[i][2] % solvent_sum %
            fit_sum;
    }

    // Print coordinates of point charges in fit_grid and their magnitudes
    for (uint i = 0; i < N_fit_surface; ++i)
        for (uint j = 0; j < N_fit_points[i]; ++j)
            cout << boost::format("%12.8f %12.8f %12.8f %12.8f\n") %
                fit_grid[i][j][0] % fit_grid[i][j][1] % fit_grid[i][j][2] %
                    fit_charges[i][j];

}
