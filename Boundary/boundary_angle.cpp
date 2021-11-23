#include <iostream>
#include <tuple>
#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::tuple;

constexpr double R2D = 180.0 / M_PI; 

struct mounting_angle {
    double mounting_roll;
    double mounting_pitch;
    double mounting_heading;

    mounting_angle(double _roll, double _pitch, double _heading) : mounting_roll(_roll), mounting_pitch(_pitch), mounting_heading(_heading){}
};

struct attitude_angle {
    double attitude_roll;
    double attitude_pitch;
    double attitude_heading;

    attitude_angle(double _roll, double _pitch, double _heading) : attitude_roll(_roll), attitude_pitch(_pitch), attitude_heading(_heading) {}
};

Eigen::Matrix3d rollMatrix(const double rollValue) {
    Eigen::Matrix3d result {};
    result << 1.0, 0.0, 0.0, 0.0, cos(rollValue), -sin(rollValue), 0.0, sin(rollValue), cos(rollValue);
    return result;
}

Eigen::Matrix3d pitchMatrix(const double pitchValue) {
    Eigen::Matrix3d result {};
    result << cos(pitchValue), 0.0, sin(pitchValue), 0.0, 1.0, 0.0, -sin(pitchValue), 0.0, cos(pitchValue);
    return result;
}

Eigen::Matrix3d headingMatrix(const double headingValue)
{
    Eigen::Matrix3d result {};
    result << cos(headingValue), -sin(headingValue), 0.0, sin(headingValue), cos(headingValue), 0.0, 0.0, 0.0, 1.0;
    return result;
}

Eigen::Matrix3d attitude2Matrix(const float roll, const float pitch = 0.0, const float heading = 0.0)
{
    return headingMatrix(heading) * pitchMatrix(pitch) * rollMatrix(roll);
}

Eigen::Matrix3d attitude2Matrix(const mounting_angle m) {
    return attitude2Matrix(m.mounting_roll, m.mounting_pitch, m.mounting_heading);
}

Eigen::Matrix3d attitude2Matrix(const attitude_angle a) {
    return attitude2Matrix(a.attitude_roll, a.attitude_pitch, a.attitude_heading);
}


tuple<double, double> boundary_angle(mounting_angle mounting_tx, attitude_angle attitude_tx,
                                     mounting_angle mounting_rx, attitude_angle attitude_rx,
                                     double txTilt) {
    Eigen::Vector3d txIdeal {1.0, 0.0, 0.0};
        Eigen::Vector3d txMount = attitude2Matrix(mounting_tx) * txIdeal;
        Eigen::Vector3d txGeo = attitude2Matrix(attitude_tx) * txMount;

        Eigen::Vector3d rxIdeal {0.0, 1.0, 0.0};
        Eigen::Vector3d rxMount = attitude2Matrix(mounting_rx) * rxIdeal;
        Eigen::Vector3d rxGeo = attitude2Matrix(attitude_rx) * rxMount;

        const double m1 = txGeo(0), n1 = txGeo(1), p1 = txGeo(2);
        const double m2 = rxGeo(0), n2 = rxGeo(1), p2 = rxGeo(2);
        double denominator = m1 * n2 - m2 * n1;

        Eigen::Vector3d directionVect = txGeo.cross(rxGeo);

        const double M1 = directionVect(0);
        const double M2 = directionVect(1);
        const double A = n2 * sin(txTilt) / denominator;
        const double B = m2 * sin(txTilt) / denominator;
        const double k = n1 / denominator;
        const double t = m1 / denominator;
        const double c = M1 * k - M2 * t;
        const double T = M1 * A - M2 * B;

        const double coeff1 = c * c - k * k - t * t;
        const double coeff2 = 2.0 * (T * c - A * k - B * t);
        const double coeff3 = T * T + 1.0 - A * A - B * B;

        const double delta = coeff2 * coeff2 - 4.0 * coeff3 * coeff1;
        const double l1 = (-coeff2 + sqrt(delta)) / (2.0 * coeff1);
        const double l2 = (-coeff2 - sqrt(delta)) / (2.0 * coeff1);


        return std::make_tuple(asin(l2) * R2D, asin(l1) * R2D);
}


int main(int argv, char **argc) {
    const double txRoll = 1.48 * M_PI / 180;
    const double txPitch = -0.32 * M_PI / 180;
    const double txHeading = 0.81 * M_PI / 180;
    const double txRollBias = -0.32 * M_PI / 180;
    const double txPitchBias = 0.14 * M_PI / 180;
    const double txHeadingBias = -0.05 * M_PI / 180;
    const double txTilt = 0.5 * M_PI / 180;
    const double rxRoll = 1.51 * M_PI / 180;
    const double rxPitch = -0.31 * M_PI / 180;
    const double rxHeading = 0.76 * M_PI / 180;
    const double rxRollBias = 0.21 * M_PI / 180;
    const double rxPitchBias = 0 * M_PI / 180;
    const double rxHeadingBias = -0.19 * M_PI / 180;

    const mounting_angle mouting_tx(txRollBias, txPitchBias, txHeadingBias);
    const mounting_angle mouting_rx(rxRollBias, rxPitchBias, rxHeadingBias);
    const attitude_angle attitude_tx(txRoll, txPitch, txHeading);
    const attitude_angle attitude_rx(rxRoll, rxPitch, rxHeading);

    tuple<double, double> angle = boundary_angle(mouting_tx, attitude_tx, mouting_rx, attitude_rx, txTilt);
    cout << std::get<0>(angle) << endl;
    cout << std::get<1>(angle) << endl;
}
