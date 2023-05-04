#include "MagneticSand.h"
#include "Constants.h"
#include "IO.h"

#include <cmath>
#include <format>
#include <numbers>
#include <random>

#include "eigen3/unsupported/Eigen/src/KroneckerProduct/KroneckerTensorProduct.h"
#include <Eigen/src/IterativeLinearSolvers/ConjugateGradient.h>

namespace PhysX {

    template class UniformSource<2>;
    template class UniformSource<3>;
    template class DipoleSource<2>;
    template class DipoleSource<3>;
    template class MagneticSand<2>;
    template class MagneticSand<3>;

    template<int Dim>
    inline UniformSource<Dim>::VectorDd UniformSource<Dim>::field_at_position(const VectorDd & pos) const {
        return _mag * _dir;
    }

    template<int Dim>
    inline DipoleSource<Dim>::VectorDd DipoleSource<Dim>::field_at_position(const VectorDd & pos) const {
        VectorDd r     = pos - _pos;
        VectorDd r_hat = r / r.norm();

        return (3 * _dipole.dot(r_hat) * r_hat - _dipole) / pow(r.norm(), 3);
    }

    template<int Dim>
    inline void MagneticSand<Dim>::writeDescription(YAML::Node & root) const {
        { // Description of particles.
            YAML::Node node;
            node["name"]                    = "particles";
            node["data_mode"]               = "dynamic";
            node["primitive_type"]          = "point_list";
            node["indexed"]                 = false;
            node["color_map"]["enabled"]    = true;
            node["color_map"]["normalized"] = true;
            root["objects"].push_back(node);
        }

        //{ // Description of magnetic dipole.
        //    YAML::Node node;
        //    node["name"]                 = "dipole";
        //    node["data_mode"]            = "dynamic";
        //    node["primitive_type"]       = "line_list";
        //    node["indexed"]              = false;
        //    node["color_map"]["enabled"] = true;
        //    root["objects"].push_back(node);
        //}

        {
            YAML::Node node;
            node["name"]                    = "history";
            node["data_mode"]               = "dynamic";
            node["primitive_type"]          = "line_list";
            node["indexed"]                 = false;
            node["color_map"]["enabled"]    = true;
            node["color_map"]["normalized"] = true;
            root["objects"].push_back(node);
        }
    }

    template<int Dim>
    void MagneticSand<Dim>::writeFrame(const std::string & frameDir, const bool staticDraw) const {
        if constexpr (Dim == 3) { // Write particles.
            std::ofstream fout(frameDir + "/particles.mesh", std::ios::binary);
            IO::writeValue(fout, uint(_spos.size()));
            for (const auto & pos : _spos) {
                // IO::writeValue(fout, Vector2d(_pos[i][0], _pos[i][1]).template cast<float>().eval());
                IO::writeValue(fout, pos.template cast<float>().eval());
                // std::cout << pos << std::endl;
            };
            if constexpr (Dim == 3) {
                for (int i = 0; i < _spos.size(); ++i) {
                    IO::writeValue(fout, VectorDf::Unit(2).eval());
                };
            }
            for (int i = 0; i < _spos.size(); ++i) {
                IO::writeValue(fout, (float) 0.);
                // IO::writeValue(fout, (float) _M[i].norm());
            };
        }

        // if constexpr (Dim == 3) { // Write particles.
        //     std::ofstream fout("C:\\Users\\zhuyu\\Desktop\\ETH\\plotting\\two multi-chain\\particle.txt");
        //     for (int i = 0; i < 100; ++i) {
        //         // IO::writeValue(fout, Vector2d(_pos[i][0], _pos[i][1]).template cast<float>().eval());
        //         fout << "v " << _pos[i][0] << ' ' << _pos[i][1] << ' ' << _pos[i][2] << std::endl;
        //     };
        // }

        //{ // Write the magnetic dipole.
        //    std::ofstream fout(frameDir + "/dipole.mesh", std::ios::binary);
        //    IO::writeValue(fout, (uint(2 * _pos.size())));
        //    for (int i = 0; i < _n; ++i) {
        //        IO::writeValue(fout, _pos[i].template cast<float>().eval());
        //        IO::writeValue(fout, (_pos[i] + 10000 * _M[i]).template cast<float>().eval());
        //    }
        //    for (int i = 0; i < _n; ++i) {
        //        const float v = float(_M[i].norm());
        //        IO::writeValue(fout, v);
        //        IO::writeValue(fout, v);
        //    }
        //}

        if constexpr (Dim == 3) {
            double maximum = -1;
            {
#pragma omp parallel for
                for (int i = 0; i < _sample_pos.size(); ++i) {
                    for (int j = 0; j < _history[i].size() - 1; ++j)
                        if (_history[i][j] == _history[i][j + 1])
                            break;
                        else {
                            maximum = std::max(maximum, cal_field_norm_at_position(_history[i][j]));
                        }
                    maximum = std::max(maximum, cal_field_norm_at_position(_history[i][_history[i].size() - 1]));
                }
            }
            std::cout << maximum << std::endl;

            int cnt = 0;
            for (int i = 0; i < _sample_pos.size(); ++i)
                for (int j = 0; j < _history[i].size() - 1; ++j)
                    if (_history[i][j] == _history[i][j + 1])
                        break;
                    else
                        cnt += 2;

            std::ofstream fout(frameDir + "/history.mesh", std::ios::binary);
            IO::writeValue(fout, uint(cnt));

            for (int i = 0; i < _sample_pos.size(); ++i) {
                for (int j = 0; j < _history[i].size() - 1; ++j)
                    if (_history[i][j] == _history[i][j + 1])
                        break;
                    else {
                        IO::writeValue(fout, _history[i][j].template cast<float>().eval());
                        IO::writeValue(fout, _history[i][j + 1].template cast<float>().eval());
                    }
            }
            if constexpr (Dim == 3) {
                for (int i = 0; i < _sample_pos.size(); ++i) {
                    for (int j = 0; j < _history[i].size() - 1; ++j)
                        if (_history[i][j] == _history[i][j + 1])
                            break;
                        else {
                            IO::writeValue(fout, VectorDf::Unit(2).eval());
                            IO::writeValue(fout, VectorDf::Unit(2).eval());
                        }
                }
            }
            for (int i = 0; i < _sample_pos.size(); ++i) {
                double last_B = cal_field_norm_at_position(_history[i][0]);
                for (int j = 0; j < _history[i].size() - 1; ++j)
                    if (_history[i][j] == _history[i][j + 1])
                        break;
                    else {
                        IO::writeValue(fout, float(last_B / maximum));
                        last_B = cal_field_norm_at_position(_history[i][j + 1]);
                        IO::writeValue(fout, float(last_B / maximum));
                    }
            }
        }
        if constexpr (Dim == 3) {
            {
                static int    frame = -1;
                int           base  = 0;
                std::ofstream fout("C:\\Users\\zhuyu\\Desktop\\bd\\bd.obj");
                std::cout << _sample_pos.size();
                for (int i = 0; i < _sample_pos.size(); ++i) {
                    int cnt = 0;

                    VectorDd lastpos = VectorDd::Zero();

                    for (const auto & pos : _history[i]) {
                        if ((pos - lastpos).norm() < 1e-10)
                            break;
                        lastpos = pos.eval();
                        ++cnt;
                    }

                    if (cnt < 2)
                        continue;

                    cnt = 0;
                    for (const auto & pos : _history[i]) {
                        if ((pos - lastpos).norm() < 1e-10)
                            break;
                        VectorDd mag = cal_field_at_position(pos);
                        fout << "v " << pos[0] << " " << pos[1] << " " << pos[2] << " " << mag[0] << " " << mag[1] << " " << mag[2] << " " << 1 << "\n";
                        lastpos = pos.eval();
                        ++cnt;
                    }

                    for (int j = 0; j < cnt - 1; ++j) {
                        fout << "l " << base + j + 1 << " " << base + j + 2 << "\n";
                    }
                    base += cnt;
                }
                ++frame;
            }
        }
    }

    template<int Dim>
    void MagneticSand<Dim>::saveFrame(const std::string & frameDir) const {
        std::ofstream fout(frameDir + "/sand.sav", std::ios::binary);
        IO::writeValue(fout, _pos);
    }

    template<int Dim>
    void MagneticSand<Dim>::loadFrame(const std::string & frameDir) {
        std::ifstream fin(frameDir + "/sand.sav", std::ios::binary);
        IO::readValue(fin, _pos);
    }

    template<int Dim>
    void MagneticSand<Dim>::initialize() {
        // double degs[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90 };
        ////double degs[] = { 90 };
        // for (auto & deg : degs)
        //     deg = deg / 180 * kPi;

        // static double torqueo[40];
        // static double torquet[40];

        if constexpr (Dim == 3) {
            std::cout << _n << std::endl;
            for (int degs = 0; degs <= 0; degs += 5) {
                double   deg  = double(degs) / 180 * kPi;
                VectorDd axis = VectorDd(std::cos(deg), std::sin(deg), 0);
                MatrixDd rotate;
                rotate << std::cos(deg), -std::sin(deg), 0,
                    std::sin(deg), std::cos(deg), 0,
                    0, 0, 1;

                for (int i = 0; i < _n; ++i)
                    _pos[i] = _spos[i];

                // VectorDd diff1, diff2, difft;
                for (int degt = degs + 90; degt <= degs + 90; degt += 5) {
                    double   degg = double(degt) / 180 * kPi;
                    VectorDd axis = VectorDd(std::cos(degg), std::sin(degg), 0);
                    // VectorDd                            axis   = VectorDd(std::cos(deg + _eps), std::sin(deg + _eps), 0);
                    std::unique_ptr<UniformSource<Dim>> source = std::make_unique<UniformSource<Dim>>(axis, 10);
                    _sources.clear();
                    _sources.push_back(std::move(source));

                    // cal_mag_linearly(_n / 2);

                    // VectorDd torquea = VectorDd::Zero();

                    // for (int i = 0; i < _n / 2; ++i) {
                    //     VectorDd field = VectorDd::Zero();
                    //     for (const auto & source : _sources)
                    //         field += source->field_at_position(_pos[i]);

                    //    torquea += _M[i].cross(field);
                    //}

                    // std::cout << "Torque from the magnetic field of the source: " << torque.norm() << std::endl;

                    // cal_mag_linearly(_n);
                    cal_mag_linearly(_n);

                    // VectorDd       torquea   = VectorDd::Zero();
                    // VectorDd       torqueb   = VectorDd::Zero();
                    // double         force     = 0;
                    // const VectorDd centerhat = _ballCenter[0] / _ballCenter[0].norm();
                    // for (int i = 0; i < _n; ++i) {
                    //     // torqueb += _M[i].cross(cal_field_at_position(_pos[i]));
                    //     VectorDd F = cal_force(i);
                    //     torquea += (_pos[i] - _ballCenter[1]).cross(F);
                    //     torqueb += _ballCenter[1].cross(F);
                    //     force += F.dot(centerhat);
                    // }

                    //// torqueo[degt / 5] = torquea[2];
                    //// torquet[degt / 5] = torqueb[2];

                    // torque1[(degt - degs) / 5][degs / 15] = torquea[2];
                    // torque2[(degt - degs) / 5][degs / 15] = torqueb[2];
                    // forcerr[(degt - degs) / 5][degs / 15] = force;
                    //// torqua1[(degt % 180) / 5][degs / 15]  = torquea[2];
                    //// torqua2[(degt % 180) / 5][degs / 15]  = torqueb[2];
                    //// forcear[(degt % 180) / 5][degs / 15]  = force;
                }

                // std::cout << degs << '\t' << diff1.norm() / (2 * _eps) << '\t' << diff2.norm() / (2 * _eps) << '\t' << difft.norm() / (2 * _eps) << std::endl;
            }
            // printTorque();

            generateSampleParticle();
        }
    }

    template<int Dim>
    void MagneticSand<Dim>::init_M() {
        _M   = std::vector<VectorDd>(_n, VectorDd::Zero());
        _pos = std::vector<VectorDd>(_n, VectorDd::Zero());
    }

    template<int Dim>
    void MagneticSand<Dim>::cal_mag_linearly(int number) {
        std::cout << "Begin calculate\n";
        int row = Dim * number;
        A       = MatrixXd::Identity(row, row);

        // double coef = pow(_r, 3) * 100;
        double coef = (_mu_r - 1) / (_mu_r + 2) * std::pow(_r, 3);

#pragma omp parallel for
        for (int i = 0; i < number; ++i)
            for (int j = 0; j < number; ++j)
                if (i != j) {
                    VectorDd r     = _pos[j] - _pos[i];
                    VectorDd r_hat = r / r.norm();

                    MatrixDd m = 3 * Eigen::kroneckerProduct(r_hat, r_hat.transpose()) - MatrixDd::Identity();

                    m = m * coef / std::pow(r.norm(), 3);

                    A.block(i * Dim, j * Dim, Dim, Dim) = -m;
                }

        b = VectorXd::Zero(row);

        for (int i = 0; i < number; ++i) {
            for (const auto & source : _sources) {
                b.block(i * Dim, 0, Dim, 1) += coef * source->field_at_position(_pos[i]);
            }
        }

        Eigen::ConjugateGradient<MatrixXd, Eigen::Lower | Eigen::Upper> cg;
        cg.compute(A);
        VectorXd x = cg.solve(b);

        for (int i = 0; i < number; ++i)
            _M[i] = x.block(i * Dim, 0, Dim, 1);
        std::cout << "End calculate\n";
    }

    template<int Dim>
    void MagneticSand<Dim>::advance(const double dt) {
        // use RK4 to advect particle's position
        // dx/dt = alpha M / |M|

        static bool   flag  = false;
        static int    frame = 0;
        static double t     = 0;
        t += dt;
        if (flag)
            if constexpr (Dim == 3) {
                double                              deg    = double(frame) / 180 * kPi;
                VectorDd                            axis   = VectorDd(std::cos(deg), std::sin(deg), 0);
                std::unique_ptr<UniformSource<Dim>> source = std::make_unique<UniformSource<Dim>>(axis, 10);

                _sources.clear();
                _sources.push_back(std::move(source));
                cal_mag_linearly(_n);

                generateSampleParticle();
                flag = false;
            }

        double alpha = 20;

        for (int i = 0; i < _sample_pos.size(); ++i, alpha *= -1)
            if (! _fixed[i]) {
                auto & pos = _sample_pos[i];
                if (! collide2particle(pos)) {
                    VectorDd k_1 = alpha * cal_field_dir_at_position(pos);
                    VectorDd k_2 = alpha * cal_field_dir_at_position(pos + 0.5 * dt * k_1);
                    VectorDd k_3 = alpha * cal_field_dir_at_position(pos + 0.5 * dt * k_2);
                    VectorDd k_4 = alpha * cal_field_dir_at_position(pos + dt * k_3);

                    pos += (k_1 + 2 * k_2 + 2 * k_3 + k_4) * dt / 6;
                } else
                    _fixed[i] = true;
                // pos[2] = 0;
            }

        for (int i = 0; i < _sample_pos.size(); ++i) {
            _history[i].push_back(_sample_pos[i]);
        }

        if (t >= 1.0) {
            t = 0;
            frame += 15;
            std::cout << "FRAME: " << frame << std::endl;
            flag = true;
        }
    }

    template<int Dim>
    MagneticSand<Dim>::VectorDd MagneticSand<Dim>::cal_field_dir_at_position(const VectorDd & pos) const {
        VectorDd M = VectorDd::Zero();
        // for (const auto & source : _sources) {
        //     M += source->field_at_position(pos);
        // }

        for (int i = 0; i < _n; ++i) {
            VectorDd r     = pos - _pos[i];
            VectorDd r_hat = r / r.norm();

            M += (3 * _M[i].dot(r_hat) * r_hat - _M[i]) / pow(r.norm(), 3);
        }

        if (M.norm() > _eps)
            return M / M.norm();
        else
            return VectorDd::Zero();
    }

    template<int Dim>
    MagneticSand<Dim>::VectorDd MagneticSand<Dim>::cal_field_at_position(const VectorDd & pos) const {
        VectorDd M = VectorDd::Zero();
        // for (const auto & source : _sources) {
        //     M += source->field_at_position(pos);
        // }

        for (int i = 0; i < _n; ++i) {
            VectorDd r     = pos - _pos[i];
            VectorDd r_hat = r / r.norm();

            M += (3 * _M[i].dot(r_hat) * r_hat - _M[i]) / pow(r.norm(), 3);
        }
        return M;
    }

    template<int Dim>
    double MagneticSand<Dim>::cal_field_norm_at_position(const VectorDd & pos) const {
        VectorDd M = VectorDd::Zero();
        // for (const auto & source : _sources) {
        //     M += source->field_at_position(pos);
        // }
        for (int i = 0; i < _n; ++i) {
            VectorDd r     = pos - _pos[i];
            VectorDd r_hat = r / r.norm();

            M += (3 * _M[i].dot(r_hat) * r_hat - _M[i]) / pow(r.norm(), 3);
        }

        return M.norm();
    }

    template<int Dim>
    MagneticSand<Dim>::VectorDd MagneticSand<Dim>::cal_force(int i) const {
        VectorDd F = VectorDd::Zero();
        for (int j = 0; j < _n; ++j)
            if (i != j) {
                VectorDd r     = _pos[i] - _pos[j];
                VectorDd r_hat = r / r.norm();

                F += 3 * (_M[i].dot(_M[j]) * r_hat + _M[i].dot(r_hat) * _M[j] + _M[j].dot(r_hat) * _M[i] - 5 * _M[i].dot(r_hat) * _M[j].dot(r_hat) * r_hat) / pow(r.norm(), 4);
            }
        return F;
    }

    template<int Dim>
    bool MagneticSand<Dim>::collide2particle(const VectorDd & pos) const {
        // if (std::abs(pos[0]) > 10 || std::abs(pos[1]) > 5 || std::abs(pos[2]) > 5)
        //     return true;
        if (std::abs(pos[0] - 0.5) > 5 || std::abs(pos[1] - 0.5) > 5)
            return true;

        if (0 < pos[0] && pos[0] < 1 && 0 < pos[1] && pos[1] < 1)
            return image[int(pos[0] * 1000)][int(pos[1] * 1000)];
        return false;
        for (const auto & center : _ballCenter)
            if ((pos - center).norm() < 1)
                return true;

        return false;
    }

    template<int Dim>
    void MagneticSand<Dim>::printTorque() {
        // just need one ball
        if (false) {
            for (int i = 0; i < 37; ++i)
                std::cout << torque1[i][0] << std::endl;
        }
        // no rotation for ball
        if (true) {
            for (int i = 0; i < 37; ++i)
                std::cout << torque1[i][0] << std::endl;
            std::cout << std::endl;
            for (int i = 0; i < 37; ++i)
                std::cout << torque2[i][0] << std::endl;
            std::cout << std::endl;
            for (int i = 0; i < 37; ++i)
                std::cout << forcerr[i][0] << std::endl;
        }
        // with rotation, so need some alignment
        if (false) {
            for (int i = 0; i < 37; ++i) {
                std::cout << i * 5;
                for (int j = 0; j < 13; ++j)
                    std::cout << "\t" << torque1[i][j];
                std::cout << std::endl;
            }
            std::cout << std::endl;
            for (int i = 0; i < 37; ++i) {
                std::cout << i * 5;
                for (int j = 0; j < 13; ++j)
                    std::cout << "\t" << torque2[i][j];
                std::cout << std::endl;
            }
            std::cout << std::endl;
            for (int i = 0; i < 37; ++i) {
                std::cout << i * 5;
                for (int j = 0; j < 13; ++j)
                    std::cout << "\t" << forcerr[i][j];
                std::cout << std::endl;
            }
            std::cout << std::endl;
            for (int i = 0; i < 37; ++i) {
                std::cout << i * 5;
                for (int j = 0; j < 13; ++j)
                    std::cout << "\t" << torqua1[i][j];
                std::cout << std::endl;
            }
            std::cout << std::endl;
            for (int i = 0; i < 37; ++i) {
                std::cout << i * 5;
                for (int j = 0; j < 13; ++j)
                    std::cout << "\t" << torqua2[i][j];
                std::cout << std::endl;
            }
            std::cout << std::endl;
            for (int i = 0; i < 37; ++i) {
                std::cout << i * 5;
                for (int j = 0; j < 13; ++j)
                    std::cout << "\t" << forcear[i][j];
                std::cout << std::endl;
            }
        }
    }

    template<int Dim>
    void MagneticSand<Dim>::generateSampleParticle() {
        std::cout << "Generating" << std::endl;

        if constexpr (Dim == 3) {
            static std::mt19937                           rng(std::random_device {}());
            static std::uniform_real_distribution<double> nd(0, 1);

            double r = 1.05;
            // double maxB = 0;
            // int    cnt  = 0;
            //  for (int i = 0; i < 11; ++i) {
            //      double theta = std::acos(double(i) / 20. * 2 - 1);
            //      int    n     = i == 0 ? 1 : 6;
            //      for (int j = 0; j < n; ++j) {
            //          double   phi = j * 2 * kPi / 20.;
            //          VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
            //          if (! collide2particle(pos)) {
            //              double B = cal_field_at_position(pos).norm();
            //              maxB     = std::max(maxB, B);
            //          }
            //      }
            //  }
            //  std::cout << "maxB = " << maxB << std::endl;

            // for (int i = 0; i < 11; ++i) {
            //     double theta = std::acos(double(i) / 20. * 2 - 1);
            //     int    n     = i == 0 ? 1 : 6;
            //     for (int j = 0; j < n; ++j) {
            //         double   phi = j * 2 * kPi / 20.;
            //         VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
            //         if (! collide2particle(pos) && nd(rng) <= cal_field_at_position(pos).norm() / maxB) {
            //             _sample_pos.push_back(pos);
            //             _sample_pos.push_back(pos);

            //            std::vector<VectorDd> temp;
            //            temp.push_back(pos);

            //            _history.push_back(temp);
            //            _history.push_back(temp);
            //            ++cnt;
            //        }
            //    }
            //}

            std::vector<VectorDd> list;

            double maxB = 0;
            int    cnt  = 0;

            for (int i = 0; i < 51; ++i)
                for (int j = 0; j < 51; ++j) {
                    VectorDd pos = VectorDd(1.0 * i / 50.0, 1.0 * j / 50.0, 0);
                    if (! image[int(pos[0] * 1000)][int(pos[1] * 1000)]) {
                        double B = cal_field_at_position(pos).norm();
                        maxB     = std::max(maxB, B);
                        list.push_back(pos);
                    }
                }

            // for (int i = 0; i < 21; ++i) {
            //     double theta = std::acos(double(i) / 40. * 2 - 1);
            //     int    n     = i == 0 ? 1 : 21;
            //     for (int j = 0; j < n; ++j) {
            //         double   phi = j * 2 * kPi / 40.;
            //         VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[0];
            //         if (! collide2particle(pos)) {
            //             double B = cal_field_at_position(pos).norm();
            //             maxB     = std::max(maxB, B);
            //             list.push_back(pos);
            //         }
            //     }
            // }
            //  for (int i = 0; i < 11; ++i) {
            //      double theta = std::acos(double(i) / 20. * 2 - 1);
            //      int    n     = i == 0 ? 1 : 11;
            //      for (int j = 0; j < n; ++j) {
            //          double   phi = j * 2 * kPi / 20.;
            //          VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[0];
            //          if (! collide2particle(pos)) {
            //              double B = cal_field_at_position(pos).norm();
            //              maxB     = std::max(maxB, B);
            //          }
            //      }
            //  }
            //  for (int i = 0; i < 11; ++i) {
            //      double theta = std::acos(double(i) / 20. * 2 - 1);
            //      int    n     = i == 0 ? 1 : 20;
            //      for (int j = 0; j < n; ++j) {
            //          double   phi = j * 2 * kPi / 20.;
            //          VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[1];
            //          if (! collide2particle(pos)) {
            //              double B = cal_field_at_position(pos).norm();
            //              maxB     = std::max(maxB, B);
            //          }
            //      }
            //  }

            std::cout << "maxB = " << maxB << std::endl;

            while (cnt < 500) {
                int      index = list.size() * nd(rng);
                VectorDd pos;
                if (index < list.size())
                    pos = list[index];
                if (nd(rng) <= cal_field_at_position(pos).norm() / maxB) {
                    _sample_pos.push_back(pos);
                    _sample_pos.push_back(pos);

                    std::vector<VectorDd> temp;
                    temp.push_back(pos);

                    _history.push_back(temp);
                    _history.push_back(temp);
                    ++cnt;
                }
            }

            // for (int i = 0; i < 21; ++i) {
            //     double theta = std::acos(double(i) / 40. * 2 - 1);
            //     int    n     = i == 0 ? 1 : 21;
            //     for (int j = 0; j < n; ++j) {
            //         double   phi = j * 2 * kPi / 40.;
            //         VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[0];
            //         if (! collide2particle(pos) && nd(rng) <= cal_field_at_position(pos).norm() / maxB) {
            //             _sample_pos.push_back(pos);
            //             _sample_pos.push_back(pos);

            //            std::vector<VectorDd> temp;
            //            temp.push_back(pos);

            //            _history.push_back(temp);
            //            _history.push_back(temp);
            //            ++cnt;
            //        }
            //    }
            //}

            // for (int i = 0; i < 11; ++i) {
            //     double theta = std::acos(double(i) / 20. * 2 - 1);
            //     int    n     = i == 0 ? 1 : 20;
            //     for (int j = 0; j < n; ++j) {
            //         double   phi = j * 2 * kPi / 20.;
            //         VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[1];
            //         if (! collide2particle(pos) && nd(rng) <= cal_field_at_position(pos).norm() / maxB) {
            //             _sample_pos.push_back(pos);
            //             _sample_pos.push_back(pos);

            //            std::vector<VectorDd> temp;
            //            temp.push_back(pos);

            //            _history.push_back(temp);
            //            _history.push_back(temp);
            //            ++cnt;
            //        }
            //    }
            // }

            // MatrixDd p1, p2, p3, p4;
            // p1 << -1, 0, 0, 0, 1, 0, 0, 0, 1;
            // p2 << 1, 0, 0, 0, -1, 0, 0, 0, 1;
            // p3 << 1, 0, 0, 0, 1, 0, 0, 0, -1;
            // p4 << 0, -1, 0, -1, 0, 0, 0, 0, 1;

            // for (int i = 0; i < cnt; ++i) {
            //     _sample_pos.push_back(p1 * _sample_pos[i * 2]);
            //     _sample_pos.push_back(p1 * _sample_pos[i * 2]);

            //    std::vector<VectorDd> temp;
            //    temp.push_back(p1 * _sample_pos[i * 2]);

            //    _history.push_back(temp);
            //    _history.push_back(temp);
            //}
            // cnt *= 2;

            // for (int i = 0; i < cnt; ++i) {
            //     _sample_pos.push_back(p2 * _sample_pos[i * 2]);
            //     _sample_pos.push_back(p2 * _sample_pos[i * 2]);

            //    std::vector<VectorDd> temp;
            //    temp.push_back(p2 * _sample_pos[i * 2]);

            //    _history.push_back(temp);
            //    _history.push_back(temp);
            //}
            // cnt *= 2;

            // for (int i = 0; i < cnt; ++i) {
            //     _sample_pos.push_back(p3 * _sample_pos[i * 2]);
            //     _sample_pos.push_back(p3 * _sample_pos[i * 2]);

            //    std::vector<VectorDd> temp;
            //    temp.push_back(p3 * _sample_pos[i * 2]);

            //    _history.push_back(temp);
            //    _history.push_back(temp);
            //}
            // cnt *= 2;

            // for (const auto & center : _ballCenter) {
            //     _sample_pos.push_back(VectorDd(1, 0, 0) + center);
            //     _sample_pos.push_back(VectorDd(1, 0, 0) + center);

            //    std::vector<VectorDd> temp;
            //    temp.push_back(VectorDd(1, 0, 0) + center);

            //    _history.push_back(temp);
            //    _history.push_back(temp);
            //    ++cnt;
            //}
            // for (const auto & center : _ballCenter) {
            //    _sample_pos.push_back(VectorDd(-1, 0, 0) + center);
            //    _sample_pos.push_back(VectorDd(-1, 0, 0) + center);

            //    std::vector<VectorDd> temp;
            //    temp.push_back(VectorDd(-1, 0, 0) + center);

            //    _history.push_back(temp);
            //    _history.push_back(temp);
            //    ++cnt;
            //}

            // for (int i = 0; i < cnt; ++i) {
            //     _sample_pos.push_back(p4 * _sample_pos[i * 2]);
            //     _sample_pos.push_back(p4 * _sample_pos[i * 2]);

            //    std::vector<VectorDd> temp;
            //    temp.push_back(p4 * _sample_pos[i * 2]);

            //    _history.push_back(temp);
            //    _history.push_back(temp);
            //}
            // cnt *= 2;

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // double theta = std::acos(0.0);
            // int    N     = 50;
            // int    n     = N;
            // double r     = 1.05;
            // double maxB  = 0.686012;

            // MatrixDd m1;
            // m1 << 1, 0, 0, 0, -1, 0, 0, 0, 1;

            // MatrixDd m2;
            // m2 << -1, 0, 0, 0, 1, 0, 0, 0, 1;

            // for (int j = 0; j < n; ++j) {
            //     double   phi = double(j) * kPi / 2. / n;
            //     VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
            //     if (! collide2particle(pos)) {
            //         double B = cal_field_at_position(pos).norm();
            //         maxB     = std::max(maxB, B);
            //     }
            // }

            // n = 2 * N;
            // for (int j = 0; j < n; ++j) {
            //     double   phi = double(j) * kPi / n - kPi / 2;
            //     VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[0];
            //     if (! collide2particle(pos)) {
            //         double B = cal_field_at_position(pos).norm();
            //         maxB     = std::max(maxB, B);
            //     }
            // }

            // for (int j = 0; j < n; ++j) {
            //     double   phi = double(j) * kPi / n - kPi / 2;
            //     VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[1];
            //     if (! collide2particle(pos)) {
            //         double B = cal_field_at_position(pos).norm();
            //         maxB     = std::max(maxB, B);
            //     }
            // }

            // n = 4 * N;
            // for (int j = 0; j < n; ++j) {
            //     double   phi = double(j) * 2 * kPi / n;
            //     VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[0];
            //     if (! collide2particle(pos)) {
            //         double B = cal_field_at_position(pos).norm();
            //         maxB     = std::max(maxB, B);
            //     }
            // }

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // int cnt = 0;
            //// n       = N;
            //// for (int j = 0; j < n; ++j) {
            ////     double   phi = double(j) * kPi / 2. / n;
            ////     VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
            ////     if (! collide2particle(pos) && nd(rng) <= cal_field_at_position(pos).norm() / maxB) {
            ////         _sample_pos.push_back(pos);
            ////         _sample_pos.push_back(pos);

            ////        std::vector<VectorDd> temp;
            ////        temp.push_back(pos);

            ////        _history.push_back(temp);
            ////        _history.push_back(temp);
            ////        cnt += 1;
            ////    }
            ////}

            // n = 2 * N;
            // for (int j = 0; j < n; ++j) {
            //     double   phi = double(j) * kPi / n - kPi / 2;
            //     VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[0];
            //     if (! collide2particle(pos) && nd(rng) <= cal_field_at_position(pos).norm() / maxB) {
            //         _sample_pos.push_back(pos);
            //         _sample_pos.push_back(pos);

            //        std::vector<VectorDd> temp;
            //        temp.push_back(pos);

            //        _history.push_back(temp);
            //        _history.push_back(temp);
            //        cnt += 1;
            //    }
            //}

            //// for (int j = 0; j < n; ++j) {
            ////     double   phi = double(j) * kPi / n;
            ////     VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[1];
            ////     if (! collide2particle(pos) && nd(rng) <= cal_field_at_position(pos).norm() / maxB) {
            ////         _sample_pos.push_back(pos);
            ////         _sample_pos.push_back(pos);

            ////        std::vector<VectorDd> temp;
            ////        temp.push_back(pos);

            ////        _history.push_back(temp);
            ////        _history.push_back(temp);
            ////        cnt += 1;
            ////    }
            ////}

            //// n = 4 * N;
            ////     for (int j = 0; j < n; ++j) {
            ////         double   phi = double(j) * 2 * kPi / n;
            ////         VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[2];
            ////         if (! collide2particle(pos) && nd(rng) <= cal_field_at_position(pos).norm() / maxB) {
            ////             {
            ////                 _sample_pos.push_back(pos);
            ////                 _sample_pos.push_back(pos);

            ////            std::vector<VectorDd> temp;
            ////            temp.push_back(pos);

            ////            _history.push_back(temp);
            ////            _history.push_back(temp);
            ////            cnt += 1;
            ////        }
            ////        {
            ////            _sample_pos.push_back(m2 * pos - _ballCenter[2] + _ballCenter[3]);
            ////            _sample_pos.push_back(m2 * pos - _ballCenter[2] + _ballCenter[3]);

            ////            std::vector<VectorDd> temp;
            ////            temp.push_back(m2 * pos - _ballCenter[2] + _ballCenter[3]);

            ////            _history.push_back(temp);
            ////            _history.push_back(temp);
            ////            cnt += 1;
            ////        }
            ////    }
            ////}
            //// for (int j = 0; j < n; ++j) {
            ////    double   phi = double(j) * 2 * kPi / n;
            ////    VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[4];
            ////    if (! collide2particle(pos) && nd(rng) <= cal_field_at_position(pos).norm() / maxB) {
            ////        _sample_pos.push_back(pos);
            ////        _sample_pos.push_back(pos);

            ////        std::vector<VectorDd> temp;
            ////        temp.push_back(pos);

            ////        _history.push_back(temp);
            ////        _history.push_back(temp);
            ////        cnt += 1;
            ////    }
            ////}
            //// for (int j = 0; j < n; ++j) {
            ////    double   phi = double(j) * 2 * kPi / n;
            ////    VectorDd pos = r * VectorDd(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)) + _ballCenter[0];
            ////    if (! collide2particle(pos) && nd(rng) <= cal_field_at_position(pos).norm() / maxB) {
            ////        _sample_pos.push_back(pos);
            ////        _sample_pos.push_back(pos);

            ////        std::vector<VectorDd> temp;
            ////        temp.push_back(pos);

            ////        _history.push_back(temp);
            ////        _history.push_back(temp);
            ////        cnt += 1;
            ////    }
            ////}

            // for (int i = 0; i < cnt; ++i) {
            //     _sample_pos.push_back(m1 * _sample_pos[i * 2]);
            //     _sample_pos.push_back(m1 * _sample_pos[i * 2]);

            //    std::vector<VectorDd> temp;
            //    temp.push_back(m1 * _sample_pos[i * 2]);

            //    _history.push_back(temp);
            //    _history.push_back(temp);
            //}
            // cnt *= 2;

            // for (int i = 0; i < cnt; ++i) {
            //     _sample_pos.push_back(m2 * _sample_pos[i * 2]);
            //     _sample_pos.push_back(m2 * _sample_pos[i * 2]);

            //    std::vector<VectorDd> temp;
            //    temp.push_back(m2 * _sample_pos[i * 2]);

            //    _history.push_back(temp);
            //    _history.push_back(temp);
            //}
            // cnt *= 2;
        }

        _fixed = std::vector<bool>(_sample_pos.size(), false);
    }

} // namespace PhysX
