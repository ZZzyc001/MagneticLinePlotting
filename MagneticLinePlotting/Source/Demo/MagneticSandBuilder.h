#pragma once

#include "Constants.h"
#include "MagneticSand.h"
#include <iostream>
#include <memory>
#include <random>

namespace PhysX {

    class MagneticSandBuilder {
    public:
        template<int Dim>
        static std::unique_ptr<MagneticSand<Dim>> build(const int option, double ba) {
            auto magneticSand = std::make_unique<MagneticSand<Dim>>();

            switch (option) {
            case 0:
                buildCase0(magneticSand.get());
                break;
            case 1:
                buildCase1(magneticSand.get());
                break;
            case 2:
                buildCase2(magneticSand.get());
                break;
            case 3:
                buildCase3(magneticSand.get());
                break;
            case 4:
                buildCase4(magneticSand.get());
                break;
            case 5:
                buildCase5(magneticSand.get());
                break;
            case 6:
                buildCase6(magneticSand.get());
                break;
            case 7:
                buildCase7(magneticSand.get());
                break;
            case 8:
                buildCase8(magneticSand.get());
                break;
            case 9:
                buildCase9(magneticSand.get(), ba);
                break;
            case 10:
                buildCase10(magneticSand.get());
                break;
            case 11:
                buildCase11(magneticSand.get());
                break;
            case 12:
                buildCase12(magneticSand.get());
                break;
            case 13:
                buildCase13(magneticSand.get());
                break;
            case 14:
                buildCase14(magneticSand.get());
                break;
            case 15:
                buildCase15(magneticSand.get());
                break;
            case 16:
                buildCase16(magneticSand.get());
                break;
            case 17:
                buildCase17(magneticSand.get());
                break;
            case 18:
                buildCase18(magneticSand.get());
                break;
            case 19:
                buildCase19(magneticSand.get());
                break;
            case 20:
                buildCase20(magneticSand.get());
                break;
            case 21:
                buildCase21(magneticSand.get());
                break;
            case 22:
                buildCase22(magneticSand.get());
                break;
            case 23:
                buildCase23(magneticSand.get());
                break;
            case 24:
                buildCase24(magneticSand.get());
                break;
            case 25:
                buildCase25(magneticSand.get());
                break;
            case 26:
                buildCase26(magneticSand.get());
                break;
            case 27:
                buildCase27(magneticSand.get());
                break;
            case 28:
                buildCase28(magneticSand.get());
                break;
            case 29:
                buildCase29(magneticSand.get());
                break;
            case 30:
                buildCase30(magneticSand.get());
                break;
            case 31:
                buildCase31(magneticSand.get());
                break;
            case 32:
                buildCasebd(magneticSand.get());
                break;
            default:
                return nullptr;
            }

            return magneticSand;
        }

    protected:
        template<int Dim>
        static void buildCase0(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            int xScale = 20;
            int yScale = 10;
            int zScale = 10;

            static std::mt19937                           rng(42);
            static std::uniform_real_distribution<double> nd(0, 2 * kPi);

            if constexpr (Dim == 2) {
                for (int i = 0; i < xScale; ++i)
                    for (int j = 3; j < 3 + yScale; ++j) {
                        double   theta = nd(rng);
                        VectorDd pos(i * 4. + .1 * cos(theta), j * 4. + .1 * sin(theta));
                        magneticSand->_pos.push_back(pos);
                        ++magneticSand->_n;
                    }

                // clang-format off
                std::unique_ptr<DipoleSource<Dim>> source = std::make_unique<DipoleSource<Dim>>(
						    VectorDd(2 * xScale, -0.5 * yScale) , 
						    VectorDd(0, 1) * 1e3);
                magneticSand->_sources.push_back(std::move(source));
                // clang-format on

                magneticSand->init_M();
            } else if constexpr (Dim == 3) {
                for (int i = 0; i < xScale; ++i)
                    for (int j = 0; j < yScale; ++j)
                        for (int k = 0; k < zScale; ++k) {
                            double theta = nd(rng);
                            double phi   = nd(rng);
                            // VectorDd pos(i * 4., j * 4., k * 4.);
                            VectorDd pos(i * 4. + .1 * sin(theta) * cos(phi), j * 4. + .1 * cos(theta), k * 4. + .1 * sin(theta) * sin(phi));
                            magneticSand->_pos.push_back(pos);
                            ++magneticSand->_n;
                        }

                // clang-format off
                std::unique_ptr<DipoleSource<Dim>> source = std::make_unique<DipoleSource<Dim>>(
				    VectorDd(2 * xScale, -1.0 * yScale, 2 * zScale) , 
				    VectorDd(0, 1, 0) * 1e2);
                magneticSand->_sources.push_back(std::move(source));
                // clang-format on

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase1(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                static std::mt19937                           rng(std::random_device {}());
                static std::uniform_real_distribution<double> nd_r_1_1(-0.5, 0.5);
                static std::uniform_real_distribution<double> nd_theta_pi(0, kPi);
                static std::uniform_real_distribution<double> nd_theta_2pi(0, 2 * kPi / 3);

                const double biases[] = { 0, 6, -6, -14, 14 };
                const double varies[] = { 0, 2 * kPi / 3, 4 * kPi / 3 };

                double coef = 2 / std::sqrt(3);

                for (const auto & bias : biases) {
                    double   theta  = nd_theta_pi(rng);
                    double   r      = nd_r_1_1(rng);
                    VectorDd origin = VectorDd(0, r * std::cos(theta) + bias, r * std::sin(theta));
                    for (int i = 0; (origin + i * 2 * VectorDd::Unit(0)).norm() + 1 < 20; ++i) {
                        double theta = nd_theta_2pi(rng);
                        for (const auto & vary : varies) {
                            VectorDd pos = VectorDd(2 * i, coef * std::cos(theta + vary), coef * std::sin(theta + vary)) + origin;

                            if (pos.norm() + 1 > 20)
                                continue;

                            // bool flag = true;
                            // for (const auto & rpos : magneticSand->_pos)
                            //     if ((pos - rpos).norm() < 2 + magneticSand->_eps) {
                            //         flag = false;
                            //         break;
                            //     }
                            // if (flag) {
                            magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                            if (i)
                                magneticSand->_spos.push_back(VectorDd(-2 * i, coef * std::cos(theta + vary), coef * std::sin(theta + vary)) + origin), ++magneticSand->_n;
                            //}
                        }
                    }
                }

                magneticSand->_ballCenter.push_back(VectorDd(0, 0, 0));

                const double   rotate_angle = kPi / 6;
                const VectorDd axis         = VectorDd(std::cos(rotate_angle), std::sin(rotate_angle), 0);

                // MatrixDd rotate;
                // rotate << std::cos(rotate_angle), -std::sin(rotate_angle), 0,
                //     std::sin(rotate_angle), std::cos(rotate_angle), 0,
                //     0, 0, 1;

                // for (auto & pos : magneticSand->_pos)
                //     pos = rotate * pos;

                std::unique_ptr<UniformSource<Dim>> source = std::make_unique<UniformSource<Dim>>(axis, 10);
                magneticSand->_sources.push_back(std::move(source));

                magneticSand->init_M();

                // for (int i = 0; i < 200; ++i) {
                //     double theta = i * 2 * kPi / 200;

                //    VectorDd pos = VectorDd(20 * std::sin(theta), 20 * std::cos(theta), 0);
                //    magneticSand->_sample_pos.push_back(pos);
                //    magneticSand->_sample_pos.push_back(pos);

                //    std::vector<VectorDd> temp;

                //    temp.push_back(pos);
                //    magneticSand->_history.push_back(temp);
                //    magneticSand->_history.push_back(temp);
                //}
            }
        }

        template<int Dim>
        static void buildCase2(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                static std::mt19937                           rng(std::random_device {}());
                static std::uniform_real_distribution<double> nd_r_1_1(-0.5, 0.5);
                static std::uniform_real_distribution<double> nd_theta_pi(0, kPi);
                static std::uniform_real_distribution<double> nd_theta_2pi(0, 2 * kPi / 3);

                magneticSand->_r = 1.0 / 30.0;

                const double biases[] = { 0, 14, -14, 25, -25 };
                const double varies[] = { 0, kPi / 3, 2 * kPi / 3, kPi, 4 * kPi / 3, 5 * kPi / 3 };
                const double starts[] = { 0, 0, kPi / 6, 0, kPi / 6 - std::atan(1.0 / (3.0 * std::sqrt(3))), kPi / 6 + std::atan(1.0 / (3.0 * std::sqrt(3))) };
                const double radius[] = { 2,
                                          4,
                                          2 * std::sqrt(3),
                                          6,
                                          2 * std::sqrt(7),
                                          2 * std::sqrt(7) };
                const int    counts[] = { 6, 3, 3, 1, 1 };

                for (int b = 0; b < 5; ++b) {
                    // double   theta  = nd_theta_pi(rng);
                    // double   r      = nd_r_1_1(rng);
                    VectorDd origin = VectorDd(0, biases[b], 0);
                    // VectorDd origin = VectorDd(0, r * std::cos(theta) + bias, r * std::sin(theta));
                    for (int i = -29; i <= 29; i += 2) {
                        VectorDd pos = VectorDd(i, 0, 0) + origin;
                        if (pos.norm() + 1 <= 30)
                            magneticSand->_spos.push_back(pos * magneticSand->_r), ++magneticSand->_n;
                        double theta = nd_theta_2pi(rng);
                        for (int j = 0; j < counts[b]; ++j)
                            for (const auto & vary : varies) {
                                VectorDd pos = VectorDd(i, radius[j] * std::cos(theta + vary + starts[j]), radius[j] * std::sin(theta + vary + starts[j])) + origin;

                                if (pos.norm() + 1 <= 30)
                                    magneticSand->_spos.push_back(pos * magneticSand->_r), ++magneticSand->_n;

                                // bool flag = true;
                                // for (const auto & rpos : magneticSand->_pos)
                                //     if ((pos - rpos).norm() < 2 + magneticSand->_eps) {
                                //         flag = false;
                                //         break;
                                //     }
                                // if (flag) {
                                // magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                                //}
                            }
                    }
                }

                // const double   rotate_angle = 0 * kPi / 6;
                // const VectorDd axis = VectorDd(std::cos(rotate_angle), std::sin(rotate_angle), 0);

                // MatrixDd rotate;
                // rotate << std::cos(rotate_angle), -std::sin(rotate_angle), 0,
                //     std::sin(rotate_angle), std::cos(rotate_angle), 0,
                //     0, 0, 1;

                // for (auto & pos : magneticSand->_pos)
                //     pos = rotate * pos;

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + 1.1 * VectorDd::Unit(0));
                    magneticSand->_spos[i] -= 1.1 * VectorDd::Unit(0);
                }
                magneticSand->_n *= 2;

                magneticSand->_ballCenter.push_back(VectorDd(-1.1, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(1.1, 0, 0));

                // std::unique_ptr<UniformSource<Dim>> source = std::make_unique<UniformSource<Dim>>(axis, 10);
                // magneticSand->_sources.push_back(std::move(source));

                magneticSand->init_M();

                // for (int i = 0; i < 21; ++i) {
                //     double theta = std::acos(double(i) / 20 * 2 - 1);
                //     int    n     = (i == 0 || i == 20) ? 1 : 20;
                //     for (int j = 0; j < n; ++j) {
                //         double   phi  = j * 2 * kPi / 20;
                //         VectorDd pos  = VectorDd(20 * std::sin(theta) * std::cos(phi) - 21, 20 * std::sin(theta) * std::sin(phi), 20 * std::cos(theta));
                //         VectorDd pos1 = VectorDd(20 * std::sin(theta) * std::cos(phi) + 21, 20 * std::sin(theta) * std::sin(phi), 20 * std::cos(theta));
                //         magneticSand->_sample_pos.push_back(pos);
                //         magneticSand->_sample_pos.push_back(pos);
                //         magneticSand->_sample_pos.push_back(pos1);
                //         magneticSand->_sample_pos.push_back(pos1);

                //        std::vector<VectorDd> temp;
                //        temp.push_back(pos);

                //        std::vector<VectorDd> temp1;
                //        temp1.push_back(pos1);

                //        magneticSand->_history.push_back(temp);
                //        magneticSand->_history.push_back(temp);
                //        magneticSand->_history.push_back(temp1);
                //        magneticSand->_history.push_back(temp1);
                //    }
                //}
            }
        }

        template<int Dim>
        static void buildCase3(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                const double biases[] = { 0, 10, -10 };

                for (const auto & bias : biases)
                    for (int i = -20; i <= 20; i += 2)
                        for (int j = -20; j <= 20; j += 2) {
                            VectorDd pos = VectorDd(i, bias, j);

                            if (pos.norm() + 1 > 20)
                                continue;

                            magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                        }

                magneticSand->_ballCenter.push_back(VectorDd(0, 0, 0));

                const double   rotate_angle = 1 * kPi / 6;
                const VectorDd axis         = VectorDd(std::cos(rotate_angle), std::sin(rotate_angle), 0);

                // MatrixDd rotate;
                // rotate << std::cos(rotate_angle), -std::sin(rotate_angle), 0,
                //     std::sin(rotate_angle), std::cos(rotate_angle), 0,
                //     0, 0, 1;

                // for (auto & pos : magneticSand->_pos)
                //     pos = rotate * pos;

                std::unique_ptr<UniformSource<Dim>> source = std::make_unique<UniformSource<Dim>>(axis, 10);
                magneticSand->_sources.push_back(std::move(source));

                magneticSand->init_M();

                // for (int i = 0; i < 200; ++i) {
                //     double theta = i * 2 * kPi / 200;

                //    VectorDd pos = VectorDd(20 * std::sin(theta), 20 * std::cos(theta), 0);
                //    magneticSand->_sample_pos.push_back(pos);
                //    magneticSand->_sample_pos.push_back(pos);

                //    std::vector<VectorDd> temp;

                //    temp.push_back(pos);
                //    magneticSand->_history.push_back(temp);
                //    magneticSand->_history.push_back(temp);
                //}
            }
        }

        template<int Dim>
        static void buildCase4(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r      = 1.0 / 30.0;
                const double biases[] = { 0, 20, -20 };
                const int    counts[] = { 2, 1, 1 };
                const double e        = 0.99;

                for (int b = 0; b < 3; ++b)
                    for (int layer = -counts[b]; layer <= counts[b]; ++layer)
                        for (int i = -29 + (layer & 1); i <= 29 - (layer & 1); i += 2)
                            for (int j = -29 + (layer & 1); j <= 29 - (layer & 1); j += 2) {
                                VectorDd pos = VectorDd(i, j, biases[b] + layer * std::sqrt(2));

                                if (std::pow(pos[0], 2) / (29 * 29 - std::pow(pos[2], 2)) + std::pow(pos[1], 2) / ((1 - e * e) * (29 * 29 - std::pow(pos[2], 2))) <= 1)
                                    magneticSand->_spos.push_back(pos * magneticSand->_r), ++magneticSand->_n;
                            }

                // const double   rotate_angle = 0 * kPi / 6;
                // const VectorDd axis         = VectorDd(std::cos(rotate_angle), std::sin(rotate_angle), 0);

                // MatrixDd rotate;
                // rotate << std::cos(rotate_angle), -std::sin(rotate_angle), 0,
                //     std::sin(rotate_angle), std::cos(rotate_angle), 0,
                //     0, 0, 1;

                // for (auto & pos : magneticSand->_pos)
                //     pos = rotate * pos;

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + 1.1 * VectorDd::Unit(0));
                    magneticSand->_spos[i] -= 1.1 * VectorDd::Unit(0);
                }
                magneticSand->_n *= 2;

                magneticSand->_ballCenter.push_back(VectorDd(-1.1, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(1.1, 0, 0));

                // std::unique_ptr<UniformSource<Dim>> source = std::make_unique<UniformSource<Dim>>(axis, 10);
                // magneticSand->_sources.push_back(std::move(source));

                magneticSand->init_M();

                // for (int i = 0; i < 21; ++i) {
                //     double theta = std::acos(double(i) / 20 * 2 - 1);
                //     int    n     = (i == 0 || i == 20) ? 1 : 20;
                //     for (int j = 0; j < n; ++j) {
                //         double   phi  = j * 2 * kPi / 20;
                //         VectorDd pos  = VectorDd(20 * std::sin(theta) * std::cos(phi) - 21, 20 * std::sin(theta) * std::sin(phi), 20 * std::cos(theta));
                //         VectorDd pos1 = VectorDd(20 * std::sin(theta) * std::cos(phi) + 21, 20 * std::sin(theta) * std::sin(phi), 20 * std::cos(theta));
                //         magneticSand->_sample_pos.push_back(pos);
                //         magneticSand->_sample_pos.push_back(pos);
                //         magneticSand->_sample_pos.push_back(pos1);
                //         magneticSand->_sample_pos.push_back(pos1);

                //        std::vector<VectorDd> temp;
                //        temp.push_back(pos);

                //        std::vector<VectorDd> temp1;
                //        temp1.push_back(pos1);

                //        magneticSand->_history.push_back(temp);
                //        magneticSand->_history.push_back(temp);
                //        magneticSand->_history.push_back(temp1);
                //        magneticSand->_history.push_back(temp1);
                //    }
                //}
            }
        }

        template<int Dim>
        static void buildCase5(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                for (int k = -20; k <= 20; k += 2)
                    for (int i = -20; i <= 20; i += 2)
                        for (int j = -20; j <= 20; j += 2) {
                            VectorDd pos = VectorDd(i, k, j);

                            if (pos.norm() + 1 > 20)
                                continue;

                            magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                        }

                magneticSand->_ballCenter.push_back(VectorDd(0, 0, 0));

                const double   rotate_angle = 1 * kPi / 6;
                const VectorDd axis         = VectorDd(std::cos(rotate_angle), std::sin(rotate_angle), 0);

                // for (auto & pos : magneticSand->_pos)
                //     pos = pos;

                std::unique_ptr<UniformSource<Dim>> source = std::make_unique<UniformSource<Dim>>(axis, 10);
                magneticSand->_sources.push_back(std::move(source));

                magneticSand->init_M();

                // for (int i = 0; i < 200; ++i) {
                //     double theta = i * 2 * kPi / 200;

                //    VectorDd pos = VectorDd(20 * std::sin(theta), 20 * std::cos(theta), 0);
                //    magneticSand->_sample_pos.push_back(pos);
                //    magneticSand->_sample_pos.push_back(pos);

                //    std::vector<VectorDd> temp;

                //    temp.push_back(pos);
                //    magneticSand->_history.push_back(temp);
                //    magneticSand->_history.push_back(temp);
                //}
            }
        }

        template<int Dim>
        static void buildCase6(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1.0 / 30.0;

                for (int i = -290; i <= 290; i += 29)
                    for (int j = -290; j <= 290; j += 29)
                        for (int k = -290; k <= 290; k += 29) {
                            VectorDd pos = VectorDd(i / 10., j / 10., k / 10.);

                            if (pos.norm() + 1 <= 30)
                                magneticSand->_spos.push_back(pos * magneticSand->_r), ++magneticSand->_n;
                        }

                // const double   rotate_angle = 0 * kPi / 6;
                // const VectorDd axis         = VectorDd(std::cos(rotate_angle), std::sin(rotate_angle), 0);

                // for (auto & pos : magneticSand->_pos)
                //     pos = pos;

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + 1.1 * VectorDd::Unit(0));
                    magneticSand->_spos[i] -= 1.1 * VectorDd::Unit(0);
                }
                magneticSand->_n *= 2;

                magneticSand->_ballCenter.push_back(VectorDd(-1.1, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(1.1, 0, 0));

                // std::unique_ptr<UniformSource<Dim>>
                //     source = std::make_unique<UniformSource<Dim>>(axis, 10);
                // magneticSand->_sources.push_back(std::move(source));

                magneticSand->init_M();

                // for (int i = 0; i < 21; ++i) {
                //     double theta = std::acos(double(i) / 20 * 2 - 1);
                //     int    n     = (i == 0 || i == 20) ? 1 : 20;
                //     for (int j = 0; j < n; ++j) {
                //         double   phi  = j * 2 * kPi / 20;
                //         VectorDd pos  = VectorDd(20 * std::sin(theta) * std::cos(phi) - 21, 20 * std::sin(theta) * std::sin(phi), 20 * std::cos(theta));
                //         VectorDd pos1 = VectorDd(20 * std::sin(theta) * std::cos(phi) + 21, 20 * std::sin(theta) * std::sin(phi), 20 * std::cos(theta));
                //         magneticSand->_sample_pos.push_back(pos);
                //         magneticSand->_sample_pos.push_back(pos);
                //         magneticSand->_sample_pos.push_back(pos1);
                //         magneticSand->_sample_pos.push_back(pos1);

                //        std::vector<VectorDd> temp;
                //        temp.push_back(pos);

                //        std::vector<VectorDd> temp1;
                //        temp1.push_back(pos1);

                //        magneticSand->_history.push_back(temp);
                //        magneticSand->_history.push_back(temp);
                //        magneticSand->_history.push_back(temp1);
                //        magneticSand->_history.push_back(temp1);
                //    }
                //}
            }
        }

        template<int Dim>
        static void buildCase7(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                static std::mt19937                           rng(std::random_device {}());
                static std::uniform_real_distribution<double> nd_r_1_1(-0.5, 0.5);
                static std::uniform_real_distribution<double> nd_theta_pi(0, kPi);
                static std::uniform_real_distribution<double> nd_theta_2pi(0, 2 * kPi / 3);

                magneticSand->_r = 1.0 / 30.0;

                const double biases[] = { 0, 9, -9, 14, -14 };
                const double varies[] = { 0, kPi / 3, 2 * kPi / 3, kPi, 4 * kPi / 3, 5 * kPi / 3 };
                const double starts[] = { 0, 0, kPi / 6 };
                const double radius[] = { 2, 4, 2 * std::sqrt(3) };
                const int    counts[] = { 3, 1, 1, 0, 0 };

                for (int b = 0; b < 5; ++b) {
                    VectorDd origin = VectorDd(0, biases[b], 0);
                    for (int i = -15; i <= 15; i += 2) {
                        VectorDd pos = VectorDd(i, 0, 0) + origin;
                        if (pos.norm() + 1 <= 16)
                            magneticSand->_spos.push_back((pos - VectorDd(14, 0, 0)) * magneticSand->_r), ++magneticSand->_n;
                        double theta = nd_theta_2pi(rng);
                        for (int j = 0; j < counts[b]; ++j)
                            for (const auto & vary : varies) {
                                VectorDd pos = VectorDd(i, radius[j] * std::cos(theta + vary + starts[j]), radius[j] * std::sin(theta + vary + starts[j])) + origin;

                                if (pos.norm() + 1 <= 16)
                                    magneticSand->_spos.push_back((pos - VectorDd(14, 0, 0)) * magneticSand->_r), ++magneticSand->_n;
                            }
                    }
                }

                for (int i = 0; i < magneticSand->_n; ++i) {
                    VectorDd pos = magneticSand->_spos[i] + VectorDd(28, 0, 0) * magneticSand->_r;

                    bool flag = true;
                    for (int j = 0; j < magneticSand->_n && flag; ++j)
                        if ((pos - magneticSand->_spos[j]).norm() < 2 * magneticSand->_r)
                            flag = false;

                    if (flag)
                        magneticSand->_spos.push_back(pos);
                }

                magneticSand->_n = magneticSand->_spos.size();

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + 1.1 * VectorDd::Unit(0));
                    magneticSand->_spos[i] -= 1.1 * VectorDd::Unit(0);
                }
                magneticSand->_n *= 2;

                magneticSand->_ballCenter.push_back(VectorDd(-1.1, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(1.1, 0, 0));

                magneticSand->init_M();
            }
        }

        // two multi-chain
        template<int Dim>
        static void buildCase8(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);
                MatrixDd                                      px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                const double R = 1. / 15.;

                while (magneticSand->_n < 2000) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = nd(rng) * 2 - 1;
                    double   y = (nd(rng) * 2 - 1) * R;
                    double   z = (nd(rng) * 2 - 1) * R;

                    if (partion <= 2.) {
                        pos          = VectorDd(x, y, z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > R || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else if (partion <= 2. + std::sqrt(3)) {
                        pos          = VectorDd(x, y, z + 0.5);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > R || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x, y, z - 0.5);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > R || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    }

                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }

                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                magneticSand->_ballCenter.push_back(VectorDd(-1.25, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(1.25, 0, 0));
                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(px * magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos[i] = magneticSand->_spos[i] + magneticSand->_ballCenter[1];
                }
                magneticSand->_n *= 2;

                magneticSand->init_M();
            }
        }

        // two multi-disk
        template<int Dim>
        static void buildCase9(MagneticSand<Dim> * const magneticSand, const double ba) {
            DECLARE_DIM_TYPES(Dim);

            std::cout << ba << std::endl;

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double Z = 1. / 25.;

                while (magneticSand->_n < 4000) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = nd(rng) * 2 - 1;
                    double   y = nd(rng) * 2 - 1;
                    double   z = (nd(rng) * 2 - 1) * Z;

                    if (partion <= 2.) {
                        pos          = VectorDd(x, y, z);
                        VectorDd pro = VectorDd(x, y, 0);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else if (partion <= 2. + std::sqrt(3)) {
                        pos          = VectorDd(x, y, z + 0.5);
                        VectorDd pro = VectorDd(x, y, 0);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x, y, z - 0.5);
                        VectorDd pro = VectorDd(x, y, 0);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    }

                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + 1.1 * VectorDd::Unit(0));
                    magneticSand->_spos[i] -= 1.1 * VectorDd::Unit(0);
                }
                magneticSand->_n *= 2;

                magneticSand->_ballCenter.push_back(VectorDd(-1.1, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(1.1, 0, 0));

                magneticSand->init_M();
            }
        }

        // two ball
        template<int Dim>
        static void buildCase10(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                while (magneticSand->_n < 4000) {
                    double x = nd(rng) * 2 - 1;
                    double y = nd(rng) * 2 - 1;
                    double z = nd(rng) * 2 - 1;

                    VectorDd pos = VectorDd(x, y, z);

                    if (pos.norm() + magneticSand->_r > 1.)
                        continue;
                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + 1 * VectorDd::Unit(0));
                    magneticSand->_spos[i] -= 1 * VectorDd::Unit(0);
                }
                magneticSand->_n *= 2;

                magneticSand->_ballCenter.push_back(VectorDd(-1, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(1, 0, 0));

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase11(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double Z  = 1. / 25.;
                const double ba = 0.5;

                const double r = 8. / 15.;

                while (magneticSand->_n < 2500) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = nd(rng) * 2 - 1;
                    double   y = nd(rng) * 2 - 1;
                    double   z = (nd(rng) * 2 - 1) * Z;
                    if (partion <= 2.) {
                        pos          = VectorDd(x, ba * y, z);
                        VectorDd pro = VectorDd(x, y, 0);

                        if (pro.norm() + magneticSand->_r > 1.)
                            continue;
                        if ((pos - (-1. + r) * VectorDd::Unit(0)).norm() + magneticSand->_r > r && (pos - (1. - r) * VectorDd::Unit(0)).norm() + magneticSand->_r > r)
                            continue;
                    } else if (partion <= 2. + std::sqrt(3)) {
                        pos          = VectorDd(x, ba * y, z + 0.5 * r);
                        VectorDd pro = VectorDd(x, y, 0);
                        if (pro.norm() + magneticSand->_r > 1.)
                            continue;
                        if ((pos - (-1. + r) * VectorDd::Unit(0)).norm() + magneticSand->_r > r && (pos - (1. - r) * VectorDd::Unit(0)).norm() + magneticSand->_r > r)
                            continue;
                    } else {
                        pos          = VectorDd(x, ba * y, z - 0.5 * r);
                        VectorDd pro = VectorDd(x, y, 0);
                        if (pro.norm() + magneticSand->_r > 1.)
                            continue;
                        if ((pos - (-1. + r) * VectorDd::Unit(0)).norm() + magneticSand->_r > r && (pos - (1. - r) * VectorDd::Unit(0)).norm() + magneticSand->_r > r)
                            continue;
                    }

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }

                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }
                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + 1.1 * VectorDd::Unit(0));
                    magneticSand->_spos[i] -= 1.1 * VectorDd::Unit(0);
                }
                magneticSand->_n *= 2;

                magneticSand->_ballCenter.push_back(VectorDd(-1.1, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(1.1, 0, 0));

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase12(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                MatrixDd rotate;
                rotate << 0, -1, 0, 1, 0, 0, 0, 0, 1;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 500) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = (nd(rng) * 2 - 1) * X;
                    double   y = nd(rng) * 2 - 1;
                    double   z = nd(rng) * 2 - 1;

                    bool tag = false;

                    if (partion <= 2.) {
                        pos          = VectorDd(x, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                        tag = true;
                    }

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag) {
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                        if (tag) {
                            magneticSand->_spos.push_back(px * pos), ++magneticSand->_n;
                        }
                    }
                }

                double delta_x = 1;
                double delta_y = 1 * std::sqrt(3);

                magneticSand->_ballCenter.push_back(VectorDd(2 * delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-2 * delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, 0, 0));

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos[i] = rotate * magneticSand->_spos[i];
                }

                for (int j = 0; j < 6; ++j)
                    if (magneticSand->_ballCenter[j][0] < 0)
                        for (int i = 0; i < magneticSand->_n; ++i) {
                            magneticSand->_spos.push_back(px * magneticSand->_spos[i] + magneticSand->_ballCenter[j]);
                        }
                    else
                        for (int i = 0; i < magneticSand->_n; ++i) {
                            magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[j]);
                        }

                magneticSand->_n *= 7;

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase13(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                while (magneticSand->_n < 500) {
                    double x = nd(rng) * 2 - 1;
                    double y = nd(rng) * 2 - 1;
                    double z = nd(rng) * 2 - 1;

                    VectorDd pos = VectorDd(x, y, z);

                    if (pos.norm() + magneticSand->_r > 1.)
                        continue;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                double delta_x = 1;
                double delta_y = std::sqrt(3);

                magneticSand->_ballCenter.push_back(VectorDd(2 * delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-2 * delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, 0, 0));

                for (int j = 0; j < 6; ++j)
                    if (magneticSand->_ballCenter[j][0] < 0)
                        for (int i = 0; i < magneticSand->_n; ++i) {
                            magneticSand->_spos.push_back(px * magneticSand->_spos[i] + magneticSand->_ballCenter[j]);
                        }
                    else
                        for (int i = 0; i < magneticSand->_n; ++i) {
                            magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[j]);
                        }

                magneticSand->_n *= 7;

                magneticSand->init_M();
            }
        }

        // two plane
        template<int Dim>
        static void buildCase14(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 5000) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = (nd(rng) * 2 - 1) * X;
                    double   y = nd(rng) * 2 - 1;
                    double   z = nd(rng) * 2 - 1;

                    bool tag = false;

                    if (partion <= 2.) {
                        pos          = VectorDd(x, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                        tag = true;
                    }
                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag) {
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                        if (tag) {
                            magneticSand->_spos.push_back(px * pos), ++magneticSand->_n;
                        }
                    }
                }

                double delta_x = 0;
                double delta_y = 1.05;

                magneticSand->_ballCenter.push_back(VectorDd(delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, -delta_y, 0));

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos[i] = py * magneticSand->_spos[i] + magneticSand->_ballCenter[1];
                }
                magneticSand->_n *= 2;

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase15(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 3333) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = (nd(rng) * 2 - 1) * X;
                    double   y = nd(rng) * 2 - 1;
                    double   z = nd(rng) * 2 - 1;

                    bool tag = false;

                    if (partion <= 2.) {
                        pos          = VectorDd(x, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                        tag = true;
                    }
                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag) {
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                        if (tag) {
                            magneticSand->_spos.push_back(px * pos), ++magneticSand->_n;
                        }
                    }
                }

                double delta_x = 1.05;
                double delta_y = 1.05 / std::sqrt(3);

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(delta_x, -delta_y, 0));
                    magneticSand->_spos.push_back(px * magneticSand->_spos[i] + VectorDd(-delta_x, -delta_y, 0));
                    magneticSand->_spos[i] += VectorDd(0, 2 * delta_y, 0);
                }
                magneticSand->_n *= 3;

                magneticSand->_ballCenter.push_back(VectorDd(delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, 2 * delta_y, 0));

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase16(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 2500) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = (nd(rng) * 2 - 1) * X;
                    double   y = nd(rng) * 2 - 1;
                    double   z = nd(rng) * 2 - 1;

                    bool tag = false;

                    if (partion <= 2.) {
                        pos          = VectorDd(x, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                        tag = true;
                    }
                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag) {
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                        if (tag) {
                            magneticSand->_spos.push_back(px * pos), ++magneticSand->_n;
                        }
                    }
                }

                double delta_x = 2.1 / std::sqrt(2);
                double delta_y = 2.1 / std::sqrt(2);

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(py * magneticSand->_spos[i] + VectorDd(0, -delta_y, 0));
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(delta_x, 0, 0));
                    magneticSand->_spos.push_back(px * magneticSand->_spos[i] + VectorDd(-delta_x, 0, 0));
                    magneticSand->_spos[i] += VectorDd(0, delta_y, 0);
                }
                magneticSand->_n *= 4;

                magneticSand->_ballCenter.push_back(VectorDd(0, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, delta_y, 0));

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase17(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 2000) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = (nd(rng) * 2 - 1) * X;
                    double   y = nd(rng) * 2 - 1;
                    double   z = nd(rng) * 2 - 1;

                    bool tag = false;

                    if (partion <= 2.) {
                        pos          = VectorDd(x, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                        tag = true;
                    }
                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag) {
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                        if (tag) {
                            magneticSand->_spos.push_back(px * pos), ++magneticSand->_n;
                        }
                    }
                }

                double delta_x = 1.05;
                double delta_y = 1.05 / std::sqrt(3);

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(delta_x, delta_y, 0));
                    magneticSand->_spos.push_back(px * magneticSand->_spos[i] + VectorDd(-delta_x, delta_y, 0));
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(2 * delta_x, -2 * delta_y, 0));
                    magneticSand->_spos.push_back(px * magneticSand->_spos[i] + VectorDd(-2 * delta_x, -2 * delta_y, 0));
                    magneticSand->_spos[i] += VectorDd(0, -2 * delta_y, 0);
                }
                magneticSand->_n *= 5;

                magneticSand->_ballCenter.push_back(VectorDd(delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(2 * delta_x, -2 * delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-2 * delta_x, -2 * delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, -2 * delta_y, 0));

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase18(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 1666) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = (nd(rng) * 2 - 1) * X;
                    double   y = nd(rng) * 2 - 1;
                    double   z = nd(rng) * 2 - 1;

                    bool tag = false;

                    if (partion <= 2.) {
                        pos          = VectorDd(x, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                        tag = true;
                    }
                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag) {
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                        if (tag) {
                            magneticSand->_spos.push_back(px * pos), ++magneticSand->_n;
                        }
                    }
                }

                double delta_x = 1.05;
                double delta_y = 1.05 * std::sqrt(3);

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(delta_x, 0, 0));
                    magneticSand->_spos.push_back(px * magneticSand->_spos[i] + VectorDd(-delta_x, 0, 0));
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(0, delta_y, 0));
                    magneticSand->_spos.push_back(px * magneticSand->_spos[i] + VectorDd(-2 * delta_x, -delta_y, 0));
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(2 * delta_x, -delta_y, 0));
                    magneticSand->_spos[i] += VectorDd(0, -delta_y, 0);
                }
                magneticSand->_n *= 6;

                magneticSand->_ballCenter.push_back(VectorDd(delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-2 * delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(2 * delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, -delta_y, 0));

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase19(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                MatrixDd rotate;
                rotate << std::cos(kPi / 4), -std::sin(kPi / 4), 0, std::sin(kPi / 4), std::cos(kPi / 4), 0, 0, 0, 1;

                while (magneticSand->_n < 3333) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = (nd(rng) * 2 - 1) * X;
                    double   y = nd(rng) * 2 - 1;
                    double   z = nd(rng) * 2 - 1;

                    if (partion <= 2.) {
                        pos          = rotate * VectorDd(x, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else if (partion <= 2. + std::sqrt(3)) {
                        pos          = rotate * VectorDd(x + 0.5, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = rotate * VectorDd(x - 0.5, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    }
                    pos *= 0.9;
                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                magneticSand->_ballCenter.push_back(VectorDd(-2. / 5. * std::sqrt(2), -2. / 5. * std::sqrt(2), 0));
                magneticSand->_ballCenter.push_back(VectorDd(2. / 5. * std::sqrt(2), 2. / 5. * std::sqrt(2), 0));
                magneticSand->_ballCenter.push_back(VectorDd(13. / 10. * std::sqrt(2), -13. / 10. * std::sqrt(2), 0));

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[1]);
                    magneticSand->_spos[i] += magneticSand->_ballCenter[2];
                }
                magneticSand->_n *= 3;

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase20(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double R = 1. / 15.;

                while (magneticSand->_n < 1428) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = nd(rng) * 2 - 1;
                    double   y = (nd(rng) * 2 - 1) * R;
                    double   z = (nd(rng) * 2 - 1) * R;

                    if (partion <= 2.) {
                        pos          = VectorDd(x, y, z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > R || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else if (partion <= 2. + std::sqrt(3)) {
                        pos          = VectorDd(x, y + 0.5, z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > R || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x, y - 0.5, z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > R || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    }
                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }

                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                double delta_x = 1.05;
                double delta_y = 1.05 * std::sqrt(3);

                MatrixDd p;
                p << -1, 0, 0, 0, 1, 0, 0, 0, 1;

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(2 * delta_x, 0, 0));
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(delta_x, delta_y, 0));
                    magneticSand->_spos.push_back(p * magneticSand->_spos[i] + VectorDd(-delta_x, delta_y, 0));
                    magneticSand->_spos.push_back(p * magneticSand->_spos[i] + VectorDd(-2 * delta_x, 0, 0));
                    magneticSand->_spos.push_back(p * magneticSand->_spos[i] + VectorDd(-delta_x, -delta_y, 0));
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(delta_x, -delta_y, 0));
                }
                magneticSand->_n *= 7;

                magneticSand->_ballCenter.push_back(VectorDd(2 * delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-2 * delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, 0, 0));

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase21(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                MatrixDd p;
                p << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double R = 1. / 15.;

                while (magneticSand->_n < 1428) {
                    double x = nd(rng) * 2 - 1;
                    double y = nd(rng) * 2 - 1;
                    double z = nd(rng) * 2 - 1;

                    VectorDd pos = VectorDd(x, y, z);

                    if (pos.norm() + magneticSand->_r > 1.)
                        continue;
                    pos *= 0.9;
                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                double delta_x = 1.;
                double delta_y = 1. * std::sqrt(3);

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(2 * delta_x, 0, 0));
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(delta_x, delta_y, 0));
                    magneticSand->_spos.push_back(p * magneticSand->_spos[i] + VectorDd(-delta_x, delta_y, 0));
                    magneticSand->_spos.push_back(p * magneticSand->_spos[i] + VectorDd(-2 * delta_x, 0, 0));
                    magneticSand->_spos.push_back(p * magneticSand->_spos[i] + VectorDd(-delta_x, -delta_y, 0));
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + VectorDd(delta_x, -delta_y, 0));
                }
                magneticSand->_n *= 7;

                magneticSand->_ballCenter.push_back(VectorDd(2 * delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-2 * delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, 0, 0));

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase22(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double Z = 1. / 25.;
                while (magneticSand->_n < 8000) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = nd(rng) * 2 - 1;
                    double   y = nd(rng) * 2 - 1;
                    double   z = (nd(rng) * 2 - 1) * Z;

                    if (partion <= 2.) {
                        pos          = VectorDd(x, y, z);
                        VectorDd pro = VectorDd(x, y, 0);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else if (partion <= 2. + std::sqrt(3)) {
                        pos          = VectorDd(x, y, z + 0.5);
                        VectorDd pro = VectorDd(x, y, 0);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x, y, z - 0.5);
                        VectorDd pro = VectorDd(x, y, 0);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    }

                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }
                magneticSand->_ballCenter.push_back(VectorDd(0, 0, 0));

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase23(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                double rad = kPi * 75. / 180.;

                MatrixDd rotate;
                rotate << std::cos(rad), -std::sin(rad), 0, std::sin(rad), std::cos(rad), 0, 0, 0, 1;

                while (magneticSand->_n < 3333) {
                    double   partion = nd(rng) * (2. + 2 * std::sqrt(3));
                    VectorDd pos;
                    double   x = (nd(rng) * 2 - 1) * X;
                    double   y = nd(rng) * 2 - 1;
                    double   z = nd(rng) * 2 - 1;

                    if (partion <= 2.) {
                        pos          = rotate * VectorDd(x, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else if (partion <= 2. + std::sqrt(3)) {
                        pos          = rotate * VectorDd(x + 0.5, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = rotate * VectorDd(x - 0.5, y, 0.1 * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    }
                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                VectorDd t = VectorDd(2 * std::cos(kPi * 15 / 180), -2 * std::sin(kPi * 15 / 180), 0);

                magneticSand->_ballCenter.push_back(VectorDd(0, 0, 0));
                magneticSand->_ballCenter.push_back(t + 0.8 * VectorDd(std::cos(rad), std::sin(rad), 0));
                magneticSand->_ballCenter.push_back(t - 0.8 * VectorDd(std::cos(rad), std::sin(rad), 0));

                std::cout << magneticSand->_ballCenter[0].transpose() << std::endl;
                std::cout << magneticSand->_ballCenter[1].transpose() << std::endl;
                std::cout << magneticSand->_ballCenter[2].transpose() << std::endl;

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[1]);
                    magneticSand->_spos[i] += magneticSand->_ballCenter[2];
                }
                magneticSand->_n *= 3;

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase24(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 5000) {
                    double x = nd(rng) * 2 - 1;
                    double y = nd(rng) * 2 - 1;
                    double z = nd(rng) * 2 - 1;

                    VectorDd pos = VectorDd(x, y, z);

                    if (pos.norm() + magneticSand->_r > 1.)
                        continue;
                    pos *= 0.9;
                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                double delta_x = 0;
                double delta_y = 1.05;

                magneticSand->_ballCenter.push_back(VectorDd(delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, -delta_y, 0));

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos[i] = py * magneticSand->_spos[i] + magneticSand->_ballCenter[1];
                }
                magneticSand->_n *= 2;

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase25(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 3333) {
                    double x = nd(rng) * 2 - 1;
                    double y = nd(rng) * 2 - 1;
                    double z = nd(rng) * 2 - 1;

                    VectorDd pos = VectorDd(x, y, z);

                    if (pos.norm() + magneticSand->_r > 1.)
                        continue;
                    pos *= 0.9;
                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                double delta_x = 1.05;
                double delta_y = 1.05 / std::sqrt(3);

                magneticSand->_ballCenter.push_back(VectorDd(delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, 2 * delta_y, 0));

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos.push_back(px * magneticSand->_spos[i] + magneticSand->_ballCenter[1]);
                    magneticSand->_spos[i] = magneticSand->_spos[i] + magneticSand->_ballCenter[2];
                }
                magneticSand->_n *= 3;

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase26(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 2500) {
                    double x = nd(rng) * 2 - 1;
                    double y = nd(rng) * 2 - 1;
                    double z = nd(rng) * 2 - 1;

                    VectorDd pos = VectorDd(x, y, z);

                    if (pos.norm() + magneticSand->_r > 1.)
                        continue;
                    pos *= 0.9;
                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                double delta_x = 1.05;
                double delta_y = 1.05 * std::sqrt(3);

                magneticSand->_ballCenter.push_back(VectorDd(delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, 0, 0));

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[1]);
                    magneticSand->_spos.push_back(py * magneticSand->_spos[i] + magneticSand->_ballCenter[2]);
                    magneticSand->_spos[i] = px * magneticSand->_spos[i] + magneticSand->_ballCenter[3];
                }
                magneticSand->_n *= 4;

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase27(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 2500) {
                    double x = nd(rng) * 2 - 1;
                    double y = nd(rng) * 2 - 1;
                    double z = nd(rng) * 2 - 1;

                    VectorDd pos = VectorDd(x, y, z);

                    if (pos.norm() + magneticSand->_r > 1.)
                        continue;
                    pos *= 0.9;
                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                double delta_x = 2.1 / std::sqrt(2);
                double delta_y = 2.1 / std::sqrt(2);

                magneticSand->_ballCenter.push_back(VectorDd(delta_x, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(0, -delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(-delta_x, 0, 0));

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[1]);
                    magneticSand->_spos.push_back(py * magneticSand->_spos[i] + magneticSand->_ballCenter[2]);
                    magneticSand->_spos[i] = px * magneticSand->_spos[i] + magneticSand->_ballCenter[3];
                }
                magneticSand->_n *= 4;

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase28(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 2000) {
                    double x = nd(rng) * 2 - 1;
                    double y = nd(rng) * 2 - 1;
                    double z = nd(rng) * 2 - 1;

                    VectorDd pos = VectorDd(x, y, z);

                    if (pos.norm() + magneticSand->_r > 1.)
                        continue;
                    pos *= 0.9;
                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                VectorDd center = VectorDd(2.1 / 5, (2.1 - 2.1 * std::sqrt(3)) / 5, 0);

                magneticSand->_ballCenter.push_back(VectorDd(2.1, 0, 0) - center);
                magneticSand->_ballCenter.push_back(VectorDd(0, 2.1, 0) - center);
                magneticSand->_ballCenter.push_back(VectorDd(1.05, -1.05 * std::sqrt(3), 0) - center);
                magneticSand->_ballCenter.push_back(VectorDd(-1.05, -1.05 * std::sqrt(3), 0) - center);
                magneticSand->_ballCenter.push_back(VectorDd(0, 0, 0) - center);

                for (const auto & center : magneticSand->_ballCenter)
                    std::cout << center << std::endl
                              << std::endl;
                exit(0);

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[1]);
                    magneticSand->_spos.push_back(py * magneticSand->_spos[i] + magneticSand->_ballCenter[2]);
                    magneticSand->_spos.push_back(px * py * magneticSand->_spos[i] + magneticSand->_ballCenter[3]);
                    magneticSand->_spos[i] = magneticSand->_spos[i] + magneticSand->_ballCenter[4];
                }
                magneticSand->_n *= 5;

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase29(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 1666) {
                    double x = nd(rng) * 2 - 1;
                    double y = nd(rng) * 2 - 1;
                    double z = nd(rng) * 2 - 1;

                    VectorDd pos = VectorDd(x, y, z);

                    if (pos.norm() + magneticSand->_r > 1.)
                        continue;
                    pos *= 0.9;
                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                VectorDd center = VectorDd((-2.1 + 2.1 * std::cos(kPi / 9) + 2.1 * std::cos(2 * kPi / 9) - 2.1 * std::cos(7 * kPi / 180)) / 6, (-2.1 - 2.1 * std::sin(kPi / 9) + 2.1 * std::sin(2 * kPi / 9) + 2.1 * std::sin(7 * kPi / 180)) / 6, 0);

                magneticSand->_ballCenter.push_back(VectorDd(-2.1, 0, 0) - center);
                magneticSand->_ballCenter.push_back(VectorDd(0, -2.1, 0) - center);
                magneticSand->_ballCenter.push_back(VectorDd(2.1 * std::cos(kPi / 9), -2.1 * std::sin(kPi / 9), 0) - center);
                magneticSand->_ballCenter.push_back(VectorDd(2.1 * std::cos(2 * kPi / 9), 2.1 * std::sin(2 * kPi / 9), 0) - center);
                magneticSand->_ballCenter.push_back(VectorDd(-2.1 * std::cos(7 * kPi / 18), 2.1 * std::sin(7 * kPi / 18), 0) - center);
                magneticSand->_ballCenter.push_back(VectorDd(0, 0, 0) - center);

                for (const auto & center : magneticSand->_ballCenter)
                    std::cout << center << std::endl
                              << std::endl;
                exit(0);

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(px * magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos.push_back(py * magneticSand->_spos[i] + magneticSand->_ballCenter[1]);
                    magneticSand->_spos.push_back(py * magneticSand->_spos[i] + magneticSand->_ballCenter[2]);
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[3]);
                    magneticSand->_spos.push_back(px * magneticSand->_spos[i] + magneticSand->_ballCenter[4]);
                    magneticSand->_spos[i] = magneticSand->_spos[i] + magneticSand->_ballCenter[5];
                }
                magneticSand->_n *= 6;

                magneticSand->init_M();
            }
        }

        // single bundle
        template<int Dim>
        static void buildCase30(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3)

            {
                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);
                MatrixDd                                      px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                const double R = 1. / 5.;

                while (magneticSand->_n < 4000) {
                    VectorDd pos;
                    double   x = nd(rng) * 2 - 1;
                    double   y = (nd(rng) * 2 - 1) * R;
                    double   z = (nd(rng) * 2 - 1) * R;

                    pos          = VectorDd(x, y, z);
                    VectorDd pro = VectorDd(0, y, z);
                    if (pro.norm() + magneticSand->_r > R || pos.norm() + magneticSand->_r > 1.)
                        continue;

                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }

                    if (flag)
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                }

                magneticSand->_ballCenter.push_back(VectorDd(-1.1, 0, 0));
                magneticSand->_ballCenter.push_back(VectorDd(1.1, 0, 0));
                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(px * magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos[i] = magneticSand->_spos[i] + magneticSand->_ballCenter[1];
                }
                magneticSand->_n *= 2;

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCase31(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                MatrixDd px, py;
                px << -1, 0, 0, 0, 1, 0, 0, 0, 1;
                py << 1, 0, 0, 0, -1, 0, 0, 0, 1;

                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 25.;

                while (magneticSand->_n < 915) {
                    VectorDd pos;
                    double   x = (nd(rng) * 2 - 1) * X;
                    double   y = nd(rng) * 2 - 1;
                    double   z = nd(rng) * 2 - 1;

                    pos          = VectorDd(x, y, 0.1 * z);
                    VectorDd pro = VectorDd(0, y, z);
                    if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                        continue;

                    pos *= 0.9;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag) {
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                    }
                }

                double delta_x = 0;
                double delta_y = 1.05;

                magneticSand->_ballCenter.push_back(VectorDd(delta_x, delta_y, 0));
                magneticSand->_ballCenter.push_back(VectorDd(delta_x, -delta_y, 0));

                for (int i = 0; i < magneticSand->_n; ++i) {
                    magneticSand->_spos.push_back(magneticSand->_spos[i] + magneticSand->_ballCenter[0]);
                    magneticSand->_spos[i] = py * magneticSand->_spos[i] + magneticSand->_ballCenter[1];
                }
                magneticSand->_n *= 2;

                magneticSand->init_M();
            }
        }

        template<int Dim>
        static void buildCasebd(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            std::cout << "Begin Read.\n";

            std::ifstream fin("C:\\Users\\zhuyu\\Desktop\\bd\\a.txt");
            char          c;
            for (int j = 999; j >= 0; --j)
                for (int i = 0; i < 1000; ++i)
                    fin.get(c), magneticSand->image[i][j] = c == '1';

            std::cout << "End Read.\n";
            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 400.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double X = 1. / 50.;

                while (magneticSand->_n < 10000) {
                    VectorDd pos;
                    double   x = nd(rng);
                    double   y = nd(rng);
                    double   z = (nd(rng) * 2 - 1) * X;

                    pos = VectorDd(x, y, z);

                    if (! magneticSand->image[(int) (x * 1000)][(int) (y * 1000)])
                        continue;

                    bool flag = true;
                    for (const auto & sp : magneticSand->_spos)
                        if ((pos - sp).norm() < 2 * magneticSand->_r) {
                            flag = false;
                            break;
                        }
                    if (flag) {
                        magneticSand->_spos.push_back(pos), ++magneticSand->_n;
                    }
                }
                std::cout << "End Generate.\n";

                magneticSand->init_M();
            }
        }
    };
} // namespace PhysX
