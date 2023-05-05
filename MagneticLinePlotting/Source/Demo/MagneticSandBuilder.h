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
                buildCase1(magneticSand.get(), ba);
                break;
            case 2:
                buildCase2(magneticSand.get());
                break;
            case 3:
                buildCase3(magneticSand.get(), ba);
                break;
            case 4:
                buildCase4(magneticSand.get(), ba);
                break;
            case 5:
                buildCase5(magneticSand.get());
                break;
            case 6:
                buildCase6(magneticSand.get(), ba);
                break;
            case 7:
                buildCase7(magneticSand.get(), ba);
                break;
            case 8:
                buildCase8(magneticSand.get(), ba);
                break;
            case 9:
                buildCase9(magneticSand.get(), ba);
                break;
            case 10:
                buildCase10(magneticSand.get(), ba);
                break;
            case 11:
                buildCase11(magneticSand.get(), ba);
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
                buildCase15(magneticSand.get(), ba);
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
            default:
                return nullptr;
            }

            return magneticSand;
        }

    protected:
        // two multi-chain
        template<int Dim>
        static void buildCase0(MagneticSand<Dim> * const magneticSand) {
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
        static void buildCase1(MagneticSand<Dim> * const magneticSand, const double ba) {
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
        static void buildCase2(MagneticSand<Dim> * const magneticSand) {
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
        static void buildCase3(MagneticSand<Dim> * const magneticSand, double ba) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
                magneticSand->_r = 1. / 100.;

                static std::mt19937                           rng(42);
                static std::uniform_real_distribution<double> nd(0, 1);

                const double Z = 1. / 25.;

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
        static void buildCase4(MagneticSand<Dim> * const magneticSand, double ba) {
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
                        pos          = VectorDd(x, y, ba * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, ba * z);
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
        static void buildCase5(MagneticSand<Dim> * const magneticSand) {
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
        static void buildCase6(MagneticSand<Dim> * const magneticSand, double ba) {
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
                        pos          = VectorDd(x, y, ba * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, ba * z);
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
        static void buildCase7(MagneticSand<Dim> * const magneticSand, double ba) {
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
                        pos          = VectorDd(x, y, ba * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, ba * z);
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
        static void buildCase8(MagneticSand<Dim> * const magneticSand, double ba) {
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
                        pos          = VectorDd(x, y, ba * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, ba * z);
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
        static void buildCase9(MagneticSand<Dim> * const magneticSand, double ba) {
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
                        pos          = VectorDd(x, y, ba * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, ba * z);
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
        static void buildCase10(MagneticSand<Dim> * const magneticSand, double ba) {
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
                        pos          = VectorDd(x, y, ba * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = VectorDd(x + 0.5, y, ba * z);
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
        static void buildCase11(MagneticSand<Dim> * const magneticSand, double ba) {
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
                        pos          = rotate * VectorDd(x, y, ba * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else if (partion <= 2. + std::sqrt(3)) {
                        pos          = rotate * VectorDd(x + 0.5, y, ba * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = rotate * VectorDd(x - 0.5, y, ba * z);
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
        static void buildCase12(MagneticSand<Dim> * const magneticSand) {
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
        static void buildCase13(MagneticSand<Dim> * const magneticSand) {
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
        static void buildCase14(MagneticSand<Dim> * const magneticSand) {
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
        static void buildCase15(MagneticSand<Dim> * const magneticSand, double ba) {
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
                        pos          = rotate * VectorDd(x, y, ba * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else if (partion <= 2. + std::sqrt(3)) {
                        pos          = rotate * VectorDd(x + 0.5, y, ba * z);
                        VectorDd pro = VectorDd(0, y, z);
                        if (pro.norm() + magneticSand->_r > 1. || pos.norm() + magneticSand->_r > 1.)
                            continue;
                    } else {
                        pos          = rotate * VectorDd(x - 0.5, y, ba * z);
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
        static void buildCase16(MagneticSand<Dim> * const magneticSand) {
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
        static void buildCase17(MagneticSand<Dim> * const magneticSand) {
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
        static void buildCase18(MagneticSand<Dim> * const magneticSand) {
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
        static void buildCase19(MagneticSand<Dim> * const magneticSand) {
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
        static void buildCase20(MagneticSand<Dim> * const magneticSand) {
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
        static void buildCase21(MagneticSand<Dim> * const magneticSand) {
            DECLARE_DIM_TYPES(Dim);

            if constexpr (Dim == 3) {
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
    };
} // namespace PhysX
