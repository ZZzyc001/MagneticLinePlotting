#pragma once

#include "Simulation.h"

#include <vector>

namespace PhysX {

    template<int Dim>
    class MagSource {
    public:
        DECLARE_DIM_TYPES(Dim)

        virtual VectorDd field_at_position(const VectorDd & pos) const = 0;
    };

    template<int Dim>
    class UniformSource : public MagSource<Dim> {
    public:
        DECLARE_DIM_TYPES(Dim)

        VectorDd _dir;
        double   _mag;

        UniformSource(const VectorDd & dir, double mag):
            _dir(dir), _mag(mag) {
            if (dir.norm() < 1e-12)
                throw std::runtime_error("Direction vector is too small");

            _dir /= _dir.norm();
        }

        VectorDd field_at_position(const VectorDd & pos) const override;
    };

    template<int Dim>
    class DipoleSource : public MagSource<Dim> {
    public:
        DECLARE_DIM_TYPES(Dim)

        VectorDd _pos;
        VectorDd _dipole;

        DipoleSource() {}

        DipoleSource(const VectorDd & pos, const VectorDd & dipole):
            _pos(pos), _dipole(dipole) {}

        VectorDd field_at_position(const VectorDd & pos) const override;
    };

    template<int Dim>
    class MagneticSand : public Simulation {
    public:
        DECLARE_DIM_TYPES(Dim)
        friend class MagneticSandBuilder;

        std::vector<VectorDd>                        _pos;
        std::vector<bool>                            _fixed;
        std::vector<VectorDd>                        _spos;
        std::vector<VectorDd>                        _M;
        std::vector<std::unique_ptr<MagSource<Dim>>> _sources;

        std::vector<VectorDd> _ballCenter;

        std::vector<VectorDd>              _sample_pos;
        std::vector<VectorDd>              s_sample_pos;
        std::vector<std::vector<VectorDd>> _history;

        bool image[1010][1010] {};

        MatrixXd A;
        VectorXd b;

        double torque1[40][20] {};
        double torque2[40][20] {};
        double forcerr[40][20] {};
        double torqua1[40][20] {};
        double torqua2[40][20] {};
        double forcear[40][20] {};

        int    _n = 0;
        double _r = 0.01;

        const double _eps = 1e-6;

        // _mu_r = mu / mu_0
        const double _mu_r = 1000;

        MagneticSand() {}

        virtual int  dimension() const override { return Dim; }
        virtual void writeDescription(YAML::Node & root) const override;
        virtual void writeFrame(const std::string & frameDir, const bool staticDraw) const override;
        virtual void saveFrame(const std::string & frameDir) const override;
        virtual void loadFrame(const std::string & frameDir) override;

        virtual void initialize() override;
        virtual void advance(const double dt) override;

        void init_M();

        void cal_mag_linearly(int number);

        VectorDd cal_field_dir_at_position(const VectorDd & pos) const;
        VectorDd cal_field_at_position(const VectorDd & pos) const;
        double   cal_field_norm_at_position(const VectorDd & pos) const;
        VectorDd cal_force(int i) const;

        bool collide2particle(const VectorDd & pos) const;

        void printTorque();

        void generateSampleParticle();
    };

} // namespace PhysX
