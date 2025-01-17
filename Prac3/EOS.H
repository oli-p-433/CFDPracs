
#include <iostream>
#include <vector>
#include <array>

class EOS {
    protected:
        double gamma;

    public:
        virtual double calcSoundSpeed(const  std::array<double,3>) = 0;
        virtual std::array<double,3> primToConsv(const std::array<double,3>)const = 0;
        virtual std::array<double,3> consvToPrim(const std::array<double,3>)const = 0;

        EOS(double);

        virtual ~EOS() = default;
};

class idealGas : public EOS {

    public:
        double calcSoundSpeed(const std::array<double,3>) override;
        std::array<double,3> primToConsv(const std::array<double,3>)const override;
        std::array<double,3> consvToPrim(const std::array<double,3>)const override;

        idealGas(double);

};

class stiffenedGas : public EOS {
    private:
        double p_inf;
        double e_inf;

    public:
        double calcSoundSpeed(const std::array<double,3>) override;
        std::array<double,3> primToConsv(const std::array<double,3>)const override;
        std::array<double,3> consvToPrim(const std::array<double,3>)const override;

        stiffenedGas(double,double,double e_inf = 0);
};