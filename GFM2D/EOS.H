#include <iostream>
#include <vector>
#include <array>

class EOS {

    protected:
        double gamma;

    public:

        enum enuma{RHO = 0, XMOM, YMOM, ENE};
        enum enumb{BASE0 = 0, UX, UY, PRES};
        enum enumc{YDIR=0,XDIR=1};

        virtual double calcSoundSpeed(const  std::array<double,4>) = 0;
        virtual std::array<double,4> primToConsv(const std::array<double,4>)const = 0;
        virtual std::array<double,4> consvToPrim(const std::array<double,4>)const = 0;
        double get_gamma();

        EOS(double);

        virtual ~EOS() = default;
};

class idealGas : public EOS {

    public:
        double calcSoundSpeed(const std::array<double,4>) override;
        std::array<double,4> primToConsv(const std::array<double,4>)const override;
        std::array<double,4> consvToPrim(const std::array<double,4>)const override;

        idealGas(double);

};

class stiffenedGas : public EOS {
    private:
        double p_inf;
        double e_inf;

    public:
        double calcSoundSpeed(const std::array<double,4>) override;
        std::array<double,4> primToConsv(const std::array<double,4>)const override;
        std::array<double,4> consvToPrim(const std::array<double,4>)const override;

        stiffenedGas(double,double,double e_inf = 0);
};