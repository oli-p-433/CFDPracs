
#include <iostream>
#include <vector>
#include <array>

class EOS {
    protected:
        double gamma;
        double BX;

    public:
        enum Cons{RHO = 0, XMOM, YMOM, ZMOM, ENE, BY, BZ};
        enum Cons2{BASE0 = 0, UX, UY, UZ, PRES};

        virtual double calcSoundSpeed(const  std::array<double,7>) = 0;
        virtual std::array<double,7> primToConsv(const std::array<double,7>)const = 0;
        virtual std::array<double,7> consvToPrim(const std::array<double,7>)const = 0;

        EOS(double,double);

        virtual ~EOS() = default;
};

class idealGas : public EOS {

    public:
        double calcSoundSpeed(const std::array<double,7>) override;
        std::array<double,7> primToConsv(const std::array<double,7>)const override;
        std::array<double,7> consvToPrim(const std::array<double,7>)const override;

        idealGas(double,double);

};

class stiffenedGas : public EOS {
    private:
        double p_inf;
        double e_inf;

    public:
        double calcSoundSpeed(const std::array<double,7>) override;
        std::array<double,7> primToConsv(const std::array<double,7>)const override;
        std::array<double,7> consvToPrim(const std::array<double,7>)const override;

        stiffenedGas(double,double,double e_inf = 0, double BX = 0);
};