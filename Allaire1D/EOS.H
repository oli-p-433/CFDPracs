
#include <iostream>
#include <vector>
#include <array>

class EOS {
    protected:
        double gamma1;
        double gamma2;

    public:
        virtual double calcSoundSpeed(const  std::array<double,5>&) = 0;
        virtual std::array<double,5> primToConsv(const std::array<double,5>&)const = 0;
        virtual std::array<double,5> consvToPrim(const std::array<double,5>&)const = 0;
        std::array<double,2> get_gamma();

        virtual std::array<double,2> get_pInf() const {
            return std::array<double,2>{0.0,0.0}; // or throw std::runtime_error("pInf not available for this EOS type");
        }

        EOS(double,double);

        virtual ~EOS() = default;
};

class idealGas : public EOS {

    public:
        double calcSoundSpeed(const std::array<double,5>&) override;
        std::array<double,5> primToConsv(const std::array<double,5>&)const override;
        std::array<double,5> consvToPrim(const std::array<double,5>&)const override;

        idealGas(double,double);

};

class stiffenedGas : public EOS {
    private:
        double p_inf1, p_inf2;
        double e_inf1, e_inf2;

    public:
        double calcSoundSpeed(const std::array<double,5>&) override;
        std::array<double,5> primToConsv(const std::array<double,5>&)const override;
        std::array<double,5> consvToPrim(const std::array<double,5>&)const override;

        std::array<double,2> get_pInf() const override {
            return std::array<double,2>{p_inf1,p_inf2};
        }

        stiffenedGas(double,double,double,double,double e_inf1 = 0, double e_inf2 = 0);

        std::array<double,2> get_pInf();
};