#ifndef CHANNEL
#define CHANNEL

#include <random>
#include <vector>

class channel
{
private:
    double mean = 0.0;
    double stddev = 0.1;

    std::default_random_engine generator;
    std::normal_distribution<double> *dist;

public:
    channel(double mean = 0.0, double stddev = 0.1);
    ~channel();

    double set_SNR_dB(double SNR_dB);
    void reset_prarms(double mean, double stddev);

    std::vector<double> add_noise(std::vector<double> &singals);
};

#endif // !CHANNEL