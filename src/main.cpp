#include "include/crc.h"

#include "include/crc.h"
#include "include/channel.h"
#include "include/file_tools.h"

#include <iostream>

////////////////////////////////////////////////

const std::string berr_file = "../data/berr.txt";
const std::string derr_file = "../data/derr.txt";
const std::string t_berr_file = "../data/t_berr.txt";
const std::string t_derr_file = "../data/t_derr.txt";
const std::string t_upper_derr_file = "../data/t_upper_derr.txt";

////////////////////////////////////////////////

const double stddev = 15.0;
const int seed = 10;
const uint64_t iterations = 1000000;
const double min_db = -20.0;
const double max_db = 1.0;
const double duration = 1.0;
const double width = (min_db < 0.0) ? ((max_db < 0.0) ? -min_db + -max_db : -min_db + max_db) : min_db + max_db;

const uint64_t capacity = static_cast<uint64_t>((1.0 / duration) * width);

inline void check_file_procedure();

inline uint64_t msg_generator();

inline uint64_t count_errors(const std::bitset<_K_> &tx, const std::bitset<_K_> &rx);

inline double count_berr(const uint64_t &errors);

inline double count_derr(const uint64_t &decoder_errors);

std::vector<double> bits_to_signal(std::bitset<_K_> &bits);

std::bitset<_K_> signal_to_bits(std::vector<double> &signals);

////////////////////////////////////////////////

uint64_t real_iterations(0);

int main(/*int argc, char const *argv[]*/)
{

    check_file_procedure();

    ////////////////////////////////////////////////
    srand(seed);

    crc crc_inst;
    channel channel_inst(0.0, stddev);

    ////////////////////////////////////////////////

    std::bitset<_M_> msg(0);
    std::bitset<_K_> code_word(0);

    std::bitset<_K_> demod_word(0);
    std::pair<std::bitset<_M_>, bool> decoder_output;

    std::vector<double> signal;

    ////////////////////////////////////////////////

    std::vector<double> berr;
    std::vector<double> derr;

    std::vector<double> t_berr;
    std::vector<double> t_upper_derr;
    std::vector<double> t_derr;

    berr.reserve(capacity);
    derr.reserve(capacity);

    t_berr.reserve(capacity);
    t_upper_derr.reserve(capacity);
    t_derr.reserve(capacity);

    uint64_t err_bits(0);
    uint64_t err_decoder(0);

    ////////////////////////////////////////////////

    double berr_prob(0.0);

    for (double Ydb(min_db); Ydb < max_db; Ydb += duration)
    {
        err_bits = 0;
        err_decoder = 0;

        berr_prob = 0.0;

        t_berr.push_back
        (
            channel_inst.set_SNR_dB(Ydb)
        );

        for (real_iterations = 0; (real_iterations < iterations) || (err_decoder == 0); ++real_iterations)
        {
            msg = std::bitset<_M_>(msg_generator());

            code_word = crc_inst.code(msg.to_ullong());

            signal = bits_to_signal(code_word);
            signal = channel_inst.add_noise(signal);

            demod_word = signal_to_bits(signal);
            decoder_output = crc_inst.decode(demod_word);

            err_bits = count_errors(code_word, demod_word);
            berr_prob += count_berr(err_bits);

            if (decoder_output.second)
            {
                if ((msg ^ decoder_output.first) != 0)
                {
                    ++err_decoder;
                }
            }
        }

        /*
        if (err_decoder == 0)
        {
            err_decoder = 1;
        }
        */

        berr.push_back(berr_prob / (double) real_iterations);
        derr.push_back(count_derr(err_decoder));
    }

    for (auto &val : t_berr)
    {
        t_upper_derr.push_back(crc_inst.upper_border_err(val));
        t_derr.push_back(crc_inst.real_err(val));
    }

    write_vector_file(berr_file, berr);
    write_vector_file(derr_file, derr);
    write_vector_file(t_berr_file, t_berr);
    write_vector_file(t_derr_file, t_derr);
    write_vector_file(t_upper_derr_file, t_upper_derr);

    return 0;
}

inline void check_file_procedure()
{
    check_for_file(berr_file.c_str());
    check_for_file(derr_file.c_str());
    check_for_file(t_berr_file.c_str());
    check_for_file(t_derr_file.c_str());
    check_for_file(t_upper_derr_file.c_str());
}

inline uint64_t msg_generator()
{
    return rand() % (1 << _M_);
}

inline uint64_t count_errors(const std::bitset<_K_> &tx, const std::bitset<_K_> &rx)
{
    return (tx ^ rx).count();
}

inline double count_berr(const uint64_t &errors)
{
    return static_cast<double>(errors) / static_cast<double>(_K_);
}

inline double count_derr(const uint64_t &decoder_errors)
{
    return static_cast<double>(decoder_errors) / static_cast<double>(real_iterations);
}

std::vector<double> bits_to_signal(std::bitset<_K_> &bits)
{
    std::vector<double> signals;
    for (uint64_t i = 0; i < _K_; ++i)
    {
        double val = (bits[i] == 0) ? -1.0 : 1.0;
        signals.push_back(val);
    }

    return signals;
}

std::bitset<_K_> signal_to_bits(std::vector<double> &signals)
{
    std::bitset<_K_> bits;
    for (uint64_t i = 0; i < _K_; ++i)
    {
        bits[i] = (signals[i] < 0.0) ? 0 : 1;
    }

    return bits;
}