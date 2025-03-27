#include "../include/crc.h"

#include <cmath>

//std::string crc::poly_string(_POLY_);
//std::bitset<_K_> crc::polynome(crc::poly_string + std::string(_K_ - poly_string.size(), '0'));

crc::crc() :poly_string(_POLY_), d(_K_), distances(_K_ + 1)
{
    polynome = std::bitset<_K_>(poly_string + std::string(_K_ - poly_string.size(), '0'));

    std::bitset<_M_> msg;
    std::bitset<_K_> code_word;

    uint64_t weight(0);

    for (uint64_t i(0); i < (1 << _M_); ++i)
    {
        msg = std::bitset<_M_>(i);
        code_word = code(msg);

        weight = code_word.count();
        ++distances[weight];

        if ((i != 0) && (d > weight))
        {
            d = weight;
        }

        code_book[i] = code_word;
    }
}

crc::~crc()
{
}

std::bitset<_K_> crc::get_check_part(const std::bitset<_M_> &msg)
{
    std::bitset<_K_> remainder(msg.to_string());
    remainder <<= _R_;

    for (uint64_t i = 0; i < _M_; ++i) {
        if (remainder.test(_K_ - 1)) {
            remainder ^= polynome;
        }
        remainder <<= 1;
    }

    remainder >>= _M_;

    return remainder;
}

std::bitset<_K_> crc::get_syndrom(const std::bitset<_K_> &code_word)
{
    std::bitset<_K_> syndrome(code_word); 

    for (uint64_t i = 0; i < _M_; ++i) {
        if (syndrome.test(_K_ - 1)) {
            syndrome ^= polynome;
        }
        syndrome <<= 1;
    }

    //syndrome >>= _R_;

    return syndrome;
}

uint64_t crc::binomial_coeff(uint64_t n, uint64_t k)
{
    if (k > n)
    {
        return 0;
    }

    if (k == 0 || k == n) 
    {
        return 1;
    }

    uint64_t result = 1;
    for (uint64_t i = 1; i <= k; ++i) {
        result = result * (n - (k - i)) / i;
    }
    return result;
}

std::bitset<_K_> crc::code(const std::bitset<_M_> &msg) 
{
    std::bitset<_K_> output(msg.to_string());
    output <<= _R_;

    std::bitset<_K_> remainder = get_check_part(msg);

    output |= remainder;

    return output;
}

std::bitset<_K_> crc::code(const uint64_t index)
{
    if (index >= (1 << _M_))
    {
        return code_book[0];
    }

    return code_book[index];
}

std::pair<std::bitset<_M_>, bool> crc::decode(const std::bitset<_K_> &code_word) 
{

    std::bitset<_M_> msg;
    for (int i = 0; i < _M_; ++i) {
        msg[i] = code_word[i + _R_];
    }

    bool is_valid = false;

    if (get_syndrom(code_word) == 0)
    {
        is_valid = true;
    }

    return std::make_pair(msg, is_valid);
}

double crc::upper_border_err(double ch_prob)
{
    double sum = 0.0;

    double ez1 = 0.0;
    double ez2 = 0.0;
    double ez3 = 0.0;

    for (uint64_t i(0); i < d; ++i)
    {
        ez1 = static_cast<double>(binomial_coeff(_K_, i));
        ez2 = ez1 * std::pow(ch_prob, i);
        ez3 = ez2 * std::pow(1.0 - ch_prob, (_K_ - i));

        sum += ez3;
        //sum += (binomial_coeff(_K_, i) * std::pow(ch_prob, i) * std::pow(1.0 - ch_prob, (_K_ - i)));
    }

    return 1.0 - sum;
}

double crc::real_err(double ch_prob)
{
    double sum = 0.0;

    double ez1 = 0.0;
    double ez2 = 0.0;

    for (uint64_t i(d); i < _K_ + 1; ++i)
    {
        ez1 = static_cast<double>(distances[i]) * std::pow(ch_prob, i);
        ez2 = ez1 * std::pow(1.0 - ch_prob, (_K_ - i));

        sum += ez2;
        //sum += (distances[i] * std::pow(ch_prob, i) * std::pow(1.0 - ch_prob, (_K_ - i)));
    }

    return sum;
}
