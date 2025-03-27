#ifndef CRC
#define CRC 

#include "settings.h"

#include <bitset>
#include <vector>

class crc
{
private:
    std::string poly_string;
    std::bitset<_K_> polynome;

    std::bitset<_K_> code_book[(1 << _M_)];

    uint64_t d;
    std::vector<uint64_t> distances;

    std::bitset<_K_> get_check_part(const std::bitset<_M_> &msg);
    std::bitset<_K_> get_syndrom(const std::bitset<_K_> &code_word);

    uint64_t binomial_coeff(uint64_t n, uint64_t k);

public:
    crc();
    ~crc();

    std::bitset<_K_> code(const std::bitset<_M_> &msg);
    std::bitset<_K_> code(const uint64_t index);

    std::pair<std::bitset<_M_>, bool> decode(const std::bitset<_K_> &code_word);

    double upper_border_err(double ch_prob);
    double real_err(double ch_prob);
};

#endif // !CRC