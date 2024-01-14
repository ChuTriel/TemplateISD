#ifndef SMALLSECRETLWE_DECODING_FILE_H
#define SMALLSECRETLWE_DECODING_FILE_H
#include <stdint.h>
#include <string>

constexpr uint64_t n = 3108;
constexpr uint64_t k = 2487;
constexpr uint64_t seed = 0;
constexpr uint64_t w = 52;
constexpr const char *s = "000101111100001100101110100100100010001100001100010100010101111010100010111000100110001110010010000010000010101101111100111010111100110000111000000001100011110000010011011010110100011010011010110011001100100001110101011100000001011100010101111110000101010001101010010011001110110110011111100101000111101010110001011011100111010011110100000001001011001101111101011100110110100010000000100100100010110010001011011101011010000000111111010101010101011100000101111111011101001000010100011010011010001110011100001100100110111100100000100001000001000001001001101011100010011111011101101010000101100100110001000011110011101001000";

constexpr const char* e_own = "000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000100000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000011000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000100000000000000000000000000010000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000001000000000000000000000000000000000000000000010000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000010000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000010000000000000000000000000000000000000010000000000000000000000000000000000000000000001000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000010000000000000100000000000000000000000000000000000000000100000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000100000000000000000010000000000000000000000000000000000000001000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000001000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000";
const std::string eW = "01011010020102001210001012001010001111000100101110110110201100101100010202001011010000200100000100";
constexpr uint32_t newN = 1408;
constexpr uint32_t addRows = 44;

#endif //SMALLSECRETLWE_DECODING_FILE_H