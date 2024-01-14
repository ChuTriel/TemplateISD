#ifndef SMALLSECRETLWE_DECODING_FILE_H
#define SMALLSECRETLWE_DECODING_FILE_H
#include <stdint.h>
#include <string>

constexpr uint64_t n = 2681;
constexpr uint64_t k = 2145;
constexpr uint64_t seed = 0;
constexpr uint64_t w = 45;
constexpr const char *s = "11011110000011001001000010111010010100111000000010001001110001000010010111101010001000110011101100010010110001001000101011110001000010110100110011001001000100000000110111010010011001111101100010100111111000000110101011000010100101101101111100011011001101011110101110100110100100110111100111010110101101101000101010000100001001000110001001011110011100111010011111111101011001101110100110001111101000011101001011111011010011110011000101110001100110110000000100010101011001100010001010000101010011010000011111110110000111100010010010100100";

constexpr const char* e_own = "00000000000000000011000000000000000000000000000000000000000000000100000000000001000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000100000000000000000000000000000000000000000000000000010000000000000000000000100000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000010000000000000000000000000000000000100000000000000010000000000000000000000000000000000000000001100000000000000000000000000000000000000000000001000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000001000000100000000000000000000000000000000000000000000010000000000000000000000100000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000100000010000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000001000000000000000000000000000000000100000000000000000000000010000000000000000000000000000000000000000000";
const std::string eW = "202010100002021000011001000101010110220110000000001010110100111111000000120110011110";
constexpr uint32_t newN = 1216;
constexpr uint32_t addRows = 38;

#endif //SMALLSECRETLWE_DECODING_FILE_H