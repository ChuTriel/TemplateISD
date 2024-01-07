#ifndef SMALLSECRETLWE_DECODING_FILE_H
#define SMALLSECRETLWE_DECODING_FILE_H
#include <cstdint>
#include <string>
constexpr uint64_t n = 381;
constexpr uint64_t k = 305;
constexpr uint64_t seed = 305;
constexpr uint64_t w = 9;
constexpr const char *s = "0101100000001000110100011101000111010100000001011101100110011101110110111111";
constexpr const char *h = "100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000110101111010110001000001100000000000110011110000001010101010001000001010111101101010011010110101000110111111100110110010000001010110100011000111001011110101110010011000000110011011111101001011110111010001001010111011010101110010011110100000101111011011001011111111100111001001010111000010010111101101010011001110000011111011010101101001101010100101000010011100001010011100010101101011011000001001101011101000001010001010000101000101010100001000010111111100111000001010011110001010101010011010011101010001011000011111100000011000101111011100110101010101011101011111001100000101101101100000001000001100011010000100011010111001111111000110110010011100000010100000100000110001010101000000000111001110110011111011101111010111001110110100100010110001111000000001111100110100110111011001000111011101001011101111001111100010010000001100000111100001110111011000110000111101110110101011000100110110000010001111101000001011010111011010111111011010110100101101101101111001101101111110111010001000011100000011010101001011101100111100000100111001000000101000011000001001111001011001101111001011011101101101100011100100100001101111100101000110100111111010111000011011001010010001110001001100100000000110011001111001000101001111010010110111111101000011001001001001110000000000111100110111111001011000100100010011000011010001010000110110011010010000000010011011100010000011101001001110011101111000111011000111011100101101010110010011101110110111110010101111111010011011001011010011010101000110111011100001101011111000111000000011010111100000001001010000000010101100011101001101000001111101001100111010111111000000110101000010111101010011100110010110010010100011111011111010010101100100100011001011001110110110011111010101111011000110010010001100010001000100101011000110101000111011111110001011110100110101110101111010110111011111100010001010101111100011000110011111010110011001011010000010100110000001001001010010001011101101110010100111110001110100101001110110100000110010010010111000011001110010011101000101101111110001100101101000011110000001100100011010000000100010000000010010000000000111000011110000101100111000101100011111011100100011100110100111101110000111100010010000110111011000000010000001100011100011100100011100011111001101100011101000000100100000101100010101111111110011001010010111010110111111010001111101001111111110111100001010010010100111100001000000001111000000000010001110001101110000101110111111110010100110011100111001000100110110011000001100011111100011110010110011101111100101100110101010101110100000110110101100011100100001101100011101111111000100110011101000101110111001010001000101000101111011010100001011101101001000111000011011100111011010111100001110010111000000000100110011100000010010110010010101000001110111000000001100011010100111110100110001101111110011110110011010000100100111111001110100011011101000101111010000000000000101010000010010110100111000001001011101110110101001010100100110010001110000000011000011110011101111111101010111101011110100101011101011110111010011111110110011010110010110111111111011000100101000110010000000000010011010001101111101100111000001111100110011000000111100110111011101111001110000111011000010001111101100100101010010000010111101100000111000000000101100001010111011110111011100011010111110100111111111110110011011100010110000110110001111000011100111010111101110100011010110100101111010110000111011100100111000101010100100001100000100110110011011111000101010100111011111111010000000111111110001010111111111110000010011101100000100010111101101000100101100000011101000011000011001111110110011000101111111100001110010000110111011001001101110000101011100000011000010111010011100110100111100111011110111001110110111001001100011110101110010111110001010100100001100101100001010010010100000011011010000011010000010100110111110110000111111111011001011111100101011010001011010101110010000011100011011001111110110101111101100001000011010100111111000011100101110101011110000011011101100010001011001010111111011100000111101010001101001000110101100110010001110110100110011000111001000011100011001010001010010110101111000110010000010001001111100101100001111111001111010001101011011100110111010100101000111111000001010011001011101010010011111100100111010101010010001101001111000001001010010010110000011011111000000110110001101001010011011010011000100110110100100101000110000100011111110001101010010010101000010111111010110110001011000110011000010000000111010101110101001010101010010000100010001111011000000101010101111011010111000100001111010000000011101100001011000001010001011110001101100110010011001010101100011010000010101000100100100110110101111001111000000101010000001011111101000101011000010000111011011101001001110100111111110011011011000010110010010001010100001001100011101001011111100111000101001101100001000011101010011000001111000001110010011011100111110110100011000100111101101110000010110110101000000000001111100011001101001010110011100000001111100001010111000011110011001100100011010111111000100010101000101100111101100100100100010000110011111011110101001000010001111011010110110000001100000101000011110111101001110100011111111100001111100100001010101100101100001010101111111001011101011001111110110100111110101111111011010111111010011111101011110010110101011000100111100110110111100111000001111011011001100010001011010000100000001111001010001111011110111011010001011110000101000111011111011110100100001101000100011011110010100001100000111101101110100110110001010101001100111100011011001001001110001100010001011110100011011010000010000001000100001010001101100100011110001001001000011101000111010001010010111110001001101110000000000101110111100001011100100101001110001010001010100101100100100001011000111000001100111110101100111110111100001111100001111101100000000101000101100000101100000101111110111010000001010010101011100000110000001100110011101000100111010001000001001010010100010010011100101010100011000101110100011000010111011100000101010001011100100110011101001001111011110011100101111010101000101001011010011111101011001011001010101001110010000001010011110001100001111111001010011011110100001000000101010100110011111100101111100101000011101100011011101111001100001011111100011010001101010111011111001110011000101110110111010010011111001010101000010110100111010011000010001011110000100011010001111010100011100010111100111110011010101011010110100001111010101101110010000100111111000110000000000101100100011001111000100000010000100100100000001111111101011111011100101100111011100110001000001111111011010001111010110001110110011100001011000001000000110001001000010010000010100101001010011110010101101011000011000010001111100100011101110011001110001110111100110101111001101111110101100101100010100010011011100100000111110011001110010010110100100000101010010111101111010011010100010001000001010110000110001100110101011111111011001110010001111001000011111110111011001001111110000101101011101110001001011110010010011100100100100001000101110000110011110101000000001010001000011011010001001100011111001010101001010110110101100001110110001000101101010100111101110111011110000110100010000011101100101100001011000110111110110100100001001000011000110110111101000010110100011010100111001100001010111110101010010001111100100110111111100000111110100000101100001010110111111100100100100001110000100110001010101001100010111011001000000100111110110000101110111010001000110010011111001100010001010110110110011110101010110001001100010000011010110001001000011010110101000101101010000011001000001110100010001111100100100010100110110100010101111001101110010111101001101000001001101101010000111011101011100011111001000110110111010010000001000001111000000101110011011100110101000111011110110100000111111001000010101100110011111100100010111110001000111001101101010110110010010011000001101110001100101011100000001110101011000111100101011111101001011010010010110000000010010110010111001001000001011111000010110010000110011111010111011001000011000110011100010010010110110001001110010011110000111111101001101100011011000101010001000101010110111110111100011101010101110011000100111101100011110101100110010001110001101110011111000001001110100001100110110010010100100001110001100110001011000011111101111111000000110000111010001110101010011101011001010011000000010011110001111101011010110111101110011010110101111100011101111001111000111011101011000110011101101000111001111001111011100000111010001101010111111011011100000010110100100101010001010111001000011011101100010110001100000101011000111011101011111010000101000111010111001100001011010001001001100010110111010001100011000111000000001001110001011010101001010101011100100001010010011000110110010100001011101110111011100101011111110000011111100111001110110001101000110110010110110100101001000101111100100111001100001100000000111101000111100111100111100111110010101100001001010111000000001011011000011100000110110000011010011010001010010111100110010000001000111110110001111001101001011010101100100011011101011100001100101100001001100001100010111111001001011100001101111100001000001100001111000010110001001100101001001001010100111011000100101000000000010001001101011000100100110111001110011001100100011010100010100100101110101010110010001111000111100010000011000000010101011011001010000111100110111110110100000111111001000101001101100101011111100011010100000101110100001001101101110010101001110011000010111111011000000011010000010000001001101100111011011110000011100110100000111001001000011001000101101110010100101100011100110000111111101011011011001111101110000101000001001001001100011100000011010100101110111110011111011011111101010000000110111000111101111000101100000011010101110100111100110011010011000101001000100100100000110000110010101010101001111001011010011100101111111001000101000101110100010110110001111110110110101110010100010011001111011110110010110110110100010100111101101100010100110010001111101000000100001011100001110010001100011111001100100000101111000101000011011001100101100001010001111101010100010010100000010100000011011001111001010011000011001010011001011001101101011011010001101100110011111100000001001001011101101010110100100010101000100111010010101001010101110101010001101111011111011000001001000011111110000010010011111110110110011100001001000111000000100101100001011010110110110010100011111110111111000110010010100101000101010100100010100101010110111000001001111101011100000010100010101110001001000101000101000010011000001000110111000001010111110111110010100011110111100001011100101101110101100000100100110001001010010000111110110111011011010010010100110101110000110000101101000110010100110111001100100101110000110000011011000010011101110101101101100110101101100100111101010000010011001010011111001101100110011101111000111101000111111010101010110011010110111011011110011010001100101011000001000001000110110111101011010010011001011110010010011000010010101011111010101101010000101011111000010111010001001100000010010100010010101010011001111110110010111001101111000001111111101010100011000100100100111000110010110101010111000011001110001100101010010111011011110000100111110000010011100000101011100111111110111101101111011011110010010101100110010011111011101011101101001101010111100100010101010100101010010111111000111001100000001001110000011100101000010101000111111001100000001110111001000010110101101101111011101000100011101101000110100100111110100101011010001001100000100011100000101101011111101110101101111010110100000111111010111111010001101001110011110110010000011011010000000011011100011101111010111010000101011010000011001111100001001110010011111011011001111000111000001011010010000111101001100101000000100010001001011100110010100000101100000000110000111110010101111000100101011001001011101001010111011101010001011111100010100100100101110010110110001011000110110000011111110000000000010110001011011000001011110101001000110100001000111000101101010101110001100011101101100001000111000001100001000111010011111001001000011010010111100001110110000011101001111111011100100101000001111011000100010010101000011110111101101100110011101110100000011010111000011010001100000011111000101011010100110100000011110001101101000010110001011111010000010111000110110110010100111100100001001011001100011110011000010000011111100101100110011101111011111101110100001001110111101000100101011011111001110101010101011010111100110010111111000001010101011000011000111100010011011000111101100101000000010011110111101101101000100111111011111000110011101100101101101110001000010001001010100110011111011011011010001100100010100101000000010101101011101000000001111001011101110000001000001000100010000001001001010000110100100100000111011011000011111001111100110011011101101111101000111101001110100001100110000111001000001000001000110001011010111010100000010101100001010001110000101000001101000100100111100001001110100000100101111011100000000001011111110111010101010011010101011101010010110010000111101110100000110111111100111111111011110100000010001100100011011010110111000010011110011111000011010101000110001000011000101101110000111011100001000010011100100110101000010000000110111110001100110011101001001101111001000000010001011110110110000010001110011111011000011000011010000000111000111000101110000010100110101110101100000100001100010010100001010001100100100000110101000101001010000111001100000101010100111100100011000100011000011001000011001000000101001001110101100001001000100000000111111011111010010010011000101001001000101111111110011111101011001010001010011100011100010010001011011011001001011001011101111011011110100011000010110001001100110100001011111010100011010001000010111001010010111000101110001101100101111110011000000001000000101100010100011001111101110011010101001101110111101100100110001011101011011101110011100011110101001010010111110001111011011000111100101110101110000111001000010110010001011010000101001010110111111001110000011010101010100111101110101010110101110001101011000111011110100010111000111001011110000111010001100111011100010101010101000001010010100101011111110101000010101010011110010110000000110010100111011001101011000010100111010110001110100101100100011001101111001100110110101100100011100010100101010010110111010000111000010100100010100001000000111011000101010001001110100110100011111011110010100101010101111001101111101011100101110100100001011010111011001001010100010100111010010101010110010001011001001101000011011111100001110111010100000111000001011011100000100001010110011101011111000111011111111111110111100001010000111010001001100010111101101100101010010010110001110111000111101100000110001011111111010001100111101110011101000010101000111110101110001111110011110101101001010101111101101111000011111101001000001110011100010000010110000010111010011000101100011011110010000100111000110100011100111111110001111110100000100110100010100011101101100110101110010111011000111000011100011111000111100010111001000101101011111000101010111100000100100001110111001000100110101010000001100110100111001111001011011001110101111000110000110011000001000011100011001010101000101100011011011100111110100111010001010010000000011010000000101010010010101101100010110110001011011010100111011011001110110110100010111111010111010111000010000001011010101000110101101110000011111100110000111000100000001110110111101011100110011011100100111101100111011001111010110000001101111010111100010001111110001011000101100110100101111100000101110101010001001110000100101011101100000111110101010100000011101000100111100011100011101110011000001010110111110100100000011001010111001011100100101111000011001100010001010011010001001101100010000010011100101001100110010110110001111000001100001100011000010110011100000010001111011110011011000000010011011110101011011001111110110011101010011110011001101011000011001101100011110101010111100000101100111000100001001011000010010010100111101101111101010000111111111101010100110011100101000100100110100101001000111111001001101010001111100110000001100101011010101001000111111010000111101110010000000001011111010011000000011001101111100100010000000001010000011101111001100100011001100011110001111111100000100101111110110010000101111101000110110101110110110101110111100011010110101001111000011111100100101111000100100101101000101100101010101010111011011101001111111100110001011110001111010001110001010001011000100110100101011001000101000000110001100101000111111000011000101111100101110000101101000100001001001111001101111101100100101110110000111101000000001111111100010000110001101110011011011110000010011101010111110010001011001001011110111000000100010010011000100001111111100001101100100001100110111001100110001111101010110100111010101100010100000111000001110101010100000000101100010100100100001111111000011111010111001100000011100101011010101111110111101111010110011011010000111000001011000011100101011110011010011110001001100110001001000000101111110000100111110110000101000001100101010111101110010000101011110010100101011110010010001001111100110111101000110011001010001000110011100001110001101110000001011010010101100111110100000000100011101011110011011111101111010100110001110001111110100111111011011001100101001011000110111110101010101111111010100010011110010011011000010001101110010010110100100111110001000111000101011100100100000010000010010101111001011000111110001000110010001001001001101101110101110000001110110111001010011000001011100001011011110011101100001111110000110111011100000011011101101110000101111001110011001000010001000110011110011011101100101111001100010001010110010101001001110101101101101111110100100000001010111001011011011010111001000001010010100111100001110000010001001111001000001000011010000111011110000011100110101100001011001100101111100100010100100011111011001001001111111110111001110110001100000101001000000111101101000100100100101110100000001010111010011010110101000111010010000101111111011111011000100011011100010000001111010101100001011000101100100100011000011010000110001010111011100111011011001111010011110101110010101011101010100110000000100100100010010110101111001110110111101100011001010000100001011011101000101001101011101001011110000000111001101011011111110000000000010111001011111011100010101001101000100001001011101111111100100001010001010010111011111000110011011100010110011110010100110000111001011010000111100100100100001010101100100110011011010100000100001100001011011010011001101001110110100101011101100110001110111100000110001110101111101100101101101100111110111101110010001011111011011010111010010000001111000001110000100101110111111110110000101000100011111001101010001111100001000010110101001010100100110000010000101100011110011000101011100110110010110111110011010100000111101110010001101110010010101001110110101110100111000001001001000101110111100001111001001001101010100011100111100001010111010111111101111011110001010110101000100011001101001000000011111001010101111100011010011000010010101010011101110111110000000010101000111111000010111111011111010011101111101100111101101001011010110111001111001011101111111101011010111010010101000010111100001101111010001010011111111110100100001110011110011001100010101000100010100010110100000110110001010100010011011111101000111010001110010111011011111010101001111010010100001101101001100111101101100001101100000110011000110100111111010000010011001010010010010100001111110011101111111000111110000010111011111000111011001101000110101111100011100100000011111001011101011010100110010011010000110010111101100011100110011100010000010101111100001011000000011010011100010000000101100110001100101100010111001001001010111101011001100101110110111111011111001100100011001000101000111000001110001101000000110010011001111010000101100001010100101101011010100110010101100101111100110110110111011000101100000000000011001100010101000010001101001111111000111100000100100110111100000010101100000011001001000111010111111101101001000011110011101110011010110011110010010101001010010011001000010000000111111000100100010101110111000010111011000110000000111011011000011110101001010111110101100100111010111101011110000001101011011010011101100010111101100101000111011110111110000101011010100111101010100100001010011101011110110111101001101000001010100010110100000111011111000000011011000101101101001011001001011000000010110001110101100000101110000010101110010011100111110101011010001010000110000001001111000000000010110110101010001001101100000001110011110111010010000010111111110010001011111100001010001110011000110011001101000101101011010010100010111000111010111011111100110100101111110110001110000111010001100101000111001001001111010001011001111101011010110010001100000111000110100110011011101101010110110100000001001101100111110110010101100001100011011111101110000011110110100111110101000000110010110001100001010100010101001100000010010111110110100000110010100110011001101000010000010001110001011101111111011111111110001010110000100101000000111011111001000100100100100100000010010010111000100110001010100101110110011011100101101101100000001000111011101100010010100110000111010001001111111000100011010100010101001011100001011000110100000101110101111000011011000000001010111011001101011100001111010100011110000110001000101000011000100101101110100110001101100101000010010110100001010001001111101011110101110101111001001000100100010101011001101111111011011111011111010100101011111111100101111101010011101011100100101010101110001010000000001110101101111011111110001011110011100000100101100011000011101101011010111110011101100010111101011101010111110101011101010101100101100111011101111101011000100000101010001011000111010111001101100110110001101000010000000111100100000111111011001110001101111001010010101000010110011111011000010101000100110000100000111011110110101101000110110011011100101100101101010110101101101000100100110001101011101001000110000100010100000100011100000110111101000111111010100001010001010111000100011000110111101111100110000110110010100101110001000010010110010101010001000110011100010101111111100011001010111100001010111111111101010001010000001011001001111000001010110111011101000011110101000101111110001100100100111111010101000000100111001010111101011101000000100101011111001111111010101010111110111111011101111100000011101101000010101100000010111111111111001111001111011110100111011101001100001101100111100011000000010010011001111111010101000000100110011010010101100101011110011111110100101100111000111101000100010010100101010010011111110100011110001110001111010011111011110111010101001010101101110010100111010010100110000011100110001001000010101011000000010101001010001110100110111110000011111111100111101111010100001011110101011100000101101011110110010000001101010010000110001110100010001001111101011010110101000000001001001100011110111101011111110101010001011000100010111100000111000111001000000011101100110101110100010100001100011001110101001100100101000001011111101001000000101101001000101110101100111010110010001000010100011111000001111011000011011001000000010101100000000001010000100110011001000100001000011110010000000101010001000100010010111111011110011110111011100010111111101111011111100011101101101000001000110111111110100001101111011000100110001110111100110001001100100000010011011110111110011100100000011100110111011101111000110100011101001111111001111101100010000101011010010011011011010100010101001001010000110010010111011101011110000010100110000000000111010111011000001101010110011001110000011011100110110101100011100110110000011110101110000011000011100001101100101100101010";
// load with `mzd_t *A = mzd_from_str(k, n, h);`
const std::string eW = "100112010102"; // 381
constexpr uint32_t newN = 221;
constexpr uint32_t addRows = 7;
#endif //SMALLSECRETLWE_DECODING_FILE_H