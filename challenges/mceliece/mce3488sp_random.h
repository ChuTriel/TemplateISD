//
// Created by Sebastian on 22.08.2023.
//

#ifndef DECODING_MCE3488_RANDOM_H
#define DECODING_MCE3488_RANDOM_H

#include <cstdint>
#include <iostream>
constexpr uint64_t n = 3488;
constexpr uint64_t k = 2720;
constexpr uint64_t w = 64;

constexpr const char* s = "001010010000101100000011010010001110001001000101101011100111000111100001010111001000010110110000111100011011111111001011111000000001111101111010001000011001111001110010100010101101100111111001101110011111100010101100110010110001101001111010010001010111010001111001011010010001010110100101111101001111011110110100100111111000000001110000101000000101010101011001010011110101101001100011011010100010111010101111110001100001010011111110111011111000000001100100011011011000000101100111000001010010001111101110000010000011011011101111011100100010001101110110100000011001011100100111001010101011000010100110011011001011000110110100100011010101010000101011101100011011011011111101000011101001101101110111100111101001101101111000000000001011001110111010000010100000000101111001";

const char* eTest = "00000000000000000000000000000000000000000100000000000001000000000000010000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001010000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000100000000000000000010000000100000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000100010000000000000000000000000000000000000000010000000000000000000000000010000000000000000000010000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000100000000001000000000000000000100000000000000000000000000000000000000000000100000000000000000000000001000000000000000000000000000000000000000000000000000000001001000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000001000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000010100000000000000000000000000000000000010000000000000000000000000000000000000000000001000000000000000010000000000000000000010000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000010000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000100000000000000000000000000000000000000000000010000000000000000000000100000000000000000000000000001000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000010000";
std::string eW = "0211000003000221012112001021110201000111010010010100011000001201121000000100200020110010100010000101102110101";
constexpr uint32_t newN = 1600;
constexpr uint32_t addRows = 50;
#endif//DECODING_MCE3488_RANDOM_H