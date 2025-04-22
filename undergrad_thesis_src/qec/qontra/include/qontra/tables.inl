/*
 *  author: Suhas Vittal
 *  date:   6 January 2024
 * */

#include <limits>
#include <iostream> //to print

namespace qontra {

inline ErrorTable&
ErrorTable::operator=(const ErrorTable& other) {
    op1q = other.op1q;
    op2q = other.op2q;
    idling = other.idling;
    m1w0 = other.m1w0;
    m0w1 = other.m0w1;
    op2q_leakage_injection = other.op2q_leakage_injection;
    op2q_leakage_transport = other.op2q_leakage_transport;
    op2q_crosstalk = other.op2q_crosstalk;
    return *this;
}

inline TimeTable&
TimeTable::operator=(const TimeTable& other) {
    op1q = other.op1q;
    op2q = other.op2q;
    t1 = other.t1;
    t2 = other.t2;
    return *this;
}

namespace tables {

inline ErrorAndTiming&
ErrorAndTiming::operator*=(fp_t x) {
    e_g1q *= x;
    e_g2q *= x;
    e_m1w0 *= x;
    e_m0w1 *= x;
    e_idle *= x;
    t1 *= x <= 1e-12 ? std::numeric_limits<fp_t>::max() : 1.0/x;
    t2 *= x <= 1e-12 ? std::numeric_limits<fp_t>::max() : 1.0/x;
    return *this;
}

inline ErrorAndTiming
operator*(ErrorAndTiming et, fp_t x) {
    et *= x; return et;
}

inline ErrorAndTiming
operator*(fp_t x, ErrorAndTiming et) {
    return et * x;
}

//to print

inline void printErrorTable(const ErrorTable& et) {
    std::cout << "ErrorTable:" << std::endl;

    // Print op1q
    for (const auto& [key, inner_map] : et.op1q) {
        for (const auto& [qubit, value] : inner_map) {
            std::cout << "  op1q[" << key << "][" << qubit << "] = " << value << std::endl;
        }
    }

    // Print op2q
    for (const auto& [key, inner_map] : et.op2q) {
        for (const auto& [pair, value] : inner_map) {
            std::cout << "  op2q[" << key << "][" << pair.first << "," << pair.second << "] = " << value << std::endl;
        }
    }

    // Print other maps
    for (const auto& [qubit, value] : et.idling)
        std::cout << "  idling[" << qubit << "] = " << value << std::endl;

    for (const auto& [qubit, value] : et.m1w0)
        std::cout << "  m1w0[" << qubit << "] = " << value << std::endl;

    for (const auto& [qubit, value] : et.m0w1)
        std::cout << "  m0w1[" << qubit << "] = " << value << std::endl;
}

inline void printTimeTable(const TimeTable& tt) {
    std::cout << "TimeTable:" << std::endl;

    for (const auto& [key, inner_map] : tt.op1q) {
        for (const auto& [qubit, value] : inner_map) {
            std::cout << "  op1q[" << key << "][" << qubit << "] = " << value << std::endl;
        }
    }

    for (const auto& [key, inner_map] : tt.op2q) {
        for (const auto& [pair, value] : inner_map) {
            std::cout << "  op2q[" << key << "][" << pair.first << "," << pair.second << "] = " << value << std::endl;
        }
    }

    for (const auto& [qubit, value] : tt.t1)
        std::cout << "  t1[" << qubit << "] = " << value << std::endl;

    for (const auto& [qubit, value] : tt.t2)
        std::cout << "  t2[" << qubit << "] = " << value << std::endl;
}

inline void printErrorAndTiming(const ErrorAndTiming& et) {
    std::cout << "ErrorAndTiming:" << std::endl;
    std::cout << "  e_g1q = " << et.e_g1q << std::endl;
    std::cout << "  e_g2q = " << et.e_g2q << std::endl;
    std::cout << "  e_m1w0 = " << et.e_m1w0 << std::endl;
    std::cout << "  e_m0w1 = " << et.e_m0w1 << std::endl;
    std::cout << "  e_idle = " << et.e_idle << std::endl;
    std::cout << "  t_g1q = " << et.t_g1q << std::endl;
    std::cout << "  t_g2q = " << et.t_g2q << std::endl;
    std::cout << "  t_ro = " << et.t_ro << std::endl;
    std::cout << "  t1 = " << et.t1 << std::endl;
    std::cout << "  t2 = " << et.t2 << std::endl;
}

}   // tables

}   // qontra
