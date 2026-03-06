#include <iostream>
#include <fstream>
#include <string>
#include <nlohmann/json.hpp>
#include <regex>


using json = nlohmann::json;

std::string normalize(const std::string &in) {
    std::string s = in;
    // Convert `array([ ... ])` -> `[ ... ]`
    s = std::regex_replace(s, std::regex(R"(array\(\s*\[)"), "[");
    s = std::regex_replace(s, std::regex(R"(\]\s*\))"), "]");
    s = std::regex_replace(s, std::regex(R"(\bNone\b)"), "null");
    s = std::regex_replace(s, std::regex(R"(\bTrue\b)"), "true");
    s = std::regex_replace(s, std::regex(R"(\bFalse\b)"), "false");
    return s;
}

int main() {
    std::string pyString = "[7, array([-1.33605819, -2.44555599,  1.60127881]), 1, [0, 7, 1, 1], 0.46080430687203794, 1, 0, 0.0, [2, 2, 1]]";
    // then:
    auto js = normalize(pyString);
    json j = json::parse(js);
    auto a = j[0];        // typed: j[0].get<int>()
    auto vec = j[3].get<std::vector<int>>();
    if (vec.size() < 4) throw std::runtime_error("bad list1");
    int v0 = vec[0], v1 = vec[1], v2 = vec[2], v3 = vec[3];
    auto f5 = j[4].get<double>();
    auto f6 = j[5].get<int>();

    std::cout << "a: " << a << "\n";
    std::cout << "v2 " << v2 << "\n";
    std::cout << "f5: " << f5 << "\n";
    std::cout << "f6: " << f6 << "\n";
    return 0;
};