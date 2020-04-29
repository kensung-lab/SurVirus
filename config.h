#ifndef SURVEYOR_CONFIG_H
#define SURVEYOR_CONFIG_H

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

int MIN_CLIP_LEN = 5;

int MAX_READ_SUPPORTED = 10000;

struct config_t {
    int threads;
    int min_sc_size, max_sc_dist;
    int read_len;
};

config_t parse_config(std::string file) {
    std::unordered_map<std::string, std::string> config_params;
    std::ifstream fin(file);
    std::string name, value;
    while (fin >> name >> value) {
        config_params[name] = value;
    }
    fin.close();

    config_t config;
    config.threads = stoi(config_params["threads"]);
    config.min_sc_size = stoi(config_params["min_sc_size"]);
    config.max_sc_dist = stoi(config_params["max_sc_dist"]);
    config.read_len = stoi(config_params["read_len"]);
    return config;
};


struct stats_t {
    int max_is;
};

stats_t parse_stats(std::string file) {
    std::unordered_map<std::string, std::string> config_params;
    std::ifstream fin(file);
    std::string name, value;
    while (fin >> name >> value) {
        config_params[name] = value;
    }
    fin.close();

    stats_t stats;
    stats.max_is = stoi(config_params["max_is"]);
    return stats;
}

#endif //SURVEYOR_CONFIG_H
