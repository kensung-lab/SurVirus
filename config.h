#ifndef SURVEYOR_CONFIG_H
#define SURVEYOR_CONFIG_H

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

int MIN_MAPQ = 0;
int MIN_CLIP_LEN = 5;
int MIN_CLIP_CONSENSUS_LEN = 15;
int MIN_DC_MAPQ = 20;
double MAX_SEQ_ERROR = 0.04;

int MAX_READ_SUPPORTED = 10000;

struct config_t {
    int threads;
    std::string rmsk_fname = "", simple_rep_fname = "";
    int read_len;
    int min_sv_len, max_sc_dist;
};

struct repeat_t {
    std::string chr;
    int start, end;
    std::string type;

    repeat_t() {}
    repeat_t(std::string& line) {
        char temp[100];
        std::stringstream ss(line);

        for (int i = 0; i < 5; i++) {
            ss >> temp;
        }
        ss >> chr >> start >> end;

        ss >> temp >> temp;
        ss >> type >> temp;
        type += "-" + std::string(temp);
    }
};

struct simple_repeat_t {
    std::string chr;
    int start, end;
    int period;

    simple_repeat_t() {}
    simple_repeat_t(std::string& line) {
        char temp[100];
        std::stringstream ss(line);

        ss >> temp;
        ss >> chr >> start >> end;
        ss >> temp;
        ss >> period;
    }
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
//    if (config_params.count("rmsk")) {
//        config.rmsk_fname = config_params["rmsk"];
//    }
//    if (config_params.count("simple_rep")) {
//        config.simple_rep_fname = config_params["simple_rep"];
//    }
//    config.read_len = stoi(config_params["read_len"]);
//    config.avg_depth = stoi(config_params["avg_depth"]);
//    config.min_is = stoi(config_params["min_is"]);
    config.min_sv_len = stoi(config_params["min_sv_len"]);
    config.max_sc_dist = stoi(config_params["max_sc_dist"]);
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
