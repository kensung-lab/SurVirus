cmake_minimum_required(VERSION 3.8)
project(SurVirus)

set(CMAKE_CXX_STANDARD 11)

include_directories(/usr/local/include)

set(ISOLATE_RELEVANT_PAIRS_SRC isolate_relevant_pairs.cpp sam_utils.h config.h libs/cptl_stl.h)
add_executable(isolate_relevant_pairs ${ISOLATE_RELEVANT_PAIRS_SRC})

set(FILTER_BY_QNAME_SRC filter_by_qname.cpp sam_utils.h config.h libs/cptl_stl.h)
add_executable(filter_by_qname ${FILTER_BY_QNAME_SRC})

set(READS_CATEGORIZER_SRC reads_categorizer.cpp sam_utils.h config.h libs/cptl_stl.h)
add_executable(reads_categorizer ${READS_CATEGORIZER_SRC})

set(REMAPPER_SRC libs/ssw.c libs/ssw_cpp.cpp remapper.cpp sam_utils.h config.h libs/cptl_stl.h cluster.h)
add_executable(remapper ${REMAPPER_SRC})

find_library(HTS_LIB hts)
target_link_libraries(isolate_relevant_pairs "${HTS_LIB}")
target_link_libraries(filter_by_qname "${HTS_LIB}")
target_link_libraries(reads_categorizer "${HTS_LIB}")
target_link_libraries(remapper "${HTS_LIB}")
