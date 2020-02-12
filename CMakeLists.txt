cmake_minimum_required(VERSION 3.5)
project(SurVirus)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

include_directories(htslib-1.9/include)
include_directories(sparsehash-sparsehash-2.0.3/include)

set(ISOLATE_RELEVANT_PAIRS_SRC isolate_relevant_pairs.cpp sam_utils.h config.h libs/cptl_stl.h)
add_executable(isolate_relevant_pairs ${ISOLATE_RELEVANT_PAIRS_SRC})

set(FILTER_BY_QNAME_SRC filter_by_qname.cpp sam_utils.h config.h libs/cptl_stl.h)
add_executable(filter_by_qname ${FILTER_BY_QNAME_SRC})

set(EXTRACT_CLIPS_SRC extract_clips.cpp sam_utils.h config.h libs/cptl_stl.h)
add_executable(extract_clips ${EXTRACT_CLIPS_SRC})

set(READS_CATEGORIZER_SRC reads_categorizer.cpp sam_utils.h config.h libs/cptl_stl.h)
add_executable(reads_categorizer ${READS_CATEGORIZER_SRC})

set(REMAPPER_SRC libs/ssw.c libs/ssw_cpp.cpp remapper.cpp sam_utils.h config.h libs/cptl_stl.h utils.h)
add_executable(remapper ${REMAPPER_SRC})

set(CONSENSUS_BUILDER_SRC bp_region_consensus_builder.cpp)
add_executable(bp_region_consensus_builder ${CONSENSUS_BUILDER_SRC})

set(FILTER_SRC libs/ssw.c libs/ssw_cpp.cpp filter.cpp)
add_executable(filter ${FILTER_SRC})

set(GET_REALIGNMENTS_FOR_QNAME get-realignments-for-qname.cpp)
add_executable(get-realignments-for-qname ${GET_REALIGNMENTS_FOR_QNAME})

find_library(HTS_LIB hts htslib-1.9/lib)
target_link_libraries(isolate_relevant_pairs "${HTS_LIB}")
target_link_libraries(filter_by_qname "${HTS_LIB}")
target_link_libraries(extract_clips "${HTS_LIB}")
target_link_libraries(reads_categorizer "${HTS_LIB}")
target_link_libraries(remapper "${HTS_LIB}")
target_link_libraries(bp_region_consensus_builder "${HTS_LIB}")
target_link_libraries(filter "${HTS_LIB}")
target_link_libraries(get-realignments-for-qname "${HTS_LIB}")
