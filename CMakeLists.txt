cmake_minimum_required(VERSION 3.5)
project(SurVirus)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

include_directories(htslib-1.11/include)
include_directories(sparsehash-sparsehash-2.0.3/include)

set(ISOLATE_RELEVANT_PAIRS_SRC isolate_relevant_pairs.cpp)
add_executable(isolate_relevant_pairs ${ISOLATE_RELEVANT_PAIRS_SRC})

set(ISOLATE_RELEVANT_PAIRS_FQ_SRC isolate_relevant_pairs_fq.cpp)
add_executable(isolate_relevant_pairs_fq ${ISOLATE_RELEVANT_PAIRS_FQ_SRC})

set(FILTER_BY_QNAME_SRC filter_by_qname.cpp)
add_executable(filter_by_qname ${FILTER_BY_QNAME_SRC})

set(EXTRACT_CLIPS_SRC extract_clips.cpp)
add_executable(extract_clips ${EXTRACT_CLIPS_SRC})

set(READS_CATEGORIZER_SRC reads_categorizer.cpp)
add_executable(reads_categorizer ${READS_CATEGORIZER_SRC})

set(MERGE_RETAINED_READS merge_retained_reads.cpp)
add_executable(merge_retained_reads ${MERGE_RETAINED_READS})

set(BUILD_REGION_READS_ASSOCIATIONS libs/ssw.c libs/ssw_cpp.cpp build_region-reads_associations.cpp)
add_executable(build_region-reads_associations ${BUILD_REGION_READS_ASSOCIATIONS})

set(REMAPPER_SRC libs/ssw.c libs/ssw_cpp.cpp remapper.cpp)
add_executable(remapper ${REMAPPER_SRC})

set(CONSENSUS_BUILDER_SRC libs/ssw.c libs/ssw_cpp.cpp bp_region_consensus_builder.cpp)
add_executable(bp_region_consensus_builder ${CONSENSUS_BUILDER_SRC})

set(FILTER_SRC libs/ssw.c libs/ssw_cpp.cpp filter.cpp)
add_executable(filter ${FILTER_SRC})

set(GET_REALIGNMENTS_FOR_QNAME get-realignments-for-qname.cpp)
add_executable(get-realignments-for-qname ${GET_REALIGNMENTS_FOR_QNAME})

find_library(HTS_LIB hts PATHS htslib-1.11/lib NO_DEFAULT_PATH)
find_package(ZLIB)
target_link_libraries(isolate_relevant_pairs "${HTS_LIB}")
target_link_libraries(isolate_relevant_pairs_fq "${HTS_LIB}")
target_link_libraries(isolate_relevant_pairs_fq ZLIB::ZLIB)
target_link_libraries(filter_by_qname "${HTS_LIB}")
target_link_libraries(extract_clips "${HTS_LIB}")
target_link_libraries(reads_categorizer "${HTS_LIB}")
target_link_libraries(merge_retained_reads "${HTS_LIB}")
target_link_libraries(build_region-reads_associations "${HTS_LIB}")
target_link_libraries(remapper "${HTS_LIB}")
target_link_libraries(bp_region_consensus_builder "${HTS_LIB}")
target_link_libraries(filter "${HTS_LIB}")
target_link_libraries(get-realignments-for-qname "${HTS_LIB}")
