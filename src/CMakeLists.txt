add_library(ome ome/integration_engine_gsl.cpp
                ome/AqqQNSEven.cpp
                ome/AqqQNSOdd.cpp
                ome/AQqPS.cpp
                ome/AgqQ.cpp
                ome/AggQ.cpp
                ome/AqgQ.cpp
                ome/AqqQPS.cpp
                ome/AQg.cpp
                ome/polAqqQNSEven.cpp
                ome/polAqqQNSOdd.cpp
                ome/polAQqPS.cpp
                ome/polAgqQ.cpp
                ome/polAggQ.cpp
                ome/polAqgQ.cpp
                ome/polAqqQPS.cpp
                ome/polAQg.cpp)
set_target_properties(ome PROPERTIES SOVERSION 1)
target_include_directories(ome PUBLIC .)
target_link_libraries(ome PRIVATE GSL::gsl)
install(TARGETS ome)
install(FILES
               ome/functions.h
               ome/integration_engine.h
               ome/integration_engine_gsl.h
               ome/laurent_polynomial.h
               ome/mellin.h
               ome/ome.h
               ome/ome_type_aliases.h
               ome/piecewise.h
               ome/rpd_distribution.h
               ome/traits.h
               ome/AqqQNSEven.h
               ome/AqqQNSOdd.h
               ome/AQqPS.h
               ome/AgqQ.h
               ome/AggQ.h
               ome/AqgQ.h
               ome/AqqQPS.h
               ome/AQg.h
               ome/polAqqQNSEven.h
               ome/polAqqQNSOdd.h
               ome/polAQqPS.h
               ome/polAgqQ.h
               ome/polAggQ.h
               ome/polAqgQ.h
               ome/polAqqQPS.h
               ome/polAQg.h
               TYPE INCLUDE
              )
