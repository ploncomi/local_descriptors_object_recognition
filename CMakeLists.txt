# cmake .. -DOpenCV_DIR=c:/opencv/build

cmake_minimum_required(VERSION 3.5)

project(siftDetector VERSION 1.0 LANGUAGES CXX)

add_definitions(-D__COMPAT_CVMAT__=1)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

set(PROJECT_SOURCES
        L_RecPuntInt/L_Descriptores.cpp
        L_RecPuntInt/L_Fnes_Pato.cpp
        L_RecPuntInt/L_Fnes_Pato_optim.cpp
        L_RecPuntInt/L_GenPuntInt.cpp
        L_RecPuntInt/L_RecPuntInt.cpp
        L_RecPuntInt/siftImpl/L_CvSift.cpp
        simple_interface.cpp
        main.cpp
)


add_executable(siftDetector
    ${PROJECT_SOURCES}
)

target_link_libraries(siftDetector ${OpenCV_LIBS})

