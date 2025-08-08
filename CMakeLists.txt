cmake_minimum_required(VERSION 3.10)
project(SimpleLinearNDInterpolator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ソースファイルを設定
set(SOURCES
    src/SimpleLinearNDInterpolator.cpp
)

# 実行可能ファイルを作成
add_executable(SimpleLinearNDInterpolator ${SOURCES})

# ターゲット固有のインクルードディレクトリを設定
target_include_directories(SimpleLinearNDInterpolator PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/third_party/qhull-2020.2/install/include
)

# Qhullライブラリをリンク
target_link_libraries(SimpleLinearNDInterpolator 
    ${CMAKE_SOURCE_DIR}/third_party/qhull-2020.2/install/lib/qhullcpp.lib
    ${CMAKE_SOURCE_DIR}/third_party/qhull-2020.2/install/lib/qhull_r.lib
)
