cmake_minimum_required(VERSION 3.10)
project(SimpleLinearNDInterpolator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Windows用の設定: RuntimeLibraryの不一致を避ける
if(WIN32)
    # qhullライブラリ（リリースビルド）との互換性のため、常にリリース用RuntimeLibraryを使用
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
    # iterator debug levelを0に設定してqhullと一致させる
    add_compile_definitions(_ITERATOR_DEBUG_LEVEL=0)
    # デバッグ情報は保持しつつ、最適化なしでコンパイル
    set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi /MD")
endif()

# GoogleTestを有効化
enable_testing()

# FetchContentモジュールを含める
include(FetchContent)

# GoogleTestをダウンロードして利用可能にする
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.17.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# ソースファイルを設定
set(SOURCES
    src/SimpleLinearNDInterpolator.cpp
)

# ライブラリとして作成（実行可能ファイルではなく）
add_library(SimpleLinearNDInterpolator ${SOURCES})

# ターゲット固有のインクルードディレクトリを設定
target_include_directories(SimpleLinearNDInterpolator PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/third_party/qhull-2020.2/install/include
)

# Qhullライブラリをリンク（cmake-buildディレクトリのライブラリを使用）
target_link_libraries(SimpleLinearNDInterpolator PUBLIC
    ${CMAKE_SOURCE_DIR}/third_party/qhull-2020.2/cmake-build/Release/qhullcpp.lib
    ${CMAKE_SOURCE_DIR}/third_party/qhull-2020.2/cmake-build/Release/qhullstatic_r.lib
)

# testsサブディレクトリを追加
add_subdirectory(tests)
