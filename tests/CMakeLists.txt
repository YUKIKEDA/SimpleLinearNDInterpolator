# GoogleTestを使用したテストを設定

# ソースファイルを設定（テスト用に必要な実装を含める）
set(TEST_SOURCES
  ${CMAKE_SOURCE_DIR}/src/SimpleLinearNDInterpolator.cpp
)

# 自動生成テストが存在する場合は取り込む
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/generated_scipy_compat_test.cpp)
  list(APPEND TEST_SOURCES generated_scipy_compat_test.cpp)
endif()

# 最近傍フォールバックテストを追加
list(APPEND TEST_SOURCES test_nearest_neighbor_fallback.cpp)

# 1次元補間テストを追加
list(APPEND TEST_SOURCES test_1d_interpolation.cpp)

# テスト用の実行可能ファイルを作成
add_executable(
  test_interpolator
  ${TEST_SOURCES}
)

# テスト実行可能ファイルにライブラリをリンク
target_link_libraries(
  test_interpolator
  gtest_main
  ${CMAKE_SOURCE_DIR}/third_party/qhull-2020.2/cmake-build/Release/qhullcpp.lib
  ${CMAKE_SOURCE_DIR}/third_party/qhull-2020.2/cmake-build/Release/qhullstatic_r.lib
)

# テスト実行可能ファイルのインクルードディレクトリを設定
target_include_directories(test_interpolator PRIVATE
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/third_party/qhull-2020.2/install/include
)

# CTestにテストを登録
include(GoogleTest)
gtest_discover_tests(test_interpolator)
