#include <gtest/gtest.h>
#include "SimpleLinearNDInterpolator.h"
#include <vector>
#include <cmath>

/**
 * @brief 1次元線形補間のテストケース
 * 
 * このテストファイルは、1次元データでのQhullエラーを回避するための
 * 1次元線形補間機能をテストします。
 */

class Test1DInterpolation : public ::testing::Test {
protected:
    void SetUp() override {
        // 1次元の点群データを準備
        points_1d = {
            {0.0},   // x = 0
            {1.0},   // x = 1
            {2.0},   // x = 2
            {3.0}    // x = 3
        };
        
        // 対応するスカラー値
        values_scalar = {1.0, 2.0, 4.0, 8.0};  // y = 2^x の近似
        
        // 対応するベクトル値（2次元）
        values_vector = {
            {1.0, 0.0},   // x=0での値
            {2.0, 1.0},   // x=1での値
            {4.0, 2.0},   // x=2での値
            {8.0, 3.0}    // x=3での値
        };
        
        // ソートされていない点群（補間器が自動的にソートすることをテスト）
        points_unsorted = {
            {2.0},
            {0.0},
            {3.0},
            {1.0}
        };
        values_unsorted = {4.0, 1.0, 8.0, 2.0};
    }

    std::vector<std::vector<double>> points_1d;
    std::vector<double> values_scalar;
    std::vector<std::vector<double>> values_vector;
    std::vector<std::vector<double>> points_unsorted;
    std::vector<double> values_unsorted;
};

// 基本的な1次元補間テスト（スカラー値）
TEST_F(Test1DInterpolation, BasicScalarInterpolation) {
    SimpleLinearNDInterpolator interpolator(points_1d, values_scalar);
    
    // 内挿テスト
    auto result = interpolator.interpolate({0.5});
    EXPECT_NEAR(result[0], 1.5, 1e-10);  // (1.0 + 2.0) / 2 = 1.5
    
    result = interpolator.interpolate({1.5});
    EXPECT_NEAR(result[0], 3.0, 1e-10);  // (2.0 + 4.0) / 2 = 3.0
    
    result = interpolator.interpolate({2.5});
    EXPECT_NEAR(result[0], 6.0, 1e-10);  // (4.0 + 8.0) / 2 = 6.0
}

// 1次元補間テスト（ベクトル値）
TEST_F(Test1DInterpolation, BasicVectorInterpolation) {
    SimpleLinearNDInterpolator interpolator(points_1d, values_vector);
    
    // 内挿テスト
    auto result = interpolator.interpolate({0.5});
    EXPECT_EQ(result.size(), 2);
    EXPECT_NEAR(result[0], 1.5, 1e-10);  // (1.0 + 2.0) / 2 = 1.5
    EXPECT_NEAR(result[1], 0.5, 1e-10);  // (0.0 + 1.0) / 2 = 0.5
    
    result = interpolator.interpolate({1.5});
    EXPECT_EQ(result.size(), 2);
    EXPECT_NEAR(result[0], 3.0, 1e-10);  // (2.0 + 4.0) / 2 = 3.0
    EXPECT_NEAR(result[1], 1.5, 1e-10);  // (1.0 + 2.0) / 2 = 1.5
}

// 範囲外のクエリに対してNaNが返されることをテスト
TEST_F(Test1DInterpolation, ExtrapolationReturnsNaN) {
    SimpleLinearNDInterpolator interpolator(points_1d, values_scalar);
    
    // 左側の範囲外の点
    auto result_left = interpolator.interpolate({-1.0});
    
    // 結果が空でないことを確認 (次元数が合っていれば空にはならないはず)
    ASSERT_FALSE(result_left.empty());
    // 結果の最初の要素が NaN であることを確認
    EXPECT_TRUE(std::isnan(result_left[0]));
    
    // 右側の範囲外の点
    auto result_right = interpolator.interpolate({4.0});

    ASSERT_FALSE(result_right.empty());
    EXPECT_TRUE(std::isnan(result_right[0]));
}

// 範囲外での最近傍補間テスト
TEST_F(Test1DInterpolation, NearestNeighborFallback) {
    SimpleLinearNDInterpolator interpolator(points_1d, values_scalar);
    
    // 左側最近傍
    auto result = interpolator.interpolate({-1.0}, true);
    EXPECT_NEAR(result[0], 1.0, 1e-10);  // 最近傍点 x=0 の値
    
    // 右側最近傍
    result = interpolator.interpolate({5.0}, true);
    EXPECT_NEAR(result[0], 8.0, 1e-10);  // 最近傍点 x=3 の値
}

// 複数点同時補間テスト
TEST_F(Test1DInterpolation, MultiplePointsInterpolation) {
    SimpleLinearNDInterpolator interpolator(points_1d, values_scalar);
    
    std::vector<std::vector<double>> query_points = {
        {0.5},
        {1.5},
        {2.5}
    };
    
    auto results = interpolator.interpolate(query_points);
    EXPECT_EQ(results.size(), 3);
    
    EXPECT_NEAR(results[0][0], 1.5, 1e-10);
    EXPECT_NEAR(results[1][0], 3.0, 1e-10);
    EXPECT_NEAR(results[2][0], 6.0, 1e-10);
}

// ソートされていない点群でのテスト
TEST_F(Test1DInterpolation, UnsortedPoints) {
    SimpleLinearNDInterpolator interpolator(points_unsorted, values_unsorted);
    
    // 内挿テスト（補間器が自動的にソートすることを確認）
    auto result = interpolator.interpolate({0.5});
    EXPECT_NEAR(result[0], 1.5, 1e-10);  // (1.0 + 2.0) / 2 = 1.5
    
    result = interpolator.interpolate({1.5});
    EXPECT_NEAR(result[0], 3.0, 1e-10);  // (2.0 + 4.0) / 2 = 3.0
}

// 境界値での補間テスト
TEST_F(Test1DInterpolation, BoundaryValues) {
    SimpleLinearNDInterpolator interpolator(points_1d, values_scalar);
    
    // 既存の点での補間（そのまま値が返されるべき）
    auto result = interpolator.interpolate({0.0});
    EXPECT_NEAR(result[0], 1.0, 1e-10);
    
    result = interpolator.interpolate({1.0});
    EXPECT_NEAR(result[0], 2.0, 1e-10);
    
    result = interpolator.interpolate({2.0});
    EXPECT_NEAR(result[0], 4.0, 1e-10);
    
    result = interpolator.interpolate({3.0});
    EXPECT_NEAR(result[0], 8.0, 1e-10);
}

// 同じx座標を持つ点がある場合のテスト
TEST_F(Test1DInterpolation, DuplicateXCoordinates) {
    std::vector<std::vector<double>> points_dup = {
        {0.0},
        {1.0},
        {1.0},  // 重複
        {2.0}
    };
    std::vector<double> values_dup = {1.0, 2.0, 3.0, 4.0};
    
    SimpleLinearNDInterpolator interpolator(points_dup, values_dup);
    
    // 重複点での補間（ソート後の最初の値が返される）
    auto result = interpolator.interpolate({1.0});
    // ソート後は [0.0, 1.0, 1.0, 2.0] の順番、値は [1.0, 2.0, 3.0, 4.0]
    // x=1.0での最初の値は2.0
    EXPECT_NEAR(result[0], 2.0, 1e-10);
}

// エラーケースのテスト
TEST_F(Test1DInterpolation, ErrorCases) {
    // 点数不足のテスト
    std::vector<std::vector<double>> single_point = {{0.0}};
    std::vector<double> single_value = {1.0};
    
    EXPECT_THROW(
        SimpleLinearNDInterpolator interpolator(single_point, single_value),
        std::invalid_argument
    );
    
    // 次元不一致のテスト
    SimpleLinearNDInterpolator interpolator(points_1d, values_scalar);
    EXPECT_THROW(
        interpolator.interpolate({0.0, 1.0}),  // 2次元のクエリ点
        std::invalid_argument
    );
}

// 線形関数での正確性テスト
TEST_F(Test1DInterpolation, LinearFunctionAccuracy) {
    // y = 2x + 1 の線形関数
    std::vector<std::vector<double>> points = {{0.0}, {1.0}, {2.0}};
    std::vector<double> values = {1.0, 3.0, 5.0};  // y = 2x + 1
    
    SimpleLinearNDInterpolator interpolator(points, values);
    
    // 任意の点での補間が正確であることを確認
    for (double x = 0.1; x <= 1.9; x += 0.1) {
        auto result = interpolator.interpolate({x});
        double expected = 2.0 * x + 1.0;
        EXPECT_NEAR(result[0], expected, 1e-10) << "Failed at x = " << x;
    }
}