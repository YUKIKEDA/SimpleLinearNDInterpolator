/**
 * @file test_nearest_neighbor_fallback.cpp
 * @brief 最近傍フォールバック機能のテストケース
 * 
 * このファイルは、SimpleLinearNDInterpolatorの最近傍フォールバック機能を
 * テストするためのテストケースを含みます。
 * 
 * テスト内容：
 * - 基本的な最近傍補間動作
 * - 凸包外での最近傍フォールバック
 * - フォールバックオプションの有無による動作の違い
 * - エッジケース（単一点、同一距離点など）
 * 
 * @author SimpleLinearNDInterpolator Project
 * @version 1.0
 */

#include <gtest/gtest.h>
#include "SimpleLinearNDInterpolator.h"
#include <vector>
#include <cmath>
#include <limits>

/**
 * @brief 最近傍フォールバック機能のテストクラス
 */
class NearestNeighborFallbackTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // 2次元の基本的な点群を設定（正方形の頂点）
        points_2d_ = {
            {0.0, 0.0},  // 原点
            {1.0, 0.0},  // x軸上の点
            {1.0, 1.0},  // 対角点
            {0.0, 1.0}   // y軸上の点
        };

        // 対応するスカラー値
        values_scalar_ = {1.0, 2.0, 4.0, 3.0};

        // 対応するベクトル値（2次元）
        values_vector_ = {
            {1.0, 0.0},
            {2.0, 0.0},
            {4.0, 1.0},
            {3.0, 0.5}
        };

        // 3次元の点群（四面体の頂点）
        points_3d_ = {
            {0.0, 0.0, 0.0},  // 原点
            {1.0, 0.0, 0.0},  // x軸上
            {0.0, 1.0, 0.0},  // y軸上
            {0.0, 0.0, 1.0}   // z軸上
        };

        values_3d_scalar_ = {1.0, 2.0, 3.0, 4.0};
    }

    std::vector<std::vector<double>> points_2d_;
    std::vector<double> values_scalar_;
    std::vector<std::vector<double>> values_vector_;
    
    std::vector<std::vector<double>> points_3d_;
    std::vector<double> values_3d_scalar_;

    // テスト用の許容誤差
    const double eps_ = 1e-10;
};

/**
 * @brief 基本的な最近傍フォールバック機能のテスト
 * 
 * 凸包の内部点では通常の線形補間、
 * 外部点では最近傍補間が行われることを確認
 */
TEST_F(NearestNeighborFallbackTest, BasicNearestNeighborFallback)
{
    SimpleLinearNDInterpolator interp(points_2d_, values_scalar_);

    // 凸包内部の点（線形補間が使用されるべき）
    std::vector<double> inside_point = {0.5, 0.5};
    auto result_inside_no_fallback = interp.interpolate(inside_point, false);
    auto result_inside_with_fallback = interp.interpolate(inside_point, true);
    
    // 凸包内部では、フォールバックの有無に関わらず同じ結果になるべき
    EXPECT_NEAR(result_inside_no_fallback[0], result_inside_with_fallback[0], eps_);
    
    // 凸包外部の点
    std::vector<double> outside_point = {2.0, 2.0};
    auto result_outside_no_fallback = interp.interpolate(outside_point, false);
    auto result_outside_with_fallback = interp.interpolate(outside_point, true);
    
    // フォールバックなしではNaN
    EXPECT_TRUE(std::isnan(result_outside_no_fallback[0]));
    
    // フォールバックありでは最近傍点の値（この場合は点(1,1)の値4.0）
    EXPECT_FALSE(std::isnan(result_outside_with_fallback[0]));
    EXPECT_NEAR(result_outside_with_fallback[0], 4.0, eps_);
}

/**
 * @brief 複数クエリ点での最近傍フォールバックテスト
 */
TEST_F(NearestNeighborFallbackTest, MultipleQueryPointsFallback)
{
    SimpleLinearNDInterpolator interp(points_2d_, values_scalar_);

    std::vector<std::vector<double>> query_points = {
        {0.5, 0.5},   // 凸包内部
        {2.0, 2.0},   // 凸包外部（(1,1)に近い）
        {-1.0, -1.0}, // 凸包外部（(0,0)に近い）
        {1.5, 0.5}    // 凸包外部（(1,0)または(1,1)に近い）
    };

    auto results_no_fallback = interp.interpolate(query_points, false);
    auto results_with_fallback = interp.interpolate(query_points, true);

    // 凸包内部の点は同じ結果
    EXPECT_NEAR(results_no_fallback[0][0], results_with_fallback[0][0], eps_);

    // 凸包外部の点でフォールバックなしではNaN
    EXPECT_TRUE(std::isnan(results_no_fallback[1][0]));
    EXPECT_TRUE(std::isnan(results_no_fallback[2][0]));
    EXPECT_TRUE(std::isnan(results_no_fallback[3][0]));

    // フォールバックありでは最近傍点の値
    EXPECT_FALSE(std::isnan(results_with_fallback[1][0]));
    EXPECT_FALSE(std::isnan(results_with_fallback[2][0]));
    EXPECT_FALSE(std::isnan(results_with_fallback[3][0]));
    
    // 各点の最近傍点の値を確認
    EXPECT_NEAR(results_with_fallback[1][0], 4.0, eps_); // (2,2) -> (1,1)の値
    EXPECT_NEAR(results_with_fallback[2][0], 1.0, eps_); // (-1,-1) -> (0,0)の値
}

/**
 * @brief ベクトル値での最近傍フォールバックテスト
 */
TEST_F(NearestNeighborFallbackTest, VectorValuesFallback)
{
    SimpleLinearNDInterpolator interp(points_2d_, values_vector_);

    std::vector<double> outside_point = {2.0, 2.0};
    auto result_no_fallback = interp.interpolate(outside_point, false);
    auto result_with_fallback = interp.interpolate(outside_point, true);

    // フォールバックなしではNaN
    EXPECT_TRUE(std::isnan(result_no_fallback[0]));
    EXPECT_TRUE(std::isnan(result_no_fallback[1]));

    // フォールバックありでは最近傍点のベクトル値（点(1,1)の値{4.0, 1.0}）
    EXPECT_FALSE(std::isnan(result_with_fallback[0]));
    EXPECT_FALSE(std::isnan(result_with_fallback[1]));
    EXPECT_NEAR(result_with_fallback[0], 4.0, eps_);
    EXPECT_NEAR(result_with_fallback[1], 1.0, eps_);
}

/**
 * @brief 3次元空間での最近傍フォールバックテスト
 */
TEST_F(NearestNeighborFallbackTest, ThreeDimensionalFallback)
{
    SimpleLinearNDInterpolator interp(points_3d_, values_3d_scalar_);

    // 凸包外部の点
    std::vector<double> outside_point = {2.0, 2.0, 2.0};
    auto result_no_fallback = interp.interpolate(outside_point, false);
    auto result_with_fallback = interp.interpolate(outside_point, true);

    // フォールバックなしではNaN
    EXPECT_TRUE(std::isnan(result_no_fallback[0]));

    // フォールバックありでは最近傍点の値
    EXPECT_FALSE(std::isnan(result_with_fallback[0]));
    
    // (2,2,2)に最も近い点を手動で計算して確認
    // 距離: (0,0,0)=2√3, (1,0,0)=√6, (0,1,0)=√6, (0,0,1)=√6
    // 最近傍は(1,0,0), (0,1,0), (0,0,1)のいずれか（同じ距離）
    // 実装では最初に見つかった点が返される
    double expected_value = 2.0; // points_3d_[1]の値
    EXPECT_NEAR(result_with_fallback[0], expected_value, eps_);
}

/**
 * @brief エッジケース：同一距離の複数点
 */
TEST_F(NearestNeighborFallbackTest, EqualDistancePoints)
{
    // 原点を中心とした正方形の頂点
    std::vector<std::vector<double>> square_points = {
        {-1.0, -1.0},
        { 1.0, -1.0},
        { 1.0,  1.0},
        {-1.0,  1.0}
    };
    std::vector<double> square_values = {1.0, 2.0, 3.0, 4.0};

    SimpleLinearNDInterpolator interp(square_points, square_values);

    // 原点は全ての点から等距離（全て√2の距離）
    // 距離計算: (-1,-1), (1,-1), (1,1), (-1,1) 全て √2 ≈ 1.414
    std::vector<double> origin = {0.0, 0.0};
    
    // まず、原点が凸包外部にあるかを確認（フォールバックなしでNaNが返されるか）
    auto result_no_fallback = interp.interpolate(origin, false);
    
    // 原点が凸包内部の場合は線形補間が行われる
    if (!std::isnan(result_no_fallback[0]))
    {
        // 凸包内部の場合は線形補間値をテスト
        auto result_with_fallback = interp.interpolate(origin, true);
        EXPECT_NEAR(result_with_fallback[0], result_no_fallback[0], eps_);
        return; // テスト終了
    }
    
    // 凸包外部の場合は最近傍補間をテスト
    auto result = interp.interpolate(origin, true);
    EXPECT_FALSE(std::isnan(result[0]));
    
    // 等距離の場合、線形探索で最初に見つかった点が選ばれるはず
    // 実際の結果をログ出力して確認
    std::cout << "Origin interpolation result: " << result[0] << std::endl;
    
    // 4つの等距離点から1つが選ばれるため、有効な値のいずれかであることを確認
    std::vector<double> valid_values = {1.0, 2.0, 3.0, 4.0};
    bool found_valid = false;
    for (double valid_val : valid_values)
    {
        if (std::abs(result[0] - valid_val) < eps_)
        {
            found_valid = true;
            break;
        }
    }
    EXPECT_TRUE(found_valid) << "Result " << result[0] << " is not one of the expected values";
}

/**
 * @brief エッジケース：単一点での補間
 */
TEST_F(NearestNeighborFallbackTest, SinglePointInterpolation)
{
    // 最小構成の点群（2次元では3点必要だが、エラーケースとして）
    std::vector<std::vector<double>> single_point = {{0.0, 0.0}};
    std::vector<double> single_value = {42.0};

    // 単一点では三角分割が不可能なため、コンストラクタで例外が発生するべき
    EXPECT_THROW(
        SimpleLinearNDInterpolator interp(single_point, single_value),
        std::invalid_argument
    );
}

/**
 * @brief 距離計算メソッドの直接テスト（非公開メソッドのため間接的テスト）
 */
TEST_F(NearestNeighborFallbackTest, DistanceCalculationAccuracy)
{
    SimpleLinearNDInterpolator interp(points_2d_, values_scalar_);

    // 凸包外部の点で最近傍が正しく選択されることを確認
    std::vector<std::vector<double>> test_points = {
        {-0.5, -0.5}, // (0,0)に最も近い
        {1.5, -0.5},  // (1,0)に最も近い
        {1.5, 1.5},   // (1,1)に最も近い
        {-0.5, 1.5}   // (0,1)に最も近い
    };

    std::vector<double> expected_values = {1.0, 2.0, 4.0, 3.0};

    for (size_t i = 0; i < test_points.size(); ++i)
    {
        auto result = interp.interpolate(test_points[i], true);
        EXPECT_NEAR(result[0], expected_values[i], eps_)
            << "Failed for test point " << i;
    }
}

/**
 * @brief パフォーマンステスト：大量の点での最近傍検索
 */
TEST_F(NearestNeighborFallbackTest, PerformanceTest)
{
    // 100点のランダム風点群を生成
    std::vector<std::vector<double>> large_points;
    std::vector<double> large_values;
    
    for (int i = 0; i < 100; ++i)
    {
        double x = (i % 10) * 0.1;
        double y = (i / 10) * 0.1;
        large_points.push_back({x, y});
        large_values.push_back(static_cast<double>(i));
    }

    SimpleLinearNDInterpolator interp(large_points, large_values);

    // 凸包外の点で最近傍補間を実行
    std::vector<double> far_point = {10.0, 10.0};
    auto result = interp.interpolate(far_point, true);

    // 結果が有効であることを確認
    EXPECT_FALSE(std::isnan(result[0]));
    EXPECT_GE(result[0], 0.0);
    EXPECT_LT(result[0], 100.0);
}

/**
 * @brief メイン関数
 */
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}