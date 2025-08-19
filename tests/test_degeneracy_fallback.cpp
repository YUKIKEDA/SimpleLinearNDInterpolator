/**
 * @file test_degeneracy_fallback.cpp
 * @brief 縮退判定とフォールバック補間機能のテスト
 * 
 * このファイルでは、SimpleLinearNDInterpolatorの縮退検出機能と
 * フォールバック補間処理をテストします。様々な縮退パターンに対して
 * 適切なフォールバック処理が動作することを検証します。
 */

#include <gtest/gtest.h>
#include "SimpleLinearNDInterpolator.h"
#include <vector>
#include <cmath>
#include <limits>

/**
 * @brief 縮退判定とフォールバック処理テスト用のフィクスチャクラス
 */
class DegeneracyFallbackTest : public ::testing::Test {
protected:
    void SetUp() override {
        tolerance = 1e-10;
    }

    /**
     * @brief 結果がNaNでないことを確認するヘルパー関数
     */
    bool isValidResult(const std::vector<double>& result) {
        for (double val : result) {
            if (std::isnan(val)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief 2つのベクトルが近似的に等しいかをチェック
     */
    bool isApproximatelyEqual(const std::vector<double>& a, const std::vector<double>& b, double eps = 1e-10) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) {
            if (std::abs(a[i] - b[i]) > eps) {
                return false;
            }
        }
        return true;
    }

    double tolerance;
};

/**
 * @brief 2次元：正常な三角形を形成する点群のテスト
 * 
 * 非縮退の点群では通常のDelaunay補間が動作することを確認
 */
TEST_F(DegeneracyFallbackTest, ValidTriangle2D) {
    std::vector<std::vector<double>> points = {
        {0.0, 0.0},  // 原点
        {1.0, 0.0},  // x軸上の点
        {0.0, 1.0}   // y軸上の点
    };
    
    std::vector<double> values = {0.0, 1.0, 1.0};
    
    // 正常な補間器が構築できることを確認
    SimpleLinearNDInterpolator interp(points, values, true);
    
    // 三角形内部の点で補間
    std::vector<double> query_point = {0.25, 0.25};
    auto result = interp.interpolate(query_point);
    
    EXPECT_TRUE(isValidResult(result));
    EXPECT_EQ(result.size(), 1);
    EXPECT_NEAR(result[0], 0.5, tolerance);
}

/**
 * @brief 2次元：同一直線上の点群の縮退フォールバックテスト
 * 
 * 同一直線上の点群でも例外を投げずに補間結果を返すことを確認
 */
TEST_F(DegeneracyFallbackTest, CollinearPoints2D_FallbackInterpolation) {
    std::vector<std::vector<double>> points = {
        {0.0, 0.0},  // 原点
        {1.0, 1.0},  // 対角線上の点1
        {2.0, 2.0},  // 対角線上の点2
        {3.0, 3.0}   // 対角線上の点3
    };
    
    std::vector<double> values = {0.0, 1.0, 2.0, 3.0};
    
    // 縮退した点群でも補間器が構築できることを確認（フォールバック有効）
    SimpleLinearNDInterpolator interp(points, values, true);
    
    // 直線上の点で補間（射影補間が動作するはず）
    std::vector<double> query_point = {1.5, 1.5};
    auto result = interp.interpolate(query_point);
    
    EXPECT_TRUE(isValidResult(result));
    EXPECT_EQ(result.size(), 1);
    // 直線補間の結果を期待
    EXPECT_NEAR(result[0], 1.5, tolerance);
}

/**
 * @brief 3次元：同一平面上の点群の縮退フォールバックテスト
 * 
 * 3次元空間で同一平面上の点群に対して2次元射影補間が動作することを確認
 */
TEST_F(DegeneracyFallbackTest, CoplanarPoints3D_FallbackInterpolation) {
    std::vector<std::vector<double>> points = {
        {0.0, 0.0, 0.0},  // xy平面上の点
        {1.0, 0.0, 0.0},  // xy平面上の点
        {0.0, 1.0, 0.0},  // xy平面上の点
        {1.0, 1.0, 0.0}   // xy平面上の点（全てz=0）
    };
    
    std::vector<std::vector<double>> values = {
        {0.0}, {1.0}, {1.0}, {2.0}
    };
    
    // 同一平面上の点群でも補間器が構築できることを確認（フォールバック有効）
    SimpleLinearNDInterpolator interp(points, values, true);
    
    // 平面上の点で補間
    std::vector<double> query_point = {0.5, 0.5, 0.0};
    auto result = interp.interpolate(query_point);
    
    EXPECT_TRUE(isValidResult(result));
    EXPECT_EQ(result.size(), 1);
    EXPECT_NEAR(result[0], 1.0, tolerance);
}

/**
 * @brief 同一点の重複による縮退フォールバックテスト
 * 
 * 全ての点が同一位置にある場合、最近傍補間にフォールバックすることを確認
 */
TEST_F(DegeneracyFallbackTest, DuplicatePoints_NearestNeighborFallback) {
    std::vector<std::vector<double>> points = {
        {1.0, 1.0},  // 同じ位置の点
        {1.0, 1.0},  // 同じ位置の点
        {1.0, 1.0},  // 同じ位置の点
        {1.0, 1.0}   // 同じ位置の点
    };
    
    std::vector<double> values = {10.0, 20.0, 30.0, 40.0};
    
    // 同一位置の点群でも補間器が構築できることを確認（フォールバック有効）
    SimpleLinearNDInterpolator interp(points, values, true);
    
    // 任意の点で補間（最近傍補間が動作するはず）
    std::vector<double> query_point = {1.1, 1.1};
    auto result = interp.interpolate(query_point);
    
    EXPECT_TRUE(isValidResult(result));
    EXPECT_EQ(result.size(), 1);
    // 全点が同一位置の場合、値の平均が返されるはず
    double expected_average = (10.0 + 20.0 + 30.0 + 40.0) / 4.0;
    EXPECT_NEAR(result[0], expected_average, tolerance);
}

/**
 * @brief 1次元での縮退フォールバックテスト
 */
TEST_F(DegeneracyFallbackTest, DegeneratePoints1D_FallbackInterpolation) {
    std::vector<std::vector<double>> points = {
        {5.0}, {5.0}, {5.0}, {5.0}  // 全て同じx座標
    };
    
    std::vector<double> values = {1.0, 2.0, 3.0, 4.0};
    
    // 1次元縮退でも補間器が構築できることを確認
    SimpleLinearNDInterpolator interp(points, values, true);
    
    // 任意の点で補間
    std::vector<double> query_point = {6.0};
    auto result = interp.interpolate(query_point);
    
    EXPECT_TRUE(isValidResult(result));
    EXPECT_EQ(result.size(), 1);
    // 全点が同一位置の場合、値の平均が返されるはず
    double expected_average = (1.0 + 2.0 + 3.0 + 4.0) / 4.0;
    EXPECT_NEAR(result[0], expected_average, tolerance);
}

/**
 * @brief 高次元での部分縮退フォールバックテスト
 * 
 * 4次元空間で3次元部分空間に制限された点群のテスト
 */
TEST_F(DegeneracyFallbackTest, HighDimensionalDegeneracy_ProjectionFallback) {
    // 4次元空間で最後の次元が常に0（3次元部分空間に制限）
    std::vector<std::vector<double>> points = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {1.0, 1.0, 0.0, 0.0},
        {1.0, 0.0, 1.0, 0.0}   // 全ての点でw=0
    };
    
    std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
    
    // 部分縮退でも補間器が構築できることを確認
    SimpleLinearNDInterpolator interp(points, values, true);
    
    // 部分空間内の点で補間
    std::vector<double> query_point = {0.5, 0.5, 0.0, 0.0};
    auto result = interp.interpolate(query_point);
    
    EXPECT_TRUE(isValidResult(result));
    EXPECT_EQ(result.size(), 1);
}

/**
 * @brief ベクトル値での縮退フォールバックテスト
 * 
 * ベクトル値データでも縮退フォールバック処理が正常に動作することを確認
 */
TEST_F(DegeneracyFallbackTest, VectorValuedDegeneracy_Fallback) {
    std::vector<std::vector<double>> points = {
        {0.0, 0.0},  // 同一直線上の点
        {1.0, 1.0},
        {2.0, 2.0},
        {3.0, 3.0}
    };
    
    std::vector<std::vector<double>> values = {
        {1.0, 0.0, 0.0},  // 3次元ベクトル値
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 1.0, 1.0}
    };
    
    // ベクトル値の縮退データでも補間器が構築できることを確認
    SimpleLinearNDInterpolator interp(points, values, true);
    
    // 直線上の点で補間
    std::vector<double> query_point = {1.5, 1.5};
    auto result = interp.interpolate(query_point);
    
    EXPECT_TRUE(isValidResult(result));
    EXPECT_EQ(result.size(), 3);
}

/**
 * @brief 数値誤差を含む近似的な縮退のフォールバックテスト
 */
TEST_F(DegeneracyFallbackTest, NumericallyDegeneratePoinst_Fallback) {
    const double small_noise = 1e-15;  // 数値許容誤差より小さいノイズ
    
    std::vector<std::vector<double>> points = {
        {0.0, 0.0},
        {1.0, 1.0},
        {2.0, 2.0 + small_noise},  // 微小なノイズ
        {3.0, 3.0 - small_noise}   // 微小なノイズ
    };
    
    std::vector<double> values = {0.0, 1.0, 2.0, 3.0};
    
    // 数値的に縮退している点群でも補間器が構築できることを確認
    SimpleLinearNDInterpolator interp(points, values, true);
    
    // 補間を実行（射影補間または最近傍補間が動作するはず）
    std::vector<double> query_point = {1.5, 1.5};
    auto result = interp.interpolate(query_point);
    
    EXPECT_TRUE(isValidResult(result));
    EXPECT_EQ(result.size(), 1);
}

/**
 * @brief 複数クエリ点での縮退フォールバックテスト
 */
TEST_F(DegeneracyFallbackTest, MultipileQueries_DegenerateFallback) {
    std::vector<std::vector<double>> points = {
        {0.0, 0.0},  // 同一直線上
        {1.0, 1.0},
        {2.0, 2.0}
    };
    
    std::vector<double> values = {0.0, 1.0, 2.0};
    
    SimpleLinearNDInterpolator interp(points, values, true);
    
    // 複数のクエリ点で補間
    std::vector<std::vector<double>> query_points = {
        {0.5, 0.5},   // 直線上
        {1.5, 1.5},   // 直線上
        {0.25, 0.25}, // 直線上
        {2.5, 2.5}    // 直線外（外挿）
    };
    
    auto results = interp.interpolate(query_points);
    
    EXPECT_EQ(results.size(), query_points.size());
    
    for (size_t i = 0; i < results.size(); ++i) {
        EXPECT_TRUE(isValidResult(results[i]));
        EXPECT_EQ(results[i].size(), 1);
    }
}

/**
 * @brief フォールバック階層の動作確認テスト
 * 
 * 異なる縮退レベルで適切なフォールバック方法が選択されることを確認
 */
TEST_F(DegeneracyFallbackTest, FallbackHierarchy_ProperSelection) {
    // ケース1: 有効次元 = 0（全点が同一位置）→ 最近傍補間
    {
        std::vector<std::vector<double>> points = {{1.0, 1.0}, {1.0, 1.0}};
        std::vector<double> values = {10.0, 20.0};
        
        SimpleLinearNDInterpolator interp(points, values, true);
        auto result = interp.interpolate({1.1, 1.1});
        
        EXPECT_TRUE(isValidResult(result));
        // 全点が同一位置の場合、値の平均が返されるはず
        double expected_average = (10.0 + 20.0) / 2.0;
        EXPECT_NEAR(result[0], expected_average, tolerance);
    }
    
    // ケース2: 有効次元 = 1（同一直線上）→ 射影補間
    {
        std::vector<std::vector<double>> points = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}};
        std::vector<double> values = {0.0, 1.0, 2.0};
        
        SimpleLinearNDInterpolator interp(points, values, true);
        auto result = interp.interpolate({1.5, 1.5});
        
        EXPECT_TRUE(isValidResult(result));
        EXPECT_NEAR(result[0], 1.5, tolerance);  // 線形補間の結果を期待
    }
}

/**
 * @brief 境界ケース：最小限の点数での縮退テスト
 */
TEST_F(DegeneracyFallbackTest, MinimalPointSet_BoundaryCase) {
    // 2次元で2点（最小限以下）
    std::vector<std::vector<double>> points = {
        {0.0, 0.0},
        {1.0, 1.0}
    };
    
    std::vector<double> values = {1.0, 2.0};
    
    SimpleLinearNDInterpolator interp(points, values, true);
    
    auto result = interp.interpolate({0.5, 0.5});
    EXPECT_TRUE(isValidResult(result));
}