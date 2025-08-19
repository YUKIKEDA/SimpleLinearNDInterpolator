/**
 * @file SimpleLinearNDInterpolator.cpp
 * @brief N次元線形補間器の実装
 * 
 * このファイルは、N次元空間内の散在する点群に対してDelaunay三角分割を
 * 使用した線形補間を行うクラスの実装を提供します。
 * 
 * 主な実装内容：
 * - Qhullライブラリを使用したDelaunay三角分割の構築
 * - 重心座標系を用いた線形補間の実行
 * - 効率的な単体検索アルゴリズム
 * - 線形方程式ソルバー
 * 
 * @author SimpleLinearNDInterpolator Project
 * @version 1.0
 */

#include "SimpleLinearNDInterpolator.h"
#include <stdexcept>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertex.h>
#include <libqhullcpp/QhullPoint.h>
#include <cmath>
#include <libqhullcpp/QhullVertexSet.h>
#include <limits>
#include <string>
#include <algorithm>
#include <Eigen/Dense>


/**
 * @brief 多次元値を持つ点群から補間器を構築（メインコンストラクタ）
 * 
 * 実装詳細：
 * 1. 入力データを内部形式でコピー保存
 * 2. 点群の次元数と点数を計算
 * 3. Delaunay三角分割を構築
 * 
 * データ構造の前提：
 * - points_[i][d]: i番目の点のd次元座標（point-major）
 * - values_[i][v]: i番目の点のv番目の値（point-major）
 * 
 * @note この実装では、点群データをpoint-major形式で保存
 * @note 三角分割の構築に失敗した場合は例外が投げられる
 * 
 * @example
 * std::vector<std::vector<double>> points = { // 2次元の点群（point-major）
 *     {0.0, 0.0},
 *     {1.0, 0.0},
 *     {2.0, 0.0}
 * };
 * std::vector<std::vector<double>> values = { // 各点における2成分の値（point-major）
 *     {1.0, 0.5},
 *     {2.0, 1.0},
 *     {1.5, 0.8}
 * };
 * 
 * SimpleLinearNDInterpolator interpolator(points, values);
 * 
 * std::vector<double> query_point = {0.5, 0.5}; // 補間点での値を取得
 * std::vector<double> interpolated_values = interpolator.interpolate(query_point);
 */
SimpleLinearNDInterpolator::SimpleLinearNDInterpolator(
    const std::vector<std::vector<double>> &points,
    const std::vector<std::vector<double>> &values,
    bool enable_degeneracy_fallback)
{
    // 入力検証（points）: point-major [num_points][num_dims]
    if (points.empty() || points[0].empty())
    {
        throw std::invalid_argument("points must be a non-empty 2D array");
    }
    if (!isRectangular(points))
    {
        throw std::invalid_argument("points must be rectangular");
    }
    const size_t num_points = points.size();
    const size_t num_dims = points[0].size();
    points_ = points; // point-major のまま保持
    n_points_ = static_cast<int>(num_points);
    n_dims_ = static_cast<int>(num_dims);

    // 入力検証（values）: point-major [num_points][value_dims]
    if (values.empty() || values[0].empty())
    {
        throw std::invalid_argument("values must be a non-empty 2D array");
    }
    if (!isRectangular(values))
    {
        throw std::invalid_argument("values must be rectangular");
    }
    if (values.size() != num_points)
    {
        throw std::invalid_argument("values.size() must equal points.size() (point-major)");
    }
    values_ = values; // point-major のまま保持

    if (enable_degeneracy_fallback)
    {
        // フォールバック有効：縮退状況を分析してフォールバック準備
        int effective_dims = 0;
        std::vector<std::vector<double>> projection_matrix;
        
        is_degenerate_ = !analyzeDegeneracy(points_, RankTolerance, effective_dims, projection_matrix);
        effective_dimensions_ = effective_dims;
        
        if (is_degenerate_)
        {
            // 縮退している場合：フォールバック補間の準備
            if (effective_dims <= 0)
            {
                // 全ての点が同一位置：最近傍補間のみ利用可能
                projection_matrix_ = std::nullopt;
                projected_points_ = std::nullopt;
                projected_interpolator_ = std::nullopt;
            }
            else
            {
                // 部分的縮退：低次元射影補間を準備
                projection_matrix_ = projection_matrix;
                setupProjectedInterpolation(points_, values);
            }
        }
        else
        {
            // 縮退していない：通常の処理
            projection_matrix_ = std::nullopt;
            projected_points_ = std::nullopt;
            projected_interpolator_ = std::nullopt;
        }
    }
    else
    {
        // フォールバック無効：従来の厳格なチェック
        if (num_points <= num_dims)
        {
            throw std::invalid_argument("insufficient points: need at least num_dims+1 points");
        }
        
        // 縮退チェック（フォールバック無し）
        int effective_dims = 0;
        std::vector<std::vector<double>> projection_matrix;
        bool is_valid = analyzeDegeneracy(points_, RankTolerance, effective_dims, projection_matrix);
        
        if (!is_valid)
        {
            throw std::invalid_argument(
                "Input point set is degenerate. Linear interpolation cannot be performed."
                "The point set may be concentrated on the same line, same plane, or same position."
                "To handle degenerate point sets, set enable_degeneracy_fallback to true in the constructor."
            );
        }
        
        // 縮退していない場合の設定
        is_degenerate_ = false;
        effective_dimensions_ = static_cast<int>(num_dims);
        projection_matrix_ = std::nullopt;
        projected_points_ = std::nullopt;
        projected_interpolator_ = std::nullopt;
    }
    
    // 自身が保持する点群が縮退していない場合にのみ、三角分割を構築する。
    // 縮退している場合、三角分割の役割は低次元空間で初期化された
    // 内部補間器(projected_interpolator_)が担当するため、ここでは実行しない。
    if (!is_degenerate_) 
    {
        buildTriangulation(points_);
    }

    // 1次元の場合はQhullを使用せずに点をソートするだけ
    if (n_dims_ == 1) {
        // 1次元の場合の初期化処理
        initializeFor1DInterpolation(points); 
        // 1次元の場合は三角分割不要、点のソート済みインデックスのみ準備
        return;
    }
}

/**
 * @brief スカラー値を持つ点群から補間器を構築（委譲コンストラクタ）
 * 
 * 実装詳細：
 * 1. スカラー値を1次元ベクトル値に変換
 * 2. メインコンストラクタに処理を委譲
 * 
 * @note convertTo2DVector()を使用してスカラー値をベクトル形式に統一
 * @note 内部的には全て多次元値として扱われる
 * 
 * @example
 * std::vector<std::vector<double>> points = { // 2次元空間の点群（point-major）
 *     {0.0, 0.0},
 *     {1.0, 1.0},
 *     {2.0, 0.0}
 * };
 * std::vector<double> values = {1.0, 2.0, 1.5};  // 各点でのスカラー値（例：温度）
 * 
 * SimpleLinearNDInterpolator interpolator(points, values);
 * 
 * std::vector<double> query_point = {0.5, 0.5}; // 補間点での値を取得
 * std::vector<double> interpolated_values = interpolator.interpolate(query_point); // interpolated_values[0] に補間されたスカラー値が格納される
 */
SimpleLinearNDInterpolator::SimpleLinearNDInterpolator(
    const std::vector<std::vector<double>> &points, 
    const std::vector<double> &values,
    bool enable_degeneracy_fallback
) : SimpleLinearNDInterpolator(points, convertTo2DVector(values), enable_degeneracy_fallback)
{
}

/**
 * @brief デストラクタ
 * 
 * 実装詳細：
 * - qhull_はunique_ptrなので自動的に解放される
 * - 他のSTLコンテナも自動的にメモリが解放される
 * 
 * @note 明示的なメモリ管理は不要
 */
SimpleLinearNDInterpolator::~SimpleLinearNDInterpolator()
{
}

/**
 * @brief 複数のクエリ点に対してベクトル値の補間を実行（メイン補間メソッド）
 * 
 * 実装アルゴリズム：
 * 1. 入力検証：全クエリ点の次元をチェック
 * 2. 結果配列を初期化（NaNで埋める）
 * 3. 各クエリ点について：
 *    a. 含む単体を検索し重心座標を計算
 *    b. 重心座標を重みとして線形結合で補間値を計算
 * 
 * データ構造の処理：
 * - results[i][k]: i番目のクエリ点のk番目の値成分
 * - values_[vertex][k]: vertex番目の点のk番目の値成分（point-major）
 * 
 * 補間公式：
 * v[k] = Σ(j) λ[j] * values[vertex[j]][k]
 * ここで、λ[j]は重心座標、vertex[j]は単体の頂点インデックス
 * 
 * @note 単体が見つからない点はNaNのまま残される
 * @note 計算量は O(M * S) ここでMはクエリ点数、Sは単体数
 * 
 * @example
 * SimpleLinearNDInterpolator interpolator; // 2次元空間での3つのクエリ点に対する補間
 * 
 * std::vector<std::vector<double>> query_points = { // ... データの設定 ...
 *     {1.5, 2.0},  // 1番目のクエリ点
 *     {3.0, 1.5},  // 2番目のクエリ点
 *     {2.5, 3.0}   // 3番目のクエリ点
 * };
 * 
 * auto results = interpolator.interpolate(query_points); // results[0] = 1番目のクエリ点での補間値ベクトル, results[1] = 2番目のクエリ点での補間値ベクトル, results[2] = 3番目のクエリ点での補間値ベクトル
 */
std::vector<std::vector<double>> SimpleLinearNDInterpolator::interpolate(
    const std::vector<std::vector<double>> &query_points,
    bool use_nearest_neighbor_fallback) const
{    
    // 縮退している場合のフォールバック処理
    if (is_degenerate_) {
        return interpolateWithDegenerateFallback(query_points, use_nearest_neighbor_fallback);
    }

    // 1次元の場合は専用の補間メソッドを使用
    if (n_dims_ == 1) {
        std::vector<std::vector<double>> results;
        results.reserve(query_points.size());
        for (const auto &qp : query_points) {
            results.push_back(interpolate1D(qp, use_nearest_neighbor_fallback));
        }
        return results;
    }
    
    // 入力チェック: 各クエリポイントの次元数が一致するか確認
    for (const auto &qp : query_points)
    {
        if (static_cast<int>(qp.size()) != n_dims_)
        {
            throw std::invalid_argument("query point dimension mismatch");
        }
    }

    // クエリポイントの数と値の次元数を取得
    const size_t num_queries = query_points.size();
    const size_t value_dims = values_[0].size();

    // 結果を格納する配列を初期化（初期値はNaN）
    std::vector<std::vector<double>> results(num_queries, std::vector<double>(value_dims, std::numeric_limits<double>::quiet_NaN()));

    // 各クエリポイントに対して補間を実行
    for (size_t i = 0; i < num_queries; ++i)
    {
        // 現在のクエリポイントが含まれる単体を探索し、重心座標を計算
        std::vector<double> bary;
        int simplex_idx = findSimplex(query_points[i], bary, BarycentricEpsilon);
        
        if (simplex_idx < 0)
        {
            // 単体が見つからない場合の処理
            if (use_nearest_neighbor_fallback)
            {
                // 最近傍補間を使用
                int nearest_idx = findNearestNeighbor(query_points[i]);
                if (nearest_idx >= 0)
                {
                    // 最近傍点の値をコピー
                    for (size_t k = 0; k < value_dims; ++k)
                    {
                        results[i][k] = values_[static_cast<size_t>(nearest_idx)][k];
                    }
                }
            }
            // そのポイントの結果はNaNのまま（デフォルト動作）
            continue;
        }

        // 見つかった単体の頂点インデックスを取得
        const auto &indices = simplices_.value()[static_cast<size_t>(simplex_idx)];
        
        // 重心座標を使用して線形補間を実行
        // 各値の次元に対して: v = sum_j (lambda_j * values[vertex_j][k])
        for (size_t k = 0; k < value_dims; ++k)
        {
            double acc = 0.0;  // 累積値を初期化
            
            // 単体の各頂点の値を重心座標で重み付けして合計
            for (size_t j = 0; j < indices.size(); ++j)
            {
                int vertex_index = indices[j];  // 頂点のインデックス
                acc += bary[j] * values_[static_cast<size_t>(vertex_index)][k];  // 重心座標 × 頂点値
            }
            results[i][k] = acc;  // 補間結果を保存
        }
    }

    return results;
}

/**
 * @brief 単一のクエリ点に対して補間を実行（便利メソッド）
 * 
 * 実装詳細：
 * 1. 単一点を複数点形式にラップ
 * 2. メイン補間メソッドを呼び出し
 * 3. 結果を単一ベクトルに変換
 * 
 * @note 実装の重複を避けるため、メイン補間メソッドに処理を委譲
 * @note パフォーマンス：単一点の場合でも配列操作のオーバーヘッドあり
 * 
 * @example
 * SimpleLinearNDInterpolator interpolator; // 2次元空間での単一クエリ点に対する補間
 * 
 * std::vector<double> query_point = {1.5, 2.0}; // ... データの設定 ...
 * 
 * auto result = interpolator.interpolate(query_point); // result = クエリ点での補間値ベクトル, 例：values_が2次元の場合、resultは2要素のベクトル
 * 
 * double first_component = result[0]; // 単一の値成分にアクセス
 * double second_component = result[1];
 * 
 * @example
 * std::vector<double> query_point_3d = {1.0, 2.0, 3.0}; // 3次元空間での単一クエリ点に対する補間
 * auto result_3d = interpolator.interpolate(query_point_3d); // result_3d = 3次元クエリ点での補間値ベクトル
 */
std::vector<double> SimpleLinearNDInterpolator::interpolate(
    const std::vector<double> &query_points,
    bool use_nearest_neighbor_fallback) const
{   
    // 1次元のquery_pointsを2次元ベクトルに変換
    std::vector<std::vector<double>> query_points_2d = {query_points};
    
    // メインの補間メソッドを呼び出し（フォールバック処理を含む）
    auto result_2d = interpolate(query_points_2d, use_nearest_neighbor_fallback);
    
    // 結果を1次元ベクトルに変換（クエリ1点分の値の次元ベクトル）
    if (!result_2d.empty()) {
        return result_2d[0];
    }

    // 補間が失敗した場合、値の次元数ぶんのNaNベクトルを返す
    const size_t value_dims = values_.empty() ? 0 : values_[0].size();
    return std::vector<double>(value_dims, std::numeric_limits<double>::quiet_NaN());
}

/**
 * @brief 入力点群に対してDelaunay三角分割を実行
 * 
 * 実装アルゴリズム：
 * 1. 点群データをQhullが期待する1次元配列に変換
 * 2. 次元数に応じてQhullオプションを設定
 * 3. Qhullを実行してDelaunay三角分割を計算
 * 4. 結果から単体リストを構築
 * 
 * データ変換の詳細：
 * - 内部形式：points_[i][d] (i番目の点のd次元座標、point-major)
 * - Qhull形式：[x1,y1,z1, x2,y2,z2, ...] (点順に次元を連続配置)
 * 
 * Qhullオプションの説明：
 * - d: Delaunay三角分割モード
 * - Qbb: 境界ボックスを計算してスケール
 * - Qc: 凸包の中心を保持
 * - Qz: 無限遠点を追加（数値安定性向上）
 * - Q12: エラー報告レベル設定
 * - Qt: 三角形/四面体/単体に分割
 * - Qx: 高次元（5次元以上）での数値安定性向上
 * 
 * @note 5次元以上の場合、Qxオプションを追加して数値安定性を向上
 * @note 例外処理：Qhullの内部エラーをランタイムエラーとして再投げ
 */
void SimpleLinearNDInterpolator::buildTriangulation(
    const std::vector<std::vector<double>> &points)
{
    // 点群を Qhull 期待形式 [x1,y1,..., x2,y2,...] にフラット化
    std::vector<double> flattened_points;
    flattened_points.reserve(static_cast<size_t>(n_points_ * n_dims_));
    for (int i = 0; i < n_points_; ++i) {
        for (int d = 0; d < n_dims_; ++d) {
            flattened_points.push_back(points[static_cast<size_t>(i)][static_cast<size_t>(d)]);
        }
    }

    // Qhullオプションの設定
    // Delaunay をシンプレクス分割にする（Qt）。高次元は Qx で数値安定化
    std::string qhull_options = "d Qbb Qc Qz Q12 Qt";
    if (n_dims_ >= 5) {
        qhull_options += " Qx";
    }

    try
    {
        // Qhullを初期化
        qhull_ = std::make_unique<orgQhull::Qhull>();
        qhull_->runQhull(
            "",
            n_dims_,
            n_points_,
            flattened_points.data(),
            qhull_options.c_str()
        );

        // シンプレクスリストを作成
        buildSimplexList();
    }
    catch(const std::exception& e)
    {
        throw std::runtime_error("Qhull error: " + std::string(e.what()));
    }
}

/**
 * @brief Qhullの結果から単体（simplex）のリストを構築
 * 
 * 実装詳細：
 * 1. Qhullが生成したファセット（面）を全て走査
 * 2. Delaunay三角分割の場合、上側ファセットを除外
 * 3. 非単体ファセットをスキップ（Qtオプションで単体化済み）
 * 4. 各ファセットの頂点IDリストを抽出
 * 
 * Delaunay三角分割の構造：
 * - Qhullは双対構造（convex hull of lifted points）を計算
 * - 下側ファセットがDelaunay単体に対応
 * - 上側ファセットは無限遠領域（除外対象）
 * 
 * 単体の頂点順序：
 * - Qhullが返す頂点順序をそのまま使用
 * - 重心座標計算時に順序を考慮
 * 
 * @note simplices_[i]は i番目の単体の頂点インデックスリスト
 * @note 各単体はn_dims_+1個の頂点を持つ（N次元単体）
 */
void SimpleLinearNDInterpolator::buildSimplexList()
{
    // simplices_を初期化
    simplices_ = std::vector<std::vector<int>>();

    // Qhullオブジェクトが存在しない場合は早期リターン
    if (!qhull_)
    {
        return;
    }

    // Qhull が生成したファセット（シンプレクス）を走査
    // facetList()は全てのファセットを返す
    for (orgQhull::QhullFacet facet : qhull_->facetList())
    {
        // Delaunay の場合、上側（upper）ファセットは除外（下側が Delaunay シンプレクス）
        // isUpperDelaunay()は上側ファセットかどうかを判定
        if (qhull_->isDelaunay() && facet.isUpperDelaunay())
        {
            continue; // 上側ファセットはスキップ
        }

        // 念のため非単体はスキップ（Qt 指定で単体化されているはず）
        // isSimplicial()はファセットが単体（三角形、四面体など）かどうかを判定
        if (!facet.isSimplicial())
        {
            continue; // 非単体ファセットはスキップ
        }

        // 現在のファセットの頂点インデックスを格納するベクター
        // n_dims_ + 1個の頂点を持つN次元単体のため、事前に容量を確保
        std::vector<int> simplex_indices;
        simplex_indices.reserve(static_cast<size_t>(n_dims_ + 1));

        // ファセットの頂点リストを取得
        auto vertices = facet.vertices();
        
        // 各頂点のポイントIDを抽出してsimplex_indicesに追加
        for (auto it = vertices.begin(); it != vertices.end(); ++it)
        {
            orgQhull::QhullVertex v = *it;                   // 頂点オブジェクトを取得
            int point_id = static_cast<int>(v.point().id()); // 頂点のポイントIDを取得
            simplex_indices.push_back(point_id);             // 頂点インデックスリストに追加
        }

        // 頂点インデックスが正しく抽出された場合のみ、単体リストに追加
        if (!simplex_indices.empty())
        {
            // std::moveを使用してベクターの所有権を移動（コピーを避ける）
            simplices_.value().push_back(std::move(simplex_indices));
        }
    }
}

/**
 * @brief クエリ点を含む単体を特定し、重心座標を計算
 * 
 * 実装アルゴリズム：
 * 1. 全ての単体を線形探索
 * 2. 各単体について重心座標を計算
 * 3. 重心座標の条件をチェック：
 *    - 全ての係数が非負（-ε以上）
 *    - 係数の和が1に近い（1±ε以内）
 * 
 * 重心座標の条件：
 * - λ[i] ≥ -ε (全てのi): 各係数が非負（数値誤差を許容）
 * - Σλ[i] ≈ 1 (±ε): 係数の和が1（正規化条件）
 * 
 * 数値安定性の考慮：
 * - εを使用して浮動小数点誤差を許容
 * - 境界近くの点でも正しく判定
 * 
 * パフォーマンス：
 * - 最悪計算量：O(S * N²) ここでSは単体数、Nは次元数
 * - 最適化の余地：空間分割データ構造の使用
 * 
 * @param query_point 検索対象のクエリ点座標
 * @param barycentric_coordinates 出力：計算された重心座標
 * @param eps 数値計算の許容誤差
 * @return 見つかった単体のインデックス、見つからない場合は-1
 * 
 * @note 線形探索のため、大量の単体がある場合は性能が劣化
 * @note 最初に条件を満たす単体を返す（複数該当する場合の順序は未定義）
 */
int SimpleLinearNDInterpolator::findSimplex(
    const std::vector<double> &query_point, 
    std::vector<double> &barycentric_coordinates, 
    double eps) const
{
    // 出力用の重心座標ベクトルをクリア
    barycentric_coordinates.clear();

    // 1次元の場合またはsimplices_が初期化されていない場合は-1を返す
    if (!simplices_.has_value())
    {
        return -1;
    }

    // 各シンプレクスを走査して、query_point が内点となるものを探す
    // 全単体に対して重心座標を計算し、有効な単体を見つけるまで続行
    for (size_t simplex_index = 0; simplex_index < simplices_.value().size(); ++simplex_index)
    {
        // 現在の単体での重心座標を計算
        // calculateBarycentricCoordinatesは指定された単体内でのクエリ点の重心座標を返す
        std::vector<double> lambdas = calculateBarycentricCoordinates(query_point, static_cast<int>(simplex_index));

        // 次元不一致はスキップ
        // 重心座標はn_dims_ + 1個の要素を持つべき（N次元単体はN+1個の頂点）
        if (lambdas.size() != static_cast<size_t>(n_dims_ + 1))
        {
            continue; // 次元が一致しない単体はスキップ
        }

        // すべての係数が -eps 以上か確認（非負性チェック）
        // 重心座標の各係数は理論的には非負（λ[i] >= 0）であるべき
        // 数値誤差を考慮して-eps以上の許容範囲を設定
        bool non_negative = true;
        double sum_lambda = 0.0; // 重心座標の合計を計算（正規化チェック用）
        
        for (double w : lambdas)
        {
            if (w < -eps) // 許容誤差epsを超えて負の値がある場合
            {
                non_negative = false;
                break; // 1つでも条件を満たさない場合は即座に終了
            }
            sum_lambda += w; // 合計を累積
        }

        // 非負性条件を満たさない場合は次の単体をチェック
        if (!non_negative)
        {
            continue; // 非負性を満たさない単体はスキップ
        }

        // 合計が 1 に近いか確認（正規化条件チェック）
        // 重心座標の理論的性質：Σλ[i] = 1
        // 数値誤差を考慮して1.0±epsの範囲内かチェック
        if (std::abs(sum_lambda - 1.0) <= eps)
        {
            // 条件を満たす単体が見つかった場合
            // std::moveを使用してベクターの所有権を移動（コピーを避ける）
            barycentric_coordinates = std::move(lambdas);
            // 見つかった単体のインデックスを返す
            return static_cast<int>(simplex_index);
        }
    }

    // 全ての単体をチェックしても条件を満たすものが見つからない場合
    // -1は「単体が見つからない」ことを示す特殊な値
    return -1;
}

/**
 * @brief 指定された単体内でのクエリ点の重心座標を計算
 * 
 * 実装アルゴリズム：
 * 1. 単体の最初の頂点V0を参照点として選択
 * 2. 線形方程式系 A*w = (query_point - V0) を構築
 *    - A[:, j] = V[j+1] - V0 (j番目の列は j+1番目の頂点への相対ベクトル)
 *    - b = query_point - V0 (クエリ点への相対ベクトル)
 * 3. 線形方程式を解いてw (relative coordinates) を取得
 * 4. 重心座標に変換：λ[0] = 1 - Σw[j], λ[j+1] = w[j]
 * 
 * 重心座標の定義：
 * query_point = Σ(i=0 to n) λ[i] * V[i]
 * ここで、Σλ[i] = 1, λ[i] >= 0 (全てのi)
 * 
 * 数学的背景：
 * - N次元単体はN+1個の頂点を持つ
 * - N次元の線形方程式系を解く（N個の未知数）
 * - 残り1つの係数は正規化条件から決定
 * 
 * @param query_point 重心座標を計算したいクエリ点
 * @param simplex_index 対象となる単体のインデックス
 * @return 計算された重心座標のベクトル（N+1要素）
 * 
 * @note 計算量：O(N³) （線形方程式の解法による）
 * @note 数値安定性：条件数の悪い単体では精度が劣化する可能性
 */
std::vector<double> SimpleLinearNDInterpolator::calculateBarycentricCoordinates(
    const std::vector<double> &query_point, 
    int simplex_index) const
{
    // 対象シンプレクスの頂点インデックス
    // simplices_が初期化されていない場合は空ベクトルを返す
    if (!simplices_.has_value())
    {
        return {};
    }
    
    // 単体インデックスの範囲チェック（負の値や配列サイズを超える値は無効）
    if (simplex_index < 0 || simplex_index >= static_cast<int>(simplices_.value().size()))
    {
        return {}; // 無効なインデックスの場合は空ベクトルを返す
    }
    
    // 指定された単体の頂点インデックスリストを取得
    const auto &indices = simplices_.value()[static_cast<size_t>(simplex_index)];
    
    // 頂点数の妥当性チェック（N次元単体はN+1個の頂点を持つべき）
    if (static_cast<int>(indices.size()) != n_dims_ + 1)
    {
        return {}; // 頂点数が期待値と異なる場合は空ベクトルを返す
    }

    // 参照頂点 V0 を取得（重心座標計算の基準点として使用）
    // V0は単体の最初の頂点で、他の頂点との相対位置を計算する際の原点となる
    const auto& V0 = getPointCoordinates(indices[0]);

    // 行列 A と 右辺 b を構築: A * w = (x - V0), 列 A[:,j] = V_{j+1} - V0
    // AはN×N行列、bはN次元ベクトル
    // 各列A[:,j]は頂点V_{j+1}から参照頂点V0への相対ベクトル
    std::vector<std::vector<double>> A(static_cast<size_t>(n_dims_), std::vector<double>(static_cast<size_t>(n_dims_), 0.0));
    std::vector<double> b(static_cast<size_t>(n_dims_), 0.0);

    // 右辺ベクトルbを構築: b = query_point - V0
    // クエリ点から参照頂点V0への相対ベクトルを計算
    for (int d = 0; d < n_dims_; ++d)
    {
        b[static_cast<size_t>(d)] = query_point[static_cast<size_t>(d)] - V0[static_cast<size_t>(d)];
    }

    // 係数行列Aを構築: A[:,j] = V_{j+1} - V0
    // 各列jは頂点V_{j+1}から参照頂点V0への相対ベクトル
    // これにより線形方程式系 A*w = b が構築される
    for (int j = 0; j < n_dims_; ++j)
    {
        // 頂点V_{j+1}の座標を取得（indices[0]はV0なので、j+1番目の頂点）
        const auto& Vj = getPointCoordinates(indices[static_cast<size_t>(j + 1)]);
        
        // 列jの各次元dについて相対ベクトルを計算
        for (int d = 0; d < n_dims_; ++d)
        {
            A[static_cast<size_t>(d)][static_cast<size_t>(j)] = Vj[static_cast<size_t>(d)] - V0[static_cast<size_t>(d)];
        }
    }

    // 線形方程式 A*w = b を解いてw（relative coordinates）を取得
    // wは重心座標の相対的な値で、後で正規化して重心座標λに変換する
    std::vector<double> w = solveLinearEquation(A, b);
    
    // 解の次元チェック（wはn_dims_個の要素を持つべき）
    if (static_cast<int>(w.size()) != n_dims_)
    {
        return {}; // 解が正しく得られなかった場合は空ベクトルを返す
    }

    // 相対座標wから重心座標λに変換
    // 重心座標の定義: λ[0] = 1 - Σw[j], λ[j+1] = w[j]
    // これにより Σλ[i] = 1 の正規化条件が満たされる
    std::vector<double> lambdas(static_cast<size_t>(n_dims_ + 1), 0.0);
    double sum_w = 0.0; // wの合計を計算（λ[0]の計算用）
    
    // λ[j+1] = w[j] を設定（j=0,1,...,n_dims_-1）
    for (int j = 0; j < n_dims_; ++j)
    {
        lambdas[static_cast<size_t>(j + 1)] = w[static_cast<size_t>(j)]; // λ[j+1] = w[j]
        sum_w += w[static_cast<size_t>(j)]; // 合計を累積
    }
    
    // λ[0] = 1 - Σw[j] を設定（正規化条件を満たすため）
    lambdas[0] = 1.0 - sum_w;
    
    // 計算された重心座標ベクトルを返す
    return lambdas;
}

/**
 * @brief 線形方程式 Ax = b をEigenライブラリで高速に解く
 * 
 * 実装アルゴリズム：
 * 1. std::vectorからEigen::MatrixXdとEigen::VectorXdに変換
 * 2. Eigenの高度に最適化されたPartialPivLU分解を使用
 * 3. 数値安定性と性能の両立を図る
 * 4. 結果をstd::vectorに変換して返す
 * 
 * Eigen使用の利点：
 * - SIMD最適化により高速演算（SSE、AVX等）
 * - キャッシュ効率的なメモリアクセスパターン
 * - 高度な数値安定性アルゴリズム
 * - コンパイル時最適化（小サイズ行列の場合）
 * 
 * PartialPivLU分解の特徴：
 * - 部分ピボット選択付きLU分解
 * - 一般的な正方行列に対して高い数値安定性
 * - O(N³)計算量だが高度に最適化された実装
 * - 条件数チェックによる特異行列の検出
 * 
 * @param A 係数行列（N×N）
 * @param b 右辺ベクトル（N要素）
 * @return 解ベクトル x（N要素）、解が存在しない場合は空ベクトル
 * 
 * @note Eigenの高精度算術により従来実装より数値安定
 * @note 小サイズ行列（N≤4）では固定サイズ最適化の恩恵大
 */
std::vector<double> SimpleLinearNDInterpolator::solveLinearEquation(
    const std::vector<std::vector<double>> &A, 
    const std::vector<double> &b) const
{
    // 入力行列のサイズを取得（N×N行列のN）
    const int n = static_cast<int>(A.size());
    
    // 入力の妥当性チェック
    // 1. 行列Aが空でないこと
    // 2. 行列Aが正方行列であること（行数 = 列数）
    // 3. 右辺ベクトルbのサイズが行列の行数と一致すること
    if (n == 0 || static_cast<int>(A[0].size()) != n || static_cast<int>(b.size()) != n)
    {
        return {}; // 入力が不正な場合は空ベクトルを返す
    }

    // std::vectorからEigen行列・ベクトルに変換
    // Eigen::MatrixXdは動的サイズ行列（実行時にサイズ決定）
    // 小サイズの場合はEigenが自動的に固定サイズ最適化を適用
    Eigen::MatrixXd eigenA(n, n);
    Eigen::VectorXd eigenB(n);
    
    // 係数行列Aをコピー（行メジャー順で効率的アクセス）
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            eigenA(i, j) = A[static_cast<size_t>(i)][static_cast<size_t>(j)];
        }
    }
    
    // 右辺ベクトルbをコピー
    for (int i = 0; i < n; ++i)
    {
        eigenB(i) = b[static_cast<size_t>(i)];
    }

    // PartialPivLU分解を使用して線形方程式を解く
    // PartialPivLUは部分ピボット選択付きLU分解で、一般的な行列に対して
    // 高い数値安定性と良好な性能のバランスを提供
    Eigen::PartialPivLU<Eigen::MatrixXd> solver(eigenA);
    
    // 行列の可逆性をチェック（特異行列の検出）
    // 行列式が0に近い場合は数値的に特異とみなす
    double det = eigenA.determinant();
    if (std::abs(det) < SingularMatrixEpsilon)
    {
        return {}; // 特異行列または数値的に不安定な場合は空ベクトルを返す
    }
    
    // 線形方程式を解く（Ax = bのxを求める）
    // Eigenの高度に最適化されたソルバーを使用
    Eigen::VectorXd eigenX = solver.solve(eigenB);
    
    // Eigen::VectorXdからstd::vectorに変換
    std::vector<double> x(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
    {
        x[static_cast<size_t>(i)] = eigenX(i);
    }
    
    // 計算された解ベクトルを返す
    return x;
}

/**
 * @brief 指定されたインデックスの点の座標を取得
 * 
 * 実装詳細：
 * - 内部データ構造：points_[i][d] = i番目の点のd次元座標（point-major）
 * - 出力形式：[x, y, z, ...]
 * 
 * データアクセスパターン：
 * - point-major形式のまま直接アクセス
 * - メモリ局所性は良好（連続したメモリ領域）
 * 
 * @param point_index 取得したい点のインデックス（0からn_points_-1）
 * @return 指定された点の座標ベクトルの参照（n_dims_要素）
 * 
 * @note 範囲チェックは行わない（呼び出し側で保証）
 * @note パフォーマンス：O(1) （配列の直接アクセス）
 */
const std::vector<double>& SimpleLinearNDInterpolator::getPointCoordinates(int point_index) const
{
    return points_[static_cast<size_t>(point_index)];
}

/**
 * @brief 1次元ベクトル値を2次元ベクトル形式に変換
 * 
 * 実装詳細：
 * - 入力：[v1, v2, v3, ...] （スカラー値の配列）
 * - 出力：[[v1], [v2], [v3], ...] （1要素ベクトルの配列）
 * 
 * 変換の目的：
 * - 内部処理でスカラー値とベクトル値を統一的に扱う
 * - コンストラクタの委譲コンストラクタで使用
 * - 実装の重複を避ける
 * 
 * メモリ効率：
 * - 各要素が1要素のベクトルになるため、メモリオーバーヘッドあり
 * - しかし、コードの簡潔性と保守性を優先
 * 
 * @param values 変換元の1次元ベクトル（スカラー値の配列）
 * @return 変換後の2次元ベクトル（各要素が1要素のベクトル）
 * 
 * @note パフォーマンス：O(N) （Nは値の数）
 * @note メモリ使用量：約2倍（ベクトル化によるオーバーヘッド）
 */
std::vector<std::vector<double>> SimpleLinearNDInterpolator::convertTo2DVector(
    const std::vector<double> &values
) const
{
    // スカラー値ベクトル（長さ N）を point-major ([N][1]) に変換
    std::vector<std::vector<double>> values_2d(values.size(), std::vector<double>(1, 0.0));
    for (size_t i = 0; i < values.size(); ++i)
    {
        values_2d[i][0] = values[i];
    }
    return values_2d;
}

/**
 * @brief 1次元補間用の初期化処理
 * 
 * 1次元補間用のソート済みインデックスを初期化します。
 * 
 * @param points 補間に使用する点群
 */
void SimpleLinearNDInterpolator::initializeFor1DInterpolation(const std::vector<std::vector<double>> &points)
{
    const size_t num_points = points.size();

    // n_points_が0の場合は何もしない
    if (num_points == 0) {
        sorted_indices_1d_ = std::nullopt; // 値がないことを明示
        return;
    }

    // ソート済みインデックスを計算
    std::vector<int> sorted_indices(num_points);
    for (size_t i = 0; i < num_points; ++i) {
        sorted_indices[i] = static_cast<int>(i);
    }
    
    // ラムダ式で引数の points をキャプチャする
    std::sort(sorted_indices.begin(), sorted_indices.end(), 
        [&points](int a, int b) {
            return points[static_cast<size_t>(a)][0] < points[static_cast<size_t>(b)][0];
        });
    
    // 計算結果をメンバ変数に格納
    sorted_indices_1d_ = std::move(sorted_indices);
}

/**
 * @brief 最近傍点のインデックスを見つける
 * 
 * 実装アルゴリズム：
 * 1. 全ての入力点を線形探索
 * 2. 各点とクエリ点の距離の2乗を計算
 * 3. 最小距離の点のインデックスを返す
 * 
 * 距離計算：
 * - ユークリッド距離の2乗を使用（平方根計算を省略）
 * - d² = Σ(xi - qi)² (iは次元インデックス)
 * 
 * パフォーマンス：
 * - 計算量：O(N * D) ここでNは点数、Dは次元数
 * - 最適化の余地：KD-treeなどの空間分割データ構造
 * 
 * @param query_point 検索対象のクエリ点
 * @return 最近傍点のインデックス、エラーの場合は-1
 * 
 * @note 線形探索のため、大量の点がある場合は性能が劣化
 * @note 同じ最小距離の点が複数ある場合、最初に見つかった点を返す
 */
int SimpleLinearNDInterpolator::findNearestNeighbor(const std::vector<double> &query_point) const
{
    // 入力検証：クエリ点の次元数チェック
    if (static_cast<int>(query_point.size()) != n_dims_)
    {
        return -1; // 次元不一致の場合はエラー
    }
    
    // 点群が空の場合はエラー
    if (n_points_ <= 0)
    {
        return -1;
    }
    
    // 最小距離（2乗）と対応する点のインデックスを初期化
    double min_distance_squared = std::numeric_limits<double>::max();
    int nearest_index = -1;
    
    // 全ての入力点に対して距離を計算
    for (int i = 0; i < n_points_; ++i)
    {
        // 現在の点の座標を取得
        const std::vector<double> &point = points_[static_cast<size_t>(i)];
        
        // ユークリッド距離の2乗を計算
        double distance_squared = calculateDistanceSquared(query_point, point);
        
        // より近い点が見つかった場合は更新
        if (distance_squared < min_distance_squared)
        {
            min_distance_squared = distance_squared;
            nearest_index = i;
        }
    }
    
    // 最近傍点のインデックスを返す
    return nearest_index;
}

/**
 * @brief 2点間のユークリッド距離の2乗を計算
 * 
 * 実装詳細：
 * - 平方根計算を省略してパフォーマンスを向上
 * - 距離の比較には2乗値で十分
 * - d² = Σ(xi - yi)² (iは次元インデックス)
 * 
 * パフォーマンス：
 * - 計算量：O(D) ここでDは次元数
 * - 平方根計算を省略することで高速化
 * 
 * @param point1 第1の点の座標
 * @param point2 第2の点の座標
 * @return 2点間のユークリッド距離の2乗
 * 
 * @note 次元数の一致は呼び出し側で保証される前提
 * @note 距離の大小比較には2乗値で十分
 */
double SimpleLinearNDInterpolator::calculateDistanceSquared(
    const std::vector<double> &point1, 
    const std::vector<double> &point2) const
{
    // 2点の次元数が異なる場合は無限大を返す（エラーケース）
    if (point1.size() != point2.size())
    {
        return std::numeric_limits<double>::max();
    }
    
    // ユークリッド距離の2乗を計算
    double distance_squared = 0.0;
    const size_t dims = point1.size();
    
    for (size_t d = 0; d < dims; ++d)
    {
        double diff = point1[d] - point2[d]; // 各次元での差分
        distance_squared += diff * diff;      // 差分の2乗を累積
    }
    
    return distance_squared;
}

/**
 * @brief 1次元線形補間を実行
 * 
 * 実装アルゴリズム：
 * 1. 点群をx座標でソート
 * 2. クエリ点の位置を特定
 * 3. 隣接する2点間で線形補間
 * 4. 範囲外の場合は最近傍補間または外挿
 * 
 * 線形補間公式：
 * v = v1 + (v2 - v1) * (x - x1) / (x2 - x1)
 * ここで、x1 ≤ x ≤ x2
 * 
 * @param query_point 1次元のクエリ点（1要素のベクトル）
 * @param use_nearest_neighbor_fallback 範囲外の場合に最近傍補間を使用するか
 * @return 補間された値のベクトル
 * 
 * @note 計算量：O(N log N) ソート処理による
 * @note 点群が2点未満の場合は例外を投げる
 */
std::vector<double> SimpleLinearNDInterpolator::interpolate1D(
    const std::vector<double> &query_point,
    bool use_nearest_neighbor_fallback) const
{
    // 入力検証：クエリ点の次元数チェック
    if (static_cast<int>(query_point.size()) != 1)
    {
        throw std::invalid_argument("query point must be 1-dimensional");
    }
    
    // 点群が2点未満の場合はエラー
    if (n_points_ < 2)
    {
        throw std::invalid_argument("1D interpolation requires at least 2 points");
    }

    // ソート済みインデックスが初期化されているかチェック
    if (!sorted_indices_1d_.has_value()) {
        // これは通常、点が0または1つの場合に発生する可能性がある
        if (n_points_ == 1) {
            // 点が1つしかない場合は、その点の値を返す
            return values_[0];
        }
        // それ以外は初期化エラーか、補間不可能な状態
        throw std::logic_error("1D interpolator is not properly initialized or has insufficient points.");
    }
    
    // optionalから値を取り出す (以降は .value() を付けてアクセス)
    const auto& sorted_indices = sorted_indices_1d_.value();
    
    const double query_x = query_point[0];
    const size_t value_dims = values_[0].size();    
    
    // クエリ点の位置を特定
    const double min_x = points_[sorted_indices[0]][0];
    const double max_x = points_[sorted_indices.back()][0];
    
    // 範囲外の場合の処理
    if (query_x < min_x || query_x > max_x) {
        if (use_nearest_neighbor_fallback) {
            return findNearestNeighbor1D(query_x, sorted_indices);
        } else {
            // 外挿は行わず、値の次元数ぶんのNaNを詰めたベクトルを返す
            return std::vector<double>(value_dims, std::numeric_limits<double>::quiet_NaN());
        }
    }
    
    // --- 二分探索による区間特定 ---
    // query_x 以上の最初の要素を指すイテレータを見つける
    auto it = std::lower_bound(sorted_indices.begin(), sorted_indices.end(), query_x,
        [this](int index, double value) {
            return points_[index][0] < value;
        });

    // 境界ケースの処理
    if (it == sorted_indices.begin()) {
        // query_xが最初の点と一致する場合
        return values_[*it];
    }
    if (it == sorted_indices.end()) {
        // query_xが最後の点と一致する場合 (upper_boundならあり得るが、念のため)
        return values_[sorted_indices.back()];
    }

    // 補間に使用する2点のインデックスを決定
    int idx2 = *it;
    int idx1 = *(--it); // 1つ前の要素

    double x1 = points_[idx1][0];
    double x2 = points_[idx2][0];

    // 線形補間を実行
    std::vector<double> result(value_dims);
    double denominator = x2 - x1;

    if (std::abs(denominator) < SamePositionEpsilon1D) {
        // 同じ位置の点の場合は、query_xに近い方の値を返す
        // (ここではx2側の値を採用するが、平均など他の戦略もありうる)
        return values_[idx2];
    }

    double t = (query_x - x1) / denominator;
    for (size_t k = 0; k < value_dims; ++k) {
        double v1 = values_[idx1][k];
        double v2 = values_[idx2][k];
        result[k] = v1 + t * (v2 - v1);
    }
    return result;
}

/**
 * @brief 1次元空間での最近傍補間を実行
 * 
 * 実装アルゴリズム：
 * 1. ソート済み点群から最近傍点を探索
 * 2. 距離の計算は1次元なので単純な絶対値差分
 * 3. 最近傍点の値をそのまま返す
 * 
 * 探索方法：
 * - バイナリサーチ風の効率的な探索
 * - 両端の場合は直接判定
 * - 内部の場合は隣接する点との距離比較
 * 
 * @param query_x クエリ点のx座標
 * @param sorted_indices x座標でソートされた点のインデックス
 * @return 最近傍点の値ベクトル
 * 
 * @note 計算量：O(1) 単純な距離比較のため
 * @note ソート済み配列を前提としているため効率的
 * 
 * @example
 * 1次元データでの最近傍補間の使用例
 * データ点: x=[1.0, 2.0, 3.0, 4.0], y=[10.0, 20.0, 30.0, 40.0]
 * クエリ点: x=2.7 → 最近傍点: x=3.0, y=30.0
 * 
 * std::vector<double> x_coords = {1.0, 2.0, 3.0, 4.0};
 * std::vector<std::vector<double>> values = {{10.0}, {20.0}, {30.0}, {40.0}};
 * 
 * std::vector<int> sorted_indices = {0, 1, 2, 3}; // ソート済みインデックス（x座標でソート済み）
 * 
 * double query_x = 2.7; // クエリ点での最近傍値を取得
 * std::vector<double> result = findNearestNeighbor1D(query_x, sorted_indices); // result = {30.0} (x=3.0の点の値)
 * 
 * double query_x_out = 0.5; // 範囲外のクエリ点
 * std::vector<double> result_out = findNearestNeighbor1D(query_x_out, sorted_indices); // result_out = {10.0} (x=1.0の点の値、最も近い端点)
 */
std::vector<double> SimpleLinearNDInterpolator::findNearestNeighbor1D(
    double query_x,
    const std::vector<int> &sorted_indices
) const
{
    // --- 1. エッジケースの処理 ---
    if (sorted_indices.empty()) {
        // 点群が空の場合、値の次元数に合わせてNaNベクトルを返す
        if (values_.empty()) {
            return {};
        }
        const size_t value_dims = values_[0].size();
        return std::vector<double>(value_dims, std::numeric_limits<double>::quiet_NaN());
    }

    if (sorted_indices.size() == 1) {
        return values_[static_cast<size_t>(sorted_indices[0])];
    }

    // --- 2. 範囲外のクエリの処理 ---
    // 最初の点より小さいか、最後の点より大きい場合は、最も近い端点を返す。
    const int first_idx = sorted_indices[0];
    const int last_idx = sorted_indices.back();

    if (query_x <= points_[static_cast<size_t>(first_idx)][0]) {
        return values_[static_cast<size_t>(first_idx)];
    }
    if (query_x >= points_[static_cast<size_t>(last_idx)][0]) {
        return values_[static_cast<size_t>(last_idx)];
    }

    // --- 3. 二分探索で隣接区間を特定 ---
    // query_x 以上のx座標を持つ最初の要素へのイテレータを見つける (O(log N))
    auto it = std::lower_bound(sorted_indices.begin(), sorted_indices.end(), query_x,
        [this](int index, double value) {
            // std::lower_bound のカスタム比較関数
            // points_[index][0] が value より小さい場合に true を返す
            return points_[static_cast<size_t>(index)][0] < value;
        });

    // 上記の範囲外チェックにより、itがbegin()やend()になることはないため、
    // 安全に前後の要素にアクセスできる。
    
    // it は query_x の右側にある点を指す
    const int idx2 = *it;
    // --it は query_x の左側にある点を指す
    const int idx1 = *(--it);

    // --- 4. 隣接する2点のうち、より近い方を選択 ---
    const double x1 = points_[static_cast<size_t>(idx1)][0];
    const double x2 = points_[static_cast<size_t>(idx2)][0];

    if (std::abs(query_x - x1) <= std::abs(query_x - x2)) {
        return values_[static_cast<size_t>(idx1)];
    } else {
        return values_[static_cast<size_t>(idx2)];
    }
}

/**
 * @brief 2次元配列が矩形であるかをチェック
 * 
 * 実装詳細：
 * - 空の配列は矩形ではない
 * - 各行の要素数が最初の行の要素数と一致しているかをチェック
 * 
 * @param m チェック対象の2次元配列
 * @return 矩形であれば true、そうでなければ false
 */
bool SimpleLinearNDInterpolator::isRectangular(const std::vector<std::vector<double>> &m)
{
    // 空の配列は矩形ではない（行が存在しないため）
    if (m.empty())
    {
        return false;
    }
    
    // 最初の行の要素数を基準として設定
    // 矩形配列では全ての行が同じ要素数を持つ必要がある
    const size_t cols = m[0].size();
    
    // 各行の要素数をチェック
    // 最初の行と異なる要素数を持つ行があれば、その配列は矩形ではない
    for (const auto &row : m)
    {
        if (row.size() != cols)
        {
            return false; // 要素数が一致しない行を発見した場合、早期リターン
        }
    }
    
    // 全ての行の要素数が一致している場合、配列は矩形
    return true;
}

/**
 * @brief 点群の縮退状況を詳細に分析し、フォールバック情報を取得
 * 
 * 実装詳細：
 * 1. SVD分析を実行して縮退状況を判定
 * 2. 有効な特異値に対応する右特異ベクトルから射影行列を構築
 * 3. 射影行列は元空間から有効部分空間への変換行列
 * 
 * 射影行列の構築：
 * - SVD: X = U Σ V^T
 * - 射影行列P = V[:, :effective_rank]（有効な右特異ベクトル）
 * - 射影点群: X_proj = X * P
 * 
 * フォールバック戦略：
 * - effective_dims = 0: 最近傍補間のみ
 * - effective_dims = 1: 1次元補間
 * - effective_dims = 2+: 低次元Delaunay補間
 */
bool SimpleLinearNDInterpolator::analyzeDegeneracy(
    const std::vector<std::vector<double>> &points,
    double rank_tolerance,
    int &effective_dims,
    std::vector<std::vector<double>> &projection_matrix
) const
{
    // 基本的な入力検証
    if (points.empty() || points[0].empty())
    {
        effective_dims = 0;
        projection_matrix.clear();
        return false;
    }
    
    const size_t num_points = points.size();
    const size_t num_dims = points[0].size();
    
    // 重心を計算
    std::vector<double> centroid(num_dims, 0.0);
    for (size_t i = 0; i < num_points; ++i)
    {
        for (size_t d = 0; d < num_dims; ++d)
        {
            centroid[d] += points[i][d];
        }
    }
    
    for (size_t d = 0; d < num_dims; ++d)
    {
        centroid[d] /= static_cast<double>(num_points);
    }
    
    // 中心化点群行列を構築
    Eigen::MatrixXd centered_matrix(static_cast<int>(num_points), static_cast<int>(num_dims));
    
    for (size_t i = 0; i < num_points; ++i)
    {
        for (size_t d = 0; d < num_dims; ++d)
        {
            centered_matrix(static_cast<int>(i), static_cast<int>(d)) = 
                points[i][d] - centroid[d];
        }
    }
    
    // SVD分解を実行
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        centered_matrix, 
        Eigen::ComputeThinU | Eigen::ComputeThinV
    );
    
    // 特異値と右特異ベクトルを取得
    Eigen::VectorXd singular_values = svd.singularValues();
    Eigen::MatrixXd V = svd.matrixV(); // 右特異ベクトル行列
    
    // 有効ランクを判定
    const double max_singular_value = singular_values(0);
    if (max_singular_value <= rank_tolerance)
    {
        // すべての点が重心に集中
        effective_dims = 0;
        projection_matrix.clear();
        return false;
    }
    
    // 有効な特異値の数を数える
    int valid_singular_values = 0;
    for (int i = 0; i < singular_values.size(); ++i)
    {
        double relative_singular_value = singular_values(i) / max_singular_value;
        if (relative_singular_value > rank_tolerance)
        {
            valid_singular_values++;
        }
    }
    
    effective_dims = valid_singular_values;
    
    // 射影行列を構築（有効な右特異ベクトルを使用）
    projection_matrix.assign(num_dims, std::vector<double>(valid_singular_values, 0.0));
    
    for (size_t i = 0; i < num_dims; ++i)
    {
        for (int j = 0; j < valid_singular_values; ++j)
        {
            projection_matrix[i][j] = V(static_cast<int>(i), j);
        }
    }
    
    // 縮退判定
    bool is_valid = (effective_dims >= static_cast<int>(num_dims));
    return is_valid;
}

/**
 * @brief 縮退した点群に対する射影補間器を設定
 * 
 * 実装アルゴリズム：
 * 1. 射影行列を使用して点群を低次元空間に射影
 * 2. 射影された点群で新しい補間器を構築
 * 3. 射影空間での補間結果を元空間に逆変換
 * 
 * 射影計算：
 * - 元点群: X [n_points × n_dims]
 * - 射影行列: P [n_dims × effective_dims]  
 * - 射影点群: X_proj = X * P [n_points × effective_dims]
 * 
 * 注意事項：
 * - 循環参照を避けるため、射影補間器は別インスタンス
 * - effective_dimensions_ <= 0の場合は最近傍補間のみ
 */
void SimpleLinearNDInterpolator::setupProjectedInterpolation(
    const std::vector<std::vector<double>> &points,
    const std::vector<std::vector<double>> &values
)
{
    // 射影行列が設定されていない場合は処理を中断
    if (!projection_matrix_.has_value())
    {
        return;
    }
    
    const auto &proj_matrix = projection_matrix_.value();
    const size_t num_points = points.size();
    const size_t original_dims = points[0].size();
    const size_t effective_dims = static_cast<size_t>(effective_dimensions_);
    
    // 有効次元が0の場合は射影補間不可
    if (effective_dims <= 0)
    {
        projected_points_ = std::nullopt;
        projected_interpolator_ = std::nullopt;
        return;
    }
    
    // 点群を射影空間に変換
    std::vector<std::vector<double>> projected_points(
        num_points, 
        std::vector<double>(effective_dims, 0.0)
    );
    
    for (size_t i = 0; i < num_points; ++i)
    {
        for (size_t j = 0; j < effective_dims; ++j)
        {
            double projected_coord = 0.0;
            for (size_t k = 0; k < original_dims; ++k)
            {
                projected_coord += points[i][k] * proj_matrix[k][j];
            }
            projected_points[i][j] = projected_coord;
        }
    }
    
    // 射影された点群を保存
    projected_points_ = projected_points;
    
    try 
    {
        // 射影空間での補間器を構築
        // 注意: 循環参照を避けるため、直接new演算子を使用
        projected_interpolator_ = std::make_unique<SimpleLinearNDInterpolator>(
            projected_points, 
            values,
            false // 再帰的なフォールバックを防ぐ
        );
    }
    catch (const std::exception &)
    {
        // 射影空間でも補間器の構築に失敗した場合
        // 最近傍補間にフォールバック
        projected_interpolator_ = std::nullopt;
    }
}

/**
 * @brief 縮退した点群に対するフォールバック補間を実行
 * 
 * 実装戦略：
 * 1. 有効次元数に応じて補間方法を選択
 * 2. 射影補間が利用可能な場合は射影空間で補間
 * 3. それ以外は最近傍補間にフォールバック
 * 
 * 射影補間のプロセス：
 * 1. クエリ点を射影空間に変換
 * 2. 射影空間で補間を実行  
 * 3. 結果はそのまま返す（値は変換不要）
 * 
 * フォールバック階層：
 * 1. 射影補間（effective_dimensions > 0 かつ射影補間器あり）
 * 2. 最近傍補間（その他の場合）
 */
std::vector<std::vector<double>> SimpleLinearNDInterpolator::interpolateWithDegenerateFallback(
    const std::vector<std::vector<double>> &query_points,
    bool use_nearest_neighbor_fallback
) const
{
    const size_t num_queries = query_points.size();
    const size_t value_dims = values_[0].size();
    
    // 結果を初期化（NaNで初期化）
    std::vector<std::vector<double>> results(
        num_queries, 
        std::vector<double>(value_dims, std::numeric_limits<double>::quiet_NaN())
    );
    
    // 有効次元が0の場合は最近傍補間のみ
    if (effective_dimensions_ <= 0) {
        // 全点が実質的に同一位置にあるため、全点の値の平均値を計算する
        if (n_points_ == 0) {
            return results; // 点がない場合はNaNのまま
        }

        const size_t value_dims = values_[0].size();
        std::vector<double> avg_values(value_dims, 0.0);
        for (int i = 0; i < n_points_; ++i) {
            for (size_t k = 0; k < value_dims; ++k) {
                avg_values[k] += values_[static_cast<size_t>(i)][k];
            }
        }

        for (size_t k = 0; k < value_dims; ++k) {
            avg_values[k] /= static_cast<double>(n_points_);
        }
        
        // 全てのクエリに対して同じ平均値を返す
        for (size_t i = 0; i < num_queries; ++i) {
            results[i] = avg_values;
        }
        return results;
    }
    
    // 射影補間器が利用可能な場合
    if (projected_interpolator_.has_value() && projection_matrix_.has_value()) {
        const auto &proj_matrix = projection_matrix_.value();
        const size_t original_dims = static_cast<size_t>(n_dims_);
        const size_t effective_dims = static_cast<size_t>(effective_dimensions_);
        
        // クエリ点を射影空間に変換
        std::vector<std::vector<double>> projected_queries(
            num_queries,
            std::vector<double>(effective_dims, 0.0)
        );
        
        for (size_t i = 0; i < num_queries; ++i) {
            if (query_points[i].size() != original_dims) {
                continue; // 次元不一致の場合はNaNのまま
            }
            
            for (size_t j = 0; j < effective_dims; ++j) {
                double projected_coord = 0.0;
                for (size_t k = 0; k < original_dims; ++k) {
                    projected_coord += query_points[i][k] * proj_matrix[k][j];
                }
                projected_queries[i][j] = projected_coord;
            }
        }
        
        try {
            // 射影空間で補間を実行
            auto projected_results = projected_interpolator_.value()->interpolate(
                projected_queries, 
                use_nearest_neighbor_fallback
            );
            
            // 射影補間の結果をコピー（値は変換不要）
            for (size_t i = 0; i < num_queries; ++i) {
                if (i < projected_results.size()) {
                    results[i] = projected_results[i];
                }
            }
            
            return results;
        }
        catch (const std::exception &) {
            // 射影補間に失敗した場合は最近傍補間にフォールバック
        }
    }
    
    // 最後のフォールバック：最近傍補間
    for (size_t i = 0; i < num_queries; ++i) {
        if (query_points[i].size() == static_cast<size_t>(n_dims_)) {
            int nearest_idx = findNearestNeighbor(query_points[i]);
            if (nearest_idx >= 0) {
                for (size_t k = 0; k < value_dims; ++k) {
                    results[i][k] = values_[static_cast<size_t>(nearest_idx)][k];
                }
            }
        }
    }
    
    return results;
}
