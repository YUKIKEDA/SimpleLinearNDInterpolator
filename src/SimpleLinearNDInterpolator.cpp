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
 * - 線形方程式ソルバー（ガウス消去法）
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

/**
 * @brief 多次元値を持つ点群から補間器を構築（メインコンストラクタ）
 * 
 * 実装詳細：
 * 1. 入力データを内部形式でコピー保存
 * 2. 点群の次元数と点数を計算
 * 3. Delaunay三角分割を構築
 * 
 * データ構造の前提：
 * - points_[d][i]: d次元目のi番目の点の座標
 * - values_[v][i]: v番目の値のi番目の点での値
 * 
 * @note この実装では、点群データを次元優先（dimension-major）形式で保存
 * @note 三角分割の構築に失敗した場合は例外が投げられる
 */
SimpleLinearNDInterpolator::SimpleLinearNDInterpolator(
    const std::vector<std::vector<double>> &points, 
    const std::vector<std::vector<double>> &values)
{
    points_ = points;
    values_ = values;
    n_dims_ = points.size();
    n_points_ = points[0].size();

    buildTriangulation(points_);
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
 */
SimpleLinearNDInterpolator::SimpleLinearNDInterpolator(
    const std::vector<std::vector<double>> &points, 
    const std::vector<double> &values
) : SimpleLinearNDInterpolator(points, convertTo2DVector(values))
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
 * - values_[k][vertex]: k番目の値成分のvertex番目の点での値
 * 
 * 補間公式：
 * v[k] = Σ(j) λ[j] * values[k][vertex[j]]
 * ここで、λ[j]は重心座標、vertex[j]は単体の頂点インデックス
 * 
 * @note 単体が見つからない点はNaNのまま残される
 * @note 計算量は O(M * S) ここでMはクエリ点数、Sは単体数
 */
std::vector<std::vector<double>> SimpleLinearNDInterpolator::interpolate(
    const std::vector<std::vector<double>> &query_points) const
{
    const double eps = 1e-10;

    // 入力チェック
    for (const auto &qp : query_points)
    {
        if (static_cast<int>(qp.size()) != n_dims_)
        {
            throw std::invalid_argument("query point dimension mismatch");
        }
    }

    const size_t num_queries = query_points.size();
    const size_t value_dims = values_.size();

    std::vector<std::vector<double>> results(num_queries, std::vector<double>(value_dims, std::numeric_limits<double>::quiet_NaN()));

    for (size_t i = 0; i < num_queries; ++i)
    {
        std::vector<double> bary;
        int simplex_idx = findSimplex(query_points[i], bary, eps);
        if (simplex_idx < 0)
        {
            // 単体が見つからない場合は NaN のまま
            continue;
        }

        const auto &indices = simplices_[static_cast<size_t>(simplex_idx)];
        // 値の線形結合: v = sum_j lambda_j * values(:, vertex_j)
        for (size_t k = 0; k < value_dims; ++k)
        {
            double acc = 0.0;
            for (size_t j = 0; j < indices.size(); ++j)
            {
                int vertex_index = indices[j];
                acc += bary[j] * values_[k][static_cast<size_t>(vertex_index)];
            }
            results[i][k] = acc;
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
 */
std::vector<double> SimpleLinearNDInterpolator::interpolate(
    const std::vector<double> &query_points) const
{
    // 1次元のquery_pointsを2次元ベクトルに変換
    std::vector<std::vector<double>> query_points_2d = {query_points};
    
    // 1番目のinterpolateメソッドを呼び出し
    auto result_2d = interpolate(query_points_2d);
    
    // 結果を1次元ベクトルに変換（クエリ1点分の値の次元ベクトル）
    if (!result_2d.empty()) {
        return result_2d[0];
    }
    return std::vector<double>();
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
 * - 内部形式：points_[d][i] (d次元目のi番目の点)
 * - Qhull形式：[x1,y1,z1, x2,y2,z2, ...] (点順に次元を連続配置)
 * 
 * Qhullオプションの説明：
 * - d: Delaunay三角分割モード
 * - Qbb: 境界ボックスを計算してスケール
 * - Qc: 凸包の中心を保持
 * - Qz: 無限遠点を追加（数値安定性向上）
 * - Q12: エラー報告レベル設定
 * - Qx: 高次元（5次元以上）での数値安定性向上
 * - Qt: 三角形/四面体/単体に分割
 * 
 * @note 5次元以上の場合、Qxオプションを追加して数値安定性を向上
 * @note 例外処理：Qhullの内部エラーをランタイムエラーとして再投げ
 */
void SimpleLinearNDInterpolator::buildTriangulation(
    const std::vector<std::vector<double>> &points)
{
    // 点群を1次元化
    std::vector<double> flattened_points;
    flattened_points.reserve(points.size() * points[0].size());
    for (const auto& point : points) {
        flattened_points.insert(flattened_points.end(), point.begin(), point.end());
    }

    // Qhullオプションの設定
    std::string qhull_options = "d Qbb Qc Qz Q12";
    if (n_dims_ >= 5) {
        qhull_options += " Qx";
    }
    qhull_options += " Qt";

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
    simplices_.clear();

    if (!qhull_)
    {
        return;
    }

    // Qhull が生成したファセット（シンプレクス）を走査
    for (orgQhull::QhullFacet facet : qhull_->facetList())
    {
        // Delaunay の場合、上側（upper）ファセットは除外（下側が Delaunay シンプレクス）
        if (qhull_->isDelaunay() && facet.isUpperDelaunay())
        {
            continue;
        }

        // 念のため非単体はスキップ（Qt 指定で単体化されているはず）
        if (!facet.isSimplicial())
        {
            continue;
        }

        std::vector<int> simplex_indices;
        simplex_indices.reserve(static_cast<size_t>(n_dims_ + 1));

        auto vertices = facet.vertices();
        for (auto it = vertices.begin(); it != vertices.end(); ++it)
        {
            orgQhull::QhullVertex v = *it;
            int point_id = static_cast<int>(v.point().id());
            simplex_indices.push_back(point_id);
        }

        if (!simplex_indices.empty())
        {
            simplices_.push_back(std::move(simplex_indices));
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
    barycentric_coordinates.clear();

    // 各シンプレクスを走査して、query_point が内点となるものを探す
    for (size_t simplex_index = 0; simplex_index < simplices_.size(); ++simplex_index)
    {
        // 重心座標を計算
        std::vector<double> lambdas = calculateBarycentricCoordinates(query_point, static_cast<int>(simplex_index));

        // 次元不一致はスキップ
        if (lambdas.size() != static_cast<size_t>(n_dims_ + 1))
        {
            continue;
        }

        // すべての係数が -eps 以上か確認
        bool non_negative = true;
        double sum_lambda = 0.0;
        for (double w : lambdas)
        {
            if (w < -eps)
            {
                non_negative = false;
                break;
            }
            sum_lambda += w;
        }

        if (!non_negative)
        {
            continue;
        }

        // 合計が 1 に近いか確認
        if (std::abs(sum_lambda - 1.0) <= eps)
        {
            barycentric_coordinates = std::move(lambdas);
            return static_cast<int>(simplex_index);
        }
    }

    // 見つからない場合
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
    if (simplex_index < 0 || simplex_index >= static_cast<int>(simplices_.size()))
    {
        return {};
    }
    const auto &indices = simplices_[static_cast<size_t>(simplex_index)];
    if (static_cast<int>(indices.size()) != n_dims_ + 1)
    {
        return {};
    }

    // 参照頂点 V0 を取得
    const std::vector<double> V0 = getPointCoordinates(indices[0]);

    // 行列 A と 右辺 b を構築: A * w = (x - V0), 列 A[:,j] = V_{j+1} - V0
    std::vector<std::vector<double>> A(static_cast<size_t>(n_dims_), std::vector<double>(static_cast<size_t>(n_dims_), 0.0));
    std::vector<double> b(static_cast<size_t>(n_dims_), 0.0);

    for (int d = 0; d < n_dims_; ++d)
    {
        b[static_cast<size_t>(d)] = query_point[static_cast<size_t>(d)] - V0[static_cast<size_t>(d)];
    }

    for (int j = 0; j < n_dims_; ++j)
    {
        const std::vector<double> Vj = getPointCoordinates(indices[static_cast<size_t>(j + 1)]);
        for (int d = 0; d < n_dims_; ++d)
        {
            A[static_cast<size_t>(d)][static_cast<size_t>(j)] = Vj[static_cast<size_t>(d)] - V0[static_cast<size_t>(d)];
        }
    }

    // w を解く
    std::vector<double> w = solveLinearEquation(A, b);
    if (static_cast<int>(w.size()) != n_dims_)
    {
        return {};
    }

    // Barycentric座標に変換
    std::vector<double> lambdas(static_cast<size_t>(n_dims_ + 1), 0.0);
    double sum_w = 0.0;
    for (int j = 0; j < n_dims_; ++j)
    {
        lambdas[static_cast<size_t>(j + 1)] = w[static_cast<size_t>(j)];
        sum_w += w[static_cast<size_t>(j)];
    }
    lambdas[0] = 1.0 - sum_w;
    return lambdas;
}

/**
 * @brief 線形方程式 Ax = b をガウス消去法で解く
 * 
 * 実装アルゴリズム：
 * 1. 拡大係数行列 [A|b] を構築
 * 2. 部分ピボット選択付きガウス消去法を実行
 *    - 各列について最大絶対値要素をピボットに選択
 *    - 行交換で数値安定性を向上
 * 3. 前進消去：下三角行列に変換
 * 4. 後退代入：解ベクトルを計算
 * 
 * ピボット選択の意義：
 * - 小さなピボット要素による除算を避ける
 * - 数値誤差の蓄積を抑制
 * - 特異行列（det(A)≈0）の検出
 * 
 * 計算量とメモリ：
 * - 時間計算量：O(N³)
 * - 空間計算量：O(N²) （拡大係数行列用）
 * 
 * 数値安定性の限界：
 * - 条件数が大きい行列では精度劣化
 * - 完全ピボット選択やLU分解の方が安定な場合あり
 * 
 * @param A 係数行列（N×N）
 * @param b 右辺ベクトル（N要素）
 * @return 解ベクトル x（N要素）、解が存在しない場合は空ベクトル
 * 
 * @note 行列が特異（det(A)=0）の場合は空ベクトルを返す
 * @note 入力検証：行列のサイズ整合性をチェック
 */
std::vector<double> SimpleLinearNDInterpolator::solveLinearEquation(
    const std::vector<std::vector<double>> &A, 
    const std::vector<double> &b) const
{
    const int n = static_cast<int>(A.size());
    if (n == 0 || static_cast<int>(A[0].size()) != n || static_cast<int>(b.size()) != n)
    {
        return {};
    }

    // 拡大係数行列 [A|b]
    std::vector<std::vector<double>> M(n, std::vector<double>(static_cast<size_t>(n + 1), 0.0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            M[static_cast<size_t>(i)][static_cast<size_t>(j)] = A[static_cast<size_t>(i)][static_cast<size_t>(j)];
        }
        M[static_cast<size_t>(i)][static_cast<size_t>(n)] = b[static_cast<size_t>(i)];
    }

    // ガウス消去（部分ピボット）
    for (int col = 0; col < n; ++col)
    {
        // ピボット選択
        int pivot = col;
        double max_abs = std::abs(M[static_cast<size_t>(col)][static_cast<size_t>(col)]);
        for (int i = col + 1; i < n; ++i)
        {
            double v = std::abs(M[static_cast<size_t>(i)][static_cast<size_t>(col)]);
            if (v > max_abs)
            {
                max_abs = v;
                pivot = i;
            }
        }
        if (max_abs == 0.0)
        {
            return {};
        }
        if (pivot != col)
        {
            std::swap(M[static_cast<size_t>(pivot)], M[static_cast<size_t>(col)]);
        }

        // 正規化
        double diag = M[static_cast<size_t>(col)][static_cast<size_t>(col)];
        for (int j = col; j <= n; ++j)
        {
            M[static_cast<size_t>(col)][static_cast<size_t>(j)] /= diag;
        }

        // 他行をゼロ化
        for (int i = 0; i < n; ++i)
        {
            if (i == col) continue;
            double factor = M[static_cast<size_t>(i)][static_cast<size_t>(col)];
            if (factor == 0.0) continue;
            for (int j = col; j <= n; ++j)
            {
                M[static_cast<size_t>(i)][static_cast<size_t>(j)] -= factor * M[static_cast<size_t>(col)][static_cast<size_t>(j)];
            }
        }
    }

    // 解の抽出
    std::vector<double> x(static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n; ++i)
    {
        x[static_cast<size_t>(i)] = M[static_cast<size_t>(i)][static_cast<size_t>(n)];
    }
    return x;
}

/**
 * @brief 指定されたインデックスの点の座標を取得
 * 
 * 実装詳細：
 * - 内部データ構造：points_[d][i] = d次元目のi番目の点の座標
 * - 出力形式：[x, y, z, ...] （次元順の座標ベクトル）
 * 
 * データアクセスパターン：
 * - 次元優先（dimension-major）形式から点優先（point-major）形式に変換
 * - メモリ局所性は最適ではないが、理解しやすいデータ構造
 * 
 * @param point_index 取得したい点のインデックス（0からn_points_-1）
 * @return 指定された点の座標ベクトル（n_dims_要素）
 * 
 * @note 範囲チェックは行わない（呼び出し側で保証）
 * @note パフォーマンス：O(N) （Nは次元数）
 */
std::vector<double> SimpleLinearNDInterpolator::getPointCoordinates(int point_index) const
{
    std::vector<double> p(static_cast<size_t>(n_dims_));
    for (int d = 0; d < n_dims_; ++d)
    {
        p[static_cast<size_t>(d)] = points_[static_cast<size_t>(d)][static_cast<size_t>(point_index)];
    }
    return p;
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
    std::vector<std::vector<double>> values_2d;
    values_2d.reserve(values.size());
    for (const auto& value : values) {
        values_2d.push_back({value});
    }
    return values_2d;
}
