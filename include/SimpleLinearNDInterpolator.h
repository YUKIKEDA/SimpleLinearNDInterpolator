#include <vector>
#include <memory>
#include <optional>

namespace orgQhull {
    class Qhull;
}

/**
 * @brief N次元線形補間器クラス
 * 
 * このクラスは、N次元空間内の散在する点群に対して線形補間を行います。
 * Delaunay三角分割（またはN次元単体分割）を使用して、クエリ点における
 * 値を補間します。
 * 
 * 補間アルゴリズム：
 * 1. 入力点群に対してDelaunay三角分割を実行
 * 2. クエリ点を含む単体（simplex）を特定
 * 3. 重心座標系を使用して線形補間を実行
 * 
 * @note このクラスはQhullライブラリを使用してDelaunay三角分割を実行します
 * @note スレッドセーフではありません
 * 
 * 使用例：
 * @code
 * std::vector<std::vector<double>> points = {{0,0}, {1,0}, {0,1}, {1,1}};
 * std::vector<double> values = {0.0, 1.0, 1.0, 2.0};
 * 
 * SimpleLinearNDInterpolator interp(points, values);
 * 
 * std::vector<double> result = interp.interpolate({0.5, 0.5});
 * @endcode
 */
class SimpleLinearNDInterpolator
{
public:
    /**
     * @brief 多次元値を持つ点群から補間器を構築するコンストラクタ
     * 
     * このコンストラクタは、各点に対して複数の値（ベクトル値）が定義されている
     * 場合に使用します。例えば、2次元平面上の各点に3次元ベクトルが対応する場合などです。
     * 
     * @param points 補間に使用する点群。各要素は座標を表すベクトル
     *               例: {{x1,y1}, {x2,y2}, ...} for 2D points
     * @param values 各点に対応する値のベクトル。points[i]にvalues[i]が対応
     *               例: {{vx1,vy1,vz1}, {vx2,vy2,vz2}, ...} for 3D vectors
     * @param enable_degeneracy_fallback 縮退した点群でフォールバック処理を有効にするか（デフォルト: false）
     * 
     * @throw std::invalid_argument points.size() != values.size()の場合
     * @throw std::invalid_argument pointsが空の場合
     * @throw std::invalid_argument 点の次元が不一致の場合
     * @throw std::invalid_argument 点数不足または縮退（enable_degeneracy_fallback=falseの場合）
     * @throw std::runtime_error Delaunay三角分割の構築に失敗した場合
     * 
     * @note 最低でもN+1個の点が必要です（Nは空間の次元数）
     * @note 全ての点は同じ次元を持つ必要があります
     * @note 全ての値ベクトルは同じ次元を持つ必要があります
     * 
     * @example
     * std::vector<std::vector<double>> points = { // 2次元平面上の点群から3次元ベクトル値を補間する例
     *     {0.0, 0.0},  // 原点
     *     {1.0, 0.0},  // x軸上の点
     *     {0.0, 1.0},  // y軸上の点
     *     {1.0, 1.0}   // 対角線上の点
     * };
     * 
     * std::vector<std::vector<double>> values = {
     *     {1.0, 0.0, 0.0},  // 赤色ベクトル
     *     {0.0, 1.0, 0.0},  // 緑色ベクトル
     *     {0.0, 0.0, 1.0},  // 青色ベクトル
     *     {0.5, 0.5, 0.5}   // グレー色ベクトル
     * };
     * 
     * SimpleLinearNDInterpolator interpolator(points, values);
     * 
     * std::vector<double> query_point = {0.5, 0.5}; // 補間点での値を取得
     * std::vector<double> interpolated_value = interpolator.interpolate(query_point);
     */
    SimpleLinearNDInterpolator(
        const std::vector<std::vector<double>> &points, 
        const std::vector<std::vector<double>> &values,
        bool enable_degeneracy_fallback = false
    );

    /**
     * @brief スカラー値を持つ点群から補間器を構築するコンストラクタ
     * 
     * このコンストラクタは、各点に対してスカラー値（単一の数値）が定義されている
     * 場合に使用します。内部的には1次元のベクトル値として扱われます。
     * 
     * @param points 補間に使用する点群。各要素は座標を表すベクトル
     *               例: {{x1,y1}, {x2,y2}, ...} for 2D points
     * @param values 各点に対応するスカラー値。points[i]にvalues[i]が対応
     *               例: {v1, v2, v3, ...} for scalar values
     * @param enable_degeneracy_fallback 縮退した点群でフォールバック処理を有効にするか（デフォルト: false）
     * 
     * @throw std::invalid_argument points.size() != values.size()の場合
     * @throw std::invalid_argument pointsまたはvaluesが空の場合
     * @throw std::invalid_argument 点の次元が不一致の場合
     * @throw std::invalid_argument 点数不足または縮退（enable_degeneracy_fallback=falseの場合）
     * @throw std::runtime_error Delaunay三角分割の構築に失敗した場合
     * 
     * @note 最低でもN+1個の点が必要です（Nは空間の次元数）
     * @note 全ての点は同じ次元を持つ必要があります
     * 
     * @example
     * std::vector<std::vector<double>> points = { // 2次元平面上の点群からスカラー値を補間する例（温度分布など）
     *     {0.0, 0.0},  // 原点
     *     {1.0, 0.0},  // x軸上の点
     *     {0.0, 1.0},  // y軸上の点
     *     {1.0, 1.0}   // 対角線上の点
     * };
     * 
     * std::vector<double> temperatures = {
     *     20.0,  // 原点の温度
     *     25.0,  // x軸上の温度
     *     22.0,  // y軸上の温度
     *     28.0   // 対角線上の温度
     * };
     * 
     * SimpleLinearNDInterpolator interpolator(points, temperatures);
     * 
     * std::vector<double> query_point = {0.5, 0.5}; 補間点での温度を取得
     * std::vector<double> interpolated_temp = interpolator.interpolate(query_point);
     * double temperature = interpolated_temp[0];  // スカラー値として取得
     */
    SimpleLinearNDInterpolator(
        const std::vector<std::vector<double>> &points, 
        const std::vector<double> &values,
        bool enable_degeneracy_fallback = false
    );

    /**
     * @brief デストラクタ
     * 
     * Qhullオブジェクトなどの内部リソースを適切に解放します。
     */
    ~SimpleLinearNDInterpolator();

    /**
     * @brief 複数のクエリ点に対してベクトル値の補間を実行
     * 
     * 複数のクエリ点に対して同時に補間を実行し、各点における補間値を返します。
     * 各クエリ点は、補間器構築時に使用した点群と同じ次元を持つ必要があります。
     * 
     * @param query_points 補間を行いたい点群。各要素は座標を表すベクトル
     *                     例: {{qx1,qy1}, {qx2,qy2}, ...} for 2D query points
     * @param use_nearest_neighbor_fallback 単体が見つからない場合に最近傍補間を使用するか
     * 
     * @return 各クエリ点における補間値のベクトル。
     *         query_points[i]に対してreturn_value[i]が対応。
     *         各補間値は、補間器構築時の値と同じ次元のベクトル。
     * 
     * @throw std::invalid_argument クエリ点の次元が不正な場合
     * @throw std::runtime_error 補間に失敗した場合（use_nearest_neighbor_fallback=falseかつ点が凸包の外部にある場合）
     * 
     * @note use_nearest_neighbor_fallback=falseの場合、凸包外の点はNaNが返される
     * @note use_nearest_neighbor_fallback=trueの場合、凸包外では最近傍点の値を返す
     * @note 最近傍補間の計算量はO(M * N)です。Mはクエリ点数、Nは入力点数
     * 
     * 使用例：
     * @code
     * std::vector<std::vector<double>> queries = {{0.3, 0.7}, {0.8, 0.2}};
     * auto results = interpolator.interpolate(queries, true);
     * @endcode
     */
    std::vector<std::vector<double>> interpolate(
        const std::vector<std::vector<double>> &query_points,
        bool use_nearest_neighbor_fallback = false
    ) const;

    /**
     * @brief 単一のクエリ点に対して補間を実行
     * 
     * 指定されたクエリ点における補間値を計算して返します。
     * クエリ点は、補間器構築時に使用した点群と同じ次元を持つ必要があります。
     * 
     * @param query_points 補間を行いたい点の座標
     *                     例: {qx, qy} for 2D query point
     * @param use_nearest_neighbor_fallback 単体が見つからない場合に最近傍補間を使用するか
     * 
     * @return クエリ点における補間値。
     *         スカラー値で構築された場合は1要素のベクトル、
     *         ベクトル値で構築された場合は対応する次元のベクトル。
     * 
     * @throw std::invalid_argument クエリ点の次元が不正な場合
     * @throw std::runtime_error 補間に失敗した場合（use_nearest_neighbor_fallback=falseかつ点が凸包の外部にある場合）
     * 
     * @note use_nearest_neighbor_fallback=falseの場合、凸包外の点はNaNが返される
     * @note use_nearest_neighbor_fallback=trueの場合、凸包外では最近傍点の値を返す
     * @note 最近傍補間の計算量はO(N)です。Nは入力点数
     * 
     * 使用例：
     * @code
     * std::vector<double> result = interpolator.interpolate({0.5, 0.5}, true);
     * @endcode
     */
    std::vector<double> interpolate(
        const std::vector<double> &query_points,
        bool use_nearest_neighbor_fallback = false
    ) const;

private:    
    /** @brief 補間に使用する点群の座標データ */
    std::vector<std::vector<double>> points_;
    
    /** @brief 各点に対応する値データ（ベクトル形式） */
    std::vector<std::vector<double>> values_;

    /** @brief 点群の次元数（2D、3Dなど） */
    int n_dims_;
    
    /** @brief 点群の総数 */
    int n_points_;

    /** @brief Delaunay三角分割を行うQhullオブジェクト（2次元以上の場合のみ使用） */
    std::optional<std::unique_ptr<orgQhull::Qhull>> qhull_;
    
    /** @brief 三角分割で生成された単体（simplex）のリスト。各単体は点のインデックスで表現（2次元以上の場合のみ使用） */
    std::optional<std::vector<std::vector<int>>> simplices_;

    /** @brief 点群が縮退しているかどうかのフラグ */
    bool is_degenerate_;
    
    /** @brief 縮退している場合の有効次元数（実際にデータが分散している次元数） */
    int effective_dimensions_;
    
    /** @brief 縮退している場合の射影行列（高次元から有効次元への射影用） */
    std::optional<std::vector<std::vector<double>>> projection_matrix_;
    
    /** @brief 縮退している場合の射影された点群データ */
    std::optional<std::vector<std::vector<double>>> projected_points_;
    
    /** @brief 縮退している場合の射影空間での補間器 */
    std::optional<std::unique_ptr<SimpleLinearNDInterpolator>> projected_interpolator_;

    /**
     * @brief 入力点群に対してDelaunay三角分割を実行
     * 
     * Qhullライブラリを使用して点群のDelaunay三角分割を計算し、
     * 内部のqhull_オブジェクトに結果を格納します。
     * 
     * @param points 三角分割を行う点群
     * @throw std::runtime_error 三角分割の構築に失敗した場合
     */
    void buildTriangulation(const std::vector<std::vector<double>> &points);
    
    /**
     * @brief Qhullの結果から単体（simplex）のリストを構築
     * 
     * Qhullによって生成された三角分割結果から、各単体を構成する
     * 点のインデックスリストを抽出してsimplices_に格納します。
     */
    void buildSimplexList();

    /**
     * @brief クエリ点を含む単体を特定し、重心座標を計算
     * 
     * 指定されたクエリ点を含む単体（simplex）を検索し、
     * その単体内での重心座標（barycentric coordinates）を計算します。
     * 
     * @param query_point 検索対象のクエリ点
     * @param barycentric_coordinates 計算された重心座標の出力先
     * @param eps 数値計算の許容誤差
     * @return 見つかった単体のインデックス。見つからない場合は-1
     */
    int findSimplex(
        const std::vector<double> &query_point, 
        std::vector<double>& barycentric_coordinates, 
        double eps
    ) const;

    /**
     * @brief 指定された単体内でのクエリ点の重心座標を計算
     * 
     * 与えられた単体内でのクエリ点の重心座標を線形方程式を解くことで計算します。
     * 重心座標は補間の重みとして使用されます。
     * 
     * @param query_point 重心座標を計算したいクエリ点
     * @param simplex_index 対象となる単体のインデックス
     * @return 計算された重心座標のベクトル
     */
    std::vector<double> calculateBarycentricCoordinates(
        const std::vector<double> &query_point, 
        int simplex_index
    ) const;

    /**
     * @brief 線形方程式 Ax = b を解く
     * 
     * 重心座標の計算で使用される線形方程式を解きます。
     * ガウス消去法やLU分解などの手法を使用します。
     * 
     * @param A 係数行列
     * @param b 右辺ベクトル
     * @return 解ベクトル x
     * @throw std::runtime_error 行列が特異（解が存在しない）の場合
     */
    std::vector<double> solveLinearEquation(
        const std::vector<std::vector<double>> &A, 
        const std::vector<double> &b
    ) const;

    /**
     * @brief 指定されたインデックスの点の座標を取得
     * 
     * @param point_index 取得したい点のインデックス
     * @return 指定された点の座標ベクトル
     */
    std::vector<double> getPointCoordinates(int point_index) const;

    /**
     * @brief 最近傍点のインデックスを見つける
     * 
     * @param query_point 検索対象のクエリ点
     * @return 最近傍点のインデックス
     */
    int findNearestNeighbor(const std::vector<double> &query_point) const;

    /**
     * @brief 2点間のユークリッド距離の2乗を計算
     * 
     * @param point1 第1の点の座標
     * @param point2 第2の点の座標
     * @return 2点間のユークリッド距離の2乗
     */
    double calculateDistanceSquared(const std::vector<double> &point1, const std::vector<double> &point2) const;

    /**
     * @brief 1次元ベクトル値を2次元ベクトル形式に変換
     * 
     * スカラー値のベクトルを、内部で統一的に扱うための
     * 2次元ベクトル形式（各要素が1要素のベクトル）に変換します。
     * 
     * @param values 変換元の1次元ベクトル
     * @return 変換後の2次元ベクトル
     */
    std::vector<std::vector<double>> convertTo2DVector(
        const std::vector<double> &values
    ) const;

    /**
     * @brief 点群の縮退状況を詳細に分析し、フォールバック情報を取得
     * 
     * 点群の縮退状況を分析し、有効次元数や射影行列などの
     * フォールバック補間に必要な情報を取得します。
     * 
     * @param points 分析対象の点群
     * @param rank_tolerance ランク判定の数値許容誤差
     * @param effective_dims [出力] 有効次元数
     * @param projection_matrix [出力] 射影行列（高次元→有効次元）
     * @return true: 縮退していない, false: 縮退している
     * 
     * @note 射影行列はSVDのV行列から構築される
     * @note 有効次元での補間が可能になる
     */
    bool analyzeDegeneracy(
        const std::vector<std::vector<double>> &points,
        double rank_tolerance,
        int &effective_dims,
        std::vector<std::vector<double>> &projection_matrix
    ) const;

    /**
     * @brief 縮退した点群に対する射影補間器を設定
     * 
     * 縮退した点群を有効次元の部分空間に射影し、
     * その空間での補間器を構築します。
     * 
     * @param points 元の点群
     * @param values 対応する値
     */
    void setupProjectedInterpolation(
        const std::vector<std::vector<double>> &points,
        const std::vector<std::vector<double>> &values
    );

    /**
     * @brief 縮退した点群に対するフォールバック補間を実行
     * 
     * 縮退タイプに応じて最適な補間方法を選択：
     * - 有効次元 > 0: 射影補間
     * - 有効次元 = 0: 最近傍補間
     * 
     * @param query_points 補間対象のクエリ点群
     * @param use_nearest_neighbor_fallback 最近傍フォールバックの使用
     * @return 補間結果
     */
    std::vector<std::vector<double>> interpolateWithDegenerateFallback(
        const std::vector<std::vector<double>> &query_points,
        bool use_nearest_neighbor_fallback
    ) const;

    
    /**
     * @brief 1次元線形補間を実行
     * 
     * 1次元空間での線形補間を行います。2つの隣接する点の間で線形補間を実行し、
     * 範囲外の場合は最近傍補間または外挿を行います。
     * 
     * @param query_point 1次元のクエリ点（1要素のベクトル）
     * @param use_nearest_neighbor_fallback 範囲外の場合に最近傍補間を使用するか
     * @return 補間された値のベクトル
     */
    std::vector<double> interpolate1D(
        const std::vector<double> &query_point,
        bool use_nearest_neighbor_fallback = false
    ) const;

    /**
     * @brief 1次元空間での最近傍補間を実行
     * 
     * 1次元空間でクエリ点に最も近い点を見つけ、その点の値を返します。
     * 
     * @param query_x クエリ点のx座標
     * @param sorted_indices x座標でソートされた点のインデックス
     * @return 最近傍点の値ベクトル
     */
    std::vector<double> findNearestNeighbor1D(
        double query_x,
        const std::vector<int> &sorted_indices
    ) const;

    /**
     * @brief 2次元配列が矩形であるかをチェック
     * 
     * @param m チェック対象の2次元配列
     * @return 矩形であれば true、そうでなければ false
     */
    static bool isRectangular(const std::vector<std::vector<double>> &m);
};