#include <vector>
#include <memory>

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
     * 
     * @throw std::invalid_argument points.size() != values.size()の場合
     * @throw std::invalid_argument pointsが空の場合
     * @throw std::invalid_argument 点の次元が不一致の場合
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
        const std::vector<std::vector<double>> &values
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
     * 
     * @throw std::invalid_argument points.size() != values.size()の場合
     * @throw std::invalid_argument pointsまたはvaluesが空の場合
     * @throw std::invalid_argument 点の次元が不一致の場合
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
        const std::vector<double> &values
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
     * 
     * @return 各クエリ点における補間値のベクトル。
     *         query_points[i]に対してreturn_value[i]が対応。
     *         各補間値は、補間器構築時の値と同じ次元のベクトル。
     * 
     * @throw std::invalid_argument クエリ点の次元が不正な場合
     * @throw std::runtime_error 補間に失敗した場合（例：点が凸包の外部にある場合）
     * 
     * @note クエリ点が入力点群の凸包外部にある場合、最も近い境界での外挿が行われます
     * @note 計算量はO(M * log(N))です。Mはクエリ点数、Nは入力点数
     * 
     * 使用例：
     * @code
     * std::vector<std::vector<double>> queries = {{0.3, 0.7}, {0.8, 0.2}};
     * auto results = interpolator.interpolate(queries);
     * @endcode
     */
    std::vector<std::vector<double>> interpolate(
        const std::vector<std::vector<double>> &query_points
    ) const;

    /**
     * @brief 単一のクエリ点に対して補間を実行
     * 
     * 指定されたクエリ点における補間値を計算して返します。
     * クエリ点は、補間器構築時に使用した点群と同じ次元を持つ必要があります。
     * 
     * @param query_points 補間を行いたい点の座標
     *                     例: {qx, qy} for 2D query point
     * 
     * @return クエリ点における補間値。
     *         スカラー値で構築された場合は1要素のベクトル、
     *         ベクトル値で構築された場合は対応する次元のベクトル。
     * 
     * @throw std::invalid_argument クエリ点の次元が不正な場合
     * @throw std::runtime_error 補間に失敗した場合（例：点が凸包の外部にある場合）
     * 
     * @note クエリ点が入力点群の凸包外部にある場合、最も近い境界での外挿が行われます
     * @note 計算量はO(log(N))です。Nは入力点数
     * 
     * 使用例：
     * @code
     * std::vector<double> result = interpolator.interpolate({0.5, 0.5});
     * @endcode
     */
    std::vector<double> interpolate(
        const std::vector<double> &query_points
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

    /** @brief Delaunay三角分割を行うQhullオブジェクト */
    std::unique_ptr<orgQhull::Qhull> qhull_;
    
    /** @brief 三角分割で生成された単体（simplex）のリスト。各単体は点のインデックスで表現 */
    std::vector<std::vector<int>> simplices_;

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
     * @brief 2次元配列が矩形であるかをチェック
     * 
     * @param m チェック対象の2次元配列
     * @return 矩形であれば true、そうでなければ false
     */
    static bool isRectangular(const std::vector<std::vector<double>> &m);
};