# SimpleLinearNDInterpolator

N次元空間内の散在する点群に対して線形補間を行うC++ライブラリです。Delaunay三角分割（またはN次元単体分割）を使用して、クエリ点における値を補間します。

## 概要

このライブラリは以下の特徴を持ちます：

- **N次元対応**: 2次元から任意の次元の空間での補間が可能
- **線形補間**: 重心座標系を使用した高精度な線形補間
- **Qhull統合**: 高速で信頼性の高いDelaunay三角分割
- **SciPy互換**: PythonのSciPyライブラリとの互換性テスト付き
- **シンプルなAPI**: 直感的で使いやすいインターフェース

## 機能

### 補間アルゴリズム

1. 入力点群に対してDelaunay三角分割を実行
2. クエリ点を含む単体（simplex）を特定
3. 重心座標系を使用して線形補間を実行

### サポートするデータ型

- **スカラー値**: 各点に単一の数値が対応
- **ベクトル値**: 各点に多次元ベクトルが対応（例：RGB色、速度ベクトルなど）

## 要件

- **C++17** 対応コンパイラ
- **CMake 3.10** 以上
- **Qhull 2020.2** ライブラリ（自動ダウンロード・ビルド）

## インストール

### 1. リポジトリのクローン

```bash
git clone <repository-url>
cd SimpleLinearNDInterpolator
```

### 2. ビルド

```bash
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

### 3. テストの実行

```bash
ctest --output-on-failure
```

## 使用方法

### 基本的な使用例

#### スカラー値の補間（2次元）

```cpp
#include "SimpleLinearNDInterpolator.h"

// 2次元平面上の点群とスカラー値
std::vector<std::vector<double>> points = {
    {0.0, 0.0},  // 原点
    {1.0, 0.0},  // x軸上の点
    {0.0, 1.0},  // y軸上の点
    {1.0, 1.0}   // 対角線上の点
};

std::vector<double> values = {0.0, 1.0, 1.0, 2.0};

// 補間器の構築
SimpleLinearNDInterpolator interp(points, values);

// 補間点での値を取得
std::vector<double> query_point = {0.5, 0.5};
std::vector<double> result = interp.interpolate(query_point);
// result[0] = 1.0 (期待値)
```

#### ベクトル値の補間（2次元から3次元ベクトル）

```cpp
// 2次元平面上の点群から3次元ベクトル値を補間
std::vector<std::vector<double>> points = {
    {0.0, 0.0},  // 原点
    {1.0, 0.0},  // x軸上の点
    {0.0, 1.0},  // y軸上の点
    {1.0, 1.0}   // 対角線上の点
};

std::vector<std::vector<double>> values = {
    {1.0, 0.0, 0.0},  // 赤色ベクトル
    {0.0, 1.0, 0.0},  // 緑色ベクトル
    {0.0, 0.0, 1.0},  // 青色ベクトル
    {0.5, 0.5, 0.5}   // グレー色ベクトル
};

SimpleLinearNDInterpolator interpolator(points, values);

std::vector<double> query_point = {0.5, 0.5};
std::vector<double> interpolated_value = interpolator.interpolate(query_point);
// interpolated_value = {0.25, 0.25, 0.25} (期待値)
```

#### 複数点の一括補間

```cpp
std::vector<std::vector<double>> queries = {
    {0.5, 0.5},
    {0.25, 0.75},
    {0.75, 0.25}
};

auto results = interp.interpolate(queries);
// results[i] に queries[i] での補間結果が格納される
```

## API リファレンス

### コンストラクタ

#### `SimpleLinearNDInterpolator(points, values)`

- **points**: 補間に使用する点群（各要素は座標を表すベクトル）
- **values**: 各点に対応する値（スカラー値またはベクトル値）

### メソッド

#### `interpolate(query_point)`

単一のクエリ点での補間値を返します。

- **query_point**: 補間したい点の座標
- **戻り値**: 補間された値のベクトル

#### `interpolate(queries)`

複数のクエリ点での補間値を一括で返します。

- **queries**: 補間したい点群の座標
- **戻り値**: 各クエリ点での補間値のベクトルの配列

## 制限事項

- **最小点数**: 空間の次元数Nに対して、最低でもN+1個の点が必要
- **点の次元**: 全ての点は同じ次元を持つ必要がある
- **値の次元**: 全ての値ベクトルは同じ次元を持つ必要がある
- **スレッドセーフ**: 現在の実装はスレッドセーフではありません

## テスト

### 自動生成テスト

SciPyとの互換性を確認するため、Pythonスクリプトを使用してテストケースを自動生成しています：

```bash
cd scripts
python generate_scipy_test.py
```

### テストの実行

```bash
# 全テストの実行
ctest --output-on-failure

# 特定のテストの実行
./tests/generated_scipy_compat_test
```

## パフォーマンス

- **構築時間**: O(n log n) - Qhullライブラリの効率性に依存
- **補間時間**: O(log n) - 単体探索の効率性に依存
- **メモリ使用量**: O(n) - 点群のサイズに比例

## ライセンス

このプロジェクトのライセンスについては、各ファイルのヘッダーを確認してください。

## 貢献

バグレポート、機能要求、プルリクエストを歓迎します。開発に参加する前に、以下の点を確認してください：

1. コードが既存のスタイルに従っていること
2. 適切なテストが含まれていること
3. ドキュメントが更新されていること

## トラブルシューティング

### よくある問題

#### ビルドエラー

- **C++17対応**: コンパイラがC++17をサポートしていることを確認
- **Qhullライブラリ**: サードパーティライブラリのビルドに失敗した場合、CMakeの設定を確認

#### 実行時エラー

- **点の次元**: 全ての点が同じ次元を持つことを確認
- **最小点数**: 空間の次元数+1個以上の点があることを確認

### サポート

問題が解決しない場合は、以下の情報とともにイシューを作成してください：

- 使用しているOSとコンパイラのバージョン
- エラーメッセージの詳細
- 再現可能な最小限のコード例

## 関連プロジェクト

- [Qhull](http://www.qhull.org/) - 計算幾何学ライブラリ
- [SciPy](https://scipy.org/) - Python科学計算ライブラリ
- [GoogleTest](https://github.com/google/googletest) - C++テストフレームワーク
