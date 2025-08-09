#!/usr/bin/env python3
"""
SciPy 互換テスト生成スクリプト

概要:
- SciPy の LinearNDInterpolator で期待値を計算し、
  同じデータ/クエリ/期待値を埋め込んだ C++ gtest を生成。

デフォルト動作:
- 引数なしで実行すると、いくつかの固定テストケースをまとめた
  `tests/generated_scipy_compat_test.cpp` を生成します。

任意モード（単一ケース生成）:
  python scripts/generate_scipy_test.py --single \
    --dims 2 --num-points 30 --num-queries 40 --seed 42 \
    --outfile tests/generated_scipy_compat_test.cpp

要件:
- numpy, scipy がインストールされていること
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

try:
    import numpy as np
    from scipy.spatial import Delaunay
    from scipy.interpolate import LinearNDInterpolator
except Exception as e:
    sys.stderr.write(
        "SciPy/numpy のインポートに失敗しました。\n"
        "pip でインストールしてください: pip install numpy scipy\n"
        f"詳細: {e}\n"
    )
    sys.exit(1)


def format_float_for_cpp(x: float) -> str:
    """C++のダブルリテラルとして十分な桁で整形。"""
    # NaN の場合はリテラルではなく関数を使う
    if np.isnan(x):
        return "std::numeric_limits<double>::quiet_NaN()"
    # 17桁は double の有効桁数をほぼ保持
    return f"{float(x):.17g}"


def generate_queries_inside_hull(points: np.ndarray, simplices: np.ndarray, num_queries: int, rng: np.random.Generator, include_vertices: int = 0) -> np.ndarray:
    """凸包内部のクエリ点を生成。必要なら元の頂点も一部含める。"""
    n_simplex = simplices.shape[0]
    d = points.shape[1]
    queries: list[np.ndarray] = []

    # 頂点を含めたい数
    include_vertices = int(max(0, min(include_vertices, num_queries, points.shape[0])))
    interior_needed = num_queries - include_vertices

    # まず各単体の重心を可能な限り追加（内部点）
    for s in range(min(n_simplex, interior_needed)):
        verts = points[simplices[s]]  # (d+1, d)
        centroid = verts.mean(axis=0)
        queries.append(centroid)
    # 残りの内部点はランダムに Dirichlet 合成
    while len(queries) < interior_needed:
        sidx = int(rng.integers(0, n_simplex))
        verts = points[simplices[sidx]]  # (d+1, d)
        w = rng.dirichlet(np.ones(d + 1))  # (d+1,)
        q = (w[:, None] * verts).sum(axis=0)
        queries.append(q)

    # 頂点も追加（先頭のいくつか）
    for i in range(include_vertices):
        queries.append(points[i])

    return np.stack(queries, axis=0)


def generate_edge_midpoints(points: np.ndarray, simplices: np.ndarray, num_midpoints: int, rng: np.random.Generator) -> np.ndarray:
    """ランダムな単体の2頂点の中点を生成。凸包内部（または境界上）の点。"""
    if num_midpoints <= 0:
        return np.empty((0, points.shape[1]), dtype=points.dtype)
    mids: list[np.ndarray] = []
    d = points.shape[1]
    for _ in range(num_midpoints):
        sidx = int(rng.integers(0, simplices.shape[0]))
        verts_idx = simplices[sidx]
        if verts_idx.size < 2:
            # 念のため
            mids.append(points[int(rng.integers(0, points.shape[0]))])
            continue
        # 頂点2つを選ぶ
        i1, i2 = rng.choice(verts_idx, size=2, replace=False)
        p1 = points[i1]
        p2 = points[i2]
        mids.append((p1 + p2) * 0.5)
    return np.stack(mids, axis=0).reshape(-1, d)


def generate_facet_points(points: np.ndarray, simplices: np.ndarray, num_points: int, rng: np.random.Generator) -> np.ndarray:
    """ファセット（d-1 次元面）上の点を生成（1つのバリセントリック係数を0に固定）。"""
    if num_points <= 0:
        return np.empty((0, points.shape[1]), dtype=points.dtype)
    outs: list[np.ndarray] = []
    d = points.shape[1]
    for _ in range(num_points):
        sidx = int(rng.integers(0, simplices.shape[0]))
        verts_idx = simplices[sidx]
        if verts_idx.size < d + 1:
            outs.append(points[int(rng.integers(0, points.shape[0]))])
            continue
        # 1頂点を除外し、残りの頂点で Dirichlet を作る（合計1）
        exclude = int(rng.integers(0, verts_idx.size))
        mask = [idx for j, idx in enumerate(verts_idx) if j != exclude]
        w = rng.dirichlet(np.ones(len(mask)))  # (d,)
        verts = points[np.array(mask)]  # (d, d)
        q = (w[:, None] * verts).sum(axis=0)
        outs.append(q)
    return np.stack(outs, axis=0).reshape(-1, d)


def generate_outside_queries(points: np.ndarray, tri: Delaunay, num_points: int, rng: np.random.Generator) -> np.ndarray:
    """凸包外の点を生成（バウンディングボックスを拡張してサンプリングし、find_simplexで判定）。"""
    if num_points <= 0:
        return np.empty((0, points.shape[1]), dtype=points.dtype)
    d = points.shape[1]
    pmin = points.min(axis=0)
    pmax = points.max(axis=0)
    span = pmax - pmin
    scale_list = [0.15, 0.3, 0.6, 1.0]
    collected: list[np.ndarray] = []
    for scale in scale_list:
        if len(collected) >= num_points:
            break
        low = pmin - span * scale
        high = pmax + span * scale
        for _ in range(50 * num_points):
            if len(collected) >= num_points:
                break
            x = rng.uniform(low, high)
            if tri.find_simplex(x) < 0:
                collected.append(x)
    if not collected:
        # フォールバック: 中心から大きくオフセット
        center = points.mean(axis=0)
        for _ in range(num_points):
            direction = rng.standard_normal(size=d)
            direction /= np.linalg.norm(direction) + 1e-12
            collected.append(center + direction * (np.linalg.norm(span) + 1.0))
    return np.stack(collected[:num_points], axis=0)


def build_cpp_test(points: np.ndarray, values: np.ndarray, queries: np.ndarray, expected: np.ndarray, test_name: str = "ScipyCompat_Generated") -> str:
    d = points.shape[1]
    n = points.shape[0]
    m = queries.shape[0]

    # 次元優先（dimension-major）の配列に変換
    coords_by_dim = [points[:, dim].tolist() for dim in range(d)]

    # 文字列化
    cpp_points_lines = []
    for dim_idx, coord in enumerate(coords_by_dim):
        numbers = ", ".join(format_float_for_cpp(x) for x in coord)
        cpp_points_lines.append(f"    /* dim {dim_idx} */ {{ {numbers} }}")
    cpp_points_body = ",\n".join(cpp_points_lines)

    cpp_values = ", ".join(format_float_for_cpp(x) for x in values.tolist())

    cpp_queries_lines = []
    for i in range(m):
        numbers = ", ".join(format_float_for_cpp(x) for x in queries[i].tolist())
        cpp_queries_lines.append(f"    /* q{i} */ {{ {numbers} }}")
    cpp_queries_body = ",\n".join(cpp_queries_lines)

    cpp_expected = ", ".join(format_float_for_cpp(x) for x in expected.tolist())

    code = f"""
#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits>
#include "SimpleLinearNDInterpolator.h"

// このファイルは scripts/generate_scipy_test.py により自動生成されています。

TEST(ScipyCompat, {test_name}) {{
  // {d} 次元, 入力点 {n} 個, クエリ {m} 個
  std::vector<std::vector<double>> points = {{
{cpp_points_body}
  }};
  std::vector<double> values = {{ {cpp_values} }};

  std::vector<std::vector<double>> queries = {{
{cpp_queries_body}
  }};

  SimpleLinearNDInterpolator interp(points, values);
  auto results = interp.interpolate(queries);

  ASSERT_EQ(results.size(), queries.size());

  const std::vector<double> expected = {{ {cpp_expected} }};
  const double tol = 1e-9;
  for (size_t i = 0; i < queries.size(); ++i) {{
    ASSERT_EQ(results[i].size(), static_cast<size_t>(1));
    const double got = results[i][0];
    const double expv = expected[i];
    if (std::isnan(expv)) {{
      ASSERT_TRUE(std::isnan(got));
    }} else {{
      ASSERT_NEAR(got, expv, tol);
    }}
  }}
}}
"""
    return code


def build_cpp_test_case(points: np.ndarray, values: np.ndarray, queries: np.ndarray, expected: np.ndarray, test_name: str) -> str:
    """インクルードを含まない TEST 本体だけを生成。"""
    d = points.shape[1]
    n = points.shape[0]
    m = queries.shape[0]

    coords_by_dim = [points[:, dim].tolist() for dim in range(d)]

    cpp_points_lines = []
    for dim_idx, coord in enumerate(coords_by_dim):
        numbers = ", ".join(format_float_for_cpp(x) for x in coord)
        cpp_points_lines.append(f"    /* dim {dim_idx} */ {{ {numbers} }}")
    cpp_points_body = ",\n".join(cpp_points_lines)

    cpp_values = ", ".join(format_float_for_cpp(x) for x in values.tolist())

    cpp_queries_lines = []
    for i in range(m):
        numbers = ", ".join(format_float_for_cpp(x) for x in queries[i].tolist())
        cpp_queries_lines.append(f"    /* q{i} */ {{ {numbers} }}")
    cpp_queries_body = ",\n".join(cpp_queries_lines)

    cpp_expected = ", ".join(format_float_for_cpp(x) for x in expected.tolist())

    return f"""
TEST(ScipyCompat, {test_name}) {{
  // {d} 次元, 入力点 {n} 個, クエリ {m} 個
  std::vector<std::vector<double>> points = {{
{cpp_points_body}
  }};
  std::vector<double> values = {{ {cpp_values} }};

  std::vector<std::vector<double>> queries = {{
{cpp_queries_body}
  }};

  SimpleLinearNDInterpolator interp(points, values);
  auto results = interp.interpolate(queries);

  ASSERT_EQ(results.size(), queries.size());

  const std::vector<double> expected = {{ {cpp_expected} }};
  const double tol = 1e-9;
  for (size_t i = 0; i < queries.size(); ++i) {{
    ASSERT_EQ(results[i].size(), static_cast<size_t>(1));
    const double got = results[i][0];
    const double expv = expected[i];
    if (std::isnan(expv)) {{
      ASSERT_TRUE(std::isnan(got));
    }} else {{
      ASSERT_NEAR(got, expv, tol);
    }}
  }}
}}
"""


def build_cpp_file(test_cases: list[str]) -> str:
    header = """
#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits>
#include "SimpleLinearNDInterpolator.h"

// このファイルは scripts/generate_scipy_test.py により自動生成されています。

"""
    return header + "\n\n".join(test_cases) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="SciPy互換テスト(C++)自動生成")
    parser.add_argument("--single", action="store_true", help="単一ケースのみ生成する（従来の引数指定モード）")
    parser.add_argument("--dims", type=int, default=2, help="空間次元 D（--single のとき有効）")
    parser.add_argument("--num-points", type=int, default=30, help="入力点数 N (>= D+1)（--single のとき有効）")
    parser.add_argument("--num-queries", type=int, default=40, help="クエリ点数 M（--single のとき有効）")
    parser.add_argument("--seed", type=int, default=42, help="乱数シード（--single のとき有効）")
    parser.add_argument("--outfile", type=str, default="", help="出力 C++ テストファイルパス（未指定なら tests/generated_scipy_compat_test.cpp）")
    args = parser.parse_args()

    # 出力先
    out_path = Path(args.outfile) if args.outfile else (Path(__file__).resolve().parents[1] / "tests" / "generated_scipy_compat_test.cpp")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if args.single:
        if args.dims < 1:
            raise SystemExit("--dims は 1 以上である必要があります")
        if args.num_points < args.dims + 1:
            raise SystemExit("--num-points は D+1 以上である必要があります")

        rng = np.random.default_rng(args.seed)
        points = rng.random((args.num_points, args.dims))
        values = rng.standard_normal(size=(args.num_points,))

        tri = Delaunay(points, qhull_options='Qbb Qc Qz')
        interp = LinearNDInterpolator(tri, values)
        if tri.simplices.size == 0:
            raise SystemExit("Delaunay 分割が生成できませんでした。点配置や次元を見直してください。")

        # 単一モードでもエッジケースを少量混ぜる
        quota_vertices = max(0, args.num_queries // 8)
        quota_mid = max(0, args.num_queries // 8)
        quota_facet = max(0, args.num_queries // 8)
        quota_outside = max(0, args.num_queries // 8)
        interior = max(0, args.num_queries - (quota_vertices + quota_mid + quota_facet + quota_outside))

        q_interior = generate_queries_inside_hull(points, tri.simplices, interior, rng, include_vertices=0)
        q_vertices = points[:min(args.num_points, quota_vertices)]
        q_mid = generate_edge_midpoints(points, tri.simplices, quota_mid, rng)
        q_facet = generate_facet_points(points, tri.simplices, quota_facet, rng)
        q_out = generate_outside_queries(points, tri, quota_outside, rng)
        queries = np.concatenate([q_interior, q_vertices, q_mid, q_facet, q_out], axis=0)
        expected = interp(queries)
        expected = np.asarray(expected).reshape(-1)

        code = build_cpp_test(points, values, queries, expected)
        out_path.write_text(code, encoding="utf-8")
        print(f"生成しました: {out_path}")
        return

    # 固定テストケース群（引数なしデフォルト）
    fixed_cases = [
        # (dims, num_points, num_queries, seed)
        # 2D
        (2, 20, 40, 1001),
        (2, 28, 40, 1002),
        (2, 35, 48, 1003),
        (2, 50, 60, 1004),
        # 3D
        (3, 24, 36, 2001),
        (3, 30, 36, 2002),
        (3, 45, 48, 2003),
        (3, 60, 60, 2004),
        # 4D
        (4, 24, 32, 3001),
        (4, 35, 40, 3002),
        (4, 50, 48, 3003),
        # 5D（少し小さめに）
        (5, 18, 24, 4001),
    ]

    test_blocks: list[str] = []
    for dims, num_points, num_queries, seed in fixed_cases:
        rng = np.random.default_rng(seed)
        points = rng.random((num_points, dims))
        values = rng.standard_normal(size=(num_points,))

        tri = Delaunay(points, qhull_options='Qbb Qc Qz')
        if tri.simplices.size == 0:
            # 不運により退化したらスキップ
            print(f"警告: Delaunay 生成失敗のためスキップ: D={dims}, N={num_points}, seed={seed}")
            continue
        interp = LinearNDInterpolator(tri, values)

        # 内部点 + 一部頂点 + 辺中点 + ファセット上 + 凸包外
        quota_vertices = max(0, num_queries // 6)
        quota_mid = max(0, num_queries // 6)
        quota_facet = max(0, num_queries // 6)
        quota_outside = max(0, num_queries // 6)
        interior = max(0, num_queries - (quota_vertices + quota_mid + quota_facet + quota_outside))
        include_vertices = min(points.shape[0], quota_vertices)

        q_interior = generate_queries_inside_hull(points, tri.simplices, interior, rng, include_vertices=0)
        q_vertices = points[:include_vertices]
        q_mid = generate_edge_midpoints(points, tri.simplices, quota_mid, rng)
        q_facet = generate_facet_points(points, tri.simplices, quota_facet, rng)
        q_out = generate_outside_queries(points, tri, quota_outside, rng)
        queries = np.concatenate([q_interior, q_vertices, q_mid, q_facet, q_out], axis=0)
        expected = interp(queries)
        expected = np.asarray(expected).reshape(-1)

        name = f"Fixed_D{dims}_N{num_points}_Q{num_queries}_S{seed}"
        test_blocks.append(build_cpp_test_case(points, values, queries, expected, name))

    if not test_blocks:
        raise SystemExit("固定テストの生成に失敗しました（全ケースで Delaunay 失敗）。再実行してください。")

    out_code = build_cpp_file(test_blocks)
    out_path.write_text(out_code, encoding="utf-8")
    print(f"固定テストを生成しました: {out_path}  (ケース数: {len(test_blocks)})")


if __name__ == "__main__":
    main()

