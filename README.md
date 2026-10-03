[![MIT License](https://custom-icon-badges.herokuapp.com/badge/license-MIT-8BB80A.svg?logo=law&logoColor=white)](LICENSE)
[![R](https://custom-icon-badges.herokuapp.com/badge/R-198CE7.svg?logo=R&logoColor=white)]()
[![COVIVIS](https://img.shields.io/badge/COVIVIS-v2.0-CCCCCC?link=https%3A%2F%2Fcovivis.soken.ac.jp%2F)](https://covivis.soken.ac.jp/)

# COVIVIS: 疫学データと下水ウイルス量の解析関数

<p align="center">
  <img src="COVIVIS_logo.png" alt="COVIVIS_logo" width="640">
</p>

`EpiPredFunctions20260819.R` は、感染症サーベイランスにおける発症日推定、疫学データと下水中ウイルス濃度の対応づけ、および予測精度の評価に用いるR関数集です。

主に、報告遅延を考慮した発症日数の復元、回帰モデル、排出（shedding）プロファイルモデルを提供します。

## 前提条件

Rで次のパッケージを利用します。

```r
install.packages(c("dplyr", "magrittr", "nleqslv", "nloptr", "nnls"))

library(dplyr)
library(magrittr)
library(nleqslv)
library(nloptr)
library(nnls)
```

関数を読み込むには、Rコンソールまたはスクリプトで次を実行します。

```r
source("EpiPredFunctions20260819.R")
```

## データ形式

多くの関数は、先頭列が `Date` 型の日付、2列目が数値の `data.frame` を入力として想定します。関数内で列名を付け替えるものがあるため、基本的には日付列を1列目、値の列を2列目に置いてください。

| 用途 | 想定する列 |
| --- | --- |
| 報告数・発症数・感染者数 | `date`, 件数 |
| 定点データ | `date`, 定点当たりの値 |
| 下水データ | `date`, ウイルス濃度 |
| 回帰用データ | `date`, 正の数値 |

対数変換を行う関数では、入力値は正である必要があります。欠測値は関数によって除外または処理されます。

## 関数一覧

### 1. 報告日から発症日への復元

#### `find_me_from_sd(xmean, xsd)`

報告遅延をWeibull分布で表すとき、遅延日数の平均 `xmean` と標準偏差 `xsd` から、Weibull分布の形状（shape）と尺度（scale）を数値的に求めます。

- **戻り値**: `c(shape, scale)`
- **使用先**: `generate_onset()`、`generate_onset_sentinel()`

<br>

#### `generate_onset(reporteddata, xmu, xsd)`

日別の報告件数から、Weibull分布に従う報告遅延をモンテカルロ法で逆算し、推定発症日別件数を生成します。

- `reporteddata`: 日付と報告件数の2列からなるデータ
- `xmu`, `xsd`: 報告遅延日数の平均と標準偏差
- **戻り値**: `date` と `onset`（推定発症数）のデータフレーム

<br>

#### `generate_onset_sentinel(sentineldata, xmu, xsd)`

定点サーベイランスの値を一時的に100倍または1000倍して整数的に扱い、`generate_onset()` と同様の方法で発症日別の推定値を返します。計算後に元のスケールへ戻します。

- `sentineldata`: 日付と定点値の2列からなるデータ
- **戻り値**: `date` と `onset` のデータフレーム

<br>

### 2. 時系列データの整形・区間推定

#### `lerp_data(timeseries)`

週次の定点データを日次系列へ変換します。2列目の値を7で割り、日付の間を線形補間します。

- `timeseries`: 日付と週次値の2列からなるデータ
- **戻り値**: 日付ごとの `date`, `value` を持つデータフレーム

<br>

#### `calc_ts_CI(obs, pred, n_param = 3, conf_level = 0.95)`

観測値と予測値の残差から、残差分散、元の標本サイズまたはラグ1自己相関に基づく有効標本サイズ、信頼区間（CI）と予測区間（PI）の幅を算出します。対数スケールのデータを想定しています。

- `obs`, `pred`: 同じ長さの観測値と予測値
- `n_param`: 推定したパラメータ数（既定値: 3）
- `conf_level`: 信頼水準（既定値: 0.95）
- **戻り値**: `ve`, `n_eff`, `rho1`, `ci_width`, `pi_width` などを含むリスト

<br>

#### `plot_result_cipi(dataframe, ci, pi)`

対数スケールの予測値に、CIおよびPIの上下限列を付与します。

- `dataframe`: 日付と予測値の2列からなるデータ
- **戻り値**: `date`, `predicted.log`, `CI_lower`, `CI_upper`, `PI_lower`, `PI_upper`

<br>

#### `plot_result_cipi_nonlog(dataframe, ci, pi)`

対数10スケールの予測値と区間幅を、元のスケールへ変換して可視化用データを作成します。

- `dataframe`: 日付と対数10予測値の2列からなるデータ
- **戻り値**: `date`, `predicted`, `CI_lower`, `CI_upper`, `PI_lower`, `PI_upper`

<br>

#### `weekly_prediction(weeklyobserved, dailypredicted)`

日次の予測感染者数を、観測された各週次日付までの7日間で合計します。

- `weeklyobserved`: `date` 列を含む週次観測データ
- `dailypredicted`: `date`, `estimated_cases` 列を含む日次予測データ
- **戻り値**: `date`, `pred_sum` を持つデータフレーム

<br>

### 3. 回帰と決定係数

#### `Rsquared_in_lm(xdata, ydata, pa, pb)`

あらかじめ得た回帰係数 `pa`（切片）と `pb`（傾き）を用いて、日付で対応づけた `xdata` と `ydata` の対数10スケール上の残差平方和・回帰平方和・全平方和・2種類の決定係数を計算します。

- **戻り値**: `c(SSE, SSR, SST, R2_1, R2_2)`
- `R2_1 = 1 - SSE / SST`
- `R2_2 = SSR / SST`

<br>

#### `Rsquared(observeddata, predicteddata)`

観測値と予測値から、SSE、SSR、SSTおよび2種類の決定係数を計算します。入力ベクトルの長さが異なる場合、またはSSTが0の場合はエラーになります。

- **戻り値**: `c(SSE, SSR, SST, R2_1, R2_2)`

<br>

### 4. 線形回帰モデル

#### `param_estim_by_LM(xdata, ydata)`

日付で対応づけた2系列に対し、`log10(y) ~ log10(x)` の単回帰モデルを当てはめます。予測区間の計算に必要な統計量も返します。

- **戻り値**: `c(intercept, slope, t_value, n, x_mean, Sxx, residual_variance, SSE, SSR, SST, R2_1, R2_2)`

<br>

#### `epi.prediction_by_LM(xdata, pa, pb, tval, num.d, xmean, sxx, uv, rsq)`

`param_estim_by_LM()` が返す回帰パラメータと統計量を使い、`xdata` に対する対数10スケールの予測値、95%信頼区間、95%予測区間を計算します。

- **戻り値**: `date`, `xdata`, `prediction.y`, `ci.up`, `ci.lw`, `pi.up`, `pi.lw` のデータフレーム
- 注: 現在の実装では引数 `rsq` は計算に使用されません。

<br>

### 5. Shedding profileモデル

#### `param_estim_by_SPM(epidemicdata, sewagedata)`

ガウス型の排出プロファイルを用い、疫学データ（I）から下水ウイルス量（V）を予測する順方向モデルのパラメータを推定します。週次の疫学データは必要に応じて日次へ補間されます。

- `epidemicdata`: 日付と感染者数の2列からなるデータ。週次・日次のどちらにも対応します。
- `sewagedata`: 日付とウイルス濃度の2列からなるデータ
- **戻り値**: 次の2要素を持つリスト
  - `params`: `v`, `m`, `sigma`, `ci_width`, `pi_width`, `R2`
  - `predicted`: `date`, `pred.virus` を持つ推定期間の日次予測値

<br>

#### `virus_prediction_by_SPM(epidemicdata, v, m, sigma)`

指定済みの排出プロファイルパラメータで、感染者数から日次の下水ウイルス量を順方向に予測します。返される予測値は対数10スケールです。

- `v`: 排出プロファイルの大きさ
- `m`: 排出ピークの時間的ずれ
- `sigma`: 排出プロファイルの幅
- **戻り値**: `date`, `pred.virus` のデータフレーム

<br>

#### `infected_prediction_by_SPM(sewagedata, est_v, est_m, est_sigma, lambda = 10^2)`

下水ウイルス量（V）から感染者数（I）を逆方向に推定します。排出プロファイルから設計行列を作り、対数スケールの下水ウイルス量の再現誤差と、対数スケールの感染者数系列に対する2階差分ペナルティを最小化します。感染者数は常に正となるよう推定され、返却時に有症率 `2/3` が掛けられます。

- `sewagedata`: `date`, `virus` 列を持つ下水データ
- `est_v`, `est_m`, `est_sigma`: `param_estim_by_SPM()` で得られる排出プロファイルパラメータ
- `lambda`: 対数感染者数系列の滑らかさを制御するペナルティ係数（既定値: `10^2`）。大きいほど推定系列は滑らかになります。
- **戻り値**: 次の2要素を持つリスト
  - `resI`: `date`, `estimated_cases` を持つ推定感染者数のデータフレーム
  - `resV`: 推定感染者数から再構成した下水ウイルス量の数値ベクトル

## 基本的な利用例

以下は、疫学データから下水ウイルス量を予測し、逆に下水データから感染者数を推定する流れの例です。

```r
source("EpiPredFunctions20260819.R")

# epidemicdata: date, infected
# sewagedata:   date, virus
fit <- param_estim_by_SPM(epidemicdata, sewagedata)

# 推定済みパラメータでウイルス量を順方向予測
virus_pred <- virus_prediction_by_SPM(
  epidemicdata,
  v = fit$params$v,
  m = fit$params$m,
  sigma = fit$params$sigma
)

# 下水ウイルス量から感染者数を逆方向推定
infected_fit <- infected_prediction_by_SPM(
  sewagedata,
  est_v = fit$params$v,
  est_m = fit$params$m,
  est_sigma = fit$params$sigma
)

# 推定感染者数
infected_pred <- infected_fit$resI
```

## 注意事項

- モンテカルロ法を使う発症日推定（`generate_onset()`、`generate_onset_sentinel()`）の結果には乱数によるばらつきがあります。再現可能な結果が必要な場合は、実行前に `set.seed()` を指定してください。
- 関数内で列名を変更するものがあります。入力データの列順と日付型を事前に確認してください。
- CI・PIおよび決定係数は、各関数で用いる対数スケール／元スケールに合わせて解釈してください。
- 本コードは研究・解析用途の関数集です。推定結果の疫学的解釈や実運用への利用は、対象データの品質、仮定、外部妥当性を確認したうえで行ってください。

## ライセンス

本ソースコードは [MIT License](LICENSE) のもとで公開されています。詳細は `LICENSE` ファイルを参照してください。
