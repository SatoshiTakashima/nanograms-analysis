# run12_test metadata config reference

このディレクトリの `config*.yaml` は、run12_test の解析条件をまとめた設定ファイルである。対象は次の3つ。

- `config_testpulse_fit.yaml`
- `config_interpolate_gain.yaml`
- `config_pipeline.yaml`

Python スクリプトで読む YAML の相対パスは、基本的に「その YAML ファイルが置かれているディレクトリ」を基準に解釈される。ComptonSoft/Ruby 側へ渡す `config_pipeline.yaml` は、run12_test ディレクトリで実行する前提のパスと、ComptonSoft 側で config ファイルの場所を基準に解釈するパスが混ざる。混乱しそうな場合は絶対パスを書くのが安全。

数値範囲については、コードが明示的にエラーを出す範囲と、解析上そう書くべき範囲がある。下では「有効な値」として両方をまとめて書く。

## config_testpulse_fit.yaml

テストパルス ROOT ファイルを読み、各 FEC のテストパルスピークをフィットして、ゲイン補正用 CSV を作る設定である。主に `run12_test/testpulse_analysis/run_gui_fitting.py` と `mymodule/analyze_testpulse_gain.py` が使う。

### input

`input.data_root`

: VATA/test-pulse データが置かれているルートディレクトリ。
有効な値は、存在するディレクトリへのパス。絶対パスまたは `config_testpulse_fit.yaml` からの相対パスを書ける。

`input.file_name`

: 各時刻ディレクトリ内で読む ROOT ファイル名。
有効な値は文字列。通常は `outfile00001_000.root` のような ROOT ファイル名。

`input.tree`

: ROOT ファイル内の TTree 名。
有効な値は文字列。現在の VATA データでは `vatatree`。

`input.time_ids`

: 解析するテストパルスデータの時刻 ID リスト。
有効な値は文字列のリスト。現在の推奨形式は `YYYYMMDD/HHMM_SS`。

例:

```yaml
time_ids:
  - "20251001/0124_04"
  - "20251001/0132_24"
```

日付をまたぐ場合も同じリストに追加する。`date_dirs` と短い時刻 ID を組み合わせる方式は run12_test では使わない。

### test_pulse

`test_pulse.channel`

: テストパルス信号を見る VATA channel 番号。
有効な値は整数。通常の VATA channel なら `0` から `63`。現在は `17`。

`test_pulse.ccal`

: テストパルスの校正容量に対応する値。
有効な値は正の数値。ComptonSoft 側のエネルギー較正でも同じ値を使うため、`config_pipeline.yaml` の `calibration.energy.ccal` と揃えるのが基本。

`test_pulse.fec_ids`

: フィット対象にする FEC 番号リスト。
有効な値は `0`, `1`, `2`, `3` の整数からなるリスト。通常は `[0, 1, 2, 3]`。

### fit

`fit.histogram_min_adu`

: フルレンジヒストグラムの ADU 下限。
有効な値は数値。`fit.histogram_max_adu` より小さくする。

`fit.histogram_max_adu`

: フルレンジヒストグラムの ADU 上限。
有効な値は数値。`fit.histogram_min_adu` より大きくする。

`fit.histogram_bin_width_adu`

: ヒストグラムの bin 幅。
有効な値は正の数値。単位は ADU。細かくするとピーク形状は見やすくなるが、bin ごとの統計は減る。

`fit.fit_window_adu`

: ピーク近傍をフィットするときの初期ウィンドウ半幅。
有効な値は正の数値。単位は ADU。GUI では FEC ごとに調整できる。

### quality

`quality.*` はガンマ線イベントの selection には直接使わない。テストパルスピークのフィット結果を、ゲイン補正値として採用してよいか判定するための条件である。判定結果は `fits_csv` の `alive` や `reason` に保存され、manual selection を使わない場合は `summary_csv` へ書く値の採用判定にも使われる。

`quality.min_events`

: その時刻・FEC に必要な最小イベント数。
有効な値は `0` 以上の整数。これ未満なら `too_few_events`。

`quality.min_fit_entries`

: フィット範囲内に必要な最小エントリ数。
有効な値は `0` 以上の整数。これ未満なら `too_few_fit_entries`。

`quality.min_peak_count`

: ピーク bin に必要な最小カウント。
有効な値は `0` 以上の数値。これ未満なら `weak_peak`。

`quality.peak_min_adu`

: 採用するピーク位置の ADU 下限。
有効な値は数値。`quality.peak_max_adu` 以下にする。範囲外なら `peak_out_of_range`。

`quality.peak_max_adu`

: 採用するピーク位置の ADU 上限。
有効な値は数値。`quality.peak_min_adu` 以上にする。範囲外なら `peak_out_of_range`。

`quality.sigma_min_adu`

: 採用するガウス幅 sigma の下限。
有効な値は `0` より大きい数値。`quality.sigma_max_adu` 以下にする。iminuit の sigma 下限にも使われる。

`quality.sigma_max_adu`

: 採用するガウス幅 sigma の上限。
有効な値は `quality.sigma_min_adu` 以上の数値。iminuit の sigma 上限にも使われる。

`quality.max_chi2_ndf`

: 採用するフィットの `chi2/ndf` 上限。
有効な値は `0` 以上の数値。これを超えると `large_chi2`。

### output

`output.outdir`

: テストパルス解析の出力ディレクトリ。
有効な値はパス文字列。相対パスの場合は YAML ファイルの場所を基準にする。現在は `../products` なので `run12_test/products` に出力される。

`output.fits_csv`

: 全フィット結果を保存する CSV ファイル名。
有効な値は文字列。`output.outdir` の下に作られる。

`output.summary_csv`

: ゲイン補正に使う要約 CSV ファイル名。
有効な値は文字列。通常は `testpulse_data.csv`。形式は `time_id,datetime,FEC0,FEC1,FEC2,FEC3`。

`output.plot_dir`

: バッチ処理で作る確認用プロットの保存先ディレクトリ名。
有効な値は文字列。`output.outdir` の下に作られる。

`output.make_plots`

: バッチ処理で確認用プロットを作るかどうか。
有効な値は boolean。YAML では `true` または `false`。

### review

`review.use_manual_selection`

: 手動選択 CSV を最終採用判定に使うかどうか。
有効な値は boolean。`true` なら `manual_selection.csv` の `use_fit` が最終判断になる。`false` なら `quality.*` で決まる `alive` が最終判断になる。

`review.selection_csv`

: 手動選択結果を保存する CSV ファイル名。
有効な値は文字列。`output.outdir` の下に作られる。

`review.plot_dir`

: GUI/レビュー用プロットの保存先ディレクトリ名。
有効な値は文字列。`output.outdir` の下に作られる。

## config_interpolate_gain.yaml

テストパルスから得たピーク位置を、ガンマ線データ取得時刻に対して時間内挿し、ガンマ線解析用のゲイン補正 CSV を作る設定である。`run12_test/testpulse_analysis/interpolate_gain_for_gamma.py` が使う。

### input

`input.testpulse_csv`

: `config_testpulse_fit.yaml` の `output.summary_csv` で作ったテストパルスゲイン時系列 CSV。
有効な値は CSV ファイルへのパス。相対パスの場合は YAML ファイルの場所を基準にする。

`input.gamma_csv`

: ガンマ線測定の時刻リストが書かれた CSV。
有効な値は CSV ファイルへのパス。`time` という列を含む必要がある。`time` 列の値は `YYYYMMDD/HHMM_SS` または `YYYYMMDD_HHMM_SS`。

`input.fec_ids`

: 内挿対象にする FEC 番号。
有効な値は `0`, `1`, `2`, `3` の整数からなるリスト。通常は `[0, 1, 2, 3]`。

### output

`output.csv`

: ガンマ線データ時刻に対して内挿したゲイン補正値を書き出す CSV。
有効な値は CSV ファイルへのパス。相対パスの場合は YAML ファイルの場所を基準にする。

出力 CSV の `time` は `YYYYMMDD/HHMM_SS` に正規化される。各 FEC について、ガンマ線データの時刻を挟む前後のテストパルスゲインがある場合だけ値を書く。前または後のどちらかがない場合、その FEC は空欄になる。

## config_pipeline.yaml

TPC tree から quicklook tree と hittree を作る ComptonSoft パイプラインの設定である。`run12_test/clustering_tpcdata.rb` から `NanoGRAMSHitExtraction` と `NanoGRAMSCalibration` に渡される。

### general

`general.efield_v_cm`

: 電場。
有効な値は正の整数。単位は V/cm。C++ 側では整数として読む。

`general.temperature_k`

: 液体アルゴン温度。
有効な値は正の数値。単位は K。

### light

`light.event_selection_mode`

: 光信号をイベント選別にどう使うかを指定する。
有効な文字列は `gamma_required`, `veto_only`, `disabled` の3つだけ。似た意味の別名は受け付けない。

- `gamma_required`: 電荷クラスタリングで作った候補に加えて、光信号にもガンマ線候補らしいピークを要求する。`light.general_analysis_channels` の波形を解析し、trigger 周りの ROI 内ピークが `light.light_gamma_thr_mV` を超えることを要求する。さらに cosmic/pileup 的な光信号があれば落とす。光信号が取れていないイベントは候補から外れる。
- `veto_only`: ガンマ線候補かどうかは基本的に電荷クラスタリングで決める。光信号は「落とすための条件」にだけ使う。ROI 内ピークが `light.light_gamma_thr_mV` を超えることは要求しない。一方で、`light.light_cosmic_thr_mV` を超える大きな光ピークや、ROI 外に `light.out_roi_peak_thr_mV` を超えるピークがあるイベントは cosmic/pileup 的として落とす。光信号がないイベントでも、veto 条件に引っかからなければ電荷側の判定で残る。
- `disabled`: 光信号によるイベント選別を完全に使わない。ガンマ線候補要求も cosmic/pileup veto も行わず、電荷クラスタリング側の条件だけで判断する。

`light.waveform_analysis`

: 複数の光波形 channel をどう扱うかを指定する。
有効な文字列は `average`, `each_channel` の2つだけ。似た意味の別名は受け付けない。

- `average`: `light.general_analysis_channels` や `light.pileup_analysis_channels` に指定した複数 channel の波形を、同じ時間軸にそろえて平均し、その平均波形からピークを探す。ランダムノイズを下げられる一方、同じ解析グループ内の channel は同じ `wave_compress` でなければならない。channel ごとの応答差は平均される。
- `each_channel`: 指定した channel を1本ずつ個別に解析し、見つかったピーク情報を統合する。channel ごとの `wave_compress` が違っても扱いやすく、特定 channel だけが強い場合にも拾いやすい。一方で平均によるノイズ低減はしない。

`average` を使う場合、同じ解析グループ内の channel は同じ `wave_compress` である必要がある。

`light.general_analysis_channels`

: 通常の光信号解析に使う DPP channel リスト。
有効な値は `0` から `7` の整数リスト。空リストは実質的に通常光解析なしになる。

`light.pileup_analysis_channels`

: pileup 判定に使う DPP channel リスト。
有効な値は `0` から `7` の整数リスト。空リストは実質的に pileup 光解析なしになる。

`light.light_gamma_thr_mV`

: ガンマ線候補として扱う光信号ピークのしきい値。
有効な値は `0` 以上の数値。単位は mV。

`light.light_cosmic_thr_mV`

: cosmic/veto 判定に使う光信号ピークのしきい値。
有効な値は `0` 以上の数値。単位は mV。

`light.pre_roi_window_us`

: ROI の前側に取る時間幅。
有効な値は `0` 以上の数値。単位は us。

`light.post_roi_window_us`

: ROI の後側に取る時間幅。
有効な値は `0` 以上の数値。単位は us。

`light.out_roi_peak_thr_mV`

: ROI 外に大きなピークがあるかを見るためのしきい値。
有効な値は `0` 以上の数値。単位は mV。pileup や veto 的な判定に使う。

`light.delay_counts` は使わない。光波形の trigger delay は、各 `tpctree.root` と同じディレクトリにある `config_dpp.yaml` の `savefile.listwave_delay.value` から DPP channel ごとに読む。

例:

```yaml
savefile:
  listwave_delay:
    value: [60, 60, 60, 60, 60, 60, 60, 60]
```

有効な値は8要素の整数リスト。各要素が DPP channel 0-7 に対応する。値は `0` 以上で、単位は DPP config 上の delay count。ComptonSoft 内では、channel ごとに `delay * 8 * wave_compress` 相当の時間として ROI の位置を決める。

### charge

`charge.adu_range`

: 電荷信号として採用する `ADU - CMN` の範囲。
有効な値は2要素の数値リスト `[min, max]`。`min < max` にする。

`charge.clustering_pix_range`

: クラスタとして採用する pixel 数の範囲。
有効な値は2要素の整数リスト `[min, max]`。`0 <= min <= max <= 64` が自然な範囲。

`charge.circ_thr`

: クラスタ形状の circularity に関するしきい値。
有効な値は `0` 以上の数値。

`charge.spread_thr`

: クラスタの広がりに関するしきい値。
有効な値は `0` 以上の数値。

`charge.drift_time_max_us`

: 採用する最大ドリフト時間。
有効な値は正の数値。単位は us。この値より遅い hit は解析から外れる。

`charge.noise_th`

: ノイズ判定に使うしきい値。
有効な値は `0` 以上の数値。小さくすると弱い信号を拾いやすくなるが、ノイズも増えやすくなる。

`charge.circ_min_hits`

: circularity 判定を行うために必要な最小 hit 数。
有効な値は `0` 以上の整数。

`charge.core_exclude_pix`

: FEC ごとに解析から除外する pixel を指定する。
有効なキーは FEC 番号 `0`, `1`, `2`, `3`。有効な値は pixel 番号 `0` から `63` の整数、整数リスト、または文字列 `peripheral`。似た意味の別名は受け付けない。

例:

```yaml
core_exclude_pix:
  0: [0, 63]
  1: [0, 63]
  2: [0, 63]
  3: [0, 63]
```

`peripheral` を指定すると、その FEC の周辺 pixel をまとめて除外する。整数リストの中に `peripheral` を混ぜることもできる。

`charge.cross_fec_merge_drift_time_tolerance_us`

: config に書いた場合だけ有効になる任意項目。異なる FEC にまたがるクラスタを drift time 差でまとめる許容幅。
有効な値は `0` 以上の数値。単位は us。書かなければ無効。

### calibration

`calibration.energy.gain_info_file`

: FEC/channel の gain パラメータが入った HDF5 ファイル。
有効な値はファイルパス。相対パスの場合は `config_pipeline.yaml` の場所を基準に ComptonSoft 側で解釈される。

`calibration.energy.q_to_kev_spline_file`

: 電荷量から keV へ変換する spline が入った ROOT ファイル。
有効な値はファイルパス。相対パスの場合は `config_pipeline.yaml` の場所を基準に ComptonSoft 側で解釈される。

`calibration.energy.max_time_us`

: エネルギー較正で使う最大 drift time。
有効な値は正の数値。単位は us。

`calibration.energy.tp_channel`

: テストパルス補正に使う channel 番号。
有効な値は整数。通常の VATA channel なら `0` から `63`。`config_testpulse_fit.yaml` の `test_pulse.channel` と揃えるのが基本。

`calibration.energy.ccal`

: テストパルスの校正容量に対応する値。
有効な値は正の数値。`config_testpulse_fit.yaml` の `test_pulse.ccal` と揃えるのが基本。

`calibration.position.anode_pos_z_cm`

: アノード面の z 位置。
有効な値は数値。単位は cm。再構成される hit 位置の z 座標に効く。

テストパルスゲイン補正を再計算する時間間隔は、`config_pipeline.yaml` ではなく Ruby 側の `NanoGRAMSCalibration` module parameter `gain_cache_seconds` で指定する。`run12_test/clustering_tpcdata.rb` では `gain_cache_seconds = 60.0` として `with_parameters` に渡している。単位は秒。何も渡さない場合の ComptonSoft 側デフォルトは `60.0`。`0` 以下にすると Unix 秒単位で更新する。

## よく触る項目

テストパルスの入力を変えるとき:

- `config_testpulse_fit.yaml` の `input.data_root`
- `config_testpulse_fit.yaml` の `input.time_ids`
- `config_testpulse_fit.yaml` の `input.file_name`

フィットの見え方や採用条件を変えるとき:

- `config_testpulse_fit.yaml` の `fit.*`
- `config_testpulse_fit.yaml` の `quality.*`
- `config_testpulse_fit.yaml` の `review.*`

ガンマ線データに対応するゲイン補正 CSV を作るとき:

- `config_interpolate_gain.yaml` の `input.gamma_csv`
- `config_interpolate_gain.yaml` の `output.csv`

TPC hit extraction の選別条件を変えるとき:

- `config_pipeline.yaml` の `light.*`
- `config_pipeline.yaml` の `charge.*`

エネルギー較正や位置較正を変えるとき:

- `config_pipeline.yaml` の `calibration.energy.*`
- `config_pipeline.yaml` の `calibration.position.*`
