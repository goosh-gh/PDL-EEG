# PDL::EEG::IO::NihonKohden

日本光電 EEG `.EEG`（`wfmblock` = EEG-1100 系 / `extblock` = EEG-1200 系、
マルチセグメント含む）を直接 PDL piddle として読む。これは動作している実装の
リファレンス。

このリポジトリの日本光電まわりの doc:

- **このファイル** — モジュールの API と挙動（返却契約・オプション・`.21E` 処理・
  DC 命名・`_ref` 正規化）
- [`NK_block_layouts.md`](NK_block_layouts.md) / [`.ja`](NK_block_layouts.ja.md)
  — `wfmblock` と `extblock` の対比・署名分岐・バイトレイアウト・修正済みの
  1レコードずれ（横断リファレンス）
- [`NK_extblock.md`](NK_extblock.md) — extblock（EEG-1200A/C）の深掘り仕様
- [`NK_files.md`](NK_files.md) — 1記録が書き出す `.EEG`/`.PNT`/`.21E`/`.LOG`/
  `.CN3`/`.PTN`/… サイドカー群の内容

## API

### `read_nk($eeg_path, %opts) -> \%rec`

読み込むファイルの署名からレイアウトを判定し（[分岐の詳細](NK_block_layouts.md)）、
`wfmblock` / `extblock` いずれの経路でも同じ返却契約の hashref を返す。

**オプション**（すべて任意）:

| オプション | 意味 |
|-----------|------|
| `block => N` | `wfmblock` の N 番目の波形ブロックを読む（既定 0）|
| `all_blocks => 1` | 全ブロックを1記録に連結 |
| `dc_base => 1\|3` | フロントパネル上でハード code 45 を何と呼ぶか。`1` → DC01…（EEG-1200 系）/ `3` → DC03…（EEG-1100 系）。世代判定と食い違えば croak |
| `no_events => 1` | `.LOG` があってもイベント解析を省略 |
| `label_map => \%h` | ラベルを上書き。キーは 1-based ch_idx（例 `{45=>'DC03'}`）|
| `fs => N` | サンプリング周波数を明示指定（ヘッダから取れない時）|
| `gap_samples => N` | DEPRECATED |

**返却契約**（両レイアウト共通のキー）:

| キー | 内容 |
|------|------|
| `data` | `[n_ch, n_samp]` float32、**µV**（DC は µV 換算後、STIM は生値）|
| `fs` | サンプリング周波数 (Hz) |
| `labels` | チャンネル名の配列 |
| `t_start` | 記録開始時刻文字列 |
| `events` | `[{t=>秒, label=>...}, …]`（`.LOG` 由来）|
| `gains` | `[n_ch]` µV/bit |
| `ch_indices` | ラベルごとの 1-based 電極コード |
| `n_blocks` / `n_ch_valid` / `block_idx` / `all_blocks` | 構成情報 |
| `device` | 署名（例 `"EEG-1100C V01.00"`）|
| `layout` | `'wfmblock'` \| `'extblock'` |
| `system_reference` | 収録基準（例 `"C3,C4"`、`.21E [SYSTEM_SETUP]`）|
| `last_pattern` | 収録時表示モンタージュ番号（`.21E [LASTPATTERN]`）|
| `gap_bounds` / `t_block_starts` / `block_meta` / `n_samp_per_block` | セグメント境界情報 |
| `units` | チャンネル毎の単位。`wfmblock` は全 `'uV'`、`extblock` は `'uV'`\|`'mV'`(DC)\|`'code'`(STIM) |

`extblock` は加えて `ch_hw_idx`（=`ch_indices`）、`stim_index`、`hdr_len` を返す。

### `clock_to_samp($rec, $wall_sec) -> $samp`

壁時計秒を `block_meta` の区分線形マップでデータサンプルに変換（記録中断を跨いでも
正しい）。ギャップに落ちる時刻は直前の実サンプルにクランプ。

### `nk_layout($path) -> ($sig, $layout, $how)`

署名からレイアウトを確定。判定経路は [`NK_block_layouts.md`](NK_block_layouts.md)。

### `nk_format_hint($model) -> ($sig, $layout, $note)`

機種名からの**非権威**な当たり付け（UX 用）。機種名はファイルに残らないので、
確定情報は常に `nk_layout`。

## `.21E`（電極名）処理

- キーは4桁ゼロ詰め = 0-based 電極インデックス。`[ELECTRODE]` を主に読む。
- `[SD_DEF]` は同じ数値キーを持つ montage/感度マトリクスなので、素朴に全セクション
  を読むと電極名を上書きしてしまう ── `[ELECTRODE]`（と `[REFERENCE]`）だけを見る。
- `[REFERENCE]` フォールバック: `[ELECTRODE]` に無く既定名も空（`-`）のインデックス
  だけ `[REFERENCE]` で補完。`[ELECTRODE]` 優先は不変。
- `$` 正規化: NK の `$` 付き参照名（`$A1`）は Perl／ファイル名で危険かつ電極名 `A1`
  と衝突するため、`A1_ref` の接尾辞形に正規化して返す。

電極コード `c`（1-based, `ch_indices`）→ 名前は `.21E[c-1]`。これは MNE-Python
`read_raw_nihon`（生バイトを 0-based テーブルにそのまま当てる）と一致する。詳細と
`wfmblock` の1レコードずれ修正は [`NK_block_layouts.md`](NK_block_layouts.md)。

## DC トリガ番号の世代差

DC チャンネル番号は**署名**で決まる（レイアウトではない）:

- `wfmblock` / EEG-1100 → `DC03`–`DC06`
- `extblock` / EEG-1200 → `DC01`–`DC04`

実験トリガは DC チャンネルに TTL レベルで入る。波形末尾の STIM ch にも code が出るが
これは記録マークで実験トリガではない。

## 関連モジュール

- `PDL::EEG::Derivation` — `bne()`（balanced non-cephalic 再基準）。収録基準
  `Avr(C3,C4)` は再基準の重みが 1 に和すると相殺するので、収録データから直接計算
  できる。詳細は [`NK_files.md`](NK_files.md) の参照・`-BN` 節。
- `PDL::EEG::IO::BESA::ASCII` — `write_mul()`（BESA ASCII multiplexed `.mul` 書出）。

## ライセンス

フォーマット解析は EDFbrowser `nk2edf`（GPL-2）および Brainstorm `in_fopen_nk.m`
を**仕様参照**したクリーンルーム実装（コード非コピー）。
