# Nihon Kohden .EEG フォーマット — レイアウト分岐（署名ベース）

`.EEG` の内部レイアウトは **ファイル先頭(@0x0000, 16char)の形式署名だけ**で決まる。
**物理機種名とは別物**。実機 EEG-1290（ヘッドボックス JE-92NX）で記録しても
署名は `EEG-1200A V01.00`。

実データ `JJ0090J6.EEG`（EEG-1290 / MMN, 38ch, 1000Hz）で確定。
参照: Brainstorm `in_fopen_nk.m`, NK2EDF。仕様参照の独自実装（コード非コピー）。
実装は `PDL::EEG::IO::NihonKohden`（`_read_extblock` ＋ 既存 wfmblock 経路）。

## 署名 → レイアウト

| レイアウト | 署名例 | ch情報の場所 |
|-----------|--------|-------------|
| `wfmblock`（旧） | `EEG-1100C V01.00`, `EEG-2100 V01.00`, `QI-403A V01.00`, `DAE-2100D V01.30` … | wfmblock 内（既存 read_nk）|
| `extblock`（新） | `EEG-1200A V01.00` | 拡張ブロック連鎖（ext_address 経由）|

- 分岐は署名文字列を `%FORMAT_LAYOUT` で引く。型番の新旧では判断しない
  （古い型番が extblock を書く／新しい型番が wfmblock を書くこともあり得る）。
- 物理機種名（EEG-1290 等）はファイルに残らない。署名だけが手掛かり。

## 共通ヘッダ（両レイアウト同じ）

| offset | 内容 |
|--------|------|
| 0x0000 | 形式署名 16 char |
| 0x0091 | control block 数 (u8) = 1 |
| 0x0092 | ctlblock_address (u32) = 0x0400 |
| 0x03EE | ext_address (u32) — wfmblock は 0、extblock は非0（例 0x27CF）|
| ctl+0x12 | data_block_address (u32) = 0x17FE |
| data+0x1A | u16 下位14bit = fs（例 0x3E8 = 1000 Hz）|

`ext_address` が 0 か非0 かも実質的な判別材料になるが、正式には署名テーブルで分岐。

## extblock レイアウトのチャンネル解決（拡張ブロック連鎖）

```
ext        = u32(@0x03EE)              # 0x27CF
eb2        = u32(ext + 18)             # 0x2BFF
eb3        = u32(eb2 + 20)             # 0x43FB
n_channels = u16(eb3 + 68) + 1         # 37 + 1(STIM) = 38
for i in 0..n_channels-2:
    hw[i]  = u16(eb3 + 72 + i*10) + 1  # 1-based hardware index
rec_address = eb3 + 72 + (n_channels-1)*10   # 波形開始 = 0x45B5
```

（eb3 は wfmblock 同型サブヘッダ: `01 TIME...` に加え eb3+0x14 に ASCII
フルタイムスタンプ `YYYYMMDDhhmmss`、eb3+0x1C に fs=1000 も持つ。）

## 波形データ（両レイアウト共通の符号化）

- 開始: `rec_address`
- 配置: サンプルインターリーブ `[t][ch]`、各値 uint16 LE
- 1 フレーム = `n_channels` × int16（末尾が STIM）
- サンプル数: 非格納 → `n_samp = floor((filesize - rec_address) / n_channels / 2)`
- アナログ ch: offset binary（中心 0x8000=32768）、`µV = (raw - 32768) * gain`
- STIM ch（末尾）: 生値（オフセット無し）

高速読込（PDL）:
```perl
my $u16 = zeroes(ushort, $nch, $n_samp);   # (ch, t)
${ $u16->get_dataref } = substr($bytes, 0, $nch*$n_samp*2);
$u16->upd_data;
my $data = (($u16->double - $off->slice(':,*1')) * $gain->slice(':,*1'))->float;
```

## ゲイン（µV/bit, ハード固定・ファイル非依存）

1-based ハードウェア ch index（hw の値）で決まる（Brainstorm の micro マップに準拠）:

- micro（µV）: index ∈ {1–42, 75, 76, 79–1096} → `(3199.902+3200)/(32767+32768)` = **0.09765624**
- DC/その他（±12V レンジ）: 上記以外 → `(12002.56+12002.9)/(32767+32768)` = **0.36629984**

本ファイル hw = `[1..20, 23..30, 45..48, 75..78, 100]`:
- DC01–04 = index 45–48 → DC ゲイン
- **A1_ref/A2_ref = index 77,78 → 非micro（DC/mV レンジ）扱い**（Brainstorm マップ上、
  43–74・77・78 は非micro）。参照管理chなので通常は解析対象外
- 他は µV ゲイン

## 電極名 `.21E`

- セクション構成: `[ELECTRODE]` `[REFERENCE]` `[SD_DEF]` `[SYSTEM_SETUP]` `[LASTPATTERN]`
- `[ELECTRODE]` を主に読む。`[SD_DEF]` は 1096 行の montage/感度マトリクスで
  同じ数値キーを持つため、素朴に全セクションを読むと電極名を上書きしてしまう（注意）
- キーは 4 桁ゼロ詰め（`0000=Fp1`）＝電極インデックス（0-based）→ 数値化して引く。CRLF 改行
- hw の値 `c` → `.21E[c-1]`
- **`[REFERENCE]` フォールバック**: `[ELECTRODE]` に無く、かつ組込み既定名も `-`（空）の
  インデックスだけ `[REFERENCE]` で補完（`[ELECTRODE]` 優先は不変）。
  `[REFERENCE]` は index 76 から参照ch群を定義: `76=$A1 77=$A2 78/79=$A+ 80/81=$Cz
  82/83=$AV 84/85=$BN 86/87=$Aav …`
- **`$` 正規化**: NK の `$`付き参照名は Perl/ファイル名で危険（`"$A1"` が変数展開）かつ
  電極名 `A1` と衝突するため、`$A1 → A1_ref` の接尾辞形に正規化して返す

本ファイルのラベル（38ch）:
```
Fp1 Fp2 F3 F4 C3 C4 P3 P4 O1 O2 F7 F8 T3 T4 T5 T6 Fz Cz Pz E
A1 A2 vEOG hEOG X1 nose lm rm DC01 DC02 DC03 DC04 BN1 BN2 A1_ref A2_ref COM STIM
```

## トリガ／イベント

- **実験トリガ（MMN 標準/逸脱）は DC チャンネルに TTL レベルで格納**。DC 差し口の番号は
  **フォーマット世代で変わる**:
  - `EEG-1100C`: **DC03–DC06**
  - `EEG-1200A`（本ファイル）: **DC01–DC04**（hw 45–48）
  DC ch はアナログ行としてそのまま raw.plot() に出る（矩形波）。閾値交差で onset を取る。
  `read_nk` は `.21E` の名前を素通しするので、各ファイルは正しい DC 名を返す。
- `.LOG`: セッション注釈のみ（REC START / task1–5 / A1+A2 OFF / 安静開眼閉眼…）。
  日本語は Shift-JIS。1 エントリ 45 byte: `[20:label][2:HH][2:MM][2:SS][19:label2]`。
  時刻は 6 桁 ASCII 秒。
- `.EVT`: ヘッダ `Tmu\tCode\tTriNo` のみで本ファイルは空。
- 波形末尾 STIM ch にも code(1/2/…) が出るが、これは記録マークで実験トリガではない。

## 2層構造: 機種名ヒント（可変）と 署名→レイアウト（確定）

- **機種名 → 記録フォーマットの当たり付け（非権威・UX用）**: `%DEVICE_FORMAT_HINT`。
  EEG-1290→1200A、EEG-1214→1100C。SW更新で 1290 が 1200B/C に、1214 が 1200A に
  移り得るので、これはあくまで「開く前の目安」。`nk_format_hint($model)` が
  (推定署名, 推定レイアウト, 注記) を返す。機種名はファイルに残らない。
- **署名 → レイアウト（ファイルから読む確定情報）**: `nk_layout($path)`。常にこちらが真。

## nk_layout のフォールバック連鎖

未知署名でもハードフェイルさせず、可能なら読む:

1. `%FORMAT_LAYOUT` 完全一致 → 権威（how=`table`）
2. 構造判定: `ext_address(@0x3EE)` 非0 → extblock、0 → wfmblock
   （検証済みの判別子: 1100C=0, 1200A≠0。how=`fallback:ext_address(...)`）
3. 名前系統（`EEG-12*`→extblock / `EEG-11*`→wfmblock）は**照合・警告用のみ**。
   構造と食い違えば構造（ファイルの真実）を採り WARNING を付す
4. `EEG-/QI-/DAE-` で始まらない非NK署名 → `undef`（誤読せず停止）

返り値 `($sig, $layout, $how)`。`$how` に判定経路が入るので、table 以外なら
`NK_DEBUG` で警告表示できる。

## read_nk への組込み（実装済み）

`read_nk` 冒頭で署名分岐:
```perl
my (undef, $layout, $how) = nk_layout($eeg_path);
croak "Unknown Nihon Kohden signature in $eeg_path" unless defined $layout;
return _read_extblock($eeg_path, %opts) if $layout eq 'extblock';
# else: 既存 wfmblock 経路（_read_wfm_header ほか）へフォールスルー
```
`_read_extblock` は既存 read_nk と同じ返却契約
（data[n_ch,n_samp]float32 µV / fs / labels / t_start / events / gains /
n_blocks / n_ch_valid / block_idx / device）＋ extblock 固有情報
（layout / ch_hw_idx / stim_index / t_block_starts / n_samp_per_block）。
共有ヘルパー `_read_bytes/_read_u*`, `_read_21e`, `_read_log`, `_bcd_byte`,
`@DEFAULT_LABELS` を再利用。1200B/C や更新後 1214 はフォールバック 1〜2 で
自動的に extblock と解決され、コード変更なしで読める。
