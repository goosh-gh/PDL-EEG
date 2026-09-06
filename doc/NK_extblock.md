# extblock（EEG-1200 系）レイアウト深掘り

`extblock` は EEG-1200 系（署名 `EEG-1200A V01.00` / `EEG-1200C V01.00`）の `.EEG`
内部レイアウト。チャンネル情報が拡張ブロック連鎖にあり、サンプル列が複数の時間
セグメントに分かれうる。ここは `_read_extblock` を触るときの実装リファレンス。

- 2系統（`wfmblock` vs `extblock`）の対比・署名分岐・共通ヘッダ・共通の符号化は
  [`NK_block_layouts.md`](NK_block_layouts.md)。
- モジュール全体の API・返却契約は [`NihonKohden.md`](NihonKohden.md)。

実データ `JJ0090J6.EEG`（EEG-1290 / JE-92NX / MMN, 37ch+STIM, 1000Hz）で確定。
参照: Brainstorm `in_fopen_nk.m`, NK2EDF（仕様参照・コード非コピー）。

## チャンネル解決（拡張ブロック連鎖）

```
ext        = u32(@0x03EE)              # 非0（例 0x27CF）
eb2        = u32(ext + 18)             # 例 0x2BFF
eb3        = u32(eb2 + 20)             # 例 0x43FB（チャンネル情報ブロック先頭）
n_channels = u16(eb3 + 68) + 1         # 例 37 + 1(STIM) = 38
for i in 0..n_channels-2:
    hw[i]  = u16(eb3 + 72 + i*10) + 1  # 1-based 電極コード
rec_address = eb3 + 72 + (n_channels-1)*10   # = eb3 + hdr_len、波形開始
```

`hdr_len = 72 + (n_channels-1)*10`（例 442）。生 `u16` は 0-based の `.21E`
インデックスで、`+1` して 1-based コードにし、共通の `.21E[code-1]` で名前を引く
（`wfmblock` 修正後と同一規約）。eb3 は `wfmblock` 同型サブヘッダ（`01 TIME…`、
`+0x14` に ASCII フルタイムスタンプ、`+0x1C` に fs）を持つ。

## セグメント

`extblock` は単一連続領域ではない。記録中断ごとにチャンネル情報ブロックの複製
（`eb3` のものとタイムスタンプ以外バイト一致、`hdr_len` バイト）がサンプル列に
挿入される。埋め込みヘッダはサンプル境界に整列し識別可能:

```
+0x00        0x01                      ブロック種別
+0x01..+0x04 "TIME"
+0x12,+0x13  0x02 0x02                 バージョン
+0x14..      "YYYYMMDDHHMMSS" (ASCII)  当該セグメント開始時刻
+0x44        u16 = n_channels - 1      チャンネル数（一致必須）
+0x48..      再びチャンネルテーブル
```

`_read_extblock` はこれらを検出して読み飛ばす。旧版は単一ブロックが EOF まで続くと
仮定し、埋め込みヘッダを約 `hdr_len/stride` サンプルの EEG として読み、中断ごとに
チャンネル位相が `hdr_len % stride` バイト（442 % 76 = 62 = 31ch 分）ずれていた。
修正済み。

## 波形データ（両レイアウト共通の符号化）

- 開始: `rec_address`。配置: サンプルインターリーブ `[t][ch]`、各値 uint16 LE。
- 1 フレーム = `n_channels` × int16（末尾が STIM）。
- サンプル数: 非格納 → `floor((filesize − rec_address) / n_channels / 2)`（セグメント
  境界のヘッダ分は除外）。
- アナログ ch: offset binary（中心 0x8000）、`µV = (raw − 32768) × gain`。
- STIM ch（末尾）: 生値（オフセット無し）。

PDL 高速読込:
```perl
my $u16 = zeroes(ushort, $nch, $n_samp);   # (ch, t)
${ $u16->get_dataref } = substr($bytes, 0, $nch*$n_samp*2);
$u16->upd_data;
my $data = (($u16->double - $off->slice(':,*1')) * $gain->slice(':,*1'))->float;
```

## ゲイン（µV/bit・ハード固定・ファイル非依存）

1-based 電極コード `hw` の値で決まる（Brainstorm の micro マップ準拠）:

- micro（µV）: index ∈ {1–42, 75, 76, 79–1096} →
  `(3199.902+3200)/(32767+32768)` = **0.09765624**
- DC/その他（±12V レンジ）: 上記以外 →
  `(12002.56+12002.9)/(32767+32768)` = **0.36629984**

例（`JJ0090J6`）の `hw = [1..20, 23..30, 45..48, 75..78, 100]`:
- `DC01–04` = index 45–48 → DC ゲイン
- `A1_ref/A2_ref` = index 77,78 → 非 micro（43–74・77・78 は非 micro）。参照管理 ch
  で通常は解析対象外
- 他は µV ゲイン

## `.21E`（extblock で見える要点）

`.21E` 一般の扱い（`[ELECTRODE]` 優先、`[SD_DEF]` 上書き回避、`[REFERENCE]`
フォールバック、`$`→`_ref` 正規化）は [`NihonKohden.md`](NihonKohden.md)。この記録の
`[REFERENCE]` は index 76 から参照 ch 群を定義:
`76=$A1 77=$A2 78/79=$A+ 80/81=$Cz 82/83=$AV 84/85=$BN 86/87=$Aav …`。

本ファイルのラベル（38ch）:
```
Fp1 Fp2 F3 F4 C3 C4 P3 P4 O1 O2 F7 F8 T3 T4 T5 T6 Fz Cz Pz E
A1 A2 vEOG hEOG X1 nose lm rm DC01 DC02 DC03 DC04 BN1 BN2 A1_ref A2_ref COM STIM
```

## トリガ／イベント

- 実験トリガ（MMN 標準/逸脱）は DC チャンネルに TTL レベルで格納。DC 差し口番号は
  世代で変わる（EEG-1200A は **DC01–DC04** = hw 45–48。EEG-1100C は DC03–DC06）。
  DC ch はアナログ行として `raw.plot()` に矩形波で出る。閾値交差で onset を取る。
- `.LOG`: セッション注釈のみ（REC START / task1–5 / A1+A2 OFF / 安静開眼閉眼…）。
  日本語は Shift-JIS。1 エントリ 45 byte: `[20:label][2:HH][2:MM][2:SS][19:label2]`、
  時刻は 6 桁 ASCII 秒。
- `.EVT`: ヘッダ `Tmu\tCode\tTriNo` のみ、本ファイルは空。
- 波形末尾 STIM ch の code(1/2/…) は記録マークで実験トリガではない。
