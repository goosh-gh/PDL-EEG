# PDL::EEG::IO::NihonKohden — 設計メモ

## 目的

日本光電 EEG-1100/2100 系バイナリ (*.eeg) を直接 PDL piddle として読み込む。
EDFbrowserの `nk2edf.cpp` (GPL-2) をリファレンスとしてフォーマット解析。

---

## ファイルセット (EEG-1100/2100 系)

| ファイル | 必須 | 内容 |
|---------|------|------|
| `*.eeg` | ✓    | 波形データ本体 |
| `*.pnt` | EDF+時 | Patient情報ブロック |
| `*.log` | EDF+時 | イベント/アノテーション |
| `*.21e` | オプション | カスタム電極名マッピング |

---

## .eeg バイナリ構造 (nk2edf.cpp より逆引き)

```
Offset      Size   内容
0x0000      16     Device signature ("EEG-1100A V01.00" 等)
0x0081      16     Control block signature
0x0091       1     Control block count (ctl_block_cnt)
0x0092+i*20  4     ctlblock_address[i] (LE uint32)

ctlblock_address + 17:  1     Data block count (datablock_cnt)
ctlblock_address + 18 + j*20: 4  wfmblock_address[j]

0x17FE       1     Waveform data block signature (== 0x01)

--- Waveform block (at wfmblock_address) ---
+0x0006      1     Number of channels (n_channels)
+0x000A      2     Sampling rate code → actual Hz
+0x000E      4     Number of samples per channel
+0x0022     ... (header continues)
+ofs         2*n_samples*n_ch  signed int16 LE サンプルデータ

Gain/物理単位:
  各チャネルに gain コード (1バイト) → μV換算係数テーブルで変換
  標準: 10μV/bit @ gain=0x00, 係数テーブルは nk2edf.cpp の get_chan_sensitivity() 参照
```

### サンプリングレート変換テーブル

| コード | Hz    |
|--------|-------|
| 0xA0   | 100   |
| 0xA1   | 200   |
| 0xA2   | 500   |
| 0xA3   | 1000  |
| 0xA4   | 2000  |

---

## モジュール設計

### PDL::EEG::IO::NihonKohden

```
read_nk($eeg_file, %opts)
  → { data   => $pdl,      # [n_ch, n_samples] float32, μV
      fs      => $hz,
      labels  => \@labels,
      t_start => $datetime_str,
      events  => \@events,   # [{t=>$sec, label=>$str}, ...]
      gains   => $pdl,       # [n_ch] μV/bit
    }
```

### 内部サブルーチン

```perl
_read_header($fh)         # device sig チェック、ctl block アドレス群返す
_read_wfm_block($fh, $addr)  # 1 waveform block の meta + raw int16
_read_21e_labels($path)   # 電極名ハッシュ (ch_idx => name)
_read_log($log_path)      # イベントリスト
_srate_from_code($code)   # サンプリングレート変換
_gain_uv($gain_code)      # μV/bit 係数
_int16_to_uv($raw_pdl, $gains_pdl)  # int16 → μV PDL変換
```

---

## PDL 変換の核心

```perl
# raw データ: [n_ch * n_samples] のバイト列 → PDL int16
my $raw = PDL->new_from_specification(short, $n_ch, $n_samples);
read($fh, ${$raw->get_dataref}, $n_ch * $n_samples * 2);
$raw->upd_data;

# バイトオーダー (LE保証): MacでもLinuxでも
# PDL 5.x: $raw->bswap2 は BE→LE。x86/ARM64 は LE なので通常不要。
# 念のためポータブルに:
use PDL::IO::Misc;
# または手動: unpack 'v*' → PDL->new([...])

# μV変換: ゲイン piddle でブロードキャスト
my $data_uv = $raw->double * $gains->slice(':,*1');  # [n_ch,1] broadcast
```

---

## イベント (LOGファイル)

nk2edf.cpp より:
```
logblock_address + 0x0012:  1   n_logs
logblock_address + 0x0014:  n_logs * 45 bytes
各エントリ (45 bytes):
  [0..19]   event label (Latin-1, スペースパディング)
  [20..21]  時刻 (秒, uint16 LE, 記録開始からのオフセット)
  [22..24]  日付関連
  ... (残りはサブイベント用)
```

---

## CPAN Distribution 構想

```
PDL-EEG-IO-NihonKohden/
├── lib/PDL/EEG/IO/NihonKohden.pm   # メイン
├── lib/PDL/EEG/IO/NihonKohden/
│   ├── Header.pm                   # バイナリヘッダ解析
│   ├── Waveform.pm                 # 波形データ読み込み
│   └── Events.pm                  # LOGファイル解析
├── t/
│   ├── 01_header.t
│   ├── 02_waveform.t               # テスト用合成バイナリ使用
│   └── 03_events.t
├── examples/
│   └── eeg_viewer_basic.pl        # PDL::Graphics::Cairo でのプロット
└── Makefile.PL
```

---

## 将来ビューア (PDL::EEG::Viewer)

```
App-PDL-Notebook 上での動作 OR giza-server Driver::GS 直接使用

$viewer = PDL::EEG::Viewer->new(
    data   => $data_uv,   # [n_ch, n_samples]
    fs     => $hz,
    labels => \@labels,
    events => \@events,
);
$viewer->show;  # giza_server インタラクティブウィンドウ
# → Driver::GS の show_interactive() を内部利用
# → スライダーで時間軸スクロール (GSP_MSG_SLIDER 既存実装を流用)
```

---

## 優先実装順序

1. `_read_header` + `_read_wfm_block` — ファイルが開けて meta が取れる
2. `read_nk` の最小実装 — [n_ch, n_samples] float32 pdl が返る
3. `_read_21e_labels` — 電極名
4. `_read_log` — イベント
5. テスト合成バイナリの作成 (実データなしでt/が通る)
6. ビューア (別ディストリビューション or examples/)

---

## ライセンス注意点

EDFbrowser は GPL-2。フォーマット仕様を**参照**して独自実装する分には問題なし。
コードを直接コピーすれば GPL-2 継承必要。→ **clean-room 実装**を推奨。
